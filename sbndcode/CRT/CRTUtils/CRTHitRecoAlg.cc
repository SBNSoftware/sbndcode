#include "CRTHitRecoAlg.h"

namespace sbnd{

  CRTHitRecoAlg::CRTHitRecoAlg(const fhicl::ParameterSet &p) 
    : fCRTGeoAlg(p.get<fhicl::ParameterSet>("GeoAlg", fhicl::ParameterSet()))
    , fADCThreshold(p.get<uint16_t>("ADCThreshold"))
    , fADCSaturation(p.get<uint16_t>("ADCSaturation"))
    , fQPedestal(p.get<double>("QPedestal"))
    , fQSlope(p.get<double>("QSlope"))
    , fPEAttenuation(p.get<double>("PEAttenuation"))
    , fPropDelay(p.get<double>("PropDelay"))
    , fTimeWalkNorm(p.get<double>("TimeWalkNorm"))
    , fTimeWalkShift(p.get<double>("TimeWalkShift"))
    , fTimeWalkSigma(p.get<double>("TimeWalkSigma"))
    , fTimeWalkOffset(p.get<double>("TimeWalkOffset"))
    , fHitCoincidenceRequirement(p.get<uint32_t>("HitCoincidenceRequirement"))
    , fT1Offset(p.get<uint32_t>("T1Offset", 0))
  {}

  CRTHitRecoAlg::CRTHitRecoAlg() {}

  CRTHitRecoAlg::~CRTHitRecoAlg() {}

  std::map<std::string, std::vector<std::vector<CRTStripHit>>> 
    CRTHitRecoAlg::ProduceStripHits(const std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBdataVec)
  {
    std::map<std::string, std::vector<std::vector<CRTStripHit>>> stripHits;
    // Iterate through each FEBData
    int count=0;
    for(unsigned feb_i = 0; feb_i < FEBdataVec.size(); ++feb_i)
      {
        art::Ptr<sbnd::crt::FEBData> FEBdata = FEBdataVec.at(feb_i);
        uint32_t mac5  = FEBdata->Mac5();
        uint32_t unixs = FEBdata->UnixS();
        CRTModuleGeo module    = fCRTGeoAlg.GetModule(mac5 * 32);
        std::string taggerName = module.taggerName;

        // Resize the vector of vectors to 2, one vector of hits for each
        // of the orientations
        if(stripHits.find(taggerName) == stripHits.end()) stripHits[taggerName].resize(2);

        // Correct for FEB readout cable length
        uint32_t t0 = FEBdata->Ts0() + module.cableDelayCorrection;
        uint32_t t1 = FEBdata->Ts1() + module.cableDelayCorrection;

        // Iterate via strip (2 SiPMs per strip)
        const auto &sipm_adcs = FEBdata->ADC();
        for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
        {
          // Add an offset for the SiPM channel number
          uint16_t channel = mac5 * 32 + adc_i;

          CRTStripGeo strip = fCRTGeoAlg.GetStrip(channel);
          CRTSiPMGeo sipm1  = fCRTGeoAlg.GetSiPM(channel);
          CRTSiPMGeo sipm2  = fCRTGeoAlg.GetSiPM(channel+1);

          // Subtract channel pedestals
          uint16_t adc1 = sipm1.pedestal < sipm_adcs[adc_i]   ? sipm_adcs[adc_i]   - sipm1.pedestal : 0;
          uint16_t adc2 = sipm2.pedestal < sipm_adcs[adc_i+1] ? sipm_adcs[adc_i+1] - sipm2.pedestal : 0;

          uint16_t pedestal1 = sipm1.pedestal;
          uint16_t pedestal2 = sipm2.pedestal;

          // Keep hit if both SiPMs above threshold
          if(adc1 > fADCThreshold && adc2 > fADCThreshold)
          {
            std::cout<<"ADCs: channel1: "<<adc1<<" with pedstals "<<sipm1.pedestal<<" "<< adc1+sipm1.pedestal<<" channel2:"<<adc2<<" with pedstals "<<sipm2.pedestal<<" "<< adc2+sipm2.pedestal<<std::endl;
            // Access width of strip from the geometry algorithm
            // === TO-DO === //
            // AMEND THE CODE AND IMPROVE CALCULATION
            double width = strip.width;
            double x     = width / 2. * tanh(log(1. * adc2/adc1)) + width / 2.;
            double ex    = 2.5;

            // Create hit
            count++;
            stripHits[module.taggerName][module.orientation].emplace_back(channel, t0, t1, unixs, x, ex, adc1, adc2, pedestal1, pedestal2, feb_i);
          }else{
            /*if (adc1>0) {
              stripHits[module.taggerName][module.orientation].back().channels_below_threshold[0]++;
              stripHits[module.taggerName][module.orientation].back().channels_below_threshold[1]+=adc1;
            }
            if (adc2>0) {
              stripHits[module.taggerName][module.orientation].back().channels_below_threshold[0]++;
              stripHits[module.taggerName][module.orientation].back().channels_below_threshold[1]+=adc2;
            }*/
          }
        }
      }
    return stripHits;
  }

  std::vector<sbn::crt::CRTHit> CRTHitRecoAlg::ProduceCRTHits(const std::map<std::string, std::vector<std::vector<CRTStripHit>>> &taggerStripHits)
  {
    std::vector<sbn::crt::CRTHit> crtHits;

    for(auto const &[tagger, stripHitsVec] : taggerStripHits)
    {
      if(stripHitsVec.size() != 2) {
        std::cout << "=== ERROR: Strip hits vector does not have size 2" << std::endl;
        continue;
      }

      // Split hits by orientation, we want to look for coincidences between strips in 
      // opposite orientations
      std::vector<CRTStripHit> hitsOrien0 = stripHitsVec[0];
      std::vector<CRTStripHit> hitsOrien1 = stripHitsVec[1];

      // Get all possible combinations of strips that could make coincident hits
      std::vector<std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit>> candidates = ProduceCRTHitCandidates(tagger, hitsOrien0, hitsOrien1);

      // Order by the hits with the largest reconstructed PE, should be better than random?
      std::sort(candidates.begin(), candidates.end(), 
                [](const std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit> &a, std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit> &b) 
                {
                  return a.second.peshit > b.second.peshit;
                });

      std::set<unsigned> used_i, used_j;
      for(auto const &cand : candidates)
      {
        if(used_i.find(cand.first.first) != used_i.end() || used_j.find(cand.first.second) != used_j.end())
          continue;

        crtHits.push_back(cand.second);
        //used_i.insert(cand.first.first);
        //used_j.insert(cand.first.second);
      }
    }
    return crtHits;
  }

  std::vector<std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit>> CRTHitRecoAlg::ProduceCRTHitCandidates(const std::string &tagger, const std::vector<CRTStripHit> &hitsOrien0, const std::vector<CRTStripHit> &hitsOrien1)
  {
    std::vector<std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit>> candidates;
    for(unsigned i = 0; i < hitsOrien0.size(); ++i)
      {
        const CRTStripHit hit0   = hitsOrien0[i];
        const CRTStripGeo strip0 = fCRTGeoAlg.GetStrip(hit0.channel);
        //const CRTSiPMGeo  hit0sipm0 = fCRTGeoAlg.GetSiPM(strip0.channel0);
        //const CRTSiPMGeo  hit0sipm1 = fCRTGeoAlg.GetSiPM(strip0.channel1);

        for(unsigned j = 0; j < hitsOrien1.size(); ++j)
          {
            const CRTStripHit hit1   = hitsOrien1[j];
            const CRTStripGeo strip1 = fCRTGeoAlg.GetStrip(hit1.channel);
            //const CRTSiPMGeo  hit1sipm0 = fCRTGeoAlg.GetSiPM(strip1.channel0);
            //const CRTSiPMGeo  hit1sipm1 = fCRTGeoAlg.GetSiPM(strip1.channel1);

            // Check whether the two strips responsible for these hits overlap.
            if(!fCRTGeoAlg.CheckOverlap(strip0, strip1)) continue;

            // Find overlap region between two strip hits
            std::vector<double> overlap = FindOverlap(hit0, hit1, strip0, strip1);
            // Using overlap region find centre and 'error'
            TVector3 pos, err;
            CentralPosition(overlap, pos, err);             

            // Just perform the correction to the ADC values, no PE reconstruction
            const std::array<uint16_t, 4> adc = {hit0.adc1, hit0.adc2, hit1.adc1, hit1.adc2};
            const uint16_t raw_adc0 = hit0.adc1 + hit0.pedestal1; //hit0sipm0.pedestal;
            const uint16_t raw_adc1 = hit0.adc2 + hit0.pedestal2; //hit0sipm1.pedestal;
            const uint16_t raw_adc2 = hit1.adc1 + hit1.pedestal1; //hit1sipm0.pedestal;
            const uint16_t raw_adc3 = hit1.adc2 + hit1.pedestal2; //hit1sipm1.pedestal;
            const std::array<uint16_t, 4> raw_adc = {raw_adc0, raw_adc1, raw_adc2, raw_adc3};
            std::array<uint16_t, 4> corr_adc  = {0, 0, 0, 0};
            CorrectADC(pos, hit0, hit1, corr_adc);

            // Reconstruct the number of photoelectrons from the ADC values
            double pe0, pe1;
            double uncorrected_pe0, uncorrected_pe1;
            double distance0, distance1;
            ReconstructPE(pos, hit0, hit1, pe0, pe1, uncorrected_pe0, uncorrected_pe1, distance0, distance1);


            // Correct timings to find how coincident the hits were
            uint32_t t0, t1;
            double diff;
            CorrectTimings(pos, hit0, hit1, pe0, pe1, t0, t1, diff);
            mf::LogInfo("CRTHitRecoAlg") << "at position: " << pos.X() << ", " 
                                         << pos.Y() << ", " << pos.Z() << "cm\n"
                                         << "diff: " << diff 
                                         << std::endl;
            if(std::abs(diff) > fHitCoincidenceRequirement) {
              continue;
            }
            //if(abs((int)hit0.s - (int)hit1.s) > 1) continue;
            const uint64_t unixs = std::min(hit0.s, hit1.s);
            sbn::crt::CRTHit crtHit({(uint8_t)hit0.febdataindex, (uint8_t)hit1.febdataindex},
                                    pe0+pe1,
                                    t0,
                                    (double)t1 - fT1Offset,
                                    diff,
                                    unixs,
                                    pos,
                                    err,
                                    tagger,
                                    hit0.channel,
                                    hit1.channel,
                                    raw_adc,
                                    adc,
                                    corr_adc,
                                    {(float)distance0, (float)uncorrected_pe0, (float)pe0},
                                    {(float)distance1, (float)uncorrected_pe1, (float)pe1});

            mf::LogInfo("CRTHitRecoAlg") << "\nCreating CRTHit "
                                         << "from FEBs: " << (unsigned) crtHit.feb_id[0] 
                                         << " & " << (unsigned) crtHit.feb_id[1] << '\n'
                                         << "at position: " << pos.X() << ", " 
                                         << pos.Y() << ", " << pos.Z() << "cm\n"
                                         << "and t1: " << t1 << '\n'
                                         << "with PE: " << pe0+pe1 << '\n'
                                         << "on tagger: " << tagger << std::endl;

            // Record which strip hits were used to make this candidate
            candidates.push_back({{i, j}, crtHit});
          }
      }
    return candidates;
  }

  std::vector<double> CRTHitRecoAlg::FindOverlap(const CRTStripHit &hit0, const CRTStripHit &hit1,
                                                 const CRTStripGeo &strip0, const CRTStripGeo &strip1)
  {
    const std::vector<double> hit0pos = fCRTGeoAlg.StripHit3DPos(strip0.name, hit0.x, hit0.ex);
    const std::vector<double> hit1pos = fCRTGeoAlg.StripHit3DPos(strip1.name, hit1.x, hit1.ex);

    std::vector<double> overlap({std::max(hit0pos[0], hit1pos[0]), // x max
                                 std::min(hit0pos[1], hit1pos[1]), // x min
                                 std::max(hit0pos[2], hit1pos[2]),
                                 std::min(hit0pos[3], hit1pos[3]),
                                 std::max(hit0pos[4], hit1pos[4]),
                                 std::min(hit0pos[5], hit1pos[5])});

    return overlap;
  }

  void CRTHitRecoAlg::CentralPosition(const std::vector<double> overlap, 
                                      TVector3 &pos, TVector3 &err)
  {
    pos = TVector3((overlap[0] + overlap[1])/2.,
                   (overlap[2] + overlap[3])/2.,
                   (overlap[4] + overlap[5])/2.);
    
    err = TVector3(std::abs((overlap[0] - overlap[1])/2.),
                   std::abs((overlap[2] - overlap[3])/2.),
                   std::abs((overlap[4] - overlap[5])/2.));
  }

  void CRTHitRecoAlg::ReconstructPE(const TVector3 &pos, const CRTStripHit &hit0, const CRTStripHit &hit1, 
                                    double &pe0, double &pe1, 
                                    double &uncorrected_pe0, double &uncorrected_pe1,
                                    double &dist0, double &dist1)
  {
    dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0.channel);
    dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1.channel);

    pe0 = ReconstructPE(dist0, hit0, uncorrected_pe0);
    pe1 = ReconstructPE(dist1, hit1, uncorrected_pe1);

 
  }

  double CRTHitRecoAlg::ReconstructPE(const double &dist, const CRTStripHit &hit, double &uncorrected_pe)
  {
    const double stripPE = ADCtoPE(hit.adc1 + hit.adc2);

    uncorrected_pe = stripPE;
    
    return stripPE * std::pow(dist - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);
  }



  void CRTHitRecoAlg::CorrectADC(const TVector3 &pos, const CRTStripHit &hit0, 
                                 const CRTStripHit &hit1, std::array<uint16_t, 4> &corr_adc)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0.channel);
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1.channel);

    const double corr0 = std::pow(dist0 - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);
    const double corr1 = std::pow(dist1 - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);

    // only correct ADC when it is not saturated
    if (hit0.adc1+hit0.pedestal1 > fADCSaturation) corr_adc[0] = hit0.adc1;
    else corr_adc[0] = hit0.adc1 * corr0;
    if (hit0.adc2+hit0.pedestal2 > fADCSaturation) corr_adc[1] = hit0.adc2;
    else corr_adc[1] = hit0.adc2 * corr0;
    if (hit1.adc1+hit1.pedestal1 > fADCSaturation) corr_adc[2] = hit1.adc1;
    else corr_adc[2] = hit1.adc1 * corr1;
    if (hit1.adc2+hit1.pedestal2 > fADCSaturation) corr_adc[3] = hit1.adc2;
    else corr_adc[3] = hit1.adc2 * corr1;
    
  }

  double CRTHitRecoAlg::ADCtoPE(const uint16_t &adc)
  {
    return ((double)adc - fQPedestal) / fQSlope;
  }

  void CRTHitRecoAlg::CorrectTimings(const TVector3 &pos, const CRTStripHit &hit0, 
                                     const CRTStripHit &hit1, const double &pe0, const double &pe1,
                                     uint32_t &t0, uint32_t &t1, double &diff)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0.channel);
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1.channel);

    uint32_t corr0 = TimingCorrectionOffset(dist0, pe0);
    uint32_t corr1 = TimingCorrectionOffset(dist1, pe1);

    t0 = (hit0.t0 - corr0 + hit1.t0 - corr1) / 2.;
    t1 = (hit0.t1 - corr0 + hit1.t1 - corr1) / 2.;

    diff = (double)(hit0.t1 - corr0) - (double)(hit1.t1 - corr1);

    mf::LogInfo("CRTHitRecoAlg") <<"strip: t1: "<<t1-1.7e6<<", diff: "<<diff<<std::endl;
  }

  uint32_t CRTHitRecoAlg::TimingCorrectionOffset(const double &dist, const double &pe)
  {
    // The timing correction consists of two effects
    //    - the propagation delay (dependent on the distance the light needs
    //          to travel along the fibre)
    //    - the time walk effect (dependent on the size of the pulse)

    double t_TimeWalk = fTimeWalkNorm * std::exp(- fTimeWalkShift * pe);
    //fTimeWalkNorm * std::exp(-0.5 * std::pow((pe - fTimeWalkShift) / fTimeWalkSigma, 2)) + fTimeWalkOffset; 
    //if (t_TimeWalk<0) t_TimeWalk=0;
    return (uint32_t)(dist * fPropDelay + t_TimeWalk);
    
    //return (uint32_t)(t_TimeWalk);
    /*
    mf::LogInfo("CRTHitRecoAlg") << "distance to readout: "<< dist << ", fPropDelay: " << fPropDelay << "TProp: " << dist * fPropDelay << ", walkTime: " << t_TimeWalk << ", pe: " << pe << std::endl;
    if (t_TimeWalk>0.){
      return (uint32_t)(dist * fPropDelay + t_TimeWalk);
    }else {
      return (uint32_t)(dist * fPropDelay);
    }*/
  } 
}