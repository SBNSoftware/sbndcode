#include "CRTHitRecoAlg.h"

namespace sbnd{

  CRTHitRecoAlg::CRTHitRecoAlg(const fhicl::ParameterSet &p) 
    : fCRTGeoAlg(p.get<fhicl::ParameterSet>("GeoAlg"))
    , fADCThreshold(p.get<uint16_t>("ADCThreshold"))
    , fQPedestal(p.get<double>("QPedestal"))
    , fQSlope(p.get<double>("QSlope"))
    , fPEAttenuation(p.get<double>("PEAttenuation"))
    , fPropDelay(p.get<double>("PropDelay"))
    , fTimeWalkNorm(p.get<double>("TimeWalkNorm"))
    , fTimeWalkShift(p.get<double>("TimeWalkShift"))
    , fTimeWalkSigma(p.get<double>("TimeWalkSigma"))
    , fTimeWalkOffset(p.get<double>("TimeWalkOffset"))
    , fHitCoincidenceRequirement(p.get<double>("HitCoincidenceRequirement"))
  {}

  CRTHitRecoAlg::CRTHitRecoAlg() {}

  CRTHitRecoAlg::~CRTHitRecoAlg() {}

  std::map<std::string, std::vector<std::vector<CRTStripHit>>> 
    CRTHitRecoAlg::ProduceStripHits(const std::vector<art::Ptr<sbnd::crt::FEBData>> &dataVec)
  {
    std::map<std::string, std::vector<std::vector<CRTStripHit>>> stripHits;

    // Iterate through each FEBData
    for(unsigned feb_i = 0; feb_i < dataVec.size(); ++feb_i)
      {
        art::Ptr<sbnd::crt::FEBData> data = dataVec.at(feb_i);
        uint32_t mac5 = data->Mac5();

        CRTModuleGeo module    = fCRTGeoAlg.GetModule(mac5 * 32);
        std::string taggerName = module.taggerName;

        // Resize the vector of vectors to 2, one vector of hits for each
        // of the orientations
        if(stripHits.find(taggerName) == stripHits.end())
          stripHits[taggerName].resize(2);

        // Correct for FEB readout cable length
        uint32_t t0 = data->Ts0() + module.cableDelayCorrection;
        uint32_t t1 = data->Ts1() + module.cableDelayCorrection;

        // Iterate via strip (2 SiPMs per strip)
        const auto &sipm_adcs = data->ADC();
        for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
          {
            // Add an offset for the SiPM channel number
            uint16_t channel = mac5 * 32 + adc_i;

            CRTStripGeo strip = fCRTGeoAlg.GetStrip(channel);
            CRTSiPMGeo sipm1  = fCRTGeoAlg.GetSiPM(channel);
            CRTSiPMGeo sipm2  = fCRTGeoAlg.GetSiPM(channel+1);

            // Subtract channel pedestals
            uint16_t adc1 = sipm1.pedestal < sipm_adcs[adc_i] ? sipm_adcs[adc_i]   - sipm1.pedestal : 0;
            uint16_t adc2 = sipm2.pedestal < sipm_adcs[adc_i+1] ? sipm_adcs[adc_i+1]   - sipm2.pedestal : 0;

            // Keep hit if both SiPMs above threshold
            if(adc1 > fADCThreshold && adc2 > fADCThreshold)
              {
                // Access width of strip from the geometry algorithm
                // === TO-DO === //
                // AMEND THE CODE AND IMPROVE CALCULATION
                double width = strip.width;
                double x     = width/2.;
                double ex    = width/2.;

                // Create hit
                stripHits[module.taggerName][module.orientation].emplace_back(channel, t0, t1, x, ex, adc1, adc2, feb_i);
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

        std::set<unsigned> used;

        for(unsigned i = 0; i < hitsOrien0.size(); ++i)
          {
            const CRTStripHit hit0   = hitsOrien0[i];
            const CRTStripGeo strip0 = fCRTGeoAlg.GetStrip(hit0.channel);

            for(unsigned j = 0; j < hitsOrien1.size(); ++j)
              {
                if(used.find(j) != used.end())
                  continue;

                const CRTStripHit hit1     = hitsOrien1[j];
                const CRTStripGeo strip1 = fCRTGeoAlg.GetStrip(hit1.channel);

                // Check whether the two strips responsible for these hits overlap.
                if(!fCRTGeoAlg.CheckOverlap(strip0, strip1))
                   continue;

                // Find overlap region between two strip hits
                std::vector<double> overlap = FindOverlap(hit0, hit1, strip0, strip1);

                // Using overlap region find centre and 'error'
                TVector3 pos, err;
                CentralPosition(overlap, pos, err);
                
                // Reconstruct the number of photoelectrons from the ADC values
                double pe0, pe1;
                ReconstructPE(pos, hit0, hit1, pe0, pe1);

                // Correct timings to find how coincident the hits were
                uint32_t t0, t1, diff;
                CorrectTimings(pos, hit0, hit1, pe0, pe1, t0, t1, diff);

                if(diff > fHitCoincidenceRequirement)
                  continue;

                sbn::crt::CRTHit crtHit({(uint8_t)hit0.febdataindex, (uint8_t)hit1.febdataindex},
                                        pe0+pe1,
                                        t0,
                                        t1,
                                        pos,
                                        err,
                                        tagger);

                mf::LogInfo("CRTHitRecoAlg") << "\nCreating CRTHit"
                                             << "from FEBs: " << (unsigned) crtHit.feb_id[0] 
                                             << " & " << (unsigned) crtHit.feb_id[1] << '\n'
                                             << "at position: " << pos.X() << ", " 
                                             << pos.Y() << ", " << pos.Z() << "cm\n"
                                             << "and t1: " << t1 << '\n'
                                             << "with PE: " << pe0+pe1 << '\n'
                                             << "on tagger: " << tagger << std::endl;

                used.insert(j);
                crtHits.emplace_back(crtHit);
                break;
              }
          }
      }
    return crtHits;
  }

  std::vector<double> CRTHitRecoAlg::FindOverlap(const CRTStripHit &hit0, const CRTStripHit &hit1,
                                                 const CRTStripGeo &strip0, const CRTStripGeo &strip1)
  {
    const std::vector<double> hit0pos = fCRTGeoAlg.StripHit3DPos(strip0.name, hit0.x, hit0.ex);
    const std::vector<double> hit1pos = fCRTGeoAlg.StripHit3DPos(strip1.name, hit1.x, hit1.ex);

    std::vector<double> overlap({std::max(hit0pos[0], hit1pos[1]),
                                 std::min(hit0pos[1], hit1pos[1]),
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

  void CRTHitRecoAlg::ReconstructPE(const TVector3 &pos, const CRTStripHit &hit0, 
                                    const CRTStripHit &hit1, double &pe0, double &pe1)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0.channel);
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1.channel);

    pe0 = ReconstructPE(dist0, hit0);
    pe1 = ReconstructPE(dist1, hit1);
  }

  double CRTHitRecoAlg::ReconstructPE(const double &dist, const CRTStripHit &hit)
  {
    const double stripPE = ADCtoPE(hit.adc1 + hit.adc2);
    return stripPE * std::pow(dist - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);
  }

  double CRTHitRecoAlg::ADCtoPE(const uint16_t &adc)
  {
    return ((double)adc - fQPedestal) / fQSlope;
  }

  void CRTHitRecoAlg::CorrectTimings(const TVector3 &pos, const CRTStripHit &hit0, 
                                     const CRTStripHit &hit1, const double &pe0, const double &pe1,
                                     uint32_t &t0, uint32_t &t1, uint32_t &diff)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0.channel);
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1.channel);

    const double corr0 = TimingCorrectionOffset(dist0, pe0);
    const double corr1 = TimingCorrectionOffset(dist1, pe1);

    t0 = (hit0.t0 - corr0 + hit1.t0 - corr1) / 2.;
    t1 = (hit0.t1 - corr0 + hit1.t1 - corr1) / 2.;
    
    diff = std::abs((hit0.t1 - corr0) - (hit1.t1 - corr1));
  }

  double CRTHitRecoAlg::TimingCorrectionOffset(const double &dist, const double &pe)
  {
    // The timing correction consists of two effects
    //    - the propagation delay (dependent on the distance the light needs
    //          to travel along the fibre)
    //    - the time walk effect (dependent on the size of the pulse)
    return dist * fPropDelay + fTimeWalkNorm * std::exp(-0.5 * std::pow((pe - fTimeWalkShift) / fTimeWalkSigma, 2)) + fTimeWalkOffset;
  } 
}
