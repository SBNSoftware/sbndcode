#include "CRTClusterCharacterisationAlg.h"

namespace sbnd::crt {
  
  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(const fhicl::ParameterSet& pset)
    : fCRTGeoAlg(pset.get<fhicl::ParameterSet>("GeoAlg", fhicl::ParameterSet()))
    , fUseT1(pset.get<bool>("UseT1", true))
    , fOverlapBuffer(pset.get<double>("OverlapBuffer"))
    , fPEAttenuation(pset.get<double>("PEAttenuation"))
    , fPropDelay(pset.get<double>("PropDelay"))
    , fTimeWalkNorm(pset.get<double>("TimeWalkNorm"))
    , fTimeWalkShift(pset.get<double>("TimeWalkShift"))
    , fTimeWalkSigma(pset.get<double>("TimeWalkSigma"))
    , fTimeWalkOffset(pset.get<double>("TimeWalkOffset"))
  {
  }

  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(){}
  
  CRTClusterCharacterisationAlg::~CRTClusterCharacterisationAlg(){}

  CRTSpacePoint CRTClusterCharacterisationAlg::CharacteriseSingleHitCluster(const art::Ptr<CRTCluster> &cluster, const art::Ptr<CRTStripHit> &stripHit)
  {
    const std::array<double, 6> hitPos = fCRTGeoAlg.StripHit3DPos(stripHit->Channel(), stripHit->Pos(), stripHit->Error());

    const double pe = ADCToPE(stripHit->Channel(), stripHit->ADC1(), stripHit->ADC2());
    
    TVector3 pos, err;
    CentralPosition(hitPos, pos, err);

    return CRTSpacePoint(pos, err, pe, stripHit->Ts1(), false);
  }

  bool CRTClusterCharacterisationAlg::CharacteriseDoubleHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits, CRTSpacePoint &spacepoint)
  {
    const art::Ptr<CRTStripHit> &hit0 = stripHits[0];
    const art::Ptr<CRTStripHit> &hit1 = stripHits[1];

    return TwoHitSpacePoint(hit0, hit1, spacepoint);
  }

  bool CRTClusterCharacterisationAlg::TwoHitSpacePoint(const art::Ptr<CRTStripHit> hit0, const art::Ptr<CRTStripHit> hit1, CRTSpacePoint &spacepoint)
  {
    const bool threeD = fCRTGeoAlg.ChannelToOrientation(hit0->Channel()) != fCRTGeoAlg.ChannelToOrientation(hit1->Channel());

    if(threeD)
      {
        if(fCRTGeoAlg.CheckOverlap(hit0->Channel(), hit1->Channel(), fOverlapBuffer))
          {
            std::array<double, 6> overlap = FindOverlap(hit0, hit1);
            
            TVector3 pos, err;
            CentralPosition(overlap, pos, err);

            const double pe   = ReconstructPE(hit0, hit1, pos);
            const double time = CorrectTime(hit0, hit1, pos);

            spacepoint = CRTSpacePoint(pos, err, pe, time, true);
            return true;
          }
        return false;
      }
    else
      {
        if(fCRTGeoAlg.AdjacentStrips(hit0->Channel(), hit1->Channel(), fOverlapBuffer))
          {
            std::array<double, 6> overlap = FindAdjacentPosition(hit0, hit1);

            TVector3 pos, err;
            CentralPosition(overlap, pos, err);

            const double pe   = ADCToPE(hit0->Channel(), hit0->ADC1(), hit0->ADC2()) + ADCToPE(hit1->Channel(), hit1->ADC1(), hit1->ADC2());
            const double time = (hit0->Ts1() + hit1->Ts1()) / 2.;

            spacepoint = CRTSpacePoint(pos, err, pe, time, false);
            return true;
          }
        return false;
      }
  }

  bool CRTClusterCharacterisationAlg::CharacteriseMultiHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits, CRTSpacePoint &spacepoint)
  {
    std::vector<CRTSpacePoint> spacepoints, complete_spacepoints;

    for(unsigned i = 0; i < stripHits.size(); ++i)
      {
        const art::Ptr<CRTStripHit> hit0 = stripHits[i];

        for(unsigned ii = i + 1; ii < stripHits.size(); ++ii)
          {
            const art::Ptr<CRTStripHit> hit1 = stripHits[ii];

            CRTSpacePoint sp;
            if(TwoHitSpacePoint(hit0, hit1, sp))
              {
                spacepoints.push_back(sp);
        
                if(sp.Complete()) 
                  complete_spacepoints.push_back(sp);
              }
          }
      }

    if(complete_spacepoints.size() == 0)
      return false;
    
    std::array<double, 6> aggregate_position({std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 
          -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()});

    double pe = 0., time = 0.;

    for(auto const &sp : complete_spacepoints)
      {
        AggregatePositions(sp.Pos(), sp.Err(), aggregate_position);
        pe   += sp.PE();
        time += sp.Time();
      }

    TVector3 pos, err;
    CentralPosition(aggregate_position, pos, err);
    pe   /= complete_spacepoints.size();
    time /= complete_spacepoints.size();
    
    spacepoint = CRTSpacePoint(pos, err, pe, time, true);
    return true;
  }

  double CRTClusterCharacterisationAlg::ADCToPE(const uint16_t channel, const uint16_t adc1, const uint16_t adc2)
  {
    return ADCToPE(channel, adc1) + ADCToPE(channel+1, adc2);
  }

  double CRTClusterCharacterisationAlg::ADCToPE(const uint16_t channel, const uint16_t adc)
  {
    return fCRTGeoAlg.GetSiPM(channel).gain * adc;
  }

  std::array<double, 6> CRTClusterCharacterisationAlg::FindOverlap(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1)
  {
    const std::array<double, 6> hit0pos = fCRTGeoAlg.StripHit3DPos(hit0->Channel(), hit0->Pos(), hit0->Error());
    const std::array<double, 6> hit1pos = fCRTGeoAlg.StripHit3DPos(hit1->Channel(), hit1->Pos(), hit1->Error());

    std::array<double, 6> overlap({std::max(hit0pos[0], hit1pos[0]),
                                   std::min(hit0pos[1], hit1pos[1]),
                                   std::max(hit0pos[2], hit1pos[2]),
                                   std::min(hit0pos[3], hit1pos[3]),
                                   std::max(hit0pos[4], hit1pos[4]),
                                   std::min(hit0pos[5], hit1pos[5])});

    return overlap;
  }

  std::array<double, 6> CRTClusterCharacterisationAlg::FindAdjacentPosition(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1)
  {
    const std::array<double, 6> hit0pos = fCRTGeoAlg.StripHit3DPos(hit0->Channel(), hit0->Pos(), hit0->Error());
    const std::array<double, 6> hit1pos = fCRTGeoAlg.StripHit3DPos(hit1->Channel(), hit1->Pos(), hit1->Error());

    std::array<double, 6> overlap({std::min(hit0pos[0], hit1pos[0]),
                                   std::max(hit0pos[1], hit1pos[1]),
                                   std::min(hit0pos[2], hit1pos[2]),
                                   std::max(hit0pos[3], hit1pos[3]),
                                   std::min(hit0pos[4], hit1pos[4]),
                                   std::max(hit0pos[5], hit1pos[5])});

    return overlap;
  }

  void CRTClusterCharacterisationAlg::CentralPosition(const std::array<double, 6> overlap, 
                                                      TVector3 &pos, TVector3 &err)
  {
    pos = TVector3((overlap[0] + overlap[1])/2.,
                   (overlap[2] + overlap[3])/2.,
                   (overlap[4] + overlap[5])/2.);
    
    err = TVector3(std::abs((overlap[0] - overlap[1])/2.),
                   std::abs((overlap[2] - overlap[3])/2.),
                   std::abs((overlap[4] - overlap[5])/2.));
  }

  double CRTClusterCharacterisationAlg::ReconstructPE(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1, const TVector3 &pos)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0->Channel());
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1->Channel());

    return ReconstructPE(hit0, dist0) + ReconstructPE(hit1, dist1);
  }

  double CRTClusterCharacterisationAlg::ReconstructPE(const art::Ptr<CRTStripHit> &hit, const double dist)
  {
    const double pe         = ADCToPE(hit->Channel(), hit->ADC1(), hit->ADC2());
    const double correction = std::pow(dist - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);

    return pe * correction;
  }

  double CRTClusterCharacterisationAlg::CorrectTime(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1, const TVector3 &pos)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0->Channel());
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1->Channel());

    const double pe0 = ReconstructPE(hit0, dist0);
    const double pe1 = ReconstructPE(hit1, dist1);

    const double corr0 = TimingCorrectionOffset(dist0, pe0);
    const double corr1 = TimingCorrectionOffset(dist1, pe1);

    if(fUseT1)
      return (hit0->Ts1() - corr0 + hit1->Ts1() - corr1) / 2.;
    else
      return (hit0->Ts0() - corr0 + hit1->Ts0() - corr1) / 2.;
  }

  double CRTClusterCharacterisationAlg::TimingCorrectionOffset(const double &dist, const double &pe)
  {
    return dist * fPropDelay + fTimeWalkNorm * std::exp(-0.5 * std::pow((pe - fTimeWalkShift) / fTimeWalkSigma, 2)) + fTimeWalkOffset;
  }

  void CRTClusterCharacterisationAlg::AggregatePositions(const TVector3 &pos, const TVector3 &err, std::array<double, 6> &agg)
  {
    agg[0] = std::min(pos.X() - err.X(), agg[0]);
    agg[1] = std::max(pos.X() + err.X(), agg[1]);
    agg[2] = std::min(pos.Y() - err.Y(), agg[2]);
    agg[3] = std::max(pos.Y() + err.Y(), agg[3]);
    agg[4] = std::min(pos.Z() - err.Z(), agg[4]);
    agg[5] = std::max(pos.Z() + err.Z(), agg[5]);
  }
}
