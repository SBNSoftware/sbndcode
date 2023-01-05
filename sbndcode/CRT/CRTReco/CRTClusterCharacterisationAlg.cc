
#include "CRTClusterCharacterisationAlg.h"

namespace sbnd::crt {
  
  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(const fhicl::ParameterSet& pset)
    : fCRTGeoAlg(pset.get<fhicl::ParameterSet>("GeoAlg", fhicl::ParameterSet()))
    , fUseT1(pset.get<bool>("UseT1", true))
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

    return CRTSpacePoint(pos, err, pe, stripHit->Ts1());
  }

  CRTSpacePoint CRTClusterCharacterisationAlg::CharacteriseDoubleHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits)
  {
    if(cluster->ThreeD())
      {
        const art::Ptr<CRTStripHit> &hit0 = stripHits[0];
        const art::Ptr<CRTStripHit> &hit1 = stripHits[1];

        if(fCRTGeoAlg.CheckOverlap(hit0->Channel(), hit1->Channel()))
          {
            std::array<double, 6> overlap   = FindOverlap(hit0, hit1);
            
            TVector3 pos, err;
            CentralPosition(overlap, pos, err);

            const double pe   = ReconstructPE(hit0, hit1, pos);
            const double time = CorrectTime(hit0, hit1, pos);

	    std::cout << "Both (" << cluster.key() << " " << cluster->Ts1() << ")" << std::endl;
            return CRTSpacePoint(pos, err, pe, time);
          }
	std::cout << "ThreeD but no overlap (" << cluster.key() << " " << cluster->Ts1() << ")" << std::endl;
	return CRTSpacePoint();
      }
    std::cout << "Not threeD (" << cluster.key() << " " << cluster->Ts1() << ")" << std::endl;
    return CRTSpacePoint();
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
}
