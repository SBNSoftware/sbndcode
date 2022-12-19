#include "CRTClusterCharacterisationAlg.h"

namespace sbnd::crt {
  
  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(const fhicl::ParameterSet& pset)
    : fCRTGeoAlg(pset.get<fhicl::ParameterSet>("GeoAlg", fhicl::ParameterSet()))
  {
  }

  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(){}
  
  CRTClusterCharacterisationAlg::~CRTClusterCharacterisationAlg(){}

  CRTSpacePoint CRTClusterCharacterisationAlg::CharacteriseSingleHitCluster(const art::Ptr<CRTCluster> &cluster, const art::Ptr<CRTStripHit> &stripHit)
  {
    const std::array<double, 6> hitPos = fCRTGeoAlg.StripHit3DPos(stripHit->Channel(), stripHit->Pos(), stripHit->Error());

    const double pe = ADCToPE(stripHit->Channel(), stripHit->ADC1()) + 
      ADCToPE(stripHit->Channel() + 1, stripHit->ADC2());
    
    const double x  = (hitPos[0] + hitPos[1])/2.;
    const double ex = std::abs(hitPos[0] - hitPos[1])/2.;
    const double y  = (hitPos[2] + hitPos[3])/2.;
    const double ey = std::abs(hitPos[2] - hitPos[3])/2.;
    const double z  = (hitPos[4] + hitPos[5])/2.;
    const double ez = std::abs(hitPos[4] - hitPos[5])/2.;

    return CRTSpacePoint(x, ex, y, ey, z, ez, pe);
  }

  double CRTClusterCharacterisationAlg::ADCToPE(const uint16_t channel, const uint16_t adc)
  {
    return fCRTGeoAlg.GetSiPM(channel).gain * adc;
  }
}
