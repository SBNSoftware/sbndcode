#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN

#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbnobj/SBND/CRT/FEBData.hh"

namespace sbnd{

  struct CRTStripHit {
    CRTStripHit(uint32_t _channel, uint32_t _t0, uint32_t _t1, double _x,
		double _ex, uint16_t _adc1, uint16_t _adc2) :
      channel(_channel)
      , t0   (_t0)
      , t1   (_t1)
      , x    (_x)
      , ex   (_ex)
      , adc1 (_adc1)
      , adc2 (_adc2)
    {};

    uint32_t channel;  // Channel ID for 1st SiPM
    uint32_t t0;       // T0 counter [ns]
    uint32_t t1;       // T1 counter [ns]
    double   x;        // Lateral position within strip [cm]
    double   ex;       // Error on lateral position [cm]
    uint16_t adc1;     // ADC 1st SiPM
    uint16_t adc2;     // ADC 2nd SiPM
  };


  class CRTHitRecoAlg {
  public:

    CRTHitRecoAlg(const fhicl::ParameterSet &p);

    CRTHitRecoAlg();

    ~CRTHitRecoAlg();

    std::vector<CRTStripHit> ProduceStripHits(std::vector<art::Ptr<sbnd::crt::FEBData>> &dataVec);
    
  private:
    
    CRTGeoAlg fCRTGeoAlg;
    uint16_t fADCThreshold;
  };

}

#endif
