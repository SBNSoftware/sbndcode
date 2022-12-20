#ifndef CRTBACKTRACKERALG_H_SEEN
#define CRTBACKTRACKERALG_H_SEEN

///////////////////////////////////////////////
// CRTClusterCharacterisationAlg.h
//
// Truth Matching Utilities for CRT analysis
// Henry Lay (h.lay@lancaster.ac.uk)
// November 2022
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// sbnobj
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

namespace sbnd::crt {
  
  class CRTClusterCharacterisationAlg {
  public:
    
    CRTClusterCharacterisationAlg();
    
    CRTClusterCharacterisationAlg(const fhicl::ParameterSet& pset);
    
    ~CRTClusterCharacterisationAlg();

    CRTSpacePoint CharacteriseSingleHitCluster(const art::Ptr<CRTCluster> &cluster, const art::Ptr<CRTStripHit> &stripHit);

    CRTSpacePoint CharacteriseDoubleHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits);

    double ADCToPE(const uint16_t channel, const uint16_t adc1, const uint16_t adc2);

    double ADCToPE(const uint16_t channel, const uint16_t adc);

    std::array<double, 6> FindOverlap(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1);

    void CentralPosition(const std::array<double, 6> overlap, TVector3 &pos, TVector3 &err);

    double ReconstructPE(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1, const TVector3 &pos);

    double ReconstructPE(const art::Ptr<CRTStripHit> &hit, const double dist);

  private:

    CRTGeoAlg fCRTGeoAlg;

    double fPEAttenuation;
  };
}

#endif
