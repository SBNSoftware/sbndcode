///////////////////////////////////////////////////////
// BlipUtils.h
//
// Helper functions and algs. Based on RecoUtils and
// functions in AnalysisTree.
//
//////////////////////////////////////////////////////
#ifndef BLIPUTIL_H_SEEN
#define BLIPUTIL_H_SEEN

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// c++
#include <vector>
#include <map>

namespace detinfo { class DetectorClocksData; }

// Helper templates for initializing arrays
namespace{  
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, ITER to, TYPE value) 
    { std::fill(from, to, value); }
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }
}

namespace BlipUtils{

  void   HitsPurity(//detinfo::DetectorClocksData const& clockData,
                    std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, float& purity, double& maxe);
  void   HitTruth(//detinfo::DetectorClocksData const& clockData, 
                    art::Ptr<recob::Hit> const& hit, int& truthid, float& frac, float& energy, float& numElectrons);
  bool   HitTruthId( //detinfo::DetectorClocksData const& clockData, 
                    art::Ptr<recob::Hit> const& hit, int& mcid);
  bool   TrackIdToMCTruth( int const trkID, art::Ptr<simb::MCTruth>& mctruth );
  bool   DoesHitHaveSimChannel( std::vector<const sim::SimChannel*> chans, art::Ptr<recob::Hit> const& hit);

}




#endif
