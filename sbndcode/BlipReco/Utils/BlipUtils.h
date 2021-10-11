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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
//namespace detinfo {
//  class DetectorClocksData;
//  class DetectorPropertiesData;
//}

// c++
#include <vector>
#include <map>

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

typedef std::vector<int> vi_t;
typedef std::set<int> si_t;
typedef std::map<int,float> mif_t;
typedef std::map<float,float> mff_t;

namespace sbnd{

  namespace BlipUtils{

  //###################################################
  //  Data structures
  //###################################################

  struct TrueBlip {
    bool      isValid       = false;
    int       TPC           = -9;
    int       LeadG4ID      = -9;
    int       LeadG4PDG     = -9;
    float     LeadEnergy    = -9;
    float     Energy        = 0;
    float     NumElectrons  = 0; // (post-drift)
    float     Length        = 0;
    TVector3  Position;
    TVector3  StartPoint;
    TVector3  EndPoint;
    vi_t      G4IDs;
    vi_t      PDGs;
    int       ID;
  };

  struct HitInfo {
    art::Ptr<recob::Hit> hit;
    si_t  g4ids;
    int   plane   = -9;
    int   tpc     = -9;
    int   wire    = -9;
    int   trkid   = -9;
    int   shwrid  = -9;
    int   hitid   = -9;
    bool  isreal  = false;
    bool  ismatch = false;
    int   g4id    = -9;
    float g4frac  = -99;
    float g4energy= -999;
    float g4charge= -999;
    float qcoll   = -999;
    float driftTicks = -999999;
  };
  
  struct HitClust {
    art::Ptr<recob::Hit> LeadHit;
    bool    isValid         = false;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     LeadHitWire     = -999;
    int     LeadHitG4ID     = -9;
    float   LeadHitCharge   = -999;
    float   LeadHitTime     = -999;
    int     TPC             = -9;
    int     Plane           = -9;
    float   Charge          = -999;
    float   Time            = -999;
    float   WeightedTime    = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     ID              = -9;
    int     BlipID          = -9;
    si_t    HitIDs;
    si_t    G4IDs;
    si_t    Wires;
    mif_t   mapWireCharge;
    mff_t   mapTimeCharge;
  };

  struct Blip {
    bool      isValid         = false;
    int       TPC             = -9;
    int       NPlanes         = -9;
    int       Planes[3]       = {false, false, false};
    float     Charge[3]       = {-999, -999, -999};
    float     Energy          = -999;
    float     EnergyESTAR     = -999;
    float     DriftTime       = -999;
    float     MaxIntersectDiff= -9;
    int       ID              = -9;
    TVector3  Position;
    si_t      ClustIDs;
    si_t      HitIDs;
  };

  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  void      CalcTotalDep(float&,float&);
  void      CalcPartDep(int, float&,float&);
  TrueBlip  MakeTrueBlip(int);
  void      GrowTrueBlip(simb::MCParticle const&, TrueBlip&);
  void      MergeTrueBlips(std::vector<TrueBlip>&, float);
  HitClust  MakeHitClust(HitInfo const&);
  HitClust  MergeHitClusts(HitClust&,HitClust&);
  void      GrowHitClust(HitInfo const&, HitClust&);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);
  bool      DoHitClustsOverlap(HitClust const&, HitClust const&);
  bool      DoHitClustsOverlap(HitClust const&,float,float);
  bool      DoHitClustsMatch(HitClust const&, HitClust const&,float);
  //Blip      MakeBlip(detinfo::DetectorPropertiesData const&, std::vector<HitClust> const&);
  Blip      MakeBlip(std::vector<HitClust> const&);
  float     ModBoxRecomb(float,float);


  //###################################################
  // General functions 
  //###################################################
  void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  bool    G4IdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  bool    DoesHitHaveSimChannel( art::Ptr<recob::Hit> const&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);
  double  DistToBoundary(const recob::Track::Point_t&);
  void    GetGeoBoundaries(double&,double&,double&,double&,double&,double&);

  }
}

#endif
