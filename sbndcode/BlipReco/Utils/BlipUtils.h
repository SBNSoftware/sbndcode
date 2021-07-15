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
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

// c++
#include <vector>
#include <map>

typedef std::vector<int> vi_t;
typedef std::set<int> si_t;
typedef std::map<int,float> mif_t;

namespace sbnd{

  namespace BlipUtils{

  //###################################################
  //  Data structures
  //###################################################

  struct TrueBlip {
    bool      isValid           = false;
    int       LeadingG4TrackID  = -9;
    float     LeadingEnergy     = -9;
    float     Energy            = 0;
    float     NumElectrons      = 0; // (post-drift)
    float     Length            = 0;
    TVector3  Position;
    TVector3  StartPoint;
    TVector3  EndPoint;
    vi_t      vG4TrackIDs;
    vi_t      vPDGs;
  };

  struct HitInfo {
    art::Ptr<recob::Hit> hit;
    si_t  g4ids;
    int   trkid;
    int   hitid;
    int   isreal; 
    int   g4id;
    float g4frac;
    float g4energy;
    float g4charge;
    float Charge;
    float Time;
  };
  
  struct HitClust {
    art::Ptr<recob::Hit> LeadHit;
    bool    isValid         = false;
    bool    isMerged        = false;
    bool    isMatched[3]    = {false, false, false};
    int     LeadHitWire     = -999;
    int     LeadHitG4TrackID= -9;
    int     TPC             = -9;
    int     Plane           = -9;
    float   LeadHitCharge   = -999;
    float   Charge          = -999;
    float   Time            = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    si_t    HitIDs;
    si_t    G4TrackIDs;
    si_t    Wires;
    mif_t   mapWireCharge;
  };

  struct Blip { 
    bool    isValid         = false;
    int     TPC             = -9;
    float   Charge[3]       = {-999, -999, -999};
    bool    WireIntersect[3]= {false, false, false};
    int     NCrossings;
    TVector3 Position;
    float   PositionRMS;
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
  bool      DoHitClustsMatch(HitClust const&, HitClust const&);
  bool      DoHitClustsMatch(HitClust const&,float,float);
  //Blip      MakeBlip(detinfo::DetectorPropertiesData const&, std::vector<HitClust> const&);
  Blip      MakeBlip(std::vector<HitClust> const&);

  //###################################################
  // General functions 
  //###################################################
  void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  bool    TrackIdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  bool    DoesHitHaveSimChannel( art::Ptr<recob::Hit> const&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);

  }
}

#endif
