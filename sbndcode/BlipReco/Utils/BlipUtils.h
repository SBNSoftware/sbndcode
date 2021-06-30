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

typedef std::vector<int> vi_t;
typedef std::set<int> si_t;

namespace BlipUtils{

  //###################################################
  //  Data structures
  //###################################################

  struct TrueBlip {
    bool      isValid;
    vi_t      vG4TrackIDs;
    vi_t      vPDGs;
    int       LeadingG4TrackID;
    float     LeadingEnergy;
    TVector3  Location;
    TVector3  StartPoint;
    TVector3  EndPoint;
    float     Energy;
    float     NumElectrons; // after drift
    float     Length;
    TrueBlip() {
      isValid = false;
      vG4TrackIDs.clear();
      vPDGs.clear();
      LeadingG4TrackID = -9;
      LeadingEnergy = -9;
      Energy  = 0;
      Length  = 0;
      NumElectrons = 0;
    }
  };

  struct HitInfo {
    art::Ptr<recob::Hit> hit;
    si_t  g4ids;
    int   trkid;
    int   hitid;
    float Charge;
    float Time;
    int   isreal; 
    int   g4id;
    int   g4frac;
    int   g4energy;
    int   g4charge;
  };
  
  struct HitClust {
    bool    isValid;
    si_t    HitIDs;
    si_t    G4TrackIDs;
    si_t    Wires;
    std::map<int,float> mapWireCharge;
    int     LeadWire;
    int     TPC;
    int     Plane;
    float   Charge;
    float   LeadHitCharge;
    float   Time;
    float   StartTime;
    float   EndTime;
  };

  struct Blip { 
    bool  isValid;
    si_t  vG4TrackIDs;
    si_t  vHitIDs;
    float LeadingG4TrackID;
    float LeadingG4Energy;
    float G4Energy;
    int   TPC;
    bool  MatchedPlanes[3];
    bool  Is3D;
    float X,Y,Z;
    float Charge;
    float Energy;
  };

  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  float     PartEnergyDep(int, float&);
  float     PartEnergyDep(int);
  TrueBlip  MakeTrueBlip(int);
  void      GrowTrueBlip(simb::MCParticle const&, TrueBlip&);
  void      MergeBlips(std::vector<TrueBlip>&, float);
  HitClust  MakeHitClust(art::Ptr<recob::Hit> const&, HitInfo const&);
  void      GrowHitClust(art::Ptr<recob::Hit> const&, HitInfo const&, HitClust&);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);

  //###################################################
  // General functions 
  //###################################################
  void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  bool    TrackIdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  bool    DoesHitHaveSimChannel( std::vector<const sim::SimChannel*>, art::Ptr<recob::Hit> const&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);

}

#endif
