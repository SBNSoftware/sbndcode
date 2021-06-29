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

namespace BlipUtils{

  //###################################################
  //  Data structure for true energy deposition
  //###################################################
  struct TrueBlip {
    bool    isValid;
    std::vector<int> vG4TrackIDs;
    std::vector<int> vPDGs;
    int       LeadingG4TrackID;
    float     LeadingEnergy;
    TVector3  Location;
    TVector3  StartPoint;
    TVector3  EndPoint;
    float   Energy;
    float   NumElectrons; // after drift
    float   Length;
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
  
  //###################################################
  //  Data structure for hit clusters
  //###################################################
  struct HitClust {
    bool isValid;
    std::vector<int> vHitIDs;
    std::vector<int> vWires;
    int   LeadingG4TrackID;
    int   LeadingG4Energy;
    int   LeadingHitID;
    int   Wire;
    int   TPC;
    int   Plane;
    float Area;
    float Charge;
    float StartTime;
    float EndTime;
    float Time;
    float Width;
  };

  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  float     PartEnergyDep(int, float&);
  float     PartEnergyDep(int);
  TrueBlip  MakeTrueBlip(int);
  void      GrowTrueBlip(simb::MCParticle&, TrueBlip&);
  void      MergeBlips(std::vector<TrueBlip>&, float);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);

  //###################################################
  // General functions 
  //###################################################
  void    HitsPurity(std::vector< art::Ptr<recob::Hit> > const&, int&, float&, double&);
  void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  bool    HitTruthId(art::Ptr<recob::Hit> const&, int&);
  bool    TrackIdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  bool    DoesHitHaveSimChannel( std::vector<const sim::SimChannel*>, art::Ptr<recob::Hit> const&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);

}

#endif
