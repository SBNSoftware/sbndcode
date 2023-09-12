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
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"

// c++
#include <vector>
#include <map>

#include "sbndcode/BlipRecoSBND/Utils/DataTypes.h"
#include "TH1D.h"


typedef std::vector<art::Ptr<sim::SimEnergyDeposit>> SEDVec_t;

geo::View_t kViews[3]={geo::kU, geo::kV, geo::kW};

namespace BlipUtils {
 
  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  //void      InitializeDetProps();
  void      FillParticleInfo(simb::MCParticle const&, blip::ParticleInfo&, SEDVec_t&, int plane=2);
  //void      CalcPartDrift(blip::ParticleInfo&, int);
  //void      CalcTotalDep(float&,int&,float&, SEDVec_t&);
  void      MakeTrueBlips(std::vector<blip::ParticleInfo>&, std::vector<blip::TrueBlip>&);
  void      GrowTrueBlip(blip::ParticleInfo&, blip::TrueBlip&);
  void      MergeTrueBlips(std::vector<blip::TrueBlip>&, float);
  void      GrowHitClust(blip::HitInfo const&, blip::HitClust&);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);
  bool      DoHitClustsOverlap(blip::HitClust const&, blip::HitClust const&);
  bool      DoHitClustsOverlap(blip::HitClust const&,float,float);
  float     CalcHitClustsOverlap(blip::HitClust const&, blip::HitClust const&);
  float     CalcOverlap(float const&, float const&, float const&, float const&);
  bool      DoChannelsIntersect(int,int);
  bool      DoHitClustsMatch(blip::HitClust const&, blip::HitClust const&,float);
  blip::HitClust  MakeHitClust(std::vector<blip::HitInfo> const&);
  blip::Blip      MakeBlip(std::vector<blip::HitClust> const&);
  

  //###################################################
  // General functions 
  //###################################################
  //void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  //si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  //bool    G4IdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool);
  double  DistToBoundary(const recob::Track::Point_t&);
  double  DistToLine(TVector3&, TVector3&, TVector3&);
  double  DistToLine2D(TVector2&, TVector2&, TVector2&);
  bool    IsInActiveVolume(geo::Point_t const&);
  void    NormalizeHist(TH1D*);
  float   FindMedian(std::vector<float>&);
  float   FindMean(std::vector<float>&);
  
}

#endif
