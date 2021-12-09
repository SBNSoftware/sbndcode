#ifndef SHOWERUTILS_H_SEEN
#define SHOWERUTILS_H_SEEN


///////////////////////////////////////////////
// Shower.h
//
// A few reco utilities like truth matching
// D Barker: Functions to are shower specific.
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "RecoUtils.h"
// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"


namespace ShowerUtils{

  std::map<int,std::vector<int> > GetShowerMothersCandidates(std::map<int,const simb::MCParticle*>& trueParticles);

  std::pair<int,double> TrueParticleIDFromTrueChain(const detinfo::DetectorClocksData& clockData, std::map<int,std::vector<int> >& ShowersMothers,const std::vector<art::Ptr<recob::Hit> >& hits, int planeid);

  std::map<geo::PlaneID,int> NumberofWiresHitByShower(const detinfo::DetectorClocksData& clockData, std::vector<int> &TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits);

  void CutShowerMothersByE(std::map<int,std::vector<int> >& ShowersMothers, std::map<int,const simb::MCParticle*>& trueParticles, float& EnergyCut);

  void CutShowerMothersByDensity(const detinfo::DetectorClocksData& clockData,
      std::map<int,std::vector<int> >& ShowersMothers,
      std::map<int,const simb::MCParticle*>& trueParticles,
      std::vector<art::Ptr<recob::Hit> >& hits,
      float& fDensityCut);

  void RemoveNoneContainedParticles(std::map<int,std::vector<int> >&  ShowersMothers, std::map<int,const simb::MCParticle*>& trueParticles, std::map<int,float>& MCTrack_Energy_map);


}

#endif
