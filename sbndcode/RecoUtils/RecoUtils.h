#ifndef FDSELECTIONUTILS_H_SEEN
#define FDSELECTIONUTILS_H_SEEN


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
//#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Shower.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larcore/Geometry/Geometry.h"


// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"


namespace RecoUtils{
  int TrueParticleID(const art::Ptr<recob::Hit>& hit);
  int TrueParticleIDFromTotalCharge(const std::vector<art::Ptr<recob::Hit> >& hits);
  int TrueParticleIDFromTotalHits(const std::vector<art::Ptr<recob::Hit> >& hits);
  bool IsInsideTPC(TVector3 position, double distance_buffer);
}

#endif
