////////////////////////////////////////////////////////////////////////
// Class:       DecayFilter
// Plugin Type: filter (Unknown Unknown)
// File:        DecayFilterG4_module.cc
//
// Generated at Thu Oct 21 03:39:48 2021 by Francisco Nicolas-Arnaldos using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

/*
 * Filter module for Geant4 decay events
 * Filter checks if a given mother particle decays into a specific set of daughter particles
 *
 * Event is kept if mother and daughters are found
 */

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Geometry
#include "larcore/Geometry/Geometry.h"
// LArSoft includes
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include <memory>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace hypana {
  class DecayFilter;
}


class hypana::DecayFilter : public art::EDFilter {
public:
  explicit DecayFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DecayFilter(DecayFilter const&) = delete;
  DecayFilter(DecayFilter&&) = delete;
  DecayFilter& operator=(DecayFilter const&) = delete;
  DecayFilter& operator=(DecayFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  void ResetVars();
  bool allDaughtersFound();

  std::string fMCModuleLabel;
  int fMotherPDG;
  std::vector<int> fDaughterPDG;

  bool fKeepEvent;
  int fParentTackID=-1;
  std::map<int, bool> fHasDaughter;

};


hypana::DecayFilter::DecayFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
  fMCModuleLabel ( p.get<std::string>("MCModuleLabel",  "largeant") ),
  fMotherPDG ( p.get<int>("MotherPDG",  0) ),
  fDaughterPDG ( p.get<std::vector<int>>("DaughterPDG",  {0}) )
{}

bool hypana::DecayFilter::filter(art::Event& e)
{

  art::Handle< std::vector<simb::MCParticle> > mclistLARG4;
  e.getByLabel(fMCModuleLabel,mclistLARG4);
  std::vector<simb::MCParticle> const& mcpartVec(*mclistLARG4);

  //Get mother track ID

  for (auto const& mother : mcpartVec) {
    if (mother.PdgCode() != fMotherPDG) continue;
    if (mother.EndProcess() != "Decay") continue;

    fParentTackID = mother.TrackId();

    //std::cout<<"Mother "<<mother.PdgCode()<<" ID="<<fParentTackID<<" EndProcess="<<mother.EndProcess()<<std::endl;

    ResetVars();

    for (auto const& p : mcpartVec) {
        if (p.Mother() != fParentTackID) continue;
        if (p.Process() != "Decay") continue;

        if (std::find(fDaughterPDG.begin(), fDaughterPDG.end(), p.PdgCode()) != fDaughterPDG.end()) {
            fHasDaughter[p.PdgCode()] = true;
            //std::cout << "Daughter " << p.PdgCode() << " ID=" << p.TrackId() << std::endl;

            if (allDaughtersFound()){
                //std::cout<<"KEEP EVENT: "<<true<<std::endl;
                return true;
            }
        }
    }
  }

  //std::cout<<"KEEP EVENT: "<<false<<std::endl;
  return false;

}

void hypana::DecayFilter::ResetVars(){
  fHasDaughter.clear();
  for(auto & pdg:fDaughterPDG){
    fHasDaughter[pdg]=false;
  }
}

bool hypana::DecayFilter::allDaughtersFound(){
    for (auto const& [pdg, found] : fHasDaughter) {
        if (!found) return false;
    }
    return true;
}

DEFINE_ART_MODULE(hypana::DecayFilter)