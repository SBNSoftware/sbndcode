////////////////////////////////////////////////////////////////////////
// Class:       MCCaloAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        MCCaloAna_module.cc
//
// Generated at Tue Apr 29 12:33:40 2025 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art_root_io/TFileService.h"
#include "TFile.h"
#include "TTree.h"

#include <algorithm>

namespace sbnd {
  class MCCaloAna;
}


class sbnd::MCCaloAna : public art::EDAnalyzer {
public:
  explicit MCCaloAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCCaloAna(MCCaloAna const&) = delete;
  MCCaloAna(MCCaloAna&&) = delete;
  MCCaloAna& operator=(MCCaloAna const&) = delete;
  MCCaloAna& operator=(MCCaloAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  TTree* _tree;

  int _run;
  int _subrun;
  int _event;

  std::vector<double> _nu_E;
  std::vector<int> _nu_CCNC;
  std::vector<double> _true_gamma;
  std::vector<double> _true_charge;
  std::vector<double> _true_energy;

  // Declare member data here.

};


sbnd::MCCaloAna::MCCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("mccalo_tree","");
  _tree->Branch("run",             &_run,        "run/I");
  _tree->Branch("subrun",          &_subrun,     "subrun/I");
  _tree->Branch("event",           &_event,      "event/I");
  _tree->Branch("nu_E",           "std::vector<double>",  &_nu_E);
  _tree->Branch("nu_CCNC",        "std::vector<int>",     &_nu_CCNC);
  _tree->Branch("true_gamma",      "std::vector<double>",  &_true_gamma);
  _tree->Branch("true_charge",     "std::vector<double>",  &_true_charge);
  _tree->Branch("true_energy",     "std::vector<double>",  &_true_energy);
}

void sbnd::MCCaloAna::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _nu_E.clear();
  _nu_CCNC.clear();
  _true_gamma.clear();
  _true_charge.clear();
  _true_energy.clear();

  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  art::Handle<std::vector<simb::MCTruth>> mctruth_handle;
  e.getByLabel("generator", mctruth_handle);

  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  if (mctruth_handle.isValid()) 
    art::fill_ptr_vector(mctruths, mctruth_handle);

  ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
  e.getByLabel("ionandscint", energyDeps_h);
  std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
  
  if (!energyDeps_h.isValid() || energyDeps_h->empty())
    std::cout << "Don't have good SimEnergyDeposits!" << std::endl;
  else 
    art::fill_ptr_vector(energyDeps, energyDeps_h);
    
  _nu_E.resize(mctruths.size(), 0);
  _nu_CCNC.resize(mctruths.size(), 0);
  _true_gamma.resize(mctruths.size(), 0);
  _true_charge.resize(mctruths.size(), 0);
  _true_energy.resize(mctruths.size(), 0);

  for (size_t n_dep=0; n_dep < energyDeps.size(); n_dep++){
    auto energyDep = energyDeps[n_dep];
    const auto trackID = energyDep->TrackID();

    art::Ptr<simb::MCTruth> mctruth = piserv->TrackIdToMCTruth_P(trackID);
    if (mctruth->Origin() != simb::kBeamNeutrino) continue;

    auto it = std::find(mctruths.begin(), mctruths.end(), mctruth);
    if (it == mctruths.end()) {
      std::cout << "No matching MCTruth found for trackID: " << trackID << std::endl;
      continue;
    }

    // get the index of the mctruth in the set 
    auto mctruth_index = std::distance(mctruths.begin(), it);

    _true_gamma[mctruth_index]  += energyDep->NumPhotons();
    _true_charge[mctruth_index] += energyDep->NumElectrons();
    _true_energy[mctruth_index] += energyDep->Energy();
  }  
  _tree->Fill();
}

DEFINE_ART_MODULE(sbnd::MCCaloAna)
