////////////////////////////////////////////////////////////////////////
// Class:       TPCRecomAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        TPCRecomAna_module.cc
//
// Generated at Wed Jul 17 14:05:22 2024 by Thomas Junk using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"

#include "art_root_io/TFileService.h"

#include <TTree.h>

using namespace std;

namespace util {
  class TPCRecomAna;
}

class util::TPCRecomAna : public art::EDAnalyzer {
public:
  explicit TPCRecomAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCRecomAna(TPCRecomAna const&) = delete;
  TPCRecomAna(TPCRecomAna&&) = delete;
  TPCRecomAna& operator=(TPCRecomAna const&) = delete;
  TPCRecomAna& operator=(TPCRecomAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  std::vector<art::Ptr<sim::SimEnergyDeposit>> _edposim_handle;

  struct AnalysisConfig {
    std::vector<std::string> producers;
    std::string instance;

    AnalysisConfig(const fhicl::ParameterSet &param);
    AnalysisConfig() {}
  };

  AnalysisConfig _config;

  
private:

  TTree *fEventTree;
  int _run;
  int _subrun;
  int _event;
  std::vector<float> _edepo;
  std::vector<float> _nele;
  std::vector<float> _angle;
};


util::TPCRecomAna::TPCRecomAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    _config(p)
{
}

void util::TPCRecomAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

  fEventTree = tfs->make<TTree>("EventTree", "event by event info");
  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);
  fEventTree->Branch("edepo", &_edepo);
  fEventTree->Branch("nele", &_nele);
  fEventTree->Branch("angle", &_angle);
}

util::TPCRecomAna::AnalysisConfig::AnalysisConfig(const fhicl::ParameterSet &param) {
  producers = param.get<std::vector<std::string>>("edepo_sim_producers");
  instance = param.get<std::string>("priorSCE", "");
}

void util::TPCRecomAna::analyze(art::Event const& e)
{
  _edposim_handle.clear();

  _run = e.id().run();
  _subrun = e.id().subRun();
  _event = e.id().event();

  _edepo.clear();
  _nele.clear();
  _angle.clear();

  // Get time stamp
  for (const std::string &prod: _config.producers) {
    art::Handle<std::vector<sim::SimEnergyDeposit>> edposim_handle;
    if (_config.instance.size()) {
      e.getByLabel(prod, _config.instance, edposim_handle);
    }
    else {
      e.getByLabel(prod, edposim_handle);
    }

    // exit if the data isn't present
    if (!timestamp_handle.isValid()) {
      std::cerr << "Error: missing timestamps with producer (" << prod << ")" << std::endl;
      return;
    }
    art::fill_ptr_vector(_edposim_handle, edposim_handle);
  }

  // analyze timestamp distributions
  // collect timestamps
  unsigned edeposim_index = 0;
  for (auto const& edposims: _edposim_handle) {
    float this_edpo = edposims->Energy();
    _edepo.push_back(this_edpo);
    edeposim_index++;
  }

  fEventTree->Fill();
}

void util::TPCRecomAna::endJob()
{
  mf::LogVerbatim("TPCRecomAna") << "TPCRecomAna finished job";
}

DEFINE_ART_MODULE(util::TPCRecomAna)
