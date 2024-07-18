////////////////////////////////////////////////////////////////////////
// Class:       TPCChoppyAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        TPCChoppyAna_module.cc
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

#include "art_root_io/TFileService.h"

#include <TTree.h>

using namespace std;

namespace util {
  class TPCChoppyAna;
}

class util::TPCChoppyAna : public art::EDAnalyzer {
public:
  explicit TPCChoppyAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCChoppyAna(TPCChoppyAna const&) = delete;
  TPCChoppyAna(TPCChoppyAna&&) = delete;
  TPCChoppyAna& operator=(TPCChoppyAna const&) = delete;
  TPCChoppyAna& operator=(TPCChoppyAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  std::vector<art::Ptr<raw::RDTimeStamp>> _raw_timestamps_handle;

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
  long _delta_tmaxmin;
  std::vector<long> _rdts;
  
};


util::TPCChoppyAna::TPCChoppyAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    _config(p)
{
}

void util::TPCChoppyAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

  fEventTree = tfs->make<TTree>("EventTree", "event by event info");
  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);
  fEventTree->Branch("delta_tmaxmin", &_delta_tmaxmin);
  fEventTree->Branch("rdts", "std::vector<long>", &_rdts);
}

util::TPCChoppyAna::AnalysisConfig::AnalysisConfig(const fhicl::ParameterSet &param) {
  producers = param.get<std::vector<std::string>>("raw_digit_producers");
  instance = param.get<std::string>("raw_digit_instance", "");
}

void util::TPCChoppyAna::analyze(art::Event const& e)
{
  _raw_timestamps_handle.clear();

  _run = e.id().run();
  _subrun = e.id().subRun();
  _event = e.id().event();

  _rdts.clear();

  // Get time stamp
  for (const std::string &prod: _config.producers) {
    art::Handle<std::vector<raw::RDTimeStamp>> timestamp_handle;
    if (_config.instance.size()) {
      e.getByLabel(prod, _config.instance, timestamp_handle);
    }
    else {
      e.getByLabel(prod, timestamp_handle);
    }

    // exit if the data isn't present
    if (!timestamp_handle.isValid()) {
      std::cerr << "Error: missing timestamps with producer (" << prod << ")" << std::endl;
      return;
    }
    art::fill_ptr_vector(_raw_timestamps_handle, timestamp_handle);
  }

  // analyze timestamp distributions
  // collect timestamps
  unsigned timeindex = 0;
  for (auto const& timestamps: _raw_timestamps_handle) {
    long this_ts = timestamps->GetTimeStamp();
    _rdts.push_back(this_ts);
    timeindex++;
  }

  // make time stamp related metrics
  long ts_min = _rdts[0];
  long ts_max = _rdts[0];
  for (long ts : _rdts) {
    if (ts < ts_min) {
      ts_min = ts;
    }
    if (ts > ts_max) {
      ts_max = ts;
    }
  }
  _delta_tmaxmin = ts_max - ts_min;

  cout << "ts_max : " << ts_max << ", ts_min : " << ts_min << endl;
  fEventTree->Fill();
}

void util::TPCChoppyAna::endJob()
{
  mf::LogVerbatim("TPCChoppyAna") << "TPCChoppyAna finished job";
}

DEFINE_ART_MODULE(util::TPCChoppyAna)
