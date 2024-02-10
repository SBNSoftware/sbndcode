////////////////////////////////////////////////////////////////////////
// Class:       TPCEventDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        TPCEventDisplay_module.cc
//
// Generated at Thu Oct  6 09:32:09 2022 by Henry Lay using cetskelgen
// from  version .
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

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "sbndcode/TPCEventDisplay/TPCEventDisplayAlg.h"

#include "TSystem.h"

namespace sbnd::crt {
  class TPCEventDisplay;
}

class sbnd::crt::TPCEventDisplay : public art::EDAnalyzer {
public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
 
    fhicl::Table<TPCEventDisplayAlg::Config> EventDisplayConfig {
      Name("EventDisplayConfig"),
        };

    fhicl::Atom<std::string> SaveDir {
      Name("SaveDir"),
        };
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit TPCEventDisplay(Parameters const &config);

  TPCEventDisplay(TPCEventDisplay const&) = delete;
  TPCEventDisplay(TPCEventDisplay&&) = delete;
  TPCEventDisplay& operator=(TPCEventDisplay const&) = delete;
  TPCEventDisplay& operator=(TPCEventDisplay&&) = delete;

  void analyze(art::Event const& e) override;

private:

  TPCEventDisplayAlg fTPCEventDisplayAlg;
  std::string        fSaveDir;
};


sbnd::crt::TPCEventDisplay::TPCEventDisplay(Parameters const& config)
  : EDAnalyzer{config}
  , fTPCEventDisplayAlg(config().EventDisplayConfig())
  , fSaveDir(config().SaveDir())
  {
    gSystem->Exec(Form("mkdir -p %s", fSaveDir.c_str()));
  }

void sbnd::crt::TPCEventDisplay::analyze(art::Event const& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  fTPCEventDisplayAlg.Draw(clockData, e, Form("%s/tpcEventDisplayEvent%i", fSaveDir.c_str(), e.event()));
}

DEFINE_ART_MODULE(sbnd::crt::TPCEventDisplay)
