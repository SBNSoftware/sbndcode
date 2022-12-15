////////////////////////////////////////////////////////////////////////
// Class:       CRTEventDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTEventDisplay_module.cc
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
#include "sbndcode/CRT/CRTEventDisplay/CRTEventDisplayAlg.h"

namespace sbnd::crt {
  class CRTEventDisplay;
}

class sbnd::crt::CRTEventDisplay : public art::EDAnalyzer {
public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
 
    fhicl::Table<CRTEventDisplayAlg::Config> EventDisplayConfig {
      Name("EventDisplayConfig"),
        };
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRTEventDisplay(Parameters const &config);

  CRTEventDisplay(CRTEventDisplay const&) = delete;
  CRTEventDisplay(CRTEventDisplay&&) = delete;
  CRTEventDisplay& operator=(CRTEventDisplay const&) = delete;
  CRTEventDisplay& operator=(CRTEventDisplay&&) = delete;

  void analyze(art::Event const& e) override;

private:

  CRTEventDisplayAlg fCRTEventDisplayAlg;
};


sbnd::crt::CRTEventDisplay::CRTEventDisplay(Parameters const& config)
  : EDAnalyzer{config}
  , fCRTEventDisplayAlg(config().EventDisplayConfig())
  {
  }

void sbnd::crt::CRTEventDisplay::analyze(art::Event const& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  fCRTEventDisplayAlg.Draw(clockData, e);
}

DEFINE_ART_MODULE(sbnd::crt::CRTEventDisplay)
