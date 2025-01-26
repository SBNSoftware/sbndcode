////////////////////////////////////////////////////////////////////////
// Class:       CRTVisualisation
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTVisualisation_module.cc
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

#include "TSystem.h"

namespace sbnd::crt {
  class CRTVisualisation;
}

class sbnd::crt::CRTVisualisation : public art::EDAnalyzer {
public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
 
    fhicl::Table<CRTEventDisplayAlg::Config> EventDisplayConfig {
      Name("EventDisplayConfig"),
        };

    fhicl::Atom<std::string> SaveDir {
      Name("SaveDir"),
        };
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRTVisualisation(Parameters const &config);

  CRTVisualisation(CRTVisualisation const&) = delete;
  CRTVisualisation(CRTVisualisation&&) = delete;
  CRTVisualisation& operator=(CRTVisualisation const&) = delete;
  CRTVisualisation& operator=(CRTVisualisation&&) = delete;

  void analyze(art::Event const& e) override;

private:

  CRTEventDisplayAlg fCRTEventDisplayAlg;
  CRTGeoAlg          fCRTGeoAlg;
  std::string        fSaveDir;
};


sbnd::crt::CRTVisualisation::CRTVisualisation(Parameters const& config)
  : EDAnalyzer{config}
  , fCRTEventDisplayAlg(config().EventDisplayConfig())
  , fCRTGeoAlg(config().EventDisplayConfig().GeoAlgConfig())
  , fSaveDir(config().SaveDir())
  {
    gSystem->Exec(Form("mkdir -p %s", fSaveDir.c_str()));
  }

void sbnd::crt::CRTVisualisation::analyze(art::Event const& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);


  fCRTEventDisplayAlg.SetHighlightedModules({50, 51, 52, 53, 54, 55, 56, 57});
  fCRTEventDisplayAlg.SetSaveFront(true);
  fCRTEventDisplayAlg.SetSaveTop(false);
  fCRTEventDisplayAlg.SetSaveSide(false);

  fCRTEventDisplayAlg.Draw(clockData, e, Form("%s/%s", fSaveDir.c_str(), "South Inner"));
}

DEFINE_ART_MODULE(sbnd::crt::CRTVisualisation)
