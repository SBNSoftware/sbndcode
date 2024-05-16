////////////////////////////////////////////////////////////////////////
// Class:       CRTChannelMappingEventDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTChannelMappingEventDisplay_module.cc
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
  class CRTChannelMappingEventDisplay;
}

class sbnd::crt::CRTChannelMappingEventDisplay : public art::EDAnalyzer {
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

  explicit CRTChannelMappingEventDisplay(Parameters const &config);

  CRTChannelMappingEventDisplay(CRTChannelMappingEventDisplay const&) = delete;
  CRTChannelMappingEventDisplay(CRTChannelMappingEventDisplay&&) = delete;
  CRTChannelMappingEventDisplay& operator=(CRTChannelMappingEventDisplay const&) = delete;
  CRTChannelMappingEventDisplay& operator=(CRTChannelMappingEventDisplay&&) = delete;

  void analyze(art::Event const& e) override;

private:

  CRTEventDisplayAlg fCRTEventDisplayAlg;
  CRTGeoAlg          fCRTGeoAlg;
  std::string        fSaveDir;
  std::vector<int>   fChosenTaggers;
};


sbnd::crt::CRTChannelMappingEventDisplay::CRTChannelMappingEventDisplay(Parameters const& config)
  : EDAnalyzer{config}
  , fCRTEventDisplayAlg(config().EventDisplayConfig())
  , fCRTGeoAlg(config().EventDisplayConfig().GeoAlgConfig())
  , fSaveDir(config().SaveDir())
  , fChosenTaggers(config().EventDisplayConfig().ChosenTaggers())
  {
    gSystem->Exec(Form("mkdir -p %s", fSaveDir.c_str()));
  }

void sbnd::crt::CRTChannelMappingEventDisplay::analyze(art::Event const& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  for(auto const& [ name, module ] : fCRTGeoAlg.GetModules())
    {
      if(std::find(fChosenTaggers.begin(), fChosenTaggers.end(), CRTCommonUtils::GetTaggerEnum(module.taggerName)) == fChosenTaggers.end())
        continue;

      fCRTEventDisplayAlg.SetHighlightedModules({module.adID});

      fCRTEventDisplayAlg.Draw(clockData, e, Form("%s/%s", fSaveDir.c_str(), name.c_str()));
    }
}

DEFINE_ART_MODULE(sbnd::crt::CRTChannelMappingEventDisplay)
