////////////////////////////////////////////////////////////////////////
// Class:       GateTimingAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        GateTimingAnalyzer_module.cc
//
// Generated at Mon Apr 14 17:25:17 2025 by Lynn Tung using cetskelgen
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

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "art_root_io/TFileService.h"
#include "TFile.h"
#include "TTree.h"

namespace sbnd {
  class GateTimingAnalyzer;
}


class sbnd::GateTimingAnalyzer : public art::EDAnalyzer {
public:
  explicit GateTimingAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GateTimingAnalyzer(GateTimingAnalyzer const&) = delete;
  GateTimingAnalyzer(GateTimingAnalyzer&&) = delete;
  GateTimingAnalyzer& operator=(GateTimingAnalyzer const&) = delete;
  GateTimingAnalyzer& operator=(GateTimingAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
b
private:

  std::string fspectdc_module_label;
  uint32_t    fspectdc_etrig_ch;

  // Declare member data here.
  TTree* tree;
  int _run, _subrun, _event;
  uint64_t    _tdc_etrig;
  uint64_t    _tdc_bes;
  uint64_t    _tdc_rwm;
};


sbnd::GateTimingAnalyzer::GateTimingAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fspectdc_module_label = p.get<std::string>("spectdc_module_name","tdcdecoder");
  fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch",4);

  art::ServiceHandle<art::TFileService> fs;
  tree = fs->make<TTree>("tree","timingtree");
  tree->Branch("run",&_run,"run/I");
  tree->Branch("subrun",&_subrun,"subrun/I");
  tree->Branch("event",&_event,"event/I");
  tree->Branch("tdc_etrig",&_tdc_etrig,"tdc_etrig/l");
  tree->Branch("tdc_bes",&_tdc_bes,"tdc_bes/l");
  tree->Branch("tdc_rwm",&_tdc_rwm,"tdc_rwm/l");
}

void sbnd::GateTimingAnalyzer::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _tdc_etrig = 0;
  _tdc_bes = 0;
  _tdc_rwm = 0;

  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
  e.getByLabel(fspectdc_module_label,tdcHandle);

  if (!tdcHandle.isValid() || tdcHandle->size() == 0)
      std::cout << "No SPECTDC products found." << std::endl;
  else{
      const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
      for (size_t i=0; i<tdc_v.size(); i++){
          auto tdc = tdc_v[i];
          const uint32_t  ch = tdc.Channel();
          if (ch==fspectdc_etrig_ch) _tdc_etrig = tdc.Timestamp();
          if (ch==1)                 _tdc_bes   = tdc.Timestamp();
          if (ch==2)                 _tdc_rwm   = tdc.Timestamp();
      }
  }
  tree->Fill();
}

DEFINE_ART_MODULE(sbnd::GateTimingAnalyzer)
