////////////////////////////////////////////////////////////////////////
// Class:       MetricAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        MetricAnalyzer_module.cc
//
// Generated at Fri Nov 22 12:06:48 2024 by Lynn Tung using cetskelgen
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

#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "art_root_io/TFileService.h"
#include "TFile.h"
#include "TTree.h"

namespace sbndaq {
  class MetricAnalyzer;
}


class sbndaq::MetricAnalyzer : public art::EDAnalyzer {
public:
  explicit MetricAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MetricAnalyzer(MetricAnalyzer const&) = delete;
  MetricAnalyzer(MetricAnalyzer&&) = delete;
  MetricAnalyzer& operator=(MetricAnalyzer const&) = delete;
  MetricAnalyzer& operator=(MetricAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  std::string fmetric_instance_name;
  std::string fmetric_module_label;              
  
  std::string fspectdc_module_label;
  uint32_t    fspectdc_etrig_ch;

  TTree* tree;
  int _run, _subrun, _event;
  float  _flash_peakpe;
  float  _flash_peaktime; // us
  float  _trigger_dt;     // us
  int    _timing_type;
  int    _tdc_etrig;
  int    _tdc_bes;
};


sbndaq::MetricAnalyzer::MetricAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fmetric_module_label = p.get<std::string>("metric_module_label","pmtmetricproducer");
  fmetric_instance_name = p.get<std::string>("metric_instance_name","offline");
  fspectdc_module_label = p.get<std::string>("spectdc_module_name","tdcdecoder");
  fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch",4);

  art::ServiceHandle<art::TFileService> fs;
  tree = fs->make<TTree>("tree","metrictree");
  tree->Branch("run",&_run,"run/I");
  tree->Branch("subrun",&_subrun,"subrun/I");
  tree->Branch("event",&_event,"event/I");
  tree->Branch("flash_peakpe",&_flash_peakpe,"flash_peakpe/F");
  tree->Branch("flash_peaktime",&_flash_peaktime,"flash_peaktime/F");
  tree->Branch("trigger_dt",&_trigger_dt,"trigger_dt/F");
  tree->Branch("timing_type",&_timing_type,"timing_type/I");
  tree->Branch("tdc_etrig",&_tdc_etrig,"tdc_etrig/I");
  tree->Branch("tdc_bes",&_tdc_bes,"tdc_bes/I");
}

void sbndaq::MetricAnalyzer::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // initialize tree variables
  _flash_peakpe = -999999999;
  _flash_peaktime = -999999999;
  _trigger_dt = -999999999;
  _timing_type = -1;
  _tdc_etrig = -999999999;
  _tdc_bes = -999999999;

  art::Handle<std::vector<sbnd::trigger::pmtSoftwareTrigger>> metricHandle;
  e.getByLabel(fmetric_module_label,fmetric_instance_name,metricHandle);
  if (!metricHandle.isValid() || metricHandle->size() == 0){
    std::cout << "No Metrics found" << std::endl;
    tree->Fill();
    return;
  }
  else{
    const std::vector<sbnd::trigger::pmtSoftwareTrigger> metric_v(*metricHandle);
    auto metric = metric_v[0];
    _trigger_dt     = float(metric.trig_ts)/1e3;
    _flash_peakpe   = metric.peakPE;
    _flash_peaktime = metric.peaktime;
  }

  // see if etrig exists
  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
  e.getByLabel(fspectdc_module_label,tdcHandle);
  bool found_ett = false;

  if (!tdcHandle.isValid() || tdcHandle->size() == 0)
      std::cout << "No SPECTDC products found." << std::endl;
  else{
      const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
      for (size_t i=0; i<tdc_v.size(); i++){
          auto tdc = tdc_v[i];
          const uint32_t  ch = tdc.Channel();
          if (ch==fspectdc_etrig_ch) found_ett = true;
          if (ch==fspectdc_etrig_ch) _tdc_etrig = tdc.Timestamp()%uint64_t(1e9);
          if (ch==1)                 _tdc_bes   = tdc.Timestamp()%uint64_t(1e9);
      }
  }
  if (found_ett) _timing_type = 0;
  else _timing_type = 1;

  tree->Fill();
}

DEFINE_ART_MODULE(sbndaq::MetricAnalyzer)
