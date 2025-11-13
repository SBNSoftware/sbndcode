////////////////////////////////////////////////////////////////////////
// Class:       PMTMCMetricProducer
// Plugin Type: producer (Unknown Unknown)
// File:        PMTMCMetricProducer_module.cc
//
// Generated at Mon Nov 10 11:37:29 2025 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
// An MC version of the PMTMetricProducer (PMT Software Trigger).
// Creates a 10 us "flash" by summing waveforms around the beam spill. 
// Input: raw::OpDetWaveforms from optical simualtion 
// Output:: std::vector<sbnd::trigger::pmtSoftwareTrigger>>
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"

#include <algorithm>
#include <memory>

#include <TTree.h>
#include <string.h>
#include "TH1D.h"
#include "TFile.h"

namespace sbnd {
  namespace trigger {
    class PMTMCMetricProducer;
  }
}

class sbnd::trigger::PMTMCMetricProducer : public art::EDProducer {
public:
  explicit PMTMCMetricProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTMCMetricProducer(PMTMCMetricProducer const&) = delete;
  PMTMCMetricProducer(PMTMCMetricProducer&&) = delete;
  PMTMCMetricProducer& operator=(PMTMCMetricProducer const&) = delete;
  PMTMCMetricProducer& operator=(PMTMCMetricProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::ServiceHandle<art::TFileService> tfs;
  std::vector <std::string> fPDTypes;
  std::string OpDetWaveformsLabel;
  opdet::sbndPDMapAlg fPDSMap;
  float fStartTime;

  float us_to_ticks = 500.;
  float ticks_to_us = 1./us_to_ticks;

  std::vector <int> fChannelsToIgnore;
  std::vector<uint32_t> sumWvfms(const std::vector<uint32_t>& v1, const std::vector<int16_t>& v2);
  float    estimateBaseline(std::vector<uint32_t> wvfm);

  TTree* tree; 
  int _run, _subrun, _event;
  float  _flash_peakpe;
  float  _flash_peaktime;
  int _npmts; 
};


sbnd::trigger::PMTMCMetricProducer::PMTMCMetricProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fPDTypes = p.get<std::vector<std::string>>("PDTypes",{"pmt_coated", "pmt_uncoated"});
  OpDetWaveformsLabel = p.get< std::string >("OpDetWaveformsLabel","opdaq");
  fStartTime = p.get<float>("StartTime",-1);
  fChannelsToIgnore = p.get<std::vector<int>>("ChannelsToIgnore",{});

  art::ServiceHandle<art::TFileService> fs;
  tree = fs->make<TTree>("tree","metrictree");
  tree->Branch("run",&_run,"run/I");
  tree->Branch("subrun",&_subrun,"subrun/I");
  tree->Branch("event",&_event,"event/I");
  tree->Branch("flash_peakpe",&_flash_peakpe,"flash_peakpe/F");
  tree->Branch("flash_peaktime",&_flash_peaktime,"flash_peaktime/F");
  tree->Branch("npmts",&_npmts, "npmts/I");

  // Call appropriate produces<>() functions here.
  produces< std::vector<sbnd::trigger::pmtSoftwareTrigger>>(); 
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::trigger::PMTMCMetricProducer::produce(art::Event& e)
{

  std::unique_ptr<std::vector<sbnd::trigger::pmtSoftwareTrigger>> trig_metrics_v = std::make_unique<std::vector<sbnd::trigger::pmtSoftwareTrigger>>();
  sbnd::trigger::pmtSoftwareTrigger trig_metrics;

  art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
  e.getByLabel(OpDetWaveformsLabel, wvfHandle);
  if(!wvfHandle.isValid() || wvfHandle->size() == 0){
    std::cout << "RawWaveform with label " << OpDetWaveformsLabel << " not found..." << std::endl;
    e.put(std::move(trig_metrics_v));
    return;
  }

  // hardcoded to match match data version
  size_t flash_len=5000;
  std::vector<uint32_t> flash(flash_len, 0);

  float readout_start = fStartTime; 
  float readout_end   = readout_start + flash_len*ticks_to_us;

  int counter=0;
  for(auto const& wvf : (*wvfHandle)) {
    int  wvf_ch = wvf.ChannelNumber();
    auto wvf_ts = wvf.TimeStamp();
    auto wvf_end = wvf_ts+wvf.Waveform().size()*2e-3;
  
    if (wvf_end < readout_start || wvf_ts > readout_end ) continue; 
    if (std::find(fChannelsToIgnore.begin(), fChannelsToIgnore.end(), wvf_ch) != fChannelsToIgnore.end() ) continue;
    if (std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(wvf_ch) ) == fPDTypes.end() ) continue;

    // find where readout_start is
    if (wvf_ts <= readout_start && wvf_end >= readout_end){
      int offset = int((readout_start - wvf_ts)*us_to_ticks);
      std::vector<short> wvfm_snippet(wvf.Waveform().begin() + offset, wvf.Waveform().begin() + offset + flash_len);
      flash = sumWvfms(flash, wvfm_snippet);
      counter++;
    }
  }

  auto flash_baseline = estimateBaseline(flash);
  auto flash_peak_it = std::min_element(flash.begin(),flash.end());
  _flash_peakpe   = (flash_baseline-(*flash_peak_it))/12.5;
  _flash_peaktime = ((flash_peak_it-flash.begin()))*ticks_to_us + readout_start; 
  _npmts = counter;

  trig_metrics.foundBeamTrigger = true;
  trig_metrics.nAboveThreshold = _npmts;
  trig_metrics.promptPE = 0; 
  trig_metrics.prelimPE = 0;
  trig_metrics.peakPE = _flash_peakpe; 
  trig_metrics.peaktime = _flash_peaktime; 
  trig_metrics.trig_ts = 0;
  trig_metrics_v->push_back(trig_metrics);

  e.put(std::move(trig_metrics_v));

  tree->Fill();
}


std::vector<uint32_t> sbnd::trigger::PMTMCMetricProducer::sumWvfms(const std::vector<uint32_t>& v1, 
                                                                   const std::vector<int16_t>& v2) 
{
  size_t result_len = (v1.size() > v2.size()) ? v1.size() : v2.size();
  std::vector<uint32_t>  result(result_len,0);
  for (size_t i = 0; i < result_len; i++){
    auto value1 = (i < v1.size()) ? v1[i] : 1e5;
    auto value2 = (i < v2.size()) ? v2[i] : 1e5;
    result.at(i) = value1 + uint32_t(value2);
  }
  return result;
}


float sbnd::trigger::PMTMCMetricProducer::estimateBaseline(std::vector<uint32_t> wvfm){
  // use a copy of 'wvfm' because nth_element might do something weird!
  const auto median_it = wvfm.begin() + wvfm.size() / 2;
  std::nth_element(wvfm.begin(), median_it , wvfm.end());
  auto median = *median_it;
  return median;
}

DEFINE_ART_MODULE(sbnd::trigger::PMTMCMetricProducer)
