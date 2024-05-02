////////////////////////////////////////////////////////////////////////
// Class:       pmtSoftwareTriggerProducer
// Plugin Type: producer (Unknown Unknown)
// File:        pmtSoftwareTriggerProducer_module.cc
//
// Generated at Thu Feb 17 13:22:51 2022 by Patrick Green using cetskelgen
// from  version .

// Module to implement software trigger metrics to the PMT Trigger simulation
// Input: artdaq fragment output from the pmtArtdaqFragmentProducer.cc module
// Calculates various PMT metrics for every event (that passes the hardware trigger)
// Output: sbnd::trigger::pmtSoftwareTrigger data product 


// Notes: 
// in simulation, the trigger time tag **points at the actual flash trigger time**, not the end of the waveform

// More information can be found at:
// https://sbnsoftware.github.io/sbndcode_wiki/SBND_Trigger
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
#include "art_root_io/TFileService.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <algorithm>
#include <valarray>
#include <numeric>

namespace sbnd {
  namespace trigger {
    class pmtSoftwareTriggerProducer;
  }
}

class sbnd::trigger::pmtSoftwareTriggerProducer : public art::EDProducer {
public:
  explicit pmtSoftwareTriggerProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pmtSoftwareTriggerProducer(pmtSoftwareTriggerProducer const&) = delete;
  pmtSoftwareTriggerProducer(pmtSoftwareTriggerProducer&&) = delete;
  pmtSoftwareTriggerProducer& operator=(pmtSoftwareTriggerProducer const&) = delete;
  pmtSoftwareTriggerProducer& operator=(pmtSoftwareTriggerProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  // fhicl parameters
  // art::Persistable is_persistable_;
  int  fMultiplicityThreshold;

  double fTriggerTimeOffset;    // offset of trigger time, default 0.5 sec
  double fWindowLength; // beam window length after trigger time, default 1.6us
  uint32_t fWvfmLength;
  int fVerbose;
  bool fSaveSumHist; 
  bool fSavePMTHists; 

  bool fCalculateBaseline;
  bool fCountPMTs;
  bool fCalculatePEMetrics;
  bool fFindFlashInfo;
  bool fFindPulses;

  std::vector<double> fInputBaseline; 
  // double fBaselineWindow;
  double fPromptWindow;
  double fPrelimWindow; 

  double fADCThreshold;
  double fADCtoPE;       
  double fPEArea;    

  // histogram info  
  std::stringstream histname; //raw waveform hist name
  art::ServiceHandle<art::TFileService> tfs;

  // event information
  int fRun, fSubrun;
  art::EventNumber_t fEvent;

  // PD information
  opdet::sbndPDMapAlg pdMap; // photon detector map
  std::vector<unsigned int> channelList; 

  // beam window
  // set in artdaqFragment producer, in reality would be provided by event builder
  bool foundBeamTrigger;
  double ticks_to_us;
  double us_to_ticks;
  uint32_t beamWindowStart;
  uint32_t beamWindowEnd;

  // waveforms
  uint32_t fTriggerTime;
  bool fWvfmsFound;
  std::vector<std::vector<uint16_t>> fWvfmsVec;

  // pmt information 
  std::vector<sbnd::trigger::pmtInfo> fpmtInfoVec;

  void checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag);
  void analyzeCAEN1730Fragment(const artdaq::Fragment &frag);
  std::vector<uint32_t> sumWvfms(const std::vector<uint32_t>& v1, const std::vector<uint16_t>& v2);
  // void   estimateBaseline(int i_ch);
  double estimateBaseline(std::vector<uint32_t> wvfm);
  double estimateBaseline(std::vector<uint16_t> wvfm);
  void SimpleThreshAlgo(int i_ch);

  TTree* _tree; 
  int _run, _sub, _evt; 
  bool   _beam_trig; 
  double _time_trig; 
  int    _npmt;
  double _beam_promptPE, _beam_prelimPE, _beam_peakPE;  // prompt (sum) PE and preliminary (sum) PE for beam window
  double _trig_promptPE, _trig_prelimPE, _trig_peakPE; // prompt (sum) PE and preliminary (sum) PE for trigger window
  double _flash_promptPE, _flash_prelimPE , _flash_peakPE; // prompt (sum) PE and preliminary (sum) PE for flash window
  double _flash_peaktime;

  std::vector<int> _hardware_trig_pairs;
  int    _hardware_trig_idx;
  int    _hardware_trig_npairs;
  int    _hardware_trig_maxPMTs;
  int    _hardware_trig_maxPMTs_idx;

  TTree* _pulse_tree;
  int _npulses; 
  
  std::vector<int> _pulse_ch; // ch number for each pulse 
  std::vector<double> _pulse_t_start; // t_start for each pulse 
  std::vector<double> _pulse_t_end;
  std::vector<double> _pulse_t_peak;
  std::vector<double> _pulse_peak;
  std::vector<double> _pulse_area;
};


sbnd::trigger::pmtSoftwareTriggerProducer::pmtSoftwareTriggerProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fMultiplicityThreshold(p.get<int>("MultiplicityThreshold", 1)),
  fTriggerTimeOffset (p.get<double>("TriggerTimeOffset", 0.5)),
  fWindowLength      (p.get<double>("WindowLength", 1.6)),
  fWvfmLength        (p.get<uint32_t>("WvfmLength", 5120)),

  fVerbose           (p.get<int>("Verbose", 1)),
  fSaveSumHist         (p.get<bool>("SaveSumHist",true)),
  fSavePMTHists        (p.get<bool>("SavePMTHists",false)),
  
  fCalculateBaseline (p.get<bool>("CalculateBaseline",true)),
  fCountPMTs         (p.get<bool>("CountPMTs",true)),
  fCalculatePEMetrics(p.get<bool>("CalculatePEMetrics",false)),
  fFindFlashInfo     (p.get<bool>("FindFlashInfo",true)),
  fFindPulses        (p.get<bool>("FindPulses", false)),
  
  fInputBaseline     (p.get<std::vector<double>>("InputBaseline")),
  // fBaselineWindow    (p.get<double>("BaselineWindow", 0.5)),
  fPromptWindow      (p.get<double>("PromptWindow", 0.5)),
  fPrelimWindow      (p.get<double>("PrelimWindow", 0.5)),
  fADCThreshold      (p.get<double>("ADCThreshold", 7960)),
  fADCtoPE           (p.get<double>("ADCtoPE", 1)),
  fPEArea            (p.get<double>("PEArea", 66.33))

{
  // Call appropriate produces<>() functions here.
  // produces< sbnd::trigger::pmtSoftwareTrigger >("", is_persistsable_);
  produces<std::vector<sbnd::trigger::pmtSoftwareTrigger>>();

  beamWindowStart = fTriggerTimeOffset*1e9;
  beamWindowEnd = beamWindowStart + fWindowLength*1000;

  // std::cout << "beam window start: " << beamWindowStart << std::endl;
  // std::cout << "beam window end: " << beamWindowEnd << std::endl;
  // std::cout << "trigger time offset: " << fTriggerTimeOffset << std::endl;
  std::cout << "PMT Software Trigger Window length: " << fWindowLength << " us " << std::endl;

  us_to_ticks = 500.; // ticks per us 
  ticks_to_us = 1/us_to_ticks; // us per ticks

  // build PD map and channel list
  auto subsetCondition = [](auto const& i)->bool { return i["pd_type"] == "pmt_coated" || i["pd_type"] == "pmt_uncoated"; };
  auto pmtMap = pdMap.getCollectionFromCondition(subsetCondition);
  for(auto const& i:pmtMap){
    channelList.push_back(i["channel"]);
  }
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("software_metrics_tree","");
  _tree->Branch("run",       &_run,       "run/I");
  _tree->Branch("sub",       &_sub,       "sub/I");
  _tree->Branch("evt",       &_evt,       "evt/I");
  _tree->Branch("beam_trig", &_beam_trig, "beam_trig/O");
  _tree->Branch("time_trig", &_time_trig, "time_trig/D");
  _tree->Branch("npmt",      &_npmt,      "npmt/I");
  _tree->Branch("beam_promptPE",  &_beam_promptPE,  "beam_promptPE/D");
  _tree->Branch("beam_prelimPE",  &_beam_prelimPE,  "beam_prelimPE/D");
  _tree->Branch("beam_peakPE",    &_beam_peakPE,    "beam_peakPE/D");
  _tree->Branch("trig_promptPE",  &_trig_promptPE,  "trig_promptPE/D");
  _tree->Branch("trig_prelimPE",  &_trig_prelimPE,  "trig_prelimPE/D");
  _tree->Branch("trig_peakPE",    &_trig_peakPE,    "trig_peakPE/D");
  _tree->Branch("flash_promptPE", &_flash_promptPE, "flash_promptPE/D");
  _tree->Branch("flash_prelimPE", &_flash_prelimPE, "flash_prelimPE/D");
  _tree->Branch("flash_peakPE",   &_flash_peakPE,   "flash_peakPE/D");
  _tree->Branch("flash_peaktime", &_flash_peaktime, "flash_peaktime/D");

  _tree->Branch("hardware_trig_pairs",  "std::vector<int>", &_hardware_trig_pairs);
  _tree->Branch("hardware_trig_idx", &_hardware_trig_idx, "hardware_trig_idx/I");
  _tree->Branch("hardware_trig_npairs", &_hardware_trig_npairs, "hardware_trig_npairs/I");
  _tree->Branch("hardware_trig_maxPMTs", &_hardware_trig_maxPMTs, "hardware_trig_maxPMTs/I");
  _tree->Branch("hardware_trig_maxPMTs_idx", &_hardware_trig_maxPMTs_idx, "hardware_trig_maxPMTs_idx/I");

  _pulse_tree = fs->make<TTree>("pulse_tree","");
  _pulse_tree->Branch("run",       &_run,       "run/I");
  _pulse_tree->Branch("sub",       &_sub,       "sub/I");
  _pulse_tree->Branch("evt",       &_evt,       "evt/I");
  _pulse_tree->Branch("npulses",   &_npulses,   "npulses/I");
  _pulse_tree->Branch("pulse_ch",      "std::vector<int>",    &_pulse_ch);
  _pulse_tree->Branch("pulse_t_start", "std::vector<double>", &_pulse_t_start);
  _pulse_tree->Branch("pulse_t_end",   "std::vector<double>", &_pulse_t_end);
  _pulse_tree->Branch("pulse_t_peak",  "std::vector<double>", &_pulse_t_peak);
  _pulse_tree->Branch("pulse_peak",    "std::vector<double>", &_pulse_peak);
  _pulse_tree->Branch("pulse_area",    "std::vector<double>", &_pulse_area);

}

void sbnd::trigger::pmtSoftwareTriggerProducer::produce(art::Event& e)
{
  // Implementation of required member function here.
  // event information
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();
  
  _run = e.run();
  _sub = e.subRun();
  _evt = e.id().event();

  if (fVerbose>=1) std::cout << "Processing Run: " << fRun << ", Subrun: " << fSubrun << ", Event: " << fEvent << std::endl;

  // reset for this event
  foundBeamTrigger = false;
  fWvfmsFound = false;
  fWvfmsVec.clear(); fWvfmsVec.resize(15*8); // 15 pmt channels per fragment, 8 fragments per trigger
  fpmtInfoVec.clear(); fpmtInfoVec.resize(15*8); 

  _beam_trig = false;
  _time_trig = -9999;
  _npmt = -9999;
  _beam_promptPE = -9999; _beam_prelimPE = -9999; _beam_peakPE = -9999;
  _trig_promptPE = -9999; _trig_prelimPE = -9999; _trig_peakPE = -9999;
  _flash_promptPE = -9999; _flash_prelimPE = -9999; _flash_peakPE = -9999;
  _flash_peaktime = -9999;

  // access hardware trigger information
  _hardware_trig_idx = -1; _hardware_trig_npairs = -1; _hardware_trig_maxPMTs = -1; _hardware_trig_maxPMTs_idx = -1;
  _hardware_trig_pairs.clear();
  art::Handle< std::vector< sbnd::comm::pmtTrigger > > triggerHandle;
  e.getByLabel("pmttriggerproducer", triggerHandle);
  std::vector<size_t> triggerIndex; 
  for(auto const& trigger : (*triggerHandle)) {
    _hardware_trig_pairs.resize(trigger.numPassed.size());
    _hardware_trig_maxPMTs = trigger.maxPMTs; 
    for (size_t idx = 0; idx < trigger.numPassed.size(); idx++) {
      _hardware_trig_pairs[idx] = trigger.numPassed[idx];
      if (trigger.numPassed[idx] == _hardware_trig_maxPMTs)
        _hardware_trig_maxPMTs_idx = idx;
      if (trigger.numPassed[idx] >= fMultiplicityThreshold) {
        _hardware_trig_idx = idx;
        _hardware_trig_npairs = trigger.numPassed[idx];
      }
    }
  } // trigger handle loop

  // get fragment handles
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  // loop over fragment handles
  for (auto &handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;
    if (handle->front().type()==sbndaq::detail::FragmentType::CAENV1730) {
      if (fVerbose>=2)   std::cout << "Found " << handle->size() << " CAEN1730 fragments" << std::endl;
      
      // identify whether any fragments correspond to the beam spill
      // loop over fragments, in steps of 8
      size_t beamFragmentIdx = 9999;
      for (size_t fragmentIdx = 0; fragmentIdx < handle->size(); fragmentIdx += 8) {
        checkCAEN1730FragmentTimeStamp(handle->at(fragmentIdx));
        if (foundBeamTrigger) {
          beamFragmentIdx = fragmentIdx;
          if (fVerbose>=2) std::cout << "Found fragment in time with beam at index: " << beamFragmentIdx << std::endl;
          break;
        }
      }
      // if set of fragment in time with beam found, process waveforms
      if (foundBeamTrigger && beamFragmentIdx != 9999) {
        for (size_t fragmentIdx = beamFragmentIdx; fragmentIdx < beamFragmentIdx+8; fragmentIdx++) {
          analyzeCAEN1730Fragment(handle->at(fragmentIdx));
        }
        fWvfmsFound = true;
      }
    }
  } // end loop over handles

  // object to store trigger metrics in
  std::unique_ptr<std::vector<sbnd::trigger::pmtSoftwareTrigger>> trig_metrics_v = std::make_unique<std::vector<sbnd::trigger::pmtSoftwareTrigger>>();
  sbnd::trigger::pmtSoftwareTrigger trig_metrics;

  if (foundBeamTrigger && fWvfmsFound) {

    trig_metrics.foundBeamTrigger = true;
    _beam_trig = true;

    // store timestamp of trigger, relative to beam window start
    double trig_ts = fTriggerTime - beamWindowStart;
    trig_metrics.triggerTimestamp = trig_ts;
    // trig_metrics.trig_ts = trig_ts;
    _time_trig = trig_ts;
    if (fVerbose>=1) std::cout << "Saving trigger timestamp: " << trig_ts << " ns" << std::endl;

    _beam_promptPE = 0; _beam_prelimPE = 0; _beam_peakPE = 0;
    _trig_promptPE = 0; _trig_prelimPE = 0; _trig_peakPE = 0;

    _flash_promptPE = 0; _flash_prelimPE = 0; _flash_peakPE = 0;
    _flash_peaktime = 0;

    int nAboveThreshold = 0;
    // find the waveform bins that correspond to the start and end of the extended spill window (0 -> 1.8 us) within the 10 us waveform 
    // !! if the triggerTimeStamp > 1000, the beginning of the beam spill is *not* contained within the waveform 
    int beamStartBin = (trig_ts >= 1000)? 0 : int(500 - (trig_ts)/2); // units of bins 
    int beamEndBin   = (trig_ts >= 1000)? (fWindowLength*1e3 - (trig_ts-1000) )/2 : (beamStartBin + (fWindowLength*1e3)/2);
    
    // wvfm loop to calculate metrics 
    _npulses = 0;
    _pulse_ch.clear(); //_pulse_npulses.clear();
    _pulse_t_start.clear(); _pulse_t_end.clear(); _pulse_t_peak.clear();
    _pulse_peak.clear(); _pulse_area.clear();

    _pulse_ch.reserve(1000); //_pulse_npulses.reserve(1000);
    _pulse_t_start.reserve(1000); _pulse_t_end.reserve(1000); _pulse_t_peak.reserve(1000);
    _pulse_peak.reserve(1000); _pulse_area.reserve(1000);

    std::vector<uint32_t> wvfm_sum(fWvfmLength, 0);
    auto prelim_offset = (beamStartBin - fPrelimWindow*us_to_ticks) < 0 ? 0 : beamStartBin - fPrelimWindow*us_to_ticks;
    for (int i_ch = 0; i_ch < 120; ++i_ch){
      auto &pmtInfo = fpmtInfoVec.at(i_ch);
      auto wvfm = fWvfmsVec[i_ch];

      if (wvfm.begin() == wvfm.end()) continue;
      // assign channel 
      pmtInfo.channel = channelList.at(i_ch);
      // count number of PMTs above threshold within the beam window
      if (fCountPMTs){
        for (int bin = beamStartBin; bin < beamEndBin; ++bin){
          auto adc = wvfm[bin];
          if (adc < fADCThreshold){ 
            nAboveThreshold++; 
            break;
          }
        }
      }
      else {nAboveThreshold=-9999;}

      // quick estimate prompt and preliminary light, assuming sampling rate of 500 MHz (2 ns per bin)
      if (fCalculatePEMetrics){
        // calculate baseline
        if (fCalculateBaseline){
          pmtInfo.baseline = estimateBaseline(wvfm);
          pmtInfo.baselineSigma = 0;
          }
        else { 
          pmtInfo.baseline=fInputBaseline.at(0); 
          pmtInfo.baselineSigma = fInputBaseline.at(1); 
          }

        double baseline = pmtInfo.baseline;
        
        auto beam_prompt_window = std::vector<uint16_t>(wvfm.begin()+beamStartBin, wvfm.begin()+beamStartBin+fPromptWindow*us_to_ticks);
        auto beam_prelim_window = std::vector<uint16_t>(wvfm.begin()+prelim_offset, wvfm.begin()+beamStartBin);

        auto trig_prompt_window = std::vector<uint16_t>(wvfm.begin()+500, wvfm.begin()+500+fPromptWindow*us_to_ticks);
        auto trig_prelim_window = std::vector<uint16_t>(wvfm.begin(), wvfm.begin()+500);

        if (fFindPulses == false){
          double beam_peakPE_ch = (baseline-(*std::min_element(beam_prompt_window.begin(), beam_prompt_window.end())))/fADCtoPE;
          double beam_promptPE_ch = baseline*beam_prompt_window.size() - std::accumulate(beam_prompt_window.begin(), beam_prompt_window.end(), 0);
          double beam_prelimPE_ch = baseline*beam_prelim_window.size() - std::accumulate(beam_prelim_window.begin(), beam_prelim_window.end(), 0);
          _beam_peakPE += beam_peakPE_ch;
          _beam_promptPE += beam_promptPE_ch;
          _beam_prelimPE += beam_prelimPE_ch;

          double trig_peakPE_ch = (baseline-(*std::min_element(trig_prompt_window.begin(), trig_prompt_window.end())))/fADCtoPE;
          double trig_promptPE_ch = baseline*trig_prompt_window.size() - std::accumulate(trig_prompt_window.begin(), trig_prompt_window.end(), 0);
          double trig_prelimPE_ch = baseline*trig_prelim_window.size() - std::accumulate(trig_prelim_window.begin(), trig_prelim_window.end(), 0);

          _trig_peakPE += trig_peakPE_ch;
          _trig_promptPE += trig_promptPE_ch;
          _trig_prelimPE += trig_prelimPE_ch;
          }
      // pulse finder + prompt and prelim calculation with pulses 
        if (fFindPulses == true){
          int pulse_counter = 0;
          std::cout << "Finding pulses for PMT #" << channelList.at(i_ch) << std::endl;
          SimpleThreshAlgo(i_ch);

          for (auto pulse : pmtInfo.pulseVec){
            pulse_counter++; 
            // times in pulse.t_**** are in units of *bins* not actually time 
            if(pulse.t_start > beamStartBin && pulse.t_end < beamEndBin) {
              std::cout << "pulse " << pulse_counter << std::endl;
              std::cout << "pulse start, peak, end: (" << (pulse.t_start - beamStartBin)*2 << ", " 
                                                      << (pulse.t_peak  - beamStartBin)*2 << ", " 
                                                      << (pulse.t_end   - beamStartBin)*2 << ") " 
                                                      << std::endl;
              std::cout << "pulse peak ADC (baseline subtracted): " << pulse.peak << std::endl;
              std::cout << "pulse area ADC (baseline subtracted): " << pulse.area << std::endl;
            }
            // if (pulse.t_start > 500 && pulse.t_end < 550) _beam_promptPE+=pulse.pe;
            // if ((trig_ts) >= 1000){ if (pulse.t_end < 500) _beam_prelimPE+=pulse.pe; }
            // else if (trig_ts < 1000){
            //   if (pulse.t_start > (500 - abs((trig_ts-1000)/2)) && pulse.t_end < 500) _beam_prelimPE+=pulse.pe; 
            // }
            _npulses++;
            _pulse_ch.push_back(channelList.at(i_ch));
            _pulse_t_start.push_back((pulse.t_start - beamStartBin)*2);
            _pulse_t_end.push_back(  (pulse.t_end   - beamStartBin)*2);
            _pulse_t_peak.push_back( (pulse.t_peak  - beamStartBin)*2);
            _pulse_peak.push_back(pulse.peak);
            _pulse_area.push_back(pulse.area);
          } // end of pulse loop 
        } // end of pulse finding if statement
      } // end of if statement for calculating PE metrics
      if (fFindFlashInfo) wvfm_sum = sumWvfms(wvfm_sum, wvfm);
    } // end of wvfm loop 
    if (fFindFlashInfo){
      // find and remove the baseline
      double flash_baseline = 0;
      if (fCalculateBaseline) flash_baseline = estimateBaseline(wvfm_sum);
      else flash_baseline = fInputBaseline.at(0)*120;

      auto flash_prompt_window = std::vector<uint32_t>(wvfm_sum.begin()+500, wvfm_sum.begin()+500+fPromptWindow*us_to_ticks);
      auto flash_prelim_window = std::vector<uint32_t>(wvfm_sum.begin(), wvfm_sum.begin()+500);

      _flash_peakPE = (flash_baseline-(*std::min_element(flash_prompt_window.begin(), flash_prompt_window.end())))/fADCtoPE;
      _flash_promptPE = flash_baseline*flash_prompt_window.size() - std::accumulate(flash_prompt_window.begin(), flash_prompt_window.end(), 0);
      _flash_prelimPE = flash_baseline*flash_prelim_window.size() - std::accumulate(flash_prelim_window.begin(), flash_prelim_window.end(), 0);
      auto flash_peak_it = std::min_element(wvfm_sum.begin()+500,  wvfm_sum.begin()+500+fPromptWindow*us_to_ticks);
      _flash_peaktime = (std::distance(wvfm_sum.begin(), flash_peak_it))*ticks_to_us + trig_ts;
    }
    _pulse_tree->Fill();

    trig_metrics.nAboveThreshold = nAboveThreshold;
    // ! need to add fcl parameters for how we want to fill these data product entries: beam vs. trigger vs. flash
    trig_metrics.promptPE = _flash_peakPE;
    trig_metrics.prelimPE = _flash_prelimPE;
    // trig_metrics.peakPE   = _flash_peakPE;
    // trig_metrics.peaktime = _flash_peaktime;
    trig_metrics.foundBeamTrigger = false;
    trig_metrics.triggerTimestamp = trig_ts;
    // trig_metrics.trig_ts = -9999;
    trig_metrics.nAboveThreshold = -9999;
    trig_metrics.promptPE = -9999;
    trig_metrics.prelimPE = -9999;
    // trig_metrics.peakPE = -9999;
    // trig_metrics.peaktime = -9999;


    // tree variables - TTree which is only part of root 
    _npmt = nAboveThreshold;

    if (fVerbose>=1 && fCountPMTs) std::cout << "nPMTs Above Threshold: " << nAboveThreshold << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "beam prompt pe: " << _beam_promptPE << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "beam prelim pe: " << _beam_prelimPE << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "beam peak pe: " << _beam_peakPE << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "trig prompt pe: " << _trig_promptPE << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "trig prelim pe: " << _trig_prelimPE << std::endl;
    if (fVerbose>=1 && fCalculatePEMetrics) std::cout << "trig peak pe: " << _trig_peakPE << std::endl;
    if (fVerbose>=1 && fFindFlashInfo) std::cout << "flash prompt pe: " << _flash_promptPE << std::endl;
    if (fVerbose>=1 && fFindFlashInfo) std::cout << "flash prelim pe: " << _flash_prelimPE << std::endl;
    if (fVerbose>=1 && fFindFlashInfo) std::cout << "flash peak pe: " << _flash_peakPE << std::endl;
    if (fVerbose>=1 && fFindFlashInfo) std::cout << "flash peak time: " << _flash_peaktime << std::endl;

    // start histo 
    if (fSavePMTHists == true){
      int hist_id = -1; 
      for (size_t i_wvfm = 0; i_wvfm < fWvfmsVec.size(); ++i_wvfm){
        std::vector<uint16_t> wvfm = fWvfmsVec[i_wvfm];
        hist_id++;
        //if (fEvent<4){
            histname.str(std::string());
            histname << "run_" << fRun  
                    << "_subrun_" <<fSubrun
                    << "_event_" << fEvent
                    << "_pmtnum_" << channelList.at(i_wvfm);
            // assuming that we save ~1 us before the triggerTimeStamp  
            double StartTime = (trig_ts-1000)*1e-3; // us
            double EndTime   = StartTime + (fWvfmLength*2)*1e-3;

            TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), "Raw Waveform", wvfm.size(), StartTime, EndTime);
            wvfmHist->GetXaxis()->SetTitle("t (#mus)");
            for(unsigned int i = 0; i < wvfm.size(); i++) {
              wvfmHist->SetBinContent(i + 1, (double)wvfm[i]);
            }
        //} 
      } // end histo
    }
    if (fSaveSumHist == true){
      histname.str(std::string());
      histname << "run_" << fRun  
               << "_subrun_" <<fSubrun
               << "_event_" << fEvent
               << "_sum" ;
      TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), "Raw Summed Waveform", wvfm_sum.size(), 0, wvfm_sum.size()*2*1e-3);
      wvfmHist->GetXaxis()->SetTitle("t (#mus)");
      for(unsigned int i = 0; i < wvfm_sum.size(); i++) {
        wvfmHist->SetBinContent(i + 1, (double)wvfm_sum[i]);
      }
    }
  }
  else{
    if (fVerbose>=1) std::cout << "Beam and wvfms not found" << std::endl;
    trig_metrics.foundBeamTrigger = false;
    // trig_metrics.trig_ts = -9999;
    trig_metrics.triggerTimestamp = -9999;
    trig_metrics.nAboveThreshold = -9999;
    trig_metrics.promptPE = -9999;
    trig_metrics.prelimPE = -9999;
    // trig_metrics.peakPE = -9999;
    // trig_metrics.peaktime = -9999;
  }
  
  trig_metrics_v->push_back(trig_metrics);
  e.put(std::move(trig_metrics_v));    
  _tree->Fill();
}

void sbnd::trigger::pmtSoftwareTriggerProducer::checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag) {

  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp
  uint32_t timestamp = md->timeStampNSec;

  // access beam signal, in ch15 of first PMT of each fragment set
  // check entry 500 (0us), at trigger time
  const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
                 + sizeof(sbndaq::CAENV1730EventHeader));
  const uint16_t* value_ptr =  data_begin;
  uint16_t value = 0;

  size_t ch_offset = (size_t)(15*fWvfmLength);
  size_t tr_offset = fTriggerTimeOffset*1e3;

  value_ptr = data_begin + ch_offset + tr_offset; // pointer arithmetic 
  value = *(value_ptr);

  if (value == 1 && timestamp >= beamWindowStart && timestamp <= beamWindowEnd) {
    foundBeamTrigger = true;
    fTriggerTime = timestamp;
  }
}

void sbnd::trigger::pmtSoftwareTriggerProducer::analyzeCAEN1730Fragment(const artdaq::Fragment &frag) {
  
  // access fragment ID; index of fragment out of set of 8 fragments
  int fragId = static_cast<int>(frag.fragmentID()); 

  // access waveforms in fragment and save
  const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
                 + sizeof(sbndaq::CAENV1730EventHeader));
  const uint16_t* value_ptr =  data_begin;
  uint16_t value = 0;

  // channel offset
  size_t nChannels = 15; // 15 pmts per fragment
  size_t ch_offset = 0;

  // loop over channels
  for (size_t i_ch = 0; i_ch < nChannels; ++i_ch){
    fWvfmsVec[i_ch + nChannels*fragId].resize(fWvfmLength);
    ch_offset = (size_t)(i_ch * fWvfmLength);
    //--loop over waveform samples
    for(size_t i_t = 0; i_t < fWvfmLength; ++i_t){ 
      value_ptr = data_begin + ch_offset + i_t; // pointer arithmetic
      value = *(value_ptr);
      fWvfmsVec[i_ch + nChannels*fragId][i_t] = value;
    } //--end loop samples
  } //--end loop channels
}

std::vector<uint32_t> sbnd::trigger::pmtSoftwareTriggerProducer::sumWvfms(const std::vector<uint32_t>& v1, const std::vector<uint16_t>& v2) {
    std::vector<uint32_t>  result(v1.size(),0);
    for (size_t i = 0; i < v1.size(); i++)
        result.at(i) = (v1[i] + v2[i]);
    return result;
}

// void sbnd::trigger::pmtSoftwareTriggerProducer::estimateBaseline(int i_ch){
//   auto wvfm = fWvfmsVec[i_ch];
//   auto &pmtInfo = fpmtInfoVec[i_ch];
//   // assuming that the first 500 ns doesn't include peaks, find the mean of the ADC count as the baseline
//   std::vector<uint16_t> subset = std::vector<uint16_t>(wvfm.begin(),  wvfm.begin()+int(fBaselineWindow*us_to_ticks));
//   double subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size());
//   double val = 0;
//   for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
//   double subset_stddev = sqrt(val/subset.size());

//   if (subset_stddev > 3){ // make this fcl parameter?
//     val = 0; subset.clear(); subset_stddev = 0;
//     subset = std::vector<uint16_t>(wvfm.end()-int(fBaselineWindow*us_to_ticks), wvfm.end());
//     subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size());
//     for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
//     subset_stddev = sqrt(val/subset.size());
//   }
//   if (fVerbose>=3) std::cout << "Baseline mean: " << subset_mean << " ADC, with " << subset_stddev << " sigma" << std::endl;
//   pmtInfo.baseline = subset_mean;
//   pmtInfo.baselineSigma = subset_stddev;
// }

// double sbnd::trigger::pmtSoftwareTriggerProducer::estimateBaseline(std::vector<uint32_t> wvfm){
//   std::vector<uint32_t> subset = std::vector<uint32_t>(wvfm.begin(), wvfm.begin()+int(fBaselineWindow*us_to_ticks));
//   double subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size());
//   double val = 0;
//   for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
//   double subset_stddev = sqrt(val/subset.size());

//   if (subset_stddev > 3){ // make this fcl parameter?
//     val = 0; subset.clear(); subset_stddev = 0;
//     subset = std::vector<uint32_t>(wvfm.end()-int(fBaselineWindow*us_to_ticks), wvfm.end());
//     subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size());
//     for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
//     subset_stddev = sqrt(val/subset.size());
//   }
//   return subset_mean;
// }

double sbnd::trigger::pmtSoftwareTriggerProducer::estimateBaseline(std::vector<uint32_t> wvfm){
    const auto median_it = wvfm.begin() + wvfm.size() / 2;
    std::nth_element(wvfm.begin(), median_it , wvfm.end());
    auto median = *median_it;
    return median;
}

double sbnd::trigger::pmtSoftwareTriggerProducer::estimateBaseline(std::vector<uint16_t> wvfm){
    const auto median_it = wvfm.begin() + wvfm.size() / 2;
    std::nth_element(wvfm.begin(), median_it , wvfm.end());
    auto median = *median_it;
    return median;
}

/*
PE threshold algorithm
*/
void sbnd::trigger::pmtSoftwareTriggerProducer::SimpleThreshAlgo(int i_ch){
  auto wvfm = fWvfmsVec[i_ch];
  auto &pmtInfo = fpmtInfoVec[i_ch]; 
  double baseline = pmtInfo.baseline;

  bool fire = false; // bool for if pulse has been detected
  int counter = 0; // counts the bin of the waveform

  // these should be fcl parameters 
  double start_adc_thres = 5, end_adc_thres = 2; 
  auto start_threshold = baseline-start_adc_thres;
  auto end_threshold   = baseline-end_adc_thres; 

  std::vector<sbnd::trigger::pmtPulse> pulse_vec;
  sbnd::trigger::pmtPulse pulse; 
  pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
  for (auto const &adc : wvfm){
    if ( !fire && ((double)adc) <= start_threshold ){ // if its a new pulse 
      fire = true;
      pulse.t_start = counter - 1 > 0 ? counter - 1 : counter;    
    }

    else if( fire && ((double)adc) > end_threshold ){ // found end of a pulse
      fire = false;
      pulse.t_end = counter < ((int)wvfm.size())  ? counter : counter - 1;
      pulse_vec.push_back(pulse);
      pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
    }   

    else if(fire){ // if we're in a pulse 
      pulse.area += (baseline-(double)adc);
      if ((baseline-(double)adc) > pulse.peak) { // Found a new maximum
        pulse.peak = (baseline-(double)adc);
        pulse.t_peak = counter;
      }
    }
    counter++;
  }

  if(fire){ // Take care of a pulse that did not finish within the readout window.
    fire = false;
    pulse.t_end = counter - 1;
    pulse_vec.push_back(pulse);
    pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
  }

  pmtInfo.pulseVec = pulse_vec;
  // calculate PE from area 
  for (auto &pulse : pmtInfo.pulseVec){pulse.pe = pulse.area/fPEArea;}
}

DEFINE_ART_MODULE(sbnd::trigger::pmtSoftwareTriggerProducer)
