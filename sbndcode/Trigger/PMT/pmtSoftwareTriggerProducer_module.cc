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
// #include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"

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
  double fTriggerTimeOffset;    // offset of trigger time, default 0.5 sec
  double fBeamWindowLength; // beam window length after trigger time, default 1.6us
  uint32_t fWvfmLength;
  bool fVerbose;
  bool fSaveHists;

  bool fCalculateBaseline;
  bool fCountPMTs;
  bool fCalculatePEMetrics;
  bool fFindPulses;

  std::vector<double> fInputBaseline;
  int fADCThreshold;
  double fPEArea; // conversion factor from ADCxns area to PE count 

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
  void estimateBaseline(int i_ch);
  void SimpleThreshAlgo(int i_ch);

  TTree* _tree; 
  int _run, _sub, _evt; 
  bool   _beam_trig; 
  double _time_trig; 
  int    _npmt;
  double _promptPE, _prelimPE;
  std::vector<double> _ch_prelimPE, _ch_promptPE;
  std::vector<int> _ch_AboveThreshold,_ch_type,_ch_ID;


  TTree* _pulse_tree;
  int _npulses; 
  std::vector<int> _ch_npulses; // number of pulses per channel
  
  std::vector<int> _pulse_ch; // ch number for each pulse 
  std::vector<double> _pulse_t_start; // t_start for each pulse 
  std::vector<double> _pulse_t_end;
  std::vector<double> _pulse_t_peak;
  std::vector<double> _pulse_peak;
  std::vector<double> _pulse_area;
};


sbnd::trigger::pmtSoftwareTriggerProducer::pmtSoftwareTriggerProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  // is_persistable_(p.get<bool>("is_persistable", true) ? art::Persistable::Yes : art::Persistable::No),
  fTriggerTimeOffset(p.get<double>("TriggerTimeOffset", 0.5)),
  fBeamWindowLength(p.get<double>("BeamWindowLength", 1.6)), 
  fWvfmLength(p.get<uint32_t>("WvfmLength", 5120)),
  fVerbose(p.get<bool>("Verbose", false)),
  fSaveHists(p.get<bool>("SaveHists",false)),
  fCalculateBaseline(p.get<bool>("CalculateBaseline",true)),
  fCountPMTs(p.get<bool>("CountPMTs",true)),
  fCalculatePEMetrics(p.get<bool>("CalculatePEMetrics",false)),
  fFindPulses(p.get<bool>("FindPulses", false)),
  fInputBaseline(p.get<std::vector<double>>("InputBaseline")),
  fADCThreshold(p.get<double>("ADCThreshold", 7960)),
  fPEArea(p.get<double>("PEArea", 66.33))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // produces< sbnd::trigger::pmtSoftwareTrigger >("", is_persistsable_);
  produces<std::vector<sbnd::trigger::pmtSoftwareTrigger>>();

  beamWindowStart = fTriggerTimeOffset*1e9;
  beamWindowEnd = beamWindowStart + fBeamWindowLength*1000;

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
  _tree->Branch("promptPE",  &_promptPE,  "promptPE/D");
  _tree->Branch("prelimPE",  &_prelimPE,  "prelimPE/D");
  _tree->Branch("ch_prelimPE","std::vector<double>",&_ch_prelimPE);
  _tree->Branch("ch_promptPE","std::vector<double>",&_ch_promptPE);
  _tree->Branch("ch_AboveThreshold","std::vector<int>",&_ch_AboveThreshold);
  _tree->Branch("ch_ID","std::vector<int>",&_ch_ID);
  _pulse_tree = fs->make<TTree>("pulse_tree","");
  _pulse_tree->Branch("run",       &_run,       "run/I");
  _pulse_tree->Branch("sub",       &_sub,       "sub/I");
  _pulse_tree->Branch("evt",       &_evt,       "evt/I");
  _pulse_tree->Branch("ch_npulses", "std::vector<int>",&_ch_npulses); 
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

  if (fVerbose) std::cout << "Processing Run: " << fRun << ", Subrun: " << fSubrun << ", Event: " << fEvent << std::endl;

  // reset for this event
  foundBeamTrigger = false;
  fWvfmsFound = false;
  fWvfmsVec.clear(); fWvfmsVec.resize(15*8); // 15 pmt channels per fragment, 8 fragments per trigger
  fpmtInfoVec.clear(); fpmtInfoVec.resize(15*8); 

  _beam_trig = false;
  _time_trig = -9999;
  _npmt = -9999;
  _promptPE = -9999; _prelimPE = -9999;
  //Initialize size, make it fast - to be pushed into the tree vector
  std::vector<double> ch_promptPE(120,-9999), ch_prelimPE(120,-9999);
  std::vector<int> ch_AboveThreshold(120,-9999),ch_ID(120,-9999);
  // _TREE_VECTOR EVENT_VECTOR

  // get fragment handles
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  // loop over fragment handles
  for (auto &handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;
    if (handle->front().type()==sbndaq::detail::FragmentType::CAENV1730) {
      if (fVerbose)   std::cout << "Found " << handle->size() << " CAEN1730 fragments" << std::endl;
      
      // identify whether any fragments correspond to the beam spill
      // loop over fragments, in steps of 8
      size_t beamFragmentIdx = 9999;
      for (size_t fragmentIdx = 0; fragmentIdx < handle->size(); fragmentIdx += 8) {
        checkCAEN1730FragmentTimeStamp(handle->at(fragmentIdx));
        if (foundBeamTrigger) {
          beamFragmentIdx = fragmentIdx;
          if (fVerbose) std::cout << "Found fragment in time with beam at index: " << beamFragmentIdx << std::endl;
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
    double triggerTimeStamp = fTriggerTime - beamWindowStart;
    trig_metrics.trig_ts = triggerTimeStamp;
    _time_trig = triggerTimeStamp;
    if (fVerbose) std::cout << "Saving trigger timestamp: " << trig_metrics.trig_ts << " ns" << std::endl;

    double promptPE = 0;
    double prelimPE = 0; 

    int nAboveThreshold = 0;
    // find the waveform bins that correspond to the start and end of the extended spill window (0 -> 1.8 us) within the 10 us waveform 
    // !! if the triggerTimeStamp > 1000, the beginning of the beam spill is *not* contained within the waveform 
    int beamStartBin = (triggerTimeStamp >= 1000)? 0 : int(500 - (triggerTimeStamp)/2); // units of bins 
    int beamEndBin   = (triggerTimeStamp >= 1000)? (fBeamWindowLength*1e3 - (triggerTimeStamp-1000) )/2 : (beamStartBin + (fBeamWindowLength*1e3)/2);
    
    std::cout << "beamStartBin: " << beamStartBin << std::endl;
    std::cout << "beamEndBin: "   << beamEndBin << std::endl;
    // wvfm loop to calculate metrics 

    _npulses = 0;
    _pulse_ch.clear(); //_pulse_npulses.clear();
    _pulse_t_start.clear(); _pulse_t_end.clear(); _pulse_t_peak.clear();
    _pulse_peak.clear(); _pulse_area.clear();

    _pulse_ch.reserve(1000); //_pulse_npulses.reserve(1000);
    _pulse_t_start.reserve(1000); _pulse_t_end.reserve(1000); _pulse_t_peak.reserve(1000);
    _pulse_peak.reserve(1000); _pulse_area.reserve(1000);

    for (int i_ch = 0; i_ch < 120; ++i_ch){
      ch_ID[i_ch] = channelList.at(i_ch);
      auto &pmtInfo = fpmtInfoVec.at(i_ch);
      auto wvfm = fWvfmsVec[i_ch];

      // assign channel 
      pmtInfo.channel = channelList.at(i_ch);

         // calculate baseline 
      if (fCalculateBaseline) estimateBaseline(i_ch);
      else { pmtInfo.baseline=fInputBaseline.at(0); pmtInfo.baselineSigma = fInputBaseline.at(1); }

      // count number of PMTs above threshold within the beam window
      if (fCountPMTs){
        for (int bin = beamStartBin; bin < beamEndBin; ++bin){
          auto adc = wvfm[bin];
          if (adc < fADCThreshold){ 
            ch_AboveThreshold[i_ch] = 1;
            nAboveThreshold++; 
            continue; 
          } 
          else{
              ch_AboveThreshold[i_ch] = 0;
            continue;
          }
        }
      }
      else {nAboveThreshold=-9999;ch_AboveThreshold[i_ch] = -9999;}

      // quick estimate prompt and preliminary light, assuming sampling rate of 500 MHz (2 ns per bin)
      if (fCalculatePEMetrics){
        double baseline = pmtInfo.baseline;
        auto prompt_window = std::vector<uint16_t>(wvfm.begin()+500, wvfm.begin()+1000);
        auto prelim_window = std::vector<uint16_t>(wvfm.begin()+beamStartBin, wvfm.begin()+500);
        if (fFindPulses == false){
          double ch_promptPE_ = (baseline-(*std::min_element(prompt_window.begin(), prompt_window.end())))/8;
          double ch_prelimPE_ = (baseline-(*std::min_element(prelim_window.begin(), prelim_window.end())))/8;
          ch_prelimPE[i_ch] = ch_prelimPE_;
          ch_promptPE[i_ch] = ch_promptPE_;
          promptPE += ch_promptPE_;
          prelimPE += ch_prelimPE_;
        }
      }
      else {promptPE = -9999; prelimPE =-9999;ch_prelimPE[i_ch] = -9999; ch_promptPE[i_ch] = -9999;}

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
          if (pulse.t_start > 500 && pulse.t_end < 550){
             ch_promptPE[i_ch] += pulse.pe;
             promptPE+=pulse.pe;
          }
          if ((triggerTimeStamp) >= 1000){ 
            if (pulse.t_end < 500) {
              ch_prelimPE[i_ch] += pulse.pe;
              prelimPE+=pulse.pe;
            } 
          }
          else if (triggerTimeStamp < 1000){
            if (pulse.t_start > (500 - abs((triggerTimeStamp-1000)/2)) && pulse.t_end < 500) {
              ch_prelimPE[i_ch] += pulse.pe;
              prelimPE+=pulse.pe;
            } 
          }
          _npulses++;
          _pulse_ch.push_back(channelList.at(i_ch));
          _pulse_t_start.push_back((pulse.t_start - beamStartBin)*2);
          _pulse_t_end.push_back(  (pulse.t_end   - beamStartBin)*2);
          _pulse_t_peak.push_back( (pulse.t_peak  - beamStartBin)*2);
          _pulse_peak.push_back(pulse.peak);
          _pulse_area.push_back(pulse.area);
        } // end of pulse loop 
      }
    } // end of wvfm loop 

    _pulse_tree->Fill();

    //We need to fill this - data product which is a LArSoft class
    trig_metrics.nAboveThreshold = nAboveThreshold;    
    trig_metrics.promptPE = promptPE;
    trig_metrics.prelimPE = prelimPE;

    // tree variables - TTree which is only part of root 
    _npmt = nAboveThreshold;
    _promptPE = promptPE; 
    _prelimPE = prelimPE;

    if (fVerbose && fCountPMTs) std::cout << "nPMTs Above Threshold: " << trig_metrics.nAboveThreshold << std::endl;
    if (fVerbose && fCalculatePEMetrics) std::cout << "prompt pe: " << trig_metrics.promptPE << std::endl;
    if (fVerbose && fCalculatePEMetrics) std::cout << "prelim pe: " << trig_metrics.prelimPE << std::endl;

    // start histo 
    if (fSaveHists == true){
      int hist_id = -1; 
      for (size_t i_wvfm = 0; i_wvfm < fWvfmsVec.size(); ++i_wvfm){
        std::vector<uint16_t> wvfm = fWvfmsVec[i_wvfm];
        hist_id++;
        //if (fEvent<4){
            histname.str(std::string());
            histname << "run_" << fRun  
                    << "subrun_" <<fSubrun
                    << "event_" << fEvent
                    << "_pmtnum_" << channelList.at(i_wvfm);
            // assuming that we save ~1 us before the triggerTimeStamp  
            double StartTime = (triggerTimeStamp-1000)*1e-3; // us
            double EndTime   = StartTime + (fWvfmLength*2)*1e-3;

            TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), "Raw Waveform", wvfm.size(), StartTime, EndTime);
            wvfmHist->GetXaxis()->SetTitle("t (#mus)");
            for(unsigned int i = 0; i < wvfm.size(); i++) {
              wvfmHist->SetBinContent(i + 1, (double)wvfm[i]);
            }
        //} 
      } // end histo
    }
  }
  else{
    if (fVerbose) std::cout << "Beam and wvfms not found" << std::endl;
    trig_metrics.foundBeamTrigger = false;
    trig_metrics.trig_ts = -9999;
    trig_metrics.nAboveThreshold = -9999;
    trig_metrics.promptPE = -9999;
    trig_metrics.prelimPE = -9999;
    // tree variables 
    _beam_trig = false; 
    _time_trig = -9999; _npmt = -9999; _promptPE = -9999; _prelimPE = -9999;
    ch_prelimPE = std::vector<double>(120,-9999);
    ch_promptPE = std::vector<double>(120,-9999);
    ch_AboveThreshold = std::vector<int>(120,-9999);
  }

      _ch_prelimPE = (ch_prelimPE);
      _ch_promptPE = (ch_promptPE);
      _ch_ID = (ch_ID);
      _ch_AboveThreshold = (ch_AboveThreshold);
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

void sbnd::trigger::pmtSoftwareTriggerProducer::estimateBaseline(int i_ch){
  auto wvfm = fWvfmsVec[i_ch];
  auto &pmtInfo = fpmtInfoVec[i_ch]; 
  // assuming that the first 500 ns doesn't include peaks, find the mean of the ADC count as the baseline 
  std::vector<uint16_t> subset = std::vector<uint16_t>(wvfm.begin(), wvfm.begin()+250);
  double subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size()); 
  double val = 0;
  for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
  double subset_stddev = sqrt(val/subset.size()); 

  // if the first 500 ns seem to be messy, use the last 500 
  if (subset_stddev > 3){ // make this fcl parameter? 
    val = 0; subset.clear(); subset_stddev = 0;
    subset = std::vector<uint16_t>(wvfm.end()-500, wvfm.end());
    subset_mean = (std::accumulate(subset.begin(), subset.end(), 0))/(subset.size());
    for (size_t i=0; i<subset.size();i++){ val += (subset[i] - subset_mean)*(subset[i] - subset_mean);}
    subset_stddev = sqrt(val/subset.size()); 
  }
  pmtInfo.baseline = subset_mean;
  pmtInfo.baselineSigma = subset_stddev;
}

/*
PE threshold algorithm
*/
void sbnd::trigger::pmtSoftwareTriggerProducer::SimpleThreshAlgo(int i_ch){
  auto wvfm = fWvfmsVec[i_ch];
  auto &pmtInfo = fpmtInfoVec[i_ch]; 
  double baseline = pmtInfo.baseline;
  // double baseline_sigma = pmtInfo.baselineSigma;
  
  bool fire = false; // bool for if pulse has been detected
  int counter = 0; // counts the bin of the waveform

  // these should be fcl parameters 
  double start_adc_thres = 5, end_adc_thres = 2; 
  // double nsigma_start = 5, nsigma_end = 3; 
  
  // auto start_threshold = ( start_adc_thres > (nsigma_start * baseline_sigma) ? (baseline-start_adc_thres) : (baseline-(nsigma_start * baseline_sigma)));
  // auto end_threshold = ( end_adc_thres > (nsigma_end * baseline_sigma) ? (baseline - end_adc_thres) : (baseline - (nsigma_end * baseline_sigma)));
  baseline = 8000;
  auto start_threshold = baseline-start_adc_thres;
  auto end_threshold   = baseline-end_adc_thres; 

  std::vector<sbnd::trigger::pmtPulse> pulse_vec;
  sbnd::trigger::pmtPulse pulse; 
  pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
  for (auto const &adc : wvfm){
    if ( !fire && ((double)adc) <= start_threshold ){ // if its a new pulse 
      fire = true;
      //vic: i move t_start back one, this helps with porch
      pulse.t_start = counter - 1 > 0 ? counter - 1 : counter;    
    }

    else if( fire && ((double)adc) > end_threshold ){ // found end of a pulse
      fire = false;
      //vic: i move t_start forward one, this helps with tail
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

// void sbnd::trigger:pmtSoftwareTriggerProducer::SlidingThreshAlgo(){
// }

DEFINE_ART_MODULE(sbnd::trigger::pmtSoftwareTriggerProducer)
