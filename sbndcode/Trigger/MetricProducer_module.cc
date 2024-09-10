////////////////////////////////////////////////////////////////////////
// Class:       MetricProducer
// Plugin Type: producer (Unknown Unknown)
// File:        MetricProducer_module.cc
//
// Michelle Stancari
// August 2022
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

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
// #include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"
// #include "sbnobj/SBND/Trigger/CRTmetric.hh"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/CRTmetric.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <iostream>
#include <algorithm>
#include <valarray>
#include <numeric>

namespace sbndaq {
  class MetricProducer;
}

class sbndaq::MetricProducer : public art::EDProducer {
public:
  explicit MetricProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MetricProducer(MetricProducer const&) = delete;
  MetricProducer(MetricProducer&&) = delete;
  MetricProducer& operator=(MetricProducer const&) = delete;
  MetricProducer& operator=(MetricProducer&&) = delete;

  // Required functions.
  void produce(art::Event& evt) override;

private:

  //CRT Metric variables

  // fhicl parameters
  art::Persistable is_persistable_;
  int fBeamWindowStart;
  int fBeamWindowEnd;
  bool fVerbose;
  bool fCalcCRTMetrics;

  //metric variables
  int hitsperplane[7];

  //PMT Metric variables

  //fhicl parameters
  bool fCalcPMTMetrics;
  double fTriggerTimeOffset;    // offset of trigger time, default 0.5 sec
  double fBeamWindowLength; // beam window length after trigger time, default 1.6us
  uint32_t fWvfmLength;

  bool fCalculateBaseline;
  bool fCountPMTs;
  bool fCalculatePEMetrics;
  bool fFindPulses;

  std::vector<double> fInputBaseline;
  int fADCThreshold;
  double fPEArea; // conversion factor from ADCxns area to PE count 

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

  //both info
  int num_crt_frags;
  int num_pmt_frags;


  void analyze_crt_fragment(artdaq::Fragment & frag);
  void checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag);
  void analyzeCAEN1730Fragment(const artdaq::Fragment &frag);
  void estimateBaseline(int i_ch);
  void SimpleThreshAlgo(int i_ch);

  // variables to create an output root tree
  TTree* _tree;
  int _run, _sub, _evt;

  // crt tree variables
  int _crt_hitsperplane[7];
  
  // pmt tree variables 
  bool   _pmt_beam_trig; 
  double _pmt_time_trig; 
  int    _pmt_npmt;
  double _pmt_promptPE, _pmt_prelimPE; 

};


sbndaq::MetricProducer::MetricProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  is_persistable_(p.get<bool>("is_persistable", true) ? art::Persistable::Yes : art::Persistable::No),
  fBeamWindowStart(p.get<int>("BeamWindowStart",320000)),
  fBeamWindowEnd(p.get<int>("BeamWindowEnd",350000)),
  fVerbose(p.get<bool>("Verbose",false)),
  fCalcCRTMetrics(p.get<bool>("CalcCRTMetrics",true)),
  fCalcPMTMetrics(p.get<bool>("CalcPMTMetrics",true)),
  fTriggerTimeOffset(p.get<double>("TriggerTimeOffset", 0.5)),
  fBeamWindowLength(p.get<double>("BeamWindowLength", 1.6)),
  fWvfmLength(p.get<uint32_t>("WvfmLength", 5120)),
  fCalculateBaseline(p.get<bool>("CalculateBaseline",true)),
  fCountPMTs(p.get<bool>("CountPMTs",true)),
  fCalculatePEMetrics(p.get<bool>("CalculatePEMetrics",false)),
  fFindPulses(p.get<bool>("FindPulses", false)),
  fADCThreshold(p.get<double>("ADCThreshold", 7960)),
  fPEArea(p.get<double>("PEArea", 66.33))
  {
  // Call appropriate produces<>() functions here.
  produces< std::vector<sbndaq::CRTmetric >>("", is_persistable_);

  produces<std::vector<sbnd::trigger::pmtSoftwareTrigger>>("", is_persistable_);

  beamWindowStart = fTriggerTimeOffset*1e9;
  beamWindowEnd = beamWindowStart + fBeamWindowLength*1000;

  // build PD map and channel list
  auto subsetCondition = [](auto const& i)->bool { return i["pd_type"] == "pmt_coated" || i["pd_type"] == "pmt_uncoated"; };
  auto pmtMap = pdMap.getCollectionFromCondition(subsetCondition);
  for(auto const& i:pmtMap){
    channelList.push_back(i["channel"]);
  }

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("softmetrictree","");
  _tree->Branch("run",       &_run,       "run/I");
  _tree->Branch("sub",       &_sub,       "sub/I");
  _tree->Branch("evt",       &_evt,       "evt/I");
  _tree->Branch("crt_hitsperplane", &_crt_hitsperplane, "crt_hitsperplane[7]/I");
  _tree->Branch("pmt_beam_trig", &_pmt_beam_trig, "pmt_beam_trig/O");
  _tree->Branch("pmt_time_trig", &_pmt_time_trig, "pmt_time_trig/D");
  _tree->Branch("pmt_npmt",      &_pmt_npmt,      "pmt_npmt/I");
  _tree->Branch("pmt_promptPE",  &_pmt_promptPE,  "pmt_promptPE/D");
  _tree->Branch("pmt_prelimPE",  &_pmt_prelimPE,  "pmt_prelimPE/D");
}


void sbndaq::MetricProducer::produce(art::Event& evt)
{
  _run = evt.run();
  _sub = evt.subRun();
  _evt = evt.id().event();

  if (fVerbose) std::cout << "Processing: Run " << _run << ", Subrun " << _sub << ", Event " << _evt << std::endl;

  // object to store required trigger information in
  std::unique_ptr<std::vector<sbndaq::CRTmetric>> crt_metrics_v = std::make_unique<std::vector<sbndaq::CRTmetric>>();
  sbndaq::CRTmetric crt_metrics; 
  // std::unique_ptr<sbndaq::CRTmetric> CRTMetricInfo = std::make_unique<sbndaq::CRTmetric>();

  // object to store trigger metrics in
  std::unique_ptr<std::vector<sbnd::trigger::pmtSoftwareTrigger>> pmt_metrics_v = std::make_unique<std::vector<sbnd::trigger::pmtSoftwareTrigger>>();
  sbnd::trigger::pmtSoftwareTrigger pmt_metrics;
  // clear variables at the beginning of the event
  // move this to constructor??
  for (int ip=0;ip<7;++ip)  { crt_metrics.hitsperplane[ip]=0; hitsperplane[ip]=0;}
  foundBeamTrigger = false;
  fWvfmsFound = false;
  fWvfmsVec.clear(); fWvfmsVec.resize(15*8); // 15 pmt channels per fragment, 8 fragments per trigger
  fpmtInfoVec.clear(); fpmtInfoVec.resize(15*8);

  _pmt_beam_trig = false;
  _pmt_time_trig = -9999;
  _pmt_npmt = -9999;
  _pmt_promptPE = -9999; _pmt_prelimPE = -9999;

  std::fill(_crt_hitsperplane, _crt_hitsperplane+7, 0);

  // get fragment handles
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();

  num_crt_frags = 0;
  num_pmt_frags = 0;
  // loop over fragment handles
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      // container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() == sbndaq::detail::FragmentType::BERNCRTV2){
          if (fVerbose)     std::cout << "    Found " << contf.block_count() << " CRT Fragments in container " << std::endl;
          if (fCalcCRTMetrics){
              for (size_t ii = 0; ii < contf.block_count(); ++ii) analyze_crt_fragment(*contf[ii].get());
          }
        }
      }
    }
    else {
      // normal fragment
      size_t beamFragmentIdx = -1;
      for (auto frag : *handle){
        beamFragmentIdx++;
        if (frag.type()==sbndaq::detail::FragmentType::BERNCRTV2) {
          num_crt_frags++;
          if (fCalcCRTMetrics){analyze_crt_fragment(frag);}
        }

        if (frag.type()==sbndaq::detail::FragmentType::CAENV1730) {
          num_pmt_frags++;
          if (fCalcPMTMetrics){
          // identify whether any fragments correspond to the beam spill
          // loop over fragments, in steps of 8
          //size_t beamFragmentIdx = 9999;
          //for (size_t fragmentIdx = 0; fragmentIdx < handle->size(); fragmentIdx += 8) {
            if (!foundBeamTrigger){
              checkCAEN1730FragmentTimeStamp(frag);//handle->at(fragmentIdx));
              // if set of fragment in time with beam found, process waveforms
              if (foundBeamTrigger && beamFragmentIdx != 9999) {
                for (size_t fragmentIdx = beamFragmentIdx; fragmentIdx < beamFragmentIdx+8; fragmentIdx++) {
                  analyzeCAEN1730Fragment(handle->at(fragmentIdx));
                }
                fWvfmsFound = true;
              }
            }
          } //if saving pmt metrics
        } //if is pmt frag
      } //loop over frags
    }
  } // loop over frag handles

  if (fVerbose)   std::cout << "Found " << num_crt_frags << " BERNCRTV2 fragments" << std::endl;
  if (fVerbose)   std::cout << "Found " << num_pmt_frags << " CAEN1730 fragments" << std::endl;

  if (fCalcCRTMetrics){

    for (int i=0;i<7;++i) {
      crt_metrics.hitsperplane[i] = hitsperplane[i];
      _crt_hitsperplane[i] = hitsperplane[i];
      }

    if (fVerbose) {
      std::cout << "CRT hit count during beam spill ";
      for (int i=0;i<7;++i) std::cout << i << ": " << crt_metrics.hitsperplane[i] << "   " ;
      std::cout << std::endl;
    }


  }//if save crt metrics

  if (fCalcPMTMetrics){

    if (foundBeamTrigger && fWvfmsFound) {
      // store timestamp of trigger, relative to beam window start
      double triggerTimeStamp = fTriggerTime - beamWindowStart;
      pmt_metrics.foundBeamTrigger = true;
      pmt_metrics.trig_ts = triggerTimeStamp;

      _pmt_beam_trig = true;
      _pmt_time_trig = triggerTimeStamp;
      
      if (fVerbose) std::cout << "Saving PMT trigger timestamp: " << triggerTimeStamp << " ns" << std::endl;

      double promptPE = 0;
      double prelimPE = 0;

      int nAboveThreshold = 0;
      // find the waveform bins that correspond to the start and end of the extended spill window (0 -> 1.8 us) within the 10 us waveform
      // if the triggerTimeStamp < 1000, the beginning of the beam spill is *not* contained within the waveform
      int beamStartBin = (triggerTimeStamp >= 1000)? 0 : int(500 - abs((triggerTimeStamp-1000)/2));
      int beamEndBin   = (triggerTimeStamp >= 1000)? int(500 + (fBeamWindowLength*1e3 - triggerTimeStamp)/2) : (beamStartBin + (fBeamWindowLength*1e3)/2);

      // wvfm loop to calculate metrics
      for (int i_ch = 0; i_ch < 120; ++i_ch){
        auto &pmtInfo = fpmtInfoVec.at(i_ch);
        auto wvfm = fWvfmsVec[i_ch];

        // assign channel
        pmtInfo.channel = channelList.at(i_ch);

        // calculate baseline
        if (fCalculateBaseline) estimateBaseline(i_ch);
        else { pmtInfo.baseline=fInputBaseline.at(0); pmtInfo.baselineSigma = fInputBaseline.at(1); }

        // count number of PMTs above threshold
        if (fCountPMTs){
          for (int bin = beamStartBin; bin < beamEndBin; ++bin){
            auto adc = wvfm[bin];
            if (adc < fADCThreshold){ nAboveThreshold++; break; } 
          }
        }
        else nAboveThreshold=-9999;

        // quick estimate prompt and preliminary light, assuming sampling rate of 500 MHz (2 ns per bin)
        if (fCalculatePEMetrics){
          double baseline = pmtInfo.baseline;
          auto prompt_window = std::vector<uint16_t>(wvfm.begin()+500, wvfm.begin()+1000);
          auto prelim_window = std::vector<uint16_t>(wvfm.begin()+beamStartBin, wvfm.begin()+500);
          if (fFindPulses == false){
            double ch_promptPE = (baseline-(*std::min_element(prompt_window.begin(), prompt_window.end())))/8;
            double ch_prelimPE = (baseline-(*std::min_element(prelim_window.begin(), prelim_window.end())))/8;
            promptPE += ch_promptPE;
            prelimPE += ch_prelimPE;
          }
          // pulse finder + prompt and prelim calculation with pulses
          if (fFindPulses == true){
            SimpleThreshAlgo(i_ch);
            for (auto pulse : pmtInfo.pulseVec){
              if (pulse.t_start > 500 && pulse.t_end < 550) promptPE+=pulse.pe;
              if ((triggerTimeStamp) >= 1000){ if (pulse.t_end < 500) prelimPE+=pulse.pe; }
              else if (triggerTimeStamp < 1000){
                if (pulse.t_start > (500 - abs((triggerTimeStamp-1000)/2)) && pulse.t_end < 500) prelimPE+=pulse.pe;
              }
            }
          }
        } 
      } // end of wvfm loop

      pmt_metrics.nAboveThreshold = nAboveThreshold;    
      pmt_metrics.promptPE = promptPE;
      pmt_metrics.prelimPE = prelimPE;

      // tree variables 
      _pmt_npmt = nAboveThreshold;
      _pmt_promptPE = promptPE; 
      _pmt_prelimPE = prelimPE;

      if (fVerbose && fCountPMTs) std::cout << "nPMTs Above Threshold: " << nAboveThreshold << std::endl;
      if (fVerbose && fCalculatePEMetrics) std::cout << "prompt pe: " << promptPE << std::endl;
      if (fVerbose && fCalculatePEMetrics) std::cout << "prelim pe: " << prelimPE << std::endl;
    } // if beam and wvfms found 
    else{
      if (fVerbose) std::cout << "Beam and wvfms not found" << std::endl;
      pmt_metrics.foundBeamTrigger = false;
      pmt_metrics.trig_ts = -9999;
      pmt_metrics.nAboveThreshold = -9999;
      pmt_metrics.promptPE = -9999;
      pmt_metrics.prelimPE = -9999;

      // tree variables 
      _pmt_beam_trig = false; 
      _pmt_time_trig = -9999; _pmt_npmt = -9999; _pmt_promptPE = -9999; _pmt_prelimPE = -9999;
    }
  }//if save pmt metrics

  // add to event
  crt_metrics_v->push_back(crt_metrics);
  evt.put(std::move(crt_metrics_v));
  pmt_metrics_v->push_back(pmt_metrics);
  evt.put(std::move(pmt_metrics_v)); 
  _tree->Fill(); 
}//produce



void sbndaq::MetricProducer::analyze_crt_fragment(artdaq::Fragment & frag)
{

  sbndaq::BernCRTFragmentV2 bern_fragment(frag);

  // use  fragment ID to get plane information
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  //frag.sequenceID()
  auto thisone = frag.fragmentID();  uint plane = (thisone & 0x0700) >> 8;
  if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}

  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
    // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
    auto thisflag = bevt->flags;
    if (thisflag & 0x2 && !(thisflag & 0xC) ) {
      // check ts1 for beam window
      auto thistime=bevt->ts1;
      if ((int)thistime>fBeamWindowStart && (int)thistime<fBeamWindowEnd) hitsperplane[plane]++;
      //crt_metrics.hitsperplane[plane]++;
    }
  }



}//analyze crt fragments


void sbndaq::MetricProducer::checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag) {

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
}//check caen 1730 timestamp

void sbndaq::MetricProducer::analyzeCAEN1730Fragment(const artdaq::Fragment &frag) {

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
}//analyze caen 1730 fragment

void sbndaq::MetricProducer::estimateBaseline(int i_ch){
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
}//estimateBaseline

void sbndaq::MetricProducer::SimpleThreshAlgo(int i_ch){
  auto wvfm = fWvfmsVec[i_ch];
  auto &pmtInfo = fpmtInfoVec[i_ch];
  double baseline = pmtInfo.baseline;
  double baseline_sigma = pmtInfo.baselineSigma;

  bool fire = false; // bool for if pulse has been detected
  int counter = 0; // counts the bin of the waveform

  // these should be fcl parameters
  double start_adc_thres = 5, end_adc_thres = 2;
  double nsigma_start = 5, nsigma_end = 3;

  auto start_threshold = ( start_adc_thres > (nsigma_start * baseline_sigma) ? (baseline-start_adc_thres) : (baseline-(nsigma_start * baseline_sigma)));
  auto end_threshold = ( end_adc_thres > (nsigma_end * baseline_sigma) ? (baseline - end_adc_thres) : (baseline - (nsigma_end * baseline_sigma)));

  std::vector<sbnd::trigger::pmtPulse> pulse_vec;
  sbnd::trigger::pmtPulse pulse;
  pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
  for (auto const &adc : wvfm){
    if ( !fire && ((double)adc) <= start_threshold ){ // if its a new pulse
      fire = true;
      //vic: i move t_start back one, this helps with porch
      pulse.t_start = counter - 1 > 0 ? counter - 1 : counter;
    }

    if( fire && ((double)adc) > end_threshold ){ // found end of a pulse
      fire = false;
      //vic: i move t_start forward one, this helps with tail
      pulse.t_end = counter < ((int)wvfm.size())  ? counter : counter - 1;
      pulse_vec.push_back(pulse);
      pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
    }

    if(fire){ // if we're in a pulse
      pulse.area += (baseline-(double)adc);
      if (pulse.peak > (baseline-(double)adc)) { // Found a new maximum
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
}//SimpleThreshAlgo


// -------------------------------------------------

// -------------------------------------------------

DEFINE_ART_MODULE(sbndaq::MetricProducer)
