////////////////////////////////////////////////////////////////////////
// Class:       pmtSoftwareTriggerProducer
// Plugin Type: producer (Unknown Unknown)
// File:        pmtSoftwareTriggerProducer_module.cc
//
// Generated at Thu Feb 17 13:22:51 2022 by Patrick Green using cetskelgen
// from  version .
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

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"

#include <memory>
#include <algorithm>
#include <valarray>

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
  art::Persistable is_persistable_;
  double fTriggerTimeOffset;    // offset of trigger time, default 0.5 sec
  double fBeamWindowLength; // beam window length after trigger time, default 1.6us
  uint32_t fWvfmLength;
  bool fVerbose;

  std::string fBaselineAlgo;
  double fInputBaseline;

  bool fAreaToPE; // Use area to calculate number of PEs
  double fSPEArea; // If AreaToPE is true, this number is used as single PE area (in ADC counts)

  // event information
  int fRun, fSubrun;
  art::EventNumber_t fEvent;

  // beam window
  // set in artdaqFragment producer, in reality would be provided by event builder
  bool foundBeamTrigger;
  uint32_t beamWindowStart;
  uint32_t beamWindowEnd;

  // waveforms
  uint32_t fTriggerTime;
  bool fWvfmsFound;
  std::vector<std::vector<uint16_t>> fWvfmsVec;
  std::vector<double> fWvfmsBaseline;
  std::vector<double> fWvfmsBaselineSigma; 

  void checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag);
  void analyzeCAEN1730Fragment(const artdaq::Fragment &frag);
  void calculateBaseline();
  void SimpleThreshAlgo();

};


sbnd::trigger::pmtSoftwareTriggerProducer::pmtSoftwareTriggerProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  is_persistable_(p.get<bool>("is_persistable", true) ? art::Persistable::Yes : art::Persistable::No),
  fTriggerTimeOffset(p.get<double>("TriggerTimeOffset", 0.5)),
  fBeamWindowLength(p.get<double>("BeamWindowLength", 1.6)), 
  fWvfmLength(p.get<uint32_t>("WvfmLength", 5120)),
  fVerbose(p.get<bool>("Verbose", false)),
  fBaselineAlgo(p.get<std::string>("BaselineAlgo", "constant")),
  fInputBaseline(p.get<double>("InputBaseline", 8000)),
  fAreaToPE(p.get<bool>("AreaToPE", false)),
  fSPEArea(p.get<double>("SPEArea", 66.33))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< sbnd::trigger::pmtSoftwareTrigger >("", is_persistable_);

  beamWindowStart = fTriggerTimeOffset*1e9;
  beamWindowEnd = beamWindowStart + fBeamWindowLength*1000;
}

void sbnd::trigger::pmtSoftwareTriggerProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  // event information
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  if (fVerbose) std::cout << "Processing Run: " << fRun << ", Subrun: " << fSubrun << ", Event: " << fEvent << std::endl;

  // reset for this event
  foundBeamTrigger = false;
  fWvfmsFound = false;
  fWvfmsVec.clear(); fWvfmsVec.resize(15*8); // 15 pmt channels per fragment, 8 fragments per trigger
  fWvfmsBaseline.clear(); fWvfmsBaseline.resize(15*8);
  fWvfmsBaselineSigma.clear(); fWvfmsBaselineSigma.resize(15*8);


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
      // std::cout << "here 1" << std::endl;
      for (size_t fragmentIdx = 0; fragmentIdx < handle->size(); fragmentIdx += 8) {
        // std::cout << "here 2" << std::endl;
        checkCAEN1730FragmentTimeStamp(handle->at(fragmentIdx));
        // std::cout << "here 3" << std::endl;
        if (foundBeamTrigger) {
          // std::cout << "here 4" << std::endl;
          beamFragmentIdx = fragmentIdx;
          if (fVerbose) std::cout << "Found fragment in time with beam at index: " << beamFragmentIdx << std::endl;
          break;
          // std::cout << "here 5" << std::endl;
        }
      }
      // std::cout << "here 5.5" << std::endl;
      // if set of fragment in time with beam found, process waveforms
      if (foundBeamTrigger && beamFragmentIdx != 9999) {
        // std::cout << "here 6" << std::endl;
        for (size_t fragmentIdx = beamFragmentIdx; fragmentIdx < beamFragmentIdx+8; fragmentIdx++) {
          // std::cout << "here 7" << std::endl;
          analyzeCAEN1730Fragment(handle->at(fragmentIdx));
          // std::cout << "here 8" << std::endl;
        }
        // std::cout << "here 9" << std::endl;
        fWvfmsFound = true;
      }
    }
  } // end loop over handles

  // object to store trigger metrics in
  std::unique_ptr<sbnd::trigger::pmtSoftwareTrigger> pmtSoftwareTriggerMetrics = std::make_unique<sbnd::trigger::pmtSoftwareTrigger>();
<<<<<<< HEAD

  if (foundBeamTrigger && fWvfmsFound) {

    // calculate baseline 
    if (fBaselineAlgo == "constant") std::fill(fWvfmsBaseline.begin(), fWvfmsBaseline.end(), fInputBaseline);
    else if (fBaselineAlgo == "estimate") calculateBaseline();
=======

  if (foundBeamTrigger && fWvfmsFound) {
    pmtSoftwareTriggerMetrics->foundBeamTrigger = true;
>>>>>>> origin/feature/pgreen_artdaq_fragment_conversion

    // store timestamp of trigger, relative to beam window start
    pmtSoftwareTriggerMetrics->triggerTimestamp = fTriggerTime - beamWindowStart;
    if (fVerbose) std::cout << "Saving trigger timestamp: " << fTriggerTime - beamWindowStart << " ns" << std::endl;  
  }
<<<<<<< HEAD
  else{
    pmtSoftwareTriggerMetrics->triggerTimestamp = -9999;
    if (fVerbose) std::cout << "Beam and wvfms not found" << std::endl;

    // add to event
    e.put(std::move(pmtSoftwareTriggerMetrics));      
  }
=======
  else {
    pmtSoftwareTriggerMetrics->foundBeamTrigger = false;
  }

  // add to event
  e.put(std::move(pmtSoftwareTriggerMetrics));    
>>>>>>> origin/feature/pgreen_artdaq_fragment_conversion
}

void sbnd::trigger::pmtSoftwareTriggerProducer::checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag) {

  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp
  uint32_t timestamp = md->timeStampNSec;
  std::cout << "timestamp: " << timestamp << std::endl;
  std::cout << "beamWindowStart: " << beamWindowStart << std::endl;
  std::cout << "beamWindowEnd: " << beamWindowEnd << std::endl;

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

  std::cout << "value: " << value << std::endl;
  
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

void sbnd::trigger::pmtSoftwareTriggerProducer::calculateBaseline(){
  for (size_t i_wvfm = 0; i_wvfm < fWvfmsVec.size(); ++i_wvfm){
    std::vector<uint16_t> wvfm = fWvfmsVec[i_wvfm];
    // assuming that the first 500 ns doesn't include peaks, find the mean of the ADC count as the baseline 
    // this is also assuming the sampling rate of the waveform is 1 ns 
    std::vector<uint16_t> slice_vec = std::vector<uint16_t>(wvfm.begin(), wvfm.begin()+500);
    std::valarray<uint16_t> slice(slice_vec.data(), slice_vec.size());

    double slice_mean = (slice.sum())/(slice.size());
    double val = 0;
    for (size_t i=0; i<slice.size();i++){ val += (slice[i] - slice_mean)*(slice[i] - slice_mean);}
    double slice_stddev = sqrt(val/slice.size()); 
    // put in some if statement about stddev 
    fWvfmsBaseline[i_wvfm] = slice_mean;
    fWvfmsBaselineSigma[i_wvfm] = slice_stddev;
    if (fVerbose) std::cout << "baseline (ADC):" << slice_mean << ", stddev:" <<  slice_stddev << std::endl;
  }
}

void sbnd::trigger::pmtSoftwareTriggerProducer::SimpleThreshAlgo(){
  for (size_t i_wvfm = 0; i_wvfm < fWvfmsVec.size(); ++i_wvfm){
    std::vector<uint16_t> wvfm = fWvfmsVec[i_wvfm];
    double baseline = fWvfmsBaseline[i_wvfm];
    double baseline_sigma = fWvfmsBaselineSigma[i_wvfm];
    
    bool fire = false; // bool for if pulse has been detected
    int counter = 0; // counts the number of pulses 

    // these should be fcl parameters 
    double start_adc_thres = 5, end_adc_thres = 2; 
    double nsigma_start = 5, nsigma_end = 3; 
    
    auto start_threshold = ( start_adc_thres > (nsigma_start * baseline_sigma) ? start_adc_thres : (nsigma_start * baseline_sigma));
    auto end_threshold = ( end_adc_thres > (nsigma_end * baseline_sigma) ? end_adc_thres : (nsigma_end * baseline_sigma));
    
    start_threshold += baseline; end_threshold += baseline; 

    std::vector<sbnd::trigger::pmtPulse> pulse_vec;
    sbnd::trigger::pmtPulse pulse; 
    for (auto const &adc : wvfm){
      if ( !fire && ((double)adc) >= start_threshold ){ // if its a new pulse 
        fire = true;
        //vic: i move t_start back one, this helps with porch
        pulse.t_start = counter - 1 > 0 ? counter - 1 : counter;    
      }
  
      if( fire && ((double)adc) < end_threshold ){ // found end of a pulse
        fire = false;
        //vic: i move t_start forward one, this helps with tail
        pulse.t_end = counter < ((int)wvfm.size())  ? counter : counter - 1;
        pulse_vec.push_back(pulse);
        pulse.area = 0; pulse.peak = 0; pulse.t_start = 0; pulse.t_end = 0; pulse.t_peak = 0;
      }   

      if(fire){ // if we're in a pulse 
        if (fAreaToPE == true) // Add this adc count to the integral
          pulse.area += ((double)adc - baseline);
        if (pulse.peak < ((double)adc - baseline)) { // Found a new maximum
          pulse.peak = ((double)adc - baseline);
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

    // for (size_t i_pulse = 0; i_pulse < pulse_vec.size(); ++i_pulse){
    //   sbnd::trigger::pmtPulse &pulse = pulse_vec[i_pulse];
    //   if (fAreaToPE == true){ pulse.pe = pulse.area/fSPEArea;}
    // }
  }
}

// void sbnd::trigger:pmtSoftwareTriggerProducer::SlidingThreshAlgo(){
// }

DEFINE_ART_MODULE(sbnd::trigger::pmtSoftwareTriggerProducer)