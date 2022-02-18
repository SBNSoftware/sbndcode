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
  bool fVerbose;

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
  uint32_t fWvfmLength;
  bool fWvfmsFound;
  std::vector<std::vector<uint16_t>> fWvfmsVec;

  void checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag);
  void analyzeCAEN1730Fragment(const artdaq::Fragment &frag);

};


sbnd::trigger::pmtSoftwareTriggerProducer::pmtSoftwareTriggerProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  is_persistable_(p.get<bool>("is_persistable", true) ? art::Persistable::Yes : art::Persistable::No),
  fVerbose(p.get<bool>("Verbose", false))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< sbnd::trigger::pmtSoftwareTrigger >("", is_persistable_);

  // Call appropriate consumes<>() for any products to be retrieved by this module.

  beamWindowStart = 0.5*1e9; // 0.5 seconds // set as fhicl parameter
  beamWindowEnd = beamWindowStart + 1600;

  fWvfmLength = 5000; // 10 us, 2ns sampling // set as fhicl parameter (or could extract from fragment?)
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

  if (foundBeamTrigger && fWvfmsFound) {

    // object to store trigger metrics in
    std::unique_ptr<sbnd::trigger::pmtSoftwareTrigger> pmtSoftwareTriggerMetrics = std::make_unique<sbnd::trigger::pmtSoftwareTrigger>();

    // store timestamp of trigger, relative to beam window start
    pmtSoftwareTriggerMetrics->triggerTimestamp = fTriggerTime - beamWindowStart;
    if (fVerbose) std::cout << "Saving trigger timestamp: " << fTriggerTime - beamWindowStart << " ns" << std::endl;

    // add to event
    e.put(std::move(pmtSoftwareTriggerMetrics));      
  }
}

void sbnd::trigger::pmtSoftwareTriggerProducer::checkCAEN1730FragmentTimeStamp(const artdaq::Fragment &frag) {

  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp and check whether in beamspill window
  uint32_t timestamp = md->timeStampNSec;
  if (timestamp >= beamWindowStart && timestamp <= beamWindowEnd) {
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

DEFINE_ART_MODULE(sbnd::trigger::pmtSoftwareTriggerProducer)