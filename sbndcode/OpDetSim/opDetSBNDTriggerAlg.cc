#include "sbndcode/OpDetSim/opDetSBNDTriggerAlg.hh"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

namespace {
  double optical_period(detinfo::DetectorClocksData const& clockData)
  {
    return clockData.OpticalClock().TickPeriod();
  }

  raw::TimeStamp_t tick_to_timestamp(detinfo::DetectorClocksData const& clockData,
                                     raw::TimeStamp_t waveform_start,
                                     size_t waveform_index)
  {
    return waveform_start + waveform_index * optical_period(clockData);
  }
}

namespace opdet {

struct TriggerPrimitive {
  raw::TimeStamp_t start;
  raw::TimeStamp_t finish;
  raw::Channel_t channel;
};

// Local static functions
// Adds a value to a vector of trigger locations, while keeping the vector sorted
void AddTriggerLocation(std::vector<std::array<raw::TimeStamp_t, 2>> &triggers, std::array<raw::TimeStamp_t,2> range) {
  typedef std::vector<std::array<raw::TimeStamp_t,2>> TimeStamps;

  TimeStamps::iterator insert = std::lower_bound(triggers.begin(), triggers.end(), range,
    [](const auto &lhs, const auto &rhs) { return lhs[0] < rhs[0]; });
  triggers.insert(insert, range);
}

void AddTriggerPrimitive(std::vector<TriggerPrimitive> &triggers, TriggerPrimitive trigger) {
  typedef std::vector<TriggerPrimitive> TimeStamps;

  TimeStamps::iterator insert = std::lower_bound(triggers.begin(), triggers.end(), trigger,
    [](auto const &lhs, auto const &rhs) { return lhs.start < rhs.start; });

  triggers.insert(insert, trigger);
}

void AddTriggerPrimitiveFinish(std::vector<TriggerPrimitive> &triggers, TriggerPrimitive trigger) {
  typedef std::vector<TriggerPrimitive> TimeStamps;
  
  TimeStamps::iterator insert = std::upper_bound(triggers.begin(), triggers.end(), trigger,
    [](auto const &lhs, auto const &rhs) { return lhs.finish < rhs.finish; });

  triggers.insert(insert, trigger);
}

opDetSBNDTriggerAlg::opDetSBNDTriggerAlg(const Config &config):
  fConfig(config)
{
  // setup the masked channels
  fConfig.MaskedChannels(fMaskedChannels);
}

void opDetSBNDTriggerAlg::FindTriggerLocations(detinfo::DetectorClocksData const& clockData,
                                               detinfo::DetectorPropertiesData const& detProp,
                                               const raw::OpDetWaveform &waveform, raw::ADC_Count_t baseline) {
  std::vector<std::array<raw::TimeStamp_t, 2>> this_trigger_ranges;
  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  raw::Channel_t channel = waveform.ChannelNumber();
  // if (channel > (unsigned)fOpDetMap.size()) return;

  // initialize the channel in the map no matter what
  if (fTriggerRangesPerChannel.count(channel) == 0) {
    fTriggerRangesPerChannel[channel] = std::vector<std::array<raw::TimeStamp_t,2>>();
  }
  
  // get the threshold -- first check if channel is Arapuca or PMT
  std::string opdet_type = fOpDetMap.pdType(channel);
  bool is_arapuca = false;
  if( opdet_type.find("xarapuca") != std::string::npos || 
      opdet_type.find("arapuca")  != std::string::npos ){
    is_arapuca = true; }
  int threshold = is_arapuca ? fConfig.TriggerThresholdADCArapuca() : fConfig.TriggerThresholdADCPMT(); 
  int polarity = is_arapuca ? fConfig.PulsePolarityArapuca() : fConfig.PulsePolarityPMT(); 

  // find the start and end points of the trigger window in this waveform
  std::array<double, 2> trigger_window = TriggerEnableWindow(clockData, detProp);
  raw::TimeStamp_t start = tick_to_timestamp(clockData, waveform.TimeStamp(), 0);

  if (start > trigger_window[1]) return;
  size_t start_i = start > trigger_window[0] ? 0 : (size_t)((trigger_window[0] - start) / optical_period(clockData));

  // fix rounding on division if necessary
  if (!IsTriggerEnabled(clockData,
                        detProp,
                        tick_to_timestamp(clockData, waveform.TimeStamp(), start_i))) {
    start_i += 1;
  }
  assert(IsTriggerEnabled(clockData,
                          detProp,
                          tick_to_timestamp(clockData, waveform.TimeStamp(), start_i)));

  // if start is past end of waveform, we can return
  if (start_i >= adcs.size()) return;

  // get the end time
  raw::TimeStamp_t end = tick_to_timestamp(clockData, waveform.TimeStamp(), adcs.size() - 1);
  size_t end_i = end < trigger_window[1] ? adcs.size()-1 : (size_t)((trigger_window[1] - start) / optical_period(clockData));
 
  // fix rounding error...
  if (IsTriggerEnabled(clockData,
                       detProp,
                       tick_to_timestamp(clockData, waveform.TimeStamp(), end_i+1)) && end_i+1 < adcs.size()) {
    end_i += 1;
  }
  assert(end_i+1 == adcs.size() || !IsTriggerEnabled(clockData,
                                                     detProp,
                                                     tick_to_timestamp(clockData, waveform.TimeStamp(), end_i+1)));

  std::vector<std::array<raw::TimeStamp_t, 2>> this_trigger_locations; 
  bool above_threshold = false;
  raw::TimeStamp_t trigger_start;
  // find all ADC counts above threshold
  for (size_t i = start_i; i <= end_i; i++) {
    raw::ADC_Count_t val = polarity * (adcs.at(i) - baseline);
    if (!above_threshold && val > threshold) {
      // new trigger! -- get the time
      // raw::TimeStamp_t this_trigger_time 
      trigger_start = tick_to_timestamp(clockData, waveform.TimeStamp(), i);
      above_threshold = true;
    }
    else if (above_threshold && (val < threshold || i+1 == end_i)) {
      raw::TimeStamp_t trigger_finish = tick_to_timestamp(clockData, waveform.TimeStamp(), i);
      AddTriggerLocation(this_trigger_locations, {{trigger_start, trigger_finish}});
      above_threshold = false;
    }
  }

  // Add in these triggers to the channel
  //
  // Small speed optimization: if this is the first time we are setting the 
  // trigger times for the channel, just move the vector we already built
  if (fTriggerRangesPerChannel[channel].size() == 0) {
    fTriggerRangesPerChannel[channel] = std::move(this_trigger_locations);
  }
  // Otherwise, merge them in and keep things sorted in time
  else {
    for (const std::array<raw::TimeStamp_t, 2> &trigger_range: this_trigger_ranges) {
      AddTriggerLocation(fTriggerRangesPerChannel[channel], trigger_range);
    }
  }

}

bool opDetSBNDTriggerAlg::IsChannelMasked(raw::Channel_t channel) const {
  // mask by channel number
  bool in_masked_list = std::find(fMaskedChannels.begin(), fMaskedChannels.end(), channel) != fMaskedChannels.end();
  if (in_masked_list) return true;

  // mask by optical detector type
  std::string opdet_type = fOpDetMap.pdType(channel);
  if (opdet_type == "bar" && fConfig.MaskLightBars() /* RIP */) return true;
  if (opdet_type == "pmt_coated" && fConfig.MaskPMTs()) return true;
  if (opdet_type == "bar_uncoated" && fConfig.MaskBarePMTs()) return true;
  if (opdet_type == "xarapucaprime" && fConfig.MaskXArapucaPrimes()) return true;
  if (opdet_type == "xarapuca" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "xarapuca_vuv" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "xarapuca_vis" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "xarapucaT1" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "xarapucaT2" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "arapucaT1" && fConfig.MaskArapucaT1s()) return true;
  if (opdet_type == "arapucaT2" && fConfig.MaskArapucaT2s()) return true;
  return false;
}

void opDetSBNDTriggerAlg::ClearTriggerLocations() {
  fTriggerLocationsPerChannel.clear();
  fTriggerRangesPerChannel.clear();
  fTriggerLocations.clear();
}

void opDetSBNDTriggerAlg::MergeTriggerLocations() {
  // If each channel is self triggered, there is no "master" set of triggers, and 
  // we don't need to do anything here
  if (fConfig.SelfTriggerPerChannel()) {
    for (const auto &trigger_ranges: fTriggerRangesPerChannel) {
      fTriggerLocationsPerChannel[trigger_ranges.first] = std::vector<raw::TimeStamp_t>();
      for (const std::array<raw::TimeStamp_t, 2> &range: trigger_ranges.second) {
        fTriggerLocationsPerChannel[trigger_ranges.first].push_back(range[0]);
      }
    }
    return;
  }

  // Otherwise we need to merge in the triggers from each channel into 
  // a master set to be issued to each channel
  //
  // The exact mechanism for how this is done should mirror the implementation
  // in the PTB. However, no such implementation is really defined yet,
  // so we implement a small generic algorithm here. This may likely have
  // to be changed later.

  // First re-sort the trigger times to be a sorted global list of (channel, time) values
  std::vector<TriggerPrimitive> all_trigger_locations;
  for (const auto &trigger_locations_pair: fTriggerRangesPerChannel) {
    raw::Channel_t this_channel = trigger_locations_pair.first;
    // check if this channel contributes to the trigger
    if (!IsChannelMasked(this_channel)) {
      // then add in the locations
      for (std::array<raw::TimeStamp_t,2> trigger_range: trigger_locations_pair.second) {
        TriggerPrimitive trigger;
        trigger.start = trigger_range[0];
        trigger.finish = trigger_range[1];
        trigger.channel = this_channel;
        AddTriggerPrimitive(all_trigger_locations, trigger);
      }
    }
  }

  // Now merge the trigger locations we have 
  //
  // What is the merging algorithm? This will probably come from the PTB.
  // For now just require some number of channels in some time window.
  // Also allow for a trigger holdoff to be set after one trigger is issued.
  //
  // Also here is where one might start to simulate a different clock speed
  // for the trigger, and for simulating drifts between channels and/or 
  // between channel and trigger. However, for now we assume that the 
  // trigger and readout have the same clock and that they are perfectly
  // synched.
  bool was_triggering = false;

  std::vector<TriggerPrimitive> primitives;
  for (const TriggerPrimitive &primitive: all_trigger_locations) {
    AddTriggerPrimitiveFinish(primitives, primitive);
    while (primitives.back().finish < primitive.start) {
      // remove the final element
      primitives.resize(primitives.size() - 1);
    }

    bool is_triggering = primitives.size() >= fConfig.TriggerChannelCount();
    if (is_triggering && !was_triggering) {
      raw::TimeStamp_t this_trigger_time = primitive.start;
      fTriggerLocations.push_back(this_trigger_time);
    }
    was_triggering = is_triggering;
  }
}

std::array<double, 2> opDetSBNDTriggerAlg::TriggerEnableWindow(detinfo::DetectorClocksData const& clockData,
                                                               detinfo::DetectorPropertiesData const& detProp) const {
  double start, end;
  // Enable triggers over the range of the TPC readout
  if (fConfig.TriggerEnableWindowIsTPCReadoutWindow()) {
    // NOTE: assumes trigger is at t=0 in G4 time
    //
    // So currently the photon timing comes from G4 and this configuration relies
    // on the TPC clock. However -- for some reason I don't understand there is 
    // no way to convert a TPC time to a G4 time in the DetectorClock interface....
    //
    // Every reasonable configuration I have seen has setup the trigger time 
    // to be at t=0. If some configuration breaks this invariant, I am sorry, I 
    // did my best.
    //start = clockData.TriggerOffsetTPC() - 1; /* Give 1us of wiggle room*/;
    //end = start + detProp.ReadOutWindowSize() * clockData.TPCClock().TickPeriod() + 1; /* take back the wiggle room */;
    
    // Mar 2021:
    // Enable triggers for 1 full drift period prior to + throughout main drift window.
    // For some reason, drift period is not accessible from ClocksData nor DetProp... (?!)
    // Hard-coding for now (TODO: FIX!!!)
    double drift  = 1300; // 1.3 ms
    double readout= detProp.ReadOutWindowSize() * clockData.TPCClock().TickPeriod();
    start = clockData.TriggerOffsetTPC() - drift;
    end = clockData.TriggerOffsetTPC() + readout;
  }
  else {
    start = fConfig.TriggerEnableWindowStart();
    end = fConfig.TriggerEnableWindowStart() + fConfig.TriggerEnableWindowLength();
  }
    
  return {{start, end}};
}

bool opDetSBNDTriggerAlg::IsTriggerEnabled(detinfo::DetectorClocksData const& clockData,
                                           detinfo::DetectorPropertiesData const& detProp,
                                           raw::TimeStamp_t trigger_time) const {
  // otherwise check the start and end
  std::array<double, 2> trigger_range = TriggerEnableWindow(clockData, detProp);

  return trigger_time >= trigger_range[0] &&
         trigger_time <= trigger_range[1];
}

const std::vector<raw::TimeStamp_t> &opDetSBNDTriggerAlg::GetTriggerTimes(raw::Channel_t channel) const {
  if (fConfig.SelfTriggerPerChannel()) {
    return fTriggerLocationsPerChannel.at(channel); 
  }
  return fTriggerLocations;

}

double opDetSBNDTriggerAlg::ReadoutWindowPreTrigger(raw::Channel_t channel) const {
  // Allow for different channels to have different readout window lengths
  // For now, we don't use this
  (void) channel;
  return fConfig.ReadoutWindowPreTrigger();
}

double opDetSBNDTriggerAlg::ReadoutWindowPostTrigger(raw::Channel_t channel) const {
  // Allow for different channels to have different readout window lengths
  // For now, we don't use this
  (void) channel;
  return fConfig.ReadoutWindowPostTrigger();
}

double opDetSBNDTriggerAlg::ReadoutWindowPreTriggerBeam(raw::Channel_t channel) const {
  // Allow for different channels to have different readout window lengths
  // For now, we don't use this
  (void) channel;
  return fConfig.ReadoutWindowPreTrigger();
}

double opDetSBNDTriggerAlg::ReadoutWindowPostTriggerBeam(raw::Channel_t channel) const {
  // Allow for different channels to have different readout window lengths
  // For now, we don't use this
  (void) channel;
  return fConfig.ReadoutWindowPostTriggerBeam();
}

std::vector<raw::OpDetWaveform> opDetSBNDTriggerAlg::ApplyTriggerLocations(detinfo::DetectorClocksData const& clockData,
                                                                           const raw::OpDetWaveform &waveform) const {
  
  // Vector of "triggered" OpDetWaveforms
  std::vector<raw::OpDetWaveform> ret;
  
  // Get the trigger times we found earlier for this channel
  raw::Channel_t channel = waveform.ChannelNumber();
  const std::vector<raw::TimeStamp_t> &trigger_times = GetTriggerTimes(channel);
  if( trigger_times.size() == 0 ) return ret;

  // Extract waveform of raw ADC counts 
  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  raw::OpDetWaveform this_waveform;

  // Set the pre- and post-readout sizes
  // NOTE: Using same readout window size for everything
  //       as of Mar 2021. Eventually we might want to shrink
  //       the readout for out-of-time stuff.
  double    preTrig	= ReadoutWindowPreTrigger(channel);
  double    postTrig	= ReadoutWindowPostTrigger(channel);
  double    ro_length	= preTrig+postTrig;			// microseconds
  unsigned  ro_samples  = ro_length/optical_period(clockData);	// samples
  double    beamTrigTime= fConfig.BeamTriggerTime();		// should be 0
  //double preTrigBeam	= ReadoutWindowPreTriggerBeam(channel); // (not used)
  //double postTrigBeam	= ReadoutWindowPostTriggerBeam(channel);// (not used)
  //double ro_length_beam = preTrigBeam+postTrigBeam;		// (not used)
  //unsigned ro_samples_beam= ro_length_beam/optical_period(clockData); 
  
  // Are beam triggers enabled?
  bool	has_beam_trigger = fConfig.BeamTriggerEnable();

  // --------------------------------------------
  // Scan the waveform
  unsigned trigger_i	= 0;
  double  next_trig	= trigger_times[trigger_i];
  double  end_time	= tick_to_timestamp(clockData,waveform.TimeStamp(), adcs.size()-1);
  bool	  isReadingOut	= false;

  for(size_t i=0; i<adcs.size(); i++){
    double time = tick_to_timestamp(clockData, waveform.TimeStamp(), i);

    // if we're out of triggers or nearing the end of the waveform, break
    if( trigger_i >= trigger_times.size() ) break;
    if( time >= end_time - ro_length ) break;
    
    bool  isTriggering	= false;
    // if we're near a trigger time, we're triggering
    if( time >= (next_trig-preTrig) &&  time < (next_trig+postTrig))
      isTriggering = true;
    //
    // if beam trigger is enabled, then we will ALWAYS read out
    // a window around the beam time, so ensure that no triggers 
    // prior to this overlap into this region
    if( has_beam_trigger ){
      if( time >= (beamTrigTime-preTrig-ro_length) &&
	  time < beamTrigTime-preTrig ) {
	isTriggering = false; 
      } else 
      if( time >= (beamTrigTime-preTrig) && time < (beamTrigTime+postTrig)) {
	isTriggering = true;
      }
    }

    // start new readout
    if( !isReadingOut && isTriggering ) {
      this_waveform = raw::OpDetWaveform(time, channel); 
      this_waveform.push_back(adcs[i]);
      isReadingOut = true;
    }
    // if alraedy reading out, keep going!
    else if( isReadingOut && this_waveform.size() < ro_samples ){
      this_waveform.push_back(adcs[i]);
    } 
    // once we've saved a full readout window, package the 
    // waveform and reset to start looking for next trigger
    else if( isReadingOut && this_waveform.size() >= ro_samples ){
      ret.push_back(std::move(this_waveform));
      isReadingOut = false;
      // keep in mind some triggers could have happened during
      // the readout, so skip ahead to the next one coming up.
      for(size_t j=trigger_i; j<trigger_times.size(); j++){
	if( trigger_times.at(j) >= time ){
	  next_trig = trigger_times.at(j);
	  trigger_i = j;
	  break;
	}
      }
    }

  }//endloop over ADCs

  //std::cout<<"We saved a total of "<<ret.size()<<" opdetwaveforms\n";
  //for(size_t i=0; i<ret.size(); i++) std::cout<<"   - timestamp "<<ret.at(i).TimeStamp()<<"   size "<<ret.at(i).Waveform().size()<<"\n";
  
  return ret;
}

} // namespace opdet
