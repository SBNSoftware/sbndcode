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
  bool is_arapuca = false;
  std::string opdet_type = fOpDetMap.pdType(channel);
  if (opdet_type == "bar" ||
      opdet_type == "xarapucaprime" ||
      opdet_type == "xarapuca" ||
      opdet_type == "xarapucaT1" ||
      opdet_type == "xarapucaT2" ||
      opdet_type == "arapucaT1" ||
      opdet_type == "arapucaT2") {
    is_arapuca = true;
  }
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

  raw::TimeStamp_t end = tick_to_timestamp(clockData, waveform.TimeStamp(), adcs.size() - 1);
  // get the end time
  // size_t end_i = end > trigger_window[1] ? adcs.size()-1 : (size_t)((trigger_window[1] - start) / optical_period(clockData));
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
  if (opdet_type == "pmt" && fConfig.MaskPMTs()) return true;
  if (opdet_type == "barepmt" && fConfig.MaskBarePMTs()) return true;
  if (opdet_type == "xarapucaprime" && fConfig.MaskXArapucaPrimes()) return true;
  if (opdet_type == "xarapuca" && fConfig.MaskXArapucas()) return true;
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
    start = clockData.TriggerOffsetTPC() - 1 /* Give 1us of wiggle room*/;
    end = start + detProp.ReadOutWindowSize() * clockData.TPCClock().TickPeriod() + 1 /* take back the wiggle room */;
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
  std::vector<raw::OpDetWaveform> ret;
  raw::Channel_t channel = waveform.ChannelNumber();
  // if (channel > (unsigned)fOpDetMap.size()) return {};

  const std::vector<raw::TimeStamp_t> &trigger_times = GetTriggerTimes(channel);

  double readout_window_pre_trigger = ReadoutWindowPreTrigger(channel);
  double readout_window_post_trigger = ReadoutWindowPostTrigger(channel);

  double beam_readout_window_pre_trigger = ReadoutWindowPreTriggerBeam(channel);
  double beam_readout_window_post_trigger = ReadoutWindowPostTriggerBeam(channel);
  double beam_trigger_time = fConfig.BeamTriggerTime();

  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  unsigned trigger_i = 0;
  raw::OpDetWaveform this_waveform;
  bool was_triggering = false;
  for (size_t i = 0; i < adcs.size(); i++) {
    double time = tick_to_timestamp(clockData, waveform.TimeStamp(), i);

    // first, scroll to the next readout window that ends after this time
    while (trigger_i < trigger_times.size() && time >= trigger_times[trigger_i] + readout_window_post_trigger) {
      trigger_i += 1;
    } 
    // see if we are reading out
    bool has_next_trigger = trigger_i < trigger_times.size();
    bool has_beam_trigger = fConfig.BeamTriggerEnable();

    bool is_triggering = false;
    // check next trigger
    if (has_next_trigger) {
      if (time >= trigger_times[trigger_i] - readout_window_pre_trigger &&
          time < trigger_times[trigger_i] + readout_window_post_trigger) {
            is_triggering = true;
      }
    }
    // otherwise check beam trigger
    if (!is_triggering && has_beam_trigger) {
      if (time >= beam_trigger_time - beam_readout_window_pre_trigger &&
          time < beam_trigger_time + beam_readout_window_post_trigger) {
        is_triggering = true;
      }
    }
    // We were triggering and we still are! Add this adc count to the list
    if (is_triggering && was_triggering) {
      this_waveform.push_back(adcs[i]); 
    } 
    // No longer triggering -- save the waveform
    else if (!is_triggering && was_triggering) {
      ret.push_back(std::move(this_waveform));
    }
    // New Trigger! make a new OpDetWaveform
    else if (is_triggering && !was_triggering) {
      this_waveform = raw::OpDetWaveform(time, channel); 
      this_waveform.push_back(adcs[i]);
    }
    // last adc -- save the waveform
    if (is_triggering && i+1 == adcs.size()) {
      ret.push_back(std::move(this_waveform));
    }
 
    was_triggering = is_triggering;
  }
  return ret;
}

} // namespace opdet
