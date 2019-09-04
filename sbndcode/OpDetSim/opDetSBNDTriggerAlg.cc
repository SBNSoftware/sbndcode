#include "opDetSBNDTriggerAlg.h"

// Local static functions
namespace opdet {
// Adds a value to a vector of trigger locations, while keeping the vector sorted
void AddTriggerLocation(std::vector<raw::TimeStamp_t> &triggers, raw::TimeStamp_t trigger) {
  typedef std::vector<raw::TimeStamp_t> TimeStamps;

  TimeStamps::iterator insert = std::lower_bound(triggers.begin(), triggers.end(), trigger);
  triggers.insert(insert, trigger);
}
void AddTriggerLocationAndChannel(std::vector<std::pair<raw::Channel_t, raw::TimeStamp_t>> &triggers, raw::TimeStamp_t trigger, raw::Channel_t channel) {
  typedef std::vector<std::pair<raw::Channel_t, raw::TimeStamp_t>> TimeStamps;

  std::pair<raw::Channel_t, raw::TimeStamp_t> trigger_pair {trigger, channel};

  TimeStamps::iterator insert = std::lower_bound(triggers.begin(), triggers.end(), trigger_pair, 
    [](auto const &lhs, auto const &rhs) { return lhs.first < rhs.first; });
  triggers.insert(insert, trigger_pair);
}

opDetSBNDTriggerAlg::opDetSBNDTriggerAlg(const Config &config, const detinfo::DetectorClocks *detector_clocks, const detinfo::DetectorProperties *detector_properties):
  fConfig(config),
  fDetectorClocks(detector_clocks),
  fDetectorProperties(detector_properties)
{}

void opDetSBNDTriggerAlg::FindTriggerLocations(const raw::OpDetWaveform &waveform, raw::ADC_Count_t baseline) {
  std::vector<raw::TimeStamp_t> this_trigger_locations;
  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  raw::Channel_t channel = waveform.ChannelNumber();
  bool above_threshold = false;
  // find all ADC counts above threshold
  for (size_t i = 0; i < adcs.size(); i++) {
    raw::ADC_Count_t val = fConfig.PulsePolarity() * (adcs[i] - baseline);
    if (!above_threshold && val > fConfig.TriggerThresholdADC()) {
      // new trigger! -- get the time
      raw::TimeStamp_t this_trigger_time = Tick2Timestamp(waveform.TimeStamp(), i);
      if (IsTriggerEnabled(this_trigger_time)) {
        AddTriggerLocation(this_trigger_locations, this_trigger_time);
      }
      above_threshold = true;
    }
    else if (above_threshold && val < fConfig.TriggerThresholdADC()) {
      above_threshold = false;
    }
  }

  // Add in these triggers to the channel
  //
  // Small speed optimization: if this is the first time we are setting the 
  // trigger times for the channel, just move the vector we already built
  if (fTriggerLocationsPerChannel[channel].size() == 0) {
    fTriggerLocationsPerChannel[channel] = std::move(this_trigger_locations);
  }
  // Otherwise, merge them in and keep things sorted in time
  else {
    for (raw::TimeStamp_t trigger: this_trigger_locations) {
      AddTriggerLocation(fTriggerLocationsPerChannel[channel], trigger);
    }
  }

}

bool opDetSBNDTriggerAlg::IsChannelMasked(raw::Channel_t channel) const {
  std::vector<unsigned> masked = fConfig.MaskedChannels();
  // mask by channel number
  bool in_masked_list = std::find(masked.begin(), masked.end(), channel) != masked.end();
  if (in_masked_list) return true;

  // mask by optical detector type
  std::string opdet_type = fOpDetMap.pdName(channel); 
  if (opdet_type == "bar" && fConfig.MaskLightBars() /* RIP */) return true;
  if (opdet_type == "pmt" && fConfig.MaskPMTs()) return true;
  if (opdet_type == "barepmt" && fConfig.MaskBarePMTs()) return true;
  if (opdet_type == "xarapucaprime" && fConfig.MaskXArapucaPrimes()) return true;
  if (opdet_type == "xarapuca" && fConfig.MaskXArapucas()) return true;
  if (opdet_type == "arapucaT1" && fConfig.MaskArapucaT1s()) return true;
  if (opdet_type == "arapucaT2" && fConfig.MaskArapucaT2s()) return true;
  
  return false;
}

void opDetSBNDTriggerAlg::MergeTriggerLocations() {
  // If each channel is self triggered, there is no "master" set of triggers, and 
  // we don't need to do anything here
  if (fConfig.SelfTriggerPerChannel()) return;

  // Otherwise we need to merge in the triggers from each channel into 
  // a master set to be issued to each channel
  //
  // The exact mechanism for how this is done should mirror the implementation
  // in the PTB. However, no such implementation is really defined yet,
  // so we implement a small generic algorithm here. This may likely have
  // to be changed later.

  // First re-sort the trigger times to be a sorted global list of (channel, time) values
  std::vector<std::pair<raw::Channel_t, raw::TimeStamp_t>> all_trigger_locations;
  for (const auto &trigger_locations_pair: fTriggerLocationsPerChannel) {
    raw::Channel_t this_channel = trigger_locations_pair.first;
    // check if this channel contributes to the trigger
    if (!IsChannelMasked(this_channel)) {
      // then add in the locations
      for (raw::TimeStamp_t trigger_time: trigger_locations_pair.second) {
        AddTriggerLocationAndChannel(all_trigger_locations, trigger_time, this_channel);
      }
    }
  }

  // Now merge the trigger locations we have 
  //
  // What is the merging algorithm? This will probably come from the PTB.
  // For now just require some number of channels in some time window.
  // Also allow for a trigger holdoff to be set after one trigger is issued.
  bool has_triggered = false;
  double last_trigger_time;

  bool is_triggering = false; 
  double this_trigger_time;
  unsigned this_trigger_count = 0;

  for (auto const &trigger_pair: all_trigger_locations) {
    raw::TimeStamp_t time = trigger_pair.second;
    if (has_triggered && last_trigger_time + fConfig.TriggerHoldoff() > time) continue;
    if (is_triggering) {
      if (time < this_trigger_time + fConfig.TriggerCountWindow()) {
        this_trigger_count += 1;
      }
      else {
        is_triggering = false;
        this_trigger_count = 0;
      }
    } 
    else {
      is_triggering = true;
      this_trigger_count = 1;
      this_trigger_time = time;
    }
    if (is_triggering && this_trigger_count == fConfig.TriggerChannelCount()) {
      fTriggerLocations.push_back(this_trigger_time);     
      has_triggered = true;
      last_trigger_time = this_trigger_time;
      is_triggering = false;
      this_trigger_count = 0;
    }
  }
}

std::array<double, 2> opDetSBNDTriggerAlg::TriggerEnableWindow() const {
  double start, end;
  // Enable triggers over the range of the TPC readout
  if (fConfig.TriggerEnableWindowIsTPCReadoutWindow()) {
    start = fDetectorClocks->TriggerTime() + fDetectorClocks->TriggerOffsetTPC();
    end = fDetectorProperties->ReadOutWindowSize() * fDetectorClocks->TPCClock().TickPeriod();
  }
  else {
    start = fConfig.TriggerEnableWindowStart();
    end = fConfig.TriggerEnableWindowStart() + fConfig.TriggerEnableWindowLength();
  }

  return {start, end};
}

bool opDetSBNDTriggerAlg::IsTriggerEnabled(raw::TimeStamp_t trigger_time) const {
  // disable triggers during the beam trigger window
  if (fConfig.BeamTriggerEnable()) {
    if (trigger_time > fConfig.BeamTriggerTime() &&
        trigger_time < fConfig.BeamTriggerTime() + fConfig.BeamTriggerHoldoff()) {
      return false; 
    }
  }

  // Enable triggers over the range of the TPC readout
  if (fConfig.TriggerEnableWindowIsTPCReadoutWindow()) {
    raw::TimeStamp_t trigger_window_start = fDetectorClocks->TriggerTime() + fDetectorClocks->TriggerOffsetTPC();
    raw::TimeStamp_t trigger_window_length = fDetectorProperties->ReadOutWindowSize() * fDetectorClocks->TPCClock().TickPeriod();
    return trigger_time > trigger_window_start &&
           trigger_time < trigger_window_start + trigger_window_length;
  }
  // Enable triggers for custom time
  else {
    return trigger_time > fConfig.TriggerEnableWindowStart() &&
           trigger_time < fConfig.TriggerEnableWindowStart() + fConfig.TriggerEnableWindowLength();
  }
}

raw::TimeStamp_t opDetSBNDTriggerAlg::Tick2Timestamp(raw::TimeStamp_t waveform_start, size_t waveform_index) const {
  return waveform_start + waveform_index * fDetectorClocks->OpticalClock().TickPeriod();
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

std::vector<raw::OpDetWaveform> opDetSBNDTriggerAlg::ApplyTriggerLocations(const raw::OpDetWaveform &waveform) const {
  std::vector<raw::OpDetWaveform> ret;
  const std::vector<raw::TimeStamp_t> &trigger_times = GetTriggerTimes(waveform.ChannelNumber());

  double readout_window_pre_trigger = ReadoutWindowPreTrigger(waveform.ChannelNumber());
  double readout_window_post_trigger = ReadoutWindowPostTrigger(waveform.ChannelNumber());

  double beam_readout_window_pre_trigger = ReadoutWindowPreTriggerBeam(waveform.ChannelNumber());
  double beam_readout_window_post_trigger = ReadoutWindowPostTriggerBeam(waveform.ChannelNumber());
  double beam_readout_start_time = fConfig.BeamTriggerTime() - beam_readout_window_post_trigger;

  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  unsigned trigger_i = 0;
  raw::OpDetWaveform this_waveform;
  bool was_triggering = false;
  for (size_t i = 0; i < adcs.size(); i++) {
    double time = Tick2Timestamp(waveform.TimeStamp(), i);

    // first, scroll to the next readout window that ends after this time
    while (trigger_i < trigger_times.size() && time > trigger_times[i] + readout_window_post_trigger) {
      trigger_i += 1;
    } 
    // see if we are reading out
    bool has_next_trigger = trigger_i < trigger_times.size();
    bool has_beam_trigger = fConfig.BeamTriggerEnable();

    bool is_triggering = false;
    // check next trigger
    if (has_next_trigger) {
      if (time > trigger_times[trigger_i] - readout_window_pre_trigger &&
          time < trigger_times[trigger_i] + readout_window_post_trigger) {
            is_triggering = true;
      }
    }
    // otherwise check beam trigger
    if (!is_triggering && has_beam_trigger) {
      if (time > beam_readout_start_time - beam_readout_window_pre_trigger &&
          time < beam_readout_start_time + beam_readout_window_post_trigger) {
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
      this_waveform = raw::OpDetWaveform(time, waveform.ChannelNumber()); 
      this_waveform.push_back(adcs[i]);
    }
    was_triggering = is_triggering;
  }
  return std::move(ret);
}

} // namespace opdet
