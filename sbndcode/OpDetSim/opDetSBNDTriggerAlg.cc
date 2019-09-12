#include "opDetSBNDTriggerAlg.h"

namespace opdet {

// Local static functions
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
{
  // setup the masked channels
  fConfig.MaskedChannels(fMaskedChannels);
}

void opDetSBNDTriggerAlg::FindTriggerLocations(const raw::OpDetWaveform &waveform, raw::ADC_Count_t baseline) {
  std::vector<raw::TimeStamp_t> this_trigger_locations;
  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  raw::Channel_t channel = waveform.ChannelNumber();
  if (channel > (unsigned)fOpDetMap.size()) return;

  // initialize the channel in the map no matter what
  if (fTriggerLocationsPerChannel.count(channel) == 0) {
    fTriggerLocationsPerChannel[channel] = std::vector<raw::TimeStamp_t>();
  }

  // get the threshold -- first check if channel is Arapuca or PMT
  bool is_arapuca = false;
  std::string opdet_type = fOpDetMap.pdName(channel);
  if (opdet_type == "bar" ||
      opdet_type == "xarapucaprime" ||
      opdet_type == "xarapuca" ||
      opdet_type == "arapucaT1" ||
      opdet_type == "arapucaT2") {
    is_arapuca = true;
  }
  int threshold = is_arapuca ? fConfig.TriggerThresholdADCArapuca() : fConfig.TriggerThresholdADCPMT(); 

  // find the start and end points of the trigger window in this waveform
  std::array<double, 2> trigger_window = TriggerEnableWindow();
  raw::TimeStamp_t start = Tick2Timestamp(waveform.TimeStamp(), 0);

  if (start > trigger_window[1]) return;
  size_t start_i = start > trigger_window[0] ? 0 : (size_t)((trigger_window[0] - start) / OpticalPeriod());

  // fix rounding on division if necessary
  if (!IsTriggerEnabled(Tick2Timestamp(waveform.TimeStamp(), start_i))) {
    start_i += 1;
  }
  assert(IsTriggerEnabled(Tick2Timestamp(waveform.TimeStamp(), start_i)));

  // if start is past end of waveform, we can return
  if (start_i >= adcs.size()) return;

  raw::TimeStamp_t end = Tick2Timestamp(waveform.TimeStamp(), adcs.size() - 1);
  // get the end time
  // size_t end_i = end > trigger_window[1] ? adcs.size()-1 : (size_t)((trigger_window[1] - start) / OpticalPeriod());
  size_t end_i = end < trigger_window[1] ? adcs.size()-1 : (size_t)((trigger_window[1] - start) / OpticalPeriod());
  // fix rounding error...
  if (IsTriggerEnabled(Tick2Timestamp(waveform.TimeStamp(), end_i+1)) && end_i+1 < adcs.size()) {
    end_i += 1;
  }
  assert(end_i+1 == adcs.size() || !IsTriggerEnabled(Tick2Timestamp(waveform.TimeStamp(), end_i+1)));

  bool above_threshold = false;
  // find all ADC counts above threshold
  for (size_t i = start_i; i <= end_i; i++) {
    raw::ADC_Count_t val = fConfig.PulsePolarity() * (adcs.at(i) - baseline);
    if (!above_threshold && val > threshold) {
      // new trigger! -- get the time
      raw::TimeStamp_t this_trigger_time = Tick2Timestamp(waveform.TimeStamp(), i);
      AddTriggerLocation(this_trigger_locations, this_trigger_time);
      above_threshold = true;
    }
    else if (above_threshold && val < threshold) {
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
  // mask by channel number
  bool in_masked_list = std::find(fMaskedChannels.begin(), fMaskedChannels.end(), channel) != fMaskedChannels.end();
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
    // NOTE: assumes trigger is at t=0 in G4 time
    //
    // So currently the photon timing comes from G4 and this configuration relies
    // on the TPC clock. However -- for some reason I don't understand there is 
    // no way to convert a TPC time to a G4 time in the DetectorClock interface....
    //
    // Every reasonable configuration I have seen has setup the trigger time 
    // to be at t=0. If some configuration breaks this invariant, I am sorry, I 
    // did my best.
    start = fDetectorClocks->TriggerOffsetTPC(); 
    end = start + fDetectorProperties->ReadOutWindowSize() * fDetectorClocks->TPCClock().TickPeriod();
  }
  else {
    start = fConfig.TriggerEnableWindowStart();
    end = fConfig.TriggerEnableWindowStart() + fConfig.TriggerEnableWindowLength();
  }

  return {start, end};
}

bool opDetSBNDTriggerAlg::IsTriggerEnabled(raw::TimeStamp_t trigger_time) const {
  // otherwise check the start and end
  std::array<double, 2> trigger_range = TriggerEnableWindow();

  return trigger_time >= trigger_range[0] &&
         trigger_time <= trigger_range[1];
}

double opDetSBNDTriggerAlg::OpticalPeriod() const {
  //TODO: FIX!!!!
  // Currently, the OpticalClock frequency is wrong in SBND configuration
  // This is currently worked-around in the rest of the OpDet simulation by multiplying the clock
  // frequency by a hard-coded factor. However, this should really be fixed in the configuration.
  //
  // For now though, multiply by the same factor to be consistent.
  return fDetectorClocks->OpticalClock().TickPeriod() * 64. / 500.;
}

raw::TimeStamp_t opDetSBNDTriggerAlg::Tick2Timestamp(raw::TimeStamp_t waveform_start, size_t waveform_index) const {
  return waveform_start * fConfig.OpDetWaveformTimeConversion() + waveform_index * OpticalPeriod();
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
  raw::Channel_t channel = waveform.ChannelNumber();
  if (channel > (unsigned)fOpDetMap.size()) return {};

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
    double time = Tick2Timestamp(waveform.TimeStamp(), i);

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
    was_triggering = is_triggering;
  }
  return std::move(ret);
}

} // namespace opdet
