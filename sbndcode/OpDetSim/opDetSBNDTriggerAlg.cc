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
  if( opdet_type.find("arapuca") != std::string::npos ) is_arapuca = true; 
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
  bool beam_trigger_added = false;
  double t_since_last_trigger = 99999.; //[us]
  double t_deadtime = fConfig.TriggerHoldoff();
  raw::TimeStamp_t trigger_start;
  // find all ADC counts above threshold
  for (size_t i = start_i; i <= end_i; i++) {
    raw::TimeStamp_t time = tick_to_timestamp(clockData, waveform.TimeStamp(),i);
    t_since_last_trigger += optical_period(clockData);
    bool isLive = (t_since_last_trigger > t_deadtime);
    raw::ADC_Count_t val = polarity * (adcs.at(i) - baseline);
    // only open new trigger if enough deadtime has passed
    if (isLive && !above_threshold && val > threshold) {
      // new trigger! -- get the time
      trigger_start = time;
      above_threshold = true;
      t_since_last_trigger = 0;
      t_deadtime = fConfig.TriggerHoldoff();
    }
    else if (above_threshold && (val < threshold || i+1 == end_i)) {
      raw::TimeStamp_t trigger_finish = time;
      AddTriggerLocation(this_trigger_locations, {{trigger_start, trigger_finish}});
      above_threshold = false;
    }
    // add beam trigger (if enabled)
    // since the clock ticks might not sync up exactly, use the closet sample
    if( isLive && fConfig.BeamTriggerEnable() && !beam_trigger_added &&
      fabs(time-fConfig.BeamTriggerTime()) <= optical_period(clockData)/2. ){
      AddTriggerLocation(this_trigger_locations, {{time,time}});
      beam_trigger_added = true;
      t_since_last_trigger = 0;
      t_deadtime = fConfig.BeamTriggerHoldoff();
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
  if (opdet_type == "pmt_uncoated" && fConfig.MaskBarePMTs()) return true;
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

  // Enable triggers starting 1 drift period prior to the start of TPC readout,
  // and continuing over the full range of the TPC readout window
  if (fConfig.TriggerEnableWindowOneDriftBeforeTPCReadout()) {
    double tpcReadoutWindow = detProp.ReadOutWindowSize() * clockData.TPCClock().TickPeriod();
    start = clockData.TriggerOffsetTPC() - fConfig.DriftPeriod() - 1; /* Give 1us of wiggle room */
    end = clockData.TriggerOffsetTPC() + tpcReadoutWindow;
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
  return fConfig.ReadoutWindowPreTriggerBeam();
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

//  std::cout
//  <<"=======================\n"
//  <<"Found "<<trigger_times.size()<<" triggers\n";
//  for(size_t i=0; i<trigger_times.size(); i++) std::cout<<"   - "<<trigger_times[i]<<"\n";

  // Extract waveform of raw ADC counts 
  const std::vector<raw::ADC_Count_t> &adcs = waveform; // upcast to get adcs
  raw::OpDetWaveform this_waveform;

  // Set the pre- and post-readout sizes
  double    preTrig	= ReadoutWindowPreTrigger(channel);
  double    postTrig	= ReadoutWindowPostTrigger(channel);
  unsigned  ro_samples  = (preTrig+postTrig)/optical_period(clockData);	// samples

  // Are beam triggers enabled?
  double    beamTrigTime= fConfig.BeamTriggerTime();		// should be 0
  double    preTrigBeam	= ReadoutWindowPreTriggerBeam(channel);
  double    postTrigBeam	= ReadoutWindowPostTriggerBeam(channel);
  unsigned  ro_samples_beam= (preTrigBeam+postTrigBeam)/optical_period(clockData); // samples

  // --------------------------------------------
  // Scan the waveform
  unsigned  trigger_i	= 0;
  double    next_trig	= trigger_times[trigger_i];
  bool	    isReadingOut	= false;
  bool      isBeamTrigger = false;
  unsigned  min_ro_samples = ro_samples;

  for(size_t i=0; i<adcs.size(); i++){
    double time = tick_to_timestamp(clockData, waveform.TimeStamp(), i);

    // if we're out of triggers or nearing the end of the waveform, break
    if( trigger_i >= trigger_times.size() ) break;
    if( i >= adcs.size()-1-min_ro_samples ) break;

    // check if the current trigger is from the beam
    isBeamTrigger = ( fabs(next_trig-beamTrigTime)<optical_period(clockData)/2 );

    // scan ahead to the "next" trigger if we've reached the end of the previous one
    double dT = time-next_trig;
    if( trigger_i < trigger_times.size()-1 &&
        ((isBeamTrigger && dT >= postTrigBeam)||(!isBeamTrigger && dT >= postTrig))) {
      for(size_t j=trigger_i+1; j<trigger_times.size(); j++){
        double this_trig = trigger_times[j];
        isBeamTrigger = ( fabs(this_trig-beamTrigTime)<optical_period(clockData)/2 );
        double t1 = this_trig-preTrig;
        double t2 = this_trig+postTrig;
        if(isBeamTrigger){
          t1 = this_trig-preTrigBeam;
          t2 = this_trig+postTrigBeam;
        }
        if(  ( fConfig.AllowTriggerOverlap() && (t2 >= time) )
           ||(!fConfig.AllowTriggerOverlap() && (t1 >= time) ) ){
	        next_trig = this_trig;
	        trigger_i = j;
	        break;
        }
      }
    }

    // if we're within the trigger window, we're triggering
    bool  isTriggering	= false;
    if(    (!isBeamTrigger && time >= (next_trig-preTrig) &&  time < (next_trig+postTrig))
        || ( isBeamTrigger && time >= (next_trig-preTrigBeam) && time < (next_trig+postTrigBeam))
      )
      isTriggering = true;

    // if already reading out, keep going!
    if( isReadingOut && this_waveform.size() < min_ro_samples ){
      this_waveform.push_back(adcs[i]);
    }
    // once we've saved at least 1 full window, only continue adding
    // to it if we are allowing trigger overlaps. Otherwise, package up 
    // the waveform into an OpDetWaveform object and reset the readout 
    else if( isReadingOut && this_waveform.size() >= min_ro_samples ){
      if( fConfig.AllowTriggerOverlap() && isTriggering ) {
        this_waveform.push_back(adcs[i]);
      } else {
        ret.push_back(std::move(this_waveform));
        isReadingOut = false;
      }
    }

    // start new readout
    if( !isReadingOut && isTriggering ) {
      this_waveform = raw::OpDetWaveform(time, channel); 
      this_waveform.push_back(adcs[i]);
      isReadingOut = true;
      min_ro_samples = ro_samples;
      if( isBeamTrigger ) min_ro_samples = ro_samples_beam;
    }

  }//endloop over ADCs

//  std::cout<<"We saved a total of "<<ret.size()<<" opdetwaveforms\n";
//  for(size_t i=0; i<ret.size(); i++) std::cout<<"   - timestamp "<<ret.at(i).TimeStamp()<<"   size "<<ret.at(i).Waveform().size()<<"\n";
  
  return ret;
}

} // namespace opdet
