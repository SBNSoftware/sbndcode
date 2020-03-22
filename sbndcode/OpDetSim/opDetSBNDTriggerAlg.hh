////////////////////////////////////////////////////////////////////////
//// File:        opDetSBNDTriggerAlg.hh
////
//// This algorithm emulates the behavior of the SBND trigger going
//// to the photon detection system. Created by Gray Putnam
//// <grayputnam@uchicago.edu>
//////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_OPDETSBNDTRIGGERALG_HH
#define SBND_OPDETSIM_OPDETSBNDTRIGGERALG_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1D.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet {

  class opDetSBNDTriggerAlg {

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<int> PulsePolarityPMT {
        Name("PulsePolarityPMT"),
        Comment("Whether pulses go down (-1) or up (1)."),
        -1
      };

      fhicl::Atom<int> PulsePolarityArapuca {
        Name("PulsePolarityArapuca"),
        Comment("Whether pulses go down (-1) or up (1)."),
        1
      };

      fhicl::Atom<int> TriggerThresholdADCArapuca {
        Name("TriggerThresholdADCArapuca"),
        Comment("Threshold for channel to issue trigger across all arapuca channels. [ADC]")
      };

      fhicl::Atom<int> TriggerThresholdADCPMT {
        Name("TriggerThresholdADCPMT"),
        Comment("Threshold for channel to issue trigger across all pmt channels. [ADC]")
      };

      fhicl::OptionalSequence<unsigned> MaskedChannels {
        Name("MaskedChannels"),
        Comment("Channels which are ignored for issuing triggers.")
      };

      fhicl::Atom<bool> MaskLightBars {
        Name("MaskLightBars"),
        Comment("Whether to mask all light bar readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskPMTs {
        Name("MaskPMTs"),
        Comment("Whether to mask all PMT readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskBarePMTs {
        Name("MaskBarePMTs"),
        Comment("Whether to mask all Bare PMT readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskXArapucaPrimes {
        Name("MaskXArapucaPrimes"),
        Comment("Whether to mask all Arapuca Prime readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskXArapucas {
        Name("MaskXArapucas"),
        Comment("Whether to mask all X-Arapuca readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskArapucaT1s {
        Name("MaskArapucaT1s"),
        Comment("Whether to mask all Arapuca T1 readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> MaskArapucaT2s {
        Name("MaskArapucaT1s"),
        Comment("Whether to mask all Arapuca T2 readout channels from issuing triggers."),
        false
      };

      fhicl::Atom<bool> SelfTriggerPerChannel {
        Name("SelfTriggerPerChannel"),
        Comment("If true, each channel's OpDetWaveform will be 'triggered' when that channel goes above threshold (this setting ignores all mask settings). If false, a group of global triggers to be issued across all channel's will instead be collected."),
        true
      };

      fhicl::Atom<double> TriggerHoldoff {
        Name("TriggerHoldoff"),
        Comment("When collecting global triggers, sets an amount of time where one trigger will not be issued following the previous one. [us]"),
        0.
      };

      fhicl::Atom<double> TriggerCountWindow {
        Name("TriggerCountWindow"),
        Comment("Size of window to count triggers arriving from different readout channels. [us]"),
        0.
      };

      fhicl::Atom<unsigned> TriggerChannelCount {
        Name("TriggerChannelCount"),
        Comment("Number of trigger signals from different readout channels required to issue global trigger."),
        1
      };

      fhicl::Atom<bool> BeamTriggerEnable {
        Name("BeamTriggerEnable"),
        Comment("Whether to also send a beam trigger."),
        false
      };

      fhicl::Atom<double> BeamTriggerTime {
        Name("BeamTriggerTime"),
        Comment("Time at which the beam trigger will be issued. [us]"),
        0.
      };

      fhicl::Atom<double> BeamTriggerHoldoff {
        Name("BeamTriggerHoldoff"),
        Comment("Time to ignore other triggers after the beam trigger time. [us]"),
        0.
      };

      fhicl::Atom<bool> TriggerEnableWindowIsTPCReadoutWindow {
        Name("TriggerEnableWindowIsTPCReadoutWindow"),
        Comment("If true, triggers will be enabled for the length of time that the TPC is read-out (DetectorClocks and DetectorProperties are used to determine this time span. Otherwise, the config options TriggerEnableWindowStart and TriggerEnableWindowLength are checked."),
        true
      };

      fhicl::Atom<double> TriggerEnableWindowStart {
        Name("TriggerEnableWindowStart"),
        Comment("Start of the window of time for which triggers are enabled. Ignored if TriggerEnableWindowIsTPCReadoutWindow is true. [us]"),
        0.
      };

      fhicl::Atom<double> TriggerEnableWindowLength {
        Name("TriggerEnableWindowLength"),
        Comment("Length of time for which trigger are enabled. Ignored if TriggerEnableWindowIsTPCReadoutWindow is true. [us]"),
        0.
      };

      fhicl::Atom<double> ReadoutWindowPreTrigger {
        Name("ReadoutWindowPreTrigger"),
        Comment("Size of readout window prior to a non-beam trigger. [us]"),
        0.
      };

      fhicl::Atom<double> ReadoutWindowPostTrigger {
        Name("ReadoutWindowPostTrigger"),
        Comment("Size of readout window post a non-beam trigger. [us]"),
        0.
      };

      fhicl::Atom<double> ReadoutWindowPreTriggerBeam {
        Name("ReadoutWindowPreTriggerBeam"),
        Comment("Size of readout window pre a beam trigger. [us]"),
        0.
      };

      fhicl::Atom<double> ReadoutWindowPostTriggerBeam {
        Name("ReadoutWindowPostTriggerBeam"),
        Comment("Size of readout window post a beam trigger. [us]"),
        0.
      };
    };

    // construct with config or with fhicl
    opDetSBNDTriggerAlg(const Config &config, const detinfo::DetectorClocks *detector_clocks, const detinfo::DetectorProperties *detector_properties);

    opDetSBNDTriggerAlg(const fhicl::ParameterSet &pset,  const detinfo::DetectorClocks *detector_clocks, const detinfo::DetectorProperties *detector_properties) :
      opDetSBNDTriggerAlg(fhicl::Table<Config>(pset, {})(), detector_clocks, detector_properties)
    {}

    //Default destructor
    ~opDetSBNDTriggerAlg() {}

    // Clear out at the end of an event
    void ClearTriggerLocations();

    // Add in a waveform to define trigger locations
    void FindTriggerLocations(const raw::OpDetWaveform &waveform, raw::ADC_Count_t baseline);

    // Merge all of the triggers together
    void MergeTriggerLocations();

    // Apply trigger locations to an input OpDetWaveform
    std::vector<raw::OpDetWaveform> ApplyTriggerLocations(const raw::OpDetWaveform &waveform) const;

    // Returns the time range over which triggers are enabled over a range [start, end]
    std::array<double, 2> TriggerEnableWindow() const;

  private:

    // internal functions
    bool IsChannelMasked(raw::Channel_t channel) const;
    bool IsTriggerEnabled(raw::TimeStamp_t trigger_time) const;
    raw::TimeStamp_t Tick2Timestamp(raw::TimeStamp_t waveform_start, size_t waveform_index) const;
    const std::vector<raw::TimeStamp_t> &GetTriggerTimes(raw::Channel_t channel) const;
    double ReadoutWindowPreTrigger(raw::Channel_t channel) const;
    double ReadoutWindowPostTrigger(raw::Channel_t channel) const;
    double ReadoutWindowPreTriggerBeam(raw::Channel_t channel) const;
    double ReadoutWindowPostTriggerBeam(raw::Channel_t channel) const;
    double OpticalPeriod() const;

    // fhicl config
    Config fConfig;

    // pointers to services
    detinfo::DetectorClocks const *fDetectorClocks;
    detinfo::DetectorProperties const *fDetectorProperties;

    // OpDet channel map
    opdet::sbndPDMapAlg fOpDetMap;

    // keeping track of triggers
    std::map<raw::Channel_t, std::vector<std::array<raw::TimeStamp_t, 2>>> fTriggerRangesPerChannel;
    std::map<raw::Channel_t, std::vector<raw::TimeStamp_t>> fTriggerLocationsPerChannel;
    std::vector<raw::TimeStamp_t> fTriggerLocations;

    std::vector<unsigned> fMaskedChannels;

  };//class opDetSBNDTriggerAlg


} //namespace

#endif // SBND_OPDETSIM_OPDETSBNDTRIGGERALG_HH
