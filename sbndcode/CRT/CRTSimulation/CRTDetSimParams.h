#ifndef SBND_CRTDETSIMPARAMS_H
#define SBND_CRTDETSIMPARAMS_H

#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "canvas/Utilities/InputTag.h"
// #include "art/Framework/Core/EDProducer.h"

namespace sbnd
{
namespace crt
{
  struct CRTDetSimParams
  {
    template<class T> using Atom = fhicl::Atom<T>;
    template<class T> using Sequence = fhicl::Sequence<T>;
    template<class T> using Table = fhicl::Table<T>;
    using Comment  = fhicl::Comment;
    using Name     = fhicl::Name;
    using string   = std::string;
    using InputTag = art::InputTag;

    Atom<double> GlobalT0Offset {
      Name("GlobalT0Offset"),
      Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    Atom<bool> UseG4RefTimeOffset {
      Name("UseG4RefTimeOffset"),
      Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    Atom<double> TDelayNorm {
      Name("TDelayNorm"),
      Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    Atom<double> TDelayShift {
      Name("TDelayShift"),
      Comment("Time delay fit: Gaussian x shift"),
      // true
    };
    Atom<double> TDelaySigma {
      Name("TDelaySigma"),
      Comment("Time delay fit: Gaussian width"),
      // true
    };
    Atom<double> TDelayOffset {
      Name("TDelayOffset"),
      Comment("Time delay fit: Gaussian baseline offset"),
      // true
    };
    Atom<double> TDelayRMSGausNorm {
      Name("TDelayRMSGausNorm"),
      Comment("Time delay RMS fit: Gaussian normalization"),
      // true
    };
    Atom<double> TDelayRMSGausShift {
      Name("TDelayRMSGausShift"),
      Comment("Time delay RMS fit: Gaussian x shift"),
      // true
    };
    Atom<double> TDelayRMSGausSigma {
      Name("TDelayRMSGausSigma"),
      Comment("Time delay fit: Gaussian width"),
      // true
    };
    Atom<double> TDelayRMSExpNorm {
      Name("TDelayRMSExpNorm"),
      Comment("Time delay RMS fit: Exponential normalization"),
      // true
    };
    Atom<double> TDelayRMSExpShift {
      Name("TDelayRMSExpShift"),
      Comment("Time delay RMS fit: Exponential x shift"),
      // true
    };
    Atom<double> TDelayRMSExpScale {
      Name("TDelayRMSExpScale"),
      Comment("Time delay RMS fit: Exponential scale"),
      // true
    };
    Atom<uint32_t> TriggerDelay {
      Name("TriggerDelay"),
      Comment("Time between signal starts and waveform goes above threshold"),
      // true
    };
    Atom<double> ClockSpeedCRT {
      Name("ClockSpeedCRT"),
      Comment("Clock speed for the CRT system [MHz]"),
      // true
    };
    Atom<double> NpeScaleNorm {
      Name("NpeScaleNorm"),
      Comment("Npe vs. distance: 1/r^2 scale"),
      // true
    };
    Atom<double> NpeScaleShift {
      Name("NpeScaleShift"),
      Comment("Npe vs. distance: 1/r^2 x shift"),
      // true
    };
    Atom<double> Q0 {
      Name("Q0"),
      Comment("Average energy deposited for mips, for charge scaling [GeV]"),
      // true
    };
    Atom<double> QPed {
      Name("QPed"),
      Comment("ADC offset for the single-peak peak mean [ADC]"),
      // true
    };
    Atom<double> QSlope {
      Name("QSlope"),
      Comment("Slope in mean ADC / Npe [ADC]"),
      // true
    };
    Atom<double> QRMS {
      Name("QRMS"),
      Comment("ADC single-pe spectrum width [ADC]"),
      // true
    };
    Atom<double> QThreshold {
      Name("QThreshold"),
      Comment("ADC charge threshold [ADC]"),
      // true
    };
    Atom<double> TResInterpolator {
      Name("TResInterpolator"),
      Comment("Interpolator time resolution [ns]"),
      // true
    };
    Atom<double> PropDelay {
      Name("PropDelay"),
      Comment("Delay in pulse arrival time [ns/m]"),
      // true
    };
    Atom<double> PropDelayError {
      Name("PropDelayError"),
      Comment("Delay in pulse arrival time, uncertainty [ns/m]"),
      // true
    };
    Atom<double> StripCoincidenceWindow {
      Name("StripCoincidenceWindow"),
      Comment("Time window for two-fiber coincidence [ns]"),
      // true
    };
    Atom<double> TaggerPlaneCoincidenceWindow {
      Name("TaggerPlaneCoincidenceWindow"),
      Comment("Time window for two-plane coincidence [ticks]"),
      // true
    };
    Atom<double> AbsLenEff {
      Name("AbsLenEff"),
      Comment("Effective abs. length for transverse Npe scaling [cm]"),
      // true
    };
    Atom<bool> UseEdep {
      Name("UseEdep"),
      Comment("Use the true G4 energy deposited, assume mip if false"),
      // true
    };
    Atom<double> SipmTimeResponse {
      Name("SipmTimeResponse"),
      Comment("Minimum time to resolve separate energy deposits [ns]"),
      // true
    };
    Atom<uint32_t> AdcSaturation {
      Name("AdcSaturation"),
      Comment("Saturation limit per SiPM in ADC counts"),
      // true
    };
    Atom<double> DeadTime {
      Name("DeadTime"),
      Comment("Saturation limit per SiPM in ADC counts"),
      // true
    };
    fhicl::Sequence<double> WaveformX {
      Name("WaveformX"),
      Comment("SiPM waveform sampled, X")
    };
    fhicl::Sequence<double> WaveformY {
      Name("WaveformY"),
      Comment("SiPM waveform sampled, Y")
    };
    Atom<bool> DoWaveformEmulation {
      Name("DoWaveformEmulation"),
      Comment("Weather or not to perform waveform simulation"),
      true
    };
  };
}
}

#endif