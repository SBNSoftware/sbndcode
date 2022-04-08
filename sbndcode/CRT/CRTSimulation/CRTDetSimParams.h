#ifndef SBND_CRTDETSIMPARAMS_H
#define SBND_CRTDETSIMPARAMS_H

#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"

namespace sbnd
{
namespace crt
{
  struct CRTDetSimParams
  {

    fhicl::Atom<double> GlobalT0Offset {
      fhicl::Name("GlobalT0Offset"),
      fhicl::Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    fhicl::Atom<bool> UseG4RefTimeOffset {
      fhicl::Name("UseG4RefTimeOffset"),
      fhicl::Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    fhicl::Atom<double> TDelayNorm {
      fhicl::Name("TDelayNorm"),
      fhicl::Comment("Time delay fit: Gaussian normalization"),
      // true
    };
    fhicl::Atom<double> TDelayShift {
      fhicl::Name("TDelayShift"),
      fhicl::Comment("Time delay fit: Gaussian x shift"),
      // true
    };
    fhicl::Atom<double> TDelaySigma {
      fhicl::Name("TDelaySigma"),
      fhicl::Comment("Time delay fit: Gaussian width"),
      // true
    };
    fhicl::Atom<double> TDelayOffset {
      fhicl::Name("TDelayOffset"),
      fhicl::Comment("Time delay fit: Gaussian baseline offset"),
      // true
    };
    fhicl::Atom<double> TDelayRMSGausNorm {
      fhicl::Name("TDelayRMSGausNorm"),
      fhicl::Comment("Time delay RMS fit: Gaussian normalization"),
      // true
    };
    fhicl::Atom<double> TDelayRMSGausShift {
      fhicl::Name("TDelayRMSGausShift"),
      fhicl::Comment("Time delay RMS fit: Gaussian x shift"),
      // true
    };
    fhicl::Atom<double> TDelayRMSGausSigma {
      fhicl::Name("TDelayRMSGausSigma"),
      fhicl::Comment("Time delay fit: Gaussian width"),
      // true
    };
    fhicl::Atom<double> TDelayRMSExpNorm {
      fhicl::Name("TDelayRMSExpNorm"),
      fhicl::Comment("Time delay RMS fit: Exponential normalization"),
      // true
    };
    fhicl::Atom<double> TDelayRMSExpShift {
      fhicl::Name("TDelayRMSExpShift"),
      fhicl::Comment("Time delay RMS fit: Exponential x shift"),
      // true
    };
    fhicl::Atom<double> TDelayRMSExpScale {
      fhicl::Name("TDelayRMSExpScale"),
      fhicl::Comment("Time delay RMS fit: Exponential scale"),
      // true
    };
    fhicl::Atom<uint32_t> TriggerDelay {
      fhicl::Name("TriggerDelay"),
      fhicl::Comment("Time between signal starts and waveform goes above threshold"),
      // true
    };
    fhicl::Atom<bool> EqualizeSiPMTimes {
      fhicl::Name("EqualizeSiPMTimes"),
      fhicl::Comment("Makes the time simulation to the two SiPMs on a strip identical."),
      false
    };
    fhicl::Atom<double> ClockSpeedCRT {
      fhicl::Name("ClockSpeedCRT"),
      fhicl::Comment("Clock speed for the CRT system [MHz]"),
      // true
    };
    fhicl::Atom<double> NpeScaleNorm {
      fhicl::Name("NpeScaleNorm"),
      fhicl::Comment("Npe vs. distance: 1/r^2 scale"),
      // true
    };
    fhicl::Atom<double> NpeScaleShift {
      fhicl::Name("NpeScaleShift"),
      fhicl::Comment("Npe vs. distance: 1/r^2 x shift"),
      // true
    };
    fhicl::Atom<double> Q0 {
      fhicl::Name("Q0"),
      fhicl::Comment("Average energy deposited for mips, for charge scaling [GeV]"),
      // true
    };
    fhicl::Atom<double> QPed {
      fhicl::Name("QPed"),
      fhicl::Comment("ADC offset for the single-peak peak mean [ADC]"),
      // true
    };
    fhicl::Atom<double> QSlope {
      fhicl::Name("QSlope"),
      fhicl::Comment("Slope in mean ADC / Npe [ADC]"),
      // true
    };
    fhicl::Atom<double> QRMS {
      fhicl::Name("QRMS"),
      fhicl::Comment("ADC single-pe spectrum width [ADC]"),
      // true
    };
    fhicl::Atom<double> QThreshold {
      fhicl::Name("QThreshold"),
      fhicl::Comment("ADC charge threshold [ADC]"),
      // true
    };
    fhicl::Atom<double> TResInterpolator {
      fhicl::Name("TResInterpolator"),
      fhicl::Comment("Interpolator time resolution [ns]"),
      // true
    };
    fhicl::Atom<double> PropDelay {
      fhicl::Name("PropDelay"),
      fhicl::Comment("Delay in pulse arrival time [ns/m]"),
      // true
    };
    fhicl::Atom<double> PropDelayError {
      fhicl::Name("PropDelayError"),
      fhicl::Comment("Delay in pulse arrival time, uncertainty [ns/m]"),
      // true
    };
    fhicl::Atom<double> StripCoincidenceWindow {
      fhicl::Name("StripCoincidenceWindow"),
      fhicl::Comment("Time window for two-fiber coincidence [ns]"),
      // true
    };
    fhicl::Atom<double> TaggerPlaneCoincidenceWindow {
      fhicl::Name("TaggerPlaneCoincidenceWindow"),
      fhicl::Comment("Time window for two-plane coincidence [ticks]"),
      // true
    };
    fhicl::Atom<double> AbsLenEff {
      fhicl::Name("AbsLenEff"),
      fhicl::Comment("Effective abs. length for transverse Npe scaling [cm]"),
      // true
    };
    fhicl::Atom<bool> UseEdep {
      fhicl::Name("UseEdep"),
      fhicl::Comment("Use the true G4 energy deposited, assume mip if false"),
      // true
    };
    fhicl::Atom<double> SipmTimeResponse {
      fhicl::Name("SipmTimeResponse"),
      fhicl::Comment("Minimum time to resolve separate energy deposits [ns]"),
      // true
    };
    fhicl::Atom<uint32_t> AdcSaturation {
      fhicl::Name("AdcSaturation"),
      fhicl::Comment("Saturation limit per SiPM in ADC counts"),
      // true
    };
    fhicl::Atom<double> DeadTime {
      fhicl::Name("DeadTime"),
      fhicl::Comment("Saturation limit per SiPM in ADC counts"),
      // true
    };
    fhicl::Sequence<double> WaveformX {
      fhicl::Name("WaveformX"),
      fhicl::Comment("SiPM waveform sampled, X")
    };
    fhicl::Sequence<double> WaveformY {
      fhicl::Name("WaveformY"),
      fhicl::Comment("SiPM waveform sampled, Y")
    };
    fhicl::Atom<bool> DoWaveformEmulation {
      fhicl::Name("DoWaveformEmulation"),
      fhicl::Comment("Weather or not to perform waveform simulation"),
      true
    };
  };
}
}

#endif