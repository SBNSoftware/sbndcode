/**
 * \brief Class for SBND CRT detector simulation parameters
 *
 * \details This class contains all the parameters for the CRT detector simulation.
 * Note that physics parameters do not have default values,
 * and all parameters need to be initialized via fhicl.
 *
 * \author Andy Mastbaum
 * \author Marco Del Tutto
 */

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
      fhicl::Comment("The global time offset to use for the CRT times"),
    };
    fhicl::Atom<bool> UseG4RefTimeOffset {
      fhicl::Name("UseG4RefTimeOffset"),
      fhicl::Comment("Wheater or not to use the G4RefTime as GlobalT0Offset"),
    };
    fhicl::Atom<double> TDelayNorm {
      fhicl::Name("TDelayNorm"),
      fhicl::Comment("Time delay fit: Exponential normalization"),
    };
    fhicl::Atom<double> TDelayScale {
      fhicl::Name("TDelayScale"),
      fhicl::Comment("Time delay fit: Exponential x scale"),
    };
    fhicl::Atom<double> TDelayRMSGausNorm {
      fhicl::Name("TDelayRMSGausNorm"),
      fhicl::Comment("Time delay RMS fit: Gaussian normalization"),
    };
    fhicl::Atom<double> TDelayRMSGausShift {
      fhicl::Name("TDelayRMSGausShift"),
      fhicl::Comment("Time delay RMS fit: Gaussian x shift"),
    };
    fhicl::Atom<double> TDelayRMSGausSigma {
      fhicl::Name("TDelayRMSGausSigma"),
      fhicl::Comment("Time delay fit: Gaussian width"),
    };
    fhicl::Atom<double> TDelayRMSExpShift {
      fhicl::Name("TDelayRMSExpShift"),
      fhicl::Comment("Time delay RMS fit: Exponential x shift"),
    };
    fhicl::Atom<double> TDelayRMSExpScale {
      fhicl::Name("TDelayRMSExpScale"),
      fhicl::Comment("Time delay RMS fit: Exponential scale"),
    };
    fhicl::Atom<double> TDelayRMSOffSetSlope {
      fhicl::Name("TDelayRMSOffSetSlope"),
      fhicl::Comment("Time delay RMS fit: Offset slope"),
    };
    fhicl::Atom<double> TDelayRMSOffSet {
      fhicl::Name("TDelayRMSOffSet"),
      fhicl::Comment("Time delay RMS fit: Offset"),
    };
    fhicl::Atom<uint32_t> TriggerDelay {
      fhicl::Name("TriggerDelay"),
      fhicl::Comment("Time between signal starts and waveform goes above threshold"),
    };
    fhicl::Atom<bool> EqualizeSiPMTimes {
      fhicl::Name("EqualizeSiPMTimes"),
      fhicl::Comment("Makes the time simulation to the two SiPMs on a strip identical"),
      false
    };
    fhicl::Atom<double> ClockSpeedCRT {
      fhicl::Name("ClockSpeedCRT"),
      fhicl::Comment("Clock speed for the CRT system [MHz]"),
    };
    fhicl::Atom<double> NpeScaleNorm {
      fhicl::Name("NpeScaleNorm"),
      fhicl::Comment("Npe vs. distance: 1/r^2 scale"),
    };
    fhicl::Atom<double> NpeScaleShift {
      fhicl::Name("NpeScaleShift"),
      fhicl::Comment("Npe vs. distance: 1/r^2 x shift"),
    };
    fhicl::Atom<double> Q0 {
      fhicl::Name("Q0"),
      fhicl::Comment("Average energy deposited for mips, for charge scaling [GeV]"),
    };
    fhicl::Atom<double> QPed {
      fhicl::Name("QPed"),
      fhicl::Comment("ADC offset for the single-peak peak mean [ADC]"),
    };
    fhicl::Atom<double> QSlope {
      fhicl::Name("QSlope"),
      fhicl::Comment("Slope in mean ADC / Npe [ADC]"),
    };
    fhicl::Atom<double> QRMS {
      fhicl::Name("QRMS"),
      fhicl::Comment("ADC single-pe spectrum width [ADC]"),
    };
    fhicl::Atom<double> QThreshold {
      fhicl::Name("QThreshold"),
      fhicl::Comment("ADC charge threshold [ADC]"),
    };
    fhicl::Atom<double> TResInterpolator {
      fhicl::Name("TResInterpolator"),
      fhicl::Comment("Interpolator time resolution [ns]"),
    };
    fhicl::Atom<double> PropDelay {
      fhicl::Name("PropDelay"),
      fhicl::Comment("Delay in pulse arrival time [ns/m]"),
    };
    fhicl::Atom<double> PropDelayError {
      fhicl::Name("PropDelayError"),
      fhicl::Comment("Delay in pulse arrival time, uncertainty [ns/m]"),
    };
    fhicl::Atom<double> StripCoincidenceWindow {
      fhicl::Name("StripCoincidenceWindow"),
      fhicl::Comment("Time window for two-fiber coincidence [ns]"),
    };
    fhicl::Atom<double> TaggerPlaneCoincidenceWindow {
      fhicl::Name("TaggerPlaneCoincidenceWindow"),
      fhicl::Comment("Time window for two-plane coincidence [ticks]"),
    };
    fhicl::Atom<double> AbsLenEff {
      fhicl::Name("AbsLenEff"),
      fhicl::Comment("Effective abs. length for transverse Npe scaling [cm]"),
    };
    fhicl::Atom<bool> UseEdep {
      fhicl::Name("UseEdep"),
      fhicl::Comment("Use the true G4 energy deposited, assume mip if false"),
    };
    fhicl::Atom<double> SipmTimeResponse {
      fhicl::Name("SipmTimeResponse"),
      fhicl::Comment("Minimum time to resolve separate energy deposits [ns]"),
    };
    fhicl::Atom<uint32_t> AdcSaturation {
      fhicl::Name("AdcSaturation"),
      fhicl::Comment("Saturation limit per SiPM in ADC counts"),
    };
    fhicl::Atom<double> DeadTime {
      fhicl::Name("DeadTime"),
      fhicl::Comment("FEB dead time after a trigger [ns]"),
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
    };
    fhicl::Atom<bool> DebugTrigger {
      fhicl::Name("DebugTrigger"),
      fhicl::Comment("If true, prints out additional debug messages for trigger debugging"),
      false
    };
  };
}
}

#endif