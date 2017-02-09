///////////////////////////////////////////////////////////////////////////////
/// Class: CRTDetSim
/// Module Type: producer
/// File: CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
///
/// Author: mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardata/DetectorInfo/ElecClock.h"

#include <string>

namespace crt {

class CRTDetSim : public art::EDProducer {
public:
  explicit CRTDetSim(fhicl::ParameterSet const & p);

  CRTDetSim(CRTDetSim const &) = delete;
  CRTDetSim(CRTDetSim &&) = delete;
  CRTDetSim& operator = (CRTDetSim const &) = delete;
  CRTDetSim& operator = (CRTDetSim &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p) override;

  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

private:
  /**
   * Get the channel trigger time relative to the start of the MC event.
   *
   * @param engine The random number generator engine
   * @param clock The clock to count ticks on
   * @param t0 The starting time (which delay is added to)
   * @param npe Number of observed photoelectrons
   * @param r Distance between the energy deposit and strip readout end [mm]
   * @return The channel trigger time [ns]
   */
  double getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                detinfo::ElecClock& clock,
                                float t0, float npeMean, float r);

  float fTDelayNorm;  //!< Time delay fit: Gaussian normalization
  float fTDelayShift;  //!< Time delay fit: Gaussian x shift
  float fTDelaySigma;  //!< Time delay fit: Gaussian width
  float fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
  float fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
  float fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
  float fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
  float fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
  float fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
  float fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
  float fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
  float fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
  float fQ0;  //!< Average energy deposited for mips, for charge scaling [GeV]
  float fQPed;  //!< ADC offset for the single-peak peak mean [ADC]
  float fQSlope;  //!< Slope in mean ADC / Npe [ADC]
  float fQRMS;  //!< ADC single-pe spectrum width [ADC]
  float fQThreshold;  //!< ADC charge threshold [ADC]
  float fTResInterpolator;  //!< Interpolator time resolution [ns]
  float fPropDelay;  //!< Delay in pulse arrival time [ns/m]
  float fPropDelayError;  //!< Delay in pulse arrival time, uncertainty [ns/m]
  float fStripCoincidenceWindow;  //!< Time window for two-fiber coincidence [ns]
  float fAbsLenEff;  //!< Effective abs. length for transverse Npe scaling [cm]
};

}  // namespace crt

