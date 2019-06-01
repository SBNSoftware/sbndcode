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

#include "lardataalg/DetectorInfo/ElecClock.h"

#include <string>

namespace sbnd {
namespace crt {

class CRTDetSim : public art::EDProducer {
public:
  explicit CRTDetSim(fhicl::ParameterSet const & p);

  CRTDetSim(CRTDetSim const &) = delete;
  CRTDetSim(CRTDetSim &&) = delete;
  CRTDetSim& operator = (CRTDetSim const &) = delete;
  CRTDetSim& operator = (CRTDetSim &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p) ;

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
   * @return Trigger clock ticks at this true hit time
   */
  uint32_t getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                detinfo::ElecClock& clock,
                                float t0, float npeMean, float r);

  double fGlobalT0Offset;  //!< Time delay fit: Gaussian normalization
  double fTDelayNorm;  //!< Time delay fit: Gaussian normalization
  double fTDelayShift;  //!< Time delay fit: Gaussian x shift
  double fTDelaySigma;  //!< Time delay fit: Gaussian width
  double fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
  double fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
  double fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
  double fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
  double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
  double fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
  double fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
  double fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
  double fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
  double fQ0;  //!< Average energy deposited for mips, for charge scaling [GeV]
  double fQPed;  //!< ADC offset for the single-peak peak mean [ADC]
  double fQSlope;  //!< Slope in mean ADC / Npe [ADC]
  double fQRMS;  //!< ADC single-pe spectrum width [ADC]
  double fQThreshold;  //!< ADC charge threshold [ADC]
  double fTResInterpolator;  //!< Interpolator time resolution [ns]
  double fPropDelay;  //!< Delay in pulse arrival time [ns/m]
  double fPropDelayError;  //!< Delay in pulse arrival time, uncertainty [ns/m]
  double fStripCoincidenceWindow;  //!< Time window for two-fiber coincidence [ns]
  double fTaggerPlaneCoincidenceWindow;  //!< Time window for two-plane coincidence [ticks]
  double fAbsLenEff;  //!< Effective abs. length for transverse Npe scaling [cm]
  bool fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
  double fSipmTimeResponse; //!< Minimum time to resolve separate energy deposits [ns]
  short fAdcSaturation; //!< Saturation limit per SiPM in ADC counts
  CLHEP::HepRandomEngine& fEngine; //!< Reference to art-managed random-number engine
};

}  // namespace crt
}  // namespace sbnd
