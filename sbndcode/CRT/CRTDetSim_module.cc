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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/ElecClock.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TFile.h"
#include "TNtuple.h"
#include "CRTData.hh"

#include <cmath>
#include <memory>

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
  float fQ0;  // Average energy deposited for mips, for charge scaling [GeV]
  float fQPed;  // ADC offset for the single-peak peak mean [ADC]
  float fQSlope;  // Slope in mean ADC / Npe [ADC]
  float fQRMS;  // ADC single-pe spectrum width [ADC]
  float fTResInterpolator;  // Interpolator time resolution [ns]
  float fPropDelay;  // Delay in pulse arrival time [ns/m]
  float fPropDelayError;  // Delay in pulse arrival time, uncertainty [ns/m]
  float fAbsLenEff;  // Effective abs. length for transverse Npe scaling [cm]
};


void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");

  fTDelayNorm = p.get<double>("TDelayNorm");
  fTDelayShift = p.get<double>("TDelayShift");
  fTDelaySigma = p.get<double>("TDelaySigma");
  fTDelayOffset = p.get<double>("TDelayOffset");
  fTDelayRMSGausNorm = p.get<double>("TDelayRMSGausNorm");
  fTDelayRMSGausShift = p.get<double>("TDelayRMSGausShift");
  fTDelayRMSGausSigma = p.get<double>("TDelayRMSGausSigma");
  fTDelayRMSExpNorm = p.get<double>("TDelayRMSExpNorm");
  fTDelayRMSExpShift = p.get<double>("TDelayRMSExpShift");
  fTDelayRMSExpScale = p.get<double>("TDelayRMSExpScale");
  fPropDelay = p.get<double>("PropDelay");
  fPropDelayError = p.get<double>("fPropDelayError");
  fTResInterpolator = p.get<double>("TResInterpolator");
  fNpeScaleNorm = p.get<double>("NpeScaleNorm");
  fNpeScaleShift = p.get<double>("NpeScaleShift");
  fQ0 = p.get<double>("Q0");
  fQPed = p.get<double>("QPed");
  fQSlope = p.get<double>("QSlope");
  fQRMS = p.get<double>("QRMS");
  fAbsLenEff = p.get<double>("AbsLenEff");
}


CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p) {
  art::ServiceHandle<rndm::NuRandomService> seeds;
  seeds->createEngine(*this, "HepJamesRandom", "crt", p, "Seed");

  this->reconfigure(p);

  produces<std::vector<crt::CRTData> >();
}


double CRTDetSim::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                         detinfo::ElecClock& clock,
                                         float t0, float npeMean, float r) {
  // Hit timing, with smearing and NPE dependence
  double tDelayMean = \
    fTDelayNorm *
      exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
    fTDelayOffset;

  double tDelayRMS = \
    fTDelayRMSGausNorm *
      exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
    fTDelayRMSExpNorm *
      exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);

  double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

  // Time resolution of the interpolator
  tDelay += CLHEP::RandGauss::shoot(engine, 0, fTResInterpolator);

  // Propagation time
  double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;

  double t = t0 + tProp + tDelay;

  // Get clock ticks
  clock.SetTime(t);
  return clock.Ticks();
}


void CRTDetSim::produce(art::Event & e) {
  std::unique_ptr<std::vector<crt::CRTData> > crtHits(
      new std::vector<crt::CRTData>);

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;

  art::ServiceHandle<detinfo::DetectorClocksService> detClocks;

  detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine("crt");

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  // Loop through truth AD channels
  for (auto& adsc : *channels) {
    const geo::AuxDetGeo& adGeo = \
        geoService->AuxDet(adsc.AuxDetID());

    const geo::AuxDetSensitiveGeo& adsGeo = \
        adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

    // Simulate the CRT response for each hit
    for (auto ide : adsc.AuxDetIDEs()) {
      // Get the hit position in strip's local coordinates
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      double world[3] = {x, y, z};
      double svHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal);

      // Distance to the readout end ("outside") depends on module position
      // FIXME: FOR NOW ASSUME ALL THE SAME DIRECTION
      double distToReadout = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);

      // The expected number of PE
      double qr = ide.energyDeposited / fQ0;  // Scale linearly with charge
      double npeExpected = \
        fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr;

      // Put PE on channels weighted by distance
      double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
      double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
      double abs0 = exp(-d0 / fAbsLenEff);
      double abs1 = exp(-d1 / fAbsLenEff);
      double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
      double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

      // Observed PE
      long npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
      long npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);

      // Time relative to trigger
      double tTrue = (ide.entryT + ide.exitT) / 2;
      uint32_t t0 = \
        getChannelTriggerTicks(engine, trigClock, tTrue, npe0, distToReadout);
      uint32_t t1 = \
        getChannelTriggerTicks(engine, trigClock, tTrue, npe1, distToReadout);

      // Time relative to PPS: Random for now
      uint32_t ppsTicks = \
        CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

      // SiPM and ADC response: Npe to ADC counts
      short q0 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0, fQRMS * npe0);
      short q1 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe1, fQRMS * npe1);

      // Adjacent channels on a strip are numbered sequentially
      uint32_t moduleID = adsc.AuxDetID();
      uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
      uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;

      // Write AuxDetDigit for each channel
      crtHits->push_back(::crt::CRTData(channel0ID, t0, ppsTicks, q0));
      crtHits->push_back(::crt::CRTData(channel1ID, t1, ppsTicks, q1));
    }
  }

  e.put(std::move(crtHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt

