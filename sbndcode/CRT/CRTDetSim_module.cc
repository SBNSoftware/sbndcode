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
#include "larsim/RandomUtils/LArSeedService.h"

#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "lardata/DetectorInfo/ElecClock.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larsim/Simulation/AuxDetSimChannel.h"

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

  virtual ~CRTDetSim();

  CRTDetSim(CRTDetSim const &) = delete;
  CRTDetSim(CRTDetSim &&) = delete;
  CRTDetSim& operator = (CRTDetSim const &) = delete;
  CRTDetSim& operator = (CRTDetSim &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p) override;

  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

  TFile* ff;
  TNtuple* hta;
};


void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}


CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p) {
  this->reconfigure(p);

  art::ServiceHandle<sim::LArSeedService> seeds;
  seeds->createEngine(*this, "HepJamesRandom", "crt", p, "Seed");

  produces<std::vector<crt::CRTData> >();

  ff = new TFile("crt_temp.root", "recreate");
  hta = new TNtuple("hta", "", "dist:npe:ttrue:tmeas:tmean:trms:tdel:tprop");
}

CRTDetSim::~CRTDetSim() {
  ff->cd();
  hta->Write();
  ff->Close();
}

void CRTDetSim::produce(art::Event & e) {
  std::unique_ptr<std::vector<crt::CRTData> > crtHits(
      new std::vector<crt::CRTData>);

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;

  auto const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  detinfo::ElecClock trigClock = detClocks->TriggerClock();

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine("crt");

  // Handle for truth AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  // Loop through truth AD channels
  for (auto& adsc : *channels) {
    const geo::AuxDetGeo& adGeo = \
        geoService->AuxDet(adsc.AuxDetID());

    const geo::AuxDetSensitiveGeo& adsGeo = \
        adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

    std::string moduleName = adGeo.TotalVolume()->GetName();  // Strip array

    for (auto ide : adsc.AuxDetIDEs()) {
      // Adjacent channels on a strip are numbered sequentially
      uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID = 2 * stripID + 0;
      uint32_t channel1ID = 2 * stripID + 1;

      // Simulate the CRT response
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      double world[3] = {x, y, z};

      // Hit position in strip's local coordinates
      double svHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal);

      // Distance to the readout end ("outside") depends on module position
      // FOR NOW ASSUME ALL THE SAME DIRECTION
      double distToReadout = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);

      // The expected number of PE for relativistic muons, based on an
      // exponential fit to the data in Figure 4.4 of SBND Document 685-v5
      double npeExpected = CLHEP::RandFlat::shootInt(140); //43.62 * exp(-distToReadout / 652.51);

      // Put PE on channels weighted by 1/r^2
      double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]) / (2 * adsGeo.HalfWidth1());  // L
      double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]) / (2 * adsGeo.HalfWidth1());  // R
      double npeExp0 = npeExpected * (d0*d0) / (d0*d0 + d1*d1);
      double npeExp1 = npeExpected * (d1*d1) / (d0*d0 + d1*d1);
      
      // Observed PE, with Gaussian smearing
      double npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
      double npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);

      // Hit timing, with smearing and NPE dependence
      // Delay vs. NPE (terrible fit to Tech Note Figure 4.7)
      float p00 = 4125.74;
      float p01 = -300.31;
      float p02 = 90.392;
      double tDelayMean = p00*exp(-0.5*pow(((npe0+npe1)/2-p01)/p02, 2)) - 1.525;

      // Delay RMS vs. NPE (fit to Tech Note Figure 4.8)
      float p0 = 2.09138;
      float p1 = 7.23993;
      float p2 = 170.027;
      float p3 =  1.6544;
      float p4 = 75.6183;
      float p5 = 79.3543;
      double tDelayRMS = p0*exp(-pow((npe0+npe1)/2-p1,2)/p2) + p3*exp(-(x-p4)/p5);

      double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

      // Time resolution of the interpolator (cf. Figure 4.6)
      tDelay += CLHEP::RandGauss::shoot(engine, 0, 1.268);

      // Propagation time
      double tProp = CLHEP::RandGauss::shoot(0.061, 0.007) * distToReadout;

      double t = (ide.entryT + ide.exitT) / 2 + tProp + tDelay;
      trigClock.SetTime(t);
      uint32_t trigTicks = trigClock.Ticks();

      // Time relative to PPS: Random for now
      uint32_t ppsTicks = \
        CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

      std::cout << "HIT " << moduleName << "/"
                << adsGeo.TotalVolume()->GetName()
                << "(" << stripID << ")" << ": "
                << trigTicks << " " << ppsTicks << " " << distToReadout
                << std::endl;


      std::cout << "DIM "
                << "HalfHeight=" << adsGeo.HalfHeight() << ", "
                << "HalfWidth1=" << adsGeo.HalfWidth1() << ", "
                << "Length=" << adsGeo.Length() << std::endl;

      std::cout << "POS ("
                << svHitPosLocal[0] << ", "
                << svHitPosLocal[1] << ", "
                << svHitPosLocal[2]
                << ")" << std::endl;

      hta->Fill(distToReadout, npe0, (ide.entryT+ide.exitT)/2, t, tDelayMean, tDelayRMS, tDelay, tProp);

      // Made-up "Digitizer" cf. SBND Document 685-v5, Figure 3.8
      short q0 = 100.0 * (1.0 * npe0);
      short q1 = 100.0 * (1.0 * npe1);

      // Write AuxDetDigit for each channel (currently the same waveform)
      crtHits->push_back(::crt::CRTData(channel0ID, 21, 42, q0));
      crtHits->push_back(::crt::CRTData(channel1ID, 21, 42, q1));
    }
  }

  e.put(std::move(crtHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt

