////////////////////////////////////////////////////////////////////////////////
/// \file CRTDetSim_module.cc
///
/// \author mastbaum@uchicago.edu
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcore/Geometry/CryostatGeo.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "CRTData.hh"
#include "CRTDetSim.h"

#include <cmath>
#include <memory>
#include <string>

namespace crt {

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
  fQThreshold = p.get<double>("QThreshold");
  fStripCoincidenceWindow = p.get<double>("StripCoincidenceWindow");
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

  std::cout << "TIMING: t0=" << t0 << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS << ", tDelay=" << tDelay << ", tDelay(interp)=";

  // Time resolution of the interpolator
  tDelay += CLHEP::RandGauss::shoot(engine, 0, fTResInterpolator);

  std::cout << tDelay << ", tProp=";

  // Propagation time
  double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;

  std::cout << tProp << ", t=";

  double t = t0 + tProp + tDelay;

  std::cout << t << std::endl;

  // Get clock ticks
  clock.SetTime(t / 1e3);  // SetTime takes microseconds
  return clock.Ticks();
}


struct Tagger {
  std::set<unsigned> planesHit;
  std::vector<crt::CRTData> data;
};


void CRTDetSim::produce(art::Event & e) {
  // A list of hit taggers, before any coincidence requirement
  std::map<std::string, Tagger> taggers;

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;

  art::ServiceHandle<geo::AuxDetGeometry> adGeoService;
  const geo::AuxDetGeometry* adG = &(*adGeoService);
  const geo::AuxDetGeometryCore* adGeoCore = adG->GetProviderPtr();

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

      std::cout << "CRT HIT in " << adsc.AuxDetID() << "/" << adsc.AuxDetSensitiveID() << std::endl;
      std::cout << "POS " << x << " " << y << " " << z << std::endl;


      size_t adID, svID;
      const uint32_t chID = adGeoCore->PositionToAuxDetChannel(world, adID, svID);
      std::cout << "AD CHANNEL: " << chID << " AD: " << adID << " SV: " << svID << std::endl;

      // Find the path to the strip geo node, to locate it in the hierarchy
      std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
      std::vector<std::vector<TGeoNode const*> > paths = \
        geoService->FindAllVolumePaths(volNames);

      std::string path = "";
      for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
          path += "/";
        }
      }

      std::cout << "PATH: " << path << std::endl;
      TGeoManager* manager = geoService->ROOTGeoManager();
      manager->cd(path.c_str());

      TGeoNode* nodeStrip = manager->GetCurrentNode();
      std::cout << "level 0 (strip): " << nodeStrip->GetName() << std::endl;
      TGeoNode* nodeArray = manager->GetMother(1);
      std::cout << "level 1 (array): " << nodeArray->GetName() << std::endl;
      TGeoNode* nodeModule = manager->GetMother(2);
      std::cout << "level 2 (module): " << nodeModule->GetName() << std::endl;
      TGeoNode* nodeTagger = manager->GetMother(3);
      std::cout << "level 3 (tagger): " << nodeTagger->GetName() << std::endl;

      // Module position in parent (tagger) frame
      double origin[3] = {0, 0, 0};
      double modulePosMother[3];
      nodeModule->LocalToMaster(origin, modulePosMother);
      std::cout << "module pos " << modulePosMother[0] << " "
                                 << modulePosMother[1] << " "
                                 << modulePosMother[2] << " "
                                 << std::endl;

      double poss[3];
      adsGeo.LocalToWorld(origin, poss);
      std::cout << "STRIP POS " << poss[0] << " " << poss[1] << " " << poss[2] << std::endl;

      // Determine plane ID (0 for z > 0, 1 for z < 0 in local coordinates)
      unsigned planeID = (modulePosMother[2] > 0);
      std::cout << "planeID = " << planeID << std::endl;

      // Determine module orientation: which way is the top (readout end)?
      bool top = (planeID == 0) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

      // Finally, what is the distance from the hit to the readout end?
      double distToReadout;
      if (top) {
        distToReadout = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);
      }
      else {
        distToReadout = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);
      }

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

      // Time relative to PPS: Random for now! (FIXME)
      uint32_t ppsTicks = \
        CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

      // SiPM and ADC response: Npe to ADC counts
      short q0 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
      short q1 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));

      // Adjacent channels on a strip are numbered sequentially
      uint32_t moduleID = adsc.AuxDetID();
      uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
      uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;

      std::string name = geoService->AuxDet(moduleID).TotalVolume()->GetName();
      TVector3 pos = adGeoCore->AuxDetChannelToPosition(2*stripID, name);
      std::cout << "AD POS: " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;

      // Apply ADC threshold and strip-level coincidence
      std::cout << "q0: " << q0 << ", q1: " << q1 << ", dt: " << abs(t0-t1) << std::endl;
      if (q0 > fQThreshold && q1 > fQThreshold && abs(t0 - t1) < fStripCoincidenceWindow) {
        Tagger& tagger = taggers[nodeTagger->GetName()];
        tagger.planesHit.insert(planeID);
        tagger.data.push_back(crt::CRTData(channel0ID, t0, ppsTicks, q0));
        tagger.data.push_back(crt::CRTData(channel1ID, t1, ppsTicks, q1));
      }
    }
  }

  // Apply coincidence trigger requirement
  //
  // Logic: For normal taggers, require at least one hit in each perpendicular
  // plane. For the bottom tagger, all hits are read out.
  std::unique_ptr<std::vector<crt::CRTData> > triggeredCRTHits(
      new std::vector<crt::CRTData>);

  for (auto trg : taggers) {
    std::cout << trg.first << ": total " << trg.second.planesHit.size() << std::endl;
    for (auto pl : trg.second.planesHit) {
      std::cout << " " << pl;
    }
    std::cout << std::endl;

    if (trg.first.find("TaggerBot") != std::string::npos || trg.second.planesHit.size() > 1) {
      for (auto d : trg.second.data) {
        triggeredCRTHits->push_back(d);
      }
    }
  }

  std::cout << "ndatas: " << triggeredCRTHits->size() << std::endl;
  e.put(std::move(triggeredCRTHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt

