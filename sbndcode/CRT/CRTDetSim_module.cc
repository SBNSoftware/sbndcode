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
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "CRTProducts/CRTData.hh"
#include "CRTDetSim.h"

#include <cmath>
#include <memory>
#include <string>

namespace sbnd {
namespace crt {

void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");

  fGlobalT0Offset = p.get<double>("GlobalT0Offset");
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
  fUseEdep = p.get<bool>("UseEdep");
  fQ0 = p.get<double>("Q0");
  fQPed = p.get<double>("QPed");
  fQSlope = p.get<double>("QSlope");
  fQRMS = p.get<double>("QRMS");
  fQThreshold = p.get<double>("QThreshold");
  fStripCoincidenceWindow = p.get<double>("StripCoincidenceWindow");
  fTaggerPlaneCoincidenceWindow = p.get<double>("TaggerPlaneCoincidenceWindow");
  fAbsLenEff = p.get<double>("AbsLenEff");
  fSipmTimeResponse = p.get<double>("SipmTimeResponse");
}


CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p)
  : EDProducer{p}
{
  art::ServiceHandle<rndm::NuRandomService> seeds;
  seeds->createEngine(*this, "HepJamesRandom", "crt", p, "Seed");

  this->reconfigure(p);

  produces<std::vector<sbnd::crt::CRTData> >();
}


uint32_t CRTDetSim::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
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
  clock.SetTime(t / 1e3);  // SetTime takes microseconds

  mf::LogInfo("CRT")
    << "CRT TIMING: t0=" << t0
    << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
    << ", tDelay=" << tDelay << ", tDelay(interp)="
    << tDelay << ", tProp=" << tProp << ", t=" << t << ", ticks=" << clock.Ticks() << "\n";

  return clock.Ticks();
}


struct Tagger {
  std::vector<std::pair<unsigned, uint32_t> > planesHit;
  std::vector<sbnd::crt::CRTData> data;
};


void CRTDetSim::produce(art::Event & e) {
  // A list of hit taggers, before any coincidence requirement
  std::map<std::string, Tagger> taggers;

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;

  art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
  detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(art::ScheduleID::first(),
                                                   moduleDescription().moduleLabel(),
                                                   "crt");

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  // Loop through truth AD channels
  for (auto& adsc : *channels) {
    const geo::AuxDetGeo& adGeo = \
        geoService->AuxDet(adsc.AuxDetID());

    const geo::AuxDetSensitiveGeo& adsGeo = \
        adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

    // Return the vector of IDEs
    std::vector<sim::AuxDetIDE> ides = adsc.AuxDetIDEs();
    std::sort(ides.begin(), ides.end(), 
              [](const sim::AuxDetIDE & a, const sim::AuxDetIDE & b) -> bool{
                return ((a.entryT + a.exitT)/2) < ((b.entryT + b.exitT)/2);
              });

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

    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());

    TGeoNode* nodeStrip = manager->GetCurrentNode();
    TGeoNode* nodeArray = manager->GetMother(1);
    TGeoNode* nodeModule = manager->GetMother(2);
    TGeoNode* nodeTagger = manager->GetMother(3);

    // Module position in parent (tagger) frame
    double origin[3] = {0, 0, 0};
    double modulePosMother[3];
    nodeModule->LocalToMaster(origin, modulePosMother);

    // Determine plane ID (1 for z > 0, 0 for z < 0 in local coordinates)
    unsigned planeID = (modulePosMother[2] > 0);

    // Determine module orientation: which way is the top (readout end)?
    bool top = (planeID == 1) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

    // Simulate the CRT response for each hit
    for (size_t ide_i = 0; ide_i < ides.size(); ide_i++) {
      
      sim::AuxDetIDE ide = ides[ide_i];

      int trackID = ide.trackID;

      // Finally, what is the distance from the hit (centroid of the entry
      // and exit points) to the readout end?
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;

      double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
      double eDep = ide.energyDeposited;
      double maxEdep = eDep;

      //ADD UP HITS AT THE SAME TIME - 2NS DIFF IS A GUESS -VERY APPROXIMATE
      if(ide_i < ides.size() - 1){
        while(ide_i < ides.size() - 1 && std::abs(tTrue-((ides[ide_i+1].entryT + ides[ide_i+1].exitT) / 2 + fGlobalT0Offset)) < fSipmTimeResponse){
          ide_i++;
          x = (x + ((ides[ide_i].entryX + ides[ide_i].exitX) / 2))/2;
          y = (y + ((ides[ide_i].entryY + ides[ide_i].exitY) / 2))/2;
          z = (z + ((ides[ide_i].entryZ + ides[ide_i].exitZ) / 2))/2;
          eDep += ides[ide_i].energyDeposited;
          tTrue = (tTrue + ((ides[ide_i].entryT + ides[ide_i].exitT) / 2))/2;
          if(ides[ide_i].energyDeposited > maxEdep && ides[ide_i].trackID > 0){
            trackID = ides[ide_i].trackID;
            maxEdep = ides[ide_i].energyDeposited;
          }
          if(trackID < 0 && (ides[ide_i].energyDeposited > maxEdep || ides[ide_i].trackID > 0)){
            trackID = ides[ide_i].trackID;
            maxEdep = ides[ide_i].energyDeposited;
          }
        }
      }

      double world[3] = {x, y, z};
      double svHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal);

      double distToReadout;
      if (top) {
        distToReadout = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);
      }
      else {
        distToReadout = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);
      }

      // The expected number of PE, using a quadratic model for the distance
      // dependence, and scaling linearly with deposited energy.
      double qr = fUseEdep ? 1.0 * eDep / fQ0 : 1.0;

      double npeExpected = \
        fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr;

      // Put PE on channels weighted by transverse distance across the strip,
      // using an exponential model
      double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
      double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
      double abs0 = exp(-d0 / fAbsLenEff);
      double abs1 = exp(-d1 / fAbsLenEff);
      double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
      double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

      // Observed PE (Poisson-fluctuated)
      long npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
      long npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);

      // Time relative to trigger, accounting for propagation delay and 'walk'
      // for the fixed-threshold discriminator
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

      // Adjacent channels on a strip are numbered sequentially.
      //
      // In the AuxDetChannelMapAlg methods, channels are identified by an
      // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
      // module, and a channel number from 0 to 32.
      uint32_t moduleID = adsc.AuxDetID();
      uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
      uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;

      // Apply ADC threshold and strip-level coincidence (both fibers fire)
      if (q0 > fQThreshold &&
          q1 > fQThreshold &&
          util::absDiff(t0, t1) < fStripCoincidenceWindow) {
        Tagger& tagger = taggers[nodeTagger->GetName()];
        tagger.planesHit.push_back({planeID, t0});
        tagger.data.push_back(sbnd::crt::CRTData(channel0ID, t0, ppsTicks, q0, trackID));
        tagger.data.push_back(sbnd::crt::CRTData(channel1ID, t1, ppsTicks, q1, trackID));
      }

      double poss[3];
      adsGeo.LocalToWorld(origin, poss);
      mf::LogInfo("CRT")
        << "CRT HIT in " << adsc.AuxDetID() << "/" << adsc.AuxDetSensitiveID() << "\n"
        << "CRT HIT POS " << x << " " << y << " " << z << "\n"
        << "CRT STRIP POS " << poss[0] << " " << poss[1] << " " << poss[2] << "\n"
        << "CRT MODULE POS " << modulePosMother[0] << " "
                             << modulePosMother[1] << " "
                             << modulePosMother[2] << " "
                             << "\n"
        << "CRT PATH: " << path << "\n"
        << "CRT level 0 (strip): " << nodeStrip->GetName() << "\n"
        << "CRT level 1 (array): " << nodeArray->GetName() << "\n"
        << "CRT level 2 (module): " << nodeModule->GetName() << "\n"
        << "CRT level 3 (tagger): " << nodeTagger->GetName() << "\n"
        << "CRT PLANE ID: " << planeID << "\n"
        << "CRT distToReadout: " << distToReadout << " " << (top ? "top" : "bot") << "\n"
        << "CRT q0: " << q0 << ", q1: " << q1 << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << "\n";
    }
  }

  // Apply coincidence trigger requirement
  std::unique_ptr<std::vector<sbnd::crt::CRTData> > triggeredCRTHits(
      new std::vector<sbnd::crt::CRTData>);

  // Logic: For normal taggers, require at least one hit in each perpendicular
  // plane. For the bottom tagger, any hit triggers read out.
  for (auto trg : taggers) {
    bool trigger = false;

    // Loop over pairs of hits
    for (auto t1 : trg.second.planesHit) {
      for (auto t2 : trg.second.planesHit) {
        if (t1 == t2) {
          continue;
        }

        // Two hits on different planes with proximal t0 times
        if (t1.first != t2.first && util::absDiff(t1.second, t2.second) < fTaggerPlaneCoincidenceWindow) {
          trigger = true;
          break;
        }
      }
    }

    if (trigger || trg.first.find("TaggerBot") != std::string::npos) {
      // Write out all hits on a tagger when there is any coincidence
      for (auto d : trg.second.data) {
        triggeredCRTHits->push_back(d);
      }
    }
  }

  mf::LogInfo("CRT") << "CRT TRIGGERED HITS: " << triggeredCRTHits->size() << "\n";

  e.put(std::move(triggeredCRTHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt
}  // namespace sbnd
