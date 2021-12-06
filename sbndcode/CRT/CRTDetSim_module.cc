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
#include "nurandom/RandomUtils/NuRandomService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"

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
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbndcode/CRT/CRTDetSim.h"

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
  fClockSpeedCRT = p.get<double>("ClockSpeedCRT");
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
  fAdcSaturation = p.get<uint32_t>("AdcSaturation");
}


CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "crt", p, "Seed"))
{
  this->reconfigure(p);

  produces<std::vector<sbnd::crt::CRTData> >();
  produces<std::vector<sim::AuxDetIDE> >();
  produces< art::Assns<sbnd::crt::CRTData , sim::AuxDetIDE> >();
}


uint32_t CRTDetSim::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                         /*detinfo::ElecClock& clock,*/
                                         float t0, float npeMean, float r) {
  // Hit timing, with smearing and NPE dependence
  double tDelayMean =
    fTDelayNorm *
      exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
    fTDelayOffset;

  double tDelayRMS =
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
  // FIXME no clock available for CRTs, have to do it by hand
  //clock.SetTime(t / 1e3);  // SetTime takes microseconds
  int time = (t / 1e3) * fClockSpeedCRT;

  mf::LogInfo("CRT")
    << "CRT TIMING: t0=" << t0
    << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
    << ", tDelay=" << tDelay << ", tDelay(interp)="
    << tDelay << ", tProp=" << tProp << ", t=" << t << /*", ticks=" << clock.Ticks() <<*/ "\n";

  return time;//clock.Ticks();
}


struct Tagger {
  std::vector<std::pair<unsigned, uint32_t> > planesHit;
  std::vector<sbnd::crt::CRTData> data;
  std::vector<std::vector<sim::AuxDetIDE>> ides;
};


void CRTDetSim::produce(art::Event & e) {
  // A list of hit taggers, before any coincidence requirement
  std::map<std::string, Tagger> taggers;

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;

  /*art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
  detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();*/

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  if (!channels.isValid()) {std::cout << "No Channels!" << std::endl;}

  // Loop through truth AD channels
  for (auto& adsc : *channels) {
    const geo::AuxDetGeo& adGeo =
        geoService->AuxDet(adsc.AuxDetID());
    const geo::AuxDetSensitiveGeo& adsGeo =
        adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

    // Return the vector of IDEs
    std::vector<sim::AuxDetIDE> ides = adsc.AuxDetIDEs();
    std::sort(ides.begin(), ides.end(),
              [](const sim::AuxDetIDE & a, const sim::AuxDetIDE & b) -> bool{
                return ((a.entryT + a.exitT)/2) < ((b.entryT + b.exitT)/2);
              });

    // std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
    std::set<std::string> volNames = { adGeo.TotalVolume()->GetName() };
    std::vector<std::vector<TGeoNode const*> > paths =
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

    // We get the array of strips first, which is the AuxDet,
    // then from the AuxDet, we get the strip by picking the
    // daughter with the ID of the AuxDetSensitive, and finally
    // from the AuxDet, we go up and pick the module and tagger
    TGeoNode* nodeArray = manager->GetCurrentNode();
    TGeoNode* nodeStrip = nodeArray->GetDaughter(adsc.AuxDetSensitiveID());
    TGeoNode* nodeModule = manager->GetMother(1);
    TGeoNode* nodeTagger = manager->GetMother(2);

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

      // Finally, what is the distance from the hit (centroid of the entry
      // and exit points) to the readout end?
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      int nides = 1;

      double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
      double tTrueLast = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
      double eDep = ide.energyDeposited;

      std::vector<sim::AuxDetIDE> trueIdes;
      trueIdes.push_back(ide);

      //ADD UP HITS AT THE SAME TIME - FIXME 2NS DIFF IS A GUESS -VERY APPROXIMATE
      if(ide_i < ides.size() - 1){
        while(ide_i < ides.size() - 1 && std::abs(tTrueLast-((ides[ide_i+1].entryT + ides[ide_i+1].exitT) / 2 + fGlobalT0Offset)) < fSipmTimeResponse){
          ide_i++;
          x += (ides[ide_i].entryX + ides[ide_i].exitX) / 2;
          y += (ides[ide_i].entryY + ides[ide_i].exitY) / 2;
          z += (ides[ide_i].entryZ + ides[ide_i].exitZ) / 2;
          eDep += ides[ide_i].energyDeposited;
          tTrue += (ides[ide_i].entryT + ides[ide_i].exitT) / 2;
          tTrueLast = (ides[ide_i].entryT + ides[ide_i].exitT) / 2;

          nides++;

          trueIdes.push_back(ides[ide_i]);
        }
      }

      x = x/nides;
      y = y/nides;
      z = z/nides;
      tTrue = tTrue/nides;

      double world[3] = {x, y, z};
      double svHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal);

      double distToReadout;
      if (top) {
        distToReadout = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);
      }
      else {
        distToReadout = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);
      }

      // The expected number of PE, using a quadratic model for the distance
      // dependence, and scaling linearly with deposited energy.
      double qr = fUseEdep ? 1.0 * eDep / fQ0 : 1.0;

      double npeExpected =
        fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr;

      // Put PE on channels weighted by transverse distance across the strip,
      // using an exponential model
      double d0 = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);  // L
      double d1 = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);  // R
      double abs0 = exp(-d0 / fAbsLenEff);
      double abs1 = exp(-d1 / fAbsLenEff);
      double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
      double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

      // Observed PE (Poisson-fluctuated)
      long npe0 = CLHEP::RandPoisson::shoot(&fEngine, npeExp0);
      long npe1 = CLHEP::RandPoisson::shoot(&fEngine, npeExp1);

      // Time relative to trigger, accounting for propagation delay and 'walk'
      // for the fixed-threshold discriminator
      uint32_t t0 =
        getChannelTriggerTicks(&fEngine, /*trigClock,*/ tTrue, npe0, distToReadout);
      uint32_t t1 =
        getChannelTriggerTicks(&fEngine, /*trigClock,*/ tTrue, npe1, distToReadout);

      // Time relative to PPS: Random for now! (FIXME)
      uint32_t ppsTicks =
        CLHEP::RandFlat::shootInt(&fEngine, /*trigClock.Frequency()*/ fClockSpeedCRT * 1e6);

      // SiPM and ADC response: Npe to ADC counts
      uint32_t q0 =
        CLHEP::RandGauss::shoot(&fEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
      if(q0 > fAdcSaturation) q0 = fAdcSaturation;
      uint32_t q1 =
        CLHEP::RandGauss::shoot(&fEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
      if(q1 > fAdcSaturation) q1 = fAdcSaturation;

      // Adjacent channels on a strip are numbered sequentially.
      //
      // In the AuxDetChannelMapAlg methods, channels are identified by an
      // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
      // module, and a channel number from 0 to 32.
      uint32_t moduleID = adsc.AuxDetID();
      uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
      uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;

      if (moduleID >= 127) {continue;}        //Ignoring MINOS modules for now.

      // Apply ADC threshold and strip-level coincidence (both fibers fire)
      if (q0 > fQThreshold &&
          q1 > fQThreshold &&
          util::absDiff(t0, t1) < fStripCoincidenceWindow) {
        Tagger& tagger = taggers[nodeTagger->GetName()];
        tagger.planesHit.push_back({planeID, t0});
        tagger.data.push_back(sbnd::crt::CRTData(channel0ID, t0, ppsTicks, q0));
        tagger.ides.push_back(trueIdes);
        tagger.data.push_back(sbnd::crt::CRTData(channel1ID, t1, ppsTicks, q1));
        tagger.ides.push_back(trueIdes);
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
  art::PtrMaker<sbnd::crt::CRTData> makeDataPtr(e);

  std::unique_ptr<std::vector<sim::AuxDetIDE> > auxDetIdes(
      new std::vector<sim::AuxDetIDE>);
  art::PtrMaker<sim::AuxDetIDE> makeIdePtr(e);

  std::unique_ptr< art::Assns<sbnd::crt::CRTData, sim::AuxDetIDE> > Dataassn( new art::Assns<sbnd::crt::CRTData, sim::AuxDetIDE>);

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
      // Write out all hits on a tagger when there is any coincidence FIXME this reads out everything!
      //for (auto d : trg.second.data) {
      for (size_t d_i = 0; d_i < trg.second.data.size(); d_i++) {
        triggeredCRTHits->push_back(trg.second.data[d_i]);
        art::Ptr<sbnd::crt::CRTData> dataPtr = makeDataPtr(triggeredCRTHits->size()-1);
        if(trg.second.data.size() == trg.second.ides.size()){
          for (size_t i_i = 0; i_i < trg.second.ides[d_i].size(); i_i++){
            auxDetIdes->push_back(trg.second.ides[d_i][i_i]);
            art::Ptr<sim::AuxDetIDE> idePtr = makeIdePtr(auxDetIdes->size()-1);
            Dataassn->addSingle(dataPtr, idePtr);
          }
        }
      }
    }
  }

  mf::LogInfo("CRT") << "CRT TRIGGERED HITS: " << triggeredCRTHits->size() << "\n";

  e.put(std::move(triggeredCRTHits));
  e.put(std::move(auxDetIdes));
  e.put(std::move(Dataassn));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt
}  // namespace sbnd
