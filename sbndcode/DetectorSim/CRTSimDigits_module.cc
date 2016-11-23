///////////////////////////////////////////////////////////////////////////////
/// Class: CRTSimDigits
/// Module Type: producer
/// File: CRTSimDigits_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
///
/// Author: mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "sbndcode/Utilities/SignalShapingServiceT1053.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardata/RawData/AuxDetDigit.h"
#include "larsim/Simulation/sim.h"
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/AuxDetSimChannel.h"

#include "TRandom.h"
#include "TTree.h"

#include <cmath>
#include <memory>

class CRTSimDigits : public art::EDProducer {
public:
  explicit CRTSimDigits(fhicl::ParameterSet const & p);

  CRTSimDigits(CRTSimDigits const &) = delete;
  CRTSimDigits(CRTSimDigits &&) = delete;
  CRTSimDigits& operator = (CRTSimDigits const &) = delete;
  CRTSimDigits& operator = (CRTSimDigits &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;

  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

private:
  TTree* fTree;
  char* taggerName;
  int moduleID, stripID, channel0ID, channel1ID, trackID, npe;
  double energyDeposited, distToReadout;
  double entryX, entryY, entryZ, entryT, exitX, exitY, exitZ, exitT;
  double exitMomentumX, exitMomentumY, exitMomentumZ;
};

void CRTSimDigits::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}

void CRTSimDigits::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("crtana", "Simulated CRT");
  fTree->Branch("taggerName", &taggerName, "taggerName/C");
  fTree->Branch("moduleID", &moduleID, "moduleID/I");
  fTree->Branch("stripID", &stripID, "stripID/I");
  fTree->Branch("channel0ID", &channel0ID, "channel0ID/I");
  //fTree->Branch("channel1ID", &channel1ID, "channel1ID/I");
  fTree->Branch("trackID", &trackID, "trackID/I");
  fTree->Branch("energyDeposited", &energyDeposited, "energyDeposited/D");
  fTree->Branch("distToReadout", &distToReadout, "distToReadout/D");
  fTree->Branch("npe", &npe, "npe/I");
  fTree->Branch("entryX", &entryX, "entryX/D");
  fTree->Branch("entryY", &entryY, "entryY/D");
  fTree->Branch("entryZ", &entryZ, "entryZ/D");
  fTree->Branch("entryT", &entryT, "entryT/D");
  fTree->Branch("exitX", &exitX, "exitX/D");
  fTree->Branch("exitY", &exitY, "exitY/D");
  fTree->Branch("exitZ", &exitZ, "exitZ/D");
  fTree->Branch("exitT", &exitT, "exitT/D");
  fTree->Branch("exitMomentumX", &exitMomentumX, "exitMomentumX/D");
  fTree->Branch("exitMomentumY", &exitMomentumY, "exitMomentumY/D");
  fTree->Branch("exitMomentumZ", &exitMomentumZ, "exitMomentumZ/D");
}

CRTSimDigits::CRTSimDigits(fhicl::ParameterSet const & p) {
  this->reconfigure(p);
  produces<std::vector<raw::AuxDetDigit> >();
}

void CRTSimDigits::produce(art::Event & e) {
  std::unique_ptr<std::vector<raw::AuxDetDigit> > crtDigits(
      new std::vector<raw::AuxDetDigit>);

  art::ServiceHandle<geo::Geometry> adGeoService;

  art::Handle<std::vector<sim::AuxDetSimChannel> > adHandle;

  e.getByLabel(fG4ModuleLabel, adHandle);

  for (auto ad=adHandle->begin(); ad!=adHandle->end(); ++ad) {
    const geo::AuxDetGeo& adGeo = adGeoService->AuxDet(ad->AuxDetID());
    const geo::AuxDetSensitiveGeo& adSGeo = adGeo.SensitiveVolume(ad->AuxDetSensitiveID());

    taggerName = const_cast<char*>(adGeo.TotalVolume()->GetName());

    for (auto ide : ad->AuxDetIDEs()) {
      stripID = ad->AuxDetSensitiveID();
      channel0ID = stripID;

      // Fill position/energy data for the analysis tree
      entryX = ide.entryX;
      entryY = ide.entryY;
      entryZ = ide.entryZ;
      entryT = ide.entryT;
      exitX = ide.exitX;
      exitY = ide.exitY;
      exitZ = ide.exitZ;
      exitT = ide.exitT;
      energyDeposited = ide.energyDeposited;
      trackID = ide.trackID;
      fTree->Fill();

      // Simulate the CRT response
      // Distance to the readout end. Which way is "outside" depends on which
      // module was hit.
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      double world[3] = {x, y, z};

      double localSV[3];
      adSGeo.WorldToLocal(world, localSV);
      stripID = std::round(localSV[0] / 11.2 + 7.5);
      double adSCenterWorld[3];
      adSGeo.GetCenter(adSCenterWorld);
      double adSCenterLocal[3];
      adGeo.WorldToLocal(adSCenterWorld, adSCenterLocal);

      // Distance to the readout end ("outside") depends on module position
      if (adSCenterLocal[1] > 0) {
        distToReadout = abs(localSV[1] - adSGeo.HalfHeight());
      }
      else {
        distToReadout = abs(localSV[1] + adSGeo.HalfHeight());
      }

      // The expected number of PE for relativistic muons, based on an
      // exponential fit to the data in Figure 4.4 of SBND Document 685-v5
      int npeExpected = 43.62 * exp(-distToReadout/652.51);
      int npeTruth = gRandom->Poisson(npeExpected);

      // Gaussian smearing for resolution
      npe = gRandom->Gaus(npeTruth, 0.25);

      // Made-up "Digitizer" cf. SBND Document 685-v5, Figure 3.8
      short q = 100.0 * (1.0 * npe);

      // "Waveform" with one sample, that is total charge
      std::vector<short> adclist = {q};

      // Write AuxDetDigit for each channel (currently the same waveform)
      crtDigits->push_back(raw::AuxDetDigit(channel0ID, adclist, taggerName));
    }
  }

  e.put(std::move(crtDigits));
}

DEFINE_ART_MODULE(CRTSimDigits)

