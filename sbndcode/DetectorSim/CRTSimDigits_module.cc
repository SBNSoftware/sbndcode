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
};

DEFINE_ART_MODULE(CRTSimDigits)

void CRTSimDigits::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}

void CRTSimDigits::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
}

CRTSimDigits::CRTSimDigits(fhicl::ParameterSet const & p) {
  this->reconfigure(p);
  produces<std::vector<raw::AuxDetDigit> >();
}

void CRTSimDigits::produce(art::Event & e) {
  art::Handle<std::vector<sim::AuxDetSimChannel> > adHandle;

  e.getByLabel(fG4ModuleLabel, adHandle);

  std::unique_ptr<std::vector<raw::AuxDetDigit> > crtDigits(
      new std::vector<raw::AuxDetDigit>);

  if (!adHandle->empty()) {
    for (auto adIter=adHandle->begin(); adIter!=adHandle->end(); ++adIter) {
      std::vector<sim::AuxDetIDE> adSimIDEs = adIter->AuxDetIDEs();

      for (auto ideIt=adSimIDEs.begin(); ideIt!=adSimIDEs.end(); ++ideIt) {
        double crtT = gRandom->Gaus(ideIt->entryT, 2.0);
        double crtQ = gRandom->Gaus(ideIt->energyDeposited, 0.01);

        short q = 512.0 + 1000.0 * crtQ;  // "Digitizer"
        std::vector<short> adclist = {q};  // "Waveform"

        // Write AuxDetDigit for a channel number (currently AuxDetID)
        crtDigits->push_back(raw::AuxDetDigit(adIter->AuxDetID(), adclist));
      }
    }
  }

  e.put(std::move(crtDigits));
}

