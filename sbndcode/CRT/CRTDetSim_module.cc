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
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbndcode/CRT/CRTDetSim.h"

#include <cmath>
#include <memory>
#include <string>

namespace sbnd {
namespace crt {

void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}


CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "crt", p, "Seed"))
  , fG4RefTime(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob().G4ToElecTime(0) * 1e3) // ns
  , fDetAlg(p.get<fhicl::ParameterSet>("DetSimParams"), fEngine, fG4RefTime)
{
  this->reconfigure(p);

  produces<std::vector<sbnd::crt::FEBData> >();
  produces<std::vector<sim::AuxDetIDE> >();
  produces<art::Assns<sbnd::crt::FEBData , sim::AuxDetIDE> >();

  consumes<std::vector<sim::AuxDetSimChannel>>(fG4ModuleLabel);
}



void CRTDetSim::produce(art::Event & e) {

  std::unique_ptr<std::vector<sbnd::crt::FEBData> > FEBDataOut(new std::vector<sbnd::crt::FEBData>);
  art::PtrMaker<sbnd::crt::FEBData> makeDataPtr(e);

  std::unique_ptr<std::vector<sim::AuxDetIDE> > auxDetIdes(new std::vector<sim::AuxDetIDE>);
  art::PtrMaker<sim::AuxDetIDE> makeIdePtr(e);

  std::unique_ptr<art::Assns<sbnd::crt::FEBData, sim::AuxDetIDE>> Dataassn
    (new art::Assns<sbnd::crt::FEBData, sim::AuxDetIDE>);

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  if (!channels.isValid()) {
    mf::LogWarning("CRTDetSim") << "No AuxDetSimChannel..." << std::endl;
  }

  fDetAlg.ClearTaggers();


  //
  // Step 1: Construct Taggers
  //
  for(auto const& adsc : *channels) {

    if(adsc.AuxDetID() == UINT_MAX || adsc.AuxDetSensitiveID() == UINT_MAX) {
      mf::LogWarning("CRTDetSim") << "AuxDetSimChannel with ID: UINT_MAX\n"
                                  << "skipping channel...";
      continue;
    }

    if(adsc.AuxDetIDEs().size() > 0) {
      mf::LogDebug("CRTDetSim") << "Adding AuxDetSimChannel with AuxDetID " << adsc.AuxDetID()
                                << ", AuxDetSensitiveID " << adsc.AuxDetSensitiveID() << std::endl;
      fDetAlg.FillTaggers(adsc.AuxDetID(), adsc.AuxDetSensitiveID(), adsc.AuxDetIDEs());
    }

  } // loop over AuxDetSimChannels


  //
  // Step 1: Apply Coincidence, deadtime, etc.
  //

  fDetAlg.CreateData();
  std::vector<std::pair<sbnd::crt::FEBData, std::vector<sim::AuxDetIDE>>> data = fDetAlg.GetData();

  //
  // Step 3: Save output
  //

  for(auto const& dataPair : data){

    FEBDataOut->push_back(dataPair.first);
    art::Ptr<FEBData> dataPtr = makeDataPtr(FEBDataOut->size()-1);

    for(auto const& ide : dataPair.second){
      auxDetIdes->push_back(ide);
      art::Ptr<sim::AuxDetIDE> idePtr = makeIdePtr(auxDetIdes->size()-1);
      Dataassn->addSingle(dataPtr, idePtr);
    }
  }

  mf::LogInfo("CRTDetSim") << "Number of FEBData objects produced: " << FEBDataOut->size() << "\n";

  e.put(std::move(FEBDataOut));
  e.put(std::move(auxDetIdes));
  e.put(std::move(Dataassn));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt
}  // namespace sbnd
