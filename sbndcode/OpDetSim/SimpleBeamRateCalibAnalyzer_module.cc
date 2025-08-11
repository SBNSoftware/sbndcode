////////////////////////////////////////////////////////////////////////
// Class:       SimpleBeamRateCalibAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        SimpleBeamRateCalibAnalyzer_module.cc
//
// Generated at Mon Aug 11 13:01:30 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

//copy includes from BeamRatesCalib_module, Then add BeamRateCalib service
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "BeamRateCalibService.h"

class SimpleBeamRateCalibAnalyzer;


class SimpleBeamRateCalibAnalyzer : public art::EDAnalyzer {
public:
  explicit SimpleBeamRateCalibAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleBeamRateCalibAnalyzer(SimpleBeamRateCalibAnalyzer const&) = delete;
  SimpleBeamRateCalibAnalyzer(SimpleBeamRateCalibAnalyzer&&) = delete;
  SimpleBeamRateCalibAnalyzer& operator=(SimpleBeamRateCalibAnalyzer const&) = delete;
  SimpleBeamRateCalibAnalyzer& operator=(SimpleBeamRateCalibAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  calib::BeamRateCalibService *fTriggerService;
  std::string fPMTName;
  // Declare member data here.

};


SimpleBeamRateCalibAnalyzer::SimpleBeamRateCalibAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fTriggerService = art::ServiceHandle<calib::BeamRateCalibService>(p)->ProviderFrom();
  fPMTName = p.get< std::string >("PMTName" );
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SimpleBeamRateCalibAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  art::ServiceHandle<art::TFileService> tfs;
  int fEvNumber = e.id().event();
  //grab waveforms to hand service
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fPMTName, waveHandle);
  int MonThreshold=50;
  std::vector<int> *MonPulse;
  bool Saving=false;
  int FlashCounter=0;
  //Save MonPulse waveform for flash zero
  fTriggerService->ConstructMonPulse(waveHandle, MonThreshold, MonPulse, Saving, FlashCounter);
  //mark with event id and flash info
  std::stringstream histname;
  histname << "event_" << fEvNumber <<"_Mon"<<"_"<<MonThreshold << "_"<<FlashCounter << "_TriggerPulse";
  TH1D *MonHist = tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(), 
                                      MonPulse->size(), 0.0, MonPulse->size()-1); //so this just breaks
  for(unsigned int i = 0; i < MonPulse->size(); i++) {
      MonHist->SetBinContent(i + 1, (double)(*MonPulse)[i]); //Loop over waveform and set bin content
    }
}

DEFINE_ART_MODULE(SimpleBeamRateCalibAnalyzer)
