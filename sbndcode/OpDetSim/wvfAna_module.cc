////////////////////////////////////////////////////////////////////////
// Class:       wvfAna
// Module Type: analyzer
// File:        wvfAna_module.cc
//
// Analyzer to read optical waveforms
//
// Authors: L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

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

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet {

  class wvfAna;

  class wvfAna : public art::EDAnalyzer {
  public:
    explicit wvfAna(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    wvfAna(wvfAna const &) = delete;
    wvfAna(wvfAna &&) = delete;
    wvfAna & operator = (wvfAna const &) = delete;
    wvfAna & operator = (wvfAna &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    //Selected optional functions
    void beginJob() override;
    void endJob() override;

    opdet::sbndPDMapAlg pdMap; //map for photon detector types

  private:

    size_t fEvNumber;
    size_t fChNumber;
    double fSampling;
    double fSampling_Daphne;
    double fStartTime;
    double fEndTime;
    //TTree *fWaveformTree;

    // Declare member data here.
    std::string fInputModuleName;
    std::vector<std::string> fOpDetsToPlot;
    std::stringstream histname;
    std::string opdetType;
    std::string opdetElectronics;
  };


  wvfAna::wvfAna(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p)  // ,
    // More initializers here.
  {
    fInputModuleName = p.get< std::string >("InputModule" );
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampling = clockData.OpticalClock().Frequency(); // MHz
    fSampling_Daphne = p.get<double>("DaphneFrequency" );
  }

  void wvfAna::beginJob()
  {

  }

  void wvfAna::analyze(art::Event const & e)
  {
    // Implementation of required member function here.
    std::cout << "My module on event #" << e.id().event() << std::endl;

    art::ServiceHandle<art::TFileService> tfs;
    fEvNumber = e.id().event();

    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
    e.getByLabel(fInputModuleName, waveHandle);

    if(!waveHandle.isValid()) {
      std::cout << Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
    }

    // // example of usage for pdMap.getCollectionWithProperty()
    // //
    // // define a container
    // auto inBoxTwo = pdMap.getCollectionWithProperty("pds_box", 2);
    // // you can cout the whole json object
    // std::cout << "inBoxTwo:\t" << inBoxTwo << "\n";
    // // traverse its components in a loop
    // for (auto const &e: inBoxTwo) {
    //   std::cout << e["pd_type"] << " " << e["channel"] << ' ' << "\n";
    // }

    // // example of usage for pdMap.getCollectionFromCondition()
    // // define a lambda function with the conditions
    // auto subsetCondition = [](auto const& e)->bool
    //   // modify conditions as you want in the curly braces below
    //   {return e["pd_type"] == "pmt_uncoated" && e["tpc"] == 0;};
    // // get the container that satisfies the conditions
    // auto uncoatedsInTPC0 = pdMap.getCollectionFromCondition(subsetCondition);
    // std::cout << "uncoatedsInTPC0.size():\t" << uncoatedsInTPC0.size() << "\n";
    // for(auto const& e:uncoatedsInTPC0){
    //   std::cout << "e:\t" << e << "\n";
    // }

    std::cout << "Number of waveforms: " << waveHandle->size() << std::endl;

    std::cout << "fOpDetsToPlot:\t";
    for (auto const& opdet : fOpDetsToPlot){std::cout << opdet << " ";}
    std::cout << std::endl;

    int hist_id = 0;
    for(auto const& wvf : (*waveHandle)) {
      fChNumber = wvf.ChannelNumber();
      opdetType = pdMap.pdType(fChNumber);
      opdetElectronics = pdMap.electronicsType(fChNumber);
      if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
      histname.str(std::string());
      histname << "event_" << fEvNumber
               << "_opchannel_" << fChNumber
               << "_" << opdetType
               << "_" << hist_id;

      fStartTime = wvf.TimeStamp(); //in us
      if (opdetElectronics == "daphne"){
				fEndTime = double(wvf.size()) / fSampling_Daphne + fStartTime;
			} //in us
			else{
        fEndTime = double(wvf.size()) / fSampling + fStartTime;
      } //in us

      //Create a new histogram
      TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", fStartTime), wvf.size(), fStartTime, fEndTime);
      for(unsigned int i = 0; i < wvf.size(); i++) {
        wvfHist->SetBinContent(i + 1, (double)wvf[i]);
      }
      hist_id++;
    }
  }

  void wvfAna::endJob()
  {
  }

  DEFINE_ART_MODULE(opdet::wvfAna)

}
