////////////////////////////////////////////////////////////////////////
// Class:       FlashAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        FlashAnalyzer_module.cc
//
// Generated at Fri Jun  7 18:29:11 2024 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// ROOT and C++ includes
#include <TTree.h>
#include <string.h>
#include "TH1D.h"
#include "TFile.h"


namespace opdet {
  class FlashAnalyzer;
}


class opdet::FlashAnalyzer : public art::EDAnalyzer {
public:
  explicit FlashAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashAnalyzer(FlashAnalyzer const&) = delete;
  FlashAnalyzer(FlashAnalyzer&&) = delete;
  FlashAnalyzer& operator=(FlashAnalyzer const&) = delete;
  FlashAnalyzer& operator=(FlashAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  std::vector <std::string> fPDTypes;
  std::string fRawWaveformsModuleLabel;

  std::stringstream histname;

  opdet::sbndPDMapAlg fPDSMap;

};


opdet::FlashAnalyzer::FlashAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fPDTypes( p.get<std::vector<std::string>>("PDTypes") ),
  fRawWaveformsModuleLabel( p.get< std::string >("RawWaveformsModuleLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module. 
}

void opdet::FlashAnalyzer::analyze(art::Event const& e)
{
// Implementation of required member function here.
  art::ServiceHandle<art::TFileService> tfs;

 
  art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
  e.getByLabel(fRawWaveformsModuleLabel, wvfHandle);
  if(!wvfHandle.isValid()){
    std::cout << "RawWaveform with label " << fRawWaveformsModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  auto wvf_length = (*wvfHandle)[0].size();
  auto wvf_start  = (*wvfHandle)[0].TimeStamp();
  auto wvf_end    = wvf_start + wvf_length*2; 

  std::vector<double> flash(wvf_length, 0);
  
  for(auto const& wvf : (*wvfHandle)) {
    int fChNumber = wvf.ChannelNumber();
    if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(fChNumber) ) != fPDTypes.end() ){
      for(unsigned int i=0;i<wvf.size();i++){
        flash.at(i) += wvf.at(i);
      }
    }
  }

  histname.str(std::string());
  histname << "event_" << e.id().event()
            << "_PMTflash" ;
  TH1D *flashHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", wvf_start), flash.size(), wvf_start, wvf_end);
  for(unsigned int i = 0; i < flash.size(); i++) {
    flashHist->SetBinContent(i + 1, (double)flash[i]);
  }
}

DEFINE_ART_MODULE(opdet::FlashAnalyzer)
