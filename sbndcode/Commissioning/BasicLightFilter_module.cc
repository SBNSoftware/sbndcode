////////////////////////////////////////////////////////////////////////
// Class:       BasicLightFilter
// Plugin Type: filter (Unknown Unknown)
// File:        BasicLightFilter_module.cc
//
// Generated at Sat Jun 29 12:28:05 2024 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include <algorithm>
#include <memory>

#include <TTree.h>
#include <string.h>
#include "TH1D.h"
#include "TFile.h"

namespace sbnd {
  class BasicLightFilter;
}


class sbnd::BasicLightFilter : public art::EDFilter {
public:
  explicit BasicLightFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BasicLightFilter(BasicLightFilter const&) = delete;
  BasicLightFilter(BasicLightFilter&&) = delete;
  BasicLightFilter& operator=(BasicLightFilter const&) = delete;
  BasicLightFilter& operator=(BasicLightFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

art::ServiceHandle<art::TFileService> tfs;

  std::vector <std::string> fPDTypes;
  std::string fRawWaveformsModuleLabel;
  opdet::sbndPDMapAlg fPDSMap;

  std::stringstream histname;

  float fFlashADCThreshold;
  float fBeamWindowMin;
  float fBeamWindowMax;
};


sbnd::BasicLightFilter::BasicLightFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},  // ,
  // More initializers here.
  fPDTypes( p.get<std::vector<std::string>>("PDTypes",{"pmt_coated", "pmt_uncoated"}) ),
  fRawWaveformsModuleLabel( p.get< std::string >("RawWaveformsModuleLabel","pmtdecoder:PMTChannels")),
  fFlashADCThreshold( p.get<float>("FlashADCThreshold",-200)),
  fBeamWindowMin(p.get<float>("BeamWindowMin")),
  fBeamWindowMax(p.get<float>("BeamWindowMax"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool sbnd::BasicLightFilter::filter(art::Event& e)
{

  art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
  e.getByLabel(fRawWaveformsModuleLabel, wvfHandle);
  if(!wvfHandle.isValid() || wvfHandle->size() == 0){
    std::cout << "RawWaveform with label " << fRawWaveformsModuleLabel << " not found..." << std::endl;
    return false;
  }

  // for obtaining the timestamp, need to make sure we dont get the timestamp/other info from channels in board 1  
  // this assumes there is only ONE flash trigger per event !!!!! 
  std::set<uint> board1_ch{39,38,11,10,9,8,63,62,37,36,7,6,61,60,64}; // this board has a touchy clock cable... ttt is not to be trusted

  size_t wvf_length=5000;
  double wvf_start=0;
  auto wvf_end = wvf_start + wvf_length*2; 

  for(auto const& wvf : (*wvfHandle)) {
    if (auto search = board1_ch.find(wvf.ChannelNumber()); search == board1_ch.end()){
      wvf_length = (*wvfHandle)[0].size();
      wvf_start  = (*wvfHandle)[0].TimeStamp();
      wvf_end    = wvf_start + wvf_length*2; 
      break;
    }
  }

  std::vector<int> flash(wvf_length, 0);
  int npmt_counter = 0;

  for(auto const& wvf : (*wvfHandle)) {
    int fChNumber = wvf.ChannelNumber();
    if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(fChNumber) ) != fPDTypes.end() ){

      // make copy to get the median 
      std::vector<short int> wvf_copy;
      wvf_copy = wvf.Waveform();

      // get median of the waveform and remove it 
      const auto median_it = wvf_copy.begin() + int(wvf_copy.size() / 2);
      std::nth_element(wvf_copy.begin(), median_it , wvf_copy.end());
      auto median = *median_it;

      for(unsigned int i=0;i<wvf.size();i++){
        flash.at(i) += wvf.at(i) - median;
      }
      npmt_counter++;
    }
  }

  auto flash_peakval = *std::min_element(flash.begin(), flash.end());
  auto flash_peaktime = std::distance(flash.begin(), std::min_element(flash.begin(), flash.end()))*2 + wvf_start; // units of ns

  std::cout << "flash peak value: " << flash_peakval << std::endl;
  std::cout << "flash peak time: " << flash_peaktime << std::endl;

  bool inbeamwindow = false;
  bool passthreshold = false;

  if ((flash_peaktime > fBeamWindowMin) && (flash_peaktime < fBeamWindowMax)) inbeamwindow = true;
  if (flash_peakval < fFlashADCThreshold) passthreshold = true;
  
  if (inbeamwindow && passthreshold){
    histname.str(std::string());
    histname << "event_" << e.id().event()
            << "_npmts_" << npmt_counter
              << "_PMTflash" ;
    TH1D *flashHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", wvf_start), flash.size(), wvf_start, wvf_end);
    for(unsigned int i = 0; i < flash.size(); i++) {
      flashHist->SetBinContent(i + 1, (double)flash[i]);
    }
  }

  return (inbeamwindow && passthreshold);
}

DEFINE_ART_MODULE(sbnd::BasicLightFilter)
