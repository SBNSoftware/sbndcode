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

  std::vector <int> fChannelsToIgnore;
};


sbnd::BasicLightFilter::BasicLightFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},  // ,
  // More initializers here.
  fPDTypes( p.get<std::vector<std::string>>("PDTypes",{"pmt_coated", "pmt_uncoated"}) ),
  fRawWaveformsModuleLabel( p.get< std::string >("RawWaveformsModuleLabel","pmtdecoder:PMTChannels")),
  fFlashADCThreshold( p.get<float>("FlashADCThreshold",-200)),
  fBeamWindowMin(p.get<float>("BeamWindowMin")),
  fBeamWindowMax(p.get<float>("BeamWindowMax")),
  fChannelsToIgnore( p.get<std::vector<int>>("ChannelsToIgnore",{}))
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

  std::vector<int> flash_tpc0(wvf_length, 0);
  std::vector<int> flash_tpc1(wvf_length, 0);
  int npmt_counter = 0;

  for(auto const& wvf : (*wvfHandle)) {
    int wvf_ch = wvf.ChannelNumber();
    if(std::find(fChannelsToIgnore.begin(), fChannelsToIgnore.end(), wvf_ch) != fChannelsToIgnore.end() ) continue;

    if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(wvf_ch) ) != fPDTypes.end() ){

      // make copy to get the median 
      std::vector<short int> wvf_copy;
      wvf_copy = wvf.Waveform();

      // get median of the waveform and remove it 
      const auto median_it = wvf_copy.begin() + int(wvf_copy.size() / 2);
      std::nth_element(wvf_copy.begin(), median_it , wvf_copy.end());
      auto median = *median_it;

      for(unsigned int i=0;i<wvf.size();i++){
        if (wvf_ch%2==0) flash_tpc0.at(i) += wvf.at(i) - median;
        else flash_tpc1.at(i) += wvf.at(i) - median;
      }
      npmt_counter++;
    }
  }

  auto flash_tpc0_peakval = *std::min_element(flash_tpc0.begin(), flash_tpc0.end());
  auto flash_tpc0_peaktime = std::distance(flash_tpc0.begin(), std::min_element(flash_tpc0.begin(), flash_tpc0.end()))*2 + wvf_start; // units of ns

  auto flash_tpc1_peakval = *std::min_element(flash_tpc1.begin(), flash_tpc1.end());
  auto flash_tpc1_peaktime = std::distance(flash_tpc1.begin(), std::min_element(flash_tpc1.begin(), flash_tpc1.end()))*2 + wvf_start; // units of ns

  std::cout << "TPC0 flash peak value = " << flash_tpc0_peakval << " PE at " << flash_tpc0_peaktime << " ns"<< std::endl;
  std::cout << "TPC1 flash peak value = " << flash_tpc1_peakval << " PE at " << flash_tpc1_peaktime << " ns"<< std::endl;

  bool tpc0_pass = false; 
  bool tpc1_pass = false;

  if ((flash_tpc0_peaktime > fBeamWindowMin) && (flash_tpc0_peaktime < fBeamWindowMax) && (flash_tpc0_peakval < fFlashADCThreshold)) tpc0_pass = true;
  if ((flash_tpc1_peaktime > fBeamWindowMin) && (flash_tpc1_peaktime < fBeamWindowMax) && (flash_tpc1_peakval < fFlashADCThreshold)) tpc1_pass = true;
  
  if (tpc0_pass | tpc1_pass){
    histname.str(std::string());
    int ntpc = (flash_tpc0_peakval < flash_tpc1_peakval) ? 0 : 1;
    histname << "event_" << e.id().event()
             << "_npmts_" << npmt_counter
             << "_tpc_" << ntpc;

    TH1D *flashHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", wvf_start), flash_tpc0.size(), wvf_start, wvf_end);
    for(unsigned int i = 0; i < flash_tpc0.size(); i++) {
      if (ntpc == 0) flashHist->SetBinContent(i + 1, (double)flash_tpc0[i]);
      else flashHist->SetBinContent(i + 1, (double)flash_tpc1[i]);
    }
  }

  return (tpc0_pass || tpc1_pass);
}

DEFINE_ART_MODULE(sbnd::BasicLightFilter)
