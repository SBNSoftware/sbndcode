////////////////////////////////////////////////////////////////////////
// Class:       pmtArtdaqFragmentProducer
// Plugin Type: producer (Unknown Unknown)
// File:        pmtArtdaqFragmentProducer_module.cc
//
// Generated at Mon Feb  7 14:18:34 2022 by Patrick Green using cetskelgen
// from  version .

// Module to convert simulated PMT waveforms into artdaq fragment format
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// artdaq includes
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// SBN/SBND includes
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

// ROOT includes
#include "TFile.h"
#include "TH1D.h"

// random numbers
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <ctime>

namespace sbnd {
  namespace trigger {
    class pmtArtdaqFragmentProducer;
  }
}


class sbnd::trigger::pmtArtdaqFragmentProducer : public art::EDProducer {
public:
  explicit pmtArtdaqFragmentProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pmtArtdaqFragmentProducer(pmtArtdaqFragmentProducer const&) = delete;
  pmtArtdaqFragmentProducer(pmtArtdaqFragmentProducer&&) = delete;
  pmtArtdaqFragmentProducer& operator=(pmtArtdaqFragmentProducer const&) = delete;
  pmtArtdaqFragmentProducer& operator=(pmtArtdaqFragmentProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  // fhicl parameters
  std::string fInputModuleNameWvfm;
  std::string fInputModuleNameTrigger;
  int fBaseline; // baseline in simulation, default 8000 ADC
  int fThreshold; // individual PMT threshold in ADC, pass if ADC less than threshold
  bool fVerbose;

  // event information
  int fRun, fSubrun;
  art::EventNumber_t fEvent;

  // PD information
  opdet::sbndPDMapAlg pdMap; // photon detector map
  std::vector<unsigned int> channelList; 

  // waveforms
  std::vector<std::vector<short>> wvf_channel;
  std::vector<std::vector<short>> wvf_bin_channel;

  // downsampled waveforms
  std::vector<std::vector<short>> wvf_channel_down;
  std::vector<std::vector<short>> wvf_bin_channel_down;  
   
  // sampling rate
  double fSampling;

  // services
  art::ServiceHandle<art::TFileService> tfs;

  // random numbers
  CLHEP::HepRandomEngine&       fTriggerTimeEngine;

};


sbnd::trigger::pmtArtdaqFragmentProducer::pmtArtdaqFragmentProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fInputModuleNameWvfm(p.get<std::string>("InputModuleNameWvfm")),
  fInputModuleNameTrigger(p.get<std::string>("InputModuleNameTrigger")),
  fBaseline(p.get<int>("Baseline",8000)),
  fThreshold(p.get<int>("Threshold",7900)),
  fVerbose(p.get<bool>("Verbose", false)),
  fTriggerTimeEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "trigger", p, "SeedTriggerTime"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<artdaq::Fragment> >();

  // get clock
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSampling = clockData.OpticalClock().Frequency(); // MHz

  // build PD map and channel list
  auto subsetCondition = [](auto const& i)->bool { return i["pd_type"] == "pmt_coated" || i["pd_type"] == "pmt_uncoated"; };
  auto pmtMap = pdMap.getCollectionFromCondition(subsetCondition);
  if (fVerbose) std::cout << "Number of PDs selected: \t" << pmtMap.size() << "\n";
  for(auto const& i:pmtMap){
    channelList.push_back(i["channel"]);
  }
   
}

void sbnd::trigger::pmtArtdaqFragmentProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  // event information
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  if (fVerbose) std::cout << "Processing Run: " << fRun << ", Subrun: " << fSubrun << ", Event: " << fEvent << std::endl;

  // access PMT waveforms and hardware trigger information
  art::Handle< std::vector< raw::OpDetWaveform > > wvfmHandle;
  art::Handle< std::vector< sbnd::comm::pmtTrigger > > triggerHandle;
  e.getByLabel(fInputModuleNameWvfm, wvfmHandle);
  e.getByLabel(fInputModuleNameTrigger, triggerHandle);

  if(!wvfmHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any waveforms, input must contain OpDetWaveforms." << "\n";
  }
  if(!triggerHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any PMT hardware trigger object, must run PMT trigger producer before running this module." << "\n";
  }

  // create empty vectors to hold waveforms for each channel
  std::vector<short> wvf_bin_0; wvf_bin_0.reserve((int)(3000/fSampling));
  std::vector<short> wvf_0; wvf_0.reserve((int)(3000/fSampling)); 
  for (double i = -1500.0; i<1500.0+(1./fSampling); i+=(1./fSampling)){
    wvf_bin_0.push_back(0);
    wvf_0.push_back(fBaseline);
  }
  for (size_t i = 0; i < channelList.size(); i++){
    wvf_bin_channel.push_back(wvf_bin_0);
    wvf_channel.push_back(wvf_0);
  }

  // counters
  int num_pmt_wvf = 0;
  
  // loop through waveform handles
  size_t wvf_id = -1;
  size_t hist_id = -1;
  for(auto const& wvf : (*wvfmHandle)) {
    wvf_id++;
    double fChNumber = wvf.ChannelNumber();
    std::string opdetType = pdMap.pdType(fChNumber);
    
    // only look at pmts
    if (opdetType != "pmt_coated" && opdetType != "pmt_uncoated") continue;
    
    num_pmt_wvf++;

    double fStartTime = wvf.TimeStamp(); // in us
    double fEndTime = double(wvf.size()) / fSampling + fStartTime; // in us

    //create binary waveform and full waveform
    std::vector<short> wvf_bin; wvf_bin.reserve((short)(3000/fSampling));
    std::vector<short> wvf_full; wvf_full.reserve((short)(3000/fSampling));

    if (fStartTime > -1500.0){
      for (double i = fStartTime+1500.0; i>0.; i-=(1./fSampling)){
        wvf_bin.push_back(0);
        wvf_full.push_back(fBaseline);
      }
    }

    for(unsigned int i = 0; i < wvf.size(); i++) {
      if((double)wvf[i]<fThreshold) wvf_bin.push_back(1);
      else  wvf_bin.push_back(0);

      wvf_full.push_back(wvf[i]);
    }

    if (fEndTime < 1500.0){
      for (double i = 1500.0-fEndTime; i>0.; i-=(1./fSampling)){
        wvf_bin.push_back(0);
        wvf_full.push_back(fBaseline);
      }
    }

    // combine waveform with any other waveforms from same channel
    int i_ch = -1.;
    auto ich = std::find(channelList.begin(), channelList.end(), fChNumber);
    if (ich != channelList.end()){
      i_ch = ich - channelList.begin();
    }
    if (wvf_bin_channel.at(i_ch).size() < wvf_bin.size()){
      if (fVerbose) std::cout<<"Binary waveform -- Previous Channel" << fChNumber <<" Size: "<<wvf_bin_channel.at(i_ch).size()<<"New Channel" << fChNumber <<" Size: "<<wvf_bin.size()<<std::endl;
      for(unsigned int i = wvf_bin_channel.at(i_ch).size(); i < wvf_bin.size(); i++) {
        wvf_bin_channel.at(i_ch).push_back(0);
      }
    }
    if (wvf_channel.at(i_ch).size() < wvf_full.size()){
      if (fVerbose) std::cout<<"Full waveform -- Previous Channel" << fChNumber <<" Size: "<<wvf_channel.at(i_ch).size()<<"New Channel" << fChNumber <<" Size: "<<wvf_full.size()<<std::endl;
      for(unsigned int i = wvf_channel.at(i_ch).size(); i < wvf_full.size(); i++) {
        wvf_channel.at(i_ch).push_back(fBaseline);
      }
    }
    for(unsigned int i = 0; i < wvf_bin.size(); i++) {
      if(wvf_bin_channel.at(i_ch).at(i) == 1 || wvf_bin[i] == 1) wvf_bin_channel.at(i_ch)[i] = 1;
      else wvf_bin_channel.at(i_ch)[i] = 0;
    }
    for(unsigned int i = 0; i < wvf_full.size(); i++) {
       wvf_channel.at(i_ch)[i] += (wvf_full[i] - fBaseline);
    } 

    hist_id++;
    //if (hist_id < 15) {
    // histogram for testing
    std::stringstream histname;
    histname << "event_" << fEvent
                 << "_opchannel_" << fChNumber
                 << "_" << opdetType
                 << "_" << hist_id
                 << "raw";

    TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), "Raw Waveform", wvf.size(), 0, (wvf.size()*2)/1000);
    wvfHist->GetXaxis()->SetTitle("t (#mus)");
    for(unsigned int i = 0; i < wvf.size(); i++) {
      wvfHist->SetBinContent(i+1, (double)wvf[i]);
    }
    //}
    
    wvf_bin.clear();
    wvfmHandle.clear();

  } // waveform handle loop

  // loop through waveforms applying factor 4 downsampling (simulate at 2ns, hardware at 8ns)
  // binary
  int wvf_num = -1;
  for (auto wvf_bin : wvf_bin_channel){
    wvf_num++;
    int fChNumber = channelList.at(wvf_num);
    double fStartTime = -1500.0;
    double fEndTime = 1500.0;

    // downscale binary waveform by 4
    std::vector<short> wvf_bin_down; wvf_bin_down.reserve((short)(3000/fSampling/4));
    for(unsigned int i = 0; i < wvf_bin.size(); i++) {
      if(i%4==0) wvf_bin_down.push_back(wvf_bin[i]);
    }

    // save
    wvf_bin_channel_down.push_back(wvf_bin_down);

    // plot for testing
    if (wvf_num < 1){
       std::stringstream histname_binary;
       histname_binary << "event_" << fEvent
                << "_opchannel_" << fChNumber
                << "_binary";
       TH1D *wvfbHist = tfs->make< TH1D >(histname_binary.str().c_str(), "Binary Waveform", wvf_bin.size(), fStartTime, fEndTime);
       wvfbHist->GetXaxis()->SetTitle("t (#mus)");
       std::cout << "size of binary wvf: " << wvf_bin.size() << std::endl;
       for(unsigned int i = 0; i < wvf_bin.size(); i++) {
         wvfbHist->SetBinContent(i + 1, wvf_bin[i]);
       }
     
       std::stringstream histname_binary_down;
       histname_binary_down << "event_" << fEvent
                << "_opchannel_" << fChNumber
                << "_binary_down";

       TH1D *wvfbdHist = tfs->make< TH1D >(histname_binary_down.str().c_str(), "Downsampled Binary Waveform", wvf_bin_down.size(), fStartTime, fEndTime);
       wvfbdHist->GetXaxis()->SetTitle("t (#mus)");
       for(unsigned int i = 0; i < wvf_bin_down.size(); i++) {
         wvfbdHist->SetBinContent(i + 1, wvf_bin_down[i]);
       }
    }
  }
  // full
  wvf_num = -1;
  for (auto wvf : wvf_channel){
    wvf_num++;
    int fChNumber = channelList.at(wvf_num);
    double fStartTime = -1500.0;
    double fEndTime = 1500.0;

    // downscale full waveform by 4
    std::vector<short> wvf_down; wvf_down.reserve((short)(3000/fSampling/4));
    for(unsigned int i = 0; i < wvf.size(); i++) {
      if(i%4==0) wvf_down.push_back(wvf[i]);
    }

    // save
    wvf_channel_down.push_back(wvf_down);

    // plot for testing
    if (wvf_num < 1){
       std::stringstream histname_full;
       histname_full << "event_" << fEvent
                << "_opchannel_" << fChNumber
                << "_full";
       TH1D *wvfHist = tfs->make< TH1D >(histname_full.str().c_str(), "Full Waveform", wvf.size(), fStartTime, fEndTime);
       wvfHist->GetXaxis()->SetTitle("t (#mus)");
       //std::cout << "wvf" << std::endl;
       for(unsigned int i = 0; i < wvf.size(); i++) {
         wvfHist->SetBinContent(i + 1, wvf[i]);
       }
     
       std::stringstream histname_full_down;
       histname_full_down << "event_" << fEvent
                << "_opchannel_" << fChNumber
                << "_full_down";

       TH1D *wvfdHist = tfs->make< TH1D >(histname_full_down.str().c_str(), "Downsampled Full Waveform", wvf_down.size(), fStartTime, fEndTime);
       wvfdHist->GetXaxis()->SetTitle("t (#mus)");
       //std::cout << "wvf_down" << std::endl;
       for(unsigned int i = 0; i < wvf_down.size(); i++) {
         wvfdHist->SetBinContent(i + 1, wvf_down[i]);
       }
    }
  }

  // access hardware trigger information
  std::vector<size_t> triggerIndex; 
  for(auto const& trigger : (*triggerHandle)) {
    for (size_t idx = 0; idx < trigger.numPassed.size(); idx++) {
      if (trigger.numPassed[idx] > 1) {
        // save index
        triggerIndex.push_back(idx);
        // skip ahead by 5120 samples (downsampled by 4) before checking again
        idx += 1280;
      }
    }
  }

  if (fVerbose) std::cout << "Number of PMT hardware triggers found: " << triggerIndex.size() << std::endl; 
  
  // fragments vector
  std::unique_ptr<std::vector<artdaq::Fragment>> vecFrag = std::make_unique<std::vector<artdaq::Fragment>>();

  // set properties of fragment that are common to event
  uint32_t nChannelsTotal = channelList.size();
  uint32_t nChannelsFrag = 15;
  uint32_t nFrag = (uint32_t)nChannelsTotal/nChannelsFrag;
  uint32_t wfm_length = 5120; // ~10us, 2ns tick

  // fragment properties
  uint32_t timestampVal = std::time(nullptr); // current time
  uint32_t sequenceIDVal = fEvent;
        
  // create and populate metadata
  sbndaq::CAENV1730FragmentMetadata metadata;
  metadata.nChannels = nChannelsFrag;
  metadata.nSamples = wfm_length;
  metadata.timeStampSec = timestampVal;
  metadata.timeStampNSec = timestampVal*1e9;

  // fragment handle properties
  uint32_t eventCounterVal = fEvent;
  uint32_t boardIDVal = 0;
  uint32_t triggerTimeTagVal = (uint32_t)CLHEP::RandFlat::shoot(&fTriggerTimeEngine, 0, 1e9/16);
   

  // loop over PMT hardware triggers
  for (auto wvfIdx : triggerIndex) {

    // index in downsampled waveform, 8ns tick
    //size_t trigIdx_down = wvfIdx;
    //size_t startIdx_down = trigIdx_down-125; // -1us
    //size_t endIdx_down = trigIdx_down+1125; // +9us

    // index in full waveform, 2ns tick
    size_t trigIdx = wvfIdx*4;
    size_t startIdx = trigIdx-500; // -1us
    //size_t endIdx = trigIdx+4500; // +9us

    // create fragments to hold waveforms, set properties and populate
    // 15 PMTs stored per fragment, 120/15 = 8 fragments per trigger
    for (size_t counter = 0; counter < nFrag; counter++) {
      
      // create fragment
      const auto fragment_datasize_bytes = metadata.ExpectedDataSize();
      uint32_t fragmentIDVal = counter;
      auto fragment_uptr = artdaq::Fragment::FragmentBytes(fragment_datasize_bytes, sequenceIDVal, fragmentIDVal, sbndaq::detail::FragmentType::CAENV1730, metadata); // unique pointer
      fragment_uptr->setTimestamp(timestampVal);

      // populate fragment header
      auto header_ptr = reinterpret_cast<sbndaq::CAENV1730EventHeader*>(fragment_uptr->dataBeginBytes());
        
      header_ptr->eventCounter = eventCounterVal;
      header_ptr->boardID = boardIDVal;
      header_ptr->triggerTimeTag = triggerTimeTagVal + trigIdx*2;  // ns // set timetag as random value for event, but encode trigger time within this event as offset
        
      // populate waveforms  
      // populate fragment with waveform
      uint16_t* data_begin = reinterpret_cast<uint16_t*>(fragment_uptr->dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
      uint16_t* value_ptr = data_begin;
      uint16_t value = 0;
      size_t ch_offset = 0;

      // loop over channels
      for (size_t i_ch = 0; i_ch < nChannelsFrag; i_ch++) {
        ch_offset = (size_t)(i_ch*wfm_length);
        // loop over waveform
        for (size_t i_t = 0; i_t < wfm_length; i_t++) {
          // set value
          value = wvf_channel[counter*nChannelsFrag + i_ch][startIdx+i_t];
          value_ptr = data_begin + ch_offset + i_t;
          *value_ptr = value;
        }
      }

      // add fragment to vector
      vecFrag->push_back(*fragment_uptr);  
    }
  }

  // plot for testing
  for (size_t i_frag = 0; i_frag < 8; i_frag++) {

  artdaq::Fragment frag = (*vecFrag)[i_frag];

  const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
                 + sizeof(sbndaq::CAENV1730EventHeader));

  const uint16_t* value_ptr =  data_begin;
  uint16_t value = 0;
  size_t ch_offset = 0;

  // loop over channels
  for (size_t i_ch = 0; i_ch < 15; i_ch++) {
    ch_offset = (size_t)(i_ch*wfm_length);
    std::vector<uint16_t> wvf_frag; wvf_frag.reserve(wfm_length);
  
    //--loop over waveform samples
    for(size_t i_t=0; i_t<wfm_length; ++i_t){ 
      //fTicksVec[i_t] = t0*Ttt_DownSamp + i_t;   /*timestamps, event level*/
      value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
      value = *(value_ptr);
      //std::cout << value << std::endl;
      wvf_frag.push_back(value);
    }

    std::stringstream histname_frag;
    histname_frag << "event_" << fEvent
            << "_opchannel_" << channelList[i_ch+i_frag*15]
            << "_full_fromFragment";

    TH1D *wvFragHist = tfs->make< TH1D >(histname_frag.str().c_str(), "Full Waveform from fragment", wvf_frag.size(), 0, 10*1e9);
    wvFragHist->GetXaxis()->SetTitle("t (ns)");
    for(unsigned int i = 0; i < wvf_frag.size(); i++) {
      wvFragHist->SetBinContent(i + 1, wvf_frag[i]);
    }
  
  }
  }

  // add fragments to event
  e.put(std::move(vecFrag));

  
}

DEFINE_ART_MODULE(sbnd::trigger::pmtArtdaqFragmentProducer)
