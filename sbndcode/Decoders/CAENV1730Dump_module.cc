//-------------------------------------------------
//---------------------------------------

////////////////////////////////////////////////////////////////////////
// Class:       CAENV1730Dump
// Module Type: analyzer
// File:        CAENV1730Dump_module.cc
// Description: Prints out information about each event.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// #include "art/Framework/Principal/DataViewImpl.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

//#include "art/Framework/Services/Optional/TFileService.h" //before art_root_io transition
#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <bitset>

namespace sbndaq {
  class CAENV1730Dump;
}

/**************************************************************************************************/

class sbndaq::CAENV1730Dump : public art::EDAnalyzer {

public:
  struct Config {
    //--one atom for each parameter
    fhicl::Atom<art::InputTag> DataLabel {
      fhicl::Name("data_label"),
      fhicl::Comment("Tag for the input data product")
    };
    fhicl::Atom<bool> Verbose {
      fhicl::Name("verbose"), 
      fhicl::Comment("toggle for additional text output")
    };
  }; //--configuration
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CAENV1730Dump(Parameters const & pset);
  virtual ~CAENV1730Dump();

  void analyze(const art::Event& evt) override;
  void beginJob() override;
  void endJob() override;
  void analyze_caen_fragment(artdaq::Fragment & frag);


private:

  //--default values
  uint32_t nChannels;//    = 16;
  uint32_t Ttt_DownSamp;// =  4; 
 /* the trigger time resolution is 16ns when waveforms are sampled at
                               * 500MHz sampling. The trigger timestamp is
                               * sampled 4 times slower than input channels*/

  TNtuple* nt_header;
  
  TH1F*    hEventCounter;
  TH1F*    hTriggerTimeTag;
  TH1F*    h_wvfm_ev0_ch0;

  TTree* fEventTree;
  int fRun;
  art::EventNumber_t fEvent;
  int fragID;
  int seqID;
  std::vector<uint64_t>  fTicksVec;
  std::vector< std::vector<uint16_t> >  fWvfmsVec;
  
  bool firstEvt = true;

  // fcl params
  art::InputTag fDataLabel;
  bool fverbose;

}; //--class CAENV1730Dump


sbndaq::CAENV1730Dump::CAENV1730Dump(CAENV1730Dump::Parameters const& pset): art::EDAnalyzer(pset)
{
  fDataLabel = pset().DataLabel();
  fverbose = pset().Verbose();
}


void sbndaq::CAENV1730Dump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; 
  nt_header       = tfs->make<TNtuple>("nt_header","CAENV1730 Header Ntuple","art_ev:caen_ev:caenv_ev_tts"); 
  /************************************************************************************************/
  hEventCounter   = tfs->make<TH1F>("hEventCounter","Event Counter Histogram",10000,0,10000);
  hTriggerTimeTag = tfs->make<TH1F>("hTriggerTimeTag","Trigger Time Tag Histogram",10,2000000000,4500000000);
  h_wvfm_ev0_ch0  = tfs->make<TH1F>("h_wvfm_ev0_ch0","Waveform",2000,0,2000);  
  /************************************************************************************************/
  //--make tree to store the channel waveform info:
  fEventTree = tfs->make<TTree>("events","waveform tree");
  fEventTree->Branch("fRun",&fRun,"fRun/I");
  fEventTree->Branch("fEvent",&fEvent,"fEvent/I");
  fEventTree->Branch("fragID",&fragID,"fragID/I");
  fEventTree->Branch("seqID",&seqID,"seqID/I");
  fEventTree->Branch("fTicksVec",&fTicksVec);
  fEventTree->Branch("fWvfmsVec",&fWvfmsVec);
}

void sbndaq::CAENV1730Dump::endJob()
{
  std::cout << "Ending CAENV1730Dump...\n";
}


sbndaq::CAENV1730Dump::~CAENV1730Dump()
{
}


void sbndaq::CAENV1730Dump::analyze(const art::Event& evt)
{


  fRun = evt.run();
  fEvent = evt.event();
    
  /************************************************************************************************/


  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;

//   #if ART_HEX_VERSION < 0x30900
//           evt.getManyByType(fragmentHandles);
//   #else
    fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();
//   #endif

  /************************************************************************************************/
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      //Container fragment
      for (auto cont : *handle) {
	artdaq::ContainerFragment contf(cont);
	if (contf.fragment_type()==sbndaq::detail::FragmentType::CAENV1730) {
	  if (fverbose) 	  std::cout << "    Found " << contf.block_count() << " CAEN Fragments in container " << std::endl;
	  fWvfmsVec.resize(16*contf.block_count());
	  for (size_t ii = 0; ii < contf.block_count(); ++ii)
	  	analyze_caen_fragment(*contf[ii].get());
	}
      }
    }
    else {
      //normal fragment
      if (handle->front().type()==sbndaq::detail::FragmentType::CAENV1730) {
	if (fverbose) std::cout << "   found normal caen fragments " << handle->size() << std::endl;
	fWvfmsVec.resize(16*handle->size());
	for (auto frag : *handle)
	  analyze_caen_fragment(frag);
      }
    }
  } // loop over frag handles



}

void sbndaq::CAENV1730Dump::analyze_caen_fragment(artdaq::Fragment & frag)  {


    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    CAENV1730Event const* event_ptr = bb.Event();
    
    CAENV1730EventHeader header = event_ptr->Header;
    
    fragID = static_cast<int>(frag.fragmentID()); 
    seqID = static_cast<int>(frag.sequenceID()); 
    
    if (fverbose) std::cout << "\tFrom header, event counter is "  << header.eventCounter   << "\n";
    if (fverbose) std::cout << "\tFrom header, triggerTimeTag is " << header.triggerTimeTag << "\n";
    if (fverbose) std::cout << "\tFrom header, board id is "       << header.boardID        << "\n";
    if (fverbose) std::cout << "\tFrom fragment, fragment id is "  << fragID << "\n";
    if (fverbose) std::cout << "\tFrom fragment, sequence id is "  << seqID << "\n";
    if (fverbose) std::cout <<  "\tFrom fragment, timestamp is  " << frag.timestamp() << std::endl;
      
    uint32_t t0 = header.triggerTimeTag;
    hEventCounter->Fill(header.eventCounter);
    hTriggerTimeTag->Fill(t0);
    nt_header->Fill(fEvent,header.eventCounter,t0);
    nChannels = md->nChannels;
    std::cout << "\tNumber of channels: " << nChannels << "\n";
    
    //--get the number of 32-bit words (quad_bytes) from the header
    uint32_t ev_size_quad_bytes = header.eventSize;
    if (fverbose) std::cout << "Event size in quad bytes is: " << ev_size_quad_bytes << "\n";
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wfm_length = data_size_double_bytes/nChannels;
    if (fverbose) std::cout << "Channel waveform length = " << wfm_length << "\n";
    
    //--store the tick value for each acquisition 
    fTicksVec.resize(wfm_length);
    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
								   + sizeof(CAENV1730EventHeader));
    const uint16_t* value_ptr =  data_begin;
    uint16_t value = 0;
    size_t ch_offset = 0;
    //--loop over channels
    for (size_t i_ch=0; i_ch<nChannels; ++i_ch){
      fWvfmsVec[i_ch].resize(wfm_length);
      ch_offset = (size_t)(i_ch * wfm_length);
      //      std::cout << "ch" << i_ch << " offset =" << ch_offset << std::endl;
      
      //--loop over waveform samples
      for(size_t i_t=0; i_t<wfm_length; ++i_t){ 
	fTicksVec[i_t] = t0*Ttt_DownSamp + i_t;   /*timestamps, event level*/
	value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
	value = *(value_ptr);
	
	if (i_ch == 0 && firstEvt) { 
	  h_wvfm_ev0_ch0->SetBinContent(i_t,value);
	}
	
	fWvfmsVec[i_ch][i_t] = value;
	
      } //--end loop samples
      firstEvt = false;
      //      std::cout << " channel finished " << std::endl;
    } //--end loop channels
    
    fEventTree->Fill();
      
}

DEFINE_ART_MODULE(sbndaq::CAENV1730Dump)
