////////////////////////////////////////////////////////////////////////
// Class:       SBNDPTBDecoder
// Plugin Type: producer 
// File:        SBNDPTBDecoder_module.cc
//
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

#include <memory>
#include <bitset>

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"

class SBNDPTBDecoder;


class SBNDPTBDecoder : public art::EDProducer {
public:
  explicit SBNDPTBDecoder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDPTBDecoder(SBNDPTBDecoder const &) = delete;
  SBNDPTBDecoder(SBNDPTBDecoder &&) = delete;
  SBNDPTBDecoder & operator = (SBNDPTBDecoder const &) = delete;
  SBNDPTBDecoder & operator = (SBNDPTBDecoder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // Declare member data here.

  std::string fInputLabel;
  std::string fInputContainerInstance;
  std::string fInputNonContainerInstance;
  std::string fOutputInstance;
  int fDebugLevel;
  
  typedef struct ptbsv
  {
    std::vector<raw::ptb::Trigger> HLTrigs;
    std::vector<raw::ptb::Trigger> LLTrigs;
    std::vector<raw::ptb::ChStatus> ChStats;
    std::vector<raw::ptb::Feedback> Feedbacks;
    std::vector<raw::ptb::Misc> Miscs;
    std::vector<raw::ptb::WordIndex> WordIndexes; 
  } ptbsv_t;
  
  void _process_PTB_AUX(const artdaq::Fragment& frag, ptbsv_t &sout);
};


SBNDPTBDecoder::SBNDPTBDecoder(fhicl::ParameterSet const & p)
  : EDProducer{p}
    // Initialize member data here.
{
  fInputLabel = p.get<std::string>("InputLabel");
  fInputContainerInstance = p.get<std::string>("InputContainerInstance");
  fInputNonContainerInstance = p.get<std::string>("InputNonContainerInstance");
  fOutputInstance = p.get<std::string>("OutputInstance");
  fDebugLevel = p.get<int>("DebugLevel",0);
  
  produces<std::vector<raw::ptb::sbndptb> >(fOutputInstance);
}


int upBitLocation (int val){
  int upBit =-1;
  for(int i=0; i<32;i++){
    if ((val & 0x01) ==1){
      upBit =i;
    }
    val = val >> 1;
  }
  return upBit;
}

int eventNum =0;
int _run;
int _subrun;
int _event;

void SBNDPTBDecoder::produce(art::Event & evt)
{
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();

  std::cout << "Run: " << _run  << "  SubRun: " << _subrun  << "  Event: " << _event << " Listed Event: " << eventNum << std::endl; 
  eventNum++;

  // look first for container fragments and then non-container fragments

  std::vector<raw::ptb::sbndptb> sbndptbs;

  art::InputTag itag1(fInputLabel, fInputContainerInstance);
  auto cont_frags = evt.getHandle<artdaq::Fragments>(itag1);
  int numcont =0;
  if (cont_frags)
    {
      for (auto const& cont : *cont_frags)
	{
	  artdaq::ContainerFragment cont_frag(cont);
	  numcont++;
	  int numfrag =0;
	  for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	    {
	      numfrag++;
              ptbsv_t sout;  // output structures
	      _process_PTB_AUX(*cont_frag[ii], sout);
              raw::ptb::sbndptb ptbdp(sout.HLTrigs,sout.LLTrigs,sout.ChStats,sout.Feedbacks,sout.Miscs,sout.WordIndexes);
	      sbndptbs.push_back(ptbdp);
	    }

	}
    }

  art::InputTag itag2(fInputLabel, fInputNonContainerInstance);
  auto frags = evt.getHandle<artdaq::Fragments>(itag2);
  if (frags)
    {
      for(auto const& frag: *frags)
	{
          ptbsv_t sout;  // output structures
	  _process_PTB_AUX(frag, sout);
          raw::ptb::sbndptb ptbdp(sout.HLTrigs,sout.LLTrigs,sout.ChStats,sout.Feedbacks,sout.Miscs,sout.WordIndexes);
          sbndptbs.push_back(ptbdp);
	}
    }

  evt.put(std::make_unique<std::vector<raw::ptb::sbndptb>>(std::move(sbndptbs)),fOutputInstance);
}






//Copied from Lines 1382-1406 from     sbndaq-artdaq/sbndaq-artdaq/ArtModules/SBND/EventAna_module.cc

 /*******************************************
  * Note: Below for the Timestamp conversion,
  *       The PTB TS is in UTC seconds since the Unix epoch in 50MHz clock ticks.
  *       This means to recover seconds since the Unix epoch use: sec_since_epoch = TS / 50e6
  *       and of course nanosec_since_epoch = (TS / 50e6) * 1e9 = TS * 20
  *
  * Note: The `CTBFragment` constructor grabs the chunk of memory the artDAQ fragment occupies.
  *       Since we know a priori the number of PTB words (from the TCP header) and the size (128b for all words)
  *       we can loop over the memory addresses, casting each 128b into the correct word type given by each word's `type`.
  *
  * There are 5 word types with the following bit arrangement:
  *       Feeback        = | 3b Word Type | 61b Payload | 64b Timestamp | (here Payload is split into code, source, payload1, payload2)
  *       LLT, HLT       = | 3b Word Type | 61b Payload | 64b Timestamp |
  *       Channel Status = | 3b Word Type | 64b Payload | 61b Timestamp | (Larger Payload to fit all input channels)
  *
  * Word Type:
  *       - 0x0 = Feedback/Error Word -> Errors from the firmware, should abort the run
  *       - 0x1 = Low Level Trigger (LLT) -> Holds a record of any asserted LLTs
  *       - 0x2 = High Level Trigger (HLT) -> Holds a record of any asserted HLTs
  *       - 0x3 = Channel status -> Holds a record of any asserted inputs
  *       - 0x7 = Timestamp Word -> No payload just a timestamp, these are periodic
  *
  * Note: Payload AND'd with mask the size of the expected number of bits just
  *       to make sure we don't get any unexpected garbage bits since we aren't using
  *       the full bits of the variable type e.g. uint64_t, unint16_t..
  **************************************/

void SBNDPTBDecoder::_process_PTB_AUX(const artdaq::Fragment& frag, ptbsv_t &sout)
{
  int numLLT =0;
  sbndaq::CTBFragment ctbfrag(frag);   // somehow the name CTBFragment stuck

  // use the same logic in sbndaq-artdaq-core/Overlays/SBND/PTBFragment.cc: operator<<
  // but separate out the HLTs and LLTs
  if (fDebugLevel > 0)
    {
      std::cout << "SBNDPTBDecoder_module: got into aux" << std::endl;
    }
  
  for (size_t iword = 0; iword < ctbfrag.NWords(); ++iword)
    {
      if (fDebugLevel > 0)
        {
          std::cout << "SBNDPTBDecoder_module: start processing word: " << iword << std::endl;
	}
      size_t ix=0;
      uint32_t wt = 0;
      std::cout << "PTB Word type [" << ctbfrag.Word(iword)->word_type << "]  ";

      if (ctbfrag.Trigger(iword))
	{
	  raw::ptb::Trigger tstruct;
	  // tstruct.word_type = ctbfrag.Trigger(iword)->word_type;
	  tstruct.word_type = ctbfrag.Word(iword)->word_type; 
	  wt = tstruct.word_type;
	  //tstruct.trigger_word = ctbfrag.Trigger(iword)->trigger_word;
	  tstruct.trigger_word = ctbfrag.Trigger(iword)->trigger_word & 0x1FFFFFFFFFFFFFFF; //& 0x1FFFFFFFFFFF extracts the 61 bit payload
	  tstruct.timestamp = ctbfrag.Trigger(iword)->timestamp * 20; //Check note above about timestamps
	  //if (ctbfrag.Trigger(iword)->IsHLT()) 
	  if (wt ==2)  
	    {
	      ix = sout.HLTrigs.size();
	      sout.HLTrigs.push_back(tstruct);
	      std::cout << tstruct.trigger_word << "  HLT Status: " << tstruct.timestamp << std::endl;
	      if (fDebugLevel > 0)
		{
		  std::cout << "SBNDPTBDecoder_module: found HLT: " << wt << " " << ix << std::endl;
		}
	    }
	  //else if (ctbfrag.Trigger(iword)->IsLLT())
	  else if (wt ==1)
	    {
	      numLLT++;
	      ix = sout.LLTrigs.size();
	      sout.LLTrigs.push_back(tstruct);
	      std::cout << tstruct.trigger_word << "  LLT Status: " << tstruct.timestamp <<  std::endl;
	      if (fDebugLevel > 0)
		{
		  std::cout << "SBNDPTBDecoder_module: found LLT: " << wt << " " << ix << std::endl;
		}
	    }
	}
      else if (ctbfrag.ChStatus(iword))
	{
	  raw::ptb::ChStatus cstruct;
	  cstruct.timestamp = ctbfrag.ChStatus(iword)->timestamp * 20; //Check above notes about timestamps
	  cstruct.beam = ctbfrag.ChStatus(iword)->beam;
	  cstruct.crt = ctbfrag.ChStatus(iword)->crt;
	  cstruct.pds = ctbfrag.ChStatus(iword)->pds;
	  cstruct.mtca = ctbfrag.ChStatus(iword)->mtca;
	  cstruct.nim = ctbfrag.ChStatus(iword)->nim;
	  cstruct.auxpds = ctbfrag.ChStatus(iword)->auxpds;
	  cstruct.word_type = ctbfrag.ChStatus(iword)->word_type;
	  wt = cstruct.word_type;
	  std::cout <<cstruct.crt << "                  CRT Status: " << cstruct.timestamp <<  std::endl;
	  ix = sout.ChStats.size();
	  sout.ChStats.push_back(cstruct);
	  if (fDebugLevel > 0)
	    {
	      std::cout << "SBNDPTBDecoder_module: found CHStat: " << wt << " " << ix << std::endl;
	    }
	}
      else if (ctbfrag.Feedback(iword))
	{
	  raw::ptb::Feedback fstruct;
	  fstruct.timestamp = ctbfrag.Feedback(iword)->timestamp * 20;
	  fstruct.code = ctbfrag.Feedback(iword)->code;
	  fstruct.source = ctbfrag.Feedback(iword)->source;
	  fstruct.payload = ctbfrag.Feedback(iword)->payload;  // broken in two in Tereza's version
	  fstruct.word_type = ctbfrag.Feedback(iword)->word_type;
	  wt = fstruct.word_type;
	  ix = sout.Feedbacks.size();
	  sout.Feedbacks.push_back(fstruct);
	  if (fDebugLevel > 0)
	    {
	      std::cout << "SBNDPTBDecoder_module: found Feedback: " << wt << " " << ix << std::endl;
	    }
	}
      else
	{
	  raw::ptb::Misc mstruct;
	  mstruct.timestamp = ctbfrag.Word(iword)->timestamp * 20;
	  mstruct.payload = ctbfrag.Word(iword)->payload;
	  mstruct.word_type = ctbfrag.Word(iword)->word_type;
	  wt = mstruct.word_type;
	  ix = sout.Miscs.size();
	  sout.Miscs.push_back(mstruct);
	  if (fDebugLevel > 0)
	    {
	      std::cout << "SBNDPTBDecoder_module: found Misc: " << wt << " " << ix << std::endl;
	    }
	}

      raw::ptb::WordIndex wstruct;
      wstruct.word_type = wt;
      wstruct.index = ix;
      sout.WordIndexes.push_back(wstruct);
      if (fDebugLevel > 0)
	{
	  std::cout << "SBNDPTBDecoder_module: index calc: " << wt << " " << ix << std::endl;
	}
    }
}


DEFINE_ART_MODULE(SBNDPTBDecoder)
