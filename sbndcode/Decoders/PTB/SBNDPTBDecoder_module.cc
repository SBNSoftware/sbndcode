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

void SBNDPTBDecoder::produce(art::Event & evt)
{


  // look first for container fragments and then non-container fragments

  std::vector<raw::ptb::sbndptb> sbndptbs;

  art::InputTag itag1(fInputLabel, fInputContainerInstance);
  auto cont_frags = evt.getHandle<artdaq::Fragments>(itag1);
  if (cont_frags)
    {
      for (auto const& cont : *cont_frags)
	{
	  artdaq::ContainerFragment cont_frag(cont);
	  for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	    {
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

void SBNDPTBDecoder::_process_PTB_AUX(const artdaq::Fragment& frag, ptbsv_t &sout)
{
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
      if (ctbfrag.Trigger(iword))
	{
	  raw::ptb::Trigger tstruct;
	  tstruct.word_type = ctbfrag.Trigger(iword)->word_type;
	  wt = tstruct.word_type;
	  tstruct.trigger_word = ctbfrag.Trigger(iword)->trigger_word;
	  tstruct.timestamp = ctbfrag.Trigger(iword)->timestamp;
	  if (ctbfrag.Trigger(iword)->IsHLT()) 
	    {
	      ix = sout.HLTrigs.size();
	      sout.HLTrigs.push_back(tstruct);
	      if (fDebugLevel > 0)
		{
		  std::cout << "SBNDPTBDecoder_module: found HLT: " << wt << " " << ix << std::endl;
		}
	    }
	  else if (ctbfrag.Trigger(iword)->IsLLT())
	    {
	      ix = sout.LLTrigs.size();
	      sout.LLTrigs.push_back(tstruct);
	      if (fDebugLevel > 0)
		{
		  std::cout << "SBNDPTBDecoder_module: found LLT: " << wt << " " << ix << std::endl;
		}
	    }
	}
      else if (ctbfrag.ChStatus(iword))
	{
	  raw::ptb::ChStatus cstruct;
	  cstruct.timestamp = ctbfrag.ChStatus(iword)->timestamp;
	  cstruct.beam = ctbfrag.ChStatus(iword)->beam;
	  cstruct.crt = ctbfrag.ChStatus(iword)->crt;
	  cstruct.pds = ctbfrag.ChStatus(iword)->pds;
	  cstruct.mtca = ctbfrag.ChStatus(iword)->mtca;
	  cstruct.nim = ctbfrag.ChStatus(iword)->nim;
	  cstruct.auxpds = ctbfrag.ChStatus(iword)->auxpds;
	  cstruct.word_type = ctbfrag.ChStatus(iword)->word_type;
	  wt = cstruct.word_type;
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
	  fstruct.timestamp = ctbfrag.Feedback(iword)->timestamp;
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
	  mstruct.timestamp = ctbfrag.Word(iword)->timestamp;
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
