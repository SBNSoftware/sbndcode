////////////////////////////////////////////////////////////////////////
// Class:       NTBDecoder
// Plugin Type: producer (Unknown Unknown)
// File:        NTBDecoder_module.cc
//
// Generated at Thu Aug  1 11:38:21 2024 by Thomas Junk using cetskelgen
// from cetlib version 3.18.02.
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

#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndntb.h"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPCFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTB_dataFormat.hh"

class NTBDecoder;


class NTBDecoder : public art::EDProducer {
public:
  explicit NTBDecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NTBDecoder(NTBDecoder const&) = delete;
  NTBDecoder(NTBDecoder&&) = delete;
  NTBDecoder& operator=(NTBDecoder const&) = delete;
  NTBDecoder& operator=(NTBDecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string fTag;
  std::string fContainerTag;
  std::string fOutputInstance;
  int fDebugLevel;

  void process_NTB_AUX(const artdaq::Fragment& frag, std::vector<raw::ntb::sbndntb> &sbndntbs);
};


NTBDecoder::NTBDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  fTag            = p.get<std::string>("ntb_data_label","daq:NEVISTB");
  fContainerTag   = p.get<std::string>("ntb_container_data_label","daq:ContainerNEVISTB");
  fOutputInstance = p.get<std::string>("OutputInstance","");
  fDebugLevel     = p.get<int>("DebugLevel",0);

  consumes<artdaq::Fragments>(fTag);
  consumes<artdaq::Fragments>(fContainerTag);
  produces<std::vector<raw::ntb::sbndntb>>(fOutputInstance);
}

void NTBDecoder::produce(art::Event& e)
{
  // Implementation of required member function here.

    // look first for container fragments and then non-container fragments

  std::vector<raw::ntb::sbndntb> sbndntbs;

  art::InputTag itag1(fContainerTag);
  auto cont_frags = e.getHandle<artdaq::Fragments>(itag1);
  if (cont_frags)
    {
      for (auto const& cont : *cont_frags)
	{
	  artdaq::ContainerFragment cont_frag(cont);
	  for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	    {
	      process_NTB_AUX(*cont_frag[ii], sbndntbs);
	    }
	}
    }

  art::InputTag itag2(fTag);
  auto frags = e.getHandle<artdaq::Fragments>(itag2);
  if (frags)
    {
      for(auto const& frag: *frags)
	{
	  process_NTB_AUX(frag, sbndntbs);
	}
    }

  e.put(std::make_unique<std::vector<raw::ntb::sbndntb>>(std::move(sbndntbs)),fOutputInstance);
}

void NTBDecoder::process_NTB_AUX(const artdaq::Fragment& frag, std::vector<raw::ntb::sbndntb> &sbndntbs)
{
  if (fDebugLevel > 0)
    {
      std::cout << "NTBDecoder_module: got into aux" << std::endl;
    }
  sbndaq::NevisTBFragment tbfrag(frag);
  const auto header = *tbfrag.header();
  const auto fdata = *tbfrag.data();
  const auto trailer = *tbfrag.trailer();
  
  raw::ntb::sbndntb ostruct;
  ostruct.boardreader_timestamp = frag.timestamp();
  ostruct.event_number = tbfrag.metadata()->EventNumber();
  ostruct.frame_number = tbfrag.metadata()->FrameNumber();
  ostruct.sample_number = tbfrag.metadata()->SampleNumber();
  ostruct.busy = header.busy;
  ostruct.sixteen_mhz_remainder = header.remainder;
  ostruct.two_mhz_sample = header.sample;
  ostruct.frame = header.getFrame();
  ostruct.ntrig = header.getTriggerNumber();
  ostruct.pmt_trig_data = fdata.pmt_trig_data;
  ostruct.pc = fdata.pc;
  ostruct.external = fdata.external;
  ostruct.active = fdata.active;
  ostruct.gate2 = fdata.gate2;
  ostruct.gate1 = fdata.gate1;
  ostruct.veto = fdata.veto;
  ostruct.calib = fdata.calib;
  ostruct.phase = fdata.getPhase();
  ostruct.gatefake = fdata.gatefake;
  ostruct.beamfake = fdata.beamfake;
  ostruct.spare1 = fdata.spare1;
  ostruct.trailer = trailer.getTrailerWord();
  sbndntbs.push_back(ostruct);
}

DEFINE_ART_MODULE(NTBDecoder)
