////////////////////////////////////////////////////////////////////////
// Class:       SBNDTPCDecoder
// Plugin Type: producer (art v2_09_06)
// File:        SBNDTPCDecoder.cxx
//
// Generated at Thu Feb  8 16:41:18 2018 by Gray Putnam using cetskelgen
// from cetlib version v3_01_03.
// moved and adapted for sbndcode by Tom Junk, August 2023
////////////////////////////////////////////////////////////////////////

#include "SBNDTPCDecoder.h"
#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <thread>

#include "TMath.h"

#include "art/Framework/Core/ModuleMacros.h"

#include "artdaq-core/Data/Fragment.hh"

#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPCFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPC/NevisTPCTypes.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPC/NevisTPCUtilities.hh"


DEFINE_ART_MODULE(daq::SBNDTPCDecoder)

// constructs a header data object from a nevis header
// construct from a nevis header
tpcAnalysis::TPCDecodeAna daq::SBNDTPCDecoder::Fragment2TPCDecodeAna(art::Event &event, const artdaq::Fragment &frag) {
  sbndaq::NevisTPCFragment fragment(frag);

  const sbndaq::NevisTPCHeader *raw_header = fragment.header();
  tpcAnalysis::TPCDecodeAna ret;

  ret.crate = (frag.fragmentID() >> 8) & 0xF;  
  ret.slot = raw_header->getSlot();
  ret.event_number = raw_header->getEventNum();
  //std::cout << "TPC decoder frame, sample: " << raw_header->getFrameNum() << " " << raw_header->get2mhzSample() << std::endl;
  ret.checksum = raw_header->getChecksum();
  
  ret.adc_word_count = raw_header->getADCWordCount();
  // ret.trig_frame_number = raw_header->getTrigFrame();
  
  // formula for getting unix timestamp from nevis frame number:
  // timestamp = frame_number * (timesize + 1) + trigger_sample
  ret.timestamp = (raw_header->getFrameNum() * ( (ULong64_t) _config.timesize + 1) + raw_header->get2mhzSample()) * ( (double)  _config.frame_to_dt);
  //std::cout << "timestamp: " << ret.timestamp << std::endl;
  ret.index = raw_header->getSlot() - _config.min_slot_no;

  ret.framenum = raw_header->getFrameNum();
  ret.samplenum = raw_header->get2mhzSample();
  ret.fragtimestamp = (ULong64_t) frag.timestamp();

  return ret;
}

daq::SBNDTPCDecoder::SBNDTPCDecoder(fhicl::ParameterSet const & param): 
  art::EDProducer{param},
  _tag(param.get<std::string>("raw_data_label", "daq"),param.get<std::string>("fragment_type_label", "NEVISTPC")),
  _config(param)
{
  consumes<artdaq::Fragments>(_tag);
  produces<RawDigits>();
  produces<RDTimeStamps>();
  produces<RDTsAssocs>();
  if (_config.produce_header) {
    produces<std::vector<tpcAnalysis::TPCDecodeAna>>();
  }
}

daq::SBNDTPCDecoder::Config::Config(fhicl::ParameterSet const & param) {

  // whether to calcualte the pedestal (and set it in SetPedestal())
  baseline_calc = param.get<bool>("baseline_calc", true);
  // whether to put headerinfo in the art root file
  produce_header = param.get<bool>("produce_header", false);

  // nevis readout window length
  timesize = param.get<unsigned>("timesize", 1);

  // nevis tick length (for timestamp)
  // should be 1/(2MHz) = 0.5mus
  frame_to_dt = param.get<double>("frame_to_dt", 1);

  // number of channels in each slot
  channel_per_slot = param.get<unsigned>("channel_per_slot", 0);
  // index of 0th slot
  min_slot_no = param.get<unsigned>("min_slot_no", 0);
}

void daq::SBNDTPCDecoder::produce(art::Event & event)
{
  auto daq_handle = event.getHandle<artdaq::Fragments>(_tag);
  
  RDPmkr rdpm(event);
  TSPmkr tspm(event);

  // output collections
  std::unique_ptr<RawDigits> rawdigit_collection(new RawDigits);
  std::unique_ptr<RDTimeStamps> rdts_collection(new RDTimeStamps);
  std::unique_ptr<RDTsAssocs> rdtsassoc_collection(new RDTsAssocs);
  std::unique_ptr<std::vector<tpcAnalysis::TPCDecodeAna>> header_collection(new std::vector<tpcAnalysis::TPCDecodeAna>);

  if ( daq_handle.isValid() ) {
    for (auto const &rawfrag: *daq_handle) {
      process_fragment(event, rawfrag, rawdigit_collection, header_collection, rdpm, tspm, rdts_collection, rdtsassoc_collection);
    }
  }
  else
    {
      mf::LogWarning("SBNDTPCDecoder_module") <<  " Invalid fragment handle: Skipping TPC digit decoding";
    }

  
  event.put(std::move(rawdigit_collection));
  event.put(std::move(rdts_collection));
  event.put(std::move(rdtsassoc_collection));

  if (_config.produce_header) {
    event.put(std::move(header_collection));
  }

  // remove TPC fragments from the art event cache.  They can be re-read from the file
  // by downstream processes if need be, but we are
  // done with them here

  if ( daq_handle.isValid() ) daq_handle.removeProduct();
}


void daq::SBNDTPCDecoder::process_fragment(art::Event &event, const artdaq::Fragment &frag, 
					   std::unique_ptr<RawDigits> &rd_collection,
					   std::unique_ptr<std::vector<tpcAnalysis::TPCDecodeAna>> &header_collection,
					   RDPmkr &rdpm,
					   TSPmkr &tspm,
					   std::unique_ptr<RDTimeStamps> &rdts_collection,
					   std::unique_ptr<RDTsAssocs> &rdtsassoc_collection) {

  art::ServiceHandle<SBND::TPCChannelMapService> channelMap;
  
  // convert fragment to Nevis fragment
  sbndaq::NevisTPCFragment fragment(frag);

  std::unordered_map<uint16_t,sbndaq::NevisTPC_Data_t> waveform_map;
  size_t n_waveforms = fragment.decode_data(waveform_map);
  (void)n_waveforms;

  // need to retrieve the timestamp from the Nevis header and save it in the art event only on request
  
  auto header_data = Fragment2TPCDecodeAna(event, frag);
  if (_config.produce_header) {
    header_collection->push_back(header_data);
  }

  unsigned int FEMCrate = (frag.fragmentID() >> 8) & 0xF;
  unsigned int FEMSlot = fragment.header()->getSlot()-_config.min_slot_no + 1;
  
  for (auto waveform: waveform_map) {
    auto chanInfo = channelMap->GetChanInfoFromFEMElements(FEMCrate,
							   FEMSlot,
							   waveform.first); // nevis_channel_id    
    if (!chanInfo.valid) continue;

    std::vector<int16_t> raw_digits_waveform;
    raw::ChannelID_t wire_id = chanInfo.offlchan;

    for (auto digit: waveform.second) {
      raw_digits_waveform.push_back( (int16_t) digit);
    }

    float median = 0;
    float sigma = 0; 
    if (_config.baseline_calc) {
      getMedianSigma(raw_digits_waveform, median, sigma);
    }

    // construct the next RawDigit object
    rd_collection->emplace_back(wire_id, raw_digits_waveform.size(), raw_digits_waveform);
    (*rd_collection)[rd_collection->size() - 1].SetPedestal( median, sigma );

    // construct the RDTimeStamp object and make the association
    rdts_collection->emplace_back(header_data.timestamp,0);
    auto const rawdigitptr = rdpm(rd_collection->size()-1);
    auto const rdtimestampptr = tspm(rdts_collection->size()-1);
    rdtsassoc_collection->addSingle(rawdigitptr,rdtimestampptr);       
  }
}

// Computes the checksum, given a nevis tpc header
//
// Also note that this only works for uncompressed data
uint32_t daq::SBNDTPCDecoder::compute_checksum(sbndaq::NevisTPCFragment &fragment) {
  uint32_t checksum = 0;

  const sbndaq::NevisTPC_ADC_t* data_ptr = fragment.data();
  // RETURN VALUE OF getADCWordCount IS OFF BY 1
  size_t n_words = fragment.header()->getADCWordCount() + 1;

  for (size_t word_ind = 0; word_ind < n_words; word_ind++) {
    const sbndaq::NevisTPC_ADC_t* word_ptr = data_ptr + word_ind;
    checksum += *word_ptr;
  }
  // only first 6 bytes of checksum are used
  return checksum & 0xFFFFFF;

}



void daq::SBNDTPCDecoder::getMedianSigma(const std::vector<int16_t> &v_adc, float &median,
					       float &sigma) {
  size_t asiz = v_adc.size();
  int imed=0;
  if (asiz == 0) {
    median = 0;
    sigma = 0;
  }
  else {
    // the RMS includes tails from bad samples and signals and may not be the best RMS calc.

    imed = TMath::Median(asiz,v_adc.data()) + 0.01;  // add an offset to make sure the floor gets the right integer
    median = imed;
    sigma = TMath::RMS(asiz,v_adc.data());
    
    // add in a correction suggested by David Adams, May 6, 2019
    
    size_t s1 = 0;
    size_t sm = 0;
    for (size_t i = 0; i < asiz; ++i) {
      if (v_adc.at(i) < imed) s1++;
      if (v_adc.at(i) == imed) sm++;
    }
    if (sm > 0) {
      float mcorr = (-0.5 + (0.5*(float) asiz - (float) s1)/ ((float) sm) );
      //if (fDebugLevel > 0)
      //  {
      //    if (std::abs(mcorr)>1.0) std::cout << "mcorr: " << mcorr << std::endl;
      //  }
      median += mcorr;
    }
  }
}
