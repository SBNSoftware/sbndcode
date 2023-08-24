#ifndef SBNDTPCDecoder_h
#define SBNDTPCDecoder_h
////////////////////////////////////////////////////////////////////////
// Class:       SBNDTPCDecoder
// Plugin Type: producer (art v2_09_06)
// File:        SBNDTPCDecoder.h
//
// Generated at Thu Feb  8 16:41:18 2018 by Gray Putnam using cetskelgen
// from cetlib version v3_01_03.
// adapted from sbndqm for sbndcode use
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "lardataobj/RawData/RawDigit.h"
#include "artdaq-core/Data/Fragment.hh"
#include "canvas/Utilities/InputTag.h"

#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPCFragment.hh"

#include "TPCDecodeAna.h"

/*
  * The Decoder module takes as input "NevisTPCFragments" and
  * outputs raw::RawDigits. It also handles in and all issues
  * with the passed in header and fragments (or at least it will).
*/

namespace daq {
  class SBNDTPCDecoder;
}


class daq::SBNDTPCDecoder : public art::EDProducer {
public:
  explicit SBNDTPCDecoder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDTPCDecoder(SBNDTPCDecoder const &) = delete;
  SBNDTPCDecoder(SBNDTPCDecoder &&) = delete;
  SBNDTPCDecoder & operator = (SBNDTPCDecoder const &) = delete;
  SBNDTPCDecoder & operator = (SBNDTPCDecoder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // get checksum from a Nevis fragment
  static uint32_t compute_checksum(sbndaq::NevisTPCFragment &fragment);

private:
  class Config {
    public:
    bool produce_header;
    bool produce_metadata;
    bool baseline_calc;
    unsigned n_mode_skip;
    bool subtract_pedestal;

    unsigned channel_per_slot;
    unsigned min_slot_no;

    // for converting nevis frame time into timestamp
    unsigned timesize;
    double frame_to_dt;

    Config(fhicl::ParameterSet const & p);
  };

  // process an individual fragment inside an art event
  void process_fragment(art::Event &event, const artdaq::Fragment &frag,
    std::unique_ptr<std::vector<raw::RawDigit>> &product_collection,
    std::unique_ptr<std::vector<tpcAnalysis::HeaderData>> &header_collection);


  // Gets the WIRE ID of the channel. This wire id can be then passed
  // to the Lariat geometry.
  raw::ChannelID_t get_wire_id(const sbndaq::NevisTPCHeader *header, uint16_t nevis_channel_id);

  // whether the given nevis readout channel is mapped to a wire
  bool is_mapped_channel(const sbndaq::NevisTPCHeader *header, uint16_t nevis_channel_id);

  // build a HeaderData object from the Nevis Header
  tpcAnalysis::HeaderData Fragment2HeaderData(art::Event &event, const artdaq::Fragment &frag);

  art::InputTag _tag;
  Config _config;
  // keeping track of incrementing numbers
  uint32_t _last_event_number;
  uint32_t _last_trig_frame_number;

  void getMedianSigma(const std::vector<int16_t> &v_adc, float &median, float &sigma);
};

#endif /* SBNDTPCDecoder_h */
