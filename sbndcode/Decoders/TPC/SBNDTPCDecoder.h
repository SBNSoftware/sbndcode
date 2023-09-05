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
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"

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

  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef art::Assns<raw::RawDigit,raw::RDTimeStamp> RDTsAssocs;
  typedef art::PtrMaker<raw::RawDigit> RDPmkr;
  typedef art::PtrMaker<raw::RDTimeStamp> TSPmkr;
    
  // process an individual fragment inside an art event
  void process_fragment(art::Event &event,
			const artdaq::Fragment &frag,
                        std::unique_ptr<RawDigits> &rd_collection,
                        std::unique_ptr<std::vector<tpcAnalysis::TPCDecodeAna>> &header_collection,
			RDPmkr &rdpm,
			TSPmkr &tspm,
			std::unique_ptr<RDTimeStamps> &rdts_collection,
			std::unique_ptr<RDTsAssocs> &rdtsassoc_collection);


  // build a TPCDecodeAna object from the Nevis Header
  tpcAnalysis::TPCDecodeAna Fragment2TPCDecodeAna(art::Event &event, const artdaq::Fragment &frag);

  art::InputTag _tag;
  Config _config;

  void getMedianSigma(const std::vector<int16_t> &v_adc, float &median, float &sigma);

};

#endif /* SBNDTPCDecoder_h */
