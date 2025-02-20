////////////////////////////////////////////////////////////////////////
// Class:       SBNDXARAPUCADecoder
// Plugin Type: producer
// File:        SBNDXARAPUCADecoder_module.cc
// Author:      Alicia Vázquez-Ramos (aliciavr@ugr.es)
////////////////////////////////////////////////////////////////////////

/**
 * @file SBNDXARAPUCADecoder_module.cc
 * @brief Defines and implements the SBNDXARAPUCADecoder class which inherits from an art::EDProducer
 * as the decoder for V1740B digitizers, intended for the X-ARAPUCAs.
 * @details The current version of the SBND X-ARAPUCAs decoder implements the updates shown in 
 * the SBN Document 38475-v1 in the SBN Document Database.
 * @note A Python version of the binary decoding is available for testing purposes. You can find 
 * it [here: V1740 binary decoder](https://github.com/aliciavr/V1740_binary_decoder).
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/CAENV1740Fragment.hh"

#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "art_root_io/TFileService.h"
#include "TH1I.h"

#include <memory>
#include <algorithm>


/**
 * @namespace sbndaq
 * @brief Contains classes and functions for inside the sbndaq system.
 */
namespace sbndaq {
  class SBNDXARAPUCADecoder;
}


/**
 * @class SBNDXARAPUCADecoder
 * @brief A class providing the V1740B digitizer decoder. It inherits from art::EDProducer and 
 * creates a vector of raw::OpDetWaveform objects as a product of the decoding.
 */
class sbndaq::SBNDXARAPUCADecoder : public art::EDProducer {
public:
  explicit SBNDXARAPUCADecoder(fhicl::ParameterSet const& p);

  SBNDXARAPUCADecoder(SBNDXARAPUCADecoder const&) = delete;
  SBNDXARAPUCADecoder(SBNDXARAPUCADecoder&&) = delete;
  SBNDXARAPUCADecoder& operator=(SBNDXARAPUCADecoder const&) = delete;
  SBNDXARAPUCADecoder& operator=(SBNDXARAPUCADecoder&&) = delete;

  // Main class method.
  void produce(art::Event& e) override;

private:

  constexpr static uint64_t NANOSEC_IN_SEC = 1'000'000'000; /**< Number of nanoseconds in one second. */
  constexpr static uint64_t MICROSEC_IN_NANOSEC = 1'000; /**< Number of nanoseconds in one microsecond. */
  constexpr static double NANOSEC_TO_MICROSEC = 1E-3; /**< Conversion factor from nanoseconds to microseconds. */

  constexpr static uint32_t BITS_PER_WORD = 32; /**< Number of bits per word in the V1740B digitizer (32-bit word). */
  constexpr static uint32_t BITS_PER_SAMPLE = 12; /**< Bit precision of the V1740B digitizer (12-bit resolution per sample). */
  constexpr static uint32_t NUM_CHANNELS_PER_GROUP = 8; /**< Number of channels in each channel group of the V1740B digitizers. */
  constexpr static uint32_t NUM_GROUPS = 8; /**< Number of channel groups for the V1740B digitizers. */
  constexpr static uint32_t NUM_SAMPLES_PER_ROUND = 24; /**< Number of samples to be read in each channel group round. */
  constexpr static uint32_t NUM_CONSECUTIVE_SAMPLES = 3; /**< Number of consecutive samples for the same channel when decoding waveforms. */
  constexpr static uint16_t NANOSEC_PER_TICK = 8; /**< Number of nanoseconds per Trigger Time Tag (TTT) tick. */

  constexpr static uint16_t SPEC_TDC_TIMING = 0; /**< Timing reference frame: SPEC-TDC Event Trigger Timestamp (ETT). */
  constexpr static uint16_t CAEN_ONLY_TIMING = 1; /**< Timing reference frame: CAEN-only. */

  unsigned int fevent_counter; /**< Event counter. */

  unsigned int fnum_caen_boards; /**< Maximum number of CAEN boards to be considered. */
  unsigned int ffragment_id_offset; /**< Offset value to get the slot of every CAEN board. */
  std::vector<unsigned int> fboard_id_list; /**< Possible CAEN board identifiers based on the slot they are installed. */

  std::string fcaen_module_label; /**< Label identifying the module where the CAEN fragments come from. */
  std::vector<std::string> fcaen_fragment_names; /**< Fragments accepted by this decoder. */

  std::string fwaveforms_instance_name; /**< Name assigned to the product instance containing the waveforms for each channel and board generated by this art::EDProducer. */
  std::string ftiming_instance_name; /**< Name assigned to the product instance containing the timing reference information for each event generated by this art::EDProducer. */

  uint16_t ftiming_priority; /**< Timing priority configured: 0 SPEC-TDC ETT, 1 CAEN-only. */
  uint16_t factive_timing_frame; /**< Active timing frame while processing each event. */
  unsigned int fns_per_sample; /**< Number of nanoseconds per sample. */

  std::string fspectdc_product_name; /**< Name assigned to the product instance containing the SPEC-TDC decoder products. */
  uint32_t fspectdc_ftrig_ch; /** Channel assigned by the SPEC-TDC to the flash triggers. */
  uint32_t fspectdc_etrig_ch; /** Channel assigned by the SPEC-TDC to the event triggers. */

  art::ServiceHandle<art::TFileService> tfs; /**<  ServiceHandle object to store the histograms in the decoder_hist.root output file. */
  int fstore_debug_waveforms; /**< Number of waveforms to store in the ServiceHandle object for debugging purposes (0: none, -1: all, n: first n waveforms each event). */

  bool fdebug_all; /**< If `true` all debug information is printed. */
  bool fdebug_tdc_handle; /**< If `true` SPEC-TDC art::Handle information is printed. */
  bool fdebug_fragments_handle; /**< If `true` V1740B CAEN fragments art::Handle information is printed. */
  bool fdebug_timing; /**< If `true` timing data is printed. */
  bool fdebug_buffer; /**< If `true` the buffer status is printed. */
  bool fdebug_waveforms; /**< If `true` waveforms decoding data is printed. */
  bool fverbose; /**< If `true` it increases verbosity of console output for detailed processing steps. */

  // Class methods.
  void decode_fragment(uint64_t event_trigger_timestamp, std::vector<size_t> & fragment_indices, const artdaq::Fragment& fragment, std::vector <raw::OpDetWaveform>& prod_wvfms);
  void save_prod_wvfm(size_t board_idx, size_t ch, double ini_wvfm_timestamp, const std::vector <std::vector <uint16_t> > & wvfms, std::vector <raw::OpDetWaveform> & prod_wvfms);
  void save_debug_wvfm(size_t board_idx, size_t fragment_idx, int ch, double ini_wvfm_timestamp, double end_wvfm_timestamp, const std::vector <std::vector <uint16_t> > & wvfms);
  uint16_t get_sample(uint64_t buffer, uint32_t msb, uint32_t lsb);
  uint32_t read_word(const uint32_t* & data_ptr);
  unsigned int get_channel_id(unsigned int board, unsigned int board_channel);
  
};

/**
 * @brief Constructor of the class SBNDXARAPUCADecoder. Given a fhicl::ParameterSet object it gets the configuration needed by this 
 * art::EDProducer from the FHiCL file.
 * @param[in] p A set of parameters containing configuration values for the decoder from the FHiCL file.
 * @details This constructor accesses the parameters in the FHiCL configuration file and set the configuration to be run for 
 * this module. It sets settings for CAEN fragments, board information, timing, debug options, verbosity option and creates an 
 * instance of the product associated with this module.
 */
sbndaq::SBNDXARAPUCADecoder::SBNDXARAPUCADecoder(fhicl::ParameterSet const& p)
: EDProducer{p}
{
  // Sets the event counter to zero.
  fevent_counter = 0;
  
  // Gets the CAEN fragments information.
  fcaen_module_label = p.get<std::string> ("caen_module_label", "daq");
  fcaen_fragment_names = p.get<std::vector <std::string> > ("caen_fragment_names", {"CAENV1740", "ContainerCAENV1740"});

  // Gets the CAEN boards information.
  ffragment_id_offset = p.get<unsigned int> ("fragment_id_offset", 41216);
  fboard_id_list = p.get<std::vector <unsigned int> > ("board_id_list", {7, 13, 16, 19});
  fnum_caen_boards = p.get<unsigned int> ("num_caen_boards", 4);
  
  // Gets the name of the instance created by this module.
  fwaveforms_instance_name = p.get<std::string> ("waveforms_instance_name", "XARAPUCAChannels");
  ftiming_instance_name = p.get<std::string> ("timing_instance_name", "");

  // Gets timing configuration.
  ftiming_priority = p.get<unsigned int> ("timing_priority", SPEC_TDC_TIMING);
  factive_timing_frame = CAEN_ONLY_TIMING;
  fns_per_sample = p.get<unsigned int> ("ns_per_sample", 16);

  // SPEC-TDC access configuration.
  fspectdc_product_name = p.get<std::string>("spectdc_product_name", "tdcdecoder");
  fspectdc_ftrig_ch = p.get<uint32_t>("spectdc_ftrig_ch", 3);
  fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch", 4);
  
  // Gets the number of waveforms to store in the debug output file.
  fstore_debug_waveforms = p.get<int> ("store_debug_waveforms", 0);

  // Gets the debug and verbose options.
  fdebug_all = p.get<bool> ("debug_all", false);
  fdebug_tdc_handle = p.get<bool> ("debug_tdc_handle", false);
  fdebug_fragments_handle = p.get<bool> ("debug_fragments_handle", false);
  fdebug_timing = p.get<bool> ("debug_timing", false);
  fdebug_buffer = p.get<bool> ("debug_buffer", false);
  fdebug_waveforms = p.get<bool> ("debug_waveforms", false);
  fverbose = p.get<bool> ("verbose", false);

  // Creates the instance product of this module.
  produces <std::vector <raw::OpDetWaveform> > (fwaveforms_instance_name);
  produces <raw::TimingReferenceInfo> (ftiming_instance_name);
}

/**
 * @brief Main function of the art::EDProducer: the produce function analyzes every event producing waveforms after the decoding.
 * @param[in] e The event to be processed.
 * @details It is the main function of the art::EDProducer module:
 * 1. Accesses the products from the SPEC TDC decoder if the priority indicates the temporal reference frame is SPEC-TDC. If 
 * the priority does not match, it defaults to the CAEN-only temporal reference frame.
 * 2. Accesses the products from V1740 CAEN fragments and processes each of them, extracting the header, timing data, and raw waveforms.
 * 3. Dumps the products.
 */
void sbndaq::SBNDXARAPUCADecoder::produce(art::Event& e)
{
  if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: entering the produce function." << std::endl;

  // Advances the event counter.
  fevent_counter++;

  // Initializes the output instance products.
  auto prod_wvfms = std::make_unique <std::vector <raw::OpDetWaveform> > ();
  auto prod_event_timing_ref_info = std::make_unique <raw::TimingReferenceInfo> ();

  // CAEN-only timing frame as default.
  prod_event_timing_ref_info->timingType = CAEN_ONLY_TIMING;
  prod_event_timing_ref_info->timingChannel = 0; 

  if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: products initialized." << std::endl;

  uint64_t event_trigger_timestamp = 0; // By default in case SPEC-TDC Event Trigger is not found.

  // Gets the SPEC-TDC product.
  if (ftiming_priority == SPEC_TDC_TIMING) {

    art::Handle <std::vector<sbnd::timing::DAQTimestamp> > tdc_handle;
    e.getByLabel(fspectdc_product_name, tdc_handle);

    // The art::Handle object is not valid.
    if (!tdc_handle.isValid()) {
      if (fdebug_tdc_handle | fdebug_all) std::cout << "\nTDC-SPEC handle not valid for " << fspectdc_product_name << "." << std::endl;
    
    // The art::Handle object is empty.
    } else if (tdc_handle->empty()) {
      if (fdebug_tdc_handle | fdebug_all) std::cout << "\nTDC-SPEC handle is empty." << std::endl;
    
    // The art::Handle object is valid and not empty
    } else  {
      if (fdebug_tdc_handle | fdebug_all) std::cout << " \nDecoded TDC-SPEC products found: " << tdc_handle->size() << " products." << std::endl;
      
      unsigned int num_event_triggers = 0;
      unsigned int num_flash_triggers = 0;
      
      if (fdebug_tdc_handle | fdebug_all) std::cout << "\tTDC Channel \t TDC Name \t TDC Timestamp [ns] \t\t\t\t\t TDC Offset [ns]" << std::endl;
      
      for (size_t t = 0; t < tdc_handle->size(); t++) {
        const uint64_t tdc_timestamp = tdc_handle->at(t).Timestamp();   // Timestamp of the signal [ns].
        const uint32_t tdc_channel = tdc_handle->at(t).Channel();       // Hardware channel.
        
        // Counts the number of flash triggers.
        if (tdc_channel == fspectdc_ftrig_ch) {
          num_flash_triggers++;
        }

        // Counts the number of event triggers and gets the last event trigger timestamp for the SPEC-TDC reference timing frame.
        if (tdc_channel == fspectdc_etrig_ch) {
          num_event_triggers++;
          event_trigger_timestamp = tdc_timestamp;

          // Updates the Time Reference Information product and active reference timing frame.
          prod_event_timing_ref_info->timingType = SPEC_TDC_TIMING;
          prod_event_timing_ref_info->timingChannel = fspectdc_etrig_ch;
          factive_timing_frame = SPEC_TDC_TIMING;

          if (fdebug_timing | fdebug_all) {
            std::cout << "\t\t Found SPEC-TDC timestamp for Event Trigger: " << event_trigger_timestamp << " ns." << std::endl;
            std::cout << "\t\t Changing current timing type to SPEC-TDC timing frame [Timing type: " 
            << prod_event_timing_ref_info->timingType << ", channel: " << prod_event_timing_ref_info->timingChannel << "]." << std::endl;
          }
        }

        if (fdebug_tdc_handle | fdebug_all) {
          const std::string tdc_name = tdc_handle->at(t).Name();        // Name of channel input.
          const uint64_t tdc_offset = tdc_handle->at(t).Offset();       // Channel specific offset [ns].
          std::cout << "\t " << tdc_channel << "   \t\t  " << tdc_name 
                    << "\t\t " << tdc_timestamp << " (" << tdc_timestamp / NANOSEC_IN_SEC << " s " << tdc_timestamp % NANOSEC_IN_SEC << " ns." << ")"
                    << "\t  " << tdc_offset << std::endl;
        }
      } // End SPEC-TDC products extraction loop.

      if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: event trigger timestamp: " << event_trigger_timestamp << " ns." << std::endl;
      if (fdebug_tdc_handle | fdebug_all) {
        std::cout << "\t Number of event triggers: " << num_event_triggers << "." << std::endl;
        std::cout << "\t Number of flash triggers: " << num_flash_triggers << "." << std::endl;
      }
    }

    tdc_handle.removeProduct();
  }

  if (fverbose) {
    std::cout << "\n > SBNDXARAPUCADecoder::produce: active timing frame: ";

    // SPEC-TDC found.
    if (factive_timing_frame == SPEC_TDC_TIMING) {
      std::cout << "SPEC-TDC timing frame.";
    // CAEN_ONLY_TIMING found or by default.
    } else {                                      
      if (ftiming_priority == CAEN_ONLY_TIMING) {
        std::cout << "CAEN-only timing frame.";
      } else if (ftiming_priority == SPEC_TDC_TIMING) {
        std::cout << "SPEC-TDC timing frame not found. Using CAEN-only timing frame as default.";
      } else {
        std::cout << "Unknown timing frame (" << ftiming_priority << "). Using CAEN-only timing frame as default.";
      }
    }

    std::cout << " [Timing type: " << prod_event_timing_ref_info->timingType << ", channel: " << prod_event_timing_ref_info->timingChannel << "]." << std::endl;
  }

  // Flag to track if valid CAEN fragments are found.
  bool found_caen = false;

  std::vector<size_t> fragment_indices(fnum_caen_boards, 0);

  if (fverbose | fdebug_fragments_handle | fdebug_all) std::cout << "\n > SBNDXARAPUCADecoder::produce: searching for V1740 fragments..." << std::endl;
  
  for (const std::string &caen_name: fcaen_fragment_names) {
    art::Handle <std::vector <artdaq::Fragment> > fragment_handle;
    e.getByLabel(fcaen_module_label, caen_name, fragment_handle);
    
    // The art::Handle object is not valid.
    if (!fragment_handle.isValid()) { 
      if (fdebug_fragments_handle | fdebug_all) std::cout << "\tHandle not valid for " << caen_name << "." << std::endl;
    
    // The art::Handle object is empty.
    } else if (fragment_handle->size() == 0) { 
      if (fdebug_fragments_handle | fdebug_all) std::cout << "\tHandle with size " << fragment_handle->size() << " for " << caen_name << "."<< std::endl;
    
    // The art::Handle object is valid and not empty.
    } else { 
      found_caen |= true;
      size_t frag_handle_size = fragment_handle->size();

      if (fdebug_fragments_handle | fdebug_all) std::cout << "\tHandle valid for " << caen_name << " (" << frag_handle_size << " objects)." << std::endl;
    
      // It is a group of containers containing a group of CAEN V1740 fragments.
      if (fragment_handle->front().type() == artdaq::Fragment::ContainerFragmentType){
        
        // For every container,
        for (size_t c = 0; c < frag_handle_size; c++) {
          auto container = fragment_handle->at(c);
          artdaq::ContainerFragment container_fragment(container);
          
          // it searches for all CAEN V1740 fragments inside the container if it contains CAEN V1740 fragments.
          if (container_fragment.fragment_type() == sbndaq::detail::FragmentType::CAENV1740) {
            size_t num_caen_fragments = container_fragment.block_count(); 
            if (fdebug_fragments_handle | fdebug_all) std::cout << "\t\tContainerCAENV1740 "<< c << " - CAENV1740 fragments found: " << num_caen_fragments << "." << std::endl;
            
            for (size_t f = 0; f < num_caen_fragments; f++) {
              const artdaq::Fragment fragment = *container_fragment[f].get();
              decode_fragment(event_trigger_timestamp, fragment_indices, fragment, *prod_wvfms);
            } // End CAEN V1740 fragments loop.
          }
        } // End Container fragments loop.
    
      // It is a group of CAEN V1740 fragments.
      } else if (fragment_handle->front().type() == sbndaq::detail::FragmentType::CAENV1740) {
        if (fdebug_fragments_handle | fdebug_all) std::cout << "\t\tCAENV1740 fragments found: " << frag_handle_size << "." << std::endl;
        
        // It searches for all CAEN V1740 fragments.
        for (size_t f = 0; f < frag_handle_size; f++) {
          const artdaq::Fragment fragment = fragment_handle->at(f);
          decode_fragment(event_trigger_timestamp, fragment_indices, fragment, *prod_wvfms);
        } // End CAEN V1740 fragments loop.
      }
    } // End extracting CAEN V1740 fragments.
    fragment_handle.removeProduct();
  } // End CAEN fragment names loop.

  // Dumps the products of this art::EDProducer.
  e.put(std::move(prod_wvfms), fwaveforms_instance_name);
  e.put(std::move(prod_event_timing_ref_info), ftiming_instance_name);

  if (fverbose) {
    if (!found_caen) {
      std::cout << "\n > SBNDXARAPUCADecoder::produce: no CAEN V1740 fragments of any type found, pushed empty waveforms." << std::endl;
    } else {
      std::cout << "\n > SBNDXARAPUCADecoder::produce: waveforms for all valid fragments dumped." << std::endl;
    }
    std::cout << "\n > SBNDXARAPUCADecoder::produce: end of the producer function.\n\n" << std::endl;
  } 
}

/**
 * @brief This function processes a single fragment from a CAEN V1740 and stores the decoded waveforms in the provided product 
 * container.
 * 
 * @param[in,out] fragment_indices A 1D vector tracking the number of fragments processed for each board.
 * @param[in] fragment The input CAEN V1740 fragment containing raw data to be decoded.
 * @param[out] prod_wvfms Vector where the decoded waveforms will be stored as raw::OpDetWaveform objects.
 * 
 * @details 
 * - Identifies the board index corresponding to the fragment ID.
 * - Verifies the fragment ID against known boards and ensures it is within a valid range.
 * - 
 * - Decodes the fragment reading raw 32-bit words from the fragment, storing them in a buffer, 
 *    and extracting 12-bit samples. The decoding includes:
 *   - Accessing metadata for the number of channels, samples, and time information.
 *   - Binary decoding of the raw waveform. To assign efficiently each sequential sample extracted. These indices formulas are 
 *    applied:
 *      - The board channel index:
 *      \f[
 *        c = \left( \frac{S}{3} \mod 8 \right) + g \times 8
 *      \f]
 *      - The channel sample index: 
 *      \f[
 *        s = (S \mod 3) + \left( \frac{S}{24} \times 3 \right) \mod s_{w}}
 *      \f]
 *    Where the group index is computed as \f$ \frac{S}{s_{g}} \f$.
 *   - Mapping samples to corresponding channels.
 * - Populates the output vector (`prod_wvfms`) with decoded waveforms and optionally generates debug waveforms output.
 * 
 */
void sbndaq::SBNDXARAPUCADecoder::decode_fragment(uint64_t event_trigger_timestamp, std::vector<size_t> & fragment_indices, const artdaq::Fragment& fragment, std::vector <raw::OpDetWaveform>& prod_wvfms) {
  auto fragment_id = fragment.fragmentID() - ffragment_id_offset;
  auto it = std::find(fboard_id_list.begin(), fboard_id_list.end(), fragment_id);
  size_t board_idx;
  bool valid_fragment = false;

  if (it != fboard_id_list.end()) {
    board_idx = it - fboard_id_list.begin();
    if (board_idx >= fnum_caen_boards) {
      if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: fragment ID " << fragment_id << " (" << board_idx << ") is out of range. Skipping this fragment..." << std::endl;
    } else {
      valid_fragment = true;
    }
  } else {
      if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: fragment ID " << fragment_id << " is not valid. Skipping this fragment..." << std::endl;
  }

  if (valid_fragment) {
    if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::decode_fragment: decoding V1740 CAEN fragment " << fragment_indices[board_idx] << " from the board " << board_idx << " (slot " << fboard_id_list[board_idx] << "):" << std::endl;

    // ===============  Accesses Event metadata and Event header for this fragment =============== //

    CAENV1740Fragment caen_fragment(fragment);
    CAENV1740FragmentMetadata const* metadata = caen_fragment.Metadata();
    uint32_t num_channels = metadata->nChannels;
    if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: number of channels: " << num_channels << std::endl;

    // Accesses the event and header data of the CAEN fragment.
    CAENV1740Event const* event = caen_fragment.Event();
    CAENV1740EventHeader header = event->Header;
    
    // Gets the number of words of the header and the waveforms.
    uint32_t num_words_per_event = header.eventSize;
    uint32_t num_words_per_header = sizeof(CAENV1740EventHeader) / sizeof(uint32_t);
    uint32_t num_words_per_wvfms = (num_words_per_event - num_words_per_header);

    uint32_t num_bits_per_all_wvfms = num_words_per_wvfms * BITS_PER_WORD;
    uint32_t num_samples_per_all_wvfms =  num_bits_per_all_wvfms / BITS_PER_SAMPLE;
    uint32_t num_remaining_bits = num_bits_per_all_wvfms % BITS_PER_SAMPLE;
    uint32_t num_samples_per_wvfm = num_samples_per_all_wvfms / num_channels;
    uint32_t num_samples_per_group = num_samples_per_wvfm * NUM_CHANNELS_PER_GROUP;

    if (fverbose | fdebug_waveforms | fdebug_all) {
      if (metadata->nSamples == num_samples_per_wvfm) {
        std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: [NOMINAL FRAGMENT] " << num_samples_per_wvfm << " samples/waveform." << " (" << num_samples_per_group << " samples per group - 8 channels per group -)." << std::endl;
      } else {
        std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: [EXTENDED FRAGMENT] " << num_samples_per_wvfm << " samples/waveform." << " (" << num_samples_per_group << " samples per group - 8 channels per group -)." << std::endl;
      }
      std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: nominal number of samples per waveform: " << metadata->nSamples << "." << std::endl;
      std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: number of words for this fragment: " << num_words_per_event << " (Header: " << num_words_per_header << ", Waveform: " << num_words_per_wvfms << ") words." << std::endl;
    }

    if (fdebug_waveforms | fdebug_all) {
      std::cout << "\t Number of bits for all the waveforms of this fragment: " << BITS_PER_WORD << "\t" << num_bits_per_all_wvfms << std::endl;
      std::cout << "\t Number of samples for all the waveforms of this fragment: " << num_samples_per_all_wvfms << std::endl;
      std::cout << "\t Number of remaining bits for this fragment: " << num_remaining_bits << std::endl;
      std::cout << "\t Number of samples per wvfm (this fragment): " << num_samples_per_wvfm << std::endl;    
      std::cout << "\t Number of samples per group (this fragment): " << num_samples_per_group << std::endl;
    }

    // ===============  Extracts timing information for this fragment =============== //

    // Gets the timing information of the CAEN fragment.
    uint32_t TTT_end_ticks = header.triggerTime();
    int64_t TTT_end_ns = TTT_end_ticks * NANOSEC_PER_TICK;
    int64_t TTT_ini_ns = TTT_end_ns - num_samples_per_wvfm * fns_per_sample;
    double TTT_end_us = TTT_end_ns * NANOSEC_TO_MICROSEC; // us.
    double TTT_ini_us = TTT_ini_ns * NANOSEC_TO_MICROSEC; // us.

    // SPEC-TDC timestamp associated to this CAEN fragment.
    int64_t event_trigger_timestamp_ns = event_trigger_timestamp % NANOSEC_IN_SEC; // ns.
    
    // If a SPEC-TDC timestamp was found it restarts the time from it. Otherwise the CAEN time frame is assigned.
    double ini_wvfm_timestamp = (TTT_ini_ns - event_trigger_timestamp_ns) * NANOSEC_TO_MICROSEC; // us.
    double end_wvfm_timestamp = (TTT_end_ns - event_trigger_timestamp_ns) * NANOSEC_TO_MICROSEC; // us.

    if (fverbose | fdebug_timing | fdebug_all) {
      std::cout << std::fixed << std::setprecision(3);
      if (factive_timing_frame == SPEC_TDC_TIMING) {
        std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: SPEC-TDC time window of " << end_wvfm_timestamp - ini_wvfm_timestamp << " us: [" << ini_wvfm_timestamp << ", " << end_wvfm_timestamp << "] us." << std::endl;
      } else {
        std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: CAEN time window of " << TTT_end_us - TTT_ini_us << " us: [" << TTT_ini_us << ", " << TTT_end_us << "] us." << std::endl;
      }
    }

    if (fdebug_timing | fdebug_all) {
      std::cout << "\t ns/tick = " << NANOSEC_PER_TICK << ", ns/sample = " << fns_per_sample << std::endl;
      std::cout << "\t CAEN trigger timestamp (TTT) of the fragment: " << std::endl;
      std::cout << "\t\tTTT ini " << TTT_ini_ns << " ns = " << TTT_ini_us << " us." << std::endl;
      std::cout << "\t\tTTT end " << TTT_end_ticks << " ticks = " << TTT_end_ns << " ns = " << TTT_end_us << " us." << std::endl;
      if (factive_timing_frame == SPEC_TDC_TIMING) {
        std::cout << "\t SPEC-TDC timestamp of the fragment: " << std::endl;
        std::cout << "\t\tSPEC-TDC difference applied to the CAEN frame: " << TTT_ini_ns << " - " << "("<< event_trigger_timestamp / NANOSEC_IN_SEC << " s) " << event_trigger_timestamp % NANOSEC_IN_SEC << " ns." << std::endl;
      }
    }

    // ===============  Start decoding the waveforms =============== //
    if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: binary decoding of the waveforms starting... " << std::endl;
    
    std::vector <std::vector <uint16_t> > wvfms(num_channels, std::vector<uint16_t>(num_samples_per_wvfm, 0));

    // Absolute sample number [0, TOTAL_NUM_SAMPLES] where TOTAL_NUM_SAMPLES is the total number of samples stored for an event.
    uint32_t S = 0;
    // Buffer variables.
    uint64_t buffer = 0;
    uint32_t bits_in_buffer = 0;

    // Data pointer to the beggining of the waveforms stores in the event.
    const uint32_t* data_ptr = reinterpret_cast<const uint32_t*>(fragment.dataBeginBytes() + sizeof(CAENV1740EventHeader));
    // Accesses each word, stores it in the buffer and then the samples are extracted from the buffer.
    for (size_t j = 0; j < num_words_per_wvfms; j++) {
      uint64_t word = read_word(data_ptr);

      // Adds the new word to the buffer and increments the number of bits stored in it.
      if (fdebug_buffer | fdebug_all) std::cout << buffer << "[word: " << word << "]" << std::endl;
      buffer |= word << bits_in_buffer;
      bits_in_buffer += BITS_PER_WORD; // bytes * 8 bits/byte
      if (fdebug_buffer | fdebug_all) std::cout << "  +" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;

      // Obtains 12-bit sequences from the buffer and assigns each sample to the channel and channel sample it belongs to. 
      while (bits_in_buffer >= BITS_PER_SAMPLE) {
        // Computes board channel, channel sample and group channel and assigns the sample to those indices.
        uint32_t g = (S / num_samples_per_group);                                                                                       // Group index.
        uint32_t c = ((S / NUM_CONSECUTIVE_SAMPLES) % NUM_CHANNELS_PER_GROUP) + g * NUM_GROUPS;                                         // Channel index.
        uint32_t s = (S % NUM_CONSECUTIVE_SAMPLES) + ((S / NUM_SAMPLES_PER_ROUND) * NUM_CONSECUTIVE_SAMPLES) % num_samples_per_wvfm;    // Sample/channel index.
        uint16_t sample = get_sample(buffer, BITS_PER_SAMPLE - 1, 0);
        wvfms[c][s] = sample;
        if (fdebug_waveforms | fdebug_all) std::cout << "\tSample: " << sample << "\tg: " << g <<  "\tch: " << c << "\ts:" << s << "\tS: " << S << std::endl;
        
        // Updates the buffer status removing the read bits and decreasing the number of bits stored in it.
        buffer >>= BITS_PER_SAMPLE;
        bits_in_buffer -= BITS_PER_SAMPLE;
        if (fdebug_buffer | fdebug_all) std::cout << "  -" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;
        
        // Increments the absolute sample step.
        S++;
      }
    }
    
    // The decoded waveforms are dumped into two products:
    // - A xarapucadecoder-art.root file with the OpDetWaveforms as the product of this producer for further analysis.
    // - A decoder_hist.root file gathering a waveform histograms.
    if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::decode_fragment: binary decoding complete, dumping products..." << std::endl;
    
    uint32_t num_debug_wvfms;

    if (fstore_debug_waveforms == -1) {
      num_debug_wvfms = num_channels;
    } else {
      num_debug_wvfms = std::min<size_t>(num_channels, fstore_debug_waveforms);
    }

    uint32_t ch;

    for (ch = 0; ch < num_debug_wvfms; ch++) {
      save_prod_wvfm(board_idx, ch, ini_wvfm_timestamp, wvfms, prod_wvfms);
      save_debug_wvfm(board_idx, fragment_indices[board_idx], ch, ini_wvfm_timestamp, end_wvfm_timestamp, wvfms);
    }

    for (;ch < num_channels; ch++) {
      save_prod_wvfm(board_idx, ch, ini_wvfm_timestamp, wvfms, prod_wvfms);
    }

    fragment_indices[board_idx]++;
  }
}

/**
 * @brief Converts a production waveform from raw ADC data into a `raw::OpDetWaveform` object,
 * assigning it a global channel ID and timestamp, and appends it to the output collection.
 *
 * @param[in] b The board index (position in the list of boards).
 * @param[in] ch The channel index (channel number from which the waveform is extracted).
 * @param[in] ini_wvfm_timestamp The initial timestamp of the waveform in microseconds.
 * @param[in] wvfms A 2D vector containing the waveforms. Each inner vector corresponds to a channel,
 * and each element of the inner vector represents a sample (ADC count).
 * @param[out] prod_wvfms A reference to the vector where the produced `raw::OpDetWaveform` objects
 * will be stored.
 *
 * @details
 * The function performs the following steps:
 * 1. Determines the global channel ID for the specified board and channel using `get_channel_id`.
 * 2. Creates a `raw::OpDetWaveform` object with the initial timestamp, global channel ID, and waveform data.
 * 3. Appends the `raw::OpDetWaveform` object to the provided `prod_wvfms` collection.
 * 
 * @pre ini_wvfm_timestamp is asumed to be given in microseconds.
 *
 * @see raw::OpDetWaveform
 */
void sbndaq::SBNDXARAPUCADecoder::save_prod_wvfm(size_t board_idx, size_t ch, double ini_wvfm_timestamp, const std::vector <std::vector <uint16_t> > & wvfms, std::vector <raw::OpDetWaveform> & prod_wvfms) {
  unsigned int channel_id = get_channel_id(board_idx, ch);
  raw::OpDetWaveform waveform(ini_wvfm_timestamp, channel_id, wvfms[ch]);
  if (fdebug_waveforms | fdebug_all) std::cout << "Pushing waveform from board " << board_idx << " (slot " << fboard_id_list[board_idx] << ") channel " << ch << " (ch_id " << channel_id << ")" << std::endl;
  prod_wvfms.push_back(waveform);
}

/**
 * @brief Saves debug waveform as a histogram for debugging purposes, based on the waveform data provided. 
 * The histogram is labeled with event number, fragment, board slot ID and channel information.
 *
 * @param[in] board_idx The board index (position in the list of boards).
 * @param[in] fragment_idx The fragment index (order in decoding).
 * @param[in] ch The channel index (channel number from which the waveform is extracted).
 * @param[in] ini_wvfm_timestamp The initial timestamp of the waveform in microseconds.
 * @param[in] end_wvfm_timestamp The final timestamp of the waveform in microseconds.
 * @param[in] wvfms A 2D vector containing the waveforms. Each inner vector corresponds to a channel,
 * and each element of the inner vector represents a sample (ADC count).
 *
 * @details
 * The function generates a histogram with:
 * - X-axis representing time in microseconds, scaled between `TTT_ini_us` and `TTT_end_us`.
 * - Y-axis representing ADC counts.
 * The histogram is stored using ROOT's `TFileService`.
 * 
 * @pre ini_wvfm_timestamp and end_wvfm_timestamp are asumed to be given in microseconds.
 *
 * @see TH1I, TFileService
 */

void sbndaq::SBNDXARAPUCADecoder::save_debug_wvfm(size_t board_idx, size_t fragment_idx, int ch, double ini_wvfm_timestamp, double end_wvfm_timestamp, const std::vector <std::vector <uint16_t> > & wvfms) {
  std::stringstream hist_name("");
  hist_name << "Event " << fevent_counter << " CH " << ch << " [frag " << fragment_idx << ", board " << board_idx << " (slot " << fboard_id_list[board_idx] << ") ] waveform";
  TH1I* hist = tfs->make<TH1I>(hist_name.str().c_str(), hist_name.str().c_str(), wvfms[ch].size(), ini_wvfm_timestamp, end_wvfm_timestamp);
  hist->GetYaxis()->SetTitle("ADCs");
  hist->GetXaxis()->SetTitle("Time [us]");
  for (size_t i = 0; i < wvfms[ch].size(); i++) {
    hist->SetBinContent(i+1, wvfms[ch][i]);
  }
}

/**
 * @brief Extract a sample from a 64-bit buffer using the specified bit positions.
 *
 * @param[in] buffer An unsigned 64-bit integer which represents a temporal buffer for the read words and where the samples are extracted from. 
 * @param[in] msb An unsigned 32-bit integer representing the most significative bit (MSB) where the readout from the buffer paramter.
 * @param[in] lsb An unsigned 32-bit integer representing the less significative bit (LSB) from we end read
 *
 * @details The function shifts the buffer to the right by the number of positions specified by `lsb` so that the least significant bit of the 
 * sample aligns with bit 0. It then applies a mask to isolate the bits between `lsb` and `msb`, inclusive.
 *
 * @return The extracted sample as a 16-bit unsigned integer.
 */
uint16_t sbndaq::SBNDXARAPUCADecoder::get_sample(uint64_t buffer, uint32_t msb, uint32_t lsb) {
  uint64_t mask = (1U << (msb - lsb + 1)) - 1;
  uint64_t sample = buffer >> lsb;
  return sample & mask;
}

/**
 * @brief Read a 32-bit word from the data pointer and advances the pointer. 
 *
 * @param[in, out] data_ptr A reference to a pointer pointing to the current position in the data.
 *
 * @details This function retrieves a 32-bit word from the memory location pointed to by `data_ptr`. After reading, it advances `data_ptr` to 
 * the next 32-bit word location.
 *
 * @return The 32-bit word read from the location pointed to by `data_ptr`.
 */
uint32_t sbndaq::SBNDXARAPUCADecoder::read_word(const uint32_t* & data_ptr) {
  uint32_t word = *data_ptr;
  data_ptr += 1;
  return word;
}

/**
 * @brief Generates a unique global channel identifier using the board slot and the channel number of that board.
 *
 * @param[in] board Index of the board in `fboard_id_list` from which to derive the board slot.
 * @param[in] board_channel The specific channel number on the given board.
 *
 * @details This function computes a `channel_id` by combining the board slot and the specific 
 * channel number on that board. 
 * The unique identifier `channel_id` (\f$ CH_{ID} $\f) is computed as follows:
 * \f[
 * CH\_{ID} = B\_{ID} \times 100 + CH\_{B}
 * \f]
 *
 * Where:
 * - \f$ B\_{ID} \f$ is the fragment ID retrieved from `fboard_id_list` based on the slot.
 * - \f$ CH\_B \f$ is the channel number on that board.
 *
 * @return A unique identifier for the specified channel as an unsigned integer.
 */
unsigned int sbndaq::SBNDXARAPUCADecoder::get_channel_id(unsigned int board, unsigned int board_channel) {
  unsigned int channel_id = fboard_id_list[board] * 100 + board_channel;
  return channel_id;
}

DEFINE_ART_MODULE(sbndaq::SBNDXARAPUCADecoder)
