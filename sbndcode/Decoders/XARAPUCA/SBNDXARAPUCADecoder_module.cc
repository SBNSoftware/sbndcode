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
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1740Fragment.hh"

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

  unsigned int fevent_counter; /**< Event counter. */

  unsigned int fnum_caen_boards; /**< Maximum number of CAEN boards to be considered. */
  unsigned int ffragment_id_offset; /**< Offset value to get the slot of every CAEN board. */
  std::vector<unsigned int> fboard_id_list; /**< Possible CAEN board identifiers based on the slot they are installed. */

  std::string fcaen_module_label; /**< Label identifying the module where the CAEN fragments come from. */
  std::vector<std::string> fcaen_fragment_names; /**< Fragments accepted by this decoder. */

  std::string fproduct_instance_name; /**< Name assigned to the product instance generated by this art::EDProducer. */

  unsigned int fns_per_sample; /**< Number of nanoseconds per sample. */
  uint16_t fns_per_tick; /**< Number of nanoseconds per trigger time stamp (TTT) tick. */

  uint32_t fsample_bits; /**< Number of bits per sample. */

  art::ServiceHandle<art::TFileService> tfs; /**<  ServiceHandle object to store the histograms in the decoder_hist.root output file. */
  int fstore_debug_waveforms; /**< Number of waveforms to store in the ServiceHandle object for debugging purposes (0: none, -1: all, n: first n waveforms). */

  bool fdebug_all; /**< If `true` all debug information is printed. */
  bool fdebug_handle; /**< If `true` art::Handle object data is printed. */
  bool fdebug_timing; /**< If `true` timing data is printed. */
  bool fdebug_buffer; /**< If `true` the buffer status is printed. */
  bool fdebug_waveforms; /**< If `true` waveforms decoding data is printed. */
  bool fverbose; /**< If `true` it increases verbosity of console output for detailed processing steps. */

  // Class methods.
  void save_prod_wvfm(size_t b, size_t ch, float TTT_ini_us, const std::vector <std::vector <uint16_t> > & wvfms, std::vector <raw::OpDetWaveform> & prod_wvfms);
  void save_debug_wvfm(size_t b, size_t f, int ch, float TTT_ini_us, float TTT_end_us, const std::vector <std::vector <uint16_t> > & wvfms);
  void add_fragment(artdaq::Fragment& fragment, std::vector <std::vector <artdaq::Fragment> >& fragments);
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
  fproduct_instance_name = p.get<std::string> ("product_instance_name", "XARAPUCAChannels");

  // Gets timing information.
  fns_per_sample = p.get<unsigned int> ("ns_per_sample", 16);
  fns_per_tick = 8;

  // Sets the number of bits per sample.
  fsample_bits = 12;
  
  // Gets the number of waveforms to store in the debug output file.
  fstore_debug_waveforms = p.get<int> ("store_debug_waveforms", 0);

  // Gets the debug and verbose options.
  fdebug_all = p.get<bool> ("debug_all", false);
  fdebug_handle = p.get<bool> ("debug_handle", false);
  fdebug_timing = p.get<bool> ("debug_timing", false);
  fdebug_buffer = p.get<bool> ("debug_buffer", false);
  fdebug_waveforms = p.get<bool> ("debug_waveforms", false);
  fverbose = p.get<bool> ("verbose", false);

  // Creates the instance product of this module.
  produces <std::vector <raw::OpDetWaveform> > (fproduct_instance_name);
}

/**
 * @brief Main function of the art::EDProducer: the produce function analyzes every event producing waveforms after the decoding.
 * @param[in] e The event to be processed.
 * @details It is the main function of the art::EDProducer module:
 * 1. Iterates over `fcaen_fragment_names` to retrieve each fragment by label. It checks if the fragment handle is valid and 
 * non-empty. If valid, it updates `found_caen` and processes each fragment based on its type: either as a 
 * container of fragments or as individual CAEN V1740 fragments.
 * 2. Access the data inside the CAEN fragments and extracts header, timing data and the raw waveform.
 * 3. Binary decoding of the raw waveform. To assign efficiently each sequential sample extracted. These indices formulas are 
 * applied:
 * - The board channel index:
 * \f[
 * c = \left( \frac{S}{3} \mod 8 \right) + g \times 8
 * \f]
 * - The channel sample index: 
 * \f[
 * s = (S \mod 3) + \left( \frac{S}{24} \times 3 \right) \mod s_{w}}
 * \f]
 * Where the group index is computed as \f$ \frac{S}{s_{g}} \f$.
 * 4. The decoded waveforms are dumped into:
 * - decoder_hist.root file (waveforms in ROOT histograms format).
 * - xarapucadecoder-art.root file (EDProducer product: a vector of raw::OpDetWaverform).
 */
void sbndaq::SBNDXARAPUCADecoder::produce(art::Event& e)
{
  if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: entering the produce function." << std::endl;

  // Advances the event counter.
  fevent_counter++;

  // Initializes the output instance product.
  auto prod_wvfms = std::make_unique <std::vector <raw::OpDetWaveform> > ();
  if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: products initialized." << std::endl;

  // Creates a vector for all fragments for <fnum_caen_boards> boards
  std::vector <std::vector <artdaq::Fragment> > fragments (fnum_caen_boards);

  // Flag to track if valid CAEN fragments are found.
  bool found_caen = false;

  // This loop iterates over the valid CAEN fragments names (`fcaen_fragment_names`) to check 
  // if the art::Handle object is not empty and has valid CAEN fragments. After this, it extracts
  // the CAEN V1740 fragments into a nested std::vector, assigning each fragment to its board.
  if (fverbose | fdebug_handle | fdebug_all) std::cout << "\n > SBNDXARAPUCADecoder::produce: searching for V1740 fragments." << std::endl;
  for (const std::string &caen_name: fcaen_fragment_names) {
    art::Handle <std::vector <artdaq::Fragment> > fragment_handle;
    e.getByLabel(fcaen_module_label, caen_name, fragment_handle);
    
    // The art::Handle object is not valid.
    if (!fragment_handle.isValid()) { 
      if (fdebug_handle | fdebug_all) std::cout << "\tHandle not valid for " << caen_name << "." << std::endl;
    
    // The art::Handle object is empty.
    } else if (fragment_handle->size() == 0) { 
      if (fdebug_handle | fdebug_all) std::cout << "\tHandle with size " << fragment_handle->size() << " for " << caen_name << "."<< std::endl;
    
    // The art::Handle object is valid and not empty.
    } else { 
      found_caen |= true;
      size_t frag_handle_size = fragment_handle->size();

      if (fdebug_handle | fdebug_all) std::cout << "\tHandle valid for " << caen_name << " (" << frag_handle_size << " objects)." << std::endl;
    
      // It is a group of containers containing a group of CAEN V1740 fragments.
      if (fragment_handle->front().type() == artdaq::Fragment::ContainerFragmentType){
        
        // For every container,
        for (size_t c = 0; c < frag_handle_size; c++) {
          auto container = fragment_handle->at(c);
          artdaq::ContainerFragment container_fragment(container);
          
          // it searches for all CAEN V1740 fragments inside the container if it contains CAEN V1740 fragments.
          if (container_fragment.fragment_type() == sbndaq::detail::FragmentType::CAENV1740) {
            size_t num_caen_fragments = container_fragment.block_count(); 
            if (fdebug_handle | fdebug_all) std::cout << "\t\tContainerCAENV1740 "<< c << " - CAENV1740 fragments found: " << num_caen_fragments << "." << std::endl;
            
            for (size_t f = 0; f < num_caen_fragments; f++) {
              artdaq::Fragment fragment = *container_fragment[f].get();
              add_fragment(fragment, fragments);
            } // End CAEN V1740 fragments loop.
          }
        } // End Container fragments loop.
    
      // It is a group of CAEN V1740 fragments.
      } else if (fragment_handle->front().type() == sbndaq::detail::FragmentType::CAENV1740) {
        if (fdebug_handle | fdebug_all) std::cout << "\t\tCAENV1740 fragments found: " << frag_handle_size << "." << std::endl;
        
        // It searches for all CAEN V1740 fragments.
        for (size_t f = 0; f < frag_handle_size; f++) {
          artdaq::Fragment fragment = fragment_handle->at(f);
          add_fragment(fragment, fragments);
        } // End CAEN V1740 fragments loop.
      }
    } // End extracting CAEN V1740 fragments.
    fragment_handle.removeProduct();
  } // End CAEN fragment names loop.

  // If no valid CAEN fragments are found it pushes empty waveforms,
  if (!found_caen) {
    e.put(std::move(prod_wvfms), fproduct_instance_name);      
    if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: no CAEN V1740 fragments of any type found, pushed empty waveforms." << std::endl;

  // otherwise it starts the decoding of the fragments found.
  } else {
    if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: decoding V1740 CAEN fragments..." << std::endl;
    for (size_t b = 0; b < fragments.size(); b++) {
      for (size_t f = 0; f < fragments[b].size(); f++) {
        if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: decoding V1740 CAEN fragment " << f << " from the board " << b << " (slot " << fboard_id_list[b] << "):" << std::endl;

        // Accesses the metadata of the CAEN fragment.
        CAENV1740Fragment caen_fragment(fragments[b][f]);
        CAENV1740FragmentMetadata const* metadata = caen_fragment.Metadata();
        uint32_t num_channels = metadata->nChannels;
        uint32_t num_samples_per_wvfm = metadata->nSamples;
        uint32_t num_samples_per_group = num_samples_per_wvfm * 8;
        if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::produce: number of channels: " << num_channels << std::endl;
        if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::produce: number of samples: " << num_samples_per_wvfm << " (" << num_samples_per_group << " samples per group - 8 channels -)"<< std::endl;

        // Accesses the event and header data of the CAEN fragment.
        CAENV1740Event const* event = caen_fragment.Event();
        CAENV1740EventHeader header = event->Header;

        // Gets the TTT data of the event.
        uint32_t TTT_end_ticks = header.triggerTime();
        float TTT_end_ns = TTT_end_ticks * fns_per_tick;
        float TTT_ini_ns = TTT_end_ns - num_samples_per_wvfm * fns_per_sample;
        float TTT_end_us = TTT_end_ns * 1E-3f;
        float TTT_ini_us = TTT_ini_ns * 1E-3f;
        if (fverbose | fdebug_timing | fdebug_all) {
          std::cout << std::fixed << std::setprecision(1);
          std::cout << "  > SBNDXARAPUCADecoder::produce: time window of " << TTT_end_us - TTT_ini_us << " us: [" << TTT_ini_us << ", " << TTT_end_us << "] us." << std::endl;
        }
        if (fdebug_timing | fdebug_all) {
          std::cout << "\tTTT_end [ticks]: " << TTT_end_ticks << "\tNanoseconds per tick: " << fns_per_tick << " ns." << std::endl;
          std::cout << "\tTTT_end [ns]: " << TTT_end_ns << "\tNanoseconds per sample: " << fns_per_sample << " ns." << std::endl;
          std::cout << "\tTTT_ini [us]: " << TTT_ini_us << std::endl;
          std::cout << "\tTTT_end [us]: " << TTT_end_us << std::endl;
        } 
          
        // Gets the number of words of the header and the waveforms.
        uint32_t num_words_per_event = header.eventSize;
        uint32_t num_words_per_header = sizeof(CAENV1740EventHeader)/sizeof(uint32_t);
        uint32_t num_words_per_wvfm = (num_words_per_event - num_words_per_header);
        if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::produce: number of words per event: " << num_words_per_event << " (Header: " << num_words_per_header << ", Waveform: " << num_words_per_wvfm << ") words." << std::endl;

        // ===============  Start decoding the waveforms =============== //
        if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::produce: binary decoding of the waveforms starting... " << std::endl;
        
        std::vector <std::vector <uint16_t> > wvfms(num_channels, std::vector<uint16_t>(num_samples_per_wvfm, 0));

        // Absolute sample number [0, TOTAL_NUM_SAMPLES] where TOTAL_NUM_SAMPLES is the total number of samples stored for an event.
        uint32_t S = 0;
        // Buffer variables.
        uint64_t buffer = 0;
        uint32_t bits_in_buffer = 0;

        // Data pointer to the beggining of the waveforms stores in the event.
        const uint32_t* data_ptr = reinterpret_cast<const uint32_t*>(fragments[b][f].dataBeginBytes() + sizeof(CAENV1740EventHeader));
        // Accesses each word, stores it in the buffer and then the samples are extracted from the buffer.
        for (size_t j = 0; j < num_words_per_wvfm; j++) {
          uint64_t word = read_word(data_ptr);

          // Adds the new word to the buffer and increments the number of bits stored in it.
          if (fdebug_buffer | fdebug_all) std::cout << buffer << "[word: " << word << "]" << std::endl;
          buffer |= word << bits_in_buffer;
          bits_in_buffer += sizeof(uint32_t) * 8; // bytes * 8 bits/byte
          if (fdebug_buffer | fdebug_all) std::cout << "  +" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;

          // Obtains 12-bit sequences from the buffer and assigns each sample to the channel and channel sample it belongs to. 
          while (bits_in_buffer >= fsample_bits) {
            // Computes board channel, channel sample and group channel and assigns the sample to those indices.
            uint32_t g = (S / num_samples_per_group);                       // Group index.
            uint32_t c = ((S / 3) % 8) + g * 8;                             // Channel index.
            uint32_t s = (S % 3) + ((S / 24) * 3) % num_samples_per_wvfm;   // Sample/channel index.
            uint16_t sample = get_sample(buffer, fsample_bits - 1, 0);
            wvfms[c][s] = sample;
            if (fdebug_waveforms | fdebug_all) std::cout << "\tSample: " << sample << "\tg: " << g <<  "\tch: " << c << "\ts:" << s << "\tS: " << S << std::endl;
            
            // Updates the buffer status removing the read bits and decreasing the number of bits stored in it.
            buffer >>= fsample_bits;
            bits_in_buffer -= fsample_bits;
            if (fdebug_buffer | fdebug_all) std::cout << "  -" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;
            
            // Increments the absolute sample step.
            S++;
          }
        }
        
        // The decoded waveforms are dumped into two products:
        // - A xarapucadecoder-art.root file with the OpDetWaveforms as the product of this producer for further analysis.
        // - A decoder_hist.root file gathering a waveform histograms.
        if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::produce: binary decoding complete, dumping products..." << std::endl;
        
        uint32_t num_debug_wvfms;

        if (fstore_debug_waveforms == -1) {
          num_debug_wvfms = num_channels;
        } else {
          num_debug_wvfms = std::min<size_t>(num_channels, fstore_debug_waveforms);
        }

        uint32_t ch;

        for (ch = 0; ch < num_debug_wvfms; ch++) {
          save_prod_wvfm(b, ch, TTT_ini_us, wvfms, *prod_wvfms);
          save_debug_wvfm(b, f, ch, TTT_ini_us, TTT_end_us, wvfms);
        }

        for (;ch < num_channels; ch++) {
          save_prod_wvfm(b, ch, TTT_ini_us, wvfms, *prod_wvfms);
        }

      }
    }

    e.put(std::move(prod_wvfms), fproduct_instance_name);
    if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: waveforms for all valid fragments dumped." << std::endl;
  }

  if (fverbose) std::cout << "\n > SBNDXARAPUCADecoder::produce: end of the producer function.\n\n" << std::endl;
}

/**
 * @brief Converts a production waveform from raw ADC data into a `raw::OpDetWaveform` object,
 * assigning it a global channel ID and timestamp, and appends it to the output collection.
 *
 * @param[in] b The board index (position in the list of boards).
 * @param[in] ch The channel index (channel number from which the waveform is extracted).
 * @param[in] TTT_ini_us The initial timestamp of the waveform in microseconds.
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
 * @see raw::OpDetWaveform
 */
void sbndaq::SBNDXARAPUCADecoder::save_prod_wvfm(size_t b, size_t ch, float TTT_ini_us, const std::vector <std::vector <uint16_t> > & wvfms, std::vector <raw::OpDetWaveform> & prod_wvfms) {
  unsigned int channel_id = get_channel_id(b, ch);
  raw::OpDetWaveform waveform(TTT_ini_us, channel_id, wvfms[ch]);
  if (fdebug_waveforms | fdebug_all) std::cout << "Pushing waveform from board " << b << " (slot " << fboard_id_list[b] << ") channel " << ch << " (ch_id " << channel_id << ")" << std::endl;
  prod_wvfms.push_back(waveform);
}

/**
 * @brief Saves debug waveform as a histogram for debugging purposes, based on the waveform data provided. 
 * The histogram is labeled with event number, fragment, board slot ID and channel information.
 *
 * @param[in] b The board index (position in the list of boards).
 * @param[in] f The fragment index.
 * @param[in] ch The channel index (channel number from which the waveform is extracted).
 * @param[in] TTT_ini_us The initial timestamp of the waveform in microseconds.
 * @param[in] TTT_end_us The final timestamp of the waveform in microseconds.
 * @param[in] wvfms A 2D vector containing the waveforms. Each inner vector corresponds to a channel,
 * and each element of the inner vector represents a sample (ADC count).
 *
 * @details
 * The function generates a histogram with:
 * - X-axis representing time in microseconds, scaled between `TTT_ini_us` and `TTT_end_us`.
 * - Y-axis representing ADC counts.
 * The histogram is stored using ROOT's `TFileService`.
 *
 * @see TH1I, TFileService
 */

void sbndaq::SBNDXARAPUCADecoder::save_debug_wvfm(size_t b, size_t f, int ch, float TTT_ini_us, float TTT_end_us, const std::vector <std::vector <uint16_t> > & wvfms) {
  std::stringstream hist_name("");
  hist_name << "Event " << fevent_counter << " CH " << ch << " [frag " << f << ", board " << b << " (slot " << fboard_id_list[b] << ") ] waveform";
  TH1I* hist = tfs->make<TH1I>(hist_name.str().c_str(), hist_name.str().c_str(), wvfms[ch].size(), TTT_ini_us, TTT_end_us);
  hist->GetYaxis()->SetTitle("ADCs");
  hist->GetXaxis()->SetTitle("Time [us]");
  for (size_t i = 0; i < wvfms[ch].size(); i++) {
    hist->SetBinContent(i+1, wvfms[ch][i]);
  }
}

/**
 * @brief Add the fragment given as a parameter to the nested vector of fragments to be processed if it satisfies the ID restrictions.
 *
 * @param[in] fragment artdaq::Fragment object containing information related to a CAENV1740 fragment.
 * @param[out] fragments A nested std::vector of artdaq::Fragment objects to be decoded.
 *
 * @details The function check the validity of the fragment ID after subtracting the offset (`ffragment_id_offset`).
 * It searches for this adjusted ID in `fboard_id_list`. If found, it calculates the insertion index.
 * - If the index is within the valid range `[0, fnum_caen_boards)`, the fragment is added to `fragments` at
 *   the specified position (this position represents the board to wich it belongs to).
 * - If the index is out of range, or if the fragment ID is not found in `fboard_id_list`, the fragment is
 *   skipped.
 */
void sbndaq::SBNDXARAPUCADecoder::add_fragment(artdaq::Fragment& fragment, std::vector <std::vector <artdaq::Fragment> >& fragments) {
  auto fragment_id = fragment.fragmentID() - ffragment_id_offset;
  auto it = std::find(fboard_id_list.begin(), fboard_id_list.end(), fragment_id);

  if (it != fboard_id_list.end()) {
    unsigned int index = it - fboard_id_list.begin();
    if (index < 0 || index >= fnum_caen_boards) {
      if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::add_fragment: fragment ID " << fragment_id << " (" << index << ") is out of range. Skipping this fragment..." << std::endl;
    } else {
      if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::add_fragment: getting a CAENV1740 fragment: [" << fragment_id << " (" << index << ")" << "]" << fragment;
      fragments.at(index).push_back(fragment);
    }
  } else {
      if (fverbose) std::cout << "  > SBNDXARAPUCADecoder::add_fragment: fragment ID " << fragment_id << " is not valid. Skipping this fragment..." << std::endl;
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
