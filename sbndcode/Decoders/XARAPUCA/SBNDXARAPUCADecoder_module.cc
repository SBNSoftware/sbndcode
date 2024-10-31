////////////////////////////////////////////////////////////////////////
// Class:       SBNDXARAPUCADecoder
// Plugin Type: producer
// File:        SBNDXARAPUCADecoder_module.cc
// Author       Alicia Vázquez-Ramos (aliciavr@ugr.es)
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
#include <algorithm>

#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1740Fragment.hh"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "art_root_io/TFileService.h"
#include "TH1D.h"

namespace sbndaq {
  class SBNDXARAPUCADecoder;
}


class sbndaq::SBNDXARAPUCADecoder : public art::EDProducer {
public:
  explicit SBNDXARAPUCADecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDXARAPUCADecoder(SBNDXARAPUCADecoder const&) = delete;
  SBNDXARAPUCADecoder(SBNDXARAPUCADecoder&&) = delete;
  SBNDXARAPUCADecoder& operator=(SBNDXARAPUCADecoder const&) = delete;
  SBNDXARAPUCADecoder& operator=(SBNDXARAPUCADecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  uint event_counter;

  uint fnum_caen_boards;
  uint ffragment_id_offset;
  std::vector<uint> fboard_id_list;

  std::string fcaen_module_label;
  std::vector<std::string> fcaen_fragment_name;

  std::string fch_instance_name;

  float fns_per_sample;
  uint16_t fns_per_tick;

  std::stringstream hist_name;
  art::ServiceHandle<art::TFileService> tfs;

  // Functions.
  void add_fragment(artdaq::Fragment& fragment, std::vector <std::vector <artdaq::Fragment> >& fragments);
  uint16_t get_range(uint64_t buffer, uint32_t msb, uint32_t lsb);
  uint32_t read_word(const uint32_t* & data_ptr);
  
};


sbndaq::SBNDXARAPUCADecoder::SBNDXARAPUCADecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  std::cout << "Entering the constructor" << std::endl;
 
  event_counter = 0;
  
  fnum_caen_boards = p.get<uint> ("num_caen_boards");
  ffragment_id_offset = p.get<uint> ("fragment_id_offset");
  
  fboard_id_list = p.get<std::vector <uint> > ("board_id_list");
  
  fcaen_module_label = p.get<std::string> ("caen_module_label");
  fcaen_fragment_name = p.get<std::vector <std::string> > ("caen_fragment_name");
  
  fch_instance_name = p.get<std::string> ("xarapucaInstanceName");

  fns_per_sample = p.get<float>("ns_per_sample");
  fns_per_tick = 8;
  
  std::cout << event_counter << "\t - xarapucaInstanceName: " << fch_instance_name << "\t - caen module label: " << fcaen_module_label << std::endl;
  std::cout << "Num CAEN boards: " << fnum_caen_boards << std::endl;

  // Call appropriate produces<>() functions here.
  produces <std::vector <raw::OpDetWaveform> > (fch_instance_name);
  
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbndaq::SBNDXARAPUCADecoder::produce(art::Event& e)
{
  std::cout << "Entering the produce function" << std::endl;
  // Implementation of required member function here.
  event_counter++;

  // Initialize output product
  auto prod_wvfms = std::make_unique <std::vector <raw::OpDetWaveform> > ();

  // Create a vector for all fragments for <fnum_caen_boards> boards
  std::vector <std::vector <artdaq::Fragment> > fragments (fnum_caen_boards);

  bool found_caen = false;

  for (const std::string &caen_name: fcaen_fragment_name) {
    art::Handle <std::vector <artdaq::Fragment> > fragment_handle;
    e.getByLabel(fcaen_module_label, caen_name, fragment_handle);
    
    if (!fragment_handle.isValid()) {
      std::cout << "\tHandle not valid for " << caen_name << std::endl;
    
    } else if (fragment_handle->size() == 0) {
      std::cout << "\tHandle with size " << fragment_handle->size() << " for " << caen_name << std::endl;
    
    } else {
      found_caen |= true;
      std::cout << "\tHandle valid for " << caen_name << " with size " << fragment_handle->size() << std::endl;
    
      // It is a group of containers containing a group of CAEN V1740 fragments
      if (fragment_handle->front().type() == artdaq::Fragment::ContainerFragmentType){
        std::cout << "\t\tCONTAINER";
        
        // For every container
        for (auto container: *fragment_handle) {
          artdaq::ContainerFragment container_fragment(container);
          
          // Search for all CAEN V1740 fragments inside the container
          if (container_fragment.fragment_type() == sbndaq::detail::FragmentType::CAENV1740) {
            for (size_t i = 0; i < container_fragment.block_count(); i++) {
              std::cout << "\t\t[container_size: " << container_fragment.block_count() << "]" << std::endl;
              artdaq::Fragment fragment = *container_fragment[i].get();
              add_fragment(fragment, fragments);
            }
          }
        }
    
      // It is a group of CAEN V1740 fragments 
      } else if (fragment_handle->front().type() == sbndaq::detail::FragmentType::CAENV1740) {
        std::cout << "\t\tCAENV1740" << std::endl;
        
        // Search for all CAEN V1740 fragments
        for (size_t i = 0; i < fragment_handle->size(); i++) {
          artdaq::Fragment fragment = fragment_handle->at(i);
          add_fragment(fragment, fragments);
        }
      }
    }
    fragment_handle.removeProduct();
  }

  if (!found_caen) {
    std::cout << "\tNo CAEN V1740 fragments of any type found, pushing empty waveforms." << std::endl;
    e.put(std::move(prod_wvfms), fch_instance_name);      
  
  } else {
    std::cout << "\tDecoding fragments... " << std::endl;
    for (size_t b = 0; b < fragments.size(); b++) {
      for (size_t f = 0; f < fragments[b].size(); f++) {
        std::cout << "\t\tCAEN fragment " << f << " from board " << b << std::endl;

        CAENV1740Fragment caen_fragment(fragments[b][f]);
        CAENV1740FragmentMetadata const* metadata = caen_fragment.Metadata();
        uint32_t num_channels = metadata->nChannels;
        uint32_t num_samples_per_wvfm = metadata->nSamples;
        uint32_t num_samples_per_group = num_samples_per_wvfm * 8;
        std::cout << "\t\t\tNumber of channels: " << num_channels << std::endl;
        std::cout << "\t\t\tNumber of samples: " << num_samples_per_wvfm << "(" << num_samples_per_group << " samples per group - 8 channels -)"<< std::endl;

        CAENV1740Event const* event = caen_fragment.Event();
        CAENV1740EventHeader header = event->Header;
        float TTT_end_ns = header.triggerTime() * fns_per_tick;
        float TTT_ini_ns = TTT_end_ns - num_samples_per_wvfm * fns_per_sample;
        float TTT_end_us = TTT_end_ns * 1E-3f;
        float TTT_ini_us = TTT_ini_ns * 1E-3f;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "TTT_ini [us]: " << TTT_ini_us << "\tTTT_end [us]: " << TTT_end_us << std::endl;
          
        //std::cout << "\t\t\tBoard ID (from header): " << header.boardID << std::endl;

        uint32_t num_words_per_event = header.eventSize;
        uint32_t num_words_per_header = sizeof(CAENV1740EventHeader)/sizeof(uint32_t);
        uint32_t num_words_per_wvfm = (num_words_per_event - num_words_per_header);

        std::cout << "\t\t\tNumber of words per event: " << num_words_per_event << " [Header: " << num_words_per_header << ", Waveform: " << num_words_per_wvfm << "] words" << std::endl;

        // ===============  Start decoding the waveforms =============== //
        std::vector <std::vector <uint16_t> > wvfms(num_channels, std::vector<uint16_t>(num_samples_per_wvfm, 0));

        // Absolute sample number        
        uint32_t S = 0;
        // Buffer variables
        uint64_t buffer = 0;
        uint32_t bits_in_buffer = 0;

        uint32_t sample_bits = 12;
        
        // Data pointer to the beggining of the waveforms stores in the event.
        const uint32_t* data_ptr = reinterpret_cast<const uint32_t*>(fragments[b][f].dataBeginBytes() + sizeof(CAENV1740EventHeader));
        for (size_t j = 0; j < num_words_per_wvfm; j++) {
          uint64_t word = read_word(data_ptr);
          //std::cout << buffer << "[word: " << word << "]" << std::endl;
          buffer |= word << bits_in_buffer;
          bits_in_buffer += sizeof(uint32_t) * 8;
          //std::cout << "  +" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;

          // Get sequences of 12 bits
          while (bits_in_buffer >= sample_bits) {
            uint32_t g = (S / num_samples_per_group);                       // Group index.
            uint32_t c = ((S / 3) % 8) + g * 8;                             // Channel index.
            uint32_t s = (S % 3) + ((S / 24) * 3) % num_samples_per_wvfm;   // Sample/channel index.
            uint16_t sample = get_range(buffer, sample_bits - 1, 0);
            //std::cout << "\tSample: " << sample << "\tg: " << g <<  "\tch: " << c << "\ts:" << s << "\tS: " << S << std::endl;
            wvfms[c][s] = sample;
            buffer >>= sample_bits;
            bits_in_buffer -= sample_bits;
            //std::cout << "  -" << buffer << " [bits in buffer: "<< bits_in_buffer << "]" << std::endl;
            S++;
          }
        }
        
        for (size_t ch = 0; ch < wvfms.size(); ch++) {
          raw::OpDetWaveform waveform(TTT_ini_us, ch, wvfms[ch]);
          prod_wvfms->push_back(waveform);
          hist_name.str("");
          hist_name << "Event " << event_counter << " CH " << ch << " [frag: " << f << ", board: " << b << "] waveform";
          TH1D* hist = tfs->make<TH1D>(hist_name.str().c_str(), hist_name.str().c_str(), wvfms[ch].size(), TTT_ini_us, TTT_end_us);
          hist->GetYaxis()->SetTitle("ADCs");
          hist->GetXaxis()->SetTitle("Time [us]");
          for (size_t i = 0; i < wvfms[ch].size(); i++) {
            hist->SetBinContent(i+1, wvfms[ch][i]);
          }
        }
      }
    }

    e.put(std::move(prod_wvfms), fch_instance_name);
  }

  std::cout << "End of producer function." << std::endl;
}

void sbndaq::SBNDXARAPUCADecoder::add_fragment(artdaq::Fragment& fragment, std::vector <std::vector <artdaq::Fragment> >& fragments) {
  auto fragment_id = fragment.fragmentID() - ffragment_id_offset;
  auto it = std::find(fboard_id_list.begin(), fboard_id_list.end(), fragment_id);

  if (it != fboard_id_list.end()) {
    uint index = it - fboard_id_list.begin();
    if (index < 0 || index >= fnum_caen_boards) {
      std::cout << "\t\t\tFragment ID " << fragment_id << " (" << index << ") is out of range. Skipping this fragment..." << std::endl;
    } else {
      std::cout << "\t\t\tGetting a CAENV1740 fragment: [" << fragment_id << " (" << index << ") " << "]" << fragment;
      fragments.at(index).push_back(fragment);
    }
  } else {
      std::cout << "\t\t\tFragment ID " << fragment_id << " is not valid. Skipping this fragment..." << std::endl;
  }

}

uint16_t sbndaq::SBNDXARAPUCADecoder::get_range(uint64_t buffer, uint32_t msb, uint32_t lsb) {
  uint64_t mask = (1U << (msb - lsb + 1)) - 1;
  uint64_t sample = buffer >> lsb;
  return sample & mask;
}

uint32_t sbndaq::SBNDXARAPUCADecoder::read_word(const uint32_t* & data_ptr) {
  uint32_t word = *data_ptr;
  data_ptr += 1;
  return word;
}

DEFINE_ART_MODULE(sbndaq::SBNDXARAPUCADecoder)
