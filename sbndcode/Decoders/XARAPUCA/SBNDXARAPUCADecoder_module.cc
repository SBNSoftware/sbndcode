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

#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1740Fragment.hh"

#include "lardataobj/RawData/OpDetWaveform.h"

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

  std::string fcaen_module_label;
  std::vector<std::string> fcaen_fragment_name;

  std::string fch_instance_name;

};


sbndaq::SBNDXARAPUCADecoder::SBNDXARAPUCADecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  std::cout << "Entering the constructor" << std::endl;
 
  event_counter = 0;
  
  fnum_caen_boards = p.get<uint> ("num_caen_boards");
  
  fcaen_module_label = p.get<std::string> ("caen_module_label");
  fcaen_fragment_name = p.get<std::vector <std::string> > ("caen_fragment_name");
  
  fch_instance_name = p.get<std::string> ("xarapucaInstanceName");
  
  std::cout << event_counter << "\t - xarapucaInstanceName: " << fch_instance_name << "\t - caen module label: " << fcaen_module_label << std::endl;
  std::cout << "Num CAEN boards: " << fnum_caen_boards << std::endl;

  // Call appropriate produces<>() functions here.
  produces< std::vector <raw::OpDetWaveform> > (fch_instance_name);
  
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbndaq::SBNDXARAPUCADecoder::produce(art::Event& e)
{
  std::cout << "Entering the produce function" << std::endl;
  // Implementation of required member function here.
  event_counter++;

  // Initialize output product
  auto wvfms = std::make_unique <std::vector <raw::OpDetWaveform> > ();

  // Create a vector for all fragments for <fnum_caen_boards> boards
  std::vector< std::vector < artdaq::Fragment> > fragments (fnum_caen_boards);

  bool found_caen = false;

  for (const std::string &caen_name: fcaen_fragment_name) {
    art::Handle<std::vector <artdaq::Fragment> > fragment_handle;
    e.getByLabel(fcaen_module_label, caen_name, fragment_handle);

    if (!fragment_handle.isValid()) {
      std::cout << "\tHandle not valid for " << caen_name << std::endl;
    } else if (fragment_handle->size() == 0) {
      std::cout << "\tHandle with size " << fragment_handle->size() << " for " << caen_name << std::endl;
    } else {
      found_caen |= true;
      std::cout << "\tHandle valid for " << caen_name << " with size " << fragment_handle->size() << std::endl;
    }
  }

  if (!found_caen) {
    std::cout << "\tNo CAEN V1740 fragments of any type found, pushing empty waveforms." << std::endl;
    e.put(std::move(wvfms), fch_instance_name);      
  } else {
    e.put(std::move(wvfms), fch_instance_name);
  }

  std::cout << std::endl;
}

DEFINE_ART_MODULE(sbndaq::SBNDXARAPUCADecoder)
