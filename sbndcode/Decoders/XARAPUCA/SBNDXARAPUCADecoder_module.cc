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

};


sbndaq::SBNDXARAPUCADecoder::SBNDXARAPUCADecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  std::cout << "Entering the constructor" << std::endl;
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbndaq::SBNDXARAPUCADecoder::produce(art::Event& e)
{
  std::cout << "Entering the produce function" << std::endl;
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(sbndaq::SBNDXARAPUCADecoder)
