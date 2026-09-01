////////////////////////////////////////////////////////////////////////
// Class:       ArtDAQFragmentBlender
// Plugin Type: producer (Unknown Unknown)
// File:        ArtDAQFragmentBlender_module.cc
//
// Generated at Tue Sep  1 12:05:42 2026 by Jacob McLaughlin using cetskelgen
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

namespace FragmentScramble {
  class ArtDAQFragmentBlender;
}


class FragmentScramble::ArtDAQFragmentBlender : public art::EDProducer {
public:
  explicit ArtDAQFragmentBlender(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ArtDAQFragmentBlender(ArtDAQFragmentBlender const&) = delete;
  ArtDAQFragmentBlender(ArtDAQFragmentBlender&&) = delete;
  ArtDAQFragmentBlender& operator=(ArtDAQFragmentBlender const&) = delete;
  ArtDAQFragmentBlender& operator=(ArtDAQFragmentBlender&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


FragmentScramble::ArtDAQFragmentBlender::ArtDAQFragmentBlender(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void FragmentScramble::ArtDAQFragmentBlender::produce(art::Event& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(FragmentScramble::ArtDAQFragmentBlender)
