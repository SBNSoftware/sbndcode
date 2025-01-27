////////////////////////////////////////////////////////////////////////
// Class:       michelETagger
// Plugin Type: filter (Unknown Unknown)
// File:        michelETagger_module.cc
//
// Generated at Mon Jan 27 14:29:41 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

class michelETagger;


class michelETagger : public art::EDFilter {
public:
  explicit michelETagger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  michelETagger(michelETagger const&) = delete;
  michelETagger(michelETagger&&) = delete;
  michelETagger& operator=(michelETagger const&) = delete;
  michelETagger& operator=(michelETagger&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.

};


michelETagger::michelETagger(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool michelETagger::filter(art::Event& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(michelETagger)
