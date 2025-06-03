////////////////////////////////////////////////////////////////////////
// Class:       GaussHitFilter
// Plugin Type: producer (Unknown Unknown)
// File:        GaussHitFilter_module.cc
//
// Generated at Tue Jun  3 15:57:00 2025 by Jacob McLaughlin using cetskelgen
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

class GaussHitFilter;


class GaussHitFilter : public art::EDProducer {
public:
  explicit GaussHitFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GaussHitFilter(GaussHitFilter const&) = delete;
  GaussHitFilter(GaussHitFilter&&) = delete;
  GaussHitFilter& operator=(GaussHitFilter const&) = delete;
  GaussHitFilter& operator=(GaussHitFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


GaussHitFilter::GaussHitFilter(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void GaussHitFilter::produce(art::Event& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(GaussHitFilter)
