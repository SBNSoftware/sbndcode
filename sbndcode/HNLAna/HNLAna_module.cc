////////////////////////////////////////////////////////////////////////
// Class:       HNLAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        HNLAna_module.cc
//
// Generated at Wed Jun 29 11:37:04 2022 by Vu Chi Lan Nguyen using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace sbnd {
  class HNLAna;
}


class sbnd::HNLAna : public art::EDAnalyzer {
public:
  explicit HNLAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HNLAna(HNLAna const&) = delete;
  HNLAna(HNLAna&&) = delete;
  HNLAna& operator=(HNLAna const&) = delete;
  HNLAna& operator=(HNLAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


sbnd::HNLAna::HNLAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::HNLAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void sbnd::HNLAna::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::HNLAna::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::HNLAna)
