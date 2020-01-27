////////////////////////////////////////////////////////////////////////
// Class:       MichelElec
// Plugin Type: analyzer (art v3_04_00)
// File:        MichelElec_module.cc
//
// Generated at Mon Jan 27 10:46:07 2020 by Georgia Chisnall using cetskelgen
// from cetlib version v3_09_00.
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

namespace trigger {
  class MichelElec;
}


class trigger::MichelElec : public art::EDAnalyzer {
public:
  explicit MichelElec(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelElec(MichelElec const&) = delete;
  MichelElec(MichelElec&&) = delete;
  MichelElec& operator=(MichelElec const&) = delete;
  MichelElec& operator=(MichelElec&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


trigger::MichelElec::MichelElec(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void trigger::MichelElec::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void trigger::MichelElec::beginJob()
{
  // Implementation of optional member function here.
}

void trigger::MichelElec::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(trigger::MichelElec)
