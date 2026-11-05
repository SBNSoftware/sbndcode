////////////////////////////////////////////////////////////////////////
// Class:       WaveDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        WaveDump_module.cc
//
// Generated at Thu Nov  5 13:45:51 2026 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
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
  class WaveDump;
}


class sbnd::WaveDump : public art::EDAnalyzer {
public:
  explicit WaveDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveDump(WaveDump const&) = delete;
  WaveDump(WaveDump&&) = delete;
  WaveDump& operator=(WaveDump const&) = delete;
  WaveDump& operator=(WaveDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


sbnd::WaveDump::WaveDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::WaveDump::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(sbnd::WaveDump)
