////////////////////////////////////////////////////////////////////////
// Class:       SimpleBeamRateCalibAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        SimpleBeamRateCalibAnalyzer_module.cc
//
// Generated at Mon Aug 11 12:51:33 2025 by Jacob McLaughlin using cetskelgen
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

class SimpleBeamRateCalibAnalyzer;


class SimpleBeamRateCalibAnalyzer : public art::EDAnalyzer {
public:
  explicit SimpleBeamRateCalibAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleBeamRateCalibAnalyzer(SimpleBeamRateCalibAnalyzer const&) = delete;
  SimpleBeamRateCalibAnalyzer(SimpleBeamRateCalibAnalyzer&&) = delete;
  SimpleBeamRateCalibAnalyzer& operator=(SimpleBeamRateCalibAnalyzer const&) = delete;
  SimpleBeamRateCalibAnalyzer& operator=(SimpleBeamRateCalibAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


SimpleBeamRateCalibAnalyzer::SimpleBeamRateCalibAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SimpleBeamRateCalibAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(SimpleBeamRateCalibAnalyzer)
