////////////////////////////////////////////////////////////////////////
// Class:       TimeStampDumper
// Plugin Type: analyzer (Unknown Unknown)
// File:        TimeStampDumper_module.cc
//
// Generated at Tue Jan  7 16:11:21 2025 by Jacob McLaughlin using cetskelgen
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

class TimeStampDumper;


class TimeStampDumper : public art::EDAnalyzer {
public:
  explicit TimeStampDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TimeStampDumper(TimeStampDumper const&) = delete;
  TimeStampDumper(TimeStampDumper&&) = delete;
  TimeStampDumper& operator=(TimeStampDumper const&) = delete;
  TimeStampDumper& operator=(TimeStampDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


TimeStampDumper::TimeStampDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void TimeStampDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(TimeStampDumper)
