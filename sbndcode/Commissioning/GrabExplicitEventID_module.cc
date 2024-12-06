////////////////////////////////////////////////////////////////////////
// Class:       GrabExplicitEventID
// Plugin Type: filter (Unknown Unknown)
// File:        GrabExplicitEventID_module.cc
//
// Generated at Wed Oct 16 22:10:13 2024 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
#ifndef GrabExplicitEventID_H
#define GrabExplicitEventID_H
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

namespace TriggerWork {
  class GrabExplicitEventID;
}


class TriggerWork::GrabExplicitEventID : public art::EDFilter {
public:
  explicit GrabExplicitEventID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GrabExplicitEventID(GrabExplicitEventID const&) = delete;
  GrabExplicitEventID(GrabExplicitEventID&&) = delete;
  GrabExplicitEventID& operator=(GrabExplicitEventID const&) = delete;
  GrabExplicitEventID& operator=(GrabExplicitEventID&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  std::vector<int> fEventIDToGrab;

private:

  // Declare member data here.

};


TriggerWork::GrabExplicitEventID::GrabExplicitEventID(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fEventIDToGrab = p.get<std::vector<int>>("EventIDs");
}

bool TriggerWork::GrabExplicitEventID::filter(art::Event& e)
{
  // Implementation of required member function here.
  int fEvNumber = e.id().event();
  bool GrabThisEvent = std::any_of(fEventIDToGrab.begin(), fEventIDToGrab.end(), [fEvNumber](int x){return ( fEvNumber == x); } );
  std::cout << fEvNumber << "  " << GrabThisEvent << std::endl;
  return GrabThisEvent;
}

void TriggerWork::GrabExplicitEventID::beginJob()
{
  // Implementation of optional member function here.
}

void TriggerWork::GrabExplicitEventID::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerWork::GrabExplicitEventID)
#endif