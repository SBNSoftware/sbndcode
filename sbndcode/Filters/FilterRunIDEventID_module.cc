////////////////////////////////////////////////////////////////////////
// Class:       FilterRunIDEventID
// Plugin Type: filter (Unknown Unknown)
// File:        FilterRunIDEventID_module.cc
//
// Generated at Mon Oct  5 09:24:12 2026 by Jacob McLaughlin using cetskelgen
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

class FilterRunIDEventID;


class FilterRunIDEventID : public art::EDFilter {
public:
  explicit FilterRunIDEventID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterRunIDEventID(FilterRunIDEventID const&) = delete;
  FilterRunIDEventID(FilterRunIDEventID&&) = delete;
  FilterRunIDEventID& operator=(FilterRunIDEventID const&) = delete;
  FilterRunIDEventID& operator=(FilterRunIDEventID&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;
  std::vector<int> fRunIDs;
  std::vector<int> fEventIDs;

private:

  // Declare member data here.

};


FilterRunIDEventID::FilterRunIDEventID(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fRunIDs = p.get<std::vector<int>>("RunIDs");
  fEventIDs = p.get<std::vector<int>>("EventIDs");
}

bool FilterRunIDEventID::filter(art::Event& e)
{
  // Implementation of required member function here.
  int ThisRunID = e.run();
  int ThisEventID = e.event();
  bool Keep = false;
  for(int i=0; i<int(fRunIDs.size()); i++)
    {
      if(ThisRunID == fRunIDs[i] && ThisEventID==fEventIDs[i]) Keep=true;
    }
  return Keep;
}

DEFINE_ART_MODULE(FilterRunIDEventID)
