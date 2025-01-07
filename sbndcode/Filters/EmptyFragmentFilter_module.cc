////////////////////////////////////////////////////////////////////////
// Class:       EmptyFragmentFilter
// Plugin Type: filter (Unknown Unknown)
// File:        EmptyFragmentFilter_module.cc
//
// Generated at Tue Jan  7 09:14:30 2025 by Henry Lay using cetskelgen
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

#include "artdaq-core/Data/ContainerFragment.hh"

namespace sbnd {
  namespace raw {
    class EmptyFragmentFilter;
  }
}


class sbnd::raw::EmptyFragmentFilter : public art::EDFilter {
public:
  explicit EmptyFragmentFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EmptyFragmentFilter(EmptyFragmentFilter const&) = delete;
  EmptyFragmentFilter(EmptyFragmentFilter&&) = delete;
  EmptyFragmentFilter& operator=(EmptyFragmentFilter const&) = delete;
  EmptyFragmentFilter& operator=(EmptyFragmentFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  std::string fEmptyFragmentModuleLabel;
  std::vector<std::string> fEmptyFragmentInstanceLabels;
};


sbnd::raw::EmptyFragmentFilter::EmptyFragmentFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
  , fEmptyFragmentModuleLabel(p.get<std::string>("EmptyFragmentModuleLabel"))
  , fEmptyFragmentInstanceLabels(p.get<std::vector<std::string>>("EmptyFragmentInstanceLabels"))
{
}

bool sbnd::raw::EmptyFragmentFilter::filter(art::Event& e)
{
  for(const std::string &emptyFragmentInstanceLabel : fEmptyFragmentInstanceLabels)
    {
      art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
      e.getByLabel(fEmptyFragmentModuleLabel, emptyFragmentInstanceLabel, fragmentHandle);

      if(fragmentHandle.isValid())
        return false;
    }

  return true;
}

DEFINE_ART_MODULE(sbnd::raw::EmptyFragmentFilter)
