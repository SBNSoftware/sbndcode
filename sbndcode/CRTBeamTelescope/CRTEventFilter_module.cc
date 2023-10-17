////////////////////////////////////////////////////////////////////////
// Class:       CRTEventFilter
// Plugin Type: filter (Unknown Unknown)
// File:        CRTEventFilter_module.cc
//
// Generated at Mon Oct  2 08:44:27 2023 by Jiaoyang Li using cetskelgen
// from  version .
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
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include <memory>

class CRTEventFilter;


class CRTEventFilter : public art::EDFilter {
public:
  explicit CRTEventFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTEventFilter(CRTEventFilter const&) = delete;
  CRTEventFilter(CRTEventFilter&&) = delete;
  CRTEventFilter& operator=(CRTEventFilter const&) = delete;
  CRTEventFilter& operator=(CRTEventFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fAuxDetIDELabel;

};


CRTEventFilter::CRTEventFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fAuxDetIDELabel = p.get<art::InputTag>("AuxDetIDELabel");
}

bool CRTEventFilter::filter(art::Event& e)
{
  auto AuxIDEHandle = e.getValidHandle<std::vector<sim::AuxDetIDE>>(fAuxDetIDELabel);
  bool isEventPassedBeamTele = false;

  if (AuxIDEHandle->size()>0){
    isEventPassedBeamTele = true;
  }

  return isEventPassedBeamTele;
}

DEFINE_ART_MODULE(CRTEventFilter)
