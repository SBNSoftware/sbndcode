////////////////////////////////////////////////////////////////////////
// Class:       ImportCRTSharpData
// Plugin Type: producer (Unknown Unknown)
// File:        ImportCRTSharpData_module.cc
//
// Generated at Tue Oct 18 19:19:21 2022 by Henry Lay using cetskelgen
// from  version .
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

#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"

class ImportCRTSharpData;


class ImportCRTSharpData : public art::EDProducer {
public:
  explicit ImportCRTSharpData(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ImportCRTSharpData(ImportCRTSharpData const&) = delete;
  ImportCRTSharpData(ImportCRTSharpData&&) = delete;
  ImportCRTSharpData& operator=(ImportCRTSharpData const&) = delete;
  ImportCRTSharpData& operator=(ImportCRTSharpData&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  void FragToFEB(const artdaq::Fragment &frag);

private:

  // Declare member data here.

};


ImportCRTSharpData::ImportCRTSharpData(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ImportCRTSharpData::produce(art::Event& e)
{
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  // Loop fragment handles
  for(auto handle : fragmentHandles)
    {
      if(!handle.isValid() || handle->size() == 0)
        continue;

      // Container or standard?
      if(handle->front().type() == artdaq::Fragment::ContainerFragmentType)
	{
	  for(auto cont : *handle)
	    {
	      artdaq::ContainerFragment contf(cont);
	      if(contf.fragment_type() == sbndaq::detail::FragmentType::BERNCRTV2)
		{
		  for(unsigned i = 0; i < contf.block_count(); ++i)
		    FragToFEB(*contf[i].get());
		}
	    }
	}
    }
}

void ImportCRTSharpData::FragToFEB(const artdaq::Fragment &frag)
{
  sbndaq::BernCRTFragmentV2 bern_frag(frag);

  std::cout << bern_frag.metadata()->hits_in_fragment() << std::endl;
}

DEFINE_ART_MODULE(ImportCRTSharpData)
