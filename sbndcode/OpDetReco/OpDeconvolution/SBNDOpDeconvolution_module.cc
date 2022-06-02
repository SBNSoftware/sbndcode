////////////////////////////////////////////////////////////////////////
// Class:       SBNDOpDeconvolution
// Plugin Type: producer (art v3_06_03)
// File:        SBNDOpDeconvolution_module.cc
//
// Generated at Tue Jul 13 06:29:02 2021 by Francisco Nicolas-Arnaldos using cetskelgen
// from cetlib version v3_11_01.
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
#include "art/Utilities/make_tool.h"

#include <memory>

#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbndcode/OpDetReco/OpDeconvolution/Alg/OpDeconvolutionAlg.hh"

namespace opdet {
  class SBNDOpDeconvolution;
}


class opdet::SBNDOpDeconvolution : public art::EDProducer {
public:
  explicit SBNDOpDeconvolution(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDOpDeconvolution(SBNDOpDeconvolution const&) = delete;
  SBNDOpDeconvolution(SBNDOpDeconvolution&&) = delete;
  SBNDOpDeconvolution& operator=(SBNDOpDeconvolution const&) = delete;
  SBNDOpDeconvolution& operator=(SBNDOpDeconvolution&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:

  // Declare member data here.
  std::string fInputLabel;
  std::vector<std::string> fPDTypes;
  std::vector<std::string> fElectronics;
  //OpDecoAlg tool
  std::unique_ptr<opdet::OpDeconvolutionAlg> fOpDecoAlgPtr;
  //PDS map
  opdet::sbndPDMapAlg pdsmap;
};


opdet::SBNDOpDeconvolution::SBNDOpDeconvolution(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fInputLabel = p.get< std::string >("InputLabel");
  fPDTypes = p.get< std::vector<std::string> >("PDTypes");
  fElectronics = p.get< std::vector<std::string> >("Electronics");
  fOpDecoAlgPtr = art::make_tool<opdet::OpDeconvolutionAlg>( p.get< fhicl::ParameterSet >("OpDecoAlg") );

  produces< std::vector< raw::OpDetWaveform > >();
}

void opdet::SBNDOpDeconvolution::produce(art::Event& e)
{
  //Load the waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);
  if (!wfHandle.isValid()) {
   mf::LogError("SBNDOpDeconvolution")<<"Input waveforms with input label "<<fInputLabel<<" couldn't be loaded..."<<std::endl;
   throw cet::exception("SBNDOpDeconvolution") << "Input waveforms with input label " << fInputLabel << " not found\n";
  }

  std::vector< raw::OpDetWaveform > RawWfVector;
  RawWfVector.reserve(wfHandle->size());

  for(auto const& wf : *wfHandle){
    if((std::find(fPDTypes.begin(), fPDTypes.end(), pdsmap.pdType(wf.ChannelNumber()) ) != fPDTypes.end() ) &&
       (std::find(fElectronics.begin(), fElectronics.end(), pdsmap.electronicsType(wf.ChannelNumber()) ) != fElectronics.end()) ){
      RawWfVector.push_back(wf);
    }
  }

  std::unique_ptr< std::vector< raw::OpDetWaveform > > DecoWf_VecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
  auto & DecoWf_Vec(*DecoWf_VecPtr);
  DecoWf_Vec = fOpDecoAlgPtr->RunDeconvolution(RawWfVector);

  e.put( std::move(DecoWf_VecPtr) );

}

DEFINE_ART_MODULE(opdet::SBNDOpDeconvolution)
