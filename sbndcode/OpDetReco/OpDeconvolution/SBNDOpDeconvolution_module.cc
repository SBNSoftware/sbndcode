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
  //OpDecoAlg
  opdet::OpDeconvolutionAlg *fOpDecoAlg;
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
  auto const opdeco_pset = p.get< fhicl::ParameterSet >("OpDecoAlg");
  fOpDecoAlg = new opdet::OpDeconvolutionAlg(opdeco_pset);

  produces< std::vector< raw::OpDetWaveform > >();
}

void opdet::SBNDOpDeconvolution::produce(art::Event& e)
{
  //Load the waveforms 
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);
  if (!wfHandle.isValid()) {
   std::cout<<"Non valid waveform handle\n";
   return;
  }

  std::vector< raw::OpDetWaveform > DecoWfVector;
  DecoWfVector.reserve(wfHandle->size());

  for(auto const& wf : *wfHandle){
    if(std::find(fPDTypes.begin(), fPDTypes.end(), pdsmap.pdType(wf.ChannelNumber()) ) != fPDTypes.end() ){
      DecoWfVector.push_back(wf);
    }
  }

  DecoWfVector=fOpDecoAlg->RunDeconvolution(DecoWfVector);

  std::unique_ptr< std::vector< raw::OpDetWaveform > > DecoWf_VecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
  for (auto & wf : DecoWfVector) {
    (*DecoWf_VecPtr).emplace_back( wf.TimeStamp(), wf.ChannelNumber(),
      std::vector<short unsigned int> (wf.Waveform().begin(), wf.Waveform().end()) );
  }

  e.put( std::move(DecoWf_VecPtr) );

}

DEFINE_ART_MODULE(opdet::SBNDOpDeconvolution)
