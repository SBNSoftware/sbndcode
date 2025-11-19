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
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandFlat.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include <memory>

#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "TF1.h"
#include "TH1F.h"

namespace opdet {
  class PMTPulseOscillation;
}


class opdet::PMTPulseOscillation : public art::EDProducer {
  public:
    explicit PMTPulseOscillation(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PMTPulseOscillation(PMTPulseOscillation const&) = delete;
    PMTPulseOscillation(PMTPulseOscillation&&) = delete;
    PMTPulseOscillation& operator=(PMTPulseOscillation const&) = delete;
    PMTPulseOscillation& operator=(PMTPulseOscillation&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;
    void CreateOscillatedWaveform(const raw::OpDetWaveform& oldWaveform , raw::OpDetWaveform& newWaveform );


  private:

    // Declare member data here.
    std::string fInputLabel;
    std::vector<std::string> fPDTypes;
    std::vector<std::string> fElectronics;
    std::string fOscillationFunction;
    double fOscillationFrequency;
    double fOscillationDampingConstant;
    double fSamplingPeriod;
    double fPulseOscillationThreshold;
    int fOscillationOffset;
    //PDS map
    opdet::sbndPDMapAlg pdsmap;
    std::vector<double> fOscillationFunctionParams;
    TF1 *fOscillationTF1;
    CLHEP::HepRandomEngine& fEngine;
};


opdet::PMTPulseOscillation::PMTPulseOscillation(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fEngine(art::ServiceHandle<rndm::NuRandomService>{}->registerAndSeedEngine(
              createEngine(0, "HepJamesRandom", "pulseosc"), "HepJamesRandom", "pulseosc", p, "Seed"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fInputLabel = p.get< std::string >("InputLabel");
  fPDTypes = p.get< std::vector<std::string> >("PDTypes");
  fElectronics = p.get< std::vector<std::string> >("Electronics");
  fOscillationFunction = p.get<std::string>("OscillationFunction");
  fSamplingPeriod = p.get<double>("SamplingPeriod");
  fOscillationFrequency = p.get<double>("OscillationFrequency");
  fOscillationDampingConstant = p.get<double>("OscillationDampingConstant");
  fOscillationFunctionParams = p.get<std::vector<double>>("OscillationFunctionParams");
  fPulseOscillationThreshold = p.get<double>("PulseOscillationThreshold");
  fOscillationOffset = p.get<int>("OscillationOffset");
  fOscillationTF1 = new TF1("OscillationTemplate", fOscillationFunction.c_str());
  for(size_t i=0; i<fOscillationFunctionParams.size(); i++){
    fOscillationTF1->SetParameter(i, fOscillationFunctionParams[i]);
  }

  produces< std::vector< raw::OpDetWaveform > >();
}

void opdet::PMTPulseOscillation::produce(art::Event& e)
{
  //Load the waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);
  if (!wfHandle.isValid()) {  
   std::cout << " not finding the handles " << std::endl;
   mf::LogError("PMTPulseOscillation")<<"Input waveforms with input label "<<fInputLabel<<" couldn't be loaded..."<<std::endl;
   throw cet::exception("PMTPulseOscillation") << "Input waveforms with input label " << fInputLabel << " not found\n";
  }
  std::unique_ptr< std::vector< raw::OpDetWaveform > > OscillatedWf_VecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
  for(auto const& wf : *wfHandle){
      if(std::find(fPDTypes.begin(), fPDTypes.end(), pdsmap.pdType(wf.ChannelNumber()) ) != fPDTypes.end() ){
        raw::OpDetWaveform new_wf = wf;
        CreateOscillatedWaveform(wf, new_wf);
        OscillatedWf_VecPtr->push_back(new_wf);
      }
      else OscillatedWf_VecPtr->push_back(wf);
  }

  e.put( std::move(OscillatedWf_VecPtr) );
}


void opdet::PMTPulseOscillation::CreateOscillatedWaveform(const raw::OpDetWaveform& oldWaveform, raw::OpDetWaveform& newWaveform)
{
    //Find the idx of the maximum (or minimum) point of the waveform
    int starting_ttick = std::distance(oldWaveform.begin(), std::min_element(oldWaveform.begin(), oldWaveform.end()));
    double baseline = accumulate( oldWaveform.begin(), oldWaveform.end(), 0.0) / oldWaveform.size();
    double pulse_max = abs(oldWaveform[starting_ttick]-baseline);
    if(pulse_max<fPulseOscillationThreshold) return;
    double oscillation_amplitude = fOscillationTF1->Eval(pulse_max);
    for(size_t i=starting_ttick; i<oldWaveform.size()-fOscillationOffset; i++){
        double delta_time = (static_cast<double>(i)-static_cast<double>(starting_ttick))*fSamplingPeriod;
        double envelope = std::exp(-delta_time / fOscillationDampingConstant);
        double phase    = 2.0 * M_PI * fOscillationFrequency * delta_time;
        double new_value = static_cast<double>(newWaveform.Waveform()[i+fOscillationOffset]) + oscillation_amplitude * envelope * std::cos(phase+M_PI);
        // Add some dither to avoid quantization effects
        double dither = CLHEP::RandFlat::shoot(&fEngine,1.0) - 0.5;
        newWaveform.Waveform()[i+fOscillationOffset] = std::round(new_value + dither);
        if(envelope < 0.05) return;
    }
}

DEFINE_ART_MODULE(opdet::PMTPulseOscillation)