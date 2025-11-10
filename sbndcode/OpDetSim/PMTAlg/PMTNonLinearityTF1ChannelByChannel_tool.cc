////////////////////////////////////////////////////////////////////////
// Specific class tool for PMTNonLinearity
// File: PMTNonLinearityTF1_tool.hh
// Base class:        PMTNonLinearity.hh
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <vector>
#include <numeric>
#include <string>
#include "TF1.h"

#include "sbndcode/OpDetSim/PMTAlg/PMTNonLinearity.hh"


namespace opdet {
  class PMTNonLinearityTF1ChannelByChannel;
}


class opdet::PMTNonLinearityTF1ChannelByChannel : opdet::PMTNonLinearity {
public:

  //Configuration parameters
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<std::string> attenuationForm {
      Name("AttenuationForm"),
      Comment("Non linearity functional form")
    };
    
    fhicl::Sequence<fhicl::Sequence<double>> attenuationFormParams {
      Name("AttenuationFormParams"),
      Comment("Parameters for non linearity functional form")
    };

    fhicl::Atom<unsigned int> attenuationPreTime {
      Name("AttenuationPreTime"),
      Comment("For the non-linear attenuation, consider photons arriving given by this time window. In nanoseconds.")
    };

    fhicl::Sequence<unsigned int> nonLinearRange {
      Name("NonLinearRange"),
      Comment("Non linear range. Assume linear response/completely saturated response out of this range. In #PE.")
    };

  };

  explicit PMTNonLinearityTF1ChannelByChannel(art::ToolConfigTable<Config> const& config);

  //Returns rescaled #pe after non linearity
  double NObservedPE(int opch, size_t bin, std::vector<unsigned int> & pe_vector) override;

private:
  //Configuration parameters
  std::string fAttenuationForm;
  std::vector<std::vector<double>> fAttenuationFormParams;
  unsigned int fAttenuationPreTime;
  std::vector<unsigned int> fNonLinearRange;

  geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  //TF1 for non linearity function
  std::map<int, TF1*> fNonLinearTF1Map;

  // Vector to store non linearity attenuation values
  std::vector<std::vector<double>> fPEAttenuation_V;
  std::vector<double> fPEAttenuation;
  std::vector<int> fPESaturationValue_V;
  int fPESaturationValue;

};


opdet::PMTNonLinearityTF1ChannelByChannel::PMTNonLinearityTF1ChannelByChannel(art::ToolConfigTable<Config> const& config)
  : fAttenuationForm { config().attenuationForm() }
  , fAttenuationFormParams { config().attenuationFormParams() }
  , fAttenuationPreTime { config().attenuationPreTime() }
  , fNonLinearRange { config().nonLinearRange() }
{
  //Initialize TF1 for each channel with channel-dependent parameters
  for(size_t opch=0; opch<fWireReadout.NOpChannels(); opch++){
    fNonLinearTF1Map[opch] = new TF1(("NonLinearTF1_"+std::to_string(opch)).c_str(), fAttenuationForm.c_str());
    for(size_t k=0; k<fAttenuationFormParams[opch].size(); k++){
      fNonLinearTF1Map[opch]->SetParameter(k, fAttenuationFormParams[opch][k]);
    }
  }

  // Initialize attenuation vector
  fPEAttenuation_V.resize(fWireReadout.NOpChannels());
  for(size_t opch=0; opch<fPEAttenuation_V.size(); opch++){
    fPEAttenuation_V[opch].resize(fNonLinearRange[1], 1);
    for(size_t pe=fNonLinearRange[0]; pe<fNonLinearRange[1]; pe++){
      fPEAttenuation_V[opch][pe] = fNonLinearTF1Map[opch]->Eval(pe)/pe;
    }
  }

  fPESaturationValue_V.resize(fWireReadout.NOpChannels(), 1);
  for(size_t opch=0; opch<fWireReadout.NOpChannels(); opch++){
    fPESaturationValue_V[opch] = std::round(fNonLinearTF1Map[opch]->Eval(fNonLinearRange[1]));
  }
}

double opdet::PMTNonLinearityTF1ChannelByChannel::NObservedPE(int opch, size_t bin, std::vector<unsigned int> & pe_vector){
  size_t start_bin = bin-fAttenuationPreTime;
  if(fAttenuationPreTime<0) start_bin=0;
  unsigned int npe_acc = std::accumulate(pe_vector.begin()+start_bin,pe_vector.begin()+bin+1, 0);

  if(npe_acc<fNonLinearRange[1]) return pe_vector[bin]*fPEAttenuation_V[opch][npe_acc];
  else return fPESaturationValue_V[opch];
}

DEFINE_ART_CLASS_TOOL(opdet::PMTNonLinearityTF1ChannelByChannel)
