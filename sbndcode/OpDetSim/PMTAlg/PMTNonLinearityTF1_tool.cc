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

#include <vector>
#include <numeric>
#include <string>
#include "TF1.h"

#include "sbndcode/OpDetSim/PMTAlg/PMTNonLinearity.hh"


namespace opdet {
  class PMTNonLinearityTF1;
}


class opdet::PMTNonLinearityTF1 : opdet::PMTNonLinearity {
public:

  //Configuration parameters
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<std::string> attenuationForm {
      Name("AttenuationForm"),
      Comment("Non linearity functional form")
    };

    fhicl::Sequence<double> attenuationFormParams {
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

  explicit PMTNonLinearityTF1(art::ToolConfigTable<Config> const& config);

  //Returns rescaled #pe after non linearity
  double NObservedPE(size_t bin, std::vector<unsigned int> & pe_vector) override;

private:
  //Configuration parameters
  std::string fAttenuationForm;
  std::vector<double> fAttenuationFormParams;
  unsigned int fAttenuationPreTime;
  std::vector<unsigned int> fNonLinearRange;
  unsigned int fPEStartSat;
  unsigned int fPEMaxSat;

  //TF1 for non linearity function
  TF1 *fNonLinearTF1;

  // Vector to store non linearity attenuation values
  std::vector<double> fPEAttenuation_V;
  int fPESaturationValue;
};


opdet::PMTNonLinearityTF1::PMTNonLinearityTF1(art::ToolConfigTable<Config> const& config)
  : fAttenuationForm { config().attenuationForm() }
  , fAttenuationFormParams { config().attenuationFormParams() }
  , fAttenuationPreTime { config().attenuationPreTime() }
  , fNonLinearRange { config().nonLinearRange() }
{
  fNonLinearTF1 = new TF1("NonLinearTF1", fAttenuationForm.c_str());
  for(size_t k=0; k<fAttenuationFormParams.size(); k++){
    fNonLinearTF1->SetParameter(k, fAttenuationFormParams[k]);
  }

  // Initialize attenuation vector
  fPEAttenuation_V.resize(fNonLinearRange[1], 1);
  for(size_t pe=fNonLinearRange[0]; pe<fNonLinearRange[1]; pe++){
    fPEAttenuation_V[pe] = fNonLinearTF1->Eval(pe)/pe;
  }
  fPESaturationValue = fNonLinearTF1->Eval(fNonLinearRange[1]);

}

double opdet::PMTNonLinearityTF1::NObservedPE(size_t bin, std::vector<unsigned int> & pe_vector){

  // get first bin
  size_t start_bin = bin-fAttenuationPreTime;
  
  if(fAttenuationPreTime<0) start_bin=0;
  unsigned int npe_acc = std::accumulate(pe_vector.begin()+start_bin,pe_vector.begin()+bin+1, 0);

  if(npe_acc<fNonLinearRange[1]) return pe_vector[bin]*fPEAttenuation_V[npe_acc];
  else return fPESaturationValue;
}


DEFINE_ART_CLASS_TOOL(opdet::PMTNonLinearityTF1)
