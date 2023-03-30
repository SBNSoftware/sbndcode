////////////////////////////////////////////////////////////////////////
// Specific class tool for PMTGainFluctuations
// File: PMTGainFluctuations1Dynode_tool.hh
// Base class:        PMTGainFluctuations.hh
// Algorithm based on function
// 'multiplicationStageGain(unsigned int i /* = 1 */) const'
// in icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
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
      Comment("Pre bins to conisder")
    };

    fhicl::Sequence<unsigned int> nonLinearRange {
      Name("NonLinearRange"),
      Comment("Non linear range. Assume linear response/completely saturated response out of this range")
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
  std::cout<<" Non Linearity Parameteres: \n";
  std::cout<<fAttenuationForm<<"  PreTime"<<fAttenuationPreTime<<" PESat: "<<fNonLinearRange[1]<<"\n";
  for(size_t k=0; k<fAttenuationFormParams.size(); k++){
    fNonLinearTF1->SetParameter(k, fAttenuationFormParams[k]);
    std::cout<<" Par "<<k<<" "<<fAttenuationFormParams[k]<<std::endl;
  }

  // Initialize attenuation vector
  fPEAttenuation_V.resize(fNonLinearRange[1], 1);
  for(size_t pe=fNonLinearRange[0]; pe<fNonLinearRange[1]; pe++){
    fPEAttenuation_V[pe] = fNonLinearTF1->Eval(pe)/pe;
  }
  fPESaturationValue = fNonLinearTF1->Eval(fNonLinearRange[1]);

  std::cout<<" fAttFunc:  ";
  for(size_t pe=0; pe<fNonLinearRange[1]; pe++){
    std::cout<<pe<<":"<<pe*fPEAttenuation_V[pe]<<":"<<fPEAttenuation_V[pe]<<"  ";
  }
  std::cout<<" PESatValue: "<<fPESaturationValue<<std::endl;

}

double opdet::PMTNonLinearityTF1::NObservedPE(size_t bin, std::vector<unsigned int> & pe_vector){
  
  /*std::vector<unsigned int>::iterator it1 = pe_vector.begin();
  std::vector<unsigned int>::iterator it2 = std::next(pe_vector.begin(), -1);
  std::cout<<" ITER: "<<std::distance(pe_vector.begin(), it1)<<" "<<std::distance(pe_vector.begin(), it2)<<std::endl;
  */

  // get first bin
  size_t start_bin = bin-fAttenuationPreTime;
  
  if(fAttenuationPreTime<0) start_bin=0;
  unsigned int npe_acc = std::accumulate(pe_vector.begin()+start_bin,pe_vector.begin()+bin+1, 0);
  
  /*std::cout<<" In bin: "<<bin<<std::endl;
  for(size_t k=start_bin; k<=bin; k++){
    std::cout<<"  "<<k<<" "<<pe_vector[k];
  }
  std::cout<< "    Acc"<<npe_acc<<"  Rescales: "<<pe_vector[bin]*fPEAttenuation_V[npe_acc]<<std::endl;*/

  //std::cout<< "   Bin="<<bin<<" Acc"<<npe_acc<<"  Rescales: "<<pe_vector[bin]*fPEAttenuation_V[npe_acc]<<std::endl;
  
  if(pe_vector[bin]>100) std::cout<<"  High En: "<<pe_vector[bin]<<" Acc="<<npe_acc<<" "<<pe_vector[bin]*fPEAttenuation_V[npe_acc]<<std::endl;

  if(npe_acc<fNonLinearRange[1]) return pe_vector[bin]*fPEAttenuation_V[npe_acc];
  else return fPESaturationValue;
}


DEFINE_ART_CLASS_TOOL(opdet::PMTNonLinearityTF1)