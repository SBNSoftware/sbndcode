////////////////////////////////////////////////////////////////////////
// Specific class tool for PMTGainFluctuations
// File: PMTGaussianGainFluctuation.hh
// Base class:        PMTGainFluctuations.hh
// Algorithm based on function
// 'multiplicationStageGain(unsigned int i /* = 1 */) const'
// in icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <vector>

#include "sbndcode/OpDetSim/PMTAlg/PMTGainFluctuations.hh"


namespace opdet {
  class PMTGaussianGainFluctuation;
}


class opdet::PMTGaussianGainFluctuation : opdet::PMTGainFluctuations {
public:

  //Configuration parameters
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<double> gainFluctuation {
      Name("GainFluctuation"),
      Comment("PMT gain fluctuation")
    };

  };

  explicit PMTGaussianGainFluctuation(art::ToolConfigTable<Config> const& config);

  //Returns fluctuated factor for SPR
  double GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng) override;

private:
  //Configuration parameters
  double fGainFluctuation;
};


opdet::PMTGaussianGainFluctuation::PMTGaussianGainFluctuation(art::ToolConfigTable<Config> const& config)
  :  fGainFluctuation { config().gainFluctuation() }
{
}


double opdet::PMTGaussianGainFluctuation::GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng){
  return CLHEP::RandGaussQ::shoot(eng, npe, npe*fGainFluctuation);
}


DEFINE_ART_CLASS_TOOL(opdet::PMTGaussianGainFluctuation)
