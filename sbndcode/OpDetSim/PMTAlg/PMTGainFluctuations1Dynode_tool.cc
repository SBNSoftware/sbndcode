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
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <vector>

#include "sbndcode/OpDetSim/PMTAlg/PMTGainFluctuations.hh"


namespace opdet {
  class PMTGainFluctuations1Dynode;
}


class opdet::PMTGainFluctuations1Dynode : opdet::PMTGainFluctuations {
public:

  //Configuration parameters
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<double> dynodeK {
      Name("DynodeK"),
      Comment("Gain at each stage is proportional to pow(Voltage,DynodeK)")
    };

    fhicl::Atom<double> gain {
      Name("Gain"),
      Comment("PMT total gain")
    };

    fhicl::Sequence<double> voltageDistribution {
      Name("VoltageDistribution"),
      Comment("PMT voltage distribution ratio at each dynode")
    };

  };

  explicit PMTGainFluctuations1Dynode(art::ToolConfigTable<Config> const& config);

  //Returns fluctuated factor for SPR
  double GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng) override;

private:
  //Configuration parameters
  double fDynodeK;
  double fGain;
  std::vector<double> fVoltageDistribution;

  double fDynodeGain;

  double DynodeGain(unsigned int dynstage);
};


opdet::PMTGainFluctuations1Dynode::PMTGainFluctuations1Dynode(art::ToolConfigTable<Config> const& config)
  : fDynodeK { config().dynodeK() }
  , fGain { config().gain() }
  , fVoltageDistribution { config().voltageDistribution() }
  , fDynodeGain { DynodeGain(1) }
{
}


double opdet::PMTGainFluctuations1Dynode::DynodeGain(unsigned int dynstage){
  double prodRho = 1.0;
  for(double rho: fVoltageDistribution) prodRho *= rho;
  double const aVk = std::pow(fGain / std::pow(prodRho, fDynodeK), 1.0/static_cast<double>(fVoltageDistribution.size()));
  return aVk * std::pow(fVoltageDistribution.at(dynstage - 1), fDynodeK);
}


double opdet::PMTGainFluctuations1Dynode::GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng){
  return CLHEP::RandPoissonQ::shoot(eng, npe*fDynodeGain)/fDynodeGain;
}


DEFINE_ART_CLASS_TOOL(opdet::PMTGainFluctuations1Dynode)
