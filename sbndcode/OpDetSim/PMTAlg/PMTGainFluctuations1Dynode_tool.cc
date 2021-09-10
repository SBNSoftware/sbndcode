////////////////////////////////////////////////////////////////////////
// Specific class tool for PMTGainFluctuations
// File: PMTGainFluctuations1Dynode_tool.hh
// Base class:        PMTGainFluctuations.hh
// Algorithm based on function
// 'multiplicationStageGain(unsigned int i /* = 1 */) const'
// in icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include <vector>

#include "sbndcode/OpDetSim/PMTAlg/PMTGainFluctuations.hh"


namespace opdet {
  class PMTGainFluctuations1Dynode;
}


class opdet::PMTGainFluctuations1Dynode : opdet::PMTGainFluctuations {
public:
  explicit PMTGainFluctuations1Dynode(fhicl::ParameterSet const& p);

  ~PMTGainFluctuations1Dynode() {}

  //Returns fluctuated factor for SPR
  double GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng) override;

private:
  //Configuration parameters
  double fDynodeK;
  unsigned int fGain;
  std::vector<double> fVoltageDistribution;

  double fRefGain;

  unsigned int DynodeGain(unsigned int dynstage);
};


opdet::PMTGainFluctuations1Dynode::PMTGainFluctuations1Dynode(fhicl::ParameterSet const& p)
{
  //read fhicl paramters
  fDynodeK = p.get< double >("DynodeK");
  fGain = p.get< unsigned int >("Gain");
  fVoltageDistribution  = p.get< std::vector<double> >("VoltageDistribution");

  fRefGain = DynodeGain(1);
}


unsigned int opdet::PMTGainFluctuations1Dynode::DynodeGain(unsigned int dynstage){
  double prodRho = 1.0;
  for(double rho: fVoltageDistribution) prodRho *= rho;
  double const aVk = std::pow(fGain / std::pow(prodRho, fDynodeK), 1.0/static_cast<double>(fVoltageDistribution.size()));
  return aVk * std::pow(fVoltageDistribution.at(dynstage - 1), fDynodeK);
}


double opdet::PMTGainFluctuations1Dynode::GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng){
  return CLHEP::RandPoissonQ::shoot(eng, npe*fRefGain)/fRefGain;
}


DEFINE_ART_CLASS_TOOL(opdet::PMTGainFluctuations1Dynode)
