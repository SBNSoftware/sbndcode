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
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

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
  double GainFluctuation(int ch, unsigned int npe, CLHEP::HepRandomEngine* eng) override;

private:
  //Configuration parameters
  double fGainFluctuation;
  
  //PMTCalibrationDatabase service
  sbndDB::PMTCalibrationDatabase const* fPMTCalibrationDatabaseService;
};


opdet::PMTGaussianGainFluctuation::PMTGaussianGainFluctuation(art::ToolConfigTable<Config> const& config)
  :  fGainFluctuation { config().gainFluctuation() }
{
  fPMTCalibrationDatabaseService = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();
}

double opdet::PMTGaussianGainFluctuation::GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng){
  return CLHEP::RandGaussQ::shoot(eng, npe, npe*fGainFluctuation);
}

double opdet::PMTGaussianGainFluctuation::GainFluctuation(int ch, unsigned int npe, CLHEP::HepRandomEngine* eng){
  fPMTCalibrationDatabaseService->getSPEAmplitude(ch);
  //double ChannelGainFluctuation = fPMTCalibrationDatabaseService->getSPEAmplitudeStd(ch);

  std::vector<double> ChannelList = {6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 36, 37, 40, 41, 60, 61, 63, 64, 69, 70, 84, 88, 89, 90, 91, 93, 94, 95, 114, 116, 117, 119, 139, 142, 144, 145, 146, 147, 148, 149, 162, 163, 164, 165, 166, 167, 168, 169, 172, 173, 192, 193, 194, 195, 216, 219, 224, 225, 227, 240, 241, 242, 243, 244, 246, 247, 250, 251, 270, 271, 272, 273, 274, 275, 294, 297, 299, 300, 301, 303, 304, 305};
  std::vector<double> GainFluct = {5.193714426799699, 5.181913831483003, 5.280416474514692, 4.804029994238094, 4.9969993127103915, 5.6435365662225, 4.511098052025863, 5.162564651723171, 5.197506947691326, 4.73425046849519, 5.436316444065758, 4.6472001482533, 4.91402476279491, 4.552307883137502, 5.65773979888281, 5.107856508860017, 4.570284351105522, 5.723143065948864, 5.34170660098225, 4.796461113523314, 5.811696395155158, 4.88667416354998, 4.571048421004813, 5.547231675334098, 5.0685328710197375, 5.277629721752067, 5.2213848111771926, 4.922208616193263, 5.64881744805435, 5.5759636730740025, 4.762583854029147, 5.162482406117229, 5.264267331735896, 4.837534429019538, 5.446403008372831, 5.632829745266334, 5.180483432141108, 5.256186197457309, 5.474659782359494, 5.728937948014495, 5.064436046710759, 5.1586177249831024, 5.456721139322966, 5.247833952107098, 5.829033532820894, 4.4812108024765065, 4.816858211887511, 5.227035827304125, 5.293601913323989, 5.0043783343266295, 5.704066167007288, 5.4673569972880705, 5.45342402149174, 5.0720562954770285, 5.937589364327475, 5.619512815443453, 5.213692051782389, 4.6662997869595895, 4.693557659173738, 5.20806643605784, 5.217913027127723, 5.048422615327319, 4.839212955588604, 5.630252374498653, 5.461346072872436, 4.861031987566962, 4.953089182627699, 5.377492070165698, 5.597638847807842, 5.506704513834637, 4.7619912972968, 5.476291232678005, 5.5629656426391385, 5.1844318094418425, 5.222288381434315, 6.0074000373468275, 4.882696752068856, 5.491555496763724, 6.8789717575875935, 5.0113713941547715, 5.619061558817605, 5.011179457622508};
  double ChannelGainFluctuation;
  auto it = std::find(ChannelList.begin(), ChannelList.end(), ch);
  if (it != ChannelList.end()) {
      // Si se encuentra el canal, obtener el índice
      size_t index = std::distance(ChannelList.begin(), it);
      ChannelGainFluctuation = GainFluct[index];

  } else {
      ChannelGainFluctuation = 5.;
  }

  return CLHEP::RandGaussQ::shoot(eng, npe, npe*ChannelGainFluctuation);
}

DEFINE_ART_CLASS_TOOL(opdet::PMTGaussianGainFluctuation)
