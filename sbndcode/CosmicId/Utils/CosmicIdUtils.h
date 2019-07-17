#ifndef COSMICIDUTILS_H_SEEN
#define COSMICIDUTILS_H_SEEN


///////////////////////////////////////////////
// CosmicIdUtils.h
//
// Reco utilities for doing cosmic removal in ana modules
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// c++
#include <vector>
#include <utility>

namespace sbnd{
namespace CosmicIdUtils{

  std::pair<std::vector<double>, std::vector<double>> FakeTpcFlashes(std::vector<simb::MCParticle> particles);

  bool BeamFlash(std::vector<double> flashes, double beamTimeMin, double beamTimeMax);
  
}
}

#endif
