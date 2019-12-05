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
#include "sbndcode/OpDetSim/sbndPDMapAlg.h"

// LArSoft
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// c++
#include <vector>
#include <utility>

namespace sbnd{
namespace CosmicIdUtils{

  std::vector<double> FlashTimes(std::vector<double> optimes);

  // Create fake PDS optical flashes from true particle energy deposits
  std::pair<std::vector<double>, std::vector<double>> FakeTpcFlashes(std::vector<simb::MCParticle> particles);

  // Create real PDS optical flashes
  std::pair<std::vector<double>, std::vector<double>> OpFlashes(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle);

  // Determine if there is a PDS flash in time with the neutrino beam
  bool BeamFlash(std::vector<double> flashes, double beamTimeMin, double beamTimeMax);

  // Determine if there is a PDS flash in time with the neutrino beam
  std::pair<bool, bool> BeamFlash(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle, double beamTimeMin, double beamTimeMax);
  
}
}

#endif
