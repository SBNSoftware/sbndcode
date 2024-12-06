#include "CosmicIdUtils.h"

namespace sbnd{

// =============================== UTILITY FUNCTIONS ==============================

  // Create fake PDS optical flashes from true particle energy deposits
  std::pair<std::vector<double>, std::vector<double>> CosmicIdUtils::FakeTpcFlashes(std::vector<simb::MCParticle> particles){
    //
    TPCGeoAlg fTpcGeo;

    // Create fake flashes in each tpc
    std::vector<double> fakeTpc0Flashes;
    std::vector<double> fakeTpc1Flashes;

    // FIXME probably shouldn't be here
    //double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    //double driftTimeMuS = fTpcGeo.MaxX()/fDetectorProperties->DriftVelocity(); // [us]

    // Loop over all true particles
    for (auto const &particle: particles){

      // Get particle info
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3;

      //Check if time is in reconstructible window
      //if(time < -driftTimeMuS || time > readoutWindowMuS) continue; 
      //Check if particle is visible, electron, muon, proton, pion, kaon, photon
      if(!(pdg==13||pdg==11||pdg==22||pdg==2212||pdg==211||pdg==321||pdg==111)) continue;

      //Loop over the trajectory
      int npts = particle.NumberTrajectoryPoints();
      double TPC0Energy = 0;
      double TPC1Energy = 0;
      for(int i = 1; i < npts; i++){
        geo::Point_t pt;
        pt.SetX(particle.Vx(i)); pt.SetY(particle.Vy(i)); pt.SetZ(particle.Vz(i));
        if(!fTpcGeo.InFiducial(pt, 0, 0)) continue;
        // Add up the energy deposited in each tpc
        if(pt.X() <= 0) TPC0Energy += particle.E(i-1) - particle.E(i);
        else TPC1Energy += particle.E(i-1) - particle.E(i);
      }
      // If the total energy deposited is > 10 MeV then create fake flash
      if(TPC0Energy > 0.01) fakeTpc0Flashes.push_back(time);
      else if(TPC1Energy > 0.01) fakeTpc1Flashes.push_back(time);
    }

    std::sort(fakeTpc0Flashes.begin(), fakeTpc0Flashes.end());
    double previousTime = -99999;
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < fakeTpc0Flashes.size(); i++){
      double time = fakeTpc0Flashes[i];
      // Combine flashes within 0.01 us
      if(std::abs(time-previousTime) < 0.01){
        fakeTpc0Flashes.erase(fakeTpc0Flashes.begin()+i);
      }
      else previousTime = time;
    }

    std::sort(fakeTpc1Flashes.begin(), fakeTpc1Flashes.end());
    previousTime = -99999;
    // Loop over flashes in tpc 1
    for(size_t i = 0; i < fakeTpc1Flashes.size(); i++){
      double time = fakeTpc1Flashes[i];
      // Combine flashes within 0.01 us
      if(std::abs(time-previousTime) < 0.01){
        fakeTpc1Flashes.erase(fakeTpc1Flashes.begin()+i);
      }
      else previousTime = time;
    }

    return std::make_pair(fakeTpc0Flashes, fakeTpc1Flashes);
  }

  // Determine if there is a PDS flash in time with the neutrino beam
  bool CosmicIdUtils::BeamFlash(std::vector<double> flashes, double beamTimeMin, double beamTimeMax){
    //
    bool beamFlash = false;
    std::sort(flashes.begin(), flashes.end());
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < flashes.size(); i++){
      double time = flashes[i];
      if(time > beamTimeMin && time < beamTimeMax) beamFlash = true;
    }

    return beamFlash;
  }

}
