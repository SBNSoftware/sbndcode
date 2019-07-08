#include "CosmicRemovalUtils.h"

namespace sbnd{

// =============================== UTILITY FUNCTIONS ==============================

  bool CosmicRemovalUtils::InFiducial(geo::Point_t point, double fiducial){

    return InFiducial(point, fiducial, fiducial);

  }

  bool CosmicRemovalUtils::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){

    return InFiducial(point, fiducial, fiducial, fiducial, fiducial, fiducialTop, fiducial);

  }

  bool CosmicRemovalUtils::InFiducial(geo::Point_t point, double minXCut, double minYCut, double minZCut, double maxXCut, double maxYCut, double maxZCut){
    
    TPCGeoAlg geo;
    double xmin = geo.MinX() + minXCut;
    double xmax = geo.MaxX() - maxXCut;
    double ymin = geo.MinY() + minYCut;
    double ymax = geo.MaxY() - maxYCut;
    double zmin = geo.MinZ() + minZCut;
    double zmax = geo.MaxZ() - maxZCut;

    double x = point.X();
    double y = point.Y();
    double z = point.Z();
    if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) return true;

    return false;
  }

  int CosmicRemovalUtils::DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits){
    //
    int tpc = hits[0]->WireID().TPC;
    for(size_t i = 0; i < hits.size(); i++){
      if((int)hits[i]->WireID().TPC != tpc) return -1;
    }
    return tpc;
  }

  std::pair<std::vector<double>, std::vector<double>> CosmicRemovalUtils::FakeTpcFlashes(std::vector<simb::MCParticle> particles){
    //
    TPCGeoAlg geo;
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    detinfo::DetectorClocks const* fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

    // Create fake flashes in each tpc
    std::vector<double> fakeTpc0Flashes;
    std::vector<double> fakeTpc1Flashes;

    double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    double driftTimeMuS = geo.MaxX()/fDetectorProperties->DriftVelocity(); // [us]

    // Loop over all true particles
    for (auto const particle: particles){

      // Get particle info
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3;

      //Check if time is in reconstructible window
      if(time < -driftTimeMuS || time > readoutWindowMuS) continue; 
      //Check if particle is visible, electron, muon, proton, pion, kaon, photon
      if(!(pdg==13||pdg==11||pdg==22||pdg==2212||pdg==211||pdg==321||pdg==111)) continue;

      //Loop over the trajectory
      int npts = particle.NumberTrajectoryPoints();
      double TPC0Energy = 0;
      double TPC1Energy = 0;
      for(int i = 1; i < npts; i++){
        geo::Point_t pt;
        pt.SetX(particle.Vx(i)); pt.SetY(particle.Vy(i)); pt.SetZ(particle.Vz(i));
        if(!CosmicRemovalUtils::InFiducial(pt, 0, 0)) continue;
        // Add up the energy deposited in each tpc
        if(pt.X() < 0) TPC0Energy += particle.E(i-1) - particle.E(i);
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
      // Combine flashes within 0.1 us
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc0Flashes.erase(fakeTpc0Flashes.begin()+i);
      }
      else previousTime = time;
    }

    std::sort(fakeTpc1Flashes.begin(), fakeTpc1Flashes.end());
    previousTime = -99999;
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < fakeTpc1Flashes.size(); i++){
      double time = fakeTpc1Flashes[i];
      // Combine flashes within 0.1 us
      if(std::abs(time-previousTime) < 0.1){
        fakeTpc1Flashes.erase(fakeTpc1Flashes.begin()+i);
      }
      else previousTime = time;
    }

    return std::make_pair(fakeTpc0Flashes, fakeTpc1Flashes);
  }

  bool CosmicRemovalUtils::BeamFlash(std::vector<double> flashes, double beamTimeLimit){
    //
    bool beamFlash = false;
    std::sort(flashes.begin(), flashes.end());
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < flashes.size(); i++){
      double time = flashes[i];
      if(time > 0 && time < beamTimeLimit) beamFlash = true;
    }

    return beamFlash;
  }

}
