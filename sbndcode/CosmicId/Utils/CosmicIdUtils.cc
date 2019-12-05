#include "CosmicIdUtils.h"

namespace sbnd{

// =============================== UTILITY FUNCTIONS ==============================

  std::vector<double> CosmicIdUtils::FlashTimes(std::vector<double> optimes){
    std::vector<double> opflashes;

    // Sort PMT times by time
    std::sort(optimes.begin(), optimes.end());

    for(size_t i = 0; i < optimes.size(); i++){
      // Flash reconstruction
      double start_time = optimes[i];
      size_t nhits = 0;
      std::vector<double> times {start_time};
      while( (i+nhits)<(optimes.size()-1) && (optimes[i+nhits]-start_time) < 6){
        nhits++;
        times.push_back(optimes[i+nhits]);
      }
      if(nhits > 100){
        opflashes.push_back((std::accumulate(times.begin(), times.end(), 0.)/times.size())-2.5);
        i = i + nhits;
      }
    }

    return opflashes;
  }

  // Create real PDS optical flashes
  std::pair<std::vector<double>, std::vector<double>> CosmicIdUtils::OpFlashes(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

    opdet::sbndPDMapAlg fChannelMap; //map for photon detector types

    // Optical flash reconstruction for numuCC
    std::vector<double> optimes_tpc0;
    std::vector<double> optimes_tpc1;
    for(auto const& ophit : (*pdsHandle)){
      // Only look at PMTs
      std::string od = fChannelMap.pdName(ophit.OpChannel());
      if( od != "pmt" ) continue;
      // Work out what TPC detector is in odd = TPC1, even = TPC0
      if(ophit.OpChannel() % 2 == 0){ 
        optimes_tpc1.push_back(ophit.PeakTime());
      }
      else{ 
        optimes_tpc0.push_back(ophit.PeakTime());
      }
    }

    std::vector<double> opflashes_tpc0 = FlashTimes(optimes_tpc0);

    std::vector<double> opflashes_tpc1 = FlashTimes(optimes_tpc1);

    return std::make_pair(opflashes_tpc0, opflashes_tpc1);
  }
  
  // Create fake PDS optical flashes from true particle energy deposits
  std::pair<std::vector<double>, std::vector<double>> CosmicIdUtils::FakeTpcFlashes(std::vector<simb::MCParticle> particles){
    //
    TPCGeoAlg fTpcGeo;

    // Create fake flashes in each tpc
    std::vector<double> fakeTpc0Flashes;
    std::vector<double> fakeTpc1Flashes;

    // Loop over all true particles
    for (auto const particle: particles){

      // Get particle info
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3;

      //Check if time is in reconstructible window
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
  std::pair<bool, bool> CosmicIdUtils::BeamFlash(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle, double beamTimeMin, double beamTimeMax){

    opdet::sbndPDMapAlg fChannelMap; //map for photon detector types

    int n_hits_tpc0 = 0;
    int n_hits_tpc1 = 0;
    for(auto const& ophit : (*pdsHandle)){
      // Only look at PMTs
      std::string od = fChannelMap.pdName(ophit.OpChannel());
      if( od != "pmt" ) continue;
      // Work out what TPC detector is in odd = TPC1, even = TPC0
      if(ophit.OpChannel() % 2 == 0){ 
        // Beam activity
        if(ophit.PeakTime() >= beamTimeMin && ophit.PeakTime() <= beamTimeMax+2.5) n_hits_tpc1++;
      }
      else{ 
        // Beam activity
        if(ophit.PeakTime() >= beamTimeMin && ophit.PeakTime() <= beamTimeMax+2.5) n_hits_tpc0++;
      }
    }

    bool tpc0Flash = false;
    if(n_hits_tpc0 > 100) tpc0Flash = true;

    bool tpc1Flash = false;
    if(n_hits_tpc1 > 100) tpc1Flash = true;

    return std::make_pair(tpc0Flash, tpc1Flash);
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
