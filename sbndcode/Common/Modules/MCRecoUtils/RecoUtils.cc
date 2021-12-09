#include "RecoUtils.h"

int RecoUtils::TrueParticleID(const detinfo::DetectorClocksData& clockData, const art::Ptr<recob::Hit>& hit) {
  double particleEnergy = 0;
  int likelyTrackID = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = trackIDs.at(idIt).trackID;
    }
  }
  return likelyTrackID;
}


int RecoUtils::TrueParticleIDFromTotalTrueEnergy(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  std::map<int,double> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
    }
  }

  //Loop over the map and find the track which contributes the highest energy to the hit vector
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = trackIDToEDepMap.begin(); mapIt != trackIDToEDepMap.end(); mapIt++){
    double energy = mapIt->second;
    double trackid = mapIt->first;
    //    std::cout << "Track ID: " << mapIt->first << "has deposited 1/3 this energy: " << mapIt->second << std::endl;
    if (energy > maxenergy){
      maxenergy = energy;
      objectTrack = trackid;

      //All Electrons that deposited charge packets are given the negative of the track id they orginated from but I find the mother just in case this is not true.
      if(trackid < 0){
  simb::MCParticle motherparticle = particleInventory->TrackIdToMotherParticle(trackid);
  objectTrack = motherparticle.TrackId();
      }
    }
  }

  return objectTrack;
}



int RecoUtils::TrueParticleIDFromTotalRecoCharge(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits) {
  // Make a map of the tracks which are associated with this object and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(clockData, hit);
    trackMap[trackID] += hit->Integral();
  }

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      objectTrack  = trackIt->first;
    }
  }
  return objectTrack;
}



int RecoUtils::TrueParticleIDFromTotalRecoHits(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits) {
  // Make a map of the tracks which are associated with this object and the number of hits they are the primary contributor to
  std::map<int,int> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(clockData, hit);
    trackMap[trackID]++;
  }

  // Pick the track which is the primary contributor to the most hits as the 'true track'
  int objectTrack = -99999;
  int highestCount = -1;
  for (std::map<int,int>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCount) {
      highestCount = trackIt->second;
      objectTrack  = trackIt->first;
    }
  }
  return objectTrack;
}



bool RecoUtils::IsInsideTPC(TVector3 position, double distance_buffer){
  double vtx[3] = {position.X(), position.Y(), position.Z()};
  bool inside = false;
  art::ServiceHandle<geo::Geometry> geom;
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

  if (geom->HasTPC(idtpc))
  {
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

    for (size_t c = 0; c < geom->Ncryostats(); c++)
    {
      const geo::CryostatGeo& cryostat = geom->Cryostat(c);
      for (size_t t = 0; t < cryostat.NTPC(); t++)
      {
        const geo::TPCGeo& tpcg = cryostat.TPC(t);
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
      }
    }

    //x
    double dista = fabs(minx - position.X());
    double distb = fabs(position.X() - maxx);
    if ((position.X() > minx) && (position.X() < maxx) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    //y
    dista = fabs(maxy - position.Y());
    distb = fabs(position.Y() - miny);
    if (inside && (position.Y() > miny) && (position.Y() < maxy) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
    //z
    dista = fabs(maxz - position.Z());
    distb = fabs(position.Z() - minz);
    if (inside && (position.Z() > minz) && (position.Z() < maxz) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
  }

  return inside;

}

int RecoUtils::NumberofHitsFromTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  int HitNum = 0;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    std::map<int,float> hitEnergies;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      hitEnergies[TMath::Abs(trackIDEs.at(idIt).trackID)] += trackIDEs.at(idIt).energy;
    }

    //Find which track deposited the most energy.
    int   likelytrack = -9999;
    float MaxEnergy   = -9999;
    for(std::map<int,float>::iterator track_iter=hitEnergies.begin();track_iter!=hitEnergies.end();++track_iter){
      if(track_iter->second > MaxEnergy){
  MaxEnergy = track_iter->second;
  likelytrack = track_iter->first;
      }
    }

    if(likelytrack == TrackID){++HitNum;}

  }
  return HitNum;
}


int RecoUtils::NumberofPrimaryHitsFromTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  int HitNum = 0;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    std::map<int,float> hitEnergies;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      if(trackIDEs.at(idIt).energy < 0.1){continue;}
      hitEnergies[trackIDEs.at(idIt).trackID] += trackIDEs.at(idIt).energy;
    }

    if(hitEnergies.size() > 1){continue;}

    //Find which track deposited the most energy.
    int   likelytrack = -9999;
    float MaxEnergy   = -9999;
    for(std::map<int,float>::iterator track_iter=hitEnergies.begin();track_iter!=hitEnergies.end();++track_iter){
      if(track_iter->second > MaxEnergy){
  MaxEnergy = track_iter->second;
  likelytrack = track_iter->first;
      }
    }

    if(likelytrack == TrackID){++HitNum;}

  }
  return HitNum;
}

int RecoUtils::NumberofPrimaryHitsWithAllTracks(const detinfo::DetectorClocksData& clockData, std::vector<int>& TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  int HitNum = 0;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    std::map<int,float> hitEnergies;
    int track_ids = 0;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {

      if(trackIDEs.at(idIt).energy < 0.3){continue;}
      if(std::find(TrackIDs.begin(),TrackIDs.end(),(int)trackIDEs.at(idIt).trackID) == TrackIDs.end()){continue;}
      ++track_ids;
    }
    if(track_ids == (int) TrackIDs.size()){++HitNum;}
  }
  return HitNum;
}

std::map<geo::PlaneID,int> RecoUtils::NumberofPlaneHitsFromTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<geo::Geometry> geom;

  std::map<geo::PlaneID, int> HitNum_plane;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);

    std::map<int,float> hitEnergies;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      hitEnergies[TMath::Abs(trackIDEs.at(idIt).trackID)] += trackIDEs.at(idIt).energy;
    }

    //Find which track deposited the most energy.
    int   likelytrack = -9999;
    float MaxEnergy   = -9999;
    for(std::map<int,float>::iterator track_iter=hitEnergies.begin();track_iter!=hitEnergies.end();++track_iter){
      if(track_iter->second > MaxEnergy){
  MaxEnergy = track_iter->second;
  likelytrack = track_iter->first;
      }
    }

    if(likelytrack == TrackID){++HitNum_plane[PlaneID];}
  }
  return HitNum_plane;
}


std::map<int,std::map<geo::PlaneID,int> > RecoUtils::NumberofPlaneHitsPerTrack(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<geo::Geometry> geom;

  std::map<int, std::map<geo::PlaneID, int> > HitNum;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();


    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);

    std::map<int,float> hitEnergies;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      hitEnergies[TMath::Abs(trackIDEs.at(idIt).trackID)] += trackIDEs.at(idIt).energy;
    }

    //Find which track deposited the most energy.
    int   likelytrack = -9999;
    float MaxEnergy   = -9999;
    for(std::map<int,float>::iterator track_iter=hitEnergies.begin();track_iter!=hitEnergies.end();++track_iter){
      if(track_iter->second > MaxEnergy){
  MaxEnergy = track_iter->second;
  likelytrack = track_iter->first;
      }
    }

    ++HitNum[likelytrack][PlaneID];
  }
  return HitNum;
}



std::map<geo::PlaneID,int> RecoUtils::NumberofHitsThatContainEnergyDepositedByTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<geo::Geometry> geom;

  //Loop over the planes and create the initial PlaneMap
  std::map<geo::PlaneID,int> HitPlaneMap;
  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){HitPlaneMap[plane_id] = 0;}

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;
    //Find the plane id
    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    //Loop over the IDEs associated to the Hit
    for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      if (TMath::Abs(trackIDEs.at(idIt).trackID) == TrackID) {++HitPlaneMap[PlaneID];}
    }
  }
  return HitPlaneMap;
}


std::map<geo::PlaneID,int> RecoUtils::NumberofHitsThatContainEnergyDepositedByTracks(const detinfo::DetectorClocksData& clockData, std::vector<int> TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits){

  //Loop over the planes and create the initial PlaneMap
  art::ServiceHandle<geo::Geometry> geom;
  std::map<geo::PlaneID,int> HitPlaneMap;
  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){HitPlaneMap[plane_id] = 0;}

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;
    //Find the plane id
    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    //Loop over the IDEs associated to the Hit
    for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      if (std::find(TrackIDs.begin(), TrackIDs.end(),TMath::Abs(trackIDEs.at(idIt).trackID)) != TrackIDs.end()) {++HitPlaneMap[PlaneID];}
    }
  }
  return HitPlaneMap;
}


std::map<geo::PlaneID,int> RecoUtils::NumberofMCWiresHit(int TrackID, const std::vector<art::Ptr<sim::SimChannel> >& simchannels){


  //I don't think there is a way to go from TrackID to a list of TrackIDEs. So Instead loop through all the sim channels. Where there is a track IDE fromt th track ID. Count it.
  art::ServiceHandle<geo::Geometry> geom;
  std::map<geo::PlaneID,int> WirePlaneMap;

  int breaking_int = 0;

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){WirePlaneMap[plane_id] = 0;}

  // Loop over SimChannel
  for(size_t simch_index=0; simch_index<simchannels.size(); ++simch_index) {

    //Get the specific simchanel
    const art::Ptr<sim::SimChannel> simch_ptr = simchannels[simch_index];

    //Get the plane.
    raw::ChannelID_t Channel = simch_ptr->Channel();
    std::vector<geo::WireID> Wire = geom->ChannelToWire(Channel);
    //sbnd is one channel per wire so the vector should be size one.
    geo::PlaneID  PlaneID = (Wire[0]).planeID();

    //Get the TDCIDEMap (Charge vs Time)
    auto tdc_ide_map = simch_ptr->TDCIDEMap();

    //Loop through the map
    for(auto const& tdc_ide_pair : tdc_ide_map) {

      //Get the IDEs associated to the TDC?
      auto const& ide_v = tdc_ide_pair.second;

      //Loop over the IDEs and add the energy.
      for(auto const& ide : ide_v) {
    if(TMath::Abs(ide.trackID) == TrackID)
    {++WirePlaneMap[PlaneID];breaking_int=1;break;}
      }

      if(breaking_int==1){breaking_int=0;break;}
    }
  }
  return WirePlaneMap;
}

std::map<int,std::map<geo::PlaneID,int> > RecoUtils::NumberofMCWiresHitMap(const std::vector<art::Ptr<sim::SimChannel> >& simchannels){



  //I don't think there is a way to go from TrackID to a list of TrackIDEs. So Instead loop through all the sim channels. Where there is a track IDE fromt th track ID. Count it.
  art::ServiceHandle<geo::Geometry> geom;
  std::map<int,std::map<geo::PlaneID,int> >  TrackWireHitMap;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  const sim::ParticleList& particles = particleInventory->ParticleList();

  std::vector<float> UsedTracks;
  UsedTracks.reserve(particles.size());

  // Loop over SimChannel
  for(size_t simch_index=0; simch_index<simchannels.size(); ++simch_index) {

    //Get the specific simchanel
    const art::Ptr<sim::SimChannel> simch_ptr = simchannels[simch_index];

    //Get the plane.
    raw::ChannelID_t Channel = simch_ptr->Channel();
    std::vector<geo::WireID> Wire = geom->ChannelToWire(Channel);
    //sbnd is one channel per wire so the vector should be size one.
    geo::PlaneID  PlaneID = (Wire[0]).planeID();

    //Get the TDCIDEMap (Charge vs Time)
    auto tdc_ide_map = simch_ptr->TDCIDEMap();

    //Loop through the map
    for(auto const& tdc_ide_pair : tdc_ide_map) {

      //Get the IDEs associated to the TDC?
      auto const& ide_v = tdc_ide_pair.second;

      //Loop over the IDEs and add the energy.
      for(auto const& ide : ide_v) {
    if(std::find(UsedTracks.begin(), UsedTracks.end(),TMath::Abs(ide.trackID)) == UsedTracks.end()){
    ++TrackWireHitMap[TMath::Abs(ide.trackID)][PlaneID];
    UsedTracks.push_back(TMath::Abs(ide.trackID));
  }
      }
    }
    UsedTracks.clear();
  }
  return TrackWireHitMap;
}



//Gets the total energy deposited by a track  by looping through the MC info on each wire.
float RecoUtils::TrueEnergyDepositedFromMCTrack(int TrackID, const std::vector<art::Ptr<sim::SimChannel> > &simchannels ){

  art::ServiceHandle<geo::Geometry> geom;

  double total_energy = 0;
  // Loop over SimChannel
  for(size_t simch_index=0; simch_index<simchannels.size(); ++simch_index) {

    //Get the specific simchanel
    const art::Ptr<sim::SimChannel> simch_ptr = simchannels[simch_index];

    //Get the plane.
    raw::ChannelID_t Channel = simch_ptr->Channel();
    std::vector<geo::WireID> Wire = geom->ChannelToWire(Channel);
    //sbnd is one channel per wire so the vector should be size one.
    int  PlaneID = (Wire[0]).Plane;

    //Get the TDCIDEMap (Charge vs Time)
    auto tdc_ide_map = simch_ptr->TDCIDEMap();

    //Loop through the map
    for(auto const& tdc_ide_pair : tdc_ide_map) {

      //Get the IDEs associated to the TDC?
      auto const& ide_v = tdc_ide_pair.second;

      //Loop over the IDEs and add the energy. Only count from the collection plane.
      for(auto const& ide : ide_v) {
  if(TMath::Abs(ide.trackID) == TrackID && PlaneID == 2){
    total_energy +=  ide.energy;
  }
      }
    }
  }

  return total_energy;
}


//Gets the total energy deposited by a track  by looping through the MC info on each wire.
std::map<int,float> RecoUtils::TrueEnergyDepositedFromMCTracks(const std::vector<art::Ptr<sim::SimChannel> > &simchannels ){

  art::ServiceHandle<geo::Geometry> geom;


  std::map<int,float> total_energies;

  // Loop over SimChannel
  for(size_t simch_index=0; simch_index<simchannels.size(); ++simch_index) {

    //Get the specific simchanel
    const art::Ptr<sim::SimChannel> simch_ptr = simchannels[simch_index];

    //Get the plane.
    raw::ChannelID_t Channel = simch_ptr->Channel();
    std::vector<geo::WireID> Wire = geom->ChannelToWire(Channel);
    //sbnd is one channel per wire so the vector should be size one.
    int  PlaneID = (Wire[0]).Plane;

    //Get the TDCIDEMap (Charge vs Time)
    auto tdc_ide_map = simch_ptr->TDCIDEMap();

    //Loop through the map
    for(auto const& tdc_ide_pair : tdc_ide_map) {

      //Get the IDEs associated to the TDC?
      auto const& ide_v = tdc_ide_pair.second;

      //Loop over the IDEs and add the energy. Only count from the collection plane.
      for(auto const& ide : ide_v) {
  if(PlaneID == 2){
    total_energies[TMath::Abs(ide.trackID)] +=  ide.energy;
  }
      }
    }
  }

  return total_energies;
}



float RecoUtils::TotalEnergyDepinHits(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, int Plane){

  float DepEnergy = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  //Loop over the hits and find the IDEs
  for(std::vector< art::Ptr<recob::Hit> >::const_iterator hitIt=hits.begin(); hitIt!=hits.end(); ++hitIt){

    //Get the plane ID
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    if(PlaneID != Plane){continue;}

    //Split the Hit into its IDE for each track it associates with.
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, (*hitIt));

    for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){

      //Find the true total energy deposited in a set of hits.
      DepEnergy += trackIDEs.at(idIt).energy;

    }
  }//Hit Loop

  return DepEnergy;
}


float RecoUtils::TotalEnergyDepinHitsFromTrack(const detinfo::DetectorClocksData& clockData,
    const std::vector<art::Ptr<recob::Hit> >& hits, int TrackID, int Plane){

  float DepEnergy = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  //Loop over the hits and find the IDEs
  for(std::vector< art::Ptr<recob::Hit> >::const_iterator hitIt=hits.begin(); hitIt!=hits.end(); ++hitIt){

    //Get the plane ID
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    if(PlaneID != Plane){continue;}

    //Split the Hit into its IDE for each track it associates with.
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, (*hitIt));

    for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){

      //Add up the contribution from the trackID.
      if(TMath::Abs(trackIDEs.at(idIt).trackID) == TrackID){
        DepEnergy += trackIDEs.at(idIt).energy;
      }
    }
  }//Hit Loop

  return DepEnergy;
}


