#include "ShowerUtils.h"

std::map<int,std::vector<int> > ShowerUtils::GetShowerMothersCandidates(std::map<int,const simb::MCParticle*>& trueParticles){

  std::map<int,std::vector<int> > ShowerMotherCandidates;

  for(auto const& particle_iter: trueParticles){

    const simb::MCParticle* particle = particle_iter.second;

    //Check we are are shower particle.
    if(TMath::Abs(particle->PdgCode()) != 11 && TMath::Abs(particle->PdgCode()) != 22){continue;}

    //Check the mother is not a shower particle
    if(trueParticles.find(particle->Mother()) != trueParticles.end()){
      if(TMath::Abs(trueParticles[particle->Mother()]->PdgCode()) == 11
          || TMath::Abs(trueParticles[particle->Mother()]->PdgCode()) == 22){

        //I propose that a daughter id always come and after a mother id the daughter can be added here and be in order.
        int mother_id = particle->Mother();
        int particle_temp = mother_id;
        //Match particle with mother
        while(mother_id != 0){
          particle_temp = trueParticles[mother_id]->Mother();
          if(trueParticles.find(particle_temp) == trueParticles.end()){break;}
          if(TMath::Abs(trueParticles[particle_temp]->PdgCode()) != 11
              && TMath::Abs(trueParticles[particle_temp]->PdgCode()) != 22){break;}
          mother_id = particle_temp;
        }

        //Add to the mother chain.
        ShowerMotherCandidates[mother_id].push_back(particle->TrackId());

        continue;
      }
    }
    ShowerMotherCandidates[particle->TrackId()].push_back(particle->TrackId());
  }
  return ShowerMotherCandidates;
}

std::pair<int,double> ShowerUtils::TrueParticleIDFromTrueChain(
    const detinfo::DetectorClocksData& clockData,
    std::map<int,std::vector<int>> &ShowersMothers,
    const std::vector<art::Ptr<recob::Hit> >& hits,
    int planeid) {

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int,double> trackIDToEDepMap;
  std::map<int,double> trackIDTo3EDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDTo3EDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      if(PlaneID == planeid){trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;}
    }
  }

  //Find the energy for each showermother.
  std::map<int,double> MotherIDtoEMap;
  std::map<int,double> MotherIDto3EMap;
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
    for(std::vector<int>::iterator daughter=(showermother->second).begin(); daughter!=(showermother->second).end(); ++daughter){
      MotherIDtoEMap[showermother->first]  +=  trackIDToEDepMap[*daughter];
      MotherIDto3EMap[showermother->first] +=  trackIDTo3EDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candidate by identifying which true shower deposited the most energy in the hits.
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = MotherIDto3EMap.begin(); mapIt != MotherIDto3EMap.end(); mapIt++){
    double energy_three = mapIt->second;
    double trackid = mapIt->first;
    if (energy_three > maxenergy){
      maxenergy = energy_three;
      objectTrack = trackid;
    }
  }

  //If the none of the shower mother deposited no energy then we cannot match this.
  if(maxenergy == 0){
    return std::make_pair(-99999,-99999);
  }

  return std::make_pair(objectTrack,MotherIDtoEMap[objectTrack]);
}

std::map<geo::PlaneID,int> ShowerUtils::NumberofWiresHitByShower(
    const detinfo::DetectorClocksData& clockData,
    std::vector<int> &TrackIDs,
    const std::vector<art::Ptr<recob::Hit> >& hits){

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  std::vector<geo::WireID> WiresUsed;
  std::map<geo::PlaneID,int> HitWirePlaneMap;

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    //Find the wire and  plane id
    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    //Check to see if the wire is already been continued.
    if(std::find(WiresUsed.begin(),WiresUsed.end(),wireid) != WiresUsed.end()){continue;}

    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      if(std::find(TrackIDs.begin(), TrackIDs.end(),TMath::Abs(trackIDs.at(idIt).trackID)) != TrackIDs.end()){
        WiresUsed.push_back(wireid);
        ++HitWirePlaneMap[PlaneID];
        break;
      }
    }
  }
  return HitWirePlaneMap;
}

void ShowerUtils::CutShowerMothersByE(std::map<int,std::vector<int> >& ShowersMothers,
    std::map<int,const simb::MCParticle*>& trueParticles, float& EnergyCut){

  //Time to cut the true showers and make sure they are a shower.
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();){

    //I've read that pair production starts to dominate at around ~100 MeV so to find how many showers we expect loop over the mother particle. Pi0=143.97 MeV min gammas = 71.985 MeV which is greater than that from electrons at ~100MeV so pi0 should always shower? So cut on anything below 100MeV in energy.

    //It ain't a shower I'm interested in if it didn't start with a pi0 or electron...probably.
    const simb::MCParticle *motherparticle = trueParticles[showermother->first];

    if(motherparticle->E() < EnergyCut){
      ShowersMothers.erase(showermother++);
    }
    else{
      ++showermother;
    }
  }
  return;
}

void ShowerUtils::CutShowerMothersByDensity(
    const detinfo::DetectorClocksData& clockData,
    std::map<int,std::vector<int> >& ShowersMothers,
    std::map<int,const simb::MCParticle*>& trueParticles,
    std::vector<art::Ptr<recob::Hit> >& hits,
    float& fDensityCut){

  //Time to cut the true showers and make sure they are a shower.
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();){

    //using the RecoUtil function calculate the number of hits that see a charge deposition from the track.
    std::map<geo::PlaneID,int> Hit_num_map = RecoUtils::NumberofHitsThatContainEnergyDepositedByTracks(clockData, showermother->second, hits);

    //Calculaute the number of wires hit.
    std::map<geo::PlaneID,int> Wire_num_map = ShowerUtils::NumberofWiresHitByShower(clockData, showermother->second, hits);

    int high_density_plane=0;

    //Compare hit density on each plane;
    for(std::map<geo::PlaneID,int>::iterator Hitnum_iter=Hit_num_map.begin(); Hitnum_iter!=Hit_num_map.end(); ++Hitnum_iter){

      if(Wire_num_map[Hitnum_iter->first] == 0){continue;}
      double Hit_num = (Hitnum_iter->second);
      double Wire_num = Wire_num_map[Hitnum_iter->first];
      double Hit_Density = Hit_num/Wire_num;

      if(Hit_Density > fDensityCut){
        ++high_density_plane;
        break;
      }
    }

    //If none of the planes have a density high than the cut remove the event.
    if(high_density_plane == 0){
      ShowersMothers.erase(showermother++);
    }
    else{
      ++showermother;
    }
  }
  return;
}

void ShowerUtils::RemoveNoneContainedParticles(std::map<int,std::vector<int> >&  ShowersMothers,
    std::map<int,const simb::MCParticle*>& trueParticles, std::map<int,float>& MCTrack_Energy_map){

  //Calculate the deposited and the true energy in the tpc by the particle chain.
  for(std::map<int,std::vector<int> >::iterator ShowerMother=ShowersMothers.begin(); ShowerMother!=ShowersMothers.end();){

    //Loop over the daughters
    float deposited_energy = 0;
    float sim_energy       = trueParticles[ShowerMother->first]->E()*1000;
    for(auto const& Daughter: ShowerMother->second){
      deposited_energy += MCTrack_Energy_map[Daughter];
    }

    //If over 90% of the shower energy is seen on the wires in truth. It is contained.
    if(deposited_energy/sim_energy < 0.9){
      ShowersMothers.erase(ShowerMother++);
    }
    else{
      ++ShowerMother;
    }
  }
  return;
}
