#include "BlipUtils.h"

//########################
// Functions
//########################

// Function to check if there was a SimChannel made for a hit (useful when checking for noise hits)
bool BlipUtils::DoesHitHaveSimChannel( std::vector<const sim::SimChannel*> chans, art::Ptr<recob::Hit> const& hit){
  bool isMatch = false;
  for(size_t sc = 0; sc < chans.size(); ++sc){
    if( chans[sc]->Channel() == hit->Channel() ) { isMatch = true; break; }
  }
  return isMatch;
}

// Function to calculate purity of a collection of hits
void BlipUtils::HitsPurity(
//  detinfo::DetectorClocksData const& clockData,
  std::vector< art::Ptr<recob::Hit> > const& hits, 
  int& trackid, 
  float& purity, 
  double& maxe 
  ) {
  
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  trackid = -1;
  purity = -1;

  std::map<int,double> trkide;
  for(size_t h = 0; h < hits.size(); ++h){
    art::Ptr<recob::Hit> hit = hits[h];
    std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToEveTrackIDEs(clockData, hit);
    for(size_t e = 0; e < eveIDs.size(); ++e){
      trkide[eveIDs[e].trackID] += eveIDs[e].energy;
    }
  }

  maxe = -1;
  double tote = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      trackid = ii->first;
    }
  }

  if (tote>0) purity = maxe/tote;

}


// This function calculates the leading MCParticle ID contributing to a hit and the
// fraction of that hit's energy comes from that particle.
void  BlipUtils::HitTruth(
  //detinfo::DetectorClocksData const& clockData, 
  art::Ptr<recob::Hit> const& hit, 
	int& truthid, 
  float& truthidEnergyFrac, 
  float& energy, 
  float& numElectrons){
  

  // Get associated sim::TrackIDEs for this hit
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
  if( !trackIDEs.size() ) return;
  // Loop through and find the leading TrackIDE, and keep
  // track of the total energy of ALL IDEs.
  float maxe = 0;
  float bestfrac = 0;
  float bestid = 0;
  float ne = 0;
  for(size_t i = 0; i < trackIDEs.size(); ++i){
    ne += trackIDEs[i].numElectrons;
    if( trackIDEs[i].energy > maxe ) {
      maxe = trackIDEs[i].energy;
      bestfrac = trackIDEs[i].energyFrac;
      bestid = trackIDEs[i].trackID;
    }
  }
  // Save the results
  truthid = bestid;
  truthidEnergyFrac = bestfrac;
  energy = maxe;
  numElectrons = ne;

}
//
// helper function to get TrackID (returns false if no match found)
bool BlipUtils::HitTruthId( 
  //detinfo::DetectorClocksData const& clockData, 
  art::Ptr<recob::Hit> const& hit, 
  int& mcid) 
{
  mcid = std::numeric_limits<int>::lowest();
  float dummy1;
  float dummy2;
  float dummy3;
  //HitTruth(clockData,hit,mcid,dummy1,dummy2,dummy3);
  HitTruth(hit,mcid,dummy1,dummy2,dummy3);
  if( mcid > std::numeric_limits<int>::lowest() ) return true;
  else return false;
}


// Get MCTruth associated with TrackID using a try bracket to avoid
// fatal exceptions (return false if no match or exception thrown)
bool BlipUtils::TrackIdToMCTruth( int const trkID, art::Ptr<simb::MCTruth>& mctruth )
{ 
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  bool matchFound = false;
  try {
    mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
    matchFound = true;
  } catch(...) {
    std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
    matchFound = false;
  }
  return matchFound;
}
