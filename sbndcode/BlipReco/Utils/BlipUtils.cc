#include "BlipUtils.h"

//########################
// Functions
//########################

namespace BlipUtils {

  // Function to check if there was a SimChannel made for a hit (useful when checking for noise hits)
  bool DoesHitHaveSimChannel( std::vector<const sim::SimChannel*> chans, art::Ptr<recob::Hit> const& hit){
    bool isMatch = false;
    for(size_t sc = 0; sc < chans.size(); ++sc){
      if( chans[sc]->Channel() == hit->Channel() ) { isMatch = true; break; }
    }
    return isMatch;
  }
  
  // Function to calculate purity of a collection of hits
  void HitsPurity(
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
  void  HitTruth(
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
  bool HitTruthId( 
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
  bool TrackIdToMCTruth( int const trkID, art::Ptr<simb::MCTruth>& mctruth )
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
  
  
  // Length of reconstructed track, trajectory by trajectory.
  double PathLength(const recob::Track& track)
  {
    return track.Length();
  }
  // Length of MC particle, trajectory by trajectory.
  double PathLength(const simb::MCParticle& part, TVector3& start, TVector3& end)
  {
    // Get geometry.
    art::ServiceHandle<geo::Geometry> geom;
  
    // Get active volume boundary.
    double xmin = -2.0 * geom->DetHalfWidth() - 1e-8;
    double xmax = 2.0 * geom->DetHalfWidth() + 1e-8;
    double ymin = -geom->DetHalfHeight() -1e-8;
    double ymax = geom->DetHalfHeight() + 1e-8;
    double zmin = 0. -1e-8;
    double zmax = geom->DetLength() + 1e-8;
    
    // Get number traj points
    int n = part.NumberTrajectoryPoints();
    if( n <= 1 ) return 0.;
    
    double  L	= 0.;
    bool	  first	= true; 
  
    // Loop over points (start with 2nd)
    for(int i = 1; i < n; ++i) {
  
      TVector3 p1(part.Vx(i),part.Vy(i),part.Vz(i));
      
      if(	  p1.X() >= xmin && p1.X() <= xmax
        &&  p1.Y() >= ymin && p1.Y() <= ymax
        &&  p1.Z() >= zmin && p1.Z() <= zmax ) {
      
        TVector3 p0(part.Vx(i-1),part.Vy(i-1),part.Vz(i-1));
        
        if(first)	start = p1; 
        else L += (p1-p0).Mag();
        first = false;
        end   = p1;
      }
    }
  
    return L;
  }
  double PathLength(const simb::MCParticle& part){
    TVector3 a,b;
    return PathLength(part,a,b);
  }

}
