#include "BlipUtils.h"

const float rmsWidthFact = 1.0;

namespace sbnd::BlipUtils {
  
  //============================================================================
  // Find total visible energy deposited in the LAr, and number of electrons drifted
  // to wires, by looking at the SimChannels for all three planes
  void CalcTotalDep(float& energy, float& ne){
    energy  = 0;
    ne      = 0;
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) {
      for(auto const &tdcide : chan->TDCIDEMap() ) {
        for(const auto& ide : tdcide.second) {
          energy += ide.energy/art::ServiceHandle<geo::Geometry>()->Nplanes();
          ne += ide.numElectrons/art::ServiceHandle<geo::Geometry>()->Nplanes();
        }
      }
    }
  }
  
  //============================================================================
  // Get the total energy deposited by this particle by looking at IDEs from 3 planes.
  void CalcPartDep(int trackID, float& energy, float& ne){
    if( energy< 0 ) energy = 0;
    if( ne    < 0 ) ne = 0;
    std::cout<<"Calculating energy dep for trackID "<<trackID<<"\n";
    int nPlanes = 0;
    float totalE_particle = 0, totalne_particle = 0;
    for(const geo::View_t view : {geo::kU, geo::kV, geo::kW} ) {
      std::vector<const sim::IDE* > ides 
        = art::ServiceHandle<cheat::BackTrackerService>()->TrackIdToSimIDEs_Ps(trackID, view);
      if( ides.size() ) nPlanes++;
      for (auto ide: ides) {
        totalE_particle += ide->energy;
        totalne_particle += ide->numElectrons;
      }
    }
    if(nPlanes) { totalE_particle /= nPlanes; totalne_particle /= nPlanes; }
    energy  += totalE_particle;
    ne      += totalne_particle;
  }


  //====================================================================
  // Make a TrueBlip object from a G4 particle, grouping contiguous particles
  // in with the blip.
  TrueBlip MakeTrueBlip(int trackID){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    simb::MCParticle part = *plist.at(trackID);
    
    TrueBlip tb;

    // If this is an electron that came from another electron, it would 
    // have already been grouped as part of the contiguous "blip" previously.
    if( part.PdgCode() == 11 && (part.Process() == "eIoni" || part.Process() == "muIoni" || part.Process() == "hIoni") ) return tb;
    
    // If this is a photon or neutron, don't even bother!
    if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) return tb;

    // Only count protons < 3 MeV KE
    float edep, ne; 
    CalcPartDep(trackID,edep,ne);
    if( part.PdgCode() == 2212 && edep > 3. ) return tb;

    // Create the new blip
    GrowTrueBlip(part,tb);

    // We want to loop through any contiguous electrons (produced
    // with process "eIoni") and add the energy they deposit into this blip.
    if( part.NumberDaughters() ) {
      for(auto it = std::next(plist.find(trackID)); it != plist.end(); it++){
        simb::MCParticle p = *(it)->second;
        if( p.PdgCode() != 2112 && (p.Process() == "eIoni" || p.Process() == "muIoni" || p.Process() == "hIoni") ){
          if( IsAncestorOf(p.TrackId(),trackID,true) ) GrowTrueBlip(p,tb);
        }
      }
    }

    // Ensure this blip has at least some energy, and
    // that it doesn't extend beyond ~3 wires
//    if( tb.Energy && tb.Length < 1. ) tb.isValid = true;
    if( tb.Energy) tb.isValid = true;
    return tb;
  }

  //====================================================================
  void GrowTrueBlip(simb::MCParticle const& part, TrueBlip& tblip){
    
    
    // Skip neutrons, photons
    if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) return;

    // Get particle path length and start/end position *inside TPC*
    TVector3 start, end;
    float len = PathLength(part,start,end);

    if( !len ) return;

    // If this is a new blip, initialize
    if( tblip.G4IDs.size() == 0 ) {
      tblip.StartPoint  = start;
      //double posArray[3]= { tblip.Position.X(), tblip.Position.Y(), tblip.Position.Z() };
      //tblip.TPC         = art::ServiceHandle<geo::Geometry>()->PositionToTPC({start.X(),start.Y(),start.Z()}).ID().TPC;
    }
    
    
    float edep = 0, ne = 0;
    CalcPartDep(part.TrackId(),edep,ne);
    tblip.G4IDs       .push_back(part.TrackId());
    tblip.PDGs       .push_back(part.PdgCode());
    tblip.Energy      += edep;
    tblip.NumElectrons+= ne;
    tblip.EndPoint    = end;
    tblip.Position    = (tblip.StartPoint + end)*0.5;
    tblip.Length      += len; 
    if(edep > tblip.LeadEnergy ) {
      tblip.LeadEnergy = edep;
      tblip.LeadG4ID = part.TrackId();
      tblip.LeadG4PDG = part.PdgCode();
    }
  }
  
  //====================================================================
  // Merge blips that are close
  void MergeTrueBlips(std::vector<TrueBlip>& vtb, float dmin){
    if( dmin <= 0 ) return;
    std::vector<TrueBlip> vtb_merged;
    std::vector<bool> isGrouped(vtb.size(),false);
    for(size_t i=0; i<vtb.size(); i++){
      if( isGrouped.at(i) ) continue;
      else isGrouped.at(i) = true;
      auto& blip_i = vtb.at(i);
      for(size_t j=i+1; j<vtb.size(); j++){
        if( isGrouped.at(j) ) continue;
        auto const& blip_j = vtb.at(j);
//        if( blip_i.TPC != blip_j.TPC ) continue;
        float d = (blip_i.Position-blip_j.Position).Mag();
        if( d < dmin ) {
          isGrouped.at(j) = true;
          blip_i.Energy   += blip_j.Energy;
          blip_i.EndPoint = blip_j.EndPoint;
          blip_i.Position = (blip_i.EndPoint+blip_i.StartPoint)*0.5;
          blip_i.Length   = (blip_i.EndPoint-blip_i.StartPoint).Mag();
          for(size_t kk=0; kk<blip_j.G4IDs.size(); kk++)
            blip_i.G4IDs.push_back(blip_j.G4IDs.at(kk)); 
          for(size_t kk=0; kk<blip_j.PDGs.size(); kk++)
            blip_i.PDGs.push_back(blip_j.PDGs.at(kk));
          if( blip_j.LeadEnergy > blip_i.LeadEnergy ) {
            blip_i.LeadEnergy = blip_j.LeadEnergy;
            blip_i.LeadG4ID = blip_j.LeadG4ID;
            blip_i.LeadG4PDG = blip_j.LeadG4PDG;
          }
        }//d < dmin
      }//loop over blip_j
      vtb_merged.push_back(blip_i);
    }
    vtb.clear();
    vtb = vtb_merged;
  }
 
  //=================================================================
  HitClust MakeHitClust(HitInfo const& hitinfo){
    art::Ptr<recob::Hit> hit = hitinfo.hit;
    HitClust hc;
    hc.LeadHit      = hit;
    hc.TPC          = hit->WireID().TPC;
    hc.Plane        = hit->WireID().Plane;
    hc.LeadHitG4ID  = hitinfo.g4id;
    hc.G4IDs        .insert(hitinfo.g4ids.begin(), hitinfo.g4ids.end());
    hc.HitIDs       .insert(hitinfo.hitid);
    hc.Wires       .insert(hit->WireID().Wire);
    hc.LeadHitWire  = hit->WireID().Wire;
    hc.Charge       = hitinfo.qcoll;
    hc.LeadHitCharge= hc.Charge;
    hc.mapWireCharge[hit->WireID().Wire] = hitinfo.qcoll;
    hc.mapTimeCharge[hitinfo.driftTicks] = hitinfo.qcoll;
    hc.LeadHitTime  = hitinfo.driftTicks;
    hc.Time         = hitinfo.driftTicks;
    hc.WeightedTime = hitinfo.driftTicks;
    hc.StartTime    = hitinfo.driftTicks - hit->RMS();
    hc.EndTime      = hitinfo.driftTicks + hit->RMS();
    hc.StartWire    = hc.LeadHitWire;
    hc.EndWire      = hc.LeadHitWire;
    hc.isValid      = true;
    hc.isMerged     = false;
    hc.isMatched    = false;
    //hc.isMatchedPl[hc.Plane] = true;
    return hc;
  }

  //=================================================================
  void GrowHitClust(HitInfo const& hitinfo, HitClust& hc){
    art::Ptr<recob::Hit> hit = hitinfo.hit;
    if( (int)hit->WireID().TPC   != hc.TPC ) return;
    if( (int)hit->WireID().Plane != hc.Plane ) return;
    if( hc.HitIDs.find(hitinfo.hitid) != hc.HitIDs.end() ) return;
    hc.G4IDs .insert(hitinfo.g4ids.begin(), hitinfo.g4ids.end());
    hc.HitIDs     .insert(hitinfo.hitid);
    hc.Wires      .insert(hit->WireID().Wire);
    hc.StartWire  = *hc.Wires.begin();
    hc.EndWire    = *hc.Wires.rbegin();
    hc.Charge     += hitinfo.qcoll;
    hc.mapWireCharge[hit->WireID().Wire] += hitinfo.qcoll;
    hc.mapTimeCharge[hitinfo.driftTicks] += hitinfo.qcoll;
    if( hitinfo.qcoll > hc.LeadHitCharge ) {
      hc.LeadHit = hit;
      hc.LeadHitCharge = hitinfo.qcoll;
      hc.LeadHitG4ID = hitinfo.g4id;
      hc.LeadHitTime = hitinfo.driftTicks;
      hc.LeadHitWire = hitinfo.hit->WireID().Wire;
    }
    hc.StartTime  = std::min(hc.StartTime,hitinfo.driftTicks - hit->RMS());
    hc.EndTime    = std::max(hc.EndTime,  hitinfo.driftTicks + hit->RMS());
    hc.Time       = (hc.StartTime+hc.EndTime)/2.;
    // recalculate charge-weighted time
    hc.WeightedTime = 0;
    for(auto const& tq : hc.mapTimeCharge ) hc.WeightedTime += (tq.first * tq.second)/hc.Charge;
  }

  //=================================================================
  HitClust MergeHitClusts(HitClust& hc1, HitClust& hc2){
    HitClust hc = hc1;
    if( (hc1.TPC != hc2.TPC)||(hc1.Plane != hc2.Plane) ) return hc;
    hc1.isMerged = true;
    hc2.isMerged = true;
    hc.G4IDs.insert(hc2.G4IDs.begin(), hc2.G4IDs.end());
    hc.HitIDs   .insert(hc2.HitIDs.begin(),     hc2.HitIDs.end());
    hc.Wires    .insert(hc2.Wires.begin(),      hc2.Wires.end());
    hc.StartWire  = *hc.Wires.begin();
    hc.EndWire    = *hc.Wires.rbegin();
    hc.Charge   += hc2.Charge;
    hc.mapWireCharge.insert(hc2.mapWireCharge.begin(),hc2.mapWireCharge.end());
    hc.mapTimeCharge.insert(hc2.mapTimeCharge.begin(),hc2.mapTimeCharge.end());
    if( hc2.LeadHitCharge > hc.LeadHitCharge ) {
      hc.LeadHit = hc2.LeadHit;
      hc.LeadHitCharge = hc2.LeadHitCharge;
      hc.LeadHitWire = hc2.LeadHitWire;
      hc.LeadHitTime = hc2.LeadHitTime;
      hc.LeadHitG4ID = hc2.LeadHitG4ID;
    }
    hc.StartTime  = std::min(hc.StartTime,hc2.StartTime);
    hc.EndTime    = std::max(hc.EndTime,hc2.EndTime);
    hc.Time       = (hc.StartTime+hc.EndTime)/2.;
    // recalculate charge-weighted time
    hc.WeightedTime = 0;
    for(auto const& tq : hc.mapTimeCharge ) hc.WeightedTime += (tq.first * tq.second)/hc.Charge;
    return hc;
  }

  //=================================================================
  Blip MakeBlip( std::vector<HitClust> const& hcs){
    Blip newblip;
    
    // Must be 1-3 clusts (no more, no less!)
    if( hcs.size() > 3  || hcs.size() < 1 ) return newblip;

    // All hits must be in same TPC, and no 2 planes can be repeated
    std::set<int> planeIDs;
    for(size_t i=0; i<hcs.size(); i++) {
      planeIDs.insert(hcs[i].Plane);
      for(size_t j=i+1; j<hcs.size(); j++){
        if( hcs[i].Plane == hcs[j].Plane )  { return newblip; }
        if( hcs[i].TPC   != hcs[j].TPC )    { return newblip; }
      }
    }
    
    int nPlanes = planeIDs.size();
    int TPC = hcs[0].TPC;

    // Calculate mean time and X (assuming in-time deposition)
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    float samplingPeriod = clockData.TPCClock().TickPeriod();
    float t = 0, x = 0;
    for(auto hc : hcs ) {
      t += samplingPeriod * hc.Time/float(hcs.size());
      x += detProp.ConvertTicksToX(hc.LeadHit->PeakTime(),hc.Plane,TPC,0)/float(hcs.size());
    }
    newblip.DriftTime = t;

    // Look for valid wire intersections and calculate
    // the mean Y/Z position from these
    std::vector<TVector3> wirex;
    for(size_t i=0; i<hcs.size(); i++) {
      newblip.ClustIDs.insert(hcs[i].ID);
      newblip.HitIDs.insert(hcs[i].HitIDs.begin(), hcs[i].HitIDs.end());
      int pli = hcs[i].Plane;
      for(size_t j=i+1; j<hcs.size(); j++){
        int plj = hcs[j].Plane;
        double y,z;
        bool match3d = art::ServiceHandle<geo::Geometry>()->ChannelsIntersect(
          hcs[i].LeadHit->Channel(),hcs[j].LeadHit->Channel(),y,z);
        if( match3d ) {
//          newblip.WireIntersect[pli] = true;
//          newblip.WireIntersect[plj] = true;
          newblip.Charge[pli] = hcs[i].Charge;
          newblip.Charge[plj] = hcs[j].Charge;
          newblip.Planes[pli] = true;
          newblip.Planes[plj] = true;
          TVector3 a(x,y,z);
          wirex.push_back(a);
        }
      }
    }
    if( wirex.size() ){
      TVector3 vecmean;
      for(size_t i=0; i<wirex.size(); i++) vecmean += wirex[i] * (1./wirex.size());
      newblip.Position = vecmean;
      newblip.TPC = TPC;
      newblip.NPlanes = nPlanes;
      newblip.isValid = true;
      // find max difference
      float max = -9;
      if( wirex.size() > 1 ) {
        for(size_t i=0; i<wirex.size(); i++){
          for(size_t j=i+1; j<wirex.size(); j++){
            float d = (wirex[i]-wirex[j]).Mag();
            if( d > max ) max = d;
          }
        }
      }
      newblip.MaxIntersectDiff = max;
    }
    return newblip;
  }



  //====================================================================
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    if( particleID == ancestorID ) return true;
    if( !plist.HasParticle(ancestorID) ) return false;
    while( particleID > ancestorID ) {
      simb::MCParticle p = pi_serv->TrackIdToParticle(particleID);
      if( !plist.HasParticle(p.Mother() ) )                       { return false; }   
      simb::MCParticle pM = pi_serv->TrackIdToParticle(p.Mother());
      if      ( pM.TrackId() == ancestorID )                      { return true;  }
      else if ( breakAtPhots == true && pM.PdgCode() == 22 )      { return false; }
      else if ( pM.Process() == "primary" || pM.TrackId() == 1 )  { return false; }
      else if ( pM.Mother() == 0 )                                { return false; }
      else    { particleID = pM.TrackId(); }              
    }
    return false;
  }

  //====================================================================
  bool DoHitsOverlap(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2){
    if( hit1->WireID() != hit2->WireID() ) return false;
    float t1 = hit1->PeakTime();
    float t2 = hit2->PeakTime();
    float sig = std::max(hit1->RMS(),hit2->RMS());
    if( fabs(t1-t2) < rmsWidthFact*sig ) return true;
    else return false;
  }
  
  //====================================================================
  bool DoHitClustsOverlap(HitClust const& hc1, HitClust const& hc2){
    // only match across different wires in same TPC
    if( hc1.TPC != hc2.TPC    ) return false;
    float x1 = hc1.StartTime;
    float x2 = hc1.EndTime;
    float y1 = hc2.StartTime;
    float y2 = hc2.EndTime;
    if( x1 <= y2 && y1 <= x2 ) return true;
    else return false;
  }
  bool DoHitClustsOverlap(HitClust const& hc1, float t1, float t2 ){
    HitClust hc2;
    hc2.TPC = hc1.TPC;
    hc2.StartTime = t1;
    hc2.EndTime = t2;
    return DoHitClustsOverlap(hc1,hc2);
  }

  //====================================================================
  bool DoChannelsIntersect(int ch1, int ch2 ){
    double y,z;
    return art::ServiceHandle<geo::Geometry>()->ChannelsIntersect(ch1,ch2,y,z);
  }
  
  //====================================================================
  bool DoHitClustsMatch(HitClust const& hc1, HitClust const& hc2, float minDiffTicks = 2){
    if( fabs(hc1.WeightedTime-hc2.WeightedTime) < minDiffTicks ) return true;
    else return false;
  }


  //=====================================================================
  // Function to check if there was a SimChannel made for a hit (useful when checking for noise hits)
  bool DoesHitHaveSimChannel( art::Ptr<recob::Hit> const& hit){
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) 
      if( chan->Channel() == hit->Channel() ) return true;
    return false;
  }
  
  //===================================================================
  float ModBoxRecomb(float dEdx, float Efield){
    float B = 0.153142;
    float A = 0.93;
    float Xi = B * dEdx / Efield;
    return log(A+Xi)/Xi;
  }

  //====================================================================
  // This function calculates the leading MCParticle ID contributing to a hit and the
  // fraction of that hit's energy comes from that particle.
  void  HitTruth(art::Ptr<recob::Hit> const& hit, int& truthid, float& truthidEnergyFrac, float& energy,float& numElectrons){
    // Get associated sim::TrackIDEs for this hit
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    std::vector<sim::TrackIDE> trackIDEs 
      = art::ServiceHandle<cheat::BackTrackerService>()->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return;
    // Loop through and find the leading TrackIDE, and keep
    // track of the total energy of ALL IDEs.
    float maxe = 0;
    float bestfrac = 0;
    float bestid = 0;
    float ne = 0;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      ne += (float)trackIDEs[i].numElectrons;
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
 
  //==================================================================
  // Returns list of all G4 track IDs associated with a hit
  std::set<int> HitTruthIds( art::Ptr<recob::Hit> const& hit){
    std::set<int> ids;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(size_t i = 0; i < trackIDEs.size(); ++i) ids.insert(trackIDEs[i].trackID);
    //std::vector<int> vout(ids.begin(),ids.end());
    //return vout;
    return ids;
  }
  
  
  //=====================================================================
  // Get MCTruth associated with TrackID using a try bracket to avoid
  // fatal exceptions (return false if no match or exception thrown)
  bool G4IdToMCTruth( int const trkID, art::Ptr<simb::MCTruth>& mctruth )
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
  

  
  //=============================================================================
  // Length of reconstructed track, trajectory by trajectory.
  double PathLength(const simb::MCParticle& part, TVector3& start, TVector3& end)
  {
    // Get geo boundaries
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
   
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
        L += (p1-p0).Mag();
        if(first)	start = p1; 
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


  //=============================================================================
  // Calculate distance to boundary.
  double DistToBoundary(const recob::Track::Point_t& pos)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
    double d1 = std::min(pos.X()-xmin,xmax-pos.X());  // dist to wire plane
    double d2 = fabs(pos.X());                        // dist to cathode
    double d3 = std::min(pos.Y()-ymin,ymax-pos.Y());  // dist to closest top/bottom
    double d4 = std::min(pos.Z()-zmin,zmax-pos.Z());  // dist to closest upstream/downstream wall
    return std::min( std::min( std::min(d1,d2), d3), d4);
  }


  //===========================================================================
  void GetGeoBoundaries(double& xmin, double& xmax, double& ymin, double& ymax, double&zmin, double& zmax){
    art::ServiceHandle<geo::Geometry> geom;
    xmin = -2.0 * geom->DetHalfWidth() - 1e-8;
    xmax = 2.0 * geom->DetHalfWidth() + 1e-8;
    ymin = -geom->DetHalfHeight() -1e-8;
    ymax = geom->DetHalfHeight() + 1e-8;
    zmin = 0. -1e-8;
    zmax = geom->DetLength() + 1e-8;
  }

  //===========================================================================
  bool IsPointInAV(float x, float y, float z){
    
    // Get geo boundaries
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
      
    if(     x >= xmin && x <= xmax
        &&  y >= ymin && y <= ymax
        &&  z >= zmin && z <= zmax ) {
      return true;
    } else {
      return false;
    }
    
  }
  
  bool IsPointInAV(TVector3& v){
    return IsPointInAV(v.X(), v.Y(), v.Z());
  }

}
