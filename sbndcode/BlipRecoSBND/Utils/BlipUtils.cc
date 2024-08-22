#include "BlipUtils.h"

namespace BlipUtils {

  //============================================================================
  // Find total visible energy deposited in the LAr, and number of electrons deposited
  // and drifted to the anode.
  /*
  void CalcTotalDep(float& energy, int& ne_dep, float& ne_anode, SEDVec_t& sedvec){
   
    // energy and electrons deposited
    energy = 0; 
    ne_dep = 0;
    for(auto& sed : sedvec ) { 
      energy += sed->Energy(); 
      ne_dep += sed->NumElectrons();
    }
    
    // electrons drifted to collection plane wires
    art::ServiceHandle<geo::Geometry> geom;
    ne_anode = 0;
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) {
      if( geom->View(chan->Channel()) != geo::kW ) continue;
      for(auto const &tdcide : chan->TDCIDEMap() ) {
        for(const auto& ide : tdcide.second) ne_anode += ide.numElectrons;
      }
    }
  }
  */
 

  //===========================================================================
  // Provided a MCParticle, calculate everything we'll need for later calculations
  // and save into ParticleInfo object
  void FillParticleInfo( const simb::MCParticle& part, blip::ParticleInfo& pinfo, SEDVec_t& sedvec, int caloPlane){
    
    // Get important info and do conversions
    pinfo.particle    = part;
    pinfo.trackId     = part.TrackId();
    pinfo.isPrimary   = (int)(part.Process() == "primary");
    pinfo.mass        = /*GeV->MeV*/1e3 * part.Mass();
    pinfo.E           = /*GeV->MeV*/1e3 * part.E();
    pinfo.endE        = /*GeV->MeV*/1e3 * part.EndE();
    pinfo.KE          = /*GeV->MeV*/1e3 * (part.E()-part.Mass());
    pinfo.endKE       = /*GeV->MeV*/1e3 * (part.EndE()-part.Mass());
    pinfo.P           = /*GeV->MeV*/1e3 * part.Momentum().Vect().Mag();
    pinfo.Px          = /*GeV->MeV*/1e3 * part.Px();
    pinfo.Py          = /*GeV->MeV*/1e3 * part.Py();
    pinfo.Pz          = /*GeV->MeV*/1e3 * part.Pz();
    pinfo.time        = /*ns ->mus*/1e-3 * part.T();
    pinfo.endtime     = /*ns ->mus*/1e-3 * part.EndT();
    pinfo.numTrajPts  = part.NumberTrajectoryPoints();

    // Pathlength (in AV) and start/end point
    pinfo.pathLength  = PathLength( part, pinfo.startPoint, pinfo.endPoint);

    // Central position of trajectory
    pinfo.position    = 0.5*(pinfo.startPoint+pinfo.endPoint);

    // Energy/charge deposited by this particle, found using SimEnergyDeposits 
    pinfo.depEnergy     = 0;
    pinfo.depElectrons  = 0;
    for(auto& sed : sedvec ) {
      if( abs(sed->TrackID()) == part.TrackId() ) {
        pinfo.depEnergy     += sed->Energy();
        pinfo.depElectrons  += sed->NumElectrons();
      }
    }
    
    return;
  
  }

  //===================================================================
  // Provided a vector of all particle information for event, fill a
  // vector of true blips
  void MakeTrueBlips( std::vector<blip::ParticleInfo>& pinfo, std::vector<blip::TrueBlip>& trueblips ) {
     
    art::ServiceHandle<geo::Geometry> geom;
    auto const& detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    
    for(size_t i=0; i<pinfo.size(); i++){
      auto& part = pinfo[i].particle;
      
      //std::cout<<"Making true blip for "<<part.TrackId()<<" (PDG "<<part.PdgCode()<<", which deposited "<<pinfo[i].depEnergy<<"\n";

      // If this is a photon or neutron, don't even bother!
      if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) continue;

      // If this is an electron that came from another electron, it would 
      // have already been grouped as part of the contiguous "blip" previously.
      std::string proc = part.Process();
      if( part.PdgCode() == 11 && ( proc == "eIoni" || proc == "muIoni" || proc == "hIoni") ) continue;

      // Create the new blip
      blip::TrueBlip tb;
      GrowTrueBlip(pinfo[i],tb);
      if( !tb.Energy ) continue;  

      // We want to loop through any contiguous electrons (produced
      // with process "eIoni") and add the energy they deposit into this blip.
      if( part.NumberDaughters() ) {
        for(size_t j=0; j<pinfo.size(); j++){
          simb::MCParticle& p = pinfo[j].particle;
          std::string pr = p.Process();
          if( p.PdgCode() != 2112 && (pr == "eIoni" || pr == "muIoni" || pr == "hIoni") ){
            if( IsAncestorOf(p.TrackId(),part.TrackId(),true) ) GrowTrueBlip(pinfo[j],tb);
          }
        }
      }
      
      // Final check -- ensure there was non-negligible number 
      // of deposted ionization electrons
      if( tb.DepElectrons < 20 ) continue;

      // Calculate TPC-specific quantities
      
      // 'ConvertXToTicks' does not account for time offset of particle (i.e., it
      // assumes particle T0 = 0 with the trigger). We need to correct for that.
      //float tick_offset = (tb.Time>0) ? tb.Time/clockData.TPCClock().TickPeriod() : 0;
      auto point = geo::Point_t{tb.Position.X(),tb.Position.Y(),tb.Position.Z()};
      auto const& tpcID   = geom->FindTPCAtPosition(point);
      auto const& planeID = art::ServiceHandle<geo::Geometry>()->GetBeginPlaneID(tpcID);
      float tick_calc = (float)detProp.ConvertXToTicks(tb.Position.X(),planeID);
      tb.DriftTime = tick_calc*clockData.TPCClock().TickPeriod() + clockData.TriggerOffsetTPC();
      
      tb.ID = trueblips.size();
      trueblips.push_back(tb);

    }
    
  }
  
  
  //====================================================================
  void GrowTrueBlip( blip::ParticleInfo& pinfo, blip::TrueBlip& tblip ) {
    
    simb::MCParticle& part = pinfo.particle;

    // Skip neutrons, photons
    if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) return;
    
    // Check that path length isn't zero
    if( !pinfo.pathLength ) return;

    // If this is a new blip, initialize
    if( !tblip.G4ChargeMap.size() ) {
      tblip.Position    = pinfo.position;
      tblip.Time        = pinfo.time;
    
    // .. otherwise, check that the new particle
    // creation time is comparable to existing blip.
    // then calculate new energy-weighted position.
    } else if ( fabs(tblip.Time-pinfo.time) < 3 ) {
      float totE = tblip.Energy + pinfo.depEnergy;
      float w1 = tblip.Energy/totE;
      float w2 = pinfo.depEnergy/totE;
      tblip.Position    = w1*tblip.Position + w2*pinfo.position;
      tblip.Time        = w1*tblip.Time     + w2*pinfo.time;
      tblip.LeadCharge  = pinfo.depElectrons;
    // ... if the particle isn't a match, show's over
    } else {
      return;
    }

    tblip.Energy      += pinfo.depEnergy;
    tblip.DepElectrons+= pinfo.depElectrons;
    tblip.NumElectrons+= std::max(0.,pinfo.numElectrons);
    
    tblip.G4ChargeMap[part.TrackId()] += pinfo.depElectrons;
    tblip.G4PDGMap[part.TrackId()]    = part.PdgCode();
    if(pinfo.depElectrons > tblip.LeadCharge ) {
      tblip.LeadCharge  = pinfo.depElectrons;
      tblip.LeadG4Index = pinfo.index;
      tblip.LeadG4ID    = part.TrackId();
      tblip.LeadG4PDG   = part.PdgCode();
    }
  }

  
  //====================================================================
  // Merge blips that are close
  void MergeTrueBlips(std::vector<blip::TrueBlip>& vtb, float dmin){
    if( dmin <= 0 ) return;
    std::vector<blip::TrueBlip> vtb_merged;
    std::vector<bool> isGrouped(vtb.size(),false);
    
    for(size_t i=0; i<vtb.size(); i++){
      if( isGrouped.at(i) ) continue;
      else isGrouped.at(i) = true;
      auto& blip_i = vtb.at(i);
      for(size_t j=i+1; j<vtb.size(); j++){
        if( isGrouped.at(j) ) continue;
        auto const& blip_j = vtb.at(j);
        if( blip_i.TPC != blip_j.TPC ) continue;
        // check that the times are similar (we don't want to merge
        // together a blip that happened much later but in the same spot)
        if( fabs(blip_i.Time - blip_j.Time) > 5 ) continue;
        float d = (blip_i.Position-blip_j.Position).Mag();
        if( d < dmin ) {
          isGrouped.at(j) = true;
          //float totE = blip_i.Energy + blip_j.Energy;
          float totQ = blip_i.DepElectrons + blip_j.DepElectrons;
          float w1 = blip_i.DepElectrons/totQ;
          float w2 = blip_j.DepElectrons/totQ;
          blip_i.Energy       += blip_j.Energy;
          blip_i.Position     = w1*blip_i.Position + w2*blip_j.Position;
          blip_i.DriftTime    = w1*blip_i.DriftTime+ w2*blip_j.DriftTime; 
          blip_i.Time         = w1*blip_i.Time + w2*blip_j.Time; 
          blip_i.DepElectrons += blip_j.DepElectrons;
          if( blip_j.NumElectrons ) blip_i.NumElectrons += blip_j.NumElectrons;
          
          blip_i.G4ChargeMap.insert(blip_j.G4ChargeMap.begin(), blip_j.G4ChargeMap.end());
          blip_i.G4PDGMap.insert(blip_j.G4PDGMap.begin(), blip_j.G4PDGMap.end());
          
          if( blip_j.LeadCharge > blip_i.LeadCharge ) {
            blip_i.LeadCharge   = blip_j.LeadCharge;
            blip_i.LeadG4ID     = blip_j.LeadG4ID;
            blip_i.LeadG4Index  = blip_j.LeadG4Index;
            blip_i.LeadG4PDG    = blip_j.LeadG4PDG;
          }
        }//d < dmin
      }//loop over blip_j
      blip_i.ID = vtb_merged.size();
      vtb_merged.push_back(blip_i);
    }
    vtb.clear();
    vtb = vtb_merged;
  }

  
  //=================================================================
  blip::HitClust MakeHitClust(std::vector<blip::HitInfo> const& hitinfoVec){
    
    blip::HitClust hc;
    if( hitinfoVec.size() ) {
      int tpc   = hitinfoVec[0].tpc;
      int plane = hitinfoVec[0].plane;
      int cryo  = hitinfoVec[0].cryo;

      // check that all hits are on same plane;
      for(auto& h : hitinfoVec ) {
        if( h.cryo  != cryo   ) return hc;
        if( h.tpc   != tpc    ) return hc;
        if( h.plane != plane  ) return hc;
      }

      // initialize values
      hc.Cryostat         = cryo;
      hc.TPC              = tpc;
      hc.Plane            = plane;
      hc.ADCs             = 0;
      hc.Charge           = 0;
      hc.SigmaCharge      = 0;
      hc.Amplitude        = 0;
      hc.NPulseTrainHits  = 0;
      float startTime     = 9e9;
      float endTime       = -9e9;
      float weightedTick  = 0;
      float weightedTime  = 0;
      float weightedGOF   = 0;
      //float weightedRatio = 0;
      float qGOF          = 0;

      // store hit times, charges, and RMS
      std::vector<float> tvec;
      std::vector<float> qvec;
      std::vector<float> dqvec;
      std::vector<float> rmsvec;

      // grow our hit cluster!
      for(auto& hitinfo : hitinfoVec ) {
        if( hc.HitIDs.find(hitinfo.hitid) != hc.HitIDs.end() ) continue;
        hc.HitIDs     .insert(hitinfo.hitid);
        hc.Wires      .insert(hitinfo.wire);
        hc.Chans      .insert(hitinfo.chan);
        float q       = (hitinfo.charge > 0)? hitinfo.charge : 0;
        float integral = hitinfo.integralADC;
        float sigma   = hitinfo.sigmaintegral;
        float dq      = (integral != 0 && sigma>0)? (sigma/integral)*q : 0;
        hc.Charge     += q;
        hc.ADCs       += hitinfo.integralADC;
        hc.Amplitude  = std::max(hc.Amplitude, hitinfo.amp );
        weightedTick  += q*hitinfo.peakTime;
        weightedTime  += q*hitinfo.driftTime;
        startTime     = std::min(startTime, hitinfo.driftTime-hitinfo.rms);
        endTime       = std::max(endTime,   hitinfo.driftTime+hitinfo.rms);
        tvec          .push_back(hitinfo.driftTime);
        qvec          .push_back(q);
        dqvec         .push_back(dq);
        rmsvec        .push_back(hitinfo.rms);
        if( hitinfo.g4trkid >= 0 ) hc.G4IDs.insert(hitinfo.g4trkid);
        if( hitinfo.gof < 0 ) {
          hc.NPulseTrainHits++;
        } else { 
          weightedGOF   += q*hitinfo.gof; 
          qGOF          += q;
        }
      }//endloop over hits
      
      // mean goodness of fit
      if( qGOF ) hc.GoodnessOfFit = weightedGOF/qGOF;

      // calculate other quantities
      hc.NHits      = hc.HitIDs.size();
      hc.NWires     = hc.Wires.size();
      hc.CenterWire =(*hc.Wires.begin()+*hc.Wires.rbegin())/2.;
      hc.CenterChan =(*hc.Chans.begin()+*hc.Chans.rbegin())/2.;
      hc.StartWire  = *hc.Wires.begin();
      hc.EndWire    = *hc.Wires.rbegin();
      hc.StartTime  = startTime;
      hc.EndTime    = endTime;
      hc.Timespan   = hc.EndTime - hc.StartTime;
      hc.Time       = weightedTime / hc.Charge;
      hc.TimeTick   = weightedTick / hc.Charge;

      // overall cluster RMS and uncertainty in charge
      float sig_sumSq = 0;
      float dt_sumSq  = 0;
      float dq        = 0;
      for(size_t i=0; i<qvec.size(); i++) {
        float w = qvec[i] / hc.Charge;
        dt_sumSq  += w*pow(tvec[i]-hc.Time,2);
        sig_sumSq += pow(w*rmsvec[i],2);
        dq        += w*dqvec[i];

      }
      hc.RMS = sqrt( sig_sumSq + dt_sumSq );
      hc.SigmaCharge = dq;


    }//endif > 0 hits
  
    // mark the cluster as valid and ship it out
    hc.isValid = true;
    return hc; 
  }


  //=================================================================
  blip::Blip MakeBlip( std::vector<blip::HitClust> const& hcs,
    detinfo::DetectorPropertiesData const& detProp,
    detinfo::DetectorClocksData const& clockData ){
    
    blip::Blip  newblip;
    
    // ------------------------------------------------
    // Must be 1-3 clusts (no more, no less!)
    if( hcs.size() > 3  || hcs.size() < 1 ) return newblip;

    // ------------------------------------------------
    // All hits must be in same TPC, and no 2 planes can be repeated
    std::set<int> planeIDs;
    for(size_t i=0; i<hcs.size(); i++) {
      planeIDs.insert(hcs[i].Plane);
      for(size_t j=i+1; j<hcs.size(); j++){
        if( hcs[i].Plane == hcs[j].Plane )  { return newblip; }
        if( hcs[i].TPC   != hcs[j].TPC )    { return newblip; }
      }
    }
      
    // ------------------------------------------------
    // detector properties initialization
    //auto const& detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    //auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    float driftVelocity   = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature()); 
    float tick_to_cm      = clockData.TPCClock().TickPeriod() * driftVelocity;
    

    newblip.Cryostat = hcs[0].Cryostat;
    newblip.TPC     = hcs[0].TPC;
    newblip.NPlanes = planeIDs.size();
    
    // ------------------------------------------------
    /// Look for valid wire intersections between 
    // central-most hits in each cluster
    std::vector<TVector3> wirex;
    for(size_t i=0; i<hcs.size(); i++) {
      int pli = hcs[i].Plane;
      
      // use view with the maximal wire extent to calculate transverse (YZ) length
      if( hcs[i].NWires > newblip.MaxWireSpan ) {
        newblip.MaxWireSpan = hcs[i].NWires;
        newblip.dYZ         = hcs[i].NWires * art::ServiceHandle<geo::Geometry>()->WirePitch(kViews[pli]);
      }
  
      for(size_t j=i+1; j<hcs.size(); j++){
        int plj = hcs[j].Plane;
          
        double y,z;
        bool match3d = false;
        // If this was already calculated, use that
        if( hcs[i].IntersectLocations.count(hcs[j].ID) ) {
          match3d = true;
          y = hcs[i].IntersectLocations.find(hcs[j].ID)->second.Y();
          z = hcs[i].IntersectLocations.find(hcs[j].ID)->second.Z();
        } else {
          match3d = art::ServiceHandle<geo::Geometry>()
            ->ChannelsIntersect(hcs[i].CenterChan,hcs[j].CenterChan,y,z);
        }

        if( match3d ) {
          TVector3 a(0., y, z);
          wirex.push_back(a);
          newblip.clusters[pli] = hcs[i];
          newblip.clusters[plj] = hcs[j];
        }
      }
    }
   
    // Require some number of intersection points.
    if( !wirex.size() ) return newblip;
    
    // Loop over the intersection points and calculate average position in 
    // YZ-plane, as well as the mean difference between intersection points.
    newblip.Position.SetXYZ(0,0,0);
    if( wirex.size() == 1 ) {
      newblip.Position= wirex[0];
    } else {
      newblip.SigmaYZ = 0;
      double fact = 1./wirex.size();
      for(auto& v : wirex ) newblip.Position  += v * fact;
      for(auto& v : wirex ) newblip.SigmaYZ   += (v-newblip.Position).Mag() * fact;
      // Ensure that difference between intersection points is
      // consistent with the maximal wire extent
      if( newblip.SigmaYZ > std::max(1.,0.5*newblip.dYZ) ) return newblip;
    }
    
    // Calculate mean drift time and X-position
    // (note that the 'time' of each of the hit clusters
    // have already been corrected for plane-to-plane offsets)
    newblip.TimeTick= 0;
    newblip.Time = 0;
    newblip.dX = 0;
    float vsize = (float)hcs.size();
    for(auto hc : hcs ) {
      newblip.TimeTick  += hc.TimeTick / vsize;
      newblip.Time      += hc.Time / vsize;
      newblip.dX  = std::max((float)(hc.EndTime-hc.StartTime)*tick_to_cm, newblip.dX);
    }
   
    //auto const& tpcID = geo::TPCID(geo::CryostatID(newblip.Cryostat),newblip.TPC);
    //auto const& planeID = art::ServiceHandle<geo::Geometry>()->GetBeginPlaneID(tpcID);
    //newblip.Position  .SetX(detProp.ConvertTicksToX(newblip.TimeTick, planeID)); 
    
    // convert ticks to X
      auto const& cryostat= art::ServiceHandle<geo::Geometry>()->Cryostat(geo::CryostatID(newblip.Cryostat));
      auto const& tpcgeom = cryostat.TPC(newblip.TPC);
      auto const  xyz     = tpcgeom.Plane(0).GetCenter();
      int         dirx    = DriftDirX(tpcgeom);
      
      newblip.Position.SetX( xyz.X() + dirx * tick_to_cm * newblip.Time );

    // this should ALREADY be accounted for at the hit-processing level in BlipRecoAlg,
    // through the use of GetXTicksOffset...
    //float offset_ticks = clockData.TriggerOffsetTPC() /  clockData.TPCClock().TickPeriod();
    //float driftTicks = newblip.Time + clockData.TriggerOffsetTPC();
    //std::cout<<"blip time "<<newblip.TimeTick<<"  TPC offset "<<clockData.TriggerOffsetTPC()<<"\n";
   
    //std::cout<<"Made new blip with recoX = "<<newblip.Position.X()<<" and drift time "<<newblip.DriftTime<<" us\n";
    // OK, we made it! Flag as "valid" and ship it out.
    newblip.isValid = true;
    return newblip;
    
  }


  //====================================================================
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    if( particleID == ancestorID  )       return true;
    if( particleID < ancestorID   )       return false;
    if( !plist.HasParticle(ancestorID) )  return false;
    while( particleID > ancestorID ) {
      
      simb::MCParticle p = pi_serv->TrackIdToParticle(particleID);
      if( !plist.HasParticle(p.Mother() ) ) { return false; }
    
      simb::MCParticle pM = pi_serv->TrackIdToParticle(p.Mother());
      if      ( pM.TrackId() == ancestorID )                      { return true;  }
      else if ( breakAtPhots == true && pM.PdgCode() == 22 )      { return false; }
      else if ( pM.Process() == "primary" || pM.TrackId() == 1 )  { return false; }
      else if ( pM.Mother() == 0 )                                { return false; }
      else    { particleID = pM.TrackId(); }              
    }

    return false;
  }
  
  //===================================================================
  int DriftDirX(geo::TPCGeo const& tpcgeom) { 
    return ((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 : -1.0);
  }

  //====================================================================
  bool DoHitsOverlap(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2){
    if( hit1->WireID() != hit2->WireID() ) return false;
    float t1 = hit1->PeakTime();
    float t2 = hit2->PeakTime();
    float sig = std::max(hit1->RMS(),hit2->RMS());
    if( fabs(t1-t2) < 1.0*sig ) return true;
    else return false;
  }

  
  //====================================================================
  bool DoHitClustsOverlap(blip::HitClust const& hc1, blip::HitClust const& hc2){
    
    // only match across different wires in same TPC
    if( hc1.TPC != hc2.TPC    ) return false;

    if(     hc1.StartTime <= hc2.EndTime 
        &&  hc2.StartTime <= hc1.EndTime )  return true;
    else return false;
  }
  bool DoHitClustsOverlap(blip::HitClust const& hc1, float t1, float t2 ){
    blip::HitClust hc2;
    hc2.TPC = hc1.TPC;
    hc2.StartTime = t1;
    hc2.EndTime = t2;
    return DoHitClustsOverlap(hc1,hc2);
  }

  //====================================================================
  // Calculates the level of time overlap between two clusters
  float CalcHitClustsOverlap(blip::HitClust const& hc1, blip::HitClust const& hc2){
    return CalcOverlap(hc1.StartTime,hc1.EndTime,hc2.StartTime,hc2.EndTime);
  }

  float CalcOverlap(const float& x1, const float& x2, const float& y1, const float& y2){
    float full_range = std::max(x2,y2) - std::min(x1,y1);
    float sum        = (x2-x1) + (y2-y1);
    float overlap    = std::max(float(0), sum-full_range);
    if( overlap > 0 ) return 2. * overlap / sum;
    else              return -1;
  }

  //====================================================================
  bool DoChannelsIntersect(int ch1, int ch2 ){
    double y,z;
    return art::ServiceHandle<geo::Geometry>()->ChannelsIntersect(ch1,ch2,y,z);
  }
  
  //====================================================================
  bool DoHitClustsMatch(blip::HitClust const& hc1, blip::HitClust const& hc2, float minDiffTicks = 2){
    if( fabs(hc1.Time-hc2.Time) < minDiffTicks ) return true;
    else return false;
  }

  //====================================================================
  // This function calculates the leading MCParticle ID contributing to a hit and the
  // fraction of that hit's energy coming from that particle.
  /*
  void  HitTruth(art::Ptr<recob::Hit> const& hit, int& truthid, float& truthidEnergyFrac, float& energy,float& numElectrons){
    // Get associated sim::TrackIDEs for this hit
    std::vector<sim::TrackIDE> trackIDEs 
      = art::ServiceHandle<cheat::BackTrackerService>()->HitToTrackIDEs(hit);
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
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(hit);
    for(size_t i = 0; i < trackIDEs.size(); ++i) ids.insert(trackIDEs[i].trackID);
    return ids;
  }
  */
  
  
  //=====================================================================
  // Get MCTruth associated with TrackID using a try bracket to avoid
  // fatal exceptions (return false if no match or exception thrown)
  /*
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
  */
  

  
  //=============================================================================
  // Length of particle trajectory
  double PathLength(const simb::MCParticle& part, TVector3& start, TVector3& end)
  {
    int n = part.NumberTrajectoryPoints();
    if( n <= 1 ) return 0.;
    double  L	= 0.;
    bool	  first	= true; 
    for(int i = 1; i < n; ++i) {
      const auto& p1 = part.Position(i).Vect();
      const auto& p0 = part.Position(i-1).Vect();
      geo::Point_t pp1 { p1.X(), p1.Y(), p1.Z() };
      if( IsInActiveVolume(pp1) ) {
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
    art::ServiceHandle<geo::Geometry> geom;
    double x = pos.X();
    double y = pos.Y();
    double z = pos.Z();
    auto const& tpcid = geom->FindTPCAtPosition(geo::Point_t{x,y,z});
    if( tpcid.TPC == geo::TPCID::InvalidID )  return -9;
    auto const& tpc = geom->TPC(tpcid);
    double dx = std::min(x-tpc.MinX(),tpc.MaxX()-x);
    double dy = std::min(y-tpc.MinY(),tpc.MaxY()-y);
    double dz = std::min(z-tpc.MinZ(),tpc.MaxZ()-z);
    return std::min( std::min(dx,dy), dz );
  }

  //===========================================================================
  // Given a line with endpoints L1,L2, return shortest distance betweene the
  // line and point 'P'
  double DistToLine(TVector3& L1, TVector3& L2, TVector3& p){
    
    // general vector formulation:
    // a = point on a line
    // n = unit vector pointing along line
    // --> d = norm[ (p-a) - ((p-a) dot n) * n ]
    // In our case, 'a' = L1
    //TVector3 a      = L1;
    TVector3 n      = (L2-L1).Unit();
    TVector3 b      = (p-L1);
    double  projLen  = b.Dot(n);
    double d = -1;
    /*
    if      ( projLen < 0             ) d = (p-L1).Mag();
    else if ( projLen > (L2-L1).Mag() ) d = (p-L2).Mag();
    else                                d = (b-projLen*n).Mag();
    */
    if( projLen > 0 && projLen < (L2-L1).Mag() ) {
      d = (b-projLen*n).Mag(); 
    } else {
      d = std::min( (p-L1).Mag(), (p-L2).Mag() );
    }
    
    return d;
  }
  
  double DistToLine2D(TVector2& L1, TVector2& L2, TVector2& p){
    TVector3 newL1(L1.X(), L1.Y(), 0);
    TVector3 newL2(L2.X(), L2.Y(), 0);
    TVector3 newp(p.X(), p.Y(), 0);
    return DistToLine(newL1,newL2,newp);
  }


  //===========================================================================
  bool IsInActiveVolume(geo::Point_t const& p){
    geo::TPCGeo const* TPC = art::ServiceHandle<geo::Geometry>()->PositionToTPCptr(p);
    return TPC? TPC->ActiveBoundingBox().ContainsPosition(p): false;
  }

  
  //==========================================================================
  void NormalizeHist(TH1D* h){
    if( h->GetEntries() > 0 ) {
      h->Scale(1./h->Integral());
      h->SetBit(TH1::kIsAverage);
      h->SetOption("HIST");
    }
  }


  float FindMedian(std::vector<float>& vec){
    if( !vec.size() ) return -9;
    size_t n = vec.size() / 2;
    std::nth_element(vec.begin(),vec.begin()+n,vec.end());
    if( n % 2 != 0 ) { // odd number of elements
      return vec[n];
    }else{
      float a = vec[n]; 
      std::nth_element(vec.begin(),vec.begin()+n-1,vec.end());
      return (a + vec[n-1]) / 2.0; 
    }
  }
  
  float FindMean(std::vector<float>& vec){
    float sum = 0;
    for(auto& v : vec ) sum += v;
    return (vec.size()>0) ? sum/vec.size() : 0;
  }


  /*
  //===========================================================================
  float ConvertTicksToX(float peakTime, int plane, int tpc, int cryo) {
    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* detClock = lar::providerFrom<detinfo::DetectorClocksService>();
    double Efield   = detProp->Efield(0);
    double Temp     = detProp->Temperature();
    // The drift velocity "Fudge factor"... need to look into this more!
    //double fudgeFact = 9.832658e-1;
    double driftVel = detProp->DriftVelocity(Efield,Temp)*fudgeFact;
    double drift    = (peakTime - detProp->GetXTicksOffset(plane,tpc,cryo))*detClock->TPCClock().TickPeriod();
    double X        = drift * driftVel;
    return X;
  }
  */

}
