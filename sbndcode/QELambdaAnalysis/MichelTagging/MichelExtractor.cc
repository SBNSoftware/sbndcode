#include "MichelExtractor.hh"


namespace {

  // Get a handle and warn if the product is not available
  template <typename T>
  bool fetchHandle(art::Event const& e, std::string const& label, art::Handle<T>& handle)
  {
    handle = e.getHandle<T>(label);
    if (!handle.isValid()) {
      mf::LogWarning("MichelExtractor")
        << "Product not found (label '" << label << "'), event "
        << e.id().run() << "_" << e.id().subRun() << "_" << e.id().event();
      return false;
    }
    return true;
  }

  // Map lookup that does not insert missing keys (operator[] does)
  template <typename Map, typename Key>
  int lookupOrDefault(Map const& m, Key key, int def = -1)
  {
    auto it = m.find(key);
    return it == m.end() ? def : it->second;
  }

}


sbnd::MichelExtractor::MichelExtractor(fhicl::ParameterSet const& p, MichelClassifier& michelReader): 
    // Same order as in private
    fMichelPset(  p.get<fhicl::ParameterSet>("MichelParams") ),

    fMCTruthLabel(p.get<std::string>("MCTruthLabel")) ,
    fMCParticleLabel(p.get<std::string>("MCParticleLabel")),
    fPOTModuleLabel(p.get<std::string>("POTLabel")),

    fMCT2MCPLabel(p.get<std::string>("MCT2MCPLabel")),
    fTrack2HitLabel(p.get<std::string>("Track2HitLabel")),
    fShower2HitLabel(p.get<std::string>("Shower2HitLabel")),

    fSliceLabel(p.get<std::string>("SliceLabel")) ,
    fSliceModuleLabel(p.get<std::string>("SliceModuleLabel")) ,
    fPFParticleLabel(p.get<std::string>("PFParticleLabel")) ,
    fPFParticleMetadataLabel(p.get<std::string>("PFParticleMetadataLabel")) ,
    fVertexLabel(p.get<std::string>("VertexLabel")) ,
    fTrackLabel(p.get<std::string>("TrackLabel")) ,
    fShowerLabel(p.get<std::string>("ShowerLabel")) ,
    fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")) ,
    fHitLabel(p.get<std::string>("HitLabel")) ,
    fHit2SpacePointLabel(p.get<std::string>("Hit2SpacePointLabel")) ,
    fBlipLabel(p.get<std::string>("BlipLabel")) ,
    fTPCPMTBarycenterFMModuleLabel(p.get<std::string>("TPCPMTBarycenterFMModuleLabel")) ,

    // Reco parameters
    fMaxVtxDist(p.get< float >("MaxVtxDist")) ,
    fDistanceThreshold(p.get< float >("DistanceThreshold")) ,
    fMinMuLength(p.get< float >("MinMuLength")),
    fFV(p.get< float >("FV")),

    fmichelReader(michelReader)
  {
    // Call appropriate consumes<>() for any products to be retrieved by this module.
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  detinfo::DetectorClocksData const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fOpticalClock =  clockData.TriggerTime();
  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  geo::PlaneID PID(0,1,0);
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  fWirePlanePosition = std::abs( wireReadout.Plane(PID).GetCenter().X() );


  fhicl::ParameterSet pt = fMichelPset;
  fmichelParams.GeoPset.MaxResRange = pt.get<float>("MaxResRange");
  fmichelParams.GeoPset.SlidingWindowSize = pt.get<int>("SlidingWindowSize");
  fmichelParams.GeoPset.MaxDEDx = pt.get<float>("MaxDEDx");
  fmichelParams.RadialPset.MaxRadius = pt.get<float>("MaxRadius");
  fmichelParams.RadialPset.AngleCut = pt.get<float>("AngleCut");
  fmichelParams.RadialPset.NSubHits = pt.get<int>("NSubHits");
  fmichelParams.RadialPset.NRadialIntervals = pt.get<int>("NRadialIntervals");
  fmichelParams.AngleMaxRadius = pt.get<float>("AngleMaxRadius");
  fmichelParams.BlipMaxWireSpan = pt.get<float>("BlipMaxWireSpan");
  fmichelParams.BlipMaxDist = pt.get<float>("BlipMaxDist");
  }


void sbnd::MichelExtractor::getEventBlips(const art::Handle<std::vector<blip::Blip>>& blipHandle) {  

    fblips.clear();

    if (!blipHandle.isValid()) {
        mf::LogWarning("MichelExtractor") << "getEventBlips: invalid blip handle";
        return;
    }

    std::vector<art::Ptr<blip::Blip>> blipVec;
    art::fill_ptr_vector(blipVec, blipHandle);

    for (auto const& b : blipVec){

        MuBlip blip;

        blip.ID           = b->ID;
        blip.IsValid      = b->isValid;
        blip.Cryostat     = b->Cryostat;
        blip.TpcID        = b->TPC;
        blip.NPlanes      = b->NPlanes;
        blip.MaxWireSpan  = b->MaxWireSpan;
        blip.TimeTick     = b->TimeTick;  
        blip.Charge       = b->Charge;
        blip.Energy       = b->Energy;
        blip.EnergyESTAR  = b->EnergyESTAR;
        blip.EnergyPSTAR  = b->EnergyPSTAR;
        blip.ProxTrkDist  = b->ProxTrkDist;
        blip.ProxTrkID    = b->ProxTrkID;
        blip.InCylinder   = b->inCylinder;
        blip.X            = b->Position.X();
        blip.Y            = b->Position.Y();
        blip.Z            = b->Position.Z();
        blip.SigmaYZ      = b->SigmaYZ;
        blip.dX           = b->dX;
        blip.dYZ          = b->dYZ;

        fblips.push_back(blip);
    }
}


// === ANALYZE FUNCTIONS ===

std::vector<art::Ptr<recob::PFParticle>> sbnd::MichelExtractor::getSlicePrimPfp(const std::vector<art::Ptr<recob::PFParticle>> slcPfpVec)
{

  std::vector<art::Ptr<recob::PFParticle>> primPfpVec;

  // Loop over pfp in slice to get prim information
  for(const art::Ptr<recob::PFParticle> &pfp : slcPfpVec) {
    // Only look for primary pfps
    if (pfp->IsPrimary()) {
      primPfpVec.push_back(pfp);
    }
  }

  return primPfpVec;

}


// Update pfp hits

std::array<std::vector<MuHit>, 3> sbnd::MichelExtractor::getPfpHits(const art::Ptr<recob::PFParticle>& pfp, const std::vector<art::Ptr<recob::Hit>>& trackHitVec)
{
    std::array<std::vector<MuHit>, 3> pfpHits;

    for (const art::Ptr<recob::Hit>& h : trackHitVec) {
        MuHit hit;

        hit.PlaneID     = h->WireID().planeID().deepestIndex();
        hit.TpcID       = h->WireID().planeID().parentID().deepestIndex();
        hit.ChannelID   = h->Channel();
        hit.PeakTime    = h->PeakTime() / 2;
        hit.Integral    = h->Integral();
        hit.RMS         = h->RMS();
        hit.PfpID       = pfp->Self();
        
        // [C1] -1 if truth maps are empty (data) or the hit is not in them
        hit.TruePdg = lookupOrDefault(fslcHitID2TruePdgMap, h.key());
        hit.TrueID  = lookupOrDefault(fslcHitID2G4IDMap, h.key());

        if (hit.PlaneID < 0 || hit.PlaneID > 2) continue;   // [C2] out-of-range index

        pfpHits[hit.PlaneID].push_back(hit);
    }

    return pfpHits;
}

void sbnd::MichelExtractor::updateHits(const art::Ptr<recob::PFParticle>& pfp, const std::vector<art::Ptr<recob::Hit>>& trackHitVec) 
{
    std::array<std::vector<MuHit>, 3>  pfpHits = getPfpHits(pfp, trackHitVec);

    for (int i=0; i<3; i++) {
        fmuHits[i].insert(fmuHits[i].end(), pfpHits[i].begin(), pfpHits[i].end());
    }
}


// Update pfp calo

bool isNonDecreasing(const std::vector<MuCalo>& pfpCalos)
{
    for (size_t i = 1; i < pfpCalos.size(); ++i) {
        if (pfpCalos[i].ResRange < pfpCalos[i - 1].ResRange) {
            return false;
        }
    }
    return true;
}

void correctCaloOrdering(std::array<std::vector<MuCalo>, 3>& pfpCalos) {

    for (int i=0; i<3; i++) {
        if (!isNonDecreasing(pfpCalos[i])) {
            std::reverse(pfpCalos[i].begin(), pfpCalos[i].end());
        }
    }
}


std::optional<std::array<std::vector<MuCalo>, 3>> sbnd::MichelExtractor::getPfpCalo(
  const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
  const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
  const art::Handle<std::vector<recob::Track>> &trackHandle 
  )
{   

    std::array<std::vector<MuCalo>, 3> pfpCalos;

    // [C3] The original std::map {{0,2},{2,0}} had no entry for plane 1, so operator[]
    // inserted 0 and plane 1 was merged into plane 0. Explicit table: plane 1 -> 1.
    constexpr std::array<int, 3> correctCaloPlane = {2, 1, 0};

    if (caloVec.size() != 3) return std::nullopt;

    // Loop over calo planes
    for (unsigned i_p = 0; i_p < 3; i_p++) {
        art::Ptr<anab::Calorimetry> calo = caloVec.at(i_p);
        if (calo.isNull()) continue;

        const std::vector<size_t> &tps = calo->TpIndices();

        // [C4] The per-point vectors must be at least as long as dEdx()
        const size_t nPoints = calo->dEdx().size();
        if (tps.size() < nPoints || calo->ResidualRange().size() < nPoints || calo->XYZ().size() < nPoints) {
            mf::LogWarning("MichelExtractor") << "getPfpCalo: inconsistent calorimetry vector sizes on plane " << i_p;
            continue;
        }

        for (unsigned i_c = 0; i_c < nPoints; i_c++) {
            
            int caloTpcID = -1;
            int caloChannelID = -1;
            double caloPeakTime = -1;
            int caloTrueID = -1;
            int caloTruePdg = -1;

            for (const art::Ptr<recob::Hit> &h: trackHitVec) {
                if (h.key() == tps[i_c]) {
                    caloTpcID = h->WireID().TPC;
                    caloChannelID = h->Channel();
                    caloPeakTime = h->PeakTime() / 2;
                    caloTrueID = lookupOrDefault(fslcHitID2G4IDMap, h.key());
                    caloTruePdg = lookupOrDefault(fslcHitID2TruePdgMap, h.key());
                }
            }

            MuCalo muCalo;

            muCalo.PlaneID   = correctCaloPlane[i_p];
            muCalo.TpcID     = caloTpcID;
            muCalo.ChannelID = caloChannelID;
            muCalo.TrueID    = caloTrueID;
            muCalo.TruePdg   = caloTruePdg;
            muCalo.PeakTime  = caloPeakTime;
            muCalo.dEdx      = calo->dEdx()[i_c];
            muCalo.ResRange  = calo->ResidualRange()[i_c];
            muCalo.X         = calo->XYZ()[i_c].X();
            muCalo.Y         = calo->XYZ()[i_c].Y();
            muCalo.Z         = calo->XYZ()[i_c].Z();

            pfpCalos[muCalo.PlaneID].push_back(muCalo);
        }
    }

    correctCaloOrdering(pfpCalos);

    return pfpCalos;
}



void sbnd::MichelExtractor::updateCalo(
  const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
  const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
  const art::Handle<std::vector<recob::Track>> &trackHandle   // [H3]
  ) 
{
  auto pfpCalos = getPfpCalo(caloVec, trackHitVec, trackHandle);

  if (!pfpCalos) return;
  if (pfpCalos->size() < 3) return;

  for (int i=0; i<3; i++) {

      if ((*pfpCalos)[i].empty()) continue;

      double rangeOffset = (*pfpCalos)[i].back().ResRange;
      for (auto& calo : fmuCalo[i]) {
          calo.ResRange += rangeOffset;
      }

      fmuCalo[i].insert(fmuCalo[i].begin(), (*pfpCalos)[i].begin(), (*pfpCalos)[i].end());
  }
}


// Concatenate muon pfps

void sbnd::MichelExtractor::concatenateMuonPFPs(
  const art::Event &e,
  const art::Ptr<recob::PFParticle>& initialPfp,
  const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,  
  const art::Handle<std::vector<recob::Track>> &trackHandle      
  )
{
    if (!pfpHandle.isValid() || !trackHandle.isValid()) {
        mf::LogWarning("MichelExtractor") << "concatenateMuonPFPs: invalid PFParticle/Track handle";
        return;
    }

    art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
    art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);
    art::FindManyP<anab::Calorimetry> track2caloAssn(trackHandle, e, fCalorimetryLabel);

    if (!pfp2trackAssn.isValid() || !track2hitAssn.isValid() || !track2caloAssn.isValid()) {
        mf::LogWarning("MichelExtractor") << "concatenateMuonPFPs: PFParticle->Track, Track->Hit or Track->Calorimetry assns not available";
        return;
    }

    art::Ptr<recob::PFParticle> currentPfp = initialPfp;
    std::vector<art::Ptr<recob::PFParticle>> visited = {currentPfp};

    std::cout << " Initial pfp in michel extractor: " << initialPfp.key() << std::endl;

    while (true){

        double minDist = 9999.0;
        art::Ptr<recob::PFParticle> closestPfp;
        art::Ptr<recob::Track> closestTrack;

        // [C5] The current track does not depend on the candidate: get it once, and
        // stop if it is null (it was dereferenced unchecked inside the loop)
        const art::Ptr<recob::Track> currentTrack = pfp2trackAssn.at(currentPfp.key());
        if (currentTrack.isNull()) {
            mf::LogWarning("MichelExtractor") << "\t\t\t Current track is null, stopping concatenation";
            break;
        }

        const double currentEndX = currentTrack->End().X();
        const double currentEndY = currentTrack->End().Y();
        const double currentEndZ = currentTrack->End().Z();

        for (const art::Ptr<recob::PFParticle> &pfp : pfpVec)
        {   
            if (std::find(visited.begin(), visited.end(), pfp) != visited.end()) continue;

            const art::Ptr<recob::Track> track = pfp2trackAssn.at(pfp.key());
            if(!track.isNonnull()) {
              mf::LogWarning("MichelExtractor")  << "\t\t\t Track is null";
              continue;
            }

            double pfpStartX = track->Start().X();
            double pfpStartY = track->Start().Y();
            double pfpStartZ = track->Start().Z();

            double dist = std::hypot(pfpStartX - currentEndX, pfpStartY - currentEndY, pfpStartZ - currentEndZ);

            // std::cout << "\t\t\t Dist: " << dist << " distTh: " << fDistanceThreshold << " minDist: " << minDist << std::endl;

            if (dist > fDistanceThreshold) continue;

            if (dist < minDist) {
                minDist = dist;
                closestPfp = pfp;
                closestTrack = track;
            }
        }

        if (!closestPfp) break;

        // std::cout << "   *** Closest pfp " << closestPfp.key() << std::endl;

        visited.push_back(closestPfp);
        currentPfp = closestPfp;
        std::vector<art::Ptr<recob::Hit>> trackHitVec = track2hitAssn.at(closestTrack.key());

        // Fill current mu track information for the closest muon pfp
        MuTrack currentMuTrack;

        currentMuTrack.ID = closestPfp->Self();
        currentMuTrack.Length = closestTrack->Length();
        currentMuTrack.StartX = closestTrack->Start().X();
        currentMuTrack.StartY = closestTrack->Start().Y();
        currentMuTrack.StartZ = closestTrack->Start().Z();
        currentMuTrack.EndX = closestTrack->End().X();
        currentMuTrack.EndY = closestTrack->End().Y();
        currentMuTrack.EndZ = closestTrack->End().Z();

        const std::vector<art::Ptr<anab::Calorimetry>> caloVec = track2caloAssn.at(closestTrack.key());

        fmuTracks.push_back(currentMuTrack);
        updateHits(closestPfp, trackHitVec);
        updateCalo(caloVec, trackHitVec, trackHandle);

    }
}


// If anything required is missing, fmuTracks is left empty and Extract() stops
void sbnd::MichelExtractor::getConcatenation(
    const art::Event &e,
    const art::Ptr<recob::PFParticle>& muonPfp,
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
    const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,  
    const art::Handle<std::vector<recob::Track>> &trackHandle       
    )
{
    fmuTracks.clear();

    if (!pfpHandle.isValid() || !trackHandle.isValid()) {
        mf::LogWarning("MichelExtractor") << "getConcatenation: invalid PFParticle/Track handle";
        return;
    }
    if (muonPfp.isNull()) {
        mf::LogWarning("MichelExtractor") << "getConcatenation: null muon PFParticle";
        return;
    }

    art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
    art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);
    art::FindManyP<anab::Calorimetry> track2caloAssn(trackHandle, e, fCalorimetryLabel);

    if (!pfp2trackAssn.isValid() || !track2hitAssn.isValid() || !track2caloAssn.isValid()) {
        mf::LogWarning("MichelExtractor") << "getConcatenation: PFParticle->Track, Track->Hit or Track->Calorimetry assns not available";
        return;
    }

    MuTrack initialMuTrack;


    // std::cout << "Starting first concatenation" << std::endl;


    // std::cout << "Getting initial track from muon key: " << muonPfp.key()<< std::endl;
    const art::Ptr<recob::Track> initialTrack = pfp2trackAssn.at(muonPfp.key());
    if (initialTrack.isNull()) {   // [C6] was dereferenced unchecked
        mf::LogWarning("MichelExtractor") << "getConcatenation: muon PFParticle has no track";
        return;
    }

    // // std::cout << "Getting initial track hits" << std::endl;
    const std::vector<art::Ptr<recob::Hit>> initialHitVec = track2hitAssn.at(initialTrack.key());

    initialMuTrack.ID = muonPfp->Self();
    initialMuTrack.Length = initialTrack->Length();
    initialMuTrack.StartX = initialTrack->Start().X();
    initialMuTrack.StartY = initialTrack->Start().Y();
    initialMuTrack.StartZ = initialTrack->Start().Z();
    initialMuTrack.EndX = initialTrack->End().X();
    initialMuTrack.EndY = initialTrack->End().Y();
    initialMuTrack.EndZ = initialTrack->End().Z();

    // // std::cout << "Getting initial track calorimetry" << std::endl;
    const std::vector<art::Ptr<anab::Calorimetry>> caloVec = track2caloAssn.at(initialTrack.key());

    fmuTracks.push_back(initialMuTrack);
    updateHits(muonPfp, initialHitVec);
    updateCalo(caloVec, initialHitVec, trackHandle);


    // std::cout << "Concatenating remaining pfps" << std::endl;
    // Concatenate the rest of the muon 
    concatenateMuonPFPs(e, muonPfp, pfpVec, pfpHandle, trackHandle);
}


int sbnd::MichelExtractor::vertex2DriftTick(double vt, double vx) {
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}

std::vector<int> sbnd::MichelExtractor::getPointInWires(double x, double y, double z) {
  //Get Vertex in wires
  int tpc = 0;
  if(x > 0)  tpc = 1;

  geo::PlaneID UPlane_ID(0,tpc,0);
  geo::PlaneID VPlane_ID(0,tpc,1);
  geo::PlaneID CPlane_ID(0,tpc,2);


  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  geo::PlaneGeo const& plane =  wireReadout.Plane(UPlane_ID);
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  int recoVU = wireID.deepestIndex();
  int channelVU = wireReadout.PlaneWireToChannel(wireID);

  geo::PlaneGeo const& plane1 =  wireReadout.Plane(VPlane_ID);
  try {
    wireID = plane1.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane1.ClosestWireID(e.suggestedWireID());
  }
  int recoVV = wireID.deepestIndex();
  int channelVV = wireReadout.PlaneWireToChannel(wireID);

  geo::PlaneGeo const& plane2 =  wireReadout.Plane(CPlane_ID);
  try {
    wireID = plane2.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane2.ClosestWireID(e.suggestedWireID());
  }
  int recoVC = wireID.deepestIndex();
  int channelVC = wireReadout.PlaneWireToChannel(wireID);

  int recoTimeTick = vertex2DriftTick(0, x);

  return {recoVU, recoVV, recoVC, recoTimeTick, channelVU, channelVV, channelVC};
}


void sbnd::MichelExtractor::setSlcHitsMaps(const art::Event &e, const std::vector<art::Ptr<recob::Hit>>& slcHitVec, const art::Handle<std::vector<simb::MCParticle>>& MCPHandle) {

    fslcHitID2G4IDMap.clear();
    fslcHitID2TruePdgMap.clear();

    if (!MCPHandle.isValid()) return; 

    auto const clockDataEvt = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    std::vector<art::Ptr<simb::MCParticle>> mcpVec;
    art::fill_ptr_vector(mcpVec, MCPHandle);

    // Loop over hits in the slice
    for (auto const& h : slcHitVec) {
        
        // Get trackID of pfp associated with hit
        int trackID = TruthMatchUtils::TrueParticleID(clockDataEvt, h, true);
        if (trackID == defInt) trackID = -1;

        fslcHitID2G4IDMap[h.key()] = trackID;

        if (trackID == -1) {
            fslcHitID2TruePdgMap[h.key()] = -1;
            continue;
        }
        
        // Look for idx of trackID in MCP vector
        auto it = std::find_if(mcpVec.begin(), mcpVec.end(),
            [trackID](const art::Ptr<simb::MCParticle>& mcp) {
                return mcp->TrackId() == trackID;
            });

        int idx = -1;
        if (it != mcpVec.end()) idx = it - mcpVec.begin();

        // Store true PDG code in map
        if (idx != -1) {
            fslcHitID2TruePdgMap[h.key()] = mcpVec[idx]->PdgCode();
        } else {
            fslcHitID2TruePdgMap[h.key()] = -1;
        }
    }
}


void sbnd::MichelExtractor::getUnassociatedHits(
  const art::Event &e, 
  const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
  const std::vector<art::Ptr<recob::Hit>>& slcHitVec, 
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,   
  const art::Handle<std::vector<recob::Track>> &trackHandle)      
{
  if (!pfpHandle.isValid() || !trackHandle.isValid()) {
    mf::LogWarning("MichelExtractor") << "getUnassociatedHits: invalid PFParticle/Track handle";
    return;
  }

  art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
  art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);

  if (!pfp2trackAssn.isValid() || !track2hitAssn.isValid()) {
    mf::LogWarning("MichelExtractor") << "getUnassociatedHits: PFParticle->Track or Track->Hit assns not available";
    return;
  }

  std::unordered_set<art::Ptr<recob::Hit>> associatedHits;

  for (auto const& pfp: pfpVec){
    const art::Ptr<recob::Track> track = pfp2trackAssn.at(pfp.key());

    if (!track) continue;

    auto trackHitVec = track2hitAssn.at(track.key());

    for (auto const& hit: trackHitVec){
      associatedHits.insert(hit);
    }
  }


  for (auto const& h : slcHitVec) {

    if (associatedHits.find(h) == associatedHits.end()) {
        MuHit hit;

        hit.PlaneID     = h->WireID().planeID().deepestIndex();
        hit.TpcID       = h->WireID().planeID().parentID().deepestIndex();
        hit.ChannelID   = h->Channel();
        hit.PeakTime    = h->PeakTime() / 2;
        hit.Integral    = h->Integral();
        hit.RMS         = h->RMS();
        // [C1] with empty truth maps (data) this was 0 instead of -1 (operator[] default),
        // inconsistent with getPfpHits
        hit.TruePdg     = lookupOrDefault(fslcHitID2TruePdgMap, h.key());
        hit.TrueID      = lookupOrDefault(fslcHitID2G4IDMap, h.key());
        hit.PfpID       = -1; // No associated PFP

        if (hit.PlaneID < 0 || hit.PlaneID > 2) continue;   // [C2]

        fmuHits[hit.PlaneID].push_back(hit);
    }
  }

}

bool sbnd::MichelExtractor::isMuonFV(){
  if (fmuTracks.empty()) return false;  
  
  double trackEndX = fmuTracks.back().EndX;
  double trackEndY = fmuTracks.back().EndY;
  double trackEndZ = fmuTracks.back().EndZ;

  // std::cout << "Reco concatenated muon end position: " << trackEndX << ", " << trackEndY << ", " << trackEndZ << std::endl;

  // [C8] && instead of bitwise & (same result here, but & does not short-circuit and is error prone)
  return ((trackEndX > -200 + fFV) && (trackEndX < 200 - fFV) &&
          (trackEndY > -200 + fFV) && (trackEndY < 200 - fFV) &&
          (trackEndZ > 0 + fFV) && (trackEndZ < 500 - fFV));
}


// Main analyze function
void sbnd::MichelExtractor::Extract(
    art::Event const& e, 
    const int SliceID, 
    const art::Ptr<recob::PFParticle>& MuonPfp,
    geo::Point_t Vertex)
{

  // ==============================================================
  // === DEFAULT OUTPUT (what the caller gets if we return early)
  // ==============================================================

  fmichelFeatures.nMuonTracks = -1;
  fmichelFeatures.lowerChi2Ratio = -1;
  fmichelFeatures.angleQ80 = -1;
  fmichelFeatures.dropRadius = -1;
  fmichelFeatures.nBlips = -1;
  fmichelFeatures.score = -999;
  fmichelFeatures.muonInFV = false;   
  fmichelFeatures.featuresComputed = false;

  // [C10] Start from a clean state in case the object is reused
  for (auto& v : fmuHits) v.clear();
  for (auto& v : fmuCalo) v.clear();
  fmuTracks.clear();

  // ==============================================================
  // === LOAD HANDLES 
  // ==============================================================

  // MC (optional: not available in data)
  art::Handle<std::vector<simb::MCParticle>> MCPHandle = e.getHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);

  // Pandora and blip handles
  art::Handle<std::vector<recob::Slice>>      sliceHandle;
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  art::Handle<std::vector<recob::Track>>      trackHandle;
  art::Handle<std::vector<blip::Blip>>        blipHandle;

  // Evaluate all of them 
  bool handlesOk = true;
  handlesOk &= fetchHandle(e, fSliceLabel,      sliceHandle);
  handlesOk &= fetchHandle(e, fPFParticleLabel, pfpHandle);
  handlesOk &= fetchHandle(e, fTrackLabel,      trackHandle);
  handlesOk &= fetchHandle(e, fBlipLabel,       blipHandle);

  if (!handlesOk) {
     mf::LogWarning("MichelExtractor") << "  Missing required products, Michel features not computed";
    return;
  }

  // Get associations
  art::FindManyP<recob::PFParticle> slice2pfpAssn(sliceHandle, e, fPFParticleLabel);
  art::FindManyP<recob::Hit> slice2hitAssn(sliceHandle, e, fSliceLabel);

  if (!slice2pfpAssn.isValid() || !slice2hitAssn.isValid()) {
    mf::LogWarning("MichelExtractor") << "Slice->PFParticle or Slice->Hit assns not available, Michel features not computed";
    return;
  }

  // ==============================================================
  // === START EVENT ANALYSIS
  // ==============================================================

  // Get event blips
  getEventBlips(blipHandle);
  
  
  // Get Slice vector
  std::vector<art::Ptr<recob::Slice>> slcVec;
  art::fill_ptr_vector(slcVec, sliceHandle);

  // Loop over slices
  for (const art::Ptr<recob::Slice> &slc: slcVec) {

    if (SliceID != slc->ID()) continue; 


    // Get PFPs and hits
    std::vector<art::Ptr<recob::PFParticle>> slcPfpVec = slice2pfpAssn.at(slc.key());
    std::vector<art::Ptr<recob::Hit>> slcHitVec        = slice2hitAssn.at(slc.key());



    if(slcHitVec.empty()) return;

    // Set hits map
    setSlcHitsMaps(e, slcHitVec, MCPHandle);   // returns immediately if MCPHandle is not valid

    // Set slice vertex from FRANS
    fslcVtxX = Vertex.X();
    fslcVtxY = Vertex.Y();
    fslcVtxZ = Vertex.Z();


    
    // Concatenate pfp
    if (slcPfpVec.empty()) {
        fmichelFeatures.score = -4;
        mf::LogWarning("MichelExtractor") << "  Slice " << SliceID << " has no pfps. Skipping slice.";
        return; 
    }

    // std::cout << "  Concatenating pfps..." << std::endl;
    getConcatenation(e, MuonPfp, slcPfpVec, pfpHandle, trackHandle);

    if (fmuTracks.empty()) { 
        mf::LogWarning("MichelExtractor") << "  Muon concatenation failed. Skipping slice.";
        fmichelFeatures.score = -4;
        return;
    }

    // Get unassociated hits
    // std::cout << "  Getting unassociated hits..." << std::endl;
    getUnassociatedHits(e, slcPfpVec, slcHitVec, pfpHandle, trackHandle);

    // Check muon is in fv
    if (!isMuonFV()) {
        // std::cout << "  Slice primary muon is not in fiducial volume. Skipping slice." << std::endl;
        fmichelFeatures.muonInFV = false;
        fmichelFeatures.score = 2;
        return;
    }
    fmichelFeatures.muonInFV = true;


    // ------- Compute michel vars ------

    std::array<double, 3> recoVtx = {fslcVtxX, fslcVtxY, fslcVtxZ};

    // std::cout << "  Running michel driver..." << std::endl;
    // std::cout << "  mu tracks size: " << fmuTracks.size() << std::endl;
    MichelDriver michelDriver(fmichelParams, fmuTracks, fmuHits, fmuCalo, fblips, recoVtx, fmichelReader);
    michelDriver.Compute();
    fmichelFeatures = michelDriver.Features();

    return;
  }

  // Only reached if no slice matched (the found-slice path always returns)
  mf::LogWarning("MichelExtractor") << "Slice " << SliceID << " not found in the slice collection";
}