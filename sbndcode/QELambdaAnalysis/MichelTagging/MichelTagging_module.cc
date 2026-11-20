#include "MichelTagging_module.hh"


sbnd::MichelTagging::MichelTagging(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
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
    fmichelModelConfig(p.get<std::string> ("MichelBDTConfig"))
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

    // =========================
    //       MICHEL BDT
    // =========================
    // Find BDT weights files
    std::string michelModelPath;
    cet::search_path sp("FW_SEARCH_PATH");
    if ( !sp.find_file(fmichelModelConfig, michelModelPath) )
      throw cet::exception("FRANS producer") << "Could not locate michel BDT configuration file: " <<
          fmichelModelConfig << " not found in FW_SEARCH_PATH\n";
    std::cout << " Michel BDT weights file name: " << michelModelPath << std::endl;

    fmichelBDT.Load(michelModelPath);
  }


void sbnd::MichelTagging::respondToOpenInputFile(art::FileBlock const &fb) {
  ffileName = fb.fileName();
}

// === RUN START-END FUNCTIONS ===

void sbnd::MichelTagging::beginJob()
{

  art::ServiceHandle <art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Add branches to the TTree
  fTree->Branch("file_name", &ffileName);
  fTree->Branch("event_id", &feventID);
  fTree->Branch("run_id", &frunID);
  fTree->Branch("subrun_id", &fsubrunID);

  // Nu branches
  fTree->Branch("gen_index", &fgenIndex);
  fTree->Branch("nu_pdg", &fnuPDG);
  fTree->Branch("nu_e0", &fnuE0);
  fTree->Branch("nu_int_mode", &fnuIntMode);
  fTree->Branch("nu_int_type", &fnuIntType);
  fTree->Branch("nu_cc_nc", &fnuCCNC);

  // Slice variables
  fTree->Branch("slc_id", &fslcID);
  fTree->Branch("slc_vtx_x", &fslcVtxX);
  fTree->Branch("slc_vtx_y", &fslcVtxY);
  fTree->Branch("slc_vtx_z", &fslcVtxZ);
  fTree->Branch("slc_vtx_wire_u", &fslcVtxWireU);
  fTree->Branch("slc_vtx_wire_v", &fslcVtxWireV);
  fTree->Branch("slc_vtx_wire_c", &fslcVtxWireC);
  fTree->Branch("slc_vtx_channel_u", &fslcVtxChannelU);
  fTree->Branch("slc_vtx_channel_v", &fslcVtxChannelV);
  fTree->Branch("slc_vtx_channel_c", &fslcVtxChannelC);
  fTree->Branch("slc_vtx_time_tick", &fslcVtxTimeTick);
  
  // Slice truth matching
  fTree->Branch("slc_true_vtx_x", &fslcTrueVtxX);
  fTree->Branch("slc_true_vtx_y", &fslcTrueVtxY);
  fTree->Branch("slc_true_vtx_z", &fslcTrueVtxZ);
  fTree->Branch("slc_true_vtx_wire_u", &fslcTrueVtxWireU);
  fTree->Branch("slc_true_vtx_wire_v", &fslcTrueVtxWireV);
  fTree->Branch("slc_true_vtx_wire_c", &fslcTrueVtxWireC);
  fTree->Branch("slc_true_vtx_channel_u", &fslcTrueVtxChannelU);
  fTree->Branch("slc_true_vtx_channel_v", &fslcTrueVtxChannelV);
  fTree->Branch("slc_true_vtx_channel_c", &fslcTrueVtxChannelC);
  fTree->Branch("slc_true_vtx_time_tick", &fslcTrueVtxTimeTick);

  // Muon info
  fTree->Branch("mu_found", &fmuFound);
  fTree->Branch("mu_decay", &fmuDecay);
  fTree->Branch("prim_mu_id", &fprimMuID);
  fTree->Branch("prim_mu_len", &fprimMuLength);
  fTree->Branch("n_mu_tracks", &fnMuTracks);

  // Bragg info
  fTree->Branch("bragg_found", &fbraggFound);
  fTree->Branch("bragg_idx", &fbraggIdx);
  fTree->Branch("bragg_channel_id", &fbraggChannelID);
  fTree->Branch("bragg_peak_time", &fbraggPeakTime);
  fTree->Branch("bragg_tpc_id", &fbraggTpcID);
  fTree->Branch("bragg_dedx", &fbraggdEdx);
  fTree->Branch("bragg_res_range", &fbraggResRange);
  fTree->Branch("bragg_x", &fbraggX);
  fTree->Branch("bragg_y", &fbraggY);
  fTree->Branch("bragg_z", &fbraggZ);

  // PCA info
  fTree->Branch("pca_valid", &fpcaValid);
  fTree->Branch("pca_dir_x", &fpcaDirX);
  fTree->Branch("pca_dir_y", &fpcaDirY);

  // Hits info
  // fTree->Branch("hits_plane_id", &fhitPlaneID);
  // fTree->Branch("hits_tpc_id", &fhitTpcID);
  // fTree->Branch("hits_channel_id", &fhitChannelID);
  // fTree->Branch("hits_peak_time", &fhitPeakTime);
  // fTree->Branch("hits_integral", &fhitIntegral);
  // fTree->Branch("hits_rms", &fhitRMS);
  // fTree->Branch("hits_true_pdg", &fhitTruePdg);
  // fTree->Branch("hits_pfp_id", &fhitPfpID);
  // fTree->Branch("hits_true_id", &fhitTrueID);

  // // Calo info
  // fTree->Branch("calo_plane_id", &fcaloPlaneID);
  // fTree->Branch("calo_tpc_id", &fcaloTpcID);
  // fTree->Branch("calo_channel_id", &fcaloChannelID);
  // fTree->Branch("calo_peak_time", &fcaloPeakTime); 
  // fTree->Branch("calo_true_id", &fcaloTrueID);
  // fTree->Branch("calo_true_pdg", &fcaloTruePdg);
  // fTree->Branch("calo_dedx", &fcalodEdx);
  // fTree->Branch("calo_res_range", &fcaloResRange);
  // fTree->Branch("calo_x", &fcaloX); 
  // fTree->Branch("calo_y", &fcaloY); 
  // fTree->Branch("calo_z", &fcaloZ);

  // Angle info
  fTree->Branch("angles", &fangles);

  // Radial info
  fTree->Branch("radial_r", &fradii);
  fTree->Branch("radial_charge", &fcharge);
  fTree->Branch("radial_dqdr", &fdQdr);
  fTree->Branch("radial_smooth_dqdr", &fsmoothdQdr);
  fTree->Branch("radial_drop_idx", &fdropIdx);
  fTree->Branch("radial_min_dqdr", &fmindQdr);

  // Michel vars
  fTree->Branch("features_computed", &fmichelFeatures.featuresComputed);
  fTree->Branch("lower_chi2_ratio", &fmichelFeatures.lowerChi2Ratio);
  fTree->Branch("all_chi2_ratio", &fmichelFeatures.allChi2Ratio);
  fTree->Branch("angle_q80", &fmichelFeatures.angleQ80);
  fTree->Branch("drop_radius", &fmichelFeatures.dropRadius);
  fTree->Branch("n_blips", &fmichelFeatures.nBlips);
  fTree->Branch("michel_score", &fmichelFeatures.score);

}

void sbnd::MichelTagging::endJob()
{
  // Implementation of optional member function here.
}



// === SUBRUN START-END FUNCTIONS ===


void sbnd::MichelTagging::beginSubRun(const art::SubRun &sr) {

}

//Code ran at the end of each subrun
void sbnd::MichelTagging::endSubRun(const art::SubRun &sr) {

}

void sbnd::MichelTagging::getEventBlips(art::ValidHandle<std::vector<blip::Blip>> blipHandle) {
    std::vector<art::Ptr<blip::Blip>> blipVec;
    art::fill_ptr_vector(blipVec, blipHandle);

    fblips.clear();

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

std::vector<art::Ptr<recob::PFParticle>> sbnd::MichelTagging::getSlicePrimPfp(const std::vector<art::Ptr<recob::PFParticle>> slcPfpVec)
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


void sbnd::MichelTagging::getSliceTrueInfo(
  const art::Event &e, 
  const std::vector<art::Ptr<recob::Hit>>& slcHitVec, 
  const art::ValidHandle<std::vector<simb::MCTruth>> MCTHandle
  )
{
    art::FindManyP<simb::MCParticle> mct2mcpAssn(MCTHandle, e, fMCParticleLabel);

    // Slice mct hits map
    std::map<art::Ptr<simb::MCTruth>, int> slcMCTHitMap;
    int maxHits = -1;
    art::Ptr<simb::MCTruth> bestMCT;

    // Loop over hits in the slice
    for (auto const& hit : slcHitVec) {
        
        // Get trackID of pfp associated with hit
        int trackID = fslcHitID2G4IDMap[hit.key()];
        if (trackID == -1) continue;

        // Get mct associated with trackID
        art::Ptr<simb::MCTruth> mct;
        try {
            mct = fpiServ->TrackIdToMCTruth_P(trackID);
        }
        catch (const cet::exception& ex) {
            mf::LogWarning("MichelTagging") << "TrackIdToMCTruth_P failed: " << ex.what();
            continue;
        }

        if (!mct) continue; 

        // Add hits to slcMCTHitMap
        int& nhits = slcMCTHitMap[mct];
        nhits++;

        // Update best MCT
        if (nhits > maxHits) {
            maxHits = nhits;
            bestMCT = mct;
        }
    }

    
    if (maxHits>0){
        // Get Gen info of best MCT
        std::cout << "Getting Gen info for best MCT..." << std::endl;
        getGenInfo(bestMCT);

        // Get G4 vector associated with best MCT
        std::cout << "Getting G4 info for best MCT..." << std::endl;
        std::vector<art::Ptr<simb::MCParticle>> mcpVec = mct2mcpAssn.at(bestMCT.key());
        getG4Info(mcpVec);
    }
}


void sbnd::MichelTagging::getGenInfo(const art::Ptr<simb::MCTruth> mct)
{
  if(mct.isNonnull() && (mct->Origin() != 0)) {
    if(mct->Origin() == 1) {
      fgenIndex = 1;

      fnuPDG = mct->GetNeutrino().Nu().PdgCode();
      fnuE0 = mct->GetNeutrino().Nu().E(0);

      fnuIntMode = mct->GetNeutrino().Mode();
      fnuIntType = mct->GetNeutrino().InteractionType();
      fnuCCNC = mct->GetNeutrino().CCNC();

      fslcTrueVtxX = mct->GetParticle(0).Vx(0);
      fslcTrueVtxY = mct->GetParticle(0).Vy(0);
      fslcTrueVtxZ = mct->GetParticle(0).Vz(0);

      geo::Point_t slcTrueVertex = {fslcTrueVtxX, fslcTrueVtxY, fslcTrueVtxZ};

      if (fGeom->HasTPC(fGeom->FindTPCAtPosition(slcTrueVertex))) {
        std::vector<int> TrueStartWire = getPointInWires(fslcTrueVtxX, fslcTrueVtxY, fslcTrueVtxZ);
        fslcTrueVtxWireU = TrueStartWire[0];
        fslcTrueVtxWireV = TrueStartWire[1];
        fslcTrueVtxWireC = TrueStartWire[2];
        fslcTrueVtxTimeTick = TrueStartWire[3];
        fslcTrueVtxChannelU = TrueStartWire[4];
        fslcTrueVtxChannelV = TrueStartWire[5];
        fslcTrueVtxChannelC = TrueStartWire[6];
      }

    } else if(mct->Origin() == 2) {
      fgenIndex = 2;
    } else {
      fgenIndex = -1;
    }
  }

}


void sbnd::MichelTagging::getG4Info(const std::vector<art::Ptr<simb::MCParticle>> mcpVec)
{
    // Only look for muon decay or captured
    fmuFound = false;
    fmuDecay = false;
    fprimMuID = -1;

    for (const art::Ptr<simb::MCParticle>& mcp : mcpVec) {

        // Find primary muon
        if (abs(mcp->PdgCode()) != 13) continue;
        if (mcp->Process() != "primary") continue; 
        fmuFound = true;
        fprimMuID = mcp->TrackId();

        // Get last index of trajectory points
        const size_t numberTrajectoryPoints = mcp->NumberTrajectoryPoints();
        const int last = numberTrajectoryPoints - 1;

        // Compute track length of the muon
        double trackLength = 0;
        for(int i_t = 1;i_t <= last; i_t++) {
            TVector3 last_position(mcp->Vx(i_t - 1), mcp->Vy(i_t - 1), mcp->Vz(i_t - 1));
            TVector3 current_position(mcp->Vx(i_t), mcp->Vy(i_t), mcp->Vz(i_t));
            trackLength += (last_position - current_position).Mag();
        }
        fprimMuLength = trackLength;

        // Loop for daughter particles of the muon
        for (const art::Ptr<simb::MCParticle>& muDau : mcpVec) {
            
            if (muDau->Mother() == fprimMuID) {
                if (abs(muDau->PdgCode()) == 12) {
                    fmuDecay = true;
                    break;
                }
            }
        }

        break;
    }

    return;
}


// Update pfp hits

std::array<std::vector<MuHit>, 3> sbnd::MichelTagging::getPfpHits(const art::Ptr<recob::PFParticle>& pfp, const std::vector<art::Ptr<recob::Hit>>& trackHitVec)
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
        hit.TruePdg     = fslcHitID2TruePdgMap[h.key()];
        hit.TrueID      = fslcHitID2G4IDMap[h.key()];
        hit.PfpID       = pfp->Self();

        pfpHits[hit.PlaneID].push_back(hit);
    }

    return pfpHits;
}

void sbnd::MichelTagging::updateHits(const art::Ptr<recob::PFParticle>& pfp, const std::vector<art::Ptr<recob::Hit>>& trackHitVec) 
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


std::optional<std::array<std::vector<MuCalo>, 3>> sbnd::MichelTagging::getPfpCalo(
  const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
  const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
  const art::ValidHandle<std::vector<recob::Track>> &trackHandle
  )
{   

    std::array<std::vector<MuCalo>, 3> pfpCalos;

    std::map<int,int> correctCaloPlane = {{0,2}, {2,0}};

    if (caloVec.size() != 3) return std::nullopt;

    // Loop over calo planes
    for (unsigned i_p = 0; i_p < 3; i_p++) {
        art::Ptr<anab::Calorimetry> calo = caloVec.at(i_p);
        const std::vector<size_t> &tps = calo->TpIndices();

        for (unsigned i_c = 0; i_c < calo->dEdx().size(); i_c++) {
            
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
                    caloTrueID = fslcHitID2G4IDMap[h.key()];
                    caloTruePdg = fslcHitID2TruePdgMap[h.key()];
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



void sbnd::MichelTagging::updateCalo(
  const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
  const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
  const art::ValidHandle<std::vector<recob::Track>> &trackHandle
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

void sbnd::MichelTagging::concatenateMuonPFPs(
  const art::Event &e,
  const art::Ptr<recob::PFParticle>& initialPfp,
  const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
  const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
  const art::ValidHandle<std::vector<recob::Track>> &trackHandle,
  const art::ValidHandle<std::vector<recob::Hit>> &hitHandle
  )
{
    art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
    art::FindOneP<recob::Shower> pfp2showerAssn(pfpHandle, e, fShowerLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn(pfpHandle, e, fPFParticleMetadataLabel);

    art::FindManyP<recob::Vertex> pfp2vertexAssn(pfpHandle, e, fVertexLabel);
    art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);
    art::FindManyP<anab::Calorimetry> track2caloAssn(trackHandle, e, fCalorimetryLabel);

    art::Ptr<recob::PFParticle> currentPfp = initialPfp;
    std::vector<art::Ptr<recob::PFParticle>> visited = {currentPfp};

    while (true){

        double minDist = 9999.0;
        art::Ptr<recob::PFParticle> closestPfp;
        art::Ptr<recob::Track> closestTrack;

        for (const art::Ptr<recob::PFParticle> &pfp : pfpVec)
        {   
            if (std::find(visited.begin(), visited.end(), pfp) != visited.end()) continue;

            const art::Ptr<recob::Track> currentTrack = pfp2trackAssn.at(currentPfp.key());
            const art::Ptr<recob::Track> track = pfp2trackAssn.at(pfp.key());
            if(!track.isNonnull()) {
              std::cout << "\t\t\t Track is null" << std::endl;
              continue;
            }

            double currentEndX = currentTrack->End().X();
            double currentEndY = currentTrack->End().Y();
            double currentEndZ = currentTrack->End().Z();

            double pfpStartX = track->Start().X();
            double pfpStartY = track->Start().Y();
            double pfpStartZ = track->Start().Z();

            double dist = std::hypot(pfpStartX - currentEndX, pfpStartY - currentEndY, pfpStartZ - currentEndZ);

            std::cout << "\t\t\t Dist: " << dist << " distTh: " << fDistanceThreshold << "minDist: " << minDist << std::endl;

            if (dist > fDistanceThreshold) continue;

            if (dist < minDist) {
                minDist = dist;
                closestPfp = pfp;
                closestTrack = track;
            }
        }

        if (!closestPfp) break;
        std::cout << "   *** Closest pfp " << closestPfp.key() << std::endl;

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


int sbnd::MichelTagging::getPfpG4ID(
    const art::Event &e,
    const std::vector<art::Ptr<recob::Hit> >& hitVec 
    )
{
    int pfpG4ID = -1;
    std::map<int, int> g4IDCountMap;

    for (const art::Ptr<recob::Hit>& hit : hitVec) {
        int trackID = fslcHitID2G4IDMap[hit.key()];
        if (trackID == -1) continue;

        g4IDCountMap[trackID]++;
    }

    int maxCount = 0;
    for (const auto& pair : g4IDCountMap) {
        if (pair.second > maxCount) {
            maxCount = pair.second;
            pfpG4ID = pair.first;
        }
    }

    return pfpG4ID;
}


void sbnd::MichelTagging::analyzePFPs(
    const art::Event &e,
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
    const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle,
    const art::ValidHandle<std::vector<recob::Hit>> &hitHandle
    )
{
    art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn(pfpHandle, e, fPFParticleMetadataLabel);

    art::FindManyP<recob::Vertex> pfp2vertexAssn(pfpHandle, e, fVertexLabel);
    art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);
    art::FindManyP<anab::Calorimetry> track2caloAssn(trackHandle, e, fCalorimetryLabel);
    
    fmuTracks.clear();

    MuTrack initialMuTrack;

    double minDist = 9999.0;
    art::Ptr<recob::PFParticle> initialPfp;
    std::vector<art::Ptr<recob::Hit>> initialHitVec = {};

    //Loop over all particles in slice to look for closer pfp to reco vertex
    for (const art::Ptr<recob::PFParticle> &pfp : pfpVec)
    {   

        const art::Ptr<recob::Track> track = pfp2trackAssn.at(pfp.key());
        if(!track.isNonnull()) continue;

        std::vector<art::Ptr<recob::Hit>> trackHitVec = track2hitAssn.at(track.key());
        int pfpG4ID = getPfpG4ID(e, trackHitVec);

        if (pfpG4ID != fprimMuID) continue;
        
        double pfpStartX = track->Start().X();
        double pfpStartY = track->Start().Y();
        double pfpStartZ = track->Start().Z();

        //double pfpEndX = track->End().X();
        //double pfpEndY = track->End().Y();
        //double pfpEndZ = track->End().Z();

        double distVtx = std::hypot(pfpStartX - fslcVtxX, pfpStartY - fslcVtxY, pfpStartZ - fslcVtxZ);

        if (distVtx > fMaxVtxDist) continue;

        if (distVtx < minDist) {
            minDist = distVtx;
            initialPfp = pfp;
            initialHitVec = trackHitVec;
        }
    }

    if (!initialPfp) return;

    std::cout << "*** Initial muon pfp: " << initialPfp.key() << std::endl;


    // Fill mu track information for the initial muon pfp
    const art::Ptr<recob::Track> initialTrack = pfp2trackAssn.at(initialPfp.key());

    initialMuTrack.ID = initialPfp->Self();
    initialMuTrack.Length = initialTrack->Length();
    initialMuTrack.StartX = initialTrack->Start().X();
    initialMuTrack.StartY = initialTrack->Start().Y();
    initialMuTrack.StartZ = initialTrack->Start().Z();
    initialMuTrack.EndX = initialTrack->End().X();
    initialMuTrack.EndY = initialTrack->End().Y();
    initialMuTrack.EndZ = initialTrack->End().Z();

    const std::vector<art::Ptr<anab::Calorimetry>> caloVec = track2caloAssn.at(initialTrack.key());

    fmuTracks.push_back(initialMuTrack);
    updateHits(initialPfp, initialHitVec);
    updateCalo(caloVec, initialHitVec, trackHandle);

    // Concatenate the rest of the muon 
    concatenateMuonPFPs(e, initialPfp, pfpVec, pfpHandle, trackHandle, hitHandle);
}


int sbnd::MichelTagging::vertex2DriftTick(double vt, double vx) {
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}

std::vector<int> sbnd::MichelTagging::getPointInWires(double x, double y, double z) {
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


void sbnd::MichelTagging::clearVars()
{
  fmuTracks.clear();
  fmuHits[0].clear();
  fmuHits[1].clear();
  fmuHits[2].clear();
  fmuCalo[0].clear();
  fmuCalo[1].clear();
  fmuCalo[2].clear();

  fnMuTracks = -1;

  // Bragg
  fbraggFound = false;
  fbraggIdx = -1;
  fbraggChannelID = -999;
  fbraggPeakTime = -999;
  fbraggTpcID = -999;
  fbraggdEdx = -999;
  fbraggResRange = -999;
  fbraggX = -999;
  fbraggY = -999;
  fbraggZ = -999;

  // PCA info
  fpcaValid = false;
  fpcaDirX = -999;
  fpcaDirY = -999;

  // Hits info
  fhitPlaneID.clear();
  fhitTpcID.clear();
  fhitChannelID.clear();
  fhitPeakTime.clear();
  fhitIntegral.clear();
  fhitRMS.clear();
  fhitTruePdg.clear();
  fhitPfpID.clear();
  fhitTrueID.clear();

  // Calo info
  fcaloPlaneID.clear();
  fcaloTpcID.clear();
  fcaloChannelID.clear();
  fcaloPeakTime.clear();
  fcaloTrueID.clear();
  fcaloTruePdg.clear();
  fcalodEdx.clear();
  fcaloResRange.clear();
  fcaloX.clear();
  fcaloY.clear();
  fcaloZ.clear();

  // Angle vars
  fangles.clear();

  // Radial vars
  fradii.clear();
  fcharge.clear();
  fdQdr.clear();
  fsmoothdQdr.clear();
  fdropIdx = -1;
  fmindQdr = -999;
}


void sbnd::MichelTagging::setSlcHitsMaps(const art::Event &e, const std::vector<art::Ptr<recob::Hit>>& slcHitVec, const art::ValidHandle<std::vector<simb::MCParticle>>& MCPHandle) {

    auto const clockDataEvt = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    std::vector<art::Ptr<simb::MCParticle>> mcpVec;
    art::fill_ptr_vector(mcpVec, MCPHandle);

    fslcHitID2G4IDMap.clear();
    fslcHitID2TruePdgMap.clear();

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

void sbnd::MichelTagging::getSlcVertex(const art::Event &e, const std::vector<art::Ptr<recob::PFParticle>>& slcPfpVec, const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle) {
    
  std::vector<art::Ptr<recob::PFParticle>> primVec = getSlicePrimPfp(slcPfpVec);

  // If there is not exactly one primary particle, continue loop
  if (primVec.size() != 1) return;
  auto prim = primVec[0];

  fslcPrimID = prim->Self();
  fslcPrimPDG = prim->PdgCode();

  art::FindManyP<recob::Vertex> pfp2vertexAssn(pfpHandle, e, fVertexLabel);

  // Get slice vertex
  std::vector<art::Ptr < recob::Vertex>> pfpVertexVec = pfp2vertexAssn.at(prim.key());
  if (pfpVertexVec.size() != 1) {
      mf::LogWarning("MichelTagging") << "Slice " << fslcID << " primary PFP has " << pfpVertexVec.size() << " associated vertices." << std::endl;
      return;
  }

  for (const art::Ptr <recob::Vertex> &vertex: pfpVertexVec) {
      fslcVtxX = vertex->position().X();
      fslcVtxY = vertex->position().Y();
      fslcVtxZ = vertex->position().Z();

      if (fGeom->HasTPC(fGeom->FindTPCAtPosition(vertex->position()))) {
      std::vector<int> recoVtxWire = getPointInWires(fslcVtxX, fslcVtxY, fslcVtxZ);
      fslcVtxWireU = recoVtxWire[0];
      fslcVtxWireV = recoVtxWire[1];
      fslcVtxWireC = recoVtxWire[2];
      fslcVtxTimeTick = recoVtxWire[3];
      fslcVtxChannelU = recoVtxWire[4];
      fslcVtxChannelV = recoVtxWire[5];
      fslcVtxChannelC = recoVtxWire[6];
      }
  }
}


void sbnd::MichelTagging::getUnassociatedHits(
  const art::Event &e, 
  const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
  const std::vector<art::Ptr<recob::Hit>>& slcHitVec, 
  const art::ValidHandle<std::vector<recob::Hit>>& hitHandle,
  const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
  const art::ValidHandle<std::vector<recob::Track>> &trackHandle)
{
  art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
  art::FindManyP<recob::Hit> track2hitAssn(trackHandle, e, fTrack2HitLabel);

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
        hit.TruePdg     = fslcHitID2TruePdgMap[h.key()];
        hit.TrueID      = fslcHitID2G4IDMap[h.key()];
        hit.PfpID       = -1; // No associated PFP

        fmuHits[hit.PlaneID].push_back(hit);
    }
  }

}

bool sbnd::MichelTagging::isMuonFV(){
  double trackEndX = fmuTracks.back().EndX;
  double trackEndY = fmuTracks.back().EndY;
  double trackEndZ = fmuTracks.back().EndZ;

  std::cout << "Reco concatenated muon end position: " << trackEndX << ", " << trackEndY << ", " << trackEndZ << std::endl;

  return ((trackEndX > -200 + fFV) & (trackEndX < 200 - fFV) &
          (trackEndY > -200 + fFV) & (trackEndY < 200 - fFV) &
          (trackEndZ > 0 + fFV) & (trackEndZ < 500 - fFV));
}


// Main analyze function
void sbnd::MichelTagging::analyze(art::Event const& e)
{

  // ==============================================================
  // === EACH ENTRY OF OUTPUT TREE WILL CORRESPOND WITH A SLICE ===
  // ==============================================================


  // Fill event ID information
  feventID = e.id().event();
  frunID = e.id().run();
  fsubrunID = e.id().subRun();

  std::cout << "ANALYZING EVENT " << frunID << ":" << fsubrunID << ":" << feventID << std::endl;


  // ==============================================================
  // === LOAD HANDLES 
  // ==============================================================

  // Get MCT handle
  art::ValidHandle<std::vector<simb::MCTruth>> MCTHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);  
  art::ValidHandle<std::vector<simb::MCParticle>> MCPHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);

  // Pandora Handles
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  art::ValidHandle<std::vector<recob::Vertex>> vertexHandle = e.getValidHandle<std::vector<recob::Vertex>>(fVertexLabel);
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
  art::ValidHandle<std::vector<blip::Blip>> blipHandle = e.getValidHandle<std::vector<blip::Blip>>(fBlipLabel);

  // Get associations
  art::FindManyP<simb::MCParticle> mct2mcpAssn(MCTHandle, e, fMCParticleLabel);
  art::FindManyP<recob::PFParticle> slice2pfpAssn(sliceHandle, e, fPFParticleLabel);
  art::FindManyP<recob::Vertex> pfp2vertexAssn(pfpHandle, e, fVertexLabel);
  art::FindManyP<recob::Hit> slice2hitAssn(sliceHandle, e, fSliceLabel);
  art::FindOneP<sbn::TPCPMTBarycenterMatch> slice2bcfmAssns(sliceHandle, e, fTPCPMTBarycenterFMModuleLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn(pfpHandle, e, fPFParticleMetadataLabel);


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
      // Clear all variables
      clearVars();

      // Get Slice ID
      fslcID = slc->ID();
      std::cout << "  Analyzing slice " << fslcID << std::endl;

      // Get PFPs and hits
      std::vector<art::Ptr<recob::PFParticle>> slcPfpVec = slice2pfpAssn.at(slc.key());
      std::vector<art::Ptr<recob::Hit>> slcHitVec        = slice2hitAssn.at(slc.key());

      if(slcHitVec.empty()) continue;

      // Set hits map
      setSlcHitsMaps(e, slcHitVec, MCPHandle);

      // Get slice truthmatching
      getSliceTrueInfo(e, slcHitVec, MCTHandle);

      if (fgenIndex != 1) continue;
      if (fnuCCNC != 0) continue;

      // Get slice reco vertex
      getSlcVertex(e, slcPfpVec, pfpHandle);

      // Quality cuts
      // -----------------
      // Check if slice has a primary muon
      if (!fmuFound) {
          std::cout << "  Slice " << fslcID << " has no primary muon. Skipping slice." << std::endl;
          continue;
      }

      // Check if primary muon length is above threshold
      if (fprimMuLength < fMinMuLength) {
          std::cout << "  Slice " << fslcID << " primary muon length too short: " << fprimMuLength << " cm. Skipping slice." << std::endl;
          continue;
      }

      // Check if slice vertex is too far from true vertex
      double vtxDist = std::hypot(fslcVtxX - fslcTrueVtxX, fslcVtxY - fslcTrueVtxY, fslcVtxZ - fslcTrueVtxZ);
      if (vtxDist > fMaxVtxDist) {
          std::cout << "  Slice " << fslcID << " vertex too far from true vertex: " << vtxDist << " cm. Skipping slice." << std::endl;
          continue;
      }
      // -----------------
      
      // Analyze pfp
      if (slcPfpVec.empty()) {
          std::cout << "  Slice " << fslcID << " has no pfps. Skipping slice." << std::endl;
          continue; 
      }
      std::cout << "  Analyzing pfps..." << std::endl;
      analyzePFPs(e, slcPfpVec, pfpHandle, trackHandle, hitHandle);

      // Get unassociated hits
      std::cout << "  Getting unassociated hits..." << std::endl;
      getUnassociatedHits(e, slcPfpVec, slcHitVec, hitHandle, pfpHandle, trackHandle);

      // Check muon is in fv
      if (!isMuonFV()) {
          std::cout << "  Slice " << fslcID << " primary muon is not in fiducial volume. Skipping slice." << std::endl;
          continue;
      }

      // Save hits variables
      for (const MuHit& hit : fmuHits[2]){
        fhitPlaneID.push_back(hit.PlaneID);
        fhitTpcID.push_back(hit.TpcID);
        fhitChannelID.push_back(hit.ChannelID);
        fhitPeakTime.push_back(hit.PeakTime);
        fhitIntegral.push_back(hit.Integral);
        fhitRMS.push_back(hit.RMS);
        fhitTruePdg.push_back(hit.TruePdg);
        fhitPfpID.push_back(hit.PfpID);
        fhitTrueID.push_back(hit.TrueID);
      }

      // Save calo variables
      for (const MuCalo& calo : fmuCalo[2]){
        fcaloPlaneID.push_back(calo.PlaneID);
        fcaloTpcID.push_back(calo.TpcID);
        fcaloChannelID.push_back(calo.ChannelID);
        fcaloPeakTime.push_back(calo.PeakTime);
        fcaloTrueID.push_back(calo.TrueID);
        fcaloTruePdg.push_back(calo.TruePdg);
        fcalodEdx.push_back(calo.dEdx);
        fcaloResRange.push_back(calo.ResRange);
        fcaloX.push_back(calo.X);
        fcaloY.push_back(calo.Y);
        fcaloZ.push_back(calo.Z);
      }

      fnMuTracks = fmuTracks.size();

      // ------- Compute michel vars ------

      std::array<double, 3> recoVtx = {fslcVtxX, fslcVtxY, fslcVtxZ};

      std::cout << "  Running michel driver..." << std::endl;
      std::cout << "  mu tracks size: " << fmuTracks.size() << std::endl;
      MichelDriver michelDriver(fmichelParams, fmuTracks, fmuHits, fmuCalo, fblips, recoVtx, fmichelBDT);
      michelDriver.Compute();
      fmichelFeatures = michelDriver.Features();
      
      auto bragg = michelDriver.getBragg();
      auto pca = michelDriver.getPca();

      std::cout << "Michel features computed: " << fmichelFeatures.featuresComputed << std::endl;
      std::cout << "Michel features" << std::endl;
      fmichelFeatures.Print();
      std::cout << "Michel score: " << fmichelFeatures.score << std::endl;
      std::cout << "True muon decay: " << fmuDecay << std::endl;

      // Save bragg information
      fbraggFound = bragg.Found;
      fbraggIdx = bragg.Idx;
      fbraggChannelID = bragg.ChannelID;
      fbraggPeakTime = bragg.PeakTime;
      fbraggTpcID = bragg.TpcID;
      fbraggdEdx = bragg.dEdx;
      fbraggResRange = bragg.ResRange;
      fbraggX = bragg.X;
      fbraggY = bragg.Y;
      fbraggZ = bragg.Z;

      // Save pca information
      fpcaValid = pca.Valid;
      fpcaDirX = pca.dirX;
      fpcaDirY = pca.dirY;

      // Save angle vars
      fangles = michelDriver.getAngles();

      // Save radial vars
      fradii = michelDriver.getRadii();
      fcharge = michelDriver.getCharge();
      fdQdr = michelDriver.getdQdr();
      fsmoothdQdr = michelDriver.getSmoothdQdr();
      fdropIdx = michelDriver.getDropIdx();
      fmindQdr = michelDriver.getMindQdr();

      // Fill the tree
      fTree->Fill();
  }
}


DEFINE_ART_MODULE(sbnd::MichelTagging)
