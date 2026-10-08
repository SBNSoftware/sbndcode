#include "LightPropagationCorrection_module.hh"

#include <stdexcept>

namespace sbnd {
    class LightPropagationCorrection;
}

sbnd::LightPropagationCorrection::LightPropagationCorrection(fhicl::ParameterSet const& p)
    : EDProducer{p},
    fReco2Label( p.get<std::string>("Reco2Label") ),
    fOpT0FinderModuleLabel( p.get<std::string>("OpT0FinderModuleLabel") ),
    fTPCPMTBarycenterFMModuleLabel( p.get<std::string>("TPCPMTBarycenterFMModuleLabel") ),
    fOpFlashLabel_tpc0 ( p.get<std::string>("OpFlashLabel_tpc0") ),
    fOpFlashLabel_tpc1 ( p.get<std::string>("OpFlashLabel_tpc1") ),
    fSpacePointLabel( p.get<std::string>("SpacePointLabel") ),
    fTPCTrackLabel( p.get<std::string>("TPCTrackLabel") ),
    fParticleIDLabel( p.get<std::string>("ParticleIDLabel") ),
    fCalorimetryLabel( p.get<std::string>("CalorimetryLabel") ),
    fPandoraRangeLabel( p.get<std::string>("PandoraRangeLabel") ),
    fPandoraMCSLabel( p.get<std::string>("PandoraMCSLabel") ),
    fOpHitsModuleLabel( p.get<std::string>("OpHitsModuleLabel") ),
    fUseCheatMC( p.get<bool>("UseCheatMC") ),
    fUseMCVIS( p.get<bool>("UseMCVIS") ),
    fMCModuleLabel( p.get<std::string>("MCModuleLabel") ),
    fMCTruthModuleLabel( p.get<std::vector<std::string>>("MCTruthModuleLabel") ),
    fMCTruthInstanceLabel( p.get<std::vector<std::string>>("MCTruthInstanceLabel") ),
    fMCTruthOrigin( p.get<std::vector<int>>("MCTruthOrigin") ),
    fMCTruthPDG( p.get<std::vector<int>>("MCTruthPDG") ),
    fSimDepProducer( p.get<std::string>("SimDepProducer") ),
    fFlashMatchingTool( p.get<std::string>("FlashMatchingTool") ),
    fSaveCorrectionTree( p.get<bool>("SaveCorrectionTree") ),
    fDoPIDCorrection( p.get<bool>("DoPIDCorrection") ),
    fSpeedOfLight( p.get<double>("SpeedOfLight") ),
    fVGroupVIS( p.get<double>("VGroupVIS") ),
    fVGroupVUV( p.get<double>("VGroupVUV") ),
    fNuScoreThreshold( p.get<double>("NuScoreThreshold") ),
    fFMScoreThreshold( p.get<double>("FMScoreThreshold") ),
    fDebug( p.get<bool>("Debug", false) )
    // 
    // More initializers here.
{
    fNOpChannels = fWireReadout.NOpChannels();

    for(unsigned int opch=0; opch<fNOpChannels; opch++){
        auto pdCenter = fWireReadout.OpDetGeoFromOpChannel(opch).GetCenter();
        fOpDetID.push_back(opch);
        fOpDetX.push_back(pdCenter.X());
        fOpDetY.push_back(pdCenter.Y());
        fOpDetZ.push_back(pdCenter.Z());
    }

    auto const& tpc = art::ServiceHandle<geo::Geometry>()->TPC();
    fDriftDistance = tpc.DriftDistance();
    fKinkDistance = 0.5*fDriftDistance*(1-fVGroupVUV/fVGroupVIS);
    fVGroupVUV_I = 1./fVGroupVUV;
    fVISLightPropTime = fDriftDistance/fVGroupVIS;

    // Initialize flash algo for both TPCs 
    auto const flash_algo  = p.get<std::string>("FlashFinderAlgo");
    auto const flash_pset_tpc0 = p.get<lightana::Config_t>("AlgoConfig_tpc0");
    auto algo_ptr_tpc0 = ::lightana::FlashAlgoFactory::get().create(flash_algo,flash_algo);
    algo_ptr_tpc0->Configure(flash_pset_tpc0);
    _mgr_tpc0.SetFlashAlgo(algo_ptr_tpc0);

    auto const flash_pset_tpc1 = p.get<lightana::Config_t>("AlgoConfig_tpc1");
    auto algo_ptr_tpc1 = ::lightana::FlashAlgoFactory::get().create(flash_algo,flash_algo);
    algo_ptr_tpc1->Configure(flash_pset_tpc1);
    _mgr_tpc1.SetFlashAlgo(algo_ptr_tpc1);

    //Initialize flash geo tool
    auto const flashgeo_pset = p.get<lightana::Config_t>("FlashGeoConfig");
    _flashgeo = art::make_tool<lightana::FlashGeoBase>(flashgeo_pset);
    
    //Initialize flash t0 tool
    auto const flasht0_pset = p.get<lightana::Config_t>("FlashT0Config");
    _flasht0calculator = art::make_tool<lightana::FlashT0Base>(flasht0_pset);
    
    produces< std::vector<sbn::CorrectedOpFlashTiming> >();
    produces<art::Assns<recob::Slice, sbn::CorrectedOpFlashTiming>>();
    produces<art::Assns<recob::OpFlash, sbn::CorrectedOpFlashTiming>>();
}

void sbnd::LightPropagationCorrection::produce(art::Event & e)
{
    fEvent = e.id().event();
    fRun = e.id().run();
    fSubrun = e.id().subRun();
    _flashgeo->InitializeFlashGeoAlgo();

    if(fSaveCorrectionTree && fUseCheatMC)
    {
        SaveTrueTrajectory(e);
    }

    if(fUseCheatMC) GetMCNeutrino(e);

    std::unique_ptr< std::vector<sbn::CorrectedOpFlashTiming> > correctedOpFlashTimes (new std::vector<sbn::CorrectedOpFlashTiming>);
    art::PtrMaker<sbn::CorrectedOpFlashTiming> make_correctedopflashtime_ptr{e};
    std::unique_ptr< art::Assns<recob::Slice, sbn::CorrectedOpFlashTiming>> newCorrectedOpFlashTimingSliceAssn (new art::Assns<recob::Slice, sbn::CorrectedOpFlashTiming>);
    std::unique_ptr< art::Assns<recob::OpFlash, sbn::CorrectedOpFlashTiming>> newCorrectedOpFlashTimingOpFlashAssn (new art::Assns<recob::OpFlash, sbn::CorrectedOpFlashTiming>);

    // --- Read Recob Slice
    ::art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fReco2Label, sliceHandle);
    // Slice to OpT0Finder
    //Get the handle for making the assns later on
    ::art::Handle<std::vector<sbn::OpT0Finder>> opt0Handle;
    e.getByLabel(fOpT0FinderModuleLabel, opt0Handle);
    // Slice to TPCPMTBarycenterFM
    ::art::Handle<std::vector<sbn::TPCPMTBarycenterMatch>> tpcpmtbarycenterfmHandle;
    e.getByLabel(fTPCPMTBarycenterFMModuleLabel, tpcpmtbarycenterfmHandle);

    //Read PFPs
    ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    e.getByLabel(fReco2Label, pfpHandle);
    //Read Recob Tracks
    ::art::Handle<std::vector<recob::Track>> trackHandle;
    e.getByLabel(fTPCTrackLabel, trackHandle);
    //Read OpFlash Handle
    art::Handle< std::vector<recob::OpFlash> > opflashListHandle_tpc0;
    e.getByLabel(fOpFlashLabel_tpc0, opflashListHandle_tpc0);
    art::Handle< std::vector<recob::OpFlash> > opflashListHandle_tpc1;
    e.getByLabel(fOpFlashLabel_tpc1, opflashListHandle_tpc1);

    art::FindManyP<sbn::OpT0Finder> slice_opt0finder_assns(sliceHandle, e, fOpT0FinderModuleLabel);
    // Slice to TPCPMTBarycenterFM
    art::FindManyP<sbn::TPCPMTBarycenterMatch> slice_tpcpmtbarycentermatching_assns(sliceHandle, e, fTPCPMTBarycenterFMModuleLabel);
    // Slice to hits
    art::FindManyP<recob::Hit> slice_hit_assns (sliceHandle, e, fReco2Label);
    //Slice to PFParticles association
    art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fReco2Label);
    //PFP to vertex
    art::FindManyP<recob::Vertex> pfp_vertex_assns(pfpHandle, e, fReco2Label);
    //PFP to space points
    art::FindManyP<recob::SpacePoint> pfp_sp_assns(pfpHandle, e, fSpacePointLabel);
    //PF to track
    art::FindManyP<recob::Track> pfp_track_assns (pfpHandle, e, fTPCTrackLabel);
    // Track to PID
    art::FindManyP<anab::ParticleID> track_to_pid_assns(trackHandle, e, fParticleIDLabel);
    // Track to calo
    art::FindManyP<anab::Calorimetry> track_to_calo_assns(trackHandle, e, fCalorimetryLabel);
    // Get the rangeP track assns
    art::InputTag muon_range_tag(fPandoraRangeLabel, "muon");   // o "pion", "proton"
    art::InputTag proton_range_tag(fPandoraRangeLabel, "proton");   // o "pion", "proton"

    // Get the MCS track assns
    art::InputTag muon_MCS_tag(fPandoraMCSLabel, "muon");   // o "pion", "proton"
    art::InputTag proton_MCS_tag(fPandoraMCSLabel, "proton");   // o "pion", "proton"
  
    art::FindManyP<sbn::RangeP> track_rangeP_assns_muon(trackHandle, e, muon_range_tag);
    art::FindManyP<sbn::RangeP> track_rangeP_assns_proton(trackHandle, e, proton_range_tag);

    art::FindManyP<recob::MCSFitResult> track_MCS_assns_muon(trackHandle, e, muon_MCS_tag);
    art::FindManyP<recob::MCSFitResult> track_MCS_assns_proton(trackHandle, e, proton_MCS_tag);

    //OpFlash to OpHit
    flashToOpHitAssns_tpc0 = std::make_unique<art::FindManyP<recob::OpHit>>( opflashListHandle_tpc0, e, fOpFlashLabel_tpc0);
    flashToOpHitAssns_tpc1 = std::make_unique<art::FindManyP<recob::OpHit>>(opflashListHandle_tpc1, e, fOpFlashLabel_tpc1);

    // PFP Metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_to_metadata(pfpHandle, e, fReco2Label);

    // --- Store candidate slices
    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);

    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfpVect;
    // --- Get the candidate slices
    for(size_t ix=0; ix<sliceVect.size(); ix++){
        //std::cout << " New slice " << std::endl;
        ResetSliceInfo();
        // --- Get the slice
        auto & slice = sliceVect[ix];
        // Now I need to get all the hits associated to this flash and get the timing for all of them
        // Get the slices PFPs
        double _sliceMaxNuScore = -9999.;
        pfpVect = slice_pfp_assns.at(slice.key());
        int nPFPs = pfpVect.size();
        // Check wether there is a cathode-crosser track within the slice
        if ( nPFPs != 0 ) {
        art::FindOne<anab::T0> f1T0( {pfpVect.at(0)}, e, fReco2Label);
        }

        fNeutrinoID = -1;
        // Get the neutrino ID and vertex position
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            fTruthMatchedTrackID=-1;
            const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp.key());
            if(pfp->IsPrimary() &&( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ){
            const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp.key());
            for (auto const pfpMeta : pfpMetaVec) {
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                if (propertiesMap.count("NuScore")) _fNuScore = propertiesMap.at("NuScore");
                if(_fNuScore>_sliceMaxNuScore) _sliceMaxNuScore = _fNuScore;
            }
            std::vector< art::Ptr<recob::Vertex> > vertexVec = pfp_vertex_assns.at(pfp.key());
            for(const art::Ptr<recob::Vertex> &ver : vertexVec){
                geo::Point_t xyz_vertex = ver->position();
                fRecoVx= xyz_vertex.X();
                fRecoVy= xyz_vertex.Y();
                fRecoVz= xyz_vertex.Z();
            }
            fNeutrinoID = pfp->Self();
            }
        }
        
        // Get informaiton on the remaining PFPs
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            if (fNeutrinoID == static_cast<int>(pfp->Self())) continue;
            bool tracksuccess = false;
            //if(pfp->Self()==neutrinoID) continue; // We already got neutrino 
            // Check if it is a clear cosmic
            bool isClearCosmic = false;
            const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp.key());
            // Check if the PFP is a clear cosmic
            for (auto const pfpMeta : pfpMetaVec) {
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                if(propertiesMap.count("IsClearCosmic")){
                isClearCosmic=true;
                }
            }
            if(isClearCosmic) continue;
            //Read the tracks and store the PFParticle start/end points
            bool pfpistrack = ::lar_pandora::LArPandoraHelper::IsTrack(pfp);
            //bool pfpisshower = ::lar_pandora::LArPandoraHelper::IsShower(pfp);
            // Neutrino daughter particles
            for (auto const pfpMeta : pfpMetaVec) {
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                if(propertiesMap.count("TrackScore")){
                    double track_score=propertiesMap.at("TrackScore");
                    std::cout << " track_score " << track_score << std::endl;
                }
            }

            // For every PFP look for the true MC Particle 
            std::vector<art::Ptr<recob::SpacePoint>> PFPSpacePointsVect = pfp_sp_assns.at(pfp.key());
            //Get the SP Hit assns
            art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
            std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
            e.getByLabel(fSpacePointLabel, eventSpacePoints);
            art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);
            art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);
            std::vector<art::Ptr<recob::Hit>> PFPhits;
            for (auto const& sp : PFPSpacePointsVect) {
                std::vector<art::Ptr<recob::Hit>> hits = SPToHitAssoc.at(sp.key());
                PFPhits.insert(PFPhits.end(), hits.begin(), hits.end());
            }
            if(fUseCheatMC) GetTruthMatchedID(e, PFPhits);

            if(!fUseCheatMC)
            {
                if(fDoPIDCorrection && pfpistrack)
                {
                    // Neutrino daughter particles
                    // If it's a track then go for PID and momentum calculation to get particle propagation time
                    std::vector<art::Ptr<recob::Track>> track_v = pfp_track_assns.at(pfp.key());
                    if(track_v.size()>1)
                        throw art::Exception(art::errors::LogicError) << "Multiple tracks associated to a PFP. This is not expected.";
                    if (track_v.size() == 1) {
                        int pdg = GetTrackPID(track_to_pid_assns, track_v[0]);
                        if (pdg != -1) {
                            double initial_momentum = GetTrackMomentum( pdg, track_rangeP_assns_muon, track_rangeP_assns_proton, track_MCS_assns_muon, track_MCS_assns_proton ,  track_v[0]);
                            if (initial_momentum > 0) tracksuccess = GetParticlePropagationTime( pdg, initial_momentum, track_to_calo_assns, track_v[0], pfp->Self(), pfp->Parent() );
                        }
                    }
                }
                if(!tracksuccess)
                {
                    fAllPFPsCorrected=false;
                    // If it's a shower jsut assume speed of light propagation
                    std::vector<art::Ptr<recob::SpacePoint>> PFPSpacePointsVect = pfp_sp_assns.at(pfp.key());
                    //Get the SP Hit assns
                    art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
                    std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
                    e.getByLabel(fSpacePointLabel, eventSpacePoints);
                    art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);
                    art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);
                    GetParticlePropagationTimeLite(pfp->Self(), PFPSpacePointsVect, SPToHitAssoc);
                }
            }
            else{
                GetParticlePropagationTimeMC(e, PFPhits, pfp->Self());
            }

        }
        
        if(_sliceMaxNuScore<fNuScoreThreshold){
            ResetSliceInfo();
            continue; // Skip to the next slice if the nu score is below threshold
        }

        if(fRecoVx == -99999. || fRecoVy == -99999. || fRecoVz == -99999.) {
                ResetSliceInfo();
                continue;
        }
        
        this->GetPropagationTimeCorrectionPerChannel();
        
        // Get all the OpT0 objects associated to the slice
        std::vector<art::Ptr<recob::OpFlash>> flashFM;
        if(fFlashMatchingTool == "OpT0Finder" ){
            const std::vector< art::Ptr<sbn::OpT0Finder> > slcOpT0Finder = slice_opt0finder_assns.at( slice.key() );
            if(slcOpT0Finder.size() == 0) continue; 
            size_t OpT0Idx = HighestOpT0ScoreIdx(slcOpT0Finder);
            // Get the flash OpT0 association
            art::FindManyP<recob::OpFlash> opflash_opt0finder_assns(opt0Handle, e, fOpT0FinderModuleLabel);
            flashFM = opflash_opt0finder_assns.at( slcOpT0Finder[OpT0Idx].key() );
            if(flashFM.size() > 1){
                throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same OpT0Finder object. This is not expected.";
            }
            _fFMScore = slcOpT0Finder[OpT0Idx]->score;
            if(_fFMScore < fFMScoreThreshold){
                ResetSliceInfo();
                continue;
            }
        }
        else if(fFlashMatchingTool == "BarycenterFM")
        {
            const std::vector< art::Ptr<sbn::TPCPMTBarycenterMatch> > slcTPCPMTBarycenter = slice_tpcpmtbarycentermatching_assns.at( slice.key() );
            if(slcTPCPMTBarycenter.size() == 0) continue;
            art::FindManyP<recob::OpFlash> opflash_tpcpmtbarycenter_assns(tpcpmtbarycenterfmHandle, e, fTPCPMTBarycenterFMModuleLabel);
            size_t BFMIdx = HighestBFMScoreIdx(slcTPCPMTBarycenter);
            flashFM = opflash_tpcpmtbarycenter_assns.at( slcTPCPMTBarycenter[BFMIdx].key() );
            if(flashFM.size() > 1){
                throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same TPCPMTBarycenterFM object. This is not expected.";
            }
            if(flashFM.size() == 0){
                throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same TPCPMTBarycenterFM object. This is not expected.";
            }
            _fFMScore = slcTPCPMTBarycenter[BFMIdx]->score;
            if(_fFMScore < fFMScoreThreshold){
                ResetSliceInfo();
                continue;
            }
        }
        else throw art::Exception(art::errors::LogicError) << " Flash matching tool " <<  fFlashMatchingTool << " not supported ." << std::endl; 

        float minTimeDiff = std::numeric_limits<float>::max();
        
        int minIdx = -1;

        sbn::CorrectedOpFlashTiming correctedOpFlashTiming_tpc0;
        sbn::CorrectedOpFlashTiming correctedOpFlashTiming_tpc1;
        if(flashFM[0]->XCenter()<0)
        {
            for(size_t i=0; i<opflashListHandle_tpc1->size(); ++i){
                double timeDiff = std::abs(flashFM[0]->Time() - opflashListHandle_tpc1->at(i).Time());
                if(timeDiff < minTimeDiff){
                    minTimeDiff = timeDiff;
                    minIdx = i;
                }
            }
            // If FM flash is on TPC0 then
            CorrectOpFlash(flashFM[0], correctedOpFlashTiming_tpc0, true, 0);
            if(opflashListHandle_tpc1->size()>0)
            {
                CorrectOpFlash(art::Ptr<recob::OpFlash>(opflashListHandle_tpc1, minIdx), correctedOpFlashTiming_tpc1, false, 1);
            }

        }
        else
        {
            for(size_t i=0; i<opflashListHandle_tpc0->size(); ++i){
                double timeDiff = std::abs(flashFM[0]->Time() - opflashListHandle_tpc0->at(i).Time());
                if(timeDiff < minTimeDiff){
                    minTimeDiff = timeDiff;
                    minIdx = i;
                }
            }
            CorrectOpFlash(flashFM[0], correctedOpFlashTiming_tpc1, true, 1);
            if(opflashListHandle_tpc0->size()>0)
            {
                CorrectOpFlash(art::Ptr<recob::OpFlash>(opflashListHandle_tpc0, minIdx), correctedOpFlashTiming_tpc0, false, 0);
            }
        }
        if(minTimeDiff > 0.1){
            std::cout << " Max time diff not compatible with simultaneous flashes" << std::endl;
            continue;
        }

        if(std::isnan(correctedOpFlashTiming_tpc0.OpFlashT0Corrected) || std::isnan(correctedOpFlashTiming_tpc1.OpFlashT0Corrected))
        {
            std::cout << " NaN values found in the corrected flash timings. Skip this slice." << std::endl;
            continue;
        }

        correctedOpFlashTimes->emplace_back(std::move(correctedOpFlashTiming_tpc0));
        auto ptr0 = make_correctedopflashtime_ptr(correctedOpFlashTimes->size() - 1);

        correctedOpFlashTimes->emplace_back(std::move(correctedOpFlashTiming_tpc1));
        auto ptr1 = make_correctedopflashtime_ptr(correctedOpFlashTimes->size() - 1);

        // Slice -> ambos timings
        newCorrectedOpFlashTimingSliceAssn->addSingle(slice, ptr0);
        newCorrectedOpFlashTimingSliceAssn->addSingle(slice, ptr1);

        // Flash -> ambos timings
        newCorrectedOpFlashTimingOpFlashAssn->addSingle(flashFM[0], ptr0);
        newCorrectedOpFlashTimingOpFlashAssn->addSingle(flashFM[0], ptr1);
    }
    if(fSaveCorrectionTree) fTree->Fill();
    ResetEventVars();

    e.put(std::move(correctedOpFlashTimes));
    e.put(std::move(newCorrectedOpFlashTimingSliceAssn));
    e.put(std::move(newCorrectedOpFlashTimingOpFlashAssn));
    return;
}

void sbnd::LightPropagationCorrection::beginJob()
{
    if(fSaveCorrectionTree)
    {
        fTree = tfs->make<TTree>("PMTWaveformFilteranalyzer", "PMT Waveform Filter Analyzer Tree");
        fTree->Branch("eventID", &fEvent, "eventID/i");
        fTree->Branch("runID", &fRun, "runID/i");
        fTree->Branch("subrunID", &fSubrun, "subrunID/i");
        fTree->Branch("NuScore", &fNuScore);
        fTree->Branch("FMScore", &fFMScore);
        fTree->Branch("OpFlashTimeOld", &fOpFlashTimeOld);
        fTree->Branch("OpFlashTimeNew", &fOpFlashTimeNew);
        fTree->Branch("OpFlashXCenter", &fOpFlashXCenter);
        fTree->Branch("OpFlashYCenter", &fOpFlashYCenter);
        fTree->Branch("OpFlashZCenter", &fOpFlashZCenter);
        fTree->Branch("OpFlashPE", &fOpFlashPE);
        fTree->Branch("SliceVx", &fSliceVx);
        fTree->Branch("SliceVy", &fSliceVy);
        fTree->Branch("SliceVz", &fSliceVz);
        fTree->Branch("SliceSPX", &fSliceSPX);
        fTree->Branch("SliceSPY", &fSliceSPY);
        fTree->Branch("SliceSPZ", &fSliceSPZ);
        fTree->Branch("SliceSPT", &fSliceSPT);
        fTree->Branch("SliceSPMomentum", &fSliceSPMomentum);
        fTree->Branch("SliceSPEnergy", &fSliceSPEnergy);
        fTree->Branch("SliceSPBeta", &fSliceSPBeta);
        fTree->Branch("SliceSPID", &fSliceSPID);
        fTree->Branch("SliceSPTMID", &fSliceSPTMID);
        fTree->Branch("TrueMomentum", &fTrueMomentum);
        fTree->Branch("TruePID", &fTruePID);
        fTree->Branch("TruthMatchedX", &fTMX);
        fTree->Branch("TruthMatchedY", &fTMY);
        fTree->Branch("TruthMatchedZ", &fTMZ);
        fTree->Branch("TruthMatchedT", &fTMT);
        fTree->Branch("TruthMatchedE", &fTME);
        fTree->Branch("TruthMatchedTID", &fTMID);
        fTree->Branch("truedE_X", &truedE_X);
        fTree->Branch("truedE_Y", &truedE_Y);
        fTree->Branch("truedE_Z", &truedE_Z);
        fTree->Branch("truedE_ID", &truedE_ID);
        fTree->Branch("truedE_T", &truedE_T);
        fTree->Branch("truedE", &truedE_);
        fTree->Branch("TruthMatchedID", &fTruthMatchedID);
        fTree->Branch("OpHitOldTime", &fOpHitOldTime);
        fTree->Branch("OpHitNewTime", &fOpHitNewTime);
        fTree->Branch("OpHitPE", &fOpHitPE);
        fTree->Branch("OpHitOpCh", &fOpHitOpCh);
    }
}

void sbnd::LightPropagationCorrection::endJob()
{
}

void sbnd::LightPropagationCorrection::ResetEventVars()
{
    fTrueVt=0.0;
    fTrueVx = 0.0;
    fTrueVy = 0.0;
    fTrueVz = 0.0;

    if(fSaveCorrectionTree)
    {
        fEvent = 0;
        fRun = 0;
        fSubrun = 0;
        _fNuScore = 0.0;
        fNuScore.clear();
        fFMScore.clear();
        fOpFlashTimeOld.clear();
        fOpFlashTimeNew.clear();
        fOpFlashXCenter.clear();
        fOpFlashYCenter.clear();
        fOpFlashZCenter.clear();
        fOpHitOldTime.clear();
        fOpHitNewTime.clear();
        fOpHitPE.clear();
        fOpHitOpCh.clear();
        fOpFlashPE.clear();
        fSliceVx.clear();
        fSliceVy.clear();
        fSliceVz.clear();
        fSliceSPX.clear();
        fSliceSPY.clear();
        fSliceSPZ.clear();
        fSliceSPT.clear();
        fSliceSPMomentum.clear();
        fSliceSPEnergy.clear();
        fSliceSPBeta.clear();
        fSliceSPID.clear();
        fSliceSPTMID.clear();
        fTrueMomentum.clear();
        fTruePID.clear();
        fTMX.clear();
        fTMY.clear();
        fTMZ.clear();
        fTMT.clear();
        fTME.clear();
        fTMID.clear();
        fTruthMatchedX.clear();
        fTruthMatchedY.clear();
        fTruthMatchedZ.clear();
        fTruthMatchedT.clear();
        fTruthMatchedE.clear();
        fTruthMatchedID.clear();

        truedE.clear();
        truedE_vecX.clear();
        truedE_vecY.clear();
        truedE_vecZ.clear();
        truedE_vecID.clear();
        truedE_vecT.clear();

        truedE_.clear();
        truedE_X.clear();
        truedE_Y.clear();
        truedE_Z.clear();
        truedE_ID.clear();
        truedE_T.clear();
    }
}

size_t sbnd::LightPropagationCorrection::HighestOpT0ScoreIdx(const std::vector< art::Ptr<sbn::OpT0Finder> > slcFM)
{
    // Gets the idx of the OpT0 object with the highest score 
    double highestOpT0Score = -99999.0; // Initialize to a negative value
    size_t highestIdx = 0;
    for(size_t jx=0; jx<slcFM.size(); jx++){
        if(slcFM[jx]->score > highestOpT0Score){
            highestOpT0Score = slcFM[jx]->score;
            highestIdx = jx;
        }
    }
    return highestIdx;
}

size_t sbnd::LightPropagationCorrection::HighestBFMScoreIdx(const std::vector< art::Ptr<sbn::TPCPMTBarycenterMatch> > slcFM)
{
    // Gets the idx of the OpT0 object with the highest score 
    double highesBFM0Score = -99999.0; // Initialize to a negative value
    size_t highestIdx = 0;
    for(size_t jx=0; jx<slcFM.size(); jx++){
        if(slcFM[jx]->score > highesBFM0Score){
            highesBFM0Score = slcFM[jx]->score;
            highestIdx = jx;
        }
    }
    return highestIdx;
}

void sbnd::LightPropagationCorrection::ResetSliceInfo()
{
    _fNuScore=-99999.0;
    fRecoVx = -99999.0;
    fRecoVy = -99999.0;
    fRecoVz = -99999.0;
    fSpacePointX.clear();
    fSpacePointY.clear();
    fSpacePointZ.clear();
    fSpacePointIntegral.clear();
    fSpacePointPFPID.clear();
    fSpacePointTMID.clear();
    fSpacePointPropagationTime.clear();
    fTimeCorrectionPerChannel.assign(fNOpChannels, 0.0); // Reset the time correction vector for each channel
    fParticlePropagationTimePerChannel.assign(fNOpChannels, 0.0); // Reset the particle propagation time vector for each channel
    fPhotonPropagationTimePerChannel.assign(fNOpChannels, 0.0); // Reset the photon propagation time vector for each channel
    fPhotonPropagationDistancePerChannel.assign(fNOpChannels, 0.0); // Reset the photon propagation distance vector for each channel
    fAllPFPsCorrected = true;
}

void sbnd::LightPropagationCorrection::GetPropagationTimeCorrectionPerChannel()
{
    // Implementation
    for(size_t opdet = 0; opdet < fOpDetID.size(); opdet++) {
        double _opDetX = fOpDetX[opdet];
        double _opDetY = fOpDetY[opdet];
        double _opDetZ = fOpDetZ[opdet];
        float minPropTime = 999999999.;
        float minPartPropTime = 999999999.;
        float minLightPropTime = 999999999.;
        float minPropDistance = 999999999.;
        bool foundSP = false;

        for(size_t sp=0; sp<fSpacePointX.size(); sp++)
        {
            bool isInSameTPC = (fSpacePointX[sp] * _opDetX) > 0;
            if(!isInSameTPC) continue; // Skip points not in the same TPC
            double dx = fSpacePointX[sp] - _opDetX;
            double dy = fSpacePointY[sp] - _opDetY;
            double dz = fSpacePointZ[sp] - _opDetZ;
            double distanceToOpDet = std::sqrt(dx*dx + dy*dy + dz*dz);
            double spToCathode;
            double cathodeToOpDet;
            if(fUseMCVIS){
                spToCathode = std::sqrt( fSpacePointX[sp]*fSpacePointX[sp]); // Distance from space point to cathode in mm
                cathodeToOpDet = std::sqrt(_opDetX*_opDetX + (dy)*(dy) + (dz)*(dz)); // Distance from cathode to OpDet in mm
            }
            else{
                spToCathode = std::sqrt( fSpacePointX[sp]*fSpacePointX[sp] + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from space point to cathode in mm
                cathodeToOpDet = std::sqrt(_opDetX*_opDetX + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from cathode to OpDet in mm
            }

            float lightPropTimeVIS = spToCathode/fVGroupVUV + cathodeToOpDet/fVGroupVIS; // Speed
            float lightPropTimeVUV = distanceToOpDet / fVGroupVUV; // Speed of light in mm/ns for VUV
            float lightPropTime = 0;
            float lightPropDistance = 0;
            const std::string pdType = fPDSMap.pdType(opdet);
            if(pdType=="pmt_coated" || pdType=="xarapuca_vuv")
                lightPropTime = std::min(lightPropTimeVIS, lightPropTimeVUV);
            else if(pdType=="pmt_uncoated" || pdType=="xarapuca_vis")
                lightPropTime = lightPropTimeVIS;
            else
                throw std::runtime_error("LightPropagationCorrection: unexpected pdType '" + pdType + "' for opdet " + std::to_string(opdet));

            lightPropDistance = distanceToOpDet;
            //float fastest_partPropTime = std::sqrt((fSpacePointX[sp]-fRecoVx)*(fSpacePointX[sp]-fRecoVx) + (fSpacePointY[sp]-fRecoVy)*(fSpacePointY[sp]-fRecoVy) + (fSpacePointZ[sp]-fRecoVz)*(fSpacePointZ[sp]-fRecoVz))/fSpeedOfLight;
            float partPropTime = fSpacePointPropagationTime[sp];
            float PropTime = lightPropTime + partPropTime;
            if(PropTime < minPropTime) {
                minPropTime = PropTime;
                minPartPropTime = partPropTime;    
                minLightPropTime = lightPropTime;
                minPropDistance = lightPropDistance;
                foundSP = true;
            }
        }
        if(!foundSP) {
            fTimeCorrectionPerChannel[opdet] = 0.0;
            fParticlePropagationTimePerChannel[opdet] = 0.0;
            fPhotonPropagationTimePerChannel[opdet] = 0.0;
            fPhotonPropagationDistancePerChannel[opdet] = 0.0;
            continue;
        }

        fTimeCorrectionPerChannel[opdet] = -minPropTime;
        fParticlePropagationTimePerChannel[opdet] = minPartPropTime;
        fPhotonPropagationTimePerChannel[opdet] = minLightPropTime;
        fPhotonPropagationDistancePerChannel[opdet] = minPropDistance;
    }
}


void sbnd::LightPropagationCorrection::CorrectOpHitTime(std::vector<art::Ptr<recob::OpHit>> OldOpHitList, std::vector<recob::OpHit> & newOpHitList)
{
        int _nophits = 0;
        _nophits += OldOpHitList.size();
        for (int i = 0; i < _nophits; ++i) {
            int opCh = OldOpHitList.at(i)->OpChannel();
            double channelCorrection = fTimeCorrectionPerChannel[opCh];
            double newPeakTime = OldOpHitList.at(i)->StartTime()+OldOpHitList.at(i)->RiseTime() + channelCorrection/1000;
            double newPeakTimeAbs = OldOpHitList.at(i)->PeakTimeAbs()+ channelCorrection/1000;
            double newStartTime = OldOpHitList.at(i)->StartTime() + channelCorrection/1000;
            double riseTime = OldOpHitList.at(i)->RiseTime();
            unsigned int frame = OldOpHitList.at(i)->Frame();
            double width = OldOpHitList.at(i)->Width();
            double area = OldOpHitList.at(i)->Area();
            double amplitude = OldOpHitList.at(i)->Amplitude();
            double pe = OldOpHitList.at(i)->PE();
            recob::OpHit newOpHit = recob::OpHit(opCh, newPeakTime, newPeakTimeAbs, newStartTime, riseTime, frame, width, area, amplitude, pe, 0.0);
            newOpHitList.push_back(newOpHit);
        }
}


void sbnd::LightPropagationCorrection::FillLiteOpHit(std::vector<recob::OpHit> const& OpHitList, std::vector<::lightana::LiteOpHit_t>& LiteOpHitList)
{
    for(auto const& oph : OpHitList) {
        ::lightana::LiteOpHit_t loph;
        loph.peak_time = oph.StartTime()+oph.RiseTime();
        loph.pe = oph.PE();
        loph.channel = oph.OpChannel();
        LiteOpHitList.emplace_back(std::move(loph));
    }
}


void sbnd::LightPropagationCorrection::FillCorrectionTree(double & newFlashTime, recob::OpFlash const& flash, std::vector<recob::OpHit> const& oldOpHitList, std::vector<recob::OpHit> const& newOpHitList){
    fOpHitOldTime.push_back({});
    fOpHitNewTime.push_back({});
    fOpHitPE.push_back({});
    fOpHitOpCh.push_back({});

    if(fDebug)
    {
        for(size_t i=0; i<fSpacePointX.size(); i++){
            fSliceSPX.push_back(fSpacePointX[i]);
            fSliceSPY.push_back(fSpacePointY[i]);
            fSliceSPZ.push_back(fSpacePointZ[i]);
            fSliceSPT.push_back(fSpacePointPropagationTime[i]);
            fSliceSPID.push_back(fSpacePointPFPID[i]);
            fSliceSPTMID.push_back(fSpacePointTMID[i]);
        }
        
        for(size_t i=0; i<fTruthMatchedX.size(); i++)
        {
            fTMX.push_back(fTruthMatchedX[i]);
            fTMY.push_back(fTruthMatchedY[i]);
            fTMZ.push_back(fTruthMatchedZ[i]);
            fTMT.push_back(fTruthMatchedT[i]);
            fTME.push_back(fTruthMatchedE[i]);
            fTMID.push_back(fTruthMatchedID[i]);
        }

        for(size_t i=0; i<truedE_vecX.size(); i++)
        {
            truedE_.push_back(truedE[i]); 
            truedE_X.push_back(truedE_vecX[i]); 
            truedE_Y.push_back(truedE_vecY[i]); 
            truedE_Z.push_back(truedE_vecZ[i]);
            truedE_ID.push_back(truedE_vecID[i]);    
            truedE_T.push_back(truedE_vecT[i]);
        }
        

        for(size_t i=0; i<oldOpHitList.size(); i++){
            fOpHitOldTime.back().push_back(oldOpHitList[i].StartTime()+oldOpHitList[i].RiseTime());
            fOpHitNewTime.back().push_back(newOpHitList[i].StartTime()+newOpHitList[i].RiseTime());
            fOpHitPE.back().push_back(oldOpHitList[i].PE());
            fOpHitOpCh.back().push_back(oldOpHitList[i].OpChannel());
        }
    }

    fNuScore.push_back(_fNuScore);
    fFMScore.push_back(_fFMScore);
    fOpFlashTimeOld.push_back(flash.Time());
    fOpFlashTimeNew.push_back(newFlashTime);
    fOpFlashXCenter.push_back(flash.XCenter());
    fOpFlashYCenter.push_back(flash.YCenter());
    fOpFlashZCenter.push_back(flash.ZCenter());
    fOpFlashPE.push_back(std::accumulate(flash.PEs().begin(), flash.PEs().end(), 0.0));
    fSliceVx.push_back(fRecoVx);
    fSliceVy.push_back(fRecoVy);
    fSliceVz.push_back(fRecoVz);

}

void sbnd::LightPropagationCorrection::CorrectOpFlash(art::Ptr<recob::OpFlash> const& flash, sbn::CorrectedOpFlashTiming &correctedOpFlashTiming, bool matched, int tpc)
{
    // Get the ophits associated to the flash
    std::vector<art::Ptr<recob::OpHit>> ophitlist;
    if(abs(flash->XCenter()) > 250) return;
    if(flash->XCenter()<0)
    {
        ophitlist = flashToOpHitAssns_tpc0->at(flash.key());
        _mgr = _mgr_tpc0; // Use the TPC 0 flash finder manager
    }
    else
    {
        ophitlist = flashToOpHitAssns_tpc1->at(flash.key());
        _mgr = _mgr_tpc1; // Use the TPC 1 flash finder manager
    }
    std::vector<recob::OpHit> newOpHitList;
    std::vector<recob::OpHit> oldOpHitList;
    for(const auto& ophit : ophitlist) {
        oldOpHitList.push_back(*ophit);
    }
    // Get the list of the corrected OpHits
    this->CorrectOpHitTime(ophitlist, newOpHitList);
    // Create the list of ophit lite to be used in the flash finder
    ::lightana::LiteOpHitArray_t ophits;
    this->FillLiteOpHit(newOpHitList, ophits);
    // Create the flash manager
    auto const flash_v = _mgr.RecoFlash(ophits);
    double originalFlashTime = flash->Time();
    double newFlashTime = 0.0;
    double particlePropTime = 0.0;
    double photonPropTime = 0.0;
    double photonPropDistance = 0.0;
    int NCoatedPMTs=0;
    int NUncoatedPMTs=0; 
    for(const auto& lflash :  flash_v) {
        // Get Flash Barycenter
        double Ycenter, Zcenter, Ywidth, Zwidth;
        _flashgeo->GetFlashLocation(lflash.channel_pe, Ycenter, Zcenter, Ywidth, Zwidth);
        // Get flasht0
        double flasht0 = lflash.time;
        // Refine t0 calculation
        flasht0 = _flasht0calculator->GetFlashT0(lflash.time, GetAssociatedLiteHits(lflash, ophits));
        recob::OpFlash flash(flasht0, lflash.time_err, flasht0,
                            ( flasht0) / 1600., lflash.channel_pe,
                            0, 0, 1, // this are just default values
                            100., -1., Ycenter, Ywidth, Zcenter, Zwidth);
        newFlashTime = flasht0;
        particlePropTime = _flasht0calculator->GetAverageMagnitude(fParticlePropagationTimePerChannel)/1000;
        photonPropDistance = _flasht0calculator->GetAverageMagnitude(fPhotonPropagationDistancePerChannel)/1000;
        photonPropTime = _flasht0calculator->GetAverageMagnitude(fPhotonPropagationTimePerChannel)/1000;
        NCoatedPMTs = _flasht0calculator->GetFlashNPMTs("pmt_coated");
        NUncoatedPMTs = _flasht0calculator->GetFlashNPMTs("pmt_uncoated");
        correctedOpFlashTiming.OpFlashT0 = originalFlashTime;
        std::cout << " Neutrino interaction time is " << fTrueVt << std::endl;
        std::cout << " original flash time " << 1000*originalFlashTime-135 << " new flash time " << 1000*newFlashTime-135 <<  " in tpc " << tpc << " and is matched " << matched << std::endl;
        correctedOpFlashTiming.OpFlashPE = flash.TotalPE();
        correctedOpFlashTiming.NuToFLight = (Zcenter/fSpeedOfLight)/1000;
        correctedOpFlashTiming.NuToFCharge = (fRecoVz/fSpeedOfLight)/1000;
        correctedOpFlashTiming.OpFlashT0Corrected = newFlashTime;
        correctedOpFlashTiming.ParticlePropagationTime = particlePropTime;
        correctedOpFlashTiming.PhotonPropagationTime = photonPropTime;
        correctedOpFlashTiming.PhotonPropagationDistance = photonPropDistance;
        correctedOpFlashTiming.MatchedOpFlash = matched;
        correctedOpFlashTiming.tpc = tpc;
        correctedOpFlashTiming.NCoatedPMTs = NCoatedPMTs;
        correctedOpFlashTiming.NUncoatedPMTs = NUncoatedPMTs;
        correctedOpFlashTiming.AllPFPsCorrected = fAllPFPsCorrected;
    }
    if(fSaveCorrectionTree){
        this->FillCorrectionTree(newFlashTime, *flash, oldOpHitList, newOpHitList);
    }
}

::lightana::LiteOpHitArray_t sbnd::LightPropagationCorrection::GetAssociatedLiteHits(::lightana::LiteOpFlash_t lite_flash, ::lightana::LiteOpHitArray_t lite_hits_v)
{
    ::lightana::LiteOpHitArray_t flash_hits_v;

    for(auto const& hitidx : lite_flash.asshit_idx) {
      flash_hits_v.emplace_back(std::move(lite_hits_v.at(hitidx)));
    }

    return flash_hits_v;
}


int sbnd::LightPropagationCorrection::GetTrackPID(art::FindManyP<anab::ParticleID> track_to_pid_assns, art::Ptr<recob::Track> track )
{
    std::vector<art::Ptr<anab::ParticleID>> pidV = track_to_pid_assns.at(track.key());

    double muonChi2   = -99999.;
    double protonChi2 = -99999.;
    
    for (size_t j = 0; j < pidV.size(); ++j) {

        const auto& pid = *pidV[j];
        int plane = pid.PlaneID().Plane;
        // Use only collection plane
        if(plane!=2) continue;

        const auto& alg_score_vector = pid.ParticleIDAlgScores();
        for (const auto& alg_score : alg_score_vector) {
            if (alg_score.fAlgName != "Chi2") continue;
            switch (std::abs(alg_score.fAssumedPdg)) {
                case 13:
                    muonChi2 = alg_score.fValue;
                    break;
                case 2212:
                    protonChi2 = alg_score.fValue;
                    break;
                default:
                    break;
            }
        }
    }
    
    if(protonChi2== -99999. || muonChi2== -99999. )
        return -1;

    if(muonChi2<25 && protonChi2>95) return 13;
    else if( protonChi2<=95) return 2212;
    else return -1;

}


double sbnd::LightPropagationCorrection::GetTrackMomentum(int pdg, art::FindManyP<sbn::RangeP> track_rangeP_assns_muon,  art::FindManyP<sbn::RangeP> track_rangeP_assns_proton, art::FindManyP<recob::MCSFitResult> track_MCS_assns_muon , art::FindManyP<recob::MCSFitResult> track_MCS_assns_proton  , art::Ptr<recob::Track> track )
{
    double p_muon = -99999;
    double p_proton = -99999;

    bool isContained; 

    double track_end_x = track->End().X();
    double track_end_y = track->End().Y();
    double track_end_z = track->End().Z();

    isContained= (abs(track_end_x)<200 && abs(track_end_y)<200 && (track_end_z>0 && track_end_z<500));

    if(isContained)
    {
        auto const& muonRange = track_rangeP_assns_muon.at(track.key());
        if (!muonRange.empty()){
            p_muon = muonRange[0]->range_p;
        }
        auto const& protonRange = track_rangeP_assns_proton.at(track.key());
        if (!protonRange.empty()) {
            p_proton = protonRange[0]->range_p;
        }
    }
    else
    {
        auto const& muonMCS = track_MCS_assns_muon.at(track.key());
        if (!muonMCS.empty())
          p_muon = muonMCS[0]->fwdMomentum();
        auto const& protonMCS = track_MCS_assns_proton.at(track.key());
        if (!protonMCS.empty())
          p_proton = protonMCS[0]->fwdMomentum();
    }

    if(pdg==13) return p_muon;
    else if(pdg == 2212) return p_proton;
    else return -9999;
}

bool sbnd::LightPropagationCorrection::GetParticlePropagationTime(
    int pdg,
    double initial_momentum,
    art::FindManyP<anab::Calorimetry> track_to_calo_assns,
    art::Ptr<recob::Track> track, int PFPID ,int parentPFPID)
{
    
    bool success=false;
    double total_time=0;
    // Particle mass in GeV/c^2
    double mass = 0.0;
    double constant_dedx;

    if (std::abs(pdg) == 13)
        constant_dedx = 2.1;
    else if (std::abs(pdg) == 2212)
        constant_dedx = 4.5;


    if (std::abs(pdg) == 13)
        mass = 0.105658;
    else if (std::abs(pdg) == 2212)
        mass = 0.938272;

    // Retrieve calorimetry objects associated with the track
    std::vector<art::Ptr<anab::Calorimetry>> caloV = track_to_calo_assns.at(track.key());

    // Get the plane with the most calo points
    size_t bestPlane = 0;
    size_t maxPoints = 0;

    for (size_t j = 0; j < caloV.size(); ++j) {

        const auto& calo = *caloV[j];

        const size_t nPoints = calo.dEdx().size();

        if (nPoints > maxPoints) {
            maxPoints = nPoints;
            bestPlane = calo.PlaneID().Plane;
        }
    }

    // Loop over calorimetry objects
    for (size_t j = 0; j < caloV.size(); ++j) {

        const auto& calo = *caloV[j];

        // Only consider the collection plane
        if (calo.PlaneID().Plane != bestPlane)
            continue;

        const auto& xyz  = calo.XYZ();
        const auto& rr   = calo.ResidualRange();
        
        // Check that all vectors have the same size
        if (xyz.size() != rr.size() || rr.size()==0) {

            std::cerr << "Inconsistent calorimetry vector sizes"
                      << std::endl;
            continue;
        }

        // Create indices for the calorimetry points
        std::vector<size_t> indices(rr.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        // Sort from the beginning to the end of the track
        // (decreasing residual range)
        std::sort(indices.begin(), indices.end(),
            [&](size_t a, size_t b) {
                return rr[a] > rr[b];
            });

        // Accumulated distance along the track
        double ds = 0.0;
        double accumulated_ds = 0.0;

        double kinetic_energy;
        double momentum;
        double total_energy;
        double beta;
        double last_beta;
        double dt;

        kinetic_energy = std::sqrt(initial_momentum * initial_momentum +
                                mass * mass) - mass;

        total_energy = kinetic_energy + mass;

        momentum = initial_momentum;

        last_beta = momentum / total_energy;
        beta = last_beta;

        size_t first_idx = indices.front();
        if(parentPFPID == fNeutrinoID )
        {
            // Sum the distance from the interaction vertex to the first point in the trajectory
            double dx = xyz[first_idx].x() - fRecoVx;
            double dy = xyz[first_idx].y() - fRecoVy;
            double dz = xyz[first_idx].z() - fRecoVz;
            ds = std::sqrt(dx*dx + dy*dy + dz*dz);
            dt = ds / (last_beta * fSpeedOfLight);
            total_time += dt;
        }
        else{
            // Look at the existing spacepoints of the parent PFP and find the closest point to the trajectory.
            double minDistance=9999999;
            int minIdx=-1;
            for(size_t i=0; i<fSpacePointX.size(); i++)
            { 
                if(fSpacePointPFPID[i]==parentPFPID){
                    double dx = fSpacePointX[i] - xyz[first_idx].x();
                    double dy = fSpacePointY[i] - xyz[first_idx].y();
                    double dz = fSpacePointZ[i] - xyz[first_idx].z();
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if(distance < minDistance) {
                        minDistance = distance;
                        minIdx = i;
                    }
                }
            }
            if(minIdx!=-1)
            {
                double dx = xyz[first_idx].x() - fSpacePointX[minIdx];
                double dy = xyz[first_idx].y() - fSpacePointY[minIdx];
                double dz = xyz[first_idx].z() - fSpacePointZ[minIdx];
                ds = std::sqrt(dx*dx + dy*dy + dz*dz);
                dt = ds / (last_beta * fSpeedOfLight);
                total_time += dt + fSpacePointPropagationTime[minIdx];
            }
            else
            {
                double dx = xyz[first_idx].x() - fRecoVx;
                double dy = xyz[first_idx].y() - fRecoVy;
                double dz = xyz[first_idx].z() - fRecoVz;
                ds = std::sqrt(dx*dx + dy*dy + dz*dz);
                dt = ds / (last_beta * fSpeedOfLight);
                total_time += dt;
            }
        }

        fSpacePointX.push_back(xyz[first_idx].x());
        fSpacePointY.push_back(xyz[first_idx].y());
        fSpacePointZ.push_back(xyz[first_idx].z());
        fSpacePointPropagationTime.push_back(total_time);
        fSpacePointPFPID.push_back(PFPID);
        fSpacePointTMID.push_back(fTruthMatchedTrackID);

        // Tell whether the particle is contained

        bool isContained; 

        double track_end_x = track->End().X();
        double track_end_y = track->End().Y();
        double track_end_z = track->End().Z();

        isContained= (abs(track_end_x)<200 && abs(track_end_y)<200 && (track_end_z>0 && track_end_z<500));

        // If is contained get the energy loss through range-momentum relation
        if(isContained)
        {
            for (size_t i = 0; i < indices.size(); ++i) {

                size_t idx = indices[i];

                // First point: no previous calorimetry point
                if (i > 0) {

                    size_t prev_idx = indices[i - 1];
                    ds = std::abs(rr[idx] - rr[prev_idx]);
                    accumulated_ds+=ds;
                        momentum = tmc.GetTrackMomentum(rr[idx], pdg);        
                        kinetic_energy = std::sqrt(momentum * momentum + mass * mass) - mass;
                        total_energy = kinetic_energy + mass;
                        beta = momentum / total_energy;
                    // Avoid negative kinetic energy

                    dt = ds / (beta * fSpeedOfLight);
                    total_time += dt;
                    // Update kinetic energy after traversing the segment
                    fSpacePointX.push_back(xyz[idx].x());
                    fSpacePointY.push_back(xyz[idx].y());
                    fSpacePointZ.push_back(xyz[idx].z());
                    fSpacePointPropagationTime.push_back(total_time);
                    fSpacePointPFPID.push_back(PFPID);
                    fSpacePointTMID.push_back(fTruthMatchedTrackID);

                    /*
                    std::cout << " Space point at position X: " << xyz[idx].x()
                            << " Y: " << xyz[idx].y()
                            << " Z: " << xyz[idx].z()
                            << " t: " << total_time
                            << " momentum " << momentum
                            << " kinetic_energy " << kinetic_energy
                            << std::endl;                    
                    */       

                    
                    /*
                    std::cout << "Point " << i
                    << " ds = " << ds
                    << " dEdx = " << dedx[prev_idx]
                    << " momentum = " << momentum
                    << " beta = " << beta
                    << " dt = " << dt
                    << " accumulated time = " << total_time
                    << " acccumulated ds = " << accumulated_ds
                    << std::endl;
                    */
                }
            }
        }
        else //If not contained, assume a constant dedX through the trajectory based on pid information
        {
            for (size_t i = 0; i < indices.size(); ++i) {
                size_t idx = indices[i];

                // First point: no previous calorimetry point
                if (i > 0) {

                    size_t prev_idx = indices[i - 1];
                    ds = std::abs(rr[idx] - rr[prev_idx]);
                    
                    double dE = constant_dedx * ds / 1000;
                    kinetic_energy -= dE;
                    total_energy = kinetic_energy + mass;
                    momentum = std::sqrt(total_energy * total_energy - mass * mass);
                    beta = momentum / total_energy;
                    dt = ds / (beta * fSpeedOfLight);
                    total_time += dt;
                    
                    fSpacePointX.push_back(xyz[idx].x());
                    fSpacePointY.push_back(xyz[idx].y());
                    fSpacePointZ.push_back(xyz[idx].z());
                    fSpacePointPropagationTime.push_back(total_time);
                    fSpacePointPFPID.push_back(PFPID);
                    fSpacePointTMID.push_back(fTruthMatchedTrackID);
                    /*
                    std::cout << " Space point at position X: " << xyz[idx].x()
                            << " Y: " << xyz[idx].y()
                            << " Z: " << xyz[idx].z()
                            << " t: " << total_time
                            << " momentum " << momentum
                            << " kinetic_energy " << kinetic_energy
                            << std::endl;                    
                    */                 

                    
                }
            }
        
        }

        success=true;
    }
    return success;
}


void sbnd::LightPropagationCorrection::GetParticlePropagationTimeLite(int pfpID, std::vector<art::Ptr<recob::SpacePoint>> PFPSpacePointsVect, art::FindManyP<recob::Hit> SPToHitAssoc)
{
    for (const art::Ptr<recob::SpacePoint> &SP: PFPSpacePointsVect){
        std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());
        if (SPHit.at(0)->WireID().Plane==2){
            fSpacePointX.push_back(SP->position().X());
            fSpacePointY.push_back(SP->position().Y());
            fSpacePointZ.push_back(SP->position().Z());
            double dx = SP->position().X() - fRecoVx;
            double dy = SP->position().Y() - fRecoVy;
            double dz = SP->position().Z() - fRecoVz;
            double prop_time = std::sqrt(dx*dx + dy*dy + dz*dz) / fSpeedOfLight;
            fSpacePointPropagationTime.push_back(prop_time);
            fSpacePointPFPID.push_back(pfpID);
            fSpacePointTMID.push_back(fTruthMatchedTrackID);
        }
    }

}


void sbnd::LightPropagationCorrection::SaveTrueTrajectory(art::Event const& e)
{
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    art::Handle< std::vector<simb::MCParticle> > mclistLARG4;
    e.getByLabel(fMCModuleLabel,mclistLARG4);
    if(!mclistLARG4.isValid()){
      std::cout << " MC particles with label " << fMCModuleLabel << " not found. " << std::endl;
      throw std::exception();
    }

    std::vector<simb::MCParticle> const& mcpartVec(*mclistLARG4);
    
    for(size_t i_p=0; i_p < mcpartVec.size(); i_p++){
        const simb::MCParticle pPart = mcpartVec[i_p];
        double x = pPart.Position().X();
        double y = pPart.Position().Y();
        double z = pPart.Position().Z();
        double initial_momentum = pPart.Momentum().P();
        int pdg = pPart.PdgCode();
        if(pPart.EndT()>-10000 && pPart.EndT()<12000 && abs(x)<200 && abs(y)<200 && z<500 && z>0)
        {
            fTrueMomentum.push_back(initial_momentum);
            fTruePID.push_back(pdg);
        }


        const simb::MCTrajectory truetrack = pPart.Trajectory();
        for(size_t i_s=0; i_s < truetrack.size(); i_s++){
            //double t = pPart.Position(i_s).T();
            double x = truetrack.X(i_s);
            double y = truetrack.Y(i_s);
            double z = truetrack.Z(i_s);
            double t = truetrack.T(i_s);
            double e = truetrack.E(i_s);
            //double mom = truetrack.Momentum(i_s).P();
            if(pPart.EndT()>-10000 && pPart.EndT()<12000 && abs(x)<200 && abs(y)<200 && z<500 && z>0)
            {
                fTruthMatchedX.push_back(x);
                fTruthMatchedY.push_back(y);
                fTruthMatchedZ.push_back(z);
                fTruthMatchedT.push_back(t);
                fTruthMatchedE.push_back(e);
                fTruthMatchedID.push_back(pPart.TrackId());
            }
        }
    }
    
    // Save true energy depositions
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (e.getByLabel(fSimDepProducer, "priorSCE" ,sedHandle)){
    art::fill_ptr_vector(sedlist, sedHandle);
    }

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    for(auto& sed : sedlist ) {
        double x = sed->MidPointX();
        double y = sed->MidPointY();
        double z = sed->MidPointZ();
        double time = sed->StartT();
        art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(sed->TrackID());

        // Check if the energy deposition is from the same particle
        if( abs(time)<15000 && abs(x)<200 && abs(y)<200 && z<500 && z>0 && truth->Origin()!=2) {
            truedE.push_back(sed->Energy());
            truedE_vecX.push_back(x);
            truedE_vecY.push_back(y);
            truedE_vecZ.push_back(z);
            truedE_vecID.push_back(sed->TrackID());
            truedE_vecT.push_back(time);
            
                  // Check if the MCParticle is a cosmic in time coincidence with the G4BeamTimeWindow specified in the fhicl file

            /*
            std::cout << " Energy deposition: " << sed->Energy()
            << " at time " << sed->StartT() << " propagation time: " << (sed->StartT() - fTrueVt)
            << " and position X: " << sed->MidPointX()
            << " Y: " << sed->MidPointY()
            << " Z: " << sed->MidPointZ()
            << std::endl;            
            */

        }
    }
}

void sbnd::LightPropagationCorrection::GetParticlePropagationTimeMC(art::Event const& e, const std::vector<art::Ptr<recob::Hit> >& recoHits, int PFPID)
{
    // NEED TO REFERENCE EVERYTHING TO THE PRODUCTION TIME OF THE NEUTRINO AND THE POSITION OF THE INTERACTION VERTEX. THROUTH MCTRUTH. 

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<detinfo::DetectorClocksService> timeservice;
    auto const clockData(timeservice->DataFor(e));
    //int trkID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,recoHits,true);

    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (e.getByLabel(fSimDepProducer, "priorSCE" ,sedHandle)){
      art::fill_ptr_vector(sedlist, sedHandle);
    }

    // ============================================
    // Calculate total deposited energy per TrackID
    // ============================================

    std::map<int, double> totalEdep;

    for (auto& sed : sedlist) {
        int id = abs(sed->TrackID());
        totalEdep[id] += sed->Energy();
    }

    // ============================================
    // Store SEDs only for particles with
    // total deposited energy >= 50 MeV
    // ============================================

    for (auto& sed : sedlist) {
        int id = abs(sed->TrackID());
        // Same particle AND total Edep >= 50 MeV
        if (totalEdep[id] >= 50 && sed->Energy() > 0.1 && sed->StartT()>-10000 && sed->StartT()<12000) {

            fSpacePointX.push_back(sed->MidPointX());
            fSpacePointY.push_back(sed->MidPointY());
            fSpacePointZ.push_back(sed->MidPointZ());
            fSpacePointPropagationTime.push_back(sed->StartT() - fTrueVt);
            fSpacePointPFPID.push_back(PFPID);
            fSpacePointTMID.push_back(fTruthMatchedTrackID);
        }
    }
}


void sbnd::LightPropagationCorrection::GetMCNeutrino(art::Event const& e)
{
    if(fMCTruthModuleLabel.size()!=fMCTruthInstanceLabel.size()){
        std::cout << "MCTruthModuleLabel and MCTruthInstanceLabel vectors must have the same size..." << std::endl;
        throw std::exception();
    }

    art::Handle< std::vector<simb::MCTruth> > MCTruthListHandle;

    for (size_t s = 0; s < fMCTruthModuleLabel.size(); s++) {

        e.getByLabel(fMCTruthModuleLabel[s], fMCTruthInstanceLabel[s], MCTruthListHandle);

        if( !MCTruthListHandle.isValid() || MCTruthListHandle->empty() ) {   
            std::cout << "MCTruth with label " << fMCTruthModuleLabel[s] << " and instance " << fMCTruthInstanceLabel[s] << " not found or empty..." << std::endl;
            throw std::exception();
        }

        std::vector<art::Ptr<simb::MCTruth> > mctruth_v;
        art::fill_ptr_vector(mctruth_v, MCTruthListHandle);

        std::cout <<"Saving MCTruth from "<<fMCTruthModuleLabel[s]<<" with instance "<<fMCTruthInstanceLabel[s];
        std::cout << " with " << mctruth_v.size() << " MCTruths." << std::endl;


        for (size_t n = 0; n < mctruth_v.size(); n++) {

            art::Ptr<simb::MCTruth> evtTruth = mctruth_v[n];    
            std::cout << "  Origin: " << evtTruth->Origin() << std::endl;
            std::cout << "  We have " << evtTruth->NParticles() << " particles." << std::endl;
            std::cout << "  Mode=" << evtTruth->GetNeutrino().Mode() <<"  IntType="<<evtTruth->GetNeutrino().InteractionType();
            std::cout << "  Target=" << evtTruth->GetNeutrino().Target()<<" CCNC=" << evtTruth->GetNeutrino().CCNC()<<std::endl;

            double nu_x, nu_y, nu_z, nu_t, nu_E;

            // Loop over particles
            for (int p = 0; p < evtTruth->NParticles(); p++){

            simb::MCParticle const& par = evtTruth->GetParticle(p);

            // Only save MCTruth if the origins is specified in the fhicl list
            if( find (fMCTruthOrigin.begin(), fMCTruthOrigin.end(), evtTruth->Origin() ) != fMCTruthOrigin.end() ){

                std::cout << "    " << par.TrackId() << "  Particle PDG: " << par.PdgCode() << " E: " << par.E() << " t: " << par.T();
                std::cout << " Mother: "<<par.Mother() << " Process: "<<par.Process()<<" Status: "<<par.StatusCode()<<std::endl;
                
                // Only save vertex if the PDG is specified in the fhicl list
                if( find (fMCTruthPDG.begin(), fMCTruthPDG.end(), par.PdgCode() ) != fMCTruthPDG.end() ){
                // For BNB neutinos
                if(par.StatusCode()==0 && evtTruth->Origin()==1 && abs(par.Vx())<200 && abs(par.Vy())<200 && par.Vz()>0 && par.Vz()<500){
                    fTrueVx=par.Vx();
                    nu_x = par.Vx();
                    fTrueVy=par.Vy();
                    nu_y = par.Vy();
                    fTrueVz=par.Vz();
                    nu_z = par.Vz();
                    fTrueVt=par.T();
                    nu_t = par.T();
                    nu_E = par.E();
                }
                }

            }

            }

            std::cout << "Vertex: " << nu_x << " " << nu_y << " " << nu_z << " T: " << nu_t << " E: " << nu_E << std::endl;

        }
    }
}


void sbnd::LightPropagationCorrection::GetTruthMatchedID(art::Event const& e, const std::vector<art::Ptr<recob::Hit> >& recoHits)
{
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<detinfo::DetectorClocksService> timeservice;
    auto const clockData(timeservice->DataFor(e));
    fTruthMatchedTrackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,recoHits,true);
}



