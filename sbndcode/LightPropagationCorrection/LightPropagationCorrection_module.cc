#include "LightPropagationCorrection_module.hh"

#include <stdexcept>

namespace sbnd {
    class LightPropagationCorrection;
}

/**
 * @brief Construct the light-propagation correction producer.
 *
 * Reads the producer configuration, caches optical detector geometry and
 * PDS channel-to-box information, configures the flash reconstruction tools
 * for both TPCs, and declares the output products.
 *
 * @param p FHiCL parameter set containing input labels, reconstruction
 *          settings, detector parameters, selection thresholds and debug options.
 */
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
    fFlashMatchingTool( p.get<std::string>("FlashMatchingTool") ),
    fDoPIDCorrection( p.get<bool>("DoPIDCorrection") ),
    fUseMCVIS( p.get<bool>("UseMCVIS") ),
    fVGroupVIS( p.get<double>("VGroupVIS") ),
    fVGroupVUV( p.get<double>("VGroupVUV") ),
    fMIPdEdX( p.get<double>("MIPdEdX") ),
    fHIPdEdX( p.get<double>("HIPdEdX") ),
    fMaxRRBeta( p.get<double>("MaxRRBeta") ),
    fCherenkovRefractiveIndex( p.get<double>("CherenkovRefractiveIndex")),
    fCherenkovCosineTolerance( p.get<double>("CherenkovCosineTolerance") ),
    fCherenkovPMTRadius(p.get<double>("CherenkovPMTRadius")),
    fCherenkovPMTRings(p.get<unsigned int>("CherenkovPMTRings")),
    fCherenkovPMTAzimuthSteps( p.get<unsigned int>("CherenkovPMTAzimuthSteps")),
    fMaxCoatedUncoatedTimeDiff( p.get<double>("MaxCoatedUncoatedTimeDiff") ),
    fMinHitPE( p.get<double>("MinHitPE") ),
    fPreWindow( p.get<double>("PreWindow") ),
    fPostWindow( p.get<double>("PostWindow") ),
    fFVx( p.get<std::vector<double>>("FVx") ),
    fFVy( p.get<std::vector<double>>("FVy") ),
    fFVz( p.get<std::vector<double>>("FVz") ),
    fNuScoreThreshold( p.get<double>("NuScoreThreshold") ),
    fFMScoreThreshold( p.get<double>("FMScoreThreshold") ),
    fTrackScoreThreshold( p.get<double>("TrackScoreThreshold") ),
    fSaveCorrectionTree( p.get<bool>("SaveCorrectionTree") ),
    fDebug( p.get<bool>("Debug", false) )
    //
    // More initializers here.
{
    fNOpChannels = fWireReadout.NOpChannels();

    for (unsigned int opch=0; opch<fNOpChannels; opch++){
        auto pdCenter = fWireReadout.OpDetGeoFromOpChannel(opch).GetCenter();
        fOpDetID.push_back(opch);
        fOpDetX.push_back(pdCenter.X());
        fOpDetY.push_back(pdCenter.Y());
        fOpDetZ.push_back(pdCenter.Z());
        fChannelBoxMap[opch] = fPDSMap.pdBox(opch);
        if(fPDSMap.pdType(opch)=="pmt_uncoated")
            fUncoatedChannelMap[fPDSMap.pdBox(opch)] = opch;
    }

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

/**
 * @brief Build corrected optical-flash timing products for one event.
 *
 * Reads the reconstructed slices and their associated particles, tracks,
 * space points, optical flashes and flash-matching information. For each
 * selected neutrino slice, it estimates particle propagation times, computes
 * the scintillation and Cherenkov light-propagation corrections, corrects the
 * associated OpHits, and reconstructs the corrected flash time.
 *
 * The corrected timing object is associated with both its source Slice and
 * original OpFlash. Slices failing the configured selections, or for which a
 * valid corrected flash cannot be reconstructed, are skipped.
 *
 * @param e art event containing the input products and receiving the outputs.
 */
void sbnd::LightPropagationCorrection::produce(art::Event & e)
{
    fEvent = e.id().event();
    fRun = e.id().run();
    fSubrun = e.id().subRun();
    _flashgeo->InitializeFlashGeoAlgo();

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
    art::FindManyP<recob::OpHit> flashToOpHitAssns_tpc0(opflashListHandle_tpc0, e, fOpFlashLabel_tpc0);
    art::FindManyP<recob::OpHit> flashToOpHitAssns_tpc1(opflashListHandle_tpc1, e, fOpFlashLabel_tpc1);


    // PFP Metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_to_metadata(pfpHandle, e, fReco2Label);

    //Get the SP Hit assns
    art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
    std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
    e.getByLabel(fSpacePointLabel, eventSpacePoints);
    art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);
    art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);

    // --- Store candidate slices
    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);

    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfpVect;
    // --- Get the candidate slices
    for (size_t ix=0; ix<sliceVect.size(); ix++){
        ResetSliceInfo();
        // --- Get the slice
        auto & slice = sliceVect[ix];
        // Now I need to get all the hits associated to this flash and get the timing for all of them
        // Get the slices PFPs
        double _sliceMaxNuScore = -9999.;
        pfpVect = slice_pfp_assns.at(slice.key());
        fNeutrinoID = -1;
        // Get the neutrino ID and vertex position
        for (const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp.key());
            if (pfp->IsPrimary() &&( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ){
                for (auto const pfpMeta : pfpMetaVec) {
                    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                    if (propertiesMap.count("NuScore")) _fNuScore = propertiesMap.at("NuScore");
                    if (_fNuScore>_sliceMaxNuScore) _sliceMaxNuScore = _fNuScore;
                }
                std::vector< art::Ptr<recob::Vertex> > vertexVec = pfp_vertex_assns.at(pfp.key());
                for (const art::Ptr<recob::Vertex> &ver : vertexVec){
                    geo::Point_t xyz_vertex = ver->position();
                    fRecoVx= xyz_vertex.X();
                    fRecoVy= xyz_vertex.Y();
                    fRecoVz= xyz_vertex.Z();
                }
                fNeutrinoID = pfp->Self();
            }
        }

        // Get informaiton on the remaining PFPs
        for (const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            if (fNeutrinoID == static_cast<int>(pfp->Self())) continue;
            bool tracksuccess = false;
            // Check if it is a clear cosmic
            bool isClearCosmic = false;
            double trackScore = -1;
            const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp.key());
            // Check if the PFP is a clear cosmic
            for (auto const pfpMeta : pfpMetaVec) {
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                if (propertiesMap.count("IsClearCosmic")){
                    isClearCosmic=true;
                }
                if (propertiesMap.count("TrackScore")){
                    trackScore=propertiesMap.at("TrackScore");
                }
            }
            if (isClearCosmic) continue;
            // Neutrino daughter particles
            std::vector<art::Ptr<recob::SpacePoint>> PFPSpacePointsVect = pfp_sp_assns.at(pfp.key());
            if (fDoPIDCorrection && trackScore>fTrackScoreThreshold)
            {
                // If it's a track then go for PID and momentum calculation to get particle propagation time
                std::vector<art::Ptr<recob::Track>> track_v = pfp_track_assns.at(pfp.key());
                if (track_v.size()>1)
                    throw art::Exception(art::errors::LogicError) << "Multiple tracks associated to a PFP. This is not expected.";
                if (track_v.size() == 1) {
                    int pdg = GetTrackPID(track_to_pid_assns, track_v[0]);
                    if (pdg != -1) {
                        double initial_momentum = GetTrackMomentum(pdg, track_rangeP_assns_muon, track_rangeP_assns_proton, track_MCS_assns_muon, track_MCS_assns_proton ,  track_v[0]);
                        if (initial_momentum > 0) tracksuccess = GetParticlePropagationTime( pdg, initial_momentum, track_to_calo_assns, track_v[0], pfp->Self(), pfp->Parent() );
                    }
                }
            }
            if (!tracksuccess){
                // If it's not a track just assume speed of light propagation
                GetParticlePropagationTimeLite(pfp->Self(), PFPSpacePointsVect, SPToHitAssoc);
            }
        }

        if (_sliceMaxNuScore<fNuScoreThreshold){
            ResetSliceInfo();
            continue; // Skip to the next slice if the nu score is below threshold
        }

        if (fRecoVx == -99999. || fRecoVy == -99999. || fRecoVz == -99999.) {
            ResetSliceInfo();
            continue;
        }

        this->GetPropagationTimeCorrectionPerChannelScintillation();
        this->GetPropagationTimeCorrectionPerChannelCherenkov();

        // Get all the OpT0 objects associated to the slice
        std::vector<art::Ptr<recob::OpFlash>> flashFM;
        if (fFlashMatchingTool == "OpT0Finder" ){
            const std::vector< art::Ptr<sbn::OpT0Finder> > slcOpT0Finder = slice_opt0finder_assns.at( slice.key() );
            if (slcOpT0Finder.size() == 0) continue;
            size_t OpT0Idx = HighestOpT0ScoreIdx(slcOpT0Finder);
            // Get the flash OpT0 association
            art::FindManyP<recob::OpFlash> opflash_opt0finder_assns(opt0Handle, e, fOpT0FinderModuleLabel);
            flashFM = opflash_opt0finder_assns.at( slcOpT0Finder[OpT0Idx].key() );
            if (flashFM.size() > 1){
                throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same OpT0Finder object. This is not expected.";
            }
            if (flashFM.size() == 0){
                throw art::Exception(art::errors::LogicError) << "There are no OpFlash objects associated to the same OpT0Finder object. This is not expected.";
            }
            _fFMScore = slcOpT0Finder[OpT0Idx]->score;
            if (_fFMScore < fFMScoreThreshold){
                ResetSliceInfo();
                continue;
            }
        }
        else if (fFlashMatchingTool == "BarycenterFM")
        {
            const std::vector< art::Ptr<sbn::TPCPMTBarycenterMatch> > slcTPCPMTBarycenter = slice_tpcpmtbarycentermatching_assns.at( slice.key() );
            if (slcTPCPMTBarycenter.size() == 0) continue;
            art::FindManyP<recob::OpFlash> opflash_tpcpmtbarycenter_assns(tpcpmtbarycenterfmHandle, e, fTPCPMTBarycenterFMModuleLabel);
            size_t BFMIdx = HighestBFMScoreIdx(slcTPCPMTBarycenter);
            flashFM = opflash_tpcpmtbarycenter_assns.at( slcTPCPMTBarycenter[BFMIdx].key() );
            if (flashFM.size() > 1){
                throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same TPCPMTBarycenterFM object. This is not expected.";
            }
            if (flashFM.size() == 0){
                throw art::Exception(art::errors::LogicError) << "There are no OpFlash objects associated to the same TPCPMTBarycenterFM object. This is not expected.";
            }
            _fFMScore = slcTPCPMTBarycenter[BFMIdx]->score;
            if (_fFMScore < fFMScoreThreshold){
                ResetSliceInfo();
                continue;
            }
        }
        else throw art::Exception(art::errors::LogicError) << " Flash matching tool " <<  fFlashMatchingTool << " not supported ." << std::endl;

        sbn::CorrectedOpFlashTiming correctedOpFlashTiming;
        bool CorrectOpFlashResult = CorrectOpFlash(flashFM[0], correctedOpFlashTiming, flashToOpHitAssns_tpc0, flashToOpHitAssns_tpc1);
        if (!CorrectOpFlashResult) continue; // Skip to the next slice if the flash is not valid
        correctedOpFlashTimes->emplace_back(std::move(correctedOpFlashTiming));
        art::Ptr<sbn::CorrectedOpFlashTiming> newCorrectedOpFlashTimingPtr = make_correctedopflashtime_ptr(correctedOpFlashTimes->size()-1);
        newCorrectedOpFlashTimingSliceAssn->addSingle(slice, newCorrectedOpFlashTimingPtr);
        newCorrectedOpFlashTimingOpFlashAssn->addSingle(flashFM[0], newCorrectedOpFlashTimingPtr);
    }

    if (fSaveCorrectionTree) fTree->Fill();
    ResetEventVars();

    e.put(std::move(correctedOpFlashTimes));
    e.put(std::move(newCorrectedOpFlashTimingSliceAssn));
    e.put(std::move(newCorrectedOpFlashTimingOpFlashAssn));
    return;
}

/**
 * @brief Initialize optional diagnostic output.
 *
 * When SaveCorrectionTree is enabled, creates the diagnostic TTree and
 * registers branches used to store event, slice, flash, space-point and
 * OpHit timing information.
 */
void sbnd::LightPropagationCorrection::beginJob()
{
    if (fSaveCorrectionTree)
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
        fTree->Branch("OpHitOldTime", &fOpHitOldTime);
        fTree->Branch("OpHitNewTime", &fOpHitNewTime);
        fTree->Branch("OpHitPE", &fOpHitPE);
        fTree->Branch("OpHitOpCh", &fOpHitOpCh);
    }
}

/**
 * @brief Perform end-of-job finalization.
 *
 * No additional finalization is currently required by this producer.
 */
void sbnd::LightPropagationCorrection::endJob()
{
}

/**
 * @brief Reset variables accumulated over the current event.
 *
 * Clears event-level diagnostic vectors and resets the stored event
 * identifiers. Diagnostic variables are only reset when the correction tree
 * is enabled.
 */
void sbnd::LightPropagationCorrection::ResetEventVars()
{

    if (fSaveCorrectionTree)
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
        fSliceSPID.clear();
    }
}

/**
 * @brief Find the highest-scoring OpT0Finder match in a slice.
 *
 * @param slcFM OpT0Finder objects associated with the slice.
 * @return Index of the object with the largest score. Returns zero if the
 *         input vector is empty.
 */
size_t sbnd::LightPropagationCorrection::HighestOpT0ScoreIdx(std::vector< art::Ptr<sbn::OpT0Finder> > const & slcFM)
{
    // Gets the idx of the OpT0 object with the highest score
    double highestOpT0Score = -99999.0; // Initialize to a negative value
    size_t highestIdx = 0;
    for (size_t jx=0; jx<slcFM.size(); jx++){
        if (slcFM[jx]->score > highestOpT0Score){
            highestOpT0Score = slcFM[jx]->score;
            highestIdx = jx;
        }
    }
    return highestIdx;
}

/**
 * @brief Find the highest-scoring TPC-PMT barycenter match in a slice.
 *
 * @param slcFM TPCPMTBarycenterMatch objects associated with the slice.
 * @return Index of the object with the largest score. Returns zero if the
 *         input vector is empty.
 */
size_t sbnd::LightPropagationCorrection::HighestBFMScoreIdx(std::vector< art::Ptr<sbn::TPCPMTBarycenterMatch> > const & slcFM)
{
    // Gets the idx of the OpT0 object with the highest score
    double highesBFM0Score = -99999.0; // Initialize to a negative value
    size_t highestIdx = 0;
    for (size_t jx=0; jx<slcFM.size(); jx++){
        if (slcFM[jx]->score > highesBFM0Score){
            highesBFM0Score = slcFM[jx]->score;
            highestIdx = jx;
        }
    }
    return highestIdx;
}

/**
 * @brief Reset quantities that are specific to the current reconstructed slice.
 *
 * Clears the particle trajectory points, resets the reconstructed neutrino
 * vertex and NuScore, and initializes both per-channel light-time correction
 * vectors to zero.
 */
void sbnd::LightPropagationCorrection::ResetSliceInfo()
{
    _fNuScore=-99999.0;
    fRecoVx = -99999.0;
    fRecoVy = -99999.0;
    fRecoVz = -99999.0;
    fParticleSpacePoints.clear();
    fTimeCorrectionPerChannelScintillation.assign(fNOpChannels, 0.0); // Reset the time correction vector for each channel
    fTimeCorrectionPerChannelCherenkov.assign(fNOpChannels, 0.0); // Reset the time correction vector for each channel
}

/**
 * @brief Calculate the scintillation-light propagation correction per optical channel.
 *
 * For each optical detector, evaluates the total time for light to propagate
 * from every stored particle space point to that detector. The propagation
 * model depends on the detector type and on whether the MCVIS geometry
 * approximation is enabled. The minimum particle-plus-photon propagation time
 * is selected, and its negative is stored as the channel correction.
 *
 * Distances are in cm and propagation velocities are expected in cm/ns;
 * consequently, the resulting correction is in ns.
 */
void sbnd::LightPropagationCorrection::GetPropagationTimeCorrectionPerChannelScintillation()
{
    // Initialize one correction value per optical channel (time in ns).
    fTimeCorrectionPerChannelScintillation.assign(fNOpChannels, 0.0); // Initialize with zero or any default value
    // Implementation
    for (size_t opdet = 0; opdet < fOpDetID.size(); opdet++) {
        double _opDetX = fOpDetX[opdet];
        double _opDetY = fOpDetY[opdet];
        double _opDetZ = fOpDetZ[opdet];
        double minPropTime = 999999999.;
        bool foundSP = false;
        for (size_t sp=0; sp<fParticleSpacePoints.size(); sp++)
        {
            bool isInSameTPC = (fParticleSpacePoints[sp].x * _opDetX) > 0;
            if(!isInSameTPC) continue; // Skip points not in the same TPC
            double dx = fParticleSpacePoints[sp].x - _opDetX;
            double dy = fParticleSpacePoints[sp].y - _opDetY;
            double dz = fParticleSpacePoints[sp].z - _opDetZ;
            double distanceToOpDet = std::sqrt(dx*dx + dy*dy + dz*dz);
            double spToCathode;
            double cathodeToOpDet;
            if (fUseMCVIS){
                spToCathode = std::sqrt(fParticleSpacePoints[sp].x*fParticleSpacePoints[sp].x); // Distance from space point to cathode in cm
                cathodeToOpDet = std::sqrt(_opDetX*_opDetX + (dy)*(dy) + (dz)*(dz)); // Distance from cathode to OpDet in cm
            }
            else{
                spToCathode = std::sqrt( fParticleSpacePoints[sp].x*fParticleSpacePoints[sp].x + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from space point to cathode in cm
                cathodeToOpDet = std::sqrt(_opDetX*_opDetX + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from cathode to OpDet in cm
            }

            // VIS model: wavelength shifting at the cathode followed by visible-light propagation to the detector.
            double lightPropTimeVIS = spToCathode/fVGroupVUV + cathodeToOpDet/fVGroupVIS; // Speed
            // Direct VUV propagation from the emission point to the optical detector.
            double lightPropTimeVUV = distanceToOpDet / fVGroupVUV; // Speed of light in cm/ns for VUV
            double lightPropTime = 0;
            const std::string pdType = fPDSMap.pdType(opdet);
            if (pdType=="pmt_coated" || pdType=="xarapuca_vuv")
                lightPropTime = std::min(lightPropTimeVIS, lightPropTimeVUV);
            else if (pdType=="pmt_uncoated" || pdType=="xarapuca_vis")
            lightPropTime = lightPropTimeVIS;
            else
            throw std::runtime_error("LightPropagationCorrection: unexpected pdType '" + pdType + "' for opdet " + std::to_string(opdet));

            //float fastest_partPropTime = std::sqrt((fSpacePointX[sp]-fRecoVx)*(fSpacePointX[sp]-fRecoVx) + (fSpacePointY[sp]-fRecoVy)*(fSpacePointY[sp]-fRecoVy) + (fSpacePointZ[sp]-fRecoVz)*(fSpacePointZ[sp]-fRecoVz))/fSpeedOfLight;
            double partPropTime = fParticleSpacePoints[sp].propagationTime;
            double PropTime = lightPropTime + partPropTime;
            if (PropTime < minPropTime) {
                minPropTime = PropTime;
                foundSP = true;
            }
        }
        if (!foundSP) {
            fTimeCorrectionPerChannelScintillation[opdet] = 0.0;
            continue;
        }
        fTimeCorrectionPerChannelScintillation[opdet] = -minPropTime;
    }
}


/**
 * @brief Apply the scintillation correction and, where selected, the Cherenkov correction to OpHits.
 *
 * First applies the per-channel scintillation correction to each input OpHit.
 * It then reconstructs a preliminary flash, selects the first qualifying hit
 * in each channel around that flash time, and compares coated and uncoated
 * PMT timing within each PDS box. For boxes passing the configured timing
 * difference requirement, the Cherenkov correction is used when its magnitude
 * is larger than the scintillation correction.
 *
 * @param OldOpHitList Original OpHits associated with the flash.
 * @param newOpHitList Output vector containing the corrected OpHits.
 */
void sbnd::LightPropagationCorrection::CorrectOpHitTime(std::vector<art::Ptr<recob::OpHit>> const& OldOpHitList, std::vector<recob::OpHit>& newOpHitList)
{
    newOpHitList.clear();

    // Apply the scintillation correction first; the Cherenkov decision is made below at box level.
    for (auto const& hit : OldOpHitList) {
        const size_t opCh = hit->OpChannel();
        if (opCh >= fNOpChannels) continue;

        const double channelCorrection = fTimeCorrectionPerChannelScintillation[opCh];
        const double newPeakTime = hit->StartTime() + hit->RiseTime() + channelCorrection / 1000.;
        const double newPeakTimeAbs = hit->PeakTimeAbs() + channelCorrection / 1000.;
        const double newStartTime = hit->StartTime() + channelCorrection / 1000.;

        newOpHitList.emplace_back(opCh, newPeakTime, newPeakTimeAbs, newStartTime,
                                  hit->RiseTime(), hit->Frame(), hit->Width(),
                                  hit->Area(), hit->Amplitude(), hit->PE(), 0.0);
    }

    // Reconstruct a preliminary flash to define the time window used for the box comparison.
    ::lightana::LiteOpHitArray_t ophits;
    FillLiteOpHit(newOpHitList, ophits);

    const auto flash_v = _mgr.RecoFlash(ophits);
    if (flash_v.empty()) return;
    double flasht0 = 0.0;
    for (const auto& lflash : flash_v) flasht0 = lflash.time;

    // Get the first OpHit of each channel around the flash t0
    std::vector<recob::OpHit> first_ophits;
    FillFirstOpHit(flasht0, newOpHitList, first_ophits);
    GetBoxTimeDiffMap(first_ophits);

    // Apply the Cherenkov correction only in boxes whose coated/uncoated timing difference passes the configured cut.
    for (auto& hit : newOpHitList) {
        const size_t opCh = hit.OpChannel();
        if (opCh >= fNOpChannels) continue;

        const size_t box = fChannelBoxMap.at(opCh);
        const double scintCorrection = fTimeCorrectionPerChannelScintillation[opCh];
        const double cherenkovCorrection = fTimeCorrectionPerChannelCherenkov[opCh];

        auto it = fBoxTimeDiffMap.find(box);
        if (it == fBoxTimeDiffMap.end()) continue;
        if (it->second > fMaxCoatedUncoatedTimeDiff && (std::abs(cherenkovCorrection) > std::abs(scintCorrection))) {
            const double correction = (cherenkovCorrection - scintCorrection) / 1000.;

            const double newPeakTime = hit.PeakTime() + correction;
            const double newPeakTimeAbs = hit.PeakTimeAbs() + correction;
            const double newStartTime = hit.StartTime() + correction;

            hit = recob::OpHit(opCh, newPeakTime, newPeakTimeAbs, newStartTime,
                               hit.RiseTime(), hit.Frame(), hit.Width(),
                               hit.Area(), hit.Amplitude(), hit.PE(), 0.0);
        }
    }
}


/**
 * @brief Select the earliest qualifying OpHit in each optical channel.
 *
 * A hit is accepted when its peak time lies inside the configured window
 * around the flash time and its PE is above MinHitPE. At most one hit per
 * channel is copied to the output.
 *
 * @param flash_time Reference flash time, in the same time units as OpHit peak times.
 * @param all_opHits Input OpHits to examine.
 * @param first_ophits Output vector containing the selected hits.
 */
void sbnd::LightPropagationCorrection::FillFirstOpHit(double flash_time, std::vector<recob::OpHit> const& all_opHits, std::vector<recob::OpHit>& first_ophits)
{
    first_ophits.clear();
    std::vector<double> firstTimePerChannel(fNOpChannels, std::numeric_limits<double>::max());
    std::vector<bool> channelHasHit(fNOpChannels, false);

    for (auto const& hit : all_opHits) {
        if (hit.PeakTime() < flash_time + fPostWindow && hit.PeakTime() > flash_time - fPreWindow && hit.PE() > fMinHitPE){
            size_t channel = hit.OpChannel();
            if (channel >= fNOpChannels) continue;
            if (!channelHasHit[channel] || hit.PeakTime() < firstTimePerChannel[channel]) {
                firstTimePerChannel[channel] = hit.PeakTime();
                channelHasHit[channel] = true;
            }
        }
    }

    for (auto const& hit : all_opHits) {
        size_t channel = hit.OpChannel();
        if (channel < fNOpChannels && channelHasHit[channel] &&
            hit.PeakTime() == firstTimePerChannel[channel]) {
            first_ophits.push_back(hit);
            channelHasHit[channel] = false;
        }
    }
}

/**
 * @brief Convert recob::OpHit objects into the lightweight format used by lightana.
 *
 * @param OpHitList Input OpHits.
 * @param LiteOpHitList Output lightweight OpHit array; cleared before filling.
 */
void sbnd::LightPropagationCorrection::FillLiteOpHit(std::vector<recob::OpHit> const& OpHitList, std::vector<::lightana::LiteOpHit_t>& LiteOpHitList)
{
    LiteOpHitList.clear();
    for (auto const& oph : OpHitList) {
        ::lightana::LiteOpHit_t loph;
        loph.peak_time = oph.StartTime()+oph.RiseTime();
        loph.pe = oph.PE();
        loph.channel = oph.OpChannel();
        LiteOpHitList.emplace_back(std::move(loph));
    }
}

/**
 * @brief Retrieve the lightweight OpHits associated with a lightweight flash.
 *
 * Uses the hit indices stored in the flash to copy the corresponding entries
 * from the supplied hit array.
 *
 * @param lite_flash Flash whose associated-hit indices are used.
 * @param lite_hits_v Array containing the lightweight OpHits.
 * @return Lightweight OpHit array associated with the flash.
 */
::lightana::LiteOpHitArray_t sbnd::LightPropagationCorrection::GetAssociatedLiteHits(::lightana::LiteOpFlash_t lite_flash, ::lightana::LiteOpHitArray_t lite_hits_v)
{
    ::lightana::LiteOpHitArray_t flash_hits_v;

    for(auto const& hitidx : lite_flash.asshit_idx) {
        flash_hits_v.emplace_back(std::move(lite_hits_v.at(hitidx)));
    }

    return flash_hits_v;
}


/**
 * @brief Calculate the largest coated-to-uncoated timing difference in each PDS box.
 *
 * Separates the selected first hits into uncoated and coated channels using
 * the channel-to-box map. For each box with an uncoated hit and at least one
 * coated hit, stores the largest absolute peak-time difference.
 *
 * @param first_ophits Selected first OpHits, normally one per channel.
 */
void sbnd::LightPropagationCorrection::GetBoxTimeDiffMap( std::vector<recob::OpHit>& first_ophits)
{
    fBoxTimeDiffMap.clear();
    std::map<size_t, double> uncoatedTimePerBox;
    std::map<size_t, std::vector<double>> coatedTimesPerBox;
    size_t channel;
    size_t box;

    for (const auto& opHit : first_ophits) {
        channel = opHit.OpChannel();
        box = fChannelBoxMap.at(channel);
        double peakTime = opHit.StartTime() + opHit.RiseTime();
        auto it = fUncoatedChannelMap.find(box);
        if (it == fUncoatedChannelMap.end()) continue;
        if (channel == it->second)
            uncoatedTimePerBox[box] = peakTime;
        else
            coatedTimesPerBox[box].push_back(peakTime);
    }

    for (const auto& [box, uncoatedTime] : uncoatedTimePerBox) {
        if (!coatedTimesPerBox.count(box)) continue;
        for (double coatedTime : coatedTimesPerBox.at(box)) {
            double timeDiff = std::abs(coatedTime - uncoatedTime);
            if (!fBoxTimeDiffMap.count(box) || timeDiff > fBoxTimeDiffMap[box])
                fBoxTimeDiffMap[box] = timeDiff;
        }
    }
}



/**
 * @brief Store flash- and OpHit-level correction diagnostics.
 *
 * Records the original and corrected flash times, flash position and PE,
 * reconstructed neutrino vertex, and (when Debug is enabled) particle
 * space-point and OpHit information.
 *
 * @param newFlashTime Reconstructed flash time after applying the corrections.
 * @param flash Original OpFlash associated with the slice.
 * @param oldOpHitList Original OpHits associated with the flash.
 * @param newOpHitList Corrected OpHits used for the new flash reconstruction.
 */
void sbnd::LightPropagationCorrection::FillCorrectionTree(double & newFlashTime, recob::OpFlash const& flash, std::vector<recob::OpHit> const& oldOpHitList, std::vector<recob::OpHit> const& newOpHitList){
    fOpHitOldTime.push_back({});
    fOpHitNewTime.push_back({});
    fOpHitPE.push_back({});
    fOpHitOpCh.push_back({});

    if (fDebug)
    {
        for (size_t i=0; i<fParticleSpacePoints.size(); i++){
            fSliceSPX.push_back(fParticleSpacePoints[i].x);
            fSliceSPY.push_back(fParticleSpacePoints[i].y);
            fSliceSPZ.push_back(fParticleSpacePoints[i].z);
            fSliceSPT.push_back(fParticleSpacePoints[i].propagationTime);
            fSliceSPID.push_back(fParticleSpacePoints[i].pfpID);
        }

        for (size_t i=0; i<oldOpHitList.size(); i++){
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


/**
 * @brief Correct the timing of an OpFlash and fill its corrected timing object.
 *
 * Selects the OpHit association and flash-reconstruction manager corresponding
 * to the TPC containing the flash. It corrects the associated OpHits,
 * reconstructs a flash from them, refines its time using the configured
 * FlashT0 tool, and fills the output timing quantities.
 *
 * @param flash Original OpFlash to correct.
 * @param correctedOpFlashTiming Output object containing original and corrected timing information.
 * @param flashToOpHitAssns_tpc0 OpFlash-to-OpHit associations for TPC 0.
 * @param flashToOpHitAssns_tpc1 OpFlash-to-OpHit associations for TPC 1.
 * @return true if a single corrected flash is reconstructed; false otherwise.
 */
bool sbnd::LightPropagationCorrection::CorrectOpFlash(art::Ptr<recob::OpFlash> const& flash, sbn::CorrectedOpFlashTiming &correctedOpFlashTiming, const art::FindManyP<recob::OpHit>& flashToOpHitAssns_tpc0, const art::FindManyP<recob::OpHit>& flashToOpHitAssns_tpc1)
{
    // Get the ophits associated to the flash
    std::vector<art::Ptr<recob::OpHit>> ophitlist;
    if (abs(flash->XCenter()) > 200) return false;
    if (flash->XCenter()<0)
    {
        ophitlist = flashToOpHitAssns_tpc0.at(flash.key());
        _mgr = _mgr_tpc0; // Use the TPC 0 flash finder manager
    }
    else
    {
        ophitlist = flashToOpHitAssns_tpc1.at(flash.key());
        _mgr = _mgr_tpc1; // Use the TPC 1 flash finder manager
    }
    std::vector<recob::OpHit> newOpHitList;
    std::vector<recob::OpHit> oldOpHitList;
    for (const auto& ophit : ophitlist) {
        oldOpHitList.push_back(*ophit);
    }
    // Get the list of the corrected OpHits
    this->CorrectOpHitTime(ophitlist, newOpHitList);
    // Create the list of ophit lite to be used in the flash finder
    ::lightana::LiteOpHitArray_t ophits;
    this->FillLiteOpHit(newOpHitList, ophits);
    // Create the flash manager
    auto const flash_v = _mgr.RecoFlash(ophits);
    if (flash_v.size() != 1) return false; // Return false if no flash is reconstructed
    double originalFlashTime = flash->Time();
    double newFlashTime = 0.0;
    for (const auto& lflash :  flash_v) {
        // Get Flash Barycenter
        double Ycenter, Zcenter, Ywidth, Zwidth;
        _flashgeo->GetFlashLocation(lflash.channel_pe, Ycenter, Zcenter, Ywidth, Zwidth);
        // Get flasht0
        double flasht0 = lflash.time;
        // Refine t0 calculation
        flasht0 = _flasht0calculator->GetFlashT0(lflash.time, GetAssociatedLiteHits(lflash, ophits));
        newFlashTime = flasht0;
        correctedOpFlashTiming.OpFlashT0 = originalFlashTime;
        correctedOpFlashTiming.NuToFLight = (Zcenter/c_light)/1000;
        correctedOpFlashTiming.NuToFCharge = (fRecoVz/c_light)/1000;
        correctedOpFlashTiming.OpFlashT0Corrected = newFlashTime;
    }
    if (fSaveCorrectionTree){
        this->FillCorrectionTree(newFlashTime, *flash, oldOpHitList, newOpHitList);
    }
    return true;
}


/**
 * @brief Check whether both endpoints of a track lie inside the configured fiducial volume.
 *
 * The fiducial volume is defined by the lower and upper bounds in fFVx,
 * fFVy and fFVz.
 *
 * @param start Track start position.
 * @param end Track end position.
 * @return true only if both endpoints are strictly inside all three coordinate bounds.
 */
bool sbnd::LightPropagationCorrection::IsInFV(geo::Point_t const& start, geo::Point_t const& end){
    bool startInFV = ( (start.x() > fFVx[0]) && (start.x() < fFVx[1]) && (start.y() > fFVy[0]) && (start.y() < fFVy[1]) && (start.z() > fFVz[0]) && (start.z() < fFVz[1]));
    bool endInFV = ( (end.x() > fFVx[0]) && (end.x() < fFVx[1]) && (end.y() > fFVy[0]) && (end.y() < fFVy[1]) && (end.z() > fFVz[0]) && (end.z() < fFVz[1]));

    return startInFV && endInFV;
}


/**
 * @brief Identify a track as a muon or proton using collection-plane PID Chi2 values.
 *
 * Reads the ParticleID algorithm scores associated with the track, keeps
 * Chi2 scores from plane 2, and applies the configured hard-coded muon/proton
 * selection criteria.
 *
 * @param track_to_pid_assns Track-to-ParticleID associations.
 * @param track Track whose PID hypothesis is evaluated.
 * @return PDG code 13 for a muon, 2212 for a proton, or -1 if no hypothesis
 *         satisfies the selection or the required scores are unavailable.
 */
int sbnd::LightPropagationCorrection::GetTrackPID(art::FindManyP<anab::ParticleID>  const& track_to_pid_assns, art::Ptr<recob::Track> const& track )
{
    std::vector<art::Ptr<anab::ParticleID>> pidV = track_to_pid_assns.at(track.key());

    double muonChi2   = -99999.;
    double protonChi2 = -99999.;

    for (size_t j = 0; j < pidV.size(); ++j) {

        const auto& pid = *pidV[j];
        int plane = pid.PlaneID().Plane;
        // Use only collection plane
        if (plane!=2) continue;

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

    if (protonChi2== -99999. || muonChi2== -99999. )
        return -1;

    if (muonChi2<25 && protonChi2>95) return 13;
    else if ( protonChi2<=95) return 2212;
    else return -1;
}



/**
 * @brief Retrieve the initial momentum estimate for a track under a given PID hypothesis.
 *
 * Uses the range-based momentum estimate for tracks whose two endpoints are
 * inside the fiducial volume. Otherwise, uses the forward MCS momentum
 * estimate. The requested particle hypothesis determines which association
 * is read.
 *
 * @param pdg Particle PDG code (muon or proton).
 * @param track_rangeP_assns_muon Track-to-muon RangeP associations.
 * @param track_rangeP_assns_proton Track-to-proton RangeP associations.
 * @param track_MCS_assns_muon Track-to-muon MCS-fit associations.
 * @param track_MCS_assns_proton Track-to-proton MCS-fit associations.
 * @param track Track for which the momentum is requested.
 * @return Momentum in GeV/c, or a negative sentinel value if unavailable or unsupported.
 */
double sbnd::LightPropagationCorrection::GetTrackMomentum(int pdg, art::FindManyP<sbn::RangeP> const& track_rangeP_assns_muon,  art::FindManyP<sbn::RangeP> const& track_rangeP_assns_proton, art::FindManyP<recob::MCSFitResult> const& track_MCS_assns_muon , art::FindManyP<recob::MCSFitResult> const& track_MCS_assns_proton  , art::Ptr<recob::Track> const& track )
{
    double p_muon = -99999;
    double p_proton = -99999;

    geo::Point_t const& start = track->Vertex();
    geo::Point_t const& end = track->End();

    bool isContained = IsInFV(start, end);
    if (isContained){
        auto const& muonRange = track_rangeP_assns_muon.at(track.key());
        if (!muonRange.empty()){
            p_muon = muonRange[0]->range_p;
        }
        auto const& protonRange = track_rangeP_assns_proton.at(track.key());
        if (!protonRange.empty()) {
            p_proton = protonRange[0]->range_p;
        }
    }
    else{
        auto const& muonMCS = track_MCS_assns_muon.at(track.key());
        if (!muonMCS.empty())
            p_muon = muonMCS[0]->fwdMomentum();
        auto const& protonMCS = track_MCS_assns_proton.at(track.key());
        if (!protonMCS.empty())
            p_proton = protonMCS[0]->fwdMomentum();
    }
    
    if (pdg==13) return p_muon;
    else if (pdg == 2212) return p_proton;
    else return -9999;
}

/**
 * @brief Estimate particle propagation times along a reconstructed track.
 *
 * Uses the track PID and initial momentum to estimate the particle velocity
 * along its trajectory. Calorimetry points are ordered using decreasing
 * residual range, and the time to each point is accumulated from the
 * interaction vertex or from the closest stored point of the parent PFP.
 *
 * For contained tracks, the momentum at each point is obtained from the
 * range-momentum calculator. For non-contained tracks, the kinetic energy is
 * decreased using a constant dE/dx value selected by particle hypothesis.
 * The resulting trajectory points and their accumulated propagation times
 * are appended to fParticleSpacePoints.
 *
 * @param pdg Particle PDG code (muon or proton).
 * @param initial_momentum Initial particle momentum in GeV/c.
 * @param track_to_calo_assns Track-to-calorimetry associations.
 * @param track Track whose trajectory is used.
 * @param PFPID Self ID of the particle-flow particle associated with the track.
 * @param parentPFPID Parent PFP ID, used to connect the trajectory to its parent.
 * @return true if at least one suitable calorimetry object is processed; false otherwise.
 */
bool sbnd::LightPropagationCorrection::GetParticlePropagationTime( int pdg, double initial_momentum, art::FindManyP<anab::Calorimetry> const& track_to_calo_assns, art::Ptr<recob::Track> const& track, int PFPID ,int parentPFPID)
{

    std::vector<ParticleSpacePoint_t> particleSpacePoints; // Temporary vector to hold the space points for this track
    
    bool success=false;
    double total_time=0;
    // Particle mass in GeV/c^2
    double mass = 0.0;
    double constant_dedx;

    if (std::abs(pdg) == 13) {
        constant_dedx = fMIPdEdX;
        mass = muon_mass;
    } else if (std::abs(pdg) == 2212) {
        constant_dedx = fHIPdEdX;
        mass = proton_mass;
    }

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
        if (calo.PlaneID().Plane != bestPlane) continue;

        const auto& xyz  = calo.XYZ();
        const auto& rr   = calo.ResidualRange();

        // Check that all vectors have the same size
        if (xyz.size() != rr.size() || rr.size()==0) {
            std::cerr << "Inconsistent calorimetry vector sizes" << std::endl;
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

        double kinetic_energy;
        double momentum;
        double total_energy;
        double beta;
        double last_beta;
        double dt;

        kinetic_energy = std::sqrt(initial_momentum * initial_momentum + mass * mass) - mass;

        total_energy = kinetic_energy + mass;

        momentum = initial_momentum;

        last_beta = momentum / total_energy;
        beta = last_beta;

        size_t first_idx = indices.front();
        if (parentPFPID == fNeutrinoID )
        {
            // Sum the distance from the interaction vertex to the first point in the trajectory and assume speed of light propagation
            double dx = xyz[first_idx].x() - fRecoVx;
            double dy = xyz[first_idx].y() - fRecoVy;
            double dz = xyz[first_idx].z() - fRecoVz;
            ds = std::sqrt(dx*dx + dy*dy + dz*dz);
            dt = ds / (last_beta * c_light);
            total_time += dt;
        }
        else{
            // Look at the existing spacepoints of the parent PFP and find the closest point to the trajectory.
            double minDistance=9999999;
            int minIdx=-1;
            for (size_t i=0; i<fParticleSpacePoints.size(); i++)
            {
                if (fParticleSpacePoints[i].pfpID==parentPFPID){
                    double dx = fParticleSpacePoints[i].x - xyz[first_idx].x();
                    double dy = fParticleSpacePoints[i].y - xyz[first_idx].y();
                    double dz = fParticleSpacePoints[i].z - xyz[first_idx].z();
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (distance < minDistance) {
                        minDistance = distance;
                        minIdx = i;
                    }
                }
            }
            if (minIdx!=-1)
            {
                double dx = xyz[first_idx].x() - fParticleSpacePoints[minIdx].x;
                double dy = xyz[first_idx].y() - fParticleSpacePoints[minIdx].y;
                double dz = xyz[first_idx].z() - fParticleSpacePoints[minIdx].z;
                ds = std::sqrt(dx*dx + dy*dy + dz*dz);
                dt = ds / (last_beta * c_light);
                total_time += dt + fParticleSpacePoints[minIdx].propagationTime;
            }
            else
            {
                double dx = xyz[first_idx].x() - fRecoVx;
                double dy = xyz[first_idx].y() - fRecoVy;
                double dz = xyz[first_idx].z() - fRecoVz;
                ds = std::sqrt(dx*dx + dy*dy + dz*dz);
                dt = ds / (last_beta * c_light);
                total_time += dt;
            }
        }

        // Add first trajectory point to the spacepoint vector
        ParticleSpacePoint_t  psp;
        psp.x = xyz[first_idx].x();
        psp.y = xyz[first_idx].y();
        psp.z = xyz[first_idx].z();
        psp.propagationTime = total_time;
        psp.pfpID = PFPID;
        particleSpacePoints.push_back(psp);

        // Tell whether the particle is contained
        geo::Point_t const& start = track->Vertex();
        geo::Point_t const& end = track->End();

        bool isContained = IsInFV(start, end);

        // Contained tracks use the range-momentum relation to obtain the momentum at each residual range.
        if (isContained)
        {
            for (size_t i = 0; i < indices.size(); ++i) {
                size_t idx = indices[i];

                // First point: no previous calorimetry point
                if (i > 0) {
                    size_t prev_idx = indices[i - 1];
                    ds = std::abs(rr[idx] - rr[prev_idx]);
                    momentum = tmc.GetTrackMomentum(rr[idx], pdg);
                    kinetic_energy = std::sqrt(momentum * momentum + mass * mass) - mass;
                    total_energy = kinetic_energy + mass;
                    beta = momentum / total_energy;
                    // Avoid negative kinetic energy

                    dt = ds / (beta * c_light);
                    total_time += dt;
                    // Update kinetic energy after traversing the segment
                    ParticleSpacePoint_t  psp;
                    psp.x = xyz[idx].x();
                    psp.y = xyz[idx].y();
                    psp.z = xyz[idx].z();
                    psp.propagationTime = total_time;
                    psp.pfpID = PFPID;
                    particleSpacePoints.push_back(psp);
                }
            }
        }
        else // For non-contained tracks, approximate the energy loss with a PID-dependent constant dE/dx.
        {
            for (size_t i = 0; i < indices.size(); ++i) {
                size_t idx = indices[i];

                // First point: no previous calorimetry point
                if (i > 0) {
                    size_t prev_idx = indices[i - 1];
                    ds = std::abs(rr[idx] - rr[prev_idx]);
                    double dE = constant_dedx * ds / 1000;
                    kinetic_energy = std::max(0.0, kinetic_energy - dE);
                    if (kinetic_energy > 0.0){
                        total_energy = kinetic_energy + mass;
                        momentum = std::sqrt(total_energy * total_energy - mass * mass);
                        beta = momentum / total_energy;
                        dt = ds / (beta * c_light);
                        total_time += dt;
                    }
                    else if(rr[idx]<fMaxRRBeta){
                        beta = 1.0;
                        dt = ds / (beta * c_light);
                        total_time += dt;
                    }
                    else return false;

                    ParticleSpacePoint_t  psp;
                    psp.x = xyz[idx].x();
                    psp.y = xyz[idx].y();
                    psp.z = xyz[idx].z();
                    psp.propagationTime = total_time;
                    psp.pfpID = PFPID;
                    particleSpacePoints.push_back(psp);
                }
            }
        }

        for(const auto& psp : particleSpacePoints) {
            fParticleSpacePoints.push_back(psp);
        }
        success=true;
    }
    return success;
}


/**
 * @brief Estimate particle propagation times using reconstructed space points.
 *
 * This fallback is used when the track-based PID and calorimetry calculation
 * is unavailable. It keeps space points associated with collection-plane
 * hits (plane 2) and estimates the propagation time from the reconstructed
 * neutrino vertex assuming the particle travels at c_light.
 *
 * @param pfpID Self ID of the particle-flow particle.
 * @param PFPSpacePointsVect Space points associated with the PFP.
 * @param SPToHitAssoc SpacePoint-to-Hit associations used to identify the plane.
 */
void sbnd::LightPropagationCorrection::GetParticlePropagationTimeLite(int pfpID, std::vector<art::Ptr<recob::SpacePoint>> const& PFPSpacePointsVect, art::FindManyP<recob::Hit> const& SPToHitAssoc)
{
    for (const art::Ptr<recob::SpacePoint> &SP: PFPSpacePointsVect){
        std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());
        if (SPHit.at(0)->WireID().Plane==2){
            double dx = SP->position().X() - fRecoVx;
            double dy = SP->position().Y() - fRecoVy;
            double dz = SP->position().Z() - fRecoVz;
            double prop_time = std::sqrt(dx*dx + dy*dy + dz*dz) / c_light;
            ParticleSpacePoint_t  psp;
            psp.x = SP->position().X();
            psp.y = SP->position().Y();
            psp.z = SP->position().Z();
            psp.propagationTime = prop_time;
            psp.pfpID = pfpID;
            fParticleSpacePoints.push_back(psp);
        }
    }

}


/**
 * @brief Calculate the Cherenkov-light propagation correction per optical channel.
 *
 * Sorts the stored trajectory points by PFP and particle propagation time.
 * For each optical detector, examines adjacent points belonging to the same
 * PFP, estimates the particle beta from their spatial and temporal separation,
 * and rejects segments below the Cherenkov threshold. The local trajectory
 * direction defines the expected Cherenkov angle through
 * cos(theta_C) = 1 / (beta * n).
 *
 * The PMT photocathode is approximated as a disk in the YZ plane. The disk
 * center and a configurable set of points on concentric rings are tested
 * against the Cherenkov-cone condition. The earliest compatible
 * particle-plus-photon arrival time is stored with the same negative sign
 * convention used by the scintillation correction.
 *
 * Distances are in cm, velocities in cm/ns, and the resulting correction is
 * in ns. The cosine tolerance and PMT disk sampling are controlled by the
 * corresponding configuration parameters.
 */
void sbnd::LightPropagationCorrection::GetPropagationTimeCorrectionPerChannelCherenkov()
{
    // Reset the Cherenkov correction for all optical channels.
    fTimeCorrectionPerChannelCherenkov.assign(fNOpChannels, 0.0);

    if (fCherenkovRefractiveIndex <= 0.0 || fVGroupVIS <= 0.0 || fCherenkovPMTRadius <= 0.0 || fCherenkovPMTRings == 0 || fCherenkovPMTAzimuthSteps == 0) {
        throw std::runtime_error(
        "GetPropagationTimeCorrectionPerChannelCherenkov: "
        "invalid Cherenkov or PMT geometry parameters."
        );
    }

    const double noTime = std::numeric_limits<double>::infinity();
    const size_t nPoints = fParticleSpacePoints.size();

    // Sort spacepoints by PFP and then by particle propagation time.
    std::vector<size_t> pointIndices(nPoints);

    for (size_t i = 0; i < nPoints; ++i) {
        pointIndices[i] = i;
    }

    std::stable_sort(pointIndices.begin(), pointIndices.end(),
    [this](size_t a, size_t b) {
        if (fParticleSpacePoints[a].pfpID != fParticleSpacePoints[b].pfpID) {
            return fParticleSpacePoints[a].pfpID < fParticleSpacePoints[b].pfpID;
        }
        return fParticleSpacePoints[a].propagationTime < fParticleSpacePoints[b].propagationTime;
    }
    );
    // Loop over optical detectors.
    for (size_t opdet = 0; opdet < fOpDetID.size(); ++opdet) {

        const double opX = fOpDetX[opdet];
        const double opY = fOpDetY[opdet];
        const double opZ = fOpDetZ[opdet];

        double minTotalTime = noTime;
        bool foundEmissionPoint = false;
        // Loop over spacepoints in sorted order.
        for (size_t j = 0; j < pointIndices.size(); ++j) {

            const size_t sp = pointIndices[j];

            const double x = fParticleSpacePoints[sp].x;
            const double y = fParticleSpacePoints[sp].y;
            const double z = fParticleSpacePoints[sp].z;
            const double particleTime = fParticleSpacePoints[sp].propagationTime;

            // Ignore spacepoints outside the PMT's TPC.
            if (x * opX <= 0.0) {
                continue;
            }

            // Find an adjacent spacepoint belonging to the same PFP.
            size_t i0 = nPoints;
            size_t i1 = nPoints;

            if (j + 1 < pointIndices.size()) {
                const size_t next = pointIndices[j + 1];
                if (fParticleSpacePoints[next].pfpID == fParticleSpacePoints[sp].pfpID) {
                    i0 = sp;
                    i1 = next;
                }
            }

            // At the end of a PFP, use the preceding segment.
            if (i0 == nPoints && j > 0) {
                const size_t previous = pointIndices[j - 1];
                if (fParticleSpacePoints[previous].pfpID == fParticleSpacePoints[sp].pfpID) {
                    i0 = previous;
                    i1 = sp;
                }
            }

            // A single point cannot define a local direction.
            if (i0 == nPoints || i1 == nPoints) {
                continue;
            }
            const double sx = fParticleSpacePoints[i1].x - fParticleSpacePoints[i0].x;
            const double sy = fParticleSpacePoints[i1].y - fParticleSpacePoints[i0].y;
            const double sz = fParticleSpacePoints[i1].z - fParticleSpacePoints[i0].z;
            const double segmentLength = std::sqrt(sx*sx + sy*sy + sz*sz);

            const double deltaTime = fParticleSpacePoints[i1].propagationTime - fParticleSpacePoints[i0].propagationTime;

            if (segmentLength <= 1.0e-6 || deltaTime <= 0.0) {
                continue;
            }

            // Estimate beta from the time between adjacent spacepoints.
            double beta = segmentLength / (c_light * deltaTime);
            if (beta > 1.05 || beta <= 0.0) continue;
            beta = std::min(beta, 1.0);
            // No Cherenkov emission below threshold.
            if (beta * fCherenkovRefractiveIndex <= 1.0) {
                continue;
            }

            const double cosCherenkov = 1.0 / (beta * fCherenkovRefractiveIndex);

            // Unit vector along the local particle trajectory.
            const double ux = sx / segmentLength;
            const double uy = sy / segmentLength;
            const double uz = sz / segmentLength;

            // Test the center and concentric rings of the PMT disk.
            // Assumption: the PMT face lies in the YZ plane.
            auto testPMTPosition = [&](double qy, double qz) {


                const double rx = opX - x;
                const double ry = qy - y;
                const double rz = qz - z;

                const double photonDistance = std::sqrt(rx*rx + ry*ry + rz*rz);
                if (photonDistance <= 1.0e-6) {
                    return;
                }

                const double cosTheta = (ux*rx + uy*ry + uz*rz) / photonDistance;
                
                // Compare the photon direction with the local track direction using the Cherenkov-cone condition.
                if (std::abs(cosTheta - cosCherenkov) >=
                fCherenkovCosineTolerance) {
                    return;
                }

                // Calculate arrival time at this position on the PMT.
                const double photonTime = photonDistance / fVGroupVIS;
                const double totalTime = particleTime + photonTime;
                
                if (std::isfinite(totalTime) && totalTime < minTotalTime) {
                    minTotalTime = totalTime;
                    foundEmissionPoint = true;
                }
            };

            // Test the center of the PMT.
            testPMTPosition(opY, opZ);

            // Test concentric rings across the disk.
            for (unsigned int ir = 1; ir <= fCherenkovPMTRings; ++ir) {
                const double radius = fCherenkovPMTRadius * static_cast<double>(ir) / fCherenkovPMTRings;
                for (unsigned int ia = 0; ia < fCherenkovPMTAzimuthSteps; ++ia) {
                    const double phi = 2.0 * M_PI * static_cast<double>(ia) / fCherenkovPMTAzimuthSteps;
                    const double qy = opY + radius * std::cos(phi);
                    const double qz = opZ + radius * std::sin(phi);
                    testPMTPosition(qy, qz);
                }
            }
        }

        // Preserve the sign convention of the existing correction.
        if (foundEmissionPoint) {
            fTimeCorrectionPerChannelCherenkov[opdet] = -minTotalTime;
        }
    }
}
