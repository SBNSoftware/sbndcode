#include "LightPropagationCorrection_module.hh"

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
    fOpHitsModuleLabel( p.get<std::string>("OpHitsModuleLabel") ),
    fSPECTDCLabel( p.get<std::string>("SPECTDCLabel") ),
    fFlashMatchingTool( p.get<std::string>("FlashMatchingTool") ),
    fSaveCorrectionTree( p.get<bool>("SaveCorrectionTree") ),
    fSpeedOfLight( p.get<double>("SpeedOfLight") ),
    fVGroupVIS( p.get<double>("VGroupVIS") ),
    fVGroupVUV( p.get<double>("VGroupVUV") ),
    fNuScoreThreshold( p.get<double>("NuScoreThreshold") ),
    fFMScoreThreshold( p.get<double>("FMScoreThreshold") ),
    fDebug( p.get<bool>("Debug", false) )
    // 
    // More initializers here.
{
    // Initialize the TimeCorrectionVector PerChannel
    for(size_t i = 0; i < fWireReadout.NOpChannels(); ++i) {
        fTimeCorrectionPerChannel.push_back(0.0); // Initialize with zero or any default value
    }

    for(unsigned int opch=0; opch<fWireReadout.NOpChannels(); opch++){
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
    //OpFlash to OpHit
    art::FindManyP<recob::OpHit> flashToOpHitAssns_tpc0(opflashListHandle_tpc0, e, fOpFlashLabel_tpc0);
    art::FindManyP<recob::OpHit> flashToOpHitAssns_tpc1(opflashListHandle_tpc1, e, fOpFlashLabel_tpc1);
    // PFP Metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_to_metadata(pfpHandle, e, fReco2Label);

    // --- Store candidate slices
    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);

    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfpVect;
    // --- Get the candidate slices
    for(size_t ix=0; ix<sliceVect.size(); ix++){
        ResetSliceInfo();
        // --- Get the slice
        auto & slice = sliceVect[ix];

        // --- Get the slice nu score and check whether we want to correct for it

        // Now I need to get all the hits associated to this flash and get the timing for all of them
        // Get the slices PFPs
        double _sliceMaxNuScore = -9999.;
        pfpVect = slice_pfp_assns.at(slice.key());
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            if(pfp->IsPrimary() &&( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ){
                const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp_to_metadata.at(pfp->Self());
                for (auto const pfpMeta : pfpMetaVec) {
                    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                    if (propertiesMap.count("NuScore")) _fNuScore = propertiesMap.at("NuScore");
                    if(_fNuScore>_sliceMaxNuScore) _sliceMaxNuScore = _fNuScore;
                }
            }
            std::vector< art::Ptr<recob::Vertex> > vertexVec = pfp_vertex_assns.at(pfp.key());
            for(const art::Ptr<recob::Vertex> &ver : vertexVec){
                geo::Point_t xyz_vertex = ver->position();
                fRecoVx= xyz_vertex.X();
                fRecoVy= xyz_vertex.Y();
                fRecoVz= xyz_vertex.Z();
            }            
            // If correct light propagation time


            // Get the SP associated to the PFP and then get the hits associated to the SP. ---> Hits associated to the PFP
            //Get the spacepoints associated to the PFParticle
            std::vector<art::Ptr<recob::SpacePoint>> PFPSpacePointsVect = pfp_sp_assns.at(pfp.key());
            //Get the SP Hit assns    
            art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
            std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
            e.getByLabel(fSpacePointLabel, eventSpacePoints);
            art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);
            art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);
            for (const art::Ptr<recob::SpacePoint> &SP: PFPSpacePointsVect){
                std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());
                if (SPHit.at(0)->WireID().Plane==2){
                    fSpacePointX.push_back(SP->position().X());
                    fSpacePointY.push_back(SP->position().Y());
                    fSpacePointZ.push_back(SP->position().Z());
                    fSpacePointIntegral.push_back(SPHit.at(0)->Integral());

                    //Fill Bayrcenter Position
                    if(SP->position().X() < 0){
                        fChargeWeightX[0] += SP->position().X() * SPHit.at(0)->Integral();
                        fChargeWeightY[0] += SP->position().Y() * SPHit.at(0)->Integral();
                        fChargeWeightZ[0] += SP->position().Z() * SPHit.at(0)->Integral();
                        fChargeTotalWeight[0] += SPHit.at(0)->Integral();
                    }
                    else{
                        fChargeWeightX[1] += SP->position().X() * SPHit.at(0)->Integral();
                        fChargeWeightY[1] += SP->position().Y() * SPHit.at(0)->Integral();
                        fChargeWeightZ[1] += SP->position().Z() * SPHit.at(0)->Integral();
                        fChargeTotalWeight[1] += SPHit.at(0)->Integral();
                    }
                }
            }
        }
        //Fill TPC 0 information
        fChargeBarycenterX[0] = fChargeWeightX[0]/fChargeTotalWeight[0];
        fChargeBarycenterY[0] = fChargeWeightY[0]/fChargeTotalWeight[0];
        fChargeBarycenterZ[0] = fChargeWeightZ[0]/fChargeTotalWeight[0];
        //Fill TPC 1 information
        fChargeBarycenterX[1] = fChargeWeightX[1]/fChargeTotalWeight[1];
        fChargeBarycenterY[1] = fChargeWeightY[1]/fChargeTotalWeight[1];
        fChargeBarycenterZ[1] = fChargeWeightZ[1]/fChargeTotalWeight[1];

        if(_sliceMaxNuScore<fNuScoreThreshold){
            ResetSliceInfo();
            continue; // Skip to the next slice if the nu score is below threshold
        }
        this->GetPropagationTimeCorrectionPerChannel();
        
        // Get the SPECTDC product required to go to the RWM reference frame

        art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
        e.getByLabel(fSPECTDCLabel, tdcHandle);
        if (!tdcHandle.isValid() || tdcHandle->size() == 0){
            std::cout << "No SPECTDC products found. Skip this event." << std::endl;
            return;
        }
        else{
            const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
            for (size_t i=0; i<tdc_v.size(); i++){
                auto tdc = tdc_v[i];
                const uint32_t  ch = tdc.Channel();
                const uint64_t  ts = tdc.Timestamp();
                if(ch == 2) fRWMTime = ts%uint64_t(1e9);
                if(ch == 4) fEventTriggerTime = ts%uint64_t(1e9);
            }
        }
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
            _fFMScore = slcTPCPMTBarycenter[BFMIdx]->score;
            if(_fFMScore < fFMScoreThreshold){
                ResetSliceInfo();
                continue;
            }
        }
        else throw art::Exception(art::errors::LogicError) << " Flash matching tool " <<  fFlashMatchingTool << " not supported ." << std::endl; 


        // Get the ophits associated to the flash
        std::vector<art::Ptr<recob::OpHit>> ophitlist;
        if(flashFM[0]->XCenter()<0)
        {
            ophitlist = flashToOpHitAssns_tpc0.at(flashFM[0].key());
            _mgr = _mgr_tpc0; // Use the TPC 0 flash finder manager
        }
        else
        {
            ophitlist = flashToOpHitAssns_tpc1.at(flashFM[0].key());
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
        double originalFlashTime = flashFM[0]->Time();
        double newFlashTime = 0.0;
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
            sbn::CorrectedOpFlashTiming correctedOpFlashTiming;
            correctedOpFlashTiming.OpFlashT0 = originalFlashTime + fEventTriggerTime/1000 - fRWMTime/1000;
            correctedOpFlashTiming.NuToFLight = (Zcenter/fSpeedOfLight)/1000;
            correctedOpFlashTiming.NuToFCharge = (fRecoVz/fSpeedOfLight)/1000;
            correctedOpFlashTiming.OpFlashT0Corrected = newFlashTime + fEventTriggerTime/1000 - fRWMTime/1000;
            correctedOpFlashTimes->emplace_back(std::move(correctedOpFlashTiming));
        }

        art::Ptr<sbn::CorrectedOpFlashTiming> newCorrectedOpFlashTimingPtr = make_correctedopflashtime_ptr(correctedOpFlashTimes->size()-1);
        newCorrectedOpFlashTimingSliceAssn->addSingle(slice, newCorrectedOpFlashTimingPtr);
        newCorrectedOpFlashTimingOpFlashAssn->addSingle(flashFM[0], newCorrectedOpFlashTimingPtr);
        if(fSaveCorrectionTree){
            this->FillCorrectionTree(newFlashTime, *flashFM[0], oldOpHitList, newOpHitList);
        }
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
        fTree->Branch("RWMTime", &fRWMTime);
        fTree->Branch("EventTriggerTime", &fEventTriggerTime);
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
    if(fSaveCorrectionTree)
    {
        fEvent = 0;
        fRun = 0;
        fSubrun = 0;
        _fNuScore = 0.0;
        fRWMTime=-99999.;
        fEventTriggerTime=-99999.;
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
    fTimeCorrectionPerChannel.resize(312, 0.0); // Reset the time correction vector for each channel
    fChargeBarycenterX.assign(2, 0.0);
    fChargeBarycenterY.assign(2, 0.0);
    fChargeBarycenterZ.assign(2, 0.0);
    fChargeWeightX.assign(2, 0.0);
    fChargeWeightY.assign(2, 0.0);
    fChargeWeightZ.assign(2, 0.0);
    fChargeTotalWeight.assign(2, 0.0);
}

void sbnd::LightPropagationCorrection::GetPropagationTimeCorrectionPerChannel()
{   
    // Implementation
    for(size_t opdet = 0; opdet < fOpDetID.size(); ++opdet) {
        double _opDetX = fOpDetX[opdet];
        double _opDetY = fOpDetY[opdet];
        double _opDetZ = fOpDetZ[opdet];
        float minPropTime = 999999999.;
        for(size_t sp=0; sp<fSpacePointX.size(); sp++)
        {
            double dx = fSpacePointX[sp] - _opDetX;
            double dy = fSpacePointY[sp] - _opDetY;
            double dz = fSpacePointZ[sp] - _opDetZ;
            double distanceToOpDet = std::sqrt(dx*dx + dy*dy + dz*dz);
            //double spToCathode = abs(fSpacePointX[sp]); // Distance from space point to cathode in mm
            //double cathodeToOpDet = std::sqrt(_opDetX*_opDetX + dy*dy + dz*dz); // Distance from cathode to OpDet in mm
            //float lightPropTimeVIS = spToCathode/fVGroupVUV + cathodeToOpDet/fVGroupVIS; // Speed

            double cathodeToOpDet = std::sqrt(_opDetX*_opDetX + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from cathode to OpDet in mm
            double spToCathode = std::sqrt( fSpacePointX[sp]*fSpacePointX[sp] + (dy/2)*(dy/2) + (dz/2)*(dz/2)); // Distance from space point to cathode in mm

            float lightPropTimeVIS = spToCathode/fVGroupVUV + cathodeToOpDet/fVGroupVIS; // Speed
            float lightPropTimeVUV = distanceToOpDet / fVGroupVUV; // Speed of light in mm/ns for VUV
            float lightPropTime = std::min(lightPropTimeVIS, lightPropTimeVUV);
            float partPropTime = std::sqrt((fSpacePointX[sp]-fRecoVx)*(fSpacePointX[sp]-fRecoVx) + (fSpacePointY[sp]-fRecoVy)*(fSpacePointY[sp]-fRecoVy) + (fSpacePointZ[sp]-fRecoVz)*(fSpacePointZ[sp]-fRecoVz))/fSpeedOfLight;
            float PropTime = lightPropTime + partPropTime;
            if(PropTime < minPropTime) minPropTime = PropTime;
        }
        fTimeCorrectionPerChannel[opdet] = -minPropTime;
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

    fSliceSPX.push_back({});
    fSliceSPY.push_back({});
    fSliceSPZ.push_back({});
    fOpHitOldTime.push_back({});
    fOpHitNewTime.push_back({});
    fOpHitPE.push_back({});
    fOpHitOpCh.push_back({});

    if(fDebug)
    {
        for(size_t i=0; i<fSpacePointX.size(); i++){
            fSliceSPX.back().push_back(fSpacePointX[i]);
            fSliceSPY.back().push_back(fSpacePointY[i]);
            fSliceSPZ.back().push_back(fSpacePointZ[i]);
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


::lightana::LiteOpHitArray_t sbnd::LightPropagationCorrection::GetAssociatedLiteHits(::lightana::LiteOpFlash_t lite_flash, ::lightana::LiteOpHitArray_t lite_hits_v)
  {
    ::lightana::LiteOpHitArray_t flash_hits_v;

    for(auto const& hitidx : lite_flash.asshit_idx) {
      flash_hits_v.emplace_back(std::move(lite_hits_v.at(hitidx)));
    }

    return flash_hits_v;
  }