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
    fOpFlashNewLabel( p.get<std::string>("OpFlashNewLabel") ),
    fSpeedOfLight( p.get<double>("SpeedOfLight") ),
    fVGroupVIS( p.get<double>("VGroupVIS") ),
    fVGroupVUV( p.get<double>("VGroupVUV") ),
    fPDFraction( p.get<double>("PDFraction") ),
    fPreWindow( p.get<double>("PreWindow") ),
    fPostWindow( p.get<double>("PostWindow") ),
    fMinHitPE( p.get<double>("MinHitPE") )
    // 
    // More initializers here.
{
    // Initialize the TimeCorrectionVector PerChannel
    for(size_t i = 0; i < 312; ++i) {
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

    produces< std::vector<recob::OpFlash>>(fOpFlashNewLabel);
    produces<art::Assns<recob::OpFlash, sbn::OpT0Finder>>();
    produces< art::Assns <recob::OpFlash, recob::Slice> >();

}

void sbnd::LightPropagationCorrection::produce(art::Event & e)
{
    std::cout << " --------------> Starting With Event " << e.id().event() << std::endl;
    std::unique_ptr< std::vector<recob::OpFlash> > opflashes(new std::vector<recob::OpFlash>);
    
    std::unique_ptr< art::Assns<recob::OpFlash, recob::Slice>> newOpFlashSliceAssn (new art::Assns<recob::OpFlash, recob::Slice>);
    art::PtrMaker<recob::OpFlash> make_opflash_ptr{e, fOpFlashNewLabel};
    std::unique_ptr< art::Assns<recob::OpFlash, sbn::OpT0Finder>> newOpFlashOpT0Assn (new art::Assns<recob::OpFlash, sbn::OpT0Finder>);
    // --- Read Recob Slice
    ::art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fReco2Label, sliceHandle);
    // Slice to OpT0Finder
    //Get the handle for making the assns later on
    ::art::Handle<std::vector<sbn::OpT0Finder>> opt0Handle;
    e.getByLabel(fOpT0FinderModuleLabel, opt0Handle);
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
        // --- Get the OpT0Finder score
        // Get all the OpT0 objects associated to the slice
        const std::vector< art::Ptr<sbn::OpT0Finder> > slcOpT0Finder = slice_opt0finder_assns.at( slice.key() );
        if(slcOpT0Finder.size() == 0){
            //std::cout << " No OpT0Finder objects associated to this slice." << std::endl;
            continue; // Skip to the next slice if no OpT0Finder objects are found
        }
        size_t OpT0Idx = HighestOpT0ScoreIdx(slcOpT0Finder);
        // Get the flash OpT0 association
        art::FindManyP<recob::OpFlash> opflash_opt0finder_assns(opt0Handle, e, fOpT0FinderModuleLabel);
        auto & flashOpT0 = opflash_opt0finder_assns.at( slcOpT0Finder[OpT0Idx].key() );
        double oldFlashTime = flashOpT0[0]->Time();
        if(flashOpT0.size() > 1){
            throw art::Exception(art::errors::LogicError) << "There are multiple OpFlash objects associated to the same OpT0Finder object. This is not expected.";
        }

        // Now I need to get all the hits associated to this flash and get the timing for all of them
        // Get the slices PFPs
        pfpVect = slice_pfp_assns.at(slice.key());
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
            //std::cout<<"   ** PFParticle: "<<pfp->Self()<<"      PDG:"<<pfp->PdgCode()<<"  Primary="<<pfp->IsPrimary()<<" Mother="<<pfp->Parent()<<std::endl;
            std::vector< art::Ptr<recob::Vertex> > vertexVec = pfp_vertex_assns.at(pfp.key());
            for(const art::Ptr<recob::Vertex> &ver : vertexVec){
                geo::Point_t xyz_vertex = ver->position();
                fRecoVx= xyz_vertex.X();
                fRecoVy= xyz_vertex.Y();
                fRecoVz= xyz_vertex.Z();
                //std::cout << " Reconstructed vertex coordinates: (" << fRecoVx << ", " << fRecoVy << ", " << fRecoVz << ")" << std::endl;
            }
            
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

                    //std::cout << " Adding SpacePoint at (" 
                    //        << SP->position().X() << ", " 
                    //        << SP->position().Y() << ", " 
                    //        << SP->position().Z() << ") with integral " 
                    //        << SPHit.at(0)->Integral() << " from PFParticle " << pfp->Self() << std::endl;
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

        this->GetPropagationTimeCorrectionPerChannel();

        // Get the ophits associated to the flash
        std::vector<art::Ptr<recob::OpHit>> ophitlist;
        if(flashOpT0[0]->XCenter()<0) ophitlist = flashToOpHitAssns_tpc0.at(flashOpT0[0].key());
        else ophitlist = flashToOpHitAssns_tpc1.at(flashOpT0[0].key());
        std::vector<recob::OpHit> newOpHitList;
        int _nophits = 0;
        _nophits += ophitlist.size();
        for (int i = 0; i < _nophits; ++i) {
            int opCh = ophitlist.at(i)->OpChannel();
            double channelCorrection = fTimeCorrectionPerChannel[opCh];
            double newPeakTime = ophitlist.at(i)->StartTime()+ophitlist.at(i)->RiseTime() + channelCorrection/1000;
            double peakTimeAbs = ophitlist.at(i)->PeakTimeAbs()+ channelCorrection/1000;
            double newStartTime = ophitlist.at(i)->StartTime() + channelCorrection/1000;
            double riseTime = ophitlist.at(i)->RiseTime();
            unsigned int frame = ophitlist.at(i)->Frame();
            double width = ophitlist.at(i)->Width();
            double area = ophitlist.at(i)->Area();
            double amplitude = ophitlist.at(i)->Amplitude();
            double pe = ophitlist.at(i)->PE();
            recob::OpHit newOpHit = recob::OpHit(opCh, newPeakTime, peakTimeAbs, newStartTime, riseTime, frame, width, area, amplitude, pe, 0.0);            
            newOpHitList.push_back(newOpHit);
            //std::cout << " OpHit old time " << ophitlist.at(i)->StartTime() << " new ophit time " << newStartTime << " with PEs " << pe << std::endl;

            //newOpHitList.push_back(newOpHit);
            //std::cout << " Correcting ophit time by " << channelCorrection 
            //            << " for channel " << opCh 
            //            << " with new PeakTimeAbs: " << ophitlist.at(i)->PeakTimeAbs() 
            //            << " and StartTime: " << ophitlist.at(i)->StartTime() << std::endl;
        }

        double newFlashTime = GetFlashT0(oldFlashTime, newOpHitList);
        //std::cout << " The new flash time is " << newFlashTime << std::endl;

        double 	timewidth = flashOpT0[0]->TimeWidth();
        double abstime = flashOpT0[0]->AbsTime();
        unsigned int frame = flashOpT0[0]->Frame();
        std::vector< double > PEperOpDet = flashOpT0[0]->PEs();
        double xCenter = flashOpT0[0]->XCenter();
        double xWidth = flashOpT0[0]->XWidth();
        double yCenter = flashOpT0[0]->YCenter();
        double yWidth = flashOpT0[0]->YWidth();
        double zCenter = flashOpT0[0]->ZCenter();
        double zWidth = flashOpT0[0]->ZWidth();
        recob::OpFlash newFlash(
            newFlashTime, timewidth, abstime, frame, PEperOpDet,
            0, 0, 1,
            xCenter, xWidth, yCenter, yWidth, zCenter, zWidth);

        opflashes->emplace_back(std::move(newFlash)); 
        std::cout << " Assns between slice " << slice->ID() << " with charge barycenter on tpc 0 " << fChargeBarycenterX[0] << " Y " << fChargeBarycenterY[0] << " Z " << fChargeBarycenterZ[0]
                  << " and tpc 1, X " << fChargeBarycenterX[1] << " Y " << fChargeBarycenterY[1] << " Z " << fChargeBarycenterZ[1]
                  << " OpT0Score " << slcOpT0Finder[OpT0Idx]->score
                  << " and OpFlash at time " << newFlashTime
                  << " with XCenter " << xCenter
                  << " and YCenter " << yCenter
                  << " and ZCenter " << zCenter
                  << std::endl;

        art::Ptr<recob::OpFlash> newOpFlashPtr = make_opflash_ptr(opflashes->size()-1);
        //if(newOpFlashPtr.isNull()) {
        //    throw art::Exception(art::errors::LogicError) << "Failed to create OpFlash pointer.";
        //}
        newOpFlashOpT0Assn->addSingle(newOpFlashPtr, slcOpT0Finder[OpT0Idx]);
        newOpFlashSliceAssn->addSingle(newOpFlashPtr, slice);
    }
    e.put(std::move(opflashes),fOpFlashNewLabel);
    e.put(std::move(newOpFlashOpT0Assn));
    e.put(std::move(newOpFlashSliceAssn));

    return;
}

void sbnd::LightPropagationCorrection::beginJob()
{

}

void sbnd::LightPropagationCorrection::endJob()
{

}

void sbnd::LightPropagationCorrection::ResetEventVars()
{

}

size_t sbnd::LightPropagationCorrection::HighestOpT0ScoreIdx(const std::vector< art::Ptr<sbn::OpT0Finder> > slcOpT0Finder)
{
    // Gets the idx of the OpT0 object with the highest score 
    double highestOpT0Score = -99999.0; // Initialize to a negative value
    size_t highestIdx = 0;
    for(size_t jx=0; jx<slcOpT0Finder.size(); jx++){
        if(slcOpT0Finder[jx]->score > highestOpT0Score){
            highestOpT0Score = slcOpT0Finder[jx]->score;
            highestIdx = jx;
        }
    }
    return highestIdx;
}

void sbnd::LightPropagationCorrection::ResetSliceInfo()
{
    fRecoVx = 0.0;
    fRecoVy = 0.0;
    fRecoVz = 0.0;
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
        size_t closestSPIdx=-999999;
        double minDistance = 999999999.;
        for(size_t sp=0; sp<fSpacePointX.size(); sp++)
        {
            double dx = fSpacePointX[sp] - _opDetX;
            double dy = fSpacePointY[sp] - _opDetY;
            double dz = fSpacePointZ[sp] - _opDetZ;
            double distanceToOpDet = std::sqrt(dx*dx + dy*dy + dz*dz);
            if(distanceToOpDet < minDistance)
            {
                minDistance = distanceToOpDet;
                closestSPIdx = sp;
            }  
        }
        double LightCorrection = GetPropagationTime(fSpacePointX[closestSPIdx]); // Speed of light in mm/ns
        double ParticlePropagationCorrection = std::hypot((fSpacePointX[closestSPIdx]-fRecoVx), (fSpacePointY[closestSPIdx]-fRecoVy), (fSpacePointZ[closestSPIdx]-fRecoVz))/fSpeedOfLight;
        double correction = ParticlePropagationCorrection + LightCorrection;
        fTimeCorrectionPerChannel[opdet] = -correction;
    }
}


double sbnd::LightPropagationCorrection::GetPropagationTime(double drift){
    // drift is here the X coordinate
    // cathode: x=0 cm, PDS: x=200 cm
    //std::cout << "Drift distance: " << fDriftDistance << " X coord " << drift <<std::endl;
    if(std::abs(drift) > fKinkDistance)
      return (fDriftDistance-std::abs(drift)) * fVGroupVUV_I ;
    else
      return std::abs(drift) * fVGroupVUV_I + fVISLightPropTime;
}

double sbnd::LightPropagationCorrection::GetFlashT0(double flash_time, std::vector<recob::OpHit> ophitlist){

    std::vector< std::pair<double, double> > selected_hits;
    double pe_sum = 0;

    // fill vector with selected hits in the specified window
    for(auto const& hit : ophitlist) {
      if( hit.PeakTime()<flash_time+fPostWindow && hit.PeakTime()>flash_time-fPreWindow && hit.PE()>fMinHitPE){
        selected_hits.push_back( std::make_pair(hit.PE(), hit.PeakTime()));
        pe_sum += hit.PE();
      }
    }

    if(pe_sum>0){
      // sort vector by number of #PE (ascending order)
      std::sort( selected_hits.begin(), selected_hits.end(), std::greater< std::pair<double, double> >() );

      double flasht0_mean=0., pe_count=0.;
      int nophits=0;

      // loop over selected ophits
      for (size_t ix=0; ix<selected_hits.size(); ix++) {
        //std::cout << " Summing ophit with PE " << selected_hits[ix].first << " and time " << selected_hits[ix].second << std::endl;
        pe_count += selected_hits[ix].first;
        flasht0_mean += selected_hits[ix].second;
        nophits++;
        if( pe_count/pe_sum>fPDFraction ) break;
      }
      return flasht0_mean/nophits;
    }
    else
      return flash_time;
  }