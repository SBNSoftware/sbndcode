#include "sbndcode/OpDetAnalyzer/PDSAnalyzer/SBNDPDSAnalyzer_module.hh"


// -------- Constructor --------
opdet::SBNDPDSAnalyzer::SBNDPDSAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fSaveMCTruth( p.get<bool>("SaveMCTruth") ),
  fSaveMCParticles( p.get<bool>("SaveMCParticles") ),
  fSaveSimPhotons( p.get<bool>("SaveSimPhotons") ),
  fSaveSimPhotonsArrivalTimes( p.get<bool>("SaveSimPhotonsArrivalTimes") ),
  fSaveRawWaveforms( p.get<bool>("SaveRawWaveforms") ),
  fSaveDeconvolvedWaveforms( p.get<bool>("SaveDeconvolvedWaveforms") ),
  fSaveOpHits( p.get<bool>("SaveOpHits") ),
  fSaveOpFlashes( p.get<bool>("SaveOpFlashes") ),
  fSaveCosmicId( p.get<bool>("SaveCosmicId") ),
  fVerbosity( p.get<int>("Verbosity") ),
  fMakePerTrackTree( p.get<bool>("MakePerTrackTree") ),
  fMakePDSGeoTree( p.get<bool>("MakePDSGeoTree") ),
  fUseSimPhotonsLite( p.get<bool>("UseSimPhotonsLite") ),
  fPDTypes( p.get<std::vector<std::string>>("PDTypes") ),
  fKeepPDGCode( p.get<std::vector<int>>("KeepPDGCode") ),
  fMCTruthOrigin( p.get<std::vector<int>>("MCTruthOrigin") ),
  fMCTruthPDG( p.get<std::vector<int>>("MCTruthPDG") ),
  fMCTruthModuleLabel( p.get<std::vector<std::string>>("MCTruthModuleLabel") ),
  fMCTruthInstanceLabel( p.get<std::vector<std::string>>("MCTruthInstanceLabel") ),
  fMCModuleLabel( p.get< std::string >("MCModuleLabel") ),
  fSimPhotonsModuleLabel( p.get<std::vector<std::string>>("SimPhotonsModuleLabel") ),
  fRawWaveformsModuleLabel( p.get< std::string >("RawWaveformsModuleLabel") ),
  fDeconvolvedWaveformsModuleLabel( p.get< std::string >("DeconvolvedWaveformsModuleLabel") ),
  fOpHitsModuleLabel( p.get<std::vector<std::string>>("OpHitsModuleLabel") ),
  fOpFlashesModuleLabel( p.get<std::vector<std::string>>("OpFlashesModuleLabel") ),
  fHitsLabel( p.get<std::string>("HitsLabel") ),
  fReco2Label( p.get<std::string>("Reco2Label") ),
  fCosmicIdModuleLabel( p.get< std::string >("CosmicIdModuleLabel") ),
  fOpT0FinderModuleLabel( p.get< std::string >("OpT0FinderModuleLabel") ),
  fSimpleFlashMatchModuleLabel( p.get< std::string >("SimpleFlashMatchModuleLabel") ),
  fG4BufferBoxX( p.get<std::vector<int>>("G4BufferBoxX") ),
  fG4BufferBoxY( p.get<std::vector<int>>("G4BufferBoxY") ),
  fG4BufferBoxZ( p.get<std::vector<int>>("G4BufferBoxZ") ),
  fG4BeamWindow( p.get<std::vector<int>>("G4BeamWindow") )
{

  // MakePerTrackTree mode requires UseSimPhotonsLite false and SaveMCParticles true
  if(fMakePerTrackTree && fUseSimPhotonsLite){
    std::cout << "MakePerTrackTree mode requires UseSimPhotonsLite false." << std::endl;
    throw std::exception();
  }
  if(fMakePerTrackTree && !fSaveMCParticles){
    std::cout << "MakePerTrackTree mode requires SaveMCParticles true." << std::endl;
    throw std::exception();
  }

}


// -------- Begi job function --------
void opdet::SBNDPDSAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("OpAnaTree", "Optical Analyzer Tree");

  fTree->Branch("eventID", &_eventID, "eventID/i");
  fTree->Branch("runID", &_runID, "runID/i");
  fTree->Branch("subrunID", &_subrunID, "subrunID/i");
  
  // MCTruth
  if(fSaveMCTruth){
    fTree->Branch("nuvX","std::vector<double>", &_nuvX);
    fTree->Branch("nuvY","std::vector<double>", &_nuvY);
    fTree->Branch("nuvZ","std::vector<double>", &_nuvZ);
    fTree->Branch("nuvT","std::vector<double>", &_nuvT);
    fTree->Branch("nuvE","std::vector<double>", &_nuvE);
  }

  // MCParticles
  if(fSaveMCParticles){
    fTree->Branch("stepX",&_mc_stepX);
    fTree->Branch("stepY",&_mc_stepY);
    fTree->Branch("stepZ",&_mc_stepZ);
    fTree->Branch("stepT",&_mc_stepT);
    fTree->Branch("dE",&_mc_dE);
    fTree->Branch("energydep",&_mc_energydep);
    fTree->Branch("energydepX",&_mc_energydepX);
    fTree->Branch("energydepY",&_mc_energydepY);
    fTree->Branch("energydepZ",&_mc_energydepZ);
    fTree->Branch("E",&_mc_E);
    fTree->Branch("StartPx",&_mc_StartPx);
    fTree->Branch("StartPy",&_mc_StartPy);
    fTree->Branch("StartPz",&_mc_StartPz);
    fTree->Branch("EndPx",&_mc_EndPx);
    fTree->Branch("EndPy",&_mc_EndPy);
    fTree->Branch("EndPz",&_mc_EndPz);
    fTree->Branch("process",&_mc_process);
    fTree->Branch("trackID",&_mc_trackID);
    fTree->Branch("motherID",&_mc_motherID);
    fTree->Branch("PDGcode",&_mc_PDGcode);
    fTree->Branch("InTimeCosmics", &_mc_InTimeCosmics);
    fTree->Branch("InTimeCosmicsTime", &_mc_InTimeCosmicsTime);
    fTree->Branch("dEtpc",&_mc_dEtpc);
    fTree->Branch("dEpromx",&_mc_dEpromx);
    fTree->Branch("dEpromy",&_mc_dEpromy);
    fTree->Branch("dEpromz",&_mc_dEpromz);
    fTree->Branch("dEspreadx",&_mc_dEspreadx);
    fTree->Branch("dEspready",&_mc_dEspready);
    fTree->Branch("dEspreadz",&_mc_dEspreadz);
    fTree->Branch("dElowedges",&_mc_dElowedges);
    fTree->Branch("dEmaxedges",&_mc_dEmaxedges);
  }

  // SimPhotons
  if(fSaveSimPhotons){
    fTree->Branch("SimPhotonsperOpChVUV",&_simPhotonsperOpChVUV);
    fTree->Branch("SimPhotonsperOpChVIS",&_simPhotonsperOpChVIS);
    fTree->Branch("NPhotons",&_NPhotons);
    fTree->Branch("NPhotonsPMTCo",&_NPhotonsPMTCo);
    fTree->Branch("NPhotonsPMTUnco",&_NPhotonsPMTUnco);
    fTree->Branch("NPhotonsPMTCoVUV",&_NPhotonsPMTCoVUV);
    fTree->Branch("NPhotonsXARAPUCAVUV",&_NPhotonsXARAPUCAVUV);
    fTree->Branch("NPhotonsXARAPUCAVIS",&_NPhotonsXARAPUCAVIS);
    if(fSaveSimPhotonsArrivalTimes){
      fTree->Branch("SimPhotonsLiteVUV",&_simPhotonsLiteVUV);
      fTree->Branch("SimPhotonsLiteVIS",&_simPhotonsLiteVIS);
    }
  }

  // Raw waveforms
  if(fSaveRawWaveforms){
    fTree->Branch("SignalsDigi", "std::vector<std::vector<double>>",&_signalsDigi);
    fTree->Branch("StampTime", "std::vector<double>",&_stampTime);
    fTree->Branch("OpChDigi", "std::vector<int>",&_opChDigi);
  }

  // Deconvolved waveforms
  if(fSaveDeconvolvedWaveforms){
    fTree->Branch("SignalsDeco", "std::vector<std::vector<double>>",&_signalsDeco);
    fTree->Branch("StampTimeDeco", "std::vector<double>",&_stampTimeDeco);
    fTree->Branch("OpChDeco", "std::vector<int>",&_opChDeco);
  }

  // OpHits
  if(fSaveOpHits){
    fTree->Branch("nophits", &_nophits, "nophits/I");
    fTree->Branch("ophit_opch", "std::vector<int>",&_ophit_opch);
    fTree->Branch("ophit_peakT", "std::vector<double>",&_ophit_peakT);
    fTree->Branch("ophit_startT", "std::vector<double>",&_ophit_startT);
    fTree->Branch("ophit_riseT", "std::vector<double>",&_ophit_riseT);
    fTree->Branch("ophit_width", "std::vector<double>",&_ophit_width);
    fTree->Branch("ophit_area", "std::vector<double>",&_ophit_area);
    fTree->Branch("ophit_amplitude", "std::vector<double>",&_ophit_amplitude);
    fTree->Branch("ophit_pe", "std::vector<double>",&_ophit_pe);
  }

  // OpFlashes
  if(fSaveOpFlashes){
    fTree->Branch("nopflash", &_nopflash, "nopflash/I");
    fTree->Branch("flash_id","std::vector<int>", &_flash_id);
    fTree->Branch("flash_time","std::vector<double>", &_flash_time);
    fTree->Branch("flash_total_pe", "std::vector<double>", &_flash_total_pe);
    fTree->Branch("flash_pe_v","std::vector<std::vector<double>>", &_flash_pe_v);
    fTree->Branch("flash_tpc", "std::vector<int>", &_flash_tpc);
    fTree->Branch("flash_y","std::vector<double>", &_flash_y);
    fTree->Branch("flash_yerr", "std::vector<double>", &_flash_yerr);
    fTree->Branch("flash_z","std::vector<double>", &_flash_z);
    fTree->Branch("flash_zerr", "std::vector<double>", &_flash_zerr);
    fTree->Branch("flash_x","std::vector<double>", &_flash_x);
    fTree->Branch("flash_xerr", "std::vector<double>", &_flash_xerr);
    fTree->Branch("flash_ophit_time", "std::vector<std::vector<double>>", &_flash_ophit_time);
    fTree->Branch("flash_ophit_risetime", "std::vector<std::vector<double>>", &_flash_ophit_risetime);
    fTree->Branch("flash_ophit_starttime", "std::vector<std::vector<double>>", &_flash_ophit_starttime);
    fTree->Branch("flash_ophit_amp", "std::vector<std::vector<double>>", &_flash_ophit_amp);
    fTree->Branch("flash_ophit_area", "std::vector<std::vector<double>>", &_flash_ophit_area);
    fTree->Branch("flash_ophit_width", "std::vector<std::vector<double>>", &_flash_ophit_width);
    fTree->Branch("flash_ophit_pe", "std::vector<std::vector<double>>", &_flash_ophit_pe);
    fTree->Branch("flash_ophit_ch", "std::vector<std::vector<int>>", &_flash_ophit_ch);
  }

  // Cosmic ID
  if(fSaveCosmicId){
    fTree->Branch("CRUMBSScore","std::vector<double>", &_CRUMBSScore);
    fTree->Branch("OpT0Score","std::vector<double>", &_opT0Score);
    fTree->Branch("OpT0Chi2","std::vector<double>", &_opT0Chi2);
    fTree->Branch("OpT0Time","std::vector<double>", &_opT0Time);
    fTree->Branch("OpT0HypoPE","std::vector<double>", &_opT0HypoPE);
    fTree->Branch("OpT0MeasPE","std::vector<double>", &_opT0MeasPE);
    fTree->Branch("SliceOrigin","std::vector<int>", &_sliceOrigin);
    fTree->Branch("SliceCompleteness","std::vector<double>", &_sliceCompleteness);
    fTree->Branch("SlicePurity","std::vector<double>", &_slicePurity);
    fTree->Branch("SFMSliceOrigin","std::vector<int>", &_sFMSliceOrigin);
    fTree->Branch("SFMScore","std::vector<double>", &_sFMScore);
    fTree->Branch("SFMScoreY","std::vector<double>", &_sFMScoreY);
    fTree->Branch("SFMScoreZ","std::vector<double>", &_sFMScoreZ);
    fTree->Branch("SFMScoreRR","std::vector<double>", &_sFMScoreRR);
    fTree->Branch("SFMScoreRatio","std::vector<double>", &_sFMScoreRatio);
    fTree->Branch("SFMScoreSlope","std::vector<double>", &_sFMScoreSlope);
    fTree->Branch("SFMScorePEtoQ","std::vector<double>", &_sFMScorePEtoQ);
    fTree->Branch("SFMTime","std::vector<double>", &_sFMTime);
    fTree->Branch("SFMPE","std::vector<double>", &_sFMPE);
  }

  // Make per track ID
  if(fMakePerTrackTree){
    fPerTrackTree = tfs->make<TTree>("OpAnaPerTrackTree", "Optical Analyzer Tree (all tracks)");
    fPerTrackTree->Branch("eventID", &_eventID, "eventID/i");
    fPerTrackTree->Branch("runID", &_runID, "runID/i");
    fPerTrackTree->Branch("subrunID", &_subrunID, "subrunID/i");
    fPerTrackTree->Branch("SimPhotonsVUV",&_simPhotonsVUV);
    fPerTrackTree->Branch("SimPhotonsVIS",&_simPhotonsVIS);
    fPerTrackTree->Branch("TrackID",&_perTrackID);
  }

  if(fMakePDSGeoTree){
    fPDSMapTree = tfs->make<TTree>("PDSMapTree", "PDS Map Tree");
    fPDSMapTree->Branch("OpDetID", "std::vector<int>", &_opDetID);
    fPDSMapTree->Branch("OpDetX", "std::vector<double>", &_opDetX);
    fPDSMapTree->Branch("OpDetY", "std::vector<double>", &_opDetY);
    fPDSMapTree->Branch("OpDetZ", "std::vector<double>", &_opDetZ);
    fPDSMapTree->Branch("OpDetType", "std::vector<int>", &_opDetType);

    if(fVerbosity>0)
      std::cout << " -- Dumping SBND PDS mapping -- \n";

    for(unsigned int opch=0; opch<fGeoService->NOpChannels(); opch++){
      auto pdCenter = fGeoService->OpDetGeoFromOpChannel(opch).GetCenter();
      _opDetID.push_back(opch);
      _opDetX.push_back(pdCenter.X());
      _opDetY.push_back(pdCenter.Y());
      _opDetZ.push_back(pdCenter.Z());
      if(fPDSMap.pdType(opch)=="pmt_coated") _opDetType.push_back(kPMTCoated);
      else if(fPDSMap.pdType(opch)=="pmt_uncoated") _opDetType.push_back(kPMTUncoated);
      else if(fPDSMap.pdType(opch)=="xarapuca_vuv") _opDetType.push_back(kXARAPUCAVUV);
      else if(fPDSMap.pdType(opch)=="xarapuca_vis") _opDetType.push_back(kXARAPUCAVIS);
      else _opDetType.push_back(kPDUnknown);

      if(fVerbosity>0)
        std::cout << opch << " " << pdCenter.X() << " " << pdCenter.Y() << " " << pdCenter.Z() << " " << fPDSMap.pdType(opch) << std::endl;
    
    }

    fPDSMapTree->Fill();

  }

}


// -------- Main function --------
void opdet::SBNDPDSAnalyzer::analyze(art::Event const& e)
{

  // --- Event General Info
  _eventID = e.id().event();
  _runID = e.id().run();
  _subrunID = e.id().subRun();
  if(fVerbosity>0)
    std::cout << " -- Running SBNDPDSAnalyzer -- \n Run=" << _runID << " Subrun=" << _subrunID << " Event=" << _eventID << std::endl;
  
  // --- Services
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<detinfo::DetectorClocksService> timeservice;

  // --- Saving MCTruths
  if(fSaveMCTruth){
    _nuvX.clear(); _nuvY.clear(); _nuvZ.clear(); _nuvT.clear(); _nuvE.clear();
    FillMCTruth(e);
  }

  // --- Saving MCParticles
  if(fSaveMCParticles){
    
    art::Handle< std::vector<simb::MCParticle> > mclistLARG4;
    e.getByLabel(fMCModuleLabel,mclistLARG4);
    if(!mclistLARG4.isValid()){
      std::cout << " MC particles with label " << fMCModuleLabel << " not found. " << std::endl;
      throw std::exception();
    }
    std::vector<simb::MCParticle> const& mcpartVec(*mclistLARG4);

    // Initialize vectors and map for SimPhotons (in MakePerTrackTree mode)
    fSimPhotonsVUVMap.clear(); 
    fSimPhotonsVISMap.clear();
    if(fMakePerTrackTree){

      std::vector<std::vector<double>> fPhVTemplate;
      fPhVTemplate.resize( fGeoService->NOpChannels() );
      
      for(size_t i_p=0; i_p < mcpartVec.size(); i_p++){
        const simb::MCParticle pPart = mcpartVec[i_p];
        fSimPhotonsVUVMap.insert(std::pair< int, std::vector<std::vector<double>> >(pPart.TrackId(), fPhVTemplate));
        fSimPhotonsVISMap.insert(std::pair< int, std::vector<std::vector<double>> >(pPart.TrackId(), fPhVTemplate));
      }
      
    }

    // Reset TTree variables
    double dE_neutrinowindow=0;
    _mc_stepX.clear(); _mc_stepY.clear(); _mc_stepZ.clear(); _mc_stepT.clear();
    _mc_dE.clear(); _mc_E.clear();
    _mc_trackID.clear(); _mc_motherID.clear(); _mc_PDGcode.clear(); _mc_process.clear();
    _mc_StartPx.clear(); _mc_StartPy.clear(); _mc_StartPz.clear();
    _mc_EndPx.clear(); _mc_EndPy.clear(); _mc_EndPz.clear();
    _mc_energydep.clear(); _mc_energydepX.clear(); _mc_energydepY.clear(); _mc_energydepZ.clear();
    _mc_InTimeCosmics=0; _mc_InTimeCosmicsTime.clear();

    if(fVerbosity>0)
      std::cout << "Saving MCParticles from" << fMCModuleLabel << std::endl;

    // Loop over the handle
    for(size_t i_p=0; i_p < mcpartVec.size(); i_p++){
      
      const simb::MCParticle pPart = mcpartVec[i_p];

      // If a list of PDG codes is provided in the fhicl file, skip particles not in the list
      if(fKeepPDGCode.size()!=0 && std::find(fKeepPDGCode.begin(),fKeepPDGCode.end(),pPart.PdgCode()) == fKeepPDGCode.end())
        continue;

      // Get the MCParticle trajectory
      const simb::MCTrajectory truetrack = pPart.Trajectory();
      
      // Vectors to store the MCParticle steps
      std::vector<double> xpoints, ypoints, zpoints, tpoints;

      // Save at least the first position points (vertex) for all particles
      xpoints.push_back(pPart.Position(0).X());
      ypoints.push_back(pPart.Position(0).Y());
      zpoints.push_back(pPart.Position(0).Z());
      tpoints.push_back(pPart.Position(0).T());

      // Loop over the trajectory points
      for(size_t i_s=1; i_s < pPart.NumberTrajectoryPoints(); i_s++){

        double t = pPart.Position(i_s).T();
        double x = pPart.Position(i_s).X();
        double y = pPart.Position(i_s).Y();
        double z = pPart.Position(i_s).Z();

        // Only store trajectory points in the G4BufferBox and G4BeamTimeWindow specified in the fhicl file
        if( x < fG4BufferBoxX.at(0) || x > fG4BufferBoxX.at(1)
              || y < fG4BufferBoxY.at(0) || y > fG4BufferBoxY.at(1)
                  || z < fG4BufferBoxZ.at(0) || z > fG4BufferBoxZ.at(1)) continue;

        if( t < fG4BeamWindow.at(0) || t > fG4BeamWindow.at(1) ) continue;
        
        xpoints.push_back(x);
        ypoints.push_back(y);
        zpoints.push_back(z);
        tpoints.push_back(t);
      }

      // Save the MCParticle information
      _mc_stepX.push_back(xpoints);
      _mc_stepY.push_back(ypoints);
      _mc_stepZ.push_back(zpoints);
      _mc_stepT.push_back(tpoints);
      _mc_E.push_back(pPart.E());
      _mc_process.push_back(pPart.Process());
      _mc_trackID.push_back(pPart.TrackId());
      _mc_motherID.push_back(pPart.Mother());
      _mc_PDGcode.push_back(pPart.PdgCode());
      _mc_StartPx.push_back(pPart.Px(0));
      _mc_StartPy.push_back(pPart.Py(0));
      _mc_StartPz.push_back(pPart.Pz(0));
      _mc_EndPx.push_back(pPart.EndPx());
      _mc_EndPy.push_back(pPart.EndPy());
      _mc_EndPz.push_back(pPart.EndPz());
      
      // Save energy deposition by the MCParticle
      // Only save energy deposition step for particles G4BeamTimeWindow specified in the fhicl file
      // Also checks in time cosmics depositing eneegy in the TPC
      // Initialize variables to store energy deposition by the MCParticle
      double endep=0;
      std::vector<double> truedE; truedE.clear();
      std::vector<double>  truedE_vecX; truedE_vecX.clear();
      std::vector<double> truedE_vecY; truedE_vecY.clear();
      std::vector<double> truedE_vecZ; truedE_vecZ.clear();

      // Get the associated SimIDEs
      std::vector< const sim::IDE * > ides_v = bt_serv->TrackIdToSimIDEs_Ps (pPart.TrackId());
      for(auto *ide:ides_v){
        // need to divide by 3 to avoid triple counting (3 wire planes)
        endep+=ide->energy/3.;
        if( pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1) ){
          truedE.push_back(ide->energy/3.);
          truedE_vecX.push_back(ide->x);
          truedE_vecY.push_back(ide->y);
          truedE_vecZ.push_back(ide->z);
          dE_neutrinowindow+=ide->energy/3.;
        }
      }

      // Save energy deposition by the MCParticle
      _mc_dE.push_back(endep);
      _mc_energydep.push_back(truedE);
      _mc_energydepX.push_back(truedE_vecX);
      _mc_energydepY.push_back(truedE_vecY);
      _mc_energydepZ.push_back(truedE_vecZ);

      // Check if the MCParticle is a cosmic in time coincidence with the G4BeamTimeWindow specified in the fhicl file
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(pPart.TrackId());
      if( pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1) && truth->Origin()==2 && endep>1.){
        _mc_InTimeCosmics++;
        _mc_InTimeCosmicsTime.push_back(pPart.Position(0).T());
      }

      // Print-out MCParticle information
      if(fVerbosity>0){
        std::cout<<std::setw(4)<<i_p<<"-.-.-.Energy="<<std::setw(12)<<pPart.E()<<" PDGCODE="<<std::setw(10)<<pPart.PdgCode();
        std::cout<<" ID="<<std::setw(6)<<pPart.TrackId()<<" Mother="<<std::setw(6)<<pPart.Mother()<<" dE="<<std::setw(9)<<endep;
        std::cout<<" T="<<std::setw(10)<<pPart.T()<<" X="<<std::setw(7)<<pPart.Vx()<<" Y="<<std::setw(7)<<pPart.Vy()<<" Z="<<std::setw(7)<<pPart.Vz()<<" NPoints:"<<pPart.NumberTrajectoryPoints();
        std::cout<<" Process: "<<pPart.Process()<<" EndProcess: "<<pPart.EndProcess()<<std::endl;
      }

    }

    // Fill average energy deposition variables
    FillAverageDepositedEnergyVariables(_mc_energydep,_mc_energydepX,_mc_energydepY,_mc_energydepZ,_mc_stepT,_mc_dEtpc,_mc_dEpromx,_mc_dEpromy,_mc_dEpromz,_mc_dEspreadx,_mc_dEspready,_mc_dEspreadz,_mc_dElowedges,_mc_dEmaxedges);

    // Print-out MCParticles summary
    if(fVerbosity>0){
      std::cout<<" InTimeCosmic "<<_mc_InTimeCosmics<<std::endl;
      std::cout<<" Energy Deposition during Beam Window="<<dE_neutrinowindow<<std::endl;
      std::cout<<"----ENERGY DEPOSITIONS\n";
      std::cout<<"*** TPC0   dE="<<_mc_dEtpc[0]<<"  <x>="<<_mc_dEpromx[0]<<"  <y>="<<_mc_dEpromy[0]<<"  <z>="<<_mc_dEpromz[0]<<" SpX="<<_mc_dEspreadx[0]<<" SpY="<<_mc_dEspready[0]<<" SpZ="<<_mc_dEspreadz[0]<<"\n";
      std::cout<<"*** TPC1   dE="<<_mc_dEtpc[1]<<"  <x>="<<_mc_dEpromx[1]<<"  <y>="<<_mc_dEpromy[1]<<"  <z>="<<_mc_dEpromz[1]<<" SpX="<<_mc_dEspreadx[1]<<" SpY="<<_mc_dEspready[1]<<" SpZ="<<_mc_dEspreadz[1]<<"\n\n";
    }

  }

  // --- Saving SimPhotons
  if(fSaveSimPhotons){
    
    ResetSimPhotons();
    
    if(fUseSimPhotonsLite){
      
      std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> fLitePhotonHandle_list;
      fLitePhotonHandle_list = e.getMany<std::vector<sim::SimPhotonsLite>>();
      if( fLitePhotonHandle_list.size() == 0 ){
        std::cout << "[Error] No available SimPhotonsLite" << std::endl;
        throw std::exception();
      }
      else{
        //Fill SimPhotons variables, loop over the Handle
        FillSimPhotonsLite(fLitePhotonHandle_list);
      }

    }
    else{

      // Handle for SimPhotons
      std::vector<art::Handle<std::vector<sim::SimPhotons> >> fPhotonHandle_list;
      fPhotonHandle_list = e.getMany<std::vector<sim::SimPhotons>>();
      if( fPhotonHandle_list.size() == 0 ){
        std::cout << "[Error] No available SimPhotons" << std::endl;
        throw std::exception();
      }
      else
        FillSimPhotons(fPhotonHandle_list);
    }

    //Fill per trackID tree with SimPhotons
    if(fMakePerTrackTree && fUseSimPhotonsLite==false){
      for(auto &phTrack:fSimPhotonsVUVMap){
        int trkID = phTrack.first;
        if(fSimPhotonsVUVMap[trkID].size()!=0 || fSimPhotonsVISMap[trkID].size()!=0){
            _simPhotonsVUV.clear();
            _simPhotonsVIS.clear();
            _simPhotonsVUV = fSimPhotonsVUVMap[trkID];
            _simPhotonsVIS = fSimPhotonsVISMap[trkID];
            _perTrackID=trkID;
            fPerTrackTree->Fill();
        }
      }
    }
  }

  // --- Saving raw waveforms
  if(fSaveRawWaveforms){

    if(fVerbosity>0)
      std::cout << "Saving raw waveforms from " << fRawWaveformsModuleLabel << std::endl;

    art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
    e.getByLabel(fRawWaveformsModuleLabel, wvfHandle);
    if(!wvfHandle.isValid()){
      std::cout << "RawWaveform with label " << fRawWaveformsModuleLabel << " not found..." << std::endl;
      throw std::exception();
    }

    _signalsDigi.clear();
    _stampTime.clear();
    _opChDigi.clear();
    
    for(auto const& wvf : (*wvfHandle)) {
      int fChNumber = wvf.ChannelNumber();
      if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(fChNumber) ) != fPDTypes.end() ){
        double t0_digi = wvf.TimeStamp();
        _stampTime.push_back(t0_digi);//time stamp in us
        _opChDigi.push_back(fChNumber);
        _signalsDigi.push_back({});

        for(unsigned int i=0;i<wvf.size();i++){
          _signalsDigi[_signalsDigi.size()-1].push_back(wvf[i]);
        }

      }
    }
  }

  //Saving Deconvolved Signals
  if(fSaveDeconvolvedWaveforms){

    if(fVerbosity>0)
      std::cout << "Saving deconvolved waveforms from " << fDeconvolvedWaveformsModuleLabel << std::endl;

    art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
    e.getByLabel(fDeconvolvedWaveformsModuleLabel, wvfHandle);
    if(!wvfHandle.isValid()){
      std::cout << "RawWaveform with label " << fRawWaveformsModuleLabel << " not found..." << std::endl;
      throw std::exception();
    }

    _signalsDeco.clear();
    _stampTimeDeco.clear(); 
    _opChDeco.clear();

    for(auto const& wvf : (*wvfHandle)) {
      
      int fChNumber = wvf.ChannelNumber();
      
      if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(fChNumber) ) != fPDTypes.end() ){
        double t0_Deco = wvf.TimeStamp();
        _stampTimeDeco.push_back(t0_Deco);//time stamp in us
        _opChDeco.push_back(fChNumber);
        _signalsDeco.push_back({});
        for(unsigned int i=0;i<wvf.size();i++){
          _signalsDeco[_signalsDeco.size()-1].push_back(wvf[i]);
        }
      }

    }

  }

  // --- Saving all OpHits
  if(fSaveOpHits){
  
    _nophits = 0;
    _ophit_opch.clear();
    _ophit_peakT.clear();
    _ophit_startT.clear();
    _ophit_riseT.clear();
    _ophit_width.clear();
    _ophit_area.clear();
    _ophit_amplitude.clear();
    _ophit_pe.clear();

    art::Handle< std::vector<recob::OpHit> > ophitListHandle;
    std::vector<art::Ptr<recob::OpHit> > ophitlist;

    for (size_t s = 0; s < fOpHitsModuleLabel.size(); s++) {

      e.getByLabel(fOpHitsModuleLabel[s], ophitListHandle);
      if(!ophitListHandle.isValid()){
        std::cout << "OpHit with label " << fOpHitsModuleLabel[s] << " not found..." << std::endl;
        throw std::exception();
      }
      
      if(fVerbosity>0)
        std::cout << "Saving OpHits from " << fOpHitsModuleLabel[s] << std::endl;

      art::fill_ptr_vector(ophitlist, ophitListHandle);
      _nophits += ophitlist.size();

      for (int i = 0; i < _nophits; ++i) {
        _ophit_opch.push_back( ophitlist.at(i)->OpChannel() );
        _ophit_peakT.push_back( ophitlist.at(i)->PeakTimeAbs() );
        _ophit_startT.push_back( ophitlist.at(i)->StartTime() );
        _ophit_riseT.push_back( ophitlist.at(i)->RiseTime() );
        _ophit_width.push_back( ophitlist.at(i)->Width() );
        _ophit_area.push_back( ophitlist.at(i)->Area() );
        _ophit_amplitude.push_back(  ophitlist.at(i)->Amplitude() );
        _ophit_pe.push_back( ophitlist.at(i)->PE() );
      }
    }

  }

  // --- Saving OpFlashes
  if(fSaveOpFlashes){

    _nopflash=0;
    _flash_id.clear();
    _flash_time.clear();
    _flash_total_pe.clear();
    _flash_pe_v.clear();
    _flash_tpc.clear();
    _flash_y.clear();
    _flash_yerr.clear();
    _flash_z.clear();
    _flash_zerr.clear();
    _flash_x.clear();
    _flash_xerr.clear();
    _flash_ophit_time.clear();
    _flash_ophit_risetime.clear();
    _flash_ophit_starttime.clear();
    _flash_ophit_amp.clear();
    _flash_ophit_area.clear();
    _flash_ophit_width.clear();
    _flash_ophit_pe.clear();
    _flash_ophit_ch.clear();

    art::Handle< std::vector<recob::OpFlash> > opflashListHandle;

    // Loop over all the OpFlash labels
    for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) {
      
      e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);
      if(!opflashListHandle.isValid()){
        std::cout << "OpFlash with label " << fOpFlashesModuleLabel[s] << " not found..." << std::endl;
        throw std::exception();
      }
      art::FindManyP<recob::OpHit> flashToOpHitAssns(opflashListHandle, e, fOpFlashesModuleLabel[s]);

      if(fVerbosity>0)
        std::cout << "Saving OpFlashes from " << fOpFlashesModuleLabel[s] << std::endl;

      for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
        // Get OpFlash
        art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
        recob::OpFlash Flash = *FlashPtr;

        if(fVerbosity>0)
          std::cout << "  *  " << _nopflash << " Time [ns]=" << 1000*Flash.AbsTime() << " PE=" << Flash.TotalPE() << std::endl;

        _flash_id.push_back( _nopflash );
        _flash_time.push_back( Flash.AbsTime() );
        _flash_total_pe.push_back( Flash.TotalPE() );
        _flash_pe_v.push_back( Flash.PEs() );
        _flash_tpc.push_back( s );
        _flash_y.push_back( Flash.YCenter() );
        _flash_yerr.push_back( Flash.YWidth() );
        _flash_x.push_back( Flash.XCenter() );
        _flash_xerr.push_back( Flash.XWidth() );
        _flash_z.push_back( Flash.ZCenter() );
        _flash_zerr.push_back( Flash.ZWidth() );
        _nopflash++;

        if(fSaveOpHits){
          
          _flash_ophit_time.push_back({});
          _flash_ophit_risetime.push_back({});
          _flash_ophit_starttime.push_back({});
          _flash_ophit_amp.push_back({});
          _flash_ophit_area.push_back({});
          _flash_ophit_width.push_back({});
          _flash_ophit_pe.push_back({});
          _flash_ophit_ch.push_back({});
          
          std::vector<art::Ptr<recob::OpHit>> ophit_v = flashToOpHitAssns.at(i);
          for (auto ophit : ophit_v) {
            _flash_ophit_time[_nopflash-1].push_back(ophit->PeakTimeAbs());
            _flash_ophit_risetime[_nopflash-1].push_back(ophit->RiseTime());
            _flash_ophit_starttime[_nopflash-1].push_back(ophit->StartTime());
            _flash_ophit_amp[_nopflash-1].push_back(ophit->Amplitude());
            _flash_ophit_area[_nopflash-1].push_back(ophit->Area());
            _flash_ophit_width[_nopflash-1].push_back(ophit->Width());
            _flash_ophit_pe[_nopflash-1].push_back(ophit->PE());
            _flash_ophit_ch[_nopflash-1].push_back(ophit->OpChannel());
          }

        }

      }
    
    }

  }

  // --- Saving CosmicID
  if(fSaveCosmicId){
    
    // Reset variables
    _CRUMBSScore.clear();
    _opT0Score.clear();
    _opT0Chi2.clear();
    _opT0Time.clear();
    _opT0MeasPE.clear();
    _opT0HypoPE.clear();
    _sliceOrigin.clear();
    _sliceCompleteness.clear();
    _slicePurity.clear();
    _sFMSliceOrigin.clear();
    _sFMScore.clear();
    _sFMTime.clear();
    _sFMScoreY.clear();
    _sFMScoreZ.clear();
    _sFMScoreRR.clear();
    _sFMScoreRatio.clear();
    _sFMScorePEtoQ.clear();
    _sFMScoreSlope.clear();
    _sFMPE.clear();

    // --- Read all hits
    ::art::Handle<std::vector<recob::Hit>> hitsHandle;
    e.getByLabel(fHitsLabel, hitsHandle);
    std::vector<art::Ptr<recob::Hit>> allHits;
    art::fill_ptr_vector(allHits, hitsHandle);
    std::map<std::string, int> allHitsTruthMap = GetAllHitsTruthMatch(e, allHits);

    // --- Read Recob Slice
    ::art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fReco2Label, sliceHandle);
    // Slice to CRUMBS
    art::FindOneP<sbn::CRUMBSResult> slice_crumbs_assns(sliceHandle, e, fCosmicIdModuleLabel);
    // Slice to OpT0Finder
    art::FindManyP<sbn::OpT0Finder> slice_opt0finder_assns(sliceHandle, e, fOpT0FinderModuleLabel);
    // Slice to hits
    art::FindManyP<recob::Hit> slice_hit_assns (sliceHandle, e, fReco2Label);
    
    // --- Store candidate slices
    std::vector< art::Ptr<recob::Slice> > sliceVect;
    art::fill_ptr_vector(sliceVect, sliceHandle);

    // --- Get the candidate slices
    for(size_t ix=0; ix<sliceVect.size(); ix++){

      auto & slice = sliceVect[ix];

      // --- Get the associated hits to the clusters in the slice
      std::vector< art::Ptr<recob::Hit> > sliceRecobHits = slice_hit_assns.at( slice.key() );
      
      if(fVerbosity>0)
        std::cout<<"  --- Slice: "<<ix<<" with "<<sliceRecobHits.size()<<" hits"<<std::endl;

      // --- Get the CRUMBS score
      const sbn::CRUMBSResult *slcCRUMBS = slice_crumbs_assns.isValid()? slice_crumbs_assns.at( slice.key() ).get():nullptr;
      if(slcCRUMBS){
        if(fVerbosity>0)
          std::cout<<"    * CRUMBS: "<<slcCRUMBS->score<<std::endl;
        _CRUMBSScore.push_back(slcCRUMBS->score);
      }
      else{
        _CRUMBSScore.push_back(-1);
      }

      // --- Get the OpT0Finder score
      const std::vector< art::Ptr<sbn::OpT0Finder> > slcOpT0Finder = slice_opt0finder_assns.at( slice.key() );
      if(slcOpT0Finder.size()>0){
        
        for(size_t jx=0; jx<slcOpT0Finder.size(); jx++){
          _opT0Score.push_back(slcOpT0Finder[jx]->score);
          _opT0Time.push_back(slcOpT0Finder[jx]->time);
          _opT0Chi2.push_back(1./slcOpT0Finder[jx]->score);
          _opT0MeasPE.push_back(slcOpT0Finder[jx]->measPE);
          _opT0HypoPE.push_back(slcOpT0Finder[jx]->hypoPE);

          // --- Truth matching
          double purity = 0;
          double completeness = 0;
          simb::MCParticle mainMCParticle;
          simb::MCTruth matchedTruth = MakeTruthMatching(e, sliceRecobHits, allHitsTruthMap, slcOpT0Finder[jx]->tpc, purity, completeness, mainMCParticle);

          // --- Store the slice origin
          _sliceOrigin.push_back(matchedTruth.Origin());
          _sliceCompleteness.push_back(completeness);
          _slicePurity.push_back(purity);

          if(fVerbosity>0){
            std::cout<<"    * OpT0Finder: "<<slcOpT0Finder[jx]->score<<std::endl;
            std::cout<<"    * Truth matching: "<<matchedTruth.Origin()<<" "<<matchedTruth.GetNeutrino().CCNC()<<" "<<purity<<" "<<completeness<<std::endl;
          }
        }
      }

      // --- Save SimpleFlash
      // Slice to PFParticle
      art::FindManyP<recob::PFParticle> slice_pfps_assns (sliceHandle, e, fReco2Label);
      std::vector< art::Ptr<recob::PFParticle> > slicePFParticles = slice_pfps_assns.at( slice.key() );
      // PFP handle
      art::Handle<std::vector<recob::PFParticle>> pfpHandle;
      e.getByLabel(fReco2Label, pfpHandle);
      // PFP to SimpleFlash
      art::FindManyP<sbn::SimpleFlashMatch> pfp_simpleflash_assns (pfpHandle, e, fSimpleFlashMatchModuleLabel);

      for(size_t jx=0; jx<slicePFParticles.size(); jx++){
        std::vector< art::Ptr<sbn::SimpleFlashMatch> > pfpSimpleFlash = pfp_simpleflash_assns.at( slicePFParticles[jx].key() );
        
        if(pfpSimpleFlash.size()>0){

          if(fVerbosity>0)
            std::cout<<"    * SimpleFlashMatch  PFParticle "<<jx<<" PDG"<<slicePFParticles[jx]->PdgCode()<<" N FMatches: "<<pfpSimpleFlash.size()<<std::endl;
              
          for(size_t kx=0; kx<pfpSimpleFlash.size(); kx++){
        
            // --- Truth matching
            double purity = 0;
            double completeness = 0;
            simb::MCParticle mainMCParticle;
            int SFMTPC = 1;
              if(pfpSimpleFlash[kx]->charge.center.X()<0){
                SFMTPC = 0;
            }
            simb::MCTruth matchedTruth = MakeTruthMatching(e, sliceRecobHits, allHitsTruthMap, SFMTPC, purity, completeness, mainMCParticle);
            int sfmOrigin = matchedTruth.Origin();
            
            _sFMSliceOrigin.push_back(sfmOrigin);
            _sFMScore.push_back(pfpSimpleFlash[kx]->score.total);
            _sFMTime.push_back(pfpSimpleFlash[kx]->time);
            _sFMScoreY.push_back(pfpSimpleFlash[kx]->score.y);
            _sFMScoreZ.push_back(pfpSimpleFlash[kx]->score.z);
            _sFMScoreRR.push_back(pfpSimpleFlash[kx]->score.rr);
            _sFMScoreRatio.push_back(pfpSimpleFlash[kx]->score.ratio);
            _sFMScorePEtoQ.push_back(pfpSimpleFlash[kx]->score.petoq);
            _sFMScoreSlope.push_back(pfpSimpleFlash[kx]->score.slope);
            _sFMPE.push_back(pfpSimpleFlash[kx]->light.pe);

            if(fVerbosity>0){
              std::cout<<"      -- Truth matching: "<<sfmOrigin<<" "<<matchedTruth.GetNeutrino().CCNC()<<" "<<purity<<" "<<completeness<<std::endl;
              std::cout<<"      -- SimpleFlashMatch: "<<kx<<" Score: "<<pfpSimpleFlash[kx]->score.total<<std::endl;
            }
          }
        }
      }

    }
  }

  // --- Fill the tree
  fTree->Fill();
}


// -------- Function to fill the MCTruth information --------
void opdet::SBNDPDSAnalyzer::FillMCTruth(art::Event const& e){
  
  
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

    if(fVerbosity>0){
      std::cout <<"Saving MCTruth from "<<fMCTruthModuleLabel[s]<<" with instance "<<fMCTruthInstanceLabel[s];
      std::cout << " with " << mctruth_v.size() << " MCTruths." << std::endl;
    }

    for (size_t n = 0; n < mctruth_v.size(); n++) {

      art::Ptr<simb::MCTruth> evtTruth = mctruth_v[n];

      if(fVerbosity>0){
        std::cout << "  Origin: " << evtTruth->Origin() << std::endl;
        std::cout << "  We have " << evtTruth->NParticles() << " particles." << std::endl;
        std::cout << "  Mode=" << evtTruth->GetNeutrino().Mode() <<"  IntType="<<evtTruth->GetNeutrino().InteractionType();
        std::cout << "  Target=" << evtTruth->GetNeutrino().Target()<<" CCNC=" << evtTruth->GetNeutrino().CCNC()<<std::endl;
      }

      double nu_x, nu_y, nu_z, nu_t, nu_E;

      // Loop over particles
      for (int p = 0; p < evtTruth->NParticles(); p++){

        simb::MCParticle const& par = evtTruth->GetParticle(p);

        // Only save MCTruth if the origins is specified in the fhicl list
        if( find (fMCTruthOrigin.begin(), fMCTruthOrigin.end(), evtTruth->Origin() ) != fMCTruthOrigin.end() ){

          if(fVerbosity>0){
            std::cout << "    " << par.TrackId() << "  Particle PDG: " << par.PdgCode() << " E: " << par.E() << " t: " << par.T();
            std::cout << " Mother: "<<par.Mother() << " Process: "<<par.Process()<<" Status: "<<par.StatusCode()<<std::endl;
          }
          
          // Only save vertex if the PDG is specified in the fhicl list
          if( find (fMCTruthPDG.begin(), fMCTruthPDG.end(), par.PdgCode() ) != fMCTruthPDG.end() ){
            // For BNB neutinos
            if(par.StatusCode()==0 && evtTruth->Origin()==1){
              nu_x=par.Vx();
              nu_y=par.Vy();
              nu_z=par.Vz();
              nu_t=par.T();
              nu_E=par.E();
            }
            // For single particle gun
            else if(par.StatusCode()==1 && evtTruth->Origin()==4){
              nu_x=par.Vx();
              nu_y=par.Vy();
              nu_z=par.Vz();
              nu_t=par.T();
              nu_E=par.E();
            }
          }

        }

      }

      if(fVerbosity>0)
        std::cout << "     Vertex: " << nu_x << " " << nu_y << " " << nu_z << " T: " << nu_t << " E: " << nu_E << std::endl;

      _nuvT.push_back(nu_t);
      _nuvX.push_back(nu_x);
      _nuvY.push_back(nu_y);
      _nuvZ.push_back(nu_z);
      _nuvE.push_back(nu_E);

    }
  }
}


// -------- Function to fill the average energy deposition --------
void opdet::SBNDPDSAnalyzer::FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
  std::vector<std::vector<double>> fenergydepY, std::vector<std::vector<double>> fenergydepZ, std::vector<std::vector<double>> fstepT,
  std::vector<double> &dEtpc, std::vector<double> &dEpromx, std::vector<double> &dEpromy, std::vector<double> &dEpromz,
  std::vector<double> &dEspreadx, std::vector<double> &dEspready, std::vector<double> &dEspreadz,
  std::vector<std::vector<double>> &dElowedges, std::vector<std::vector<double>> &dEmaxedges)
{

  // Initialize variables
  dEtpc.clear(); dEpromx.clear(); dEpromy.clear(); dEpromz.clear();;
  dEtpc.resize(2, 0);
  dEpromx.resize(2, fDefaultSimIDE); dEpromy.resize(2, fDefaultSimIDE); dEpromz.resize(2, fDefaultSimIDE);
  dEspreadx.clear(); dEspready.clear(); dEspreadz.clear();
  dEspreadx.resize(2, fDefaultSimIDE); dEspready.resize(2, fDefaultSimIDE); dEspreadz.resize(2, fDefaultSimIDE);

  int ndeps_tpc0=0, ndeps_tpc1=0;
  double dEpromx_tpc0=0, dEpromy_tpc0=0, dEpromz_tpc0=0;
  double dEpromx_tpc1=0, dEpromy_tpc1=0, dEpromz_tpc1=0;
  double spreadx_tpc0=0, spreadx_tpc1=0;
  double spready_tpc0=0, spready_tpc1=0;
  double spreadz_tpc0=0, spreadz_tpc1=0;
  double minz_tpc0=1e3, maxz_tpc0=-1e3, minz_tpc1=1e3, maxz_tpc1=-1e3;
  double miny_tpc0=1e3, maxy_tpc0=-1e3, miny_tpc1=1e3, maxy_tpc1=-1e3;
  double minx_tpc0=1e3, maxx_tpc0=-1e3, minx_tpc1=1e3, maxx_tpc1=-1e3;
  double dE_tpc0=0, dE_tpc1=0;
  
  for(size_t k=0; k<fenergydep.size(); k++){
    for(size_t j=0; j<fenergydep.at(k).size(); j++){
      
      // TPC0
      if(fenergydepX.at(k).at(j)<0){

        ndeps_tpc0++;
        dE_tpc0+=fenergydep.at(k).at(j);
        dEpromx_tpc0+=fenergydep.at(k).at(j)*fenergydepX.at(k).at(j);
        dEpromy_tpc0+=fenergydep.at(k).at(j)*fenergydepY.at(k).at(j);
        dEpromz_tpc0+=fenergydep.at(k).at(j)*fenergydepZ.at(k).at(j);
        spreadx_tpc0+=fenergydep.at(k).at(j)*fenergydepX.at(k).at(j)*fenergydepX.at(k).at(j);
        spready_tpc0+=fenergydep.at(k).at(j)*fenergydepY.at(k).at(j)*fenergydepY.at(k).at(j);
        spreadz_tpc0+=fenergydep.at(k).at(j)*fenergydepZ.at(k).at(j)*fenergydepZ.at(k).at(j);
        
        // Update min/max values for TPC0
        minx_tpc0 = std::min(minx_tpc0, fenergydepX.at(k).at(j));
        maxx_tpc0 = std::max(maxx_tpc0, fenergydepX.at(k).at(j));
        miny_tpc0 = std::min(miny_tpc0, fenergydepY.at(k).at(j));
        maxy_tpc0 = std::max(maxy_tpc0, fenergydepY.at(k).at(j));
        minz_tpc0 = std::min(minz_tpc0, fenergydepZ.at(k).at(j));
        maxz_tpc0 = std::max(maxz_tpc0, fenergydepZ.at(k).at(j));
      }

      // TPC1
      else{
        ndeps_tpc1++;
        dE_tpc1+=fenergydep.at(k).at(j);
        dEpromx_tpc1+=fenergydep.at(k).at(j)*fenergydepX.at(k).at(j);
        dEpromy_tpc1+=fenergydep.at(k).at(j)*fenergydepY.at(k).at(j);
        dEpromz_tpc1+=fenergydep.at(k).at(j)*fenergydepZ.at(k).at(j);
        spreadx_tpc1+=fenergydep.at(k).at(j)*fenergydepX.at(k).at(j)*fenergydepX.at(k).at(j);
        spready_tpc1+=fenergydep.at(k).at(j)*fenergydepY.at(k).at(j)*fenergydepY.at(k).at(j);
        spreadz_tpc1+=fenergydep.at(k).at(j)*fenergydepZ.at(k).at(j)*fenergydepZ.at(k).at(j);
        
        // Update min/max values for TPC1
        minx_tpc1 = std::min(minx_tpc1, fenergydepX.at(k).at(j));
        maxx_tpc1 = std::max(maxx_tpc1, fenergydepX.at(k).at(j));
        miny_tpc1 = std::min(miny_tpc1, fenergydepY.at(k).at(j));
        maxy_tpc1 = std::max(maxy_tpc1, fenergydepY.at(k).at(j));
        minz_tpc1 = std::min(minz_tpc1, fenergydepZ.at(k).at(j));
        maxz_tpc1 = std::max(maxz_tpc1, fenergydepZ.at(k).at(j));
      }
    }
  }

  if(ndeps_tpc0!=0){
    dEpromx_tpc0=dEpromx_tpc0/dE_tpc0;
    dEpromy_tpc0=dEpromy_tpc0/dE_tpc0;
    dEpromz_tpc0=dEpromz_tpc0/dE_tpc0;
    spreadx_tpc0=sqrt( spreadx_tpc0/dE_tpc0-dEpromx_tpc0*dEpromx_tpc0 );
    spready_tpc0=std::sqrt(spready_tpc0/dE_tpc0-dEpromy_tpc0*dEpromy_tpc0);
    spreadz_tpc0=std::sqrt(spreadz_tpc0/dE_tpc0-dEpromz_tpc0*dEpromz_tpc0);
    dEtpc[0]=dE_tpc0;dEpromx[0]=dEpromx_tpc0;dEpromy[0]=dEpromy_tpc0;dEpromz[0]=dEpromz_tpc0;
    dEspreadx[0]=spreadx_tpc0;dEspready[0]=spready_tpc0;dEspreadz[0]=spreadz_tpc0;
  }
  if(ndeps_tpc1!=0){
    dEpromx_tpc1=dEpromx_tpc1/dE_tpc1;
    dEpromy_tpc1=dEpromy_tpc1/dE_tpc1;
    dEpromz_tpc1=dEpromz_tpc1/dE_tpc1;
    spreadx_tpc1=std::sqrt(spreadx_tpc1/dE_tpc1-dEpromx_tpc1*dEpromx_tpc1);
    spready_tpc1=std::sqrt(spready_tpc1/dE_tpc1-dEpromy_tpc1*dEpromy_tpc1);
    spreadz_tpc1=std::sqrt(spreadz_tpc1/dE_tpc1-dEpromz_tpc1*dEpromz_tpc1);
    dEtpc[1]=dE_tpc1;dEpromx[1]=dEpromx_tpc1;dEpromy[1]=dEpromy_tpc1;dEpromz[1]=dEpromz_tpc1;
    dEspreadx[1]=spreadx_tpc1;dEspready[1]=spready_tpc1;dEspreadz[1]=spreadz_tpc1;
  }

  dElowedges.clear(); dElowedges.resize(2);
  dEmaxedges.clear(); dEmaxedges.resize(2);
  dElowedges.at(0).push_back(minx_tpc0); dElowedges.at(0).push_back(miny_tpc0);dElowedges.at(0).push_back(minz_tpc0);
  dElowedges.at(1).push_back(minx_tpc1); dElowedges.at(1).push_back(miny_tpc1);dElowedges.at(1).push_back(minz_tpc1);
  dEmaxedges.at(0).push_back(maxx_tpc0);dEmaxedges.at(0).push_back(maxy_tpc0);dEmaxedges.at(0).push_back(maxz_tpc0);
  dEmaxedges.at(1).push_back(maxx_tpc1);dEmaxedges.at(1).push_back(maxy_tpc1);dEmaxedges.at(1).push_back(maxz_tpc1);

  return;
}

// -------- Restet SimPhotons variables --------
void opdet::SBNDPDSAnalyzer::ResetSimPhotons(){
    //Cleaning
    _simPhotonsLiteVIS.clear(); _simPhotonsLiteVIS.resize(fGeoService->NOpChannels());
    _simPhotonsLiteVUV.clear(); _simPhotonsLiteVUV.resize(fGeoService->NOpChannels());
    _simPhotonsperOpChVUV.clear(); _simPhotonsperOpChVUV.resize(fGeoService->NOpChannels(), 0);
    _simPhotonsperOpChVIS.clear(); _simPhotonsperOpChVIS.resize(fGeoService->NOpChannels(), 0);
    _NPhotons=0;
    _NPhotonsPMTCo=0; _NPhotonsPMTUnco=0; _NPhotonsPMTCoVUV=0;
    _NPhotonsXARAPUCAVUV=0; _NPhotonsXARAPUCAVIS=0;
}


// -------- Functoin to fill the SimPhotons --------
void opdet::SBNDPDSAnalyzer::FillSimPhotons(std::vector<art::Handle<std::vector<sim::SimPhotons> >> photonHandle_list){
  
  if(fVerbosity>0)
    std::cout << "Saving SimPhotons for NOpCh=" << fGeoService->NOpChannels() << std::endl;

  //Fill SimPhotons variables, loop over the Handle
  for ( const art::Handle<std::vector<sim::SimPhotons>>& photonHandle: (photonHandle_list) ){
    
    std::string spLabel = photonHandle.provenance()->moduleLabel();

    // Check if the module label is in the list of modules to save
    if(std::find(fSimPhotonsModuleLabel.begin(), fSimPhotonsModuleLabel.end(), spLabel)==fSimPhotonsModuleLabel.end()) continue;
    
    if(fVerbosity>0)
      std::cout<<"   Saving Handle: "<<photonHandle.provenance()->moduleLabel()<<" "<<photonHandle.provenance()->productInstanceName()<<std::endl;
   
    // Reflected light
    bool reflected = (photonHandle.provenance()->productInstanceName() == "Reflected");
    
    // Loop over the SimPhotons
    for ( auto const& itOpDet : (*photonHandle) ){

      int opch = itOpDet.OpChannel();
      int nphotons=0;

      std::string pd_type = fPDSMap.pdType(opch);
      // Channels not sensitive to reflected light
      if(reflected && pd_type=="xarapuca_vuv") continue;
      // Channels not sensitive to direct light
      if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;
      
      // Only save the PD types specified in the fhicl list
      if(std::find(fPDTypes.begin(), fPDTypes.end(), pd_type ) != fPDTypes.end() ){

        const sim::SimPhotons& TheHit = itOpDet;
        for (const sim::OnePhoton& Phot : TheHit) {
          
          // Save times
          if(fSaveSimPhotonsArrivalTimes){
            if(reflected){
              _simPhotonsLiteVIS.at(opch).push_back(Phot.Time);
            }
            else{
              _simPhotonsLiteVUV.at(opch).push_back(Phot.Time);
            }
          }

          // Fill SimPhotons-TrackID maps
          if(fMakePerTrackTree){
            if(reflected){
              fSimPhotonsVISMap[std::abs(Phot.MotherTrackID)].at(opch).push_back(Phot.Time);
            }
            else{
              fSimPhotonsVUVMap[std::abs(Phot.MotherTrackID)].at(opch).push_back(Phot.Time);
            }
          }

          nphotons++;
        }

        //Notice that with the new LArG4 there are 4 SimPhotons instance names, so we have to add every instance (hybrid model)
        if(reflected){
          _simPhotonsperOpChVIS.at(opch)+=nphotons;
        }
        else{
          _simPhotonsperOpChVUV.at(opch)+=nphotons;
        }

        _NPhotons+=nphotons;
        if(pd_type=="pmt_coated") {
          _NPhotonsPMTCo+=nphotons;
          if(!reflected) _NPhotonsPMTCoVUV+=nphotons;
        }
        else if(pd_type=="pmt_uncoated") _NPhotonsPMTUnco+=nphotons;
        else if(pd_type=="xarapuca_vuv") _NPhotonsXARAPUCAVUV+=nphotons;
        else if(pd_type=="xarapuca_vis") _NPhotonsXARAPUCAVIS+=nphotons;
      }
    }
  }
}


// -------- Function to fill the SimPhotonsLite --------
void opdet::SBNDPDSAnalyzer::FillSimPhotonsLite(std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> photonHandle_list){
  
  if(fVerbosity>0)
    std::cout << "Saving SimPhotonsLite for NOpCh=" << fGeoService->NOpChannels() << std::endl;
  
  //Fill SimPhotons variables, loop over the Handle
  for ( const art::Handle<std::vector<sim::SimPhotonsLite>>& litePhotonHandle: (photonHandle_list) ){

    std::string spLabel = litePhotonHandle.provenance()->moduleLabel();

    // Check if the module label is in the list of modules to save
    if(std::find(fSimPhotonsModuleLabel.begin(), fSimPhotonsModuleLabel.end(), spLabel)==fSimPhotonsModuleLabel.end()) continue;
    
    if(fVerbosity>0)
      std::cout<<"   Saving Handle: "<<litePhotonHandle.provenance()->moduleLabel()<<" "<<litePhotonHandle.provenance()->productInstanceName()<<std::endl;
   
    // Reflected light
    bool reflected = (litePhotonHandle.provenance()->productInstanceName() == "Reflected");

    // Loop over the SimPhotonsLite
    for ( auto const& fLitePhotons : (*litePhotonHandle) ){
      
      int opch=fLitePhotons.OpChannel;

      std::string pd_type = fPDSMap.pdType(opch);
      // Channels not sensitive to reflected light
      if(reflected && pd_type=="xarapuca_vuv") continue;
      // Channels not sensitive to direct light
      if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;
      
      // Only save the PD types specified in the fhicl list
      if(std::find(fPDTypes.begin(), fPDTypes.end(), pd_type ) != fPDTypes.end() ){
        
        std::map<int, int> fLitePhotons_map = fLitePhotons.DetectedPhotons;
        int nphotons=0;
        
        for(auto fphoton = fLitePhotons_map.begin(); fphoton!= fLitePhotons_map.end(); fphoton++){

          nphotons+=fphoton->second;

          if(fSaveSimPhotonsArrivalTimes){
            for(int i = 0; i < fphoton->second ; i++) {
              if(reflected){
                _simPhotonsLiteVIS.at(opch).push_back(fphoton->first);
              }
              else{
                _simPhotonsLiteVUV.at(opch).push_back(fphoton->first);
              }
            }
          }
        }

        // Fill #photons per OpChannel
        if(reflected){
          _simPhotonsperOpChVIS.at(opch)+=nphotons;
        }
        else{
          _simPhotonsperOpChVUV.at(opch)+=nphotons;
        }

        // Total number of photons
        _NPhotons+=nphotons;

        if(pd_type=="pmt_coated") {
          _NPhotonsPMTCo+=nphotons;
          if(!reflected) _NPhotonsPMTCoVUV+=nphotons;
        }
        else if(pd_type=="pmt_uncoated") _NPhotonsPMTUnco+=nphotons;
        else if(pd_type=="xarapuca_vuv") _NPhotonsXARAPUCAVUV+=nphotons;
        else if(pd_type=="xarapuca_vis") _NPhotonsXARAPUCAVIS+=nphotons;

      }
    }
  }

}


// -------- Function to make the truth matching --------
simb::MCTruth opdet::SBNDPDSAnalyzer::MakeTruthMatching(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &recoHits, std::map<std::string, int>& allHitsTruthMap, unsigned int TPC, double& purity, double& completeness, simb::MCParticle &mainMCParticle)
{
    art::ServiceHandle<detinfo::DetectorClocksService> timeservice;
    auto const clockData(timeservice->DataFor(e));
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    std::map<std::string, int> truthCounter;
    std::map<std::string, int> truthCounterTPC;
    std::map<std::string, simb::MCTruth > truthMap;
    std::map<int, int> trackIdCounter;
    simb::MCTruth truth;
    
    
    for (auto const& hit : recoHits){
      int trkID = TruthMatchUtils::TrueParticleID(clockData,hit,true);

      if ( trkID != std::numeric_limits<int>::max() && trkID != std::numeric_limits<int>::min()) {
        art::Ptr<simb::MCTruth> truthPtr = pi_serv->TrackIdToMCTruth_P(trkID);
        truth = *truthPtr;
        size_t truthKey = truthPtr.key();
        int truthOrigin = truth.Origin();
        //std::cout<<"Hit "<<hit->WireID()<<" matched to "<<trkID<<" Origin: "<<truthOrigin<<" "<<truthKey<<std::endl;
        std::string label = std::to_string(truthOrigin)+"_"+std::to_string(truthKey);
        truthCounter[label]++;
        truthMap[label] = truth;
        if(hit->WireID().TPC==TPC)
          truthCounterTPC[label]++;
        trackIdCounter[trkID]++;
      }
      
    }

    // for completeness and purity calculations
    int maxHits = 0;
    std::string labelMatched = "";
    for (auto const& [label, nHits] : truthCounter){
      if(nHits>maxHits){
        maxHits = nHits;
        labelMatched = label;
      }
    }

    if(labelMatched!=""){
      // calculate purity
      purity = (double)maxHits/(double)recoHits.size();

      // calculate completeness
      int nHitsTruth = allHitsTruthMap[labelMatched];
      completeness = (double)maxHits/(double)nHitsTruth;
    }
    else{
      purity = 0.;
      completeness = 0.;
    }


    // return truth with max number of hits in the TPC  
    maxHits = 0;
    labelMatched = "";
    for (auto const& [label, nHits] : truthCounterTPC){
      if(nHits>maxHits){
        maxHits = nHits;
        truth = truthMap[label];
        labelMatched = label;
      }
    }

    // assign the MCParticle with the most hits
    maxHits = 0;
    int trackIdMatched = -1;
    for (auto const& [trackId, nHits] : trackIdCounter){
      if(nHits>maxHits){
        maxHits = nHits;
        trackIdMatched = trackId;
      }
    }

    if(trackIdMatched!=-1){
      mainMCParticle = pi_serv->TrackIdToParticle(trackIdMatched);
    }

    return truth;
}


// -------- Function to return a map with the number of associated hits per MCTruth --------
std::map<std::string, int> opdet::SBNDPDSAnalyzer::GetAllHitsTruthMatch(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &allHits){
  
    art::ServiceHandle<detinfo::DetectorClocksService> timeservice;
    auto const clockData(timeservice->DataFor(e));
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // For completeness calculation
    std::map<std::string, int> allHitsTruthMap;
    art::Ptr<simb::MCTruth> truth;
    for (auto const& hit : allHits){
      int trkID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      if ( trkID != std::numeric_limits<int>::max() && trkID != std::numeric_limits<int>::min()) {
        truth = pi_serv->TrackIdToMCTruth_P(trkID);
        size_t truthKey = truth.key();
        int truthOrigin = truth->Origin();
        std::string label = std::to_string(truthOrigin)+"_"+std::to_string(truthKey);
        allHitsTruthMap[label]++;
      } 
    }

    return allHitsTruthMap;
}
