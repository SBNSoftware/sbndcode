#include "SBNDPDSProducer_module.hh"

#include "larcorealg/Geometry/OpDetGeo.h"

// -------- Constructor --------
opdet::SBNDPDSProducer::SBNDPDSProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fMCTruthOrigin( p.get<std::vector<int>>("MCTruthOrigin") ),
    fMCTruthPDG( p.get<std::vector<int>>("MCTruthPDG") ),
    fMCTruthModuleLabel( p.get<std::vector<std::string>>("MCTruthModuleLabel") ),
    fMCTruthInstanceLabel( p.get<std::vector<std::string>>("MCTruthInstanceLabel") ),
    fMCModuleLabel( p.get<std::string>("MCModuleLabel") ),
    fOpHitsModuleLabel( p.get<std::vector<std::string>>("OpHitsModuleLabel") ),
    fOpFlashesModuleLabel( p.get<std::vector<std::string>>("OpFlashesModuleLabel") ),
    fG4BufferBoxX( p.get<std::vector<int>>("G4BufferBoxX") ),
    fG4BufferBoxY( p.get<std::vector<int>>("G4BufferBoxY") ),
    fG4BufferBoxZ( p.get<std::vector<int>>("G4BufferBoxZ") ),
    fG4BeamWindow( p.get<std::vector<int>>("G4BeamWindow") ),
    fKeepPDGCode( p.get<std::vector<int>>("KeepPDGCode", {}) ),
    fSaveOpHits( p.get<bool>("SaveOpHits", true) ),
    fVerbosity( p.get<int>("Verbosity") ),
    fCoatedPMTMapPath( p.get<std::string>("CoatedPMTMapPath", "/exp/sbnd/app/users/svidales/larsoft_v10_06_01_develop/srcs/sbndcode/sbndcode/PosRecoCVN/pmt_maps/coated_pmt_map_realistic_flipped.csv") ),
    fUncoatedPMTMapPath( p.get<std::string>("UncoatedPMTMapPath", "/exp/sbnd/app/users/svidales/larsoft_v10_06_01_develop/srcs/sbndcode/sbndcode/PosRecoCVN/pmt_maps/uncoated_pmt_map_realistic_flipped.csv") ),
    dE_neutrinowindow( 0.0 ) 
{
    produces<PixelMapVars>();
    
    // Initialize channel dictionary
    InitializeChannelDict();
    
    // Classify channels by type and parity
    ClassifyChannels();
    
    // Load PMT maps
    LoadPMTMaps();
}


// -------- Main function --------
void opdet::SBNDPDSProducer::produce(art::Event& e)
{

  // --- Event General Info
  _eventID = e.id().event();
  _runID = e.id().run();
  _subrunID = e.id().subRun();
  if(fVerbosity>0)
    std::cout << " -- Running SBNDPDSProducer -- \n Run=" << _runID << " Subrun=" << _subrunID << " Event=" << _eventID << std::endl;

  // --- Services
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<detinfo::DetectorClocksService> timeservice;

  // --- Saving MCTruths
  _nuvT.clear(); _nuvX.clear(); _nuvY.clear(); _nuvZ.clear(); _nuvE.clear();
  FillMCTruth(e);

  // --- Load MCParticles from event
  mcpartVec.clear();
  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
  e.getByLabel(fMCModuleLabel, mcParticleHandle);
  if(mcParticleHandle.isValid()){
    mcpartVec = *mcParticleHandle;
    if(fVerbosity>0)
      std::cout << "Loaded " << mcpartVec.size() << " MCParticles from " << fMCModuleLabel << std::endl;
  } else {
    if(fVerbosity>0)
      std::cout << "MCParticles with label " << fMCModuleLabel << " not found. mcpartVec will be empty." << std::endl;
  }

  // --- Saving MCParticles
  _mc_stepX.clear(); _mc_stepY.clear(); _mc_stepZ.clear(); _mc_stepT.clear();
  _mc_dE.clear(); _mc_E.clear();
  _mc_trackID.clear(); _mc_motherID.clear(); _mc_PDGcode.clear(); _mc_process.clear();
  _mc_StartPx.clear(); _mc_StartPy.clear(); _mc_StartPz.clear();
  _mc_EndPx.clear(); _mc_EndPy.clear(); _mc_EndPz.clear();
  _mc_energydep.clear(); _mc_energydepX.clear(); _mc_energydepY.clear(); _mc_energydepZ.clear();
  _mc_InTimeCosmics=0; _mc_InTimeCosmicsTime.clear();
  dE_neutrinowindow = 0.0;

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
    if(fVerbosity>1){
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

  // --- Saving all OpHits
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
      std::string msg = "OpHit with label " + fOpHitsModuleLabel[s] + " not found in event " + std::to_string(_eventID);
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    
    if(fVerbosity>0)
      std::cout << "Saving OpHits from " << fOpHitsModuleLabel[s] << std::endl;

    art::fill_ptr_vector(ophitlist, ophitListHandle);
    _nophits += ophitlist.size();

    for (size_t i = 0; i < ophitlist.size(); ++i) {
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

  // --- Saving OpFlashes
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
      std::string msg = "OpFlash with label " + fOpFlashesModuleLabel[s] + " not found in event " + std::to_string(_eventID);
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    art::FindManyP<recob::OpHit> flashToOpHitAssns(opflashListHandle, e, fOpFlashesModuleLabel[s]);

    if(fVerbosity>1)
      std::cout << "Saving OpFlashes from " << fOpFlashesModuleLabel[s] << std::endl;

    for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
      // Get OpFlash
      art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
      recob::OpFlash Flash = *FlashPtr;

      // Filter cosmic rays: only keep flashes in beam window (0.367 - 1.9 μs)
      if(Flash.AbsTime() < 0.367 || Flash.AbsTime() > 1.9) continue;

      if(fVerbosity>1)
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
        
        // Solo guardar PE, Channel y Time para mapas del detector
        _flash_ophit_time.push_back({});
        _flash_ophit_pe.push_back({});
        _flash_ophit_ch.push_back({});
        
        std::vector<art::Ptr<recob::OpHit>> ophit_v = flashToOpHitAssns.at(i);
        for (auto ophit : ophit_v) {
          _flash_ophit_time[_nopflash-1].push_back(ophit->PeakTimeAbs());
          _flash_ophit_pe[_nopflash-1].push_back(ophit->PE());
          _flash_ophit_ch[_nopflash-1].push_back(ophit->OpChannel());
        }

      }

    }
  
  }

  // Apply first filter: only keep events where nuvT has exactly one element
  bool passFilter1 = (_nuvT.size() == 1);
  
  // Apply second filter: only keep events with at least one flash (flash_ophit_pe has at least one element)
  bool passFilter2 = false;
  if (!_flash_ophit_pe.empty()) {
    // Check if there's at least one flash with optical hits
    for (const auto& flash_pe : _flash_ophit_pe) {
      if (!flash_pe.empty()) {
        passFilter2 = true;
        break;
      }
    }
  }
  
  bool passFilter = passFilter1 && passFilter2;
  
  if(fVerbosity > 0) {
    std::cout << "Filter check: nuvT size = " << _nuvT.size() << ", pass filter 1 = " << passFilter1 << std::endl;
    std::cout << "Filter check: flash count = " << _flash_ophit_pe.size() << ", has ophits = " << passFilter2 << ", pass filter 2 = " << passFilter2 << std::endl;
    std::cout << "Overall filter result = " << passFilter << std::endl;
  }
  
  // Apply flash selection logic if event passes filters
  if(passFilter) {
    ApplyFlashSelection();
    ApplyFinalEnergyFilter();
    
    // Only create PE matrix and images if we have final filtered data
    if (!_flash_ophit_pe_final.empty()) {
      CreatePEMatrix();
      CreatePEImages();
    }
  }

  // Crear y llenar el producto PixelMapVars
  auto pixelVars = std::make_unique<PixelMapVars>();
  
  if(passFilter) {
    // Apply mask to all variables - use final filtered data after all cuts
    pixelVars->flash_ophit_pe = _flash_ophit_pe_final;
    pixelVars->flash_ophit_ch = _flash_ophit_ch_final;
    pixelVars->flash_ophit_time = _flash_ophit_time_final;
    pixelVars->nuvT = _nuvT_final;
    pixelVars->dEpromx = _mc_dEpromx_final;
    pixelVars->dEpromy = _mc_dEpromy_final;
    pixelVars->dEpromz = _mc_dEpromz_final;
    pixelVars->dEtpc = _mc_dEtpc_final;
    pixelVars->nuvZ = _nuvZ_final;
    
    if(fVerbosity > 0) {
      std::cout << "Event passed filter - data stored in PixelMapVars" << std::endl;
    }
  } else {
    // Initialize empty vectors for events that don't pass filter
    pixelVars->flash_ophit_pe.clear();
    pixelVars->flash_ophit_ch.clear();
    pixelVars->flash_ophit_time.clear();
    pixelVars->nuvT.clear();
    pixelVars->dEpromx.clear();
    pixelVars->dEpromy.clear();
    pixelVars->dEpromz.clear();
    pixelVars->dEtpc.clear();
    pixelVars->nuvZ.clear();
    
    if(fVerbosity > 0) {
      std::cout << "Event failed filter - empty data stored in PixelMapVars" << std::endl;
    }
  }
  
  // Always include the channel dictionary (independent of filter)
  pixelVars->channel_dict = fChannelDict;
  
  // Only assign PE matrix and images if they were created (event passed all filters)
  pixelVars->pe_matrix = _pe_matrix;
  pixelVars->pe_images = _pe_images;
  pixelVars->coated_pmt_map = _coated_pmt_map;
  pixelVars->uncoated_pmt_map = _uncoated_pmt_map;

  if(fVerbosity > 0) {
    std::cout << "=== PixelMapVars Summary ===" << std::endl;
    std::cout << "flash_ophit_pe size: " << pixelVars->flash_ophit_pe.size() << std::endl;
    std::cout << "pe_matrix size: " << pixelVars->pe_matrix.size() << "x" 
              << (pixelVars->pe_matrix.empty() ? 0 : pixelVars->pe_matrix[0].size()) << std::endl;
    std::cout << "pe_images size: " << pixelVars->pe_images.size() << std::endl;
    std::cout << "Only PixelMapVars should be in output ROOT file" << std::endl;
  }

  // Poner el producto en el evento
  e.put(std::move(pixelVars));

}


// -------- Create PE Matrix --------
void opdet::SBNDPDSProducer::CreatePEMatrix() {
  _pe_matrix.clear();
  
  // Use final filtered data if available, otherwise use selected data, fallback to original data
  const auto* pe_data = &_flash_ophit_pe;
  const auto* ch_data = &_flash_ophit_ch;
  
  if (!_flash_ophit_pe_final.empty() && !_flash_ophit_ch_final.empty()) {
    pe_data = &_flash_ophit_pe_final;
    ch_data = &_flash_ophit_ch_final;
    if(fVerbosity > 0) std::cout << "CreatePEMatrix: Using final filtered data" << std::endl;
  } else if (!_flash_ophit_pe_sel.empty() && !_flash_ophit_ch_sel.empty()) {
    pe_data = &_flash_ophit_pe_sel;
    ch_data = &_flash_ophit_ch_sel;
    if(fVerbosity > 0) std::cout << "CreatePEMatrix: Using selected data" << std::endl;
  } else {
    if(fVerbosity > 0) std::cout << "CreatePEMatrix: Using original data" << std::endl;
  }
  
  if (pe_data->empty() || ch_data->empty()) {
    if(fVerbosity > 0) {
      std::cout << "CreatePEMatrix: No flash data available at all" << std::endl;
    }
    // Initialize empty matrix for consistency
    _pe_matrix.resize(1, std::vector<float>(312, 0.0f));
    return;
  }
  
  // Initialize PE matrix with zeros [1 event x 312 channels]
  _pe_matrix.resize(1, std::vector<float>(312, 0.0f));
  
  // Sum PE values from all flashes for each channel into single event
  for (size_t flash = 0; flash < pe_data->size(); ++flash) {
    // Bounds checking for channel data consistency
    if (flash >= ch_data->size()) {
      if(fVerbosity > 0) {
        std::cout << "Warning: PE data has more flashes than channel data at flash " << flash << std::endl;
      }
      break;
    }
    
    size_t min_size = std::min((*pe_data)[flash].size(), (*ch_data)[flash].size());
    for (size_t j = 0; j < min_size; ++j) {
      float pe = (*pe_data)[flash][j];
      int channel = (*ch_data)[flash][j];
      
      // Ensure channel is within valid range [0, 311]
      if (channel >= 0 && channel < 312) {
        _pe_matrix[0][channel] += pe;  // Always add to event 0 (single event)
      } else if (fVerbosity > 0) {
        std::cout << "Warning: Invalid channel " << channel << " at flash " << flash << ", hit " << j << std::endl;
      }
    }
  }
  
  if(fVerbosity > 0) {
    std::cout << "CreatePEMatrix: Generated matrix with " << _pe_matrix.size() 
              << " events and 312 channels" << std::endl;
  }
}


// -------- Load PMT Maps --------
void opdet::SBNDPDSProducer::LoadPMTMaps() {
  // Load coated PMT map
  std::ifstream coated_file(fCoatedPMTMapPath);
  std::ifstream uncoated_file(fUncoatedPMTMapPath);
  
  if (!coated_file.is_open()) {
    std::string msg = "Could not open coated PMT map file: " + fCoatedPMTMapPath;
    std::cout << "Warning: " << msg << std::endl;
    if(fVerbosity > 0) throw std::runtime_error(msg);
    return;
  }
  
  if (!uncoated_file.is_open()) {
    std::string msg = "Could not open uncoated PMT map file: " + fUncoatedPMTMapPath;
    std::cout << "Warning: " << msg << std::endl;
    if(fVerbosity > 0) throw std::runtime_error(msg);
    return;
  }
  
  std::string line;
  _coated_pmt_map.clear();
  _uncoated_pmt_map.clear();
  
  // Load coated map
  while (std::getline(coated_file, line)) {
    std::vector<int> row;
    std::stringstream ss(line);
    std::string cell;
    
    while (std::getline(ss, cell, ',')) {
      row.push_back(std::stoi(cell));
    }
    _coated_pmt_map.push_back(row);
  }
  coated_file.close();
  
  // Load uncoated map
  while (std::getline(uncoated_file, line)) {
    std::vector<int> row;
    std::stringstream ss(line);
    std::string cell;
    
    while (std::getline(ss, cell, ',')) {
      row.push_back(std::stoi(cell));
    }
    _uncoated_pmt_map.push_back(row);
  }
  uncoated_file.close();
  
  if(fVerbosity > 0) {
    std::cout << "Loaded PMT maps: coated(" << _coated_pmt_map.size() << "x" << (_coated_pmt_map.empty() ? 0 : _coated_pmt_map[0].size()) 
              << "), uncoated(" << _uncoated_pmt_map.size() << "x" << (_uncoated_pmt_map.empty() ? 0 : _uncoated_pmt_map[0].size()) << ")" << std::endl;
  }
}


// -------- Select Non-Empty Half --------
std::vector<std::vector<float>> opdet::SBNDPDSProducer::SelectNonEmptyHalf(
    const std::vector<std::vector<float>>& left_half, 
    const std::vector<std::vector<float>>& right_half,
    const std::string& method) {
  
  float left_score = 0.0f, right_score = 0.0f;
  
  if (method == "max") {
    for (const auto& row : left_half) {
      for (float val : row) {
        left_score = std::max(left_score, val);
      }
    }
    for (const auto& row : right_half) {
      for (float val : row) {
        right_score = std::max(right_score, val);
      }
    }
  }
  else if (method == "sum") {
    for (const auto& row : left_half) {
      for (float val : row) {
        left_score += val;
      }
    }
    for (const auto& row : right_half) {
      for (float val : row) {
        right_score += val;
      }
    }
  }
  else if (method == "nonzero") {
    for (const auto& row : left_half) {
      for (float val : row) {
        if (val != 0.0f) left_score += 1.0f;
      }
    }
    for (const auto& row : right_half) {
      for (float val : row) {
        if (val != 0.0f) right_score += 1.0f;
      }
    }
  }
  else if (method == "mean_top") {
    std::vector<float> left_vals, right_vals;
    for (const auto& row : left_half) {
      for (float val : row) {
        left_vals.push_back(val);
      }
    }
    for (const auto& row : right_half) {
      for (float val : row) {
        right_vals.push_back(val);
      }
    }
    
    std::sort(left_vals.rbegin(), left_vals.rend());
    std::sort(right_vals.rbegin(), right_vals.rend());
    
    int top_n = std::min(5, (int)left_vals.size());
    for (int i = 0; i < top_n; ++i) {
      left_score += left_vals[i];
      right_score += right_vals[i];
    }
    left_score /= top_n;
    right_score /= top_n;
  }
  
  return left_score >= right_score ? left_half : right_half;
}


// -------- Create PE Images --------
void opdet::SBNDPDSProducer::CreatePEImages() {
  _pe_images.clear();
  
  if (_pe_matrix.empty() || _coated_pmt_map.empty() || _uncoated_pmt_map.empty()) {
    if(fVerbosity > 0) {
      std::cout << "CreatePEImages: Missing PE matrix or PMT maps" << std::endl;
    }
    return;
  }
  
  int n_flashes = _pe_matrix.size();  // This is actually number of flashes, not events
  int n_events = 1;  // For single event processing, we always have 1 event
  int ch_y = _coated_pmt_map.size();
  int ch_z = _coated_pmt_map[0].size();
  int map_count = 2; // coated and uncoated
  
  // Initialize pe_matrices_map with shape (map_count, n_events, ch_y, ch_z)
  std::vector<std::vector<std::vector<std::vector<float>>>> pe_matrices_map(
    map_count, std::vector<std::vector<std::vector<float>>>(
      n_events, std::vector<std::vector<float>>(
        ch_y, std::vector<float>(ch_z, 0.0f))));
  
  // Map pe_matrix values to pe_matrices_map for each event
  std::vector<std::vector<std::vector<int>>*> maps = {&_uncoated_pmt_map, &_coated_pmt_map};
  
  for (int idx = 0; idx < map_count; ++idx) {
    auto& map_ = *maps[idx];
    for (int y = 0; y < ch_y; ++y) {
      for (int z = 0; z < ch_z; ++z) {
        int channel = map_[y][z];
        if (channel >= 0 && channel < 312) {
          // Sum PE values from all flashes for this channel for the single event
          float total_pe = 0.0f;
          for (int flash = 0; flash < n_flashes; ++flash) {
            total_pe += _pe_matrix[flash][channel];
          }
          pe_matrices_map[idx][0][y][z] = total_pe;  // Single event at index 0
        }
      }
    }
  }
  
  // Normalize each map group
  float max_val_1 = 0.0f;
  
  // Find max values for normalization
  for (int i = 0; i < n_events; ++i) {
    for (int y = 0; y < ch_y; ++y) {
      for (int z = 0; z < ch_z; ++z) {
        max_val_1 = std::max(max_val_1, std::max(pe_matrices_map[0][i][y][z], pe_matrices_map[1][i][y][z]));
      }
    }
  }
  max_val_1 = (max_val_1 > 0) ? max_val_1 : 1.0f;
  
  // Apply normalization
  for (int idx = 0; idx < map_count; ++idx) {
    for (int i = 0; i < n_events; ++i) {
      for (int y = 0; y < ch_y; ++y) {
        for (int z = 0; z < ch_z; ++z) {
          pe_matrices_map[idx][i][y][z] /= max_val_1;
        }
      }
    }
  }
  
  // Create image with half selection: shape (n_events, ch_y/2, ch_z, map_count)
  _pe_images.resize(n_events, std::vector<std::vector<std::vector<float>>>(
    ch_y/2, std::vector<std::vector<float>>(
      ch_z, std::vector<float>(map_count, 0.0f))));
  
  for (int i = 0; i < n_events; ++i) {
    for (int idx = 0; idx < map_count; ++idx) {
      // Get the matrix for this event and map: shape (ch_y, ch_z)
      const auto& event_matrix = pe_matrices_map[idx][i];
      
      // Split into top and bottom halves
      std::vector<std::vector<float>> top_half(event_matrix.begin(), event_matrix.begin() + ch_y/2);
      std::vector<std::vector<float>> bottom_half(event_matrix.begin() + ch_y/2, event_matrix.end());
      
      // Select the non-empty half
      auto selected_half = SelectNonEmptyHalf(top_half, bottom_half, "max");
      
      // Store the selected half in the image
      for (int y = 0; y < ch_y/2; ++y) {
        for (int z = 0; z < ch_z; ++z) {
          _pe_images[i][y][z][idx] = selected_half[y][z];
        }
      }
    }
  }
  
  if(fVerbosity > 0) {
    std::cout << "CreatePEImages: Generated images with shape (" << n_events 
              << ", " << ch_y/2 << ", " << ch_z << ", " << map_count << ")" << std::endl;
  }
}


// -------- Function to fill the MCTruth information --------
void opdet::SBNDPDSProducer::FillMCTruth(art::Event const& e){
  
  
  if(fMCTruthModuleLabel.size()!=fMCTruthInstanceLabel.size()){
    std::string msg = "MCTruthModuleLabel (" + std::to_string(fMCTruthModuleLabel.size()) + 
                     ") and MCTruthInstanceLabel (" + std::to_string(fMCTruthInstanceLabel.size()) + 
                     ") vectors must have the same size";
    std::cout << msg << std::endl;
    throw std::invalid_argument(msg);
  }

  art::Handle< std::vector<simb::MCTruth> > MCTruthListHandle;

  for (size_t s = 0; s < fMCTruthModuleLabel.size(); s++) {
    
    e.getByLabel(fMCTruthModuleLabel[s], fMCTruthInstanceLabel[s], MCTruthListHandle);

    if( !MCTruthListHandle.isValid() || MCTruthListHandle->empty() ) {   
      std::string msg = "MCTruth with label " + fMCTruthModuleLabel[s] + " and instance " + fMCTruthInstanceLabel[s] + 
                       " not found or empty in event " + std::to_string(_eventID);
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
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
void opdet::SBNDPDSProducer::FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
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


// -------- Function to initialize the channel dictionary --------
void opdet::SBNDPDSProducer::InitializeChannelDict(){
  
  // Create PDS mapping algorithm instance
  opdet::sbndPDMapAlg pdsMap;
  
  // Clear the channel dictionary
  fChannelDict.clear();
  
  // Create channel dictionary mapping OpDetID (channel) to OpDetType
  for(size_t ch = 0; ch < pdsMap.size(); ++ch) {
    std::string pdType = pdsMap.pdType(ch);
    int typeInt = -1;
    
    // Convert string type to integer following SBND convention
    if (pdType == "pmt_coated") {
      typeInt = static_cast<int>(opdet::kPMTCoated);
    } else if (pdType == "pmt_uncoated") {
      typeInt = static_cast<int>(opdet::kPMTUncoated);
    } else if (pdType == "xarapuca_vuv") {
      typeInt = static_cast<int>(opdet::kXARAPUCAVUV);
    } else if (pdType == "xarapuca_vis") {
      typeInt = static_cast<int>(opdet::kXARAPUCAVIS);
    } else {
      typeInt = static_cast<int>(opdet::kPDUnknown);
    }
    
    fChannelDict[static_cast<int>(ch)] = typeInt;
    
    // Debug individual channels only at highest verbosity
    if(fVerbosity > 2) {
      std::cout << "Channel " << ch << " -> Type: " << pdType << " (" << typeInt << ")" << std::endl;
    }
  }
  
  if(fVerbosity > 0) {
    std::cout << "Initialized channel dictionary with " << fChannelDict.size() << " channels" << std::endl;
  }
}


// -------- Function to classify channels by type and parity --------
void opdet::SBNDPDSProducer::ClassifyChannels(){
  
  fPMTEven.clear(); fPMTOdd.clear();
  fXASEven.clear(); fXASOdd.clear();
  
  for(const auto& pair : fChannelDict) {
    int ch = pair.first;
    int type = pair.second;
    
    bool isPMT = (type == static_cast<int>(opdet::kPMTCoated) || 
                  type == static_cast<int>(opdet::kPMTUncoated));
    bool isXAS = (type == static_cast<int>(opdet::kXARAPUCAVUV) || 
                  type == static_cast<int>(opdet::kXARAPUCAVIS));
    bool isEven = (ch % 2 == 0);
    
    if(isPMT) {
      if(isEven) fPMTEven.insert(ch);
      else fPMTOdd.insert(ch);
    } else if(isXAS) {
      if(isEven) fXASEven.insert(ch);
      else fXASOdd.insert(ch);
    }
  }
  
  if(fVerbosity > 0) {
    std::cout << "Channel classification - PMT even: " << fPMTEven.size() 
              << ", PMT odd: " << fPMTOdd.size()
              << ", XAS even: " << fXASEven.size() 
              << ", XAS odd: " << fXASOdd.size() << std::endl;
  }
}


// -------- Function to categorize first channel of a flash --------
int opdet::SBNDPDSProducer::CategorizeFirstChannel(const std::vector<int>& channels){
  
  if(channels.empty()) return -1;
  
  int ch = channels[0];
  if(fPMTEven.count(ch)) return 0;
  if(fPMTOdd.count(ch))  return 1;
  if(fXASEven.count(ch)) return 2;
  if(fXASOdd.count(ch))  return 3;
  
  return -1; // Unclassified
}


// -------- Template function to filter arrays by mask --------
template<typename T>
std::vector<std::vector<T>> opdet::SBNDPDSProducer::FilterByMask(
    const std::vector<std::vector<T>>& array, 
    const std::vector<std::vector<bool>>& mask) {
  
  std::vector<std::vector<T>> result;
  
  for(size_t i = 0; i < array.size() && i < mask.size(); ++i) {
    std::vector<T> filteredEvent;
    for(size_t j = 0; j < array[i].size() && j < mask[i].size(); ++j) {
      if(mask[i][j]) {
        filteredEvent.push_back(array[i][j]);
      }
    }
    result.push_back(filteredEvent);
  }
  
  return result;
}


// -------- Main flash selection function --------
void opdet::SBNDPDSProducer::ApplyFlashSelection(){
  
  // Clear output vectors
  _flash_ophit_pe_sel.clear();
  _flash_ophit_ch_sel.clear(); 
  _flash_ophit_time_sel.clear();
  _categorized_flashes.clear();
  _mc_dEpromx_sel.clear();
  _mc_dEpromy_sel.clear();
  _mc_dEpromz_sel.clear();
  _mc_dEtpc_sel.clear();
  
  if(_flash_ophit_ch.empty()) {
    if(fVerbosity > 0) std::cout << "No flashes to process" << std::endl;
    return;
  }
  
  // 1. Categorize flashes based on first channel
  std::vector<int> categorized;
  for(const auto& flash_ch : _flash_ophit_ch) {
    categorized.push_back(CategorizeFirstChannel(flash_ch));
  }
  
  // 2. Calculate sum PE per flash
  std::vector<float> sum_pe;
  for(const auto& flash_pe : _flash_ophit_pe) {
    float total = 0;
    for(float pe : flash_pe) total += pe;
    sum_pe.push_back(total);
  }
  
  // 3. Create masks and calculate sums by group
  float sum_even = 0, sum_odd = 0;
  std::vector<bool> mask_even, mask_odd;
  
  for(size_t i = 0; i < categorized.size(); ++i) {
    bool is_even = (categorized[i] == 0 || categorized[i] == 2);
    bool is_odd  = (categorized[i] == 1 || categorized[i] == 3);
    
    mask_even.push_back(is_even);
    mask_odd.push_back(is_odd);
    
    if(is_even) sum_even += sum_pe[i];
    if(is_odd)  sum_odd += sum_pe[i];
  }
  
  // 4. Decision logic
  bool decision = (sum_even >= sum_odd);
  size_t n_flashes = categorized.size();
  
  std::vector<bool> selected_mask;
  if(n_flashes <= 2) {
    // Keep all flashes for <= 2 flashes
    selected_mask.assign(n_flashes, true);
  } else {
    // Select group based on decision
    selected_mask = decision ? mask_even : mask_odd;
  }
  
  // 5. Apply mask to flash data
  std::vector<std::vector<std::vector<float>>> pe_2d = {_flash_ophit_pe};
  std::vector<std::vector<std::vector<int>>> ch_2d = {_flash_ophit_ch};
  std::vector<std::vector<std::vector<float>>> time_2d = {_flash_ophit_time};
  std::vector<std::vector<bool>> mask_2d = {selected_mask};
  
  auto pe_filtered = FilterByMask(pe_2d, mask_2d);
  auto ch_filtered = FilterByMask(ch_2d, mask_2d);  
  auto time_filtered = FilterByMask(time_2d, mask_2d);
  
  if(!pe_filtered.empty()) _flash_ophit_pe_sel = pe_filtered[0];
  if(!ch_filtered.empty()) _flash_ophit_ch_sel = ch_filtered[0];
  if(!time_filtered.empty()) _flash_ophit_time_sel = time_filtered[0];
  
  // 6. TPC selection - select winning TPC based on decision
  std::vector<std::vector<double>> selector = {{decision ? 1.0 : 0.0, decision ? 0.0 : 1.0}};
  std::vector<std::vector<double>> dEpromx_2d = {_mc_dEpromx};
  std::vector<std::vector<double>> dEpromy_2d = {_mc_dEpromy};
  std::vector<std::vector<double>> dEpromz_2d = {_mc_dEpromz};
  std::vector<std::vector<double>> dEtpc_2d = {_mc_dEtpc};
  
  // Select TPC values based on decision with bounds checking
  size_t max_tpc = std::min({(size_t)2, _mc_dEpromx.size(), _mc_dEpromy.size(), _mc_dEpromz.size(), _mc_dEtpc.size()});
  for(size_t i = 0; i < max_tpc; ++i) {
    if(selector[0][i] > 0.5) {
      _mc_dEpromx_sel.push_back(_mc_dEpromx[i]);
      _mc_dEpromy_sel.push_back(_mc_dEpromy[i]);
      _mc_dEpromz_sel.push_back(_mc_dEpromz[i]);
      _mc_dEtpc_sel.push_back(_mc_dEtpc[i]);
      break; // Only select one TPC
    }
  }
  
  if(fVerbosity > 0) {
    std::cout << "Flash selection: " << n_flashes << " flashes, decision=" << decision 
              << " (even=" << sum_even << ", odd=" << sum_odd << ")" << std::endl;
    std::cout << "Selected " << _flash_ophit_pe_sel.size() << " flashes" << std::endl;
  }
}


// -------- Final energy deposition filter --------
void opdet::SBNDPDSProducer::ApplyFinalEnergyFilter(){
  
  // Clear final output vectors
  _flash_ophit_pe_final.clear();
  _flash_ophit_ch_final.clear();
  _flash_ophit_time_final.clear();
  _nuvT_final.clear();
  _nuvZ_final.clear();
  _mc_dEpromx_final.clear();
  _mc_dEpromy_final.clear();
  _mc_dEpromz_final.clear();
  _mc_dEtpc_final.clear();
  
  // Check if we have valid selected data
  if(_mc_dEpromx_sel.empty() || _mc_dEpromy_sel.empty() || 
     _mc_dEpromz_sel.empty() || _mc_dEtpc_sel.empty()) {
    if(fVerbosity > 0) {
      std::cout << "Final energy filter: No selected energy data available" << std::endl;
    }
    return;
  }
  
  // Apply energy deposition mask: 
  // (dEpromx != -999) & (dEpromy != -999) & (dEpromz != -999) & (dEtpc > 50)
  bool passEnergyFilter = false;
  
  // Check all selected TPC data (typically just one after flash selection)
  for(size_t i = 0; i < _mc_dEpromx_sel.size(); ++i) {
    bool validX = (_mc_dEpromx_sel[i] != fDefaultSimIDE);
    bool validY = (_mc_dEpromy_sel[i] != fDefaultSimIDE);  
    bool validZ = (_mc_dEpromz_sel[i] != fDefaultSimIDE);
    bool energyCut = (_mc_dEtpc_sel[i] > 50.0);
    
    if(validX && validY && validZ && energyCut) {
      passEnergyFilter = true;
      // Store the passing values
      _mc_dEpromx_final.push_back(_mc_dEpromx_sel[i]);
      _mc_dEpromy_final.push_back(_mc_dEpromy_sel[i]);
      _mc_dEpromz_final.push_back(_mc_dEpromz_sel[i]);
      _mc_dEtpc_final.push_back(_mc_dEtpc_sel[i]);
      
      if(fVerbosity > 0) {
        std::cout << "Energy filter passed: dE=" << _mc_dEtpc_sel[i] 
                  << " MeV, pos=(" << _mc_dEpromx_sel[i] << "," 
                  << _mc_dEpromy_sel[i] << "," << _mc_dEpromz_sel[i] << ")" << std::endl;
      }
      break; // Only need one passing entry
    }
  }
  
  if(passEnergyFilter) {
    // Copy flash data if energy filter passes
    _flash_ophit_pe_final = _flash_ophit_pe_sel;
    _flash_ophit_ch_final = _flash_ophit_ch_sel;
    _flash_ophit_time_final = _flash_ophit_time_sel;
    
    // Copy neutrino data
    _nuvT_final = _nuvT;
    _nuvZ_final = _nuvZ;
    
    if(fVerbosity > 0) {
      std::cout << "Final energy filter: Event PASSED - " 
                << _flash_ophit_pe_final.size() << " flashes retained" << std::endl;
    }
  } else {
    if(fVerbosity > 0) {
      std::cout << "Final energy filter: Event FAILED - no valid energy deposition (>50 MeV)" << std::endl;
    }
  }
}