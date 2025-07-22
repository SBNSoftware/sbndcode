#include "sbndcode/PosRecoCVN/module/SBNDPDSProducer_module.hh"

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
    dE_neutrinowindow( 0.0 )
{
    produces<PixelMapVars>();
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
      std::cout << "OpHit with label " << fOpHitsModuleLabel[s] << " not found..." << std::endl;
      throw std::exception();
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

  // Crear y llenar el producto PixelMapVars
  auto pixelVars = std::make_unique<PixelMapVars>();
  // Copiar solo los datos necesarios para crear mapas del detector
  pixelVars->flash_ophit_pe = _flash_ophit_pe;
  pixelVars->flash_ophit_ch = _flash_ophit_ch;
  pixelVars->flash_ophit_time = _flash_ophit_time;
  pixelVars->nuvT = _nuvT;
  pixelVars->dEpromx = _mc_dEpromx;
  pixelVars->dEpromy = _mc_dEpromy;
  pixelVars->dEpromz = _mc_dEpromz;
  pixelVars->dEtpc = _mc_dEtpc;
  pixelVars->nuvZ = _nuvZ;

  // Poner el producto en el evento
  e.put(std::move(pixelVars));

}


// -------- Function to fill the MCTruth information --------
void opdet::SBNDPDSProducer::FillMCTruth(art::Event const& e){
  
  
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