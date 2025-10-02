#include "PosRecoCVNProducer_module.hh"

#include "larcorealg/Geometry/OpDetGeo.h"
#include <chrono>
#include <unordered_map>
#include <iomanip>
#include <cstdlib>

// ============================================================================
// Constructor - Initialize CNN position reconstruction module
//   - Read FCL configuration parameters
//   - Initialize channel maps (PMT type classification)
//   - Load TensorFlow model from multiple search paths
// ============================================================================
opdet::PosRecoCVNProducer::PosRecoCVNProducer(fhicl::ParameterSet const& p)
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
    fSavePixelMapVars( p.get<bool>("SavePixelMapVars", true) ),
    fVerbosity( p.get<int>("Verbosity") ),
    fCoatedPMTMapPath("coatedPMT_map.csv"),
    fUncoatedPMTMapPath("uncoatedPMT_map.csv"),
    fModelPath("saved_model"),
    fProcessingMode( p.get<std::string>("ProcessingMode", "MC_testing") ),
    fInputNames( p.get<std::vector<std::string>>("InputNames", {}) ),
    fOutputNames( p.get<std::vector<std::string>>("OutputNames", {}) ),
    fCustomNormFactor( p.get<double>("CustomNormFactor", -1.0) ),
    fPredictionTolerance( p.get<double>("PredictionTolerance", 0.05) ),
    fSkipNeutrinoFilter( p.get<bool>("SkipNeutrinoFilter", false) ),
    fSbndcodeVersion( p.get<std::string>("SbndcodeVersion", "v10_09_00") )
{
    if(fSavePixelMapVars) produces<PixelMapVars>();

    InitializeChannelDict();
    ClassifyChannels();
    LoadPMTMaps();

    // Load TensorFlow model - search multiple paths
    if (!fModelPath.empty()) {
        std::string model_path = fModelPath;
        std::vector<std::string> search_paths = {
            model_path,
            "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/inference/tf/v0901_trained_w_165k_resnet18/" + model_path,
            std::string(getenv("MRB_INSTALL") ? getenv("MRB_INSTALL") : "") + "/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/inference/tf/v0901_trained_w_165k_resnet18/" + model_path,
            std::string(getenv("MRB_SOURCE") ? getenv("MRB_SOURCE") : "") + "/sbndcode/sbndcode/PosRecoCVN/inference/tf/v0901_trained_w_165k_resnet18/" + model_path,
            "../" + model_path,
            "./" + model_path
        };

        bool model_found = false;
        for (const auto& path : search_paths) {
            std::ifstream test_file(path + "/saved_model.pb");
            if (test_file.is_open()) {
                model_path = path;
                model_found = true;
                test_file.close();
                break;
            }
        }

        if (!model_found) {
            std::cerr << "ERROR: TensorFlow model not found: " << fModelPath << std::endl;
        } else {
            int nOutputs = fOutputNames.empty() ? -1 : fOutputNames.size();
            fTFGraph = tf::Graph::create(model_path.c_str(), fInputNames, fOutputNames, true, 1, nOutputs > 0 ? nOutputs : 3);
            if (!fTFGraph) {
                std::cerr << "ERROR: Failed to load model from " << model_path << std::endl;
            } else if(fVerbosity > 0) {
                std::cout << "TensorFlow model loaded: " << model_path << std::endl;
            }
        }
    }
}

// ============================================================================
// beginJob - Initialize output TTree with appropriate branches
//   - MC_testing: predictions + ground truth + errors
//   - MC_inference/DATA_inference: predictions only
// ============================================================================
void opdet::PosRecoCVNProducer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fInferenceTree = tfs->make<TTree>("inference_tree", "CNN Position Reconstruction");

  // Always present: event ID and predictions
  fInferenceTree->Branch("run", &fTreeData.run, "run/I");
  fInferenceTree->Branch("subrun", &fTreeData.subrun, "subrun/I");
  fInferenceTree->Branch("event", &fTreeData.event, "event/I");
  fInferenceTree->Branch("passed_filters", &fTreeData.passedFilters, "passed_filters/O");
  fInferenceTree->Branch("pred_x", &fTreeData.predX, "pred_x/D");
  fInferenceTree->Branch("pred_y", &fTreeData.predY, "pred_y/D");
  fInferenceTree->Branch("pred_z", &fTreeData.predZ, "pred_z/D");

  // MC_testing mode: add ground truth and performance metrics
  if(fProcessingMode == "MC_testing") {
    fInferenceTree->Branch("true_x", &fTreeData.trueX, "true_x/D");
    fInferenceTree->Branch("true_y", &fTreeData.trueY, "true_y/D");
    fInferenceTree->Branch("true_z", &fTreeData.trueZ, "true_z/D");
    fInferenceTree->Branch("diff_x", &fTreeData.diffX, "diff_x/D");
    fInferenceTree->Branch("diff_y", &fTreeData.diffY, "diff_y/D");
    fInferenceTree->Branch("diff_z", &fTreeData.diffZ, "diff_z/D");
    fInferenceTree->Branch("error_3d", &fTreeData.error3D, "error_3d/D");
    fInferenceTree->Branch("nuv_t", &fTreeData.nuvT, "nuv_t/D");
    fInferenceTree->Branch("nuv_z", &fTreeData.nuvZ, "nuv_z/D");
    fInferenceTree->Branch("deposited_energy", &fTreeData.dEtpc, "deposited_energy/D");
  }

  if(fVerbosity > 0) {
    std::cout << "TTree initialized (" << fProcessingMode << "): "
              << fInferenceTree->GetNbranches() << " branches" << std::endl;
  }
}

// ============================================================================
// produce - Main event processing (per-event analysis)
//   1. Load MC truth and particles (mode-dependent)
//   2. Load and filter OpFlashes/OpHits (beam window: 0.367-1.9 μs)
//   3. Apply event filters (neutrino count, optical data availability)
//   4. Process flashes: channel classification → energy filter
//   5. Generate PE images: matrix [312 ch] → images [ch_y/2, ch_z, 2]
//   6. Run CNN inference: predict (x,y,z) position
//   7. Fill output TTree and optionally save PixelMapVars
// ============================================================================
void opdet::PosRecoCVNProducer::produce(art::Event& e)
{
  auto event_start = std::chrono::high_resolution_clock::now();
  auto section_start = event_start;

  ClearEventData();
  LogTiming("ClearEventData", section_start);

  // Event identification
  fEventInfo.eventID = e.id().event();
  fEventInfo.runID = e.id().run();
  fEventInfo.subrunID = e.id().subRun();

  if(fVerbosity > 1) {
    std::cout << "Run " << fEventInfo.runID << " | Subrun " << fEventInfo.subrunID
              << " | Event " << fEventInfo.eventID << std::endl;
  }

  // Art services
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<detinfo::DetectorClocksService> timeservice;

  // --- Load MCTruth: Full data for MC_testing, minimal (nuvT only) for MC_inference, none for DATA
  bool calculateGroundTruth = (fProcessingMode == "MC_testing");

  if (calculateGroundTruth) {
    if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();
    FillMCTruth(e);
    LogTiming("FillMCTruth", section_start);
  } else {
    if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();

    // MC_inference mode: Extract only neutrino time for event filtering
    art::Handle<std::vector<simb::MCTruth>> mcListHandle;
    for (size_t s = 0; s < fMCTruthModuleLabel.size(); s++) {
      e.getByLabel(fMCTruthModuleLabel[s], fMCTruthInstanceLabel[s], mcListHandle);
      if(!mcListHandle.isValid()) continue;

      for (unsigned int i = 0; i < mcListHandle->size(); ++i) {
        art::Ptr<simb::MCTruth> mctruth(mcListHandle, i);
        if (mctruth->Origin() == simb::kBeamNeutrino) {
          simb::MCParticle nu = mctruth->GetNeutrino().Nu();
          fMCData.nuvT.push_back(nu.T());
        }
      }
    }

    if(fVerbosity > 1) std::cout << "Extracted " << fMCData.nuvT.size() << " neutrinos" << std::endl;
    LogTiming("MinimalMCTruth", section_start);
  }

  // --- Load MCParticles (only for MC_testing mode)
  if (calculateGroundTruth) {
    if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();

    fMCData.mcpartVec.clear();
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    e.getByLabel(fMCModuleLabel, mcParticleHandle);
    if(mcParticleHandle.isValid()){
      fMCData.mcpartVec = *mcParticleHandle;
      if(fVerbosity>2) std::cout << "Loaded " << fMCData.mcpartVec.size() << " MCParticles" << std::endl;
    } else {
      if(fVerbosity>0) std::cout << "WARNING: MCParticles not found (" << fMCModuleLabel << ")" << std::endl;
    }

    LogTiming("LoadMCParticles", section_start);
  }

  // --- Process MCParticles: Extract trajectories and energy depositions for ground truth
  if (calculateGroundTruth) {
    fMCData.stepX.clear(); fMCData.stepY.clear(); fMCData.stepZ.clear(); fMCData.stepT.clear();
    fMCData.dE.clear(); fMCData.E.clear();
    fMCData.trackID.clear(); fMCData.motherID.clear(); fMCData.PDGcode.clear(); fMCData.process.clear();
    fMCData.StartPx.clear(); fMCData.StartPy.clear(); fMCData.StartPz.clear();
    fMCData.EndPx.clear(); fMCData.EndPy.clear(); fMCData.EndPz.clear();
    fMCData.energydep.clear(); fMCData.energydepX.clear(); fMCData.energydepY.clear(); fMCData.energydepZ.clear();
    fMCData.InTimeCosmics=0; fMCData.InTimeCosmicsTime.clear();
    fMCData.neutrinowindow = 0.0;

  if(fVerbosity>2) std::cout << "Processing MCParticles" << std::endl;

  // Loop through all MCParticles
  for(const auto& pPart : fMCData.mcpartVec){

    // Filter by PDG code if list specified in FCL
    if(!fKeepPDGCode.empty() && std::find(fKeepPDGCode.begin(),fKeepPDGCode.end(),pPart.PdgCode()) == fKeepPDGCode.end())
      continue;

    const simb::MCTrajectory truetrack = pPart.Trajectory();
    std::vector<double> xpoints, ypoints, zpoints, tpoints;

    // Always save particle vertex (first trajectory point)
    xpoints.push_back(pPart.Position(0).X());
    ypoints.push_back(pPart.Position(0).Y());
    zpoints.push_back(pPart.Position(0).Z());
    tpoints.push_back(pPart.Position(0).T());

    // Loop trajectory points, filtering by G4BufferBox and G4BeamTimeWindow (FCL parameters)
    for(size_t i_s=1; i_s < pPart.NumberTrajectoryPoints(); i_s++){
      double t = pPart.Position(i_s).T();
      double x = pPart.Position(i_s).X();
      double y = pPart.Position(i_s).Y();
      double z = pPart.Position(i_s).Z();

      // Apply spatial and temporal cuts
      if( x < fG4BufferBoxX.at(0) || x > fG4BufferBoxX.at(1)
            || y < fG4BufferBoxY.at(0) || y > fG4BufferBoxY.at(1)
                || z < fG4BufferBoxZ.at(0) || z > fG4BufferBoxZ.at(1)) continue;
      if( t < fG4BeamWindow.at(0) || t > fG4BeamWindow.at(1) ) continue;

      xpoints.push_back(x);
      ypoints.push_back(y);
      zpoints.push_back(z);
      tpoints.push_back(t);
    }

    // Store particle trajectory and kinematics
    fMCData.stepX.push_back(xpoints);
    fMCData.stepY.push_back(ypoints);
    fMCData.stepZ.push_back(zpoints);
    fMCData.stepT.push_back(tpoints);
    fMCData.E.push_back(pPart.E());
    fMCData.process.push_back(pPart.Process());
    fMCData.trackID.push_back(pPart.TrackId());
    fMCData.motherID.push_back(pPart.Mother());
    fMCData.PDGcode.push_back(pPart.PdgCode());
    fMCData.StartPx.push_back(pPart.Px(0));
    fMCData.StartPy.push_back(pPart.Py(0));
    fMCData.StartPz.push_back(pPart.Pz(0));
    fMCData.EndPx.push_back(pPart.EndPx());
    fMCData.EndPy.push_back(pPart.EndPy());
    fMCData.EndPz.push_back(pPart.EndPz());

    // Extract energy depositions from SimIDEs
    double endep=0;
    std::vector<double> truedE, truedE_vecX, truedE_vecY, truedE_vecZ;
    std::vector<const sim::IDE*> ides_v = bt_serv->TrackIdToSimIDEs_Ps(pPart.TrackId());
    for(auto *ide:ides_v){
      endep+=ide->energy/3.;  // Divide by 3: avoid triple-counting over 3 wire planes
      if( pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1) ){
        truedE.push_back(ide->energy/3.);
        truedE_vecX.push_back(ide->x);
        truedE_vecY.push_back(ide->y);
        truedE_vecZ.push_back(ide->z);
        fMCData.neutrinowindow+=ide->energy/3.;
      }
    }

    fMCData.dE.push_back(endep);
    fMCData.energydep.push_back(truedE);
    fMCData.energydepX.push_back(truedE_vecX);
    fMCData.energydepY.push_back(truedE_vecY);
    fMCData.energydepZ.push_back(truedE_vecZ);

    // Count in-time cosmics with significant energy deposition
    art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(pPart.TrackId());
    if( pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1) && truth->Origin()==2 && endep>1.){
      fMCData.InTimeCosmics++;
      fMCData.InTimeCosmicsTime.push_back(pPart.Position(0).T());
    }

    if(fVerbosity>2){
      std::cout<<"-.-.-.E="<<std::setw(9)<<pPart.E()<<" PDG="<<std::setw(7)<<pPart.PdgCode()
               <<" ID="<<std::setw(5)<<pPart.TrackId()<<" Mother="<<std::setw(5)<<pPart.Mother()<<" dE="<<std::setw(7)<<endep
               <<" T="<<std::setw(8)<<pPart.T()<<" V=("<<pPart.Vx()<<","<<pPart.Vy()<<","<<pPart.Vz()
               <<") NP:"<<pPart.NumberTrajectoryPoints()<<" "<<pPart.Process()<<std::endl;
    }
  }

  // Calculate energy-weighted centroids (dEprom*) and other statistics
  FillAverageDepositedEnergyVariables(fMCData.energydep,fMCData.energydepX,fMCData.energydepY,fMCData.energydepZ,fMCData.stepT,fMCData.dEtpc,fMCData.dEpromx,fMCData.dEpromy,fMCData.dEpromz,fMCData.dEspreadx,fMCData.dEspready,fMCData.dEspreadz,fMCData.dElowedges,fMCData.dEmaxedges);

  if(fVerbosity>1){
    std::cout<<"InTimeCosmic: "<<fMCData.InTimeCosmics<<" | BeamWindow dE: "<<fMCData.neutrinowindow<<" MeV\n"
             <<"TPC0  dE="<<fMCData.dEtpc[0]<<"  <x,y,z>=("<<fMCData.dEpromx[0]<<","<<fMCData.dEpromy[0]<<","<<fMCData.dEpromz[0]
             <<")  Sp=("<<fMCData.dEspreadx[0]<<","<<fMCData.dEspready[0]<<","<<fMCData.dEspreadz[0]<<")\n"
             <<"TPC1  dE="<<fMCData.dEtpc[1]<<"  <x,y,z>=("<<fMCData.dEpromx[1]<<","<<fMCData.dEpromy[1]<<","<<fMCData.dEpromz[1]
             <<")  Sp=("<<fMCData.dEspreadx[1]<<","<<fMCData.dEspready[1]<<","<<fMCData.dEspreadz[1]<<")\n";
  }

    LogTiming("ProcessMCParticles", section_start);
  }

  art::Handle< std::vector<recob::OpHit> > ophitListHandle;
  std::vector<art::Ptr<recob::OpHit> > ophitlist;

  for (size_t s = 0; s < fOpHitsModuleLabel.size(); s++) {
    e.getByLabel(fOpHitsModuleLabel[s], ophitListHandle);
    if(!ophitListHandle.isValid()){
      std::string msg = "OpHit not found: " + fOpHitsModuleLabel[s] + " (event " + std::to_string(fEventInfo.eventID) + ")";
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    if(fVerbosity>2) std::cout << "Found OpHits: " << fOpHitsModuleLabel[s] << std::endl;
    art::fill_ptr_vector(ophitlist, ophitListHandle);
  }

  LogTiming("ProcessOpHits", section_start);
  section_start = std::chrono::high_resolution_clock::now();

  // --- Process OpFlashes: Filter by beam window (0.367-1.9 μs) and extract flash positions + OpHits
  art::Handle< std::vector<recob::OpFlash> > opflashListHandle;

  for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) {
    e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);
    if(!opflashListHandle.isValid()){
      std::string msg = "OpFlash not found: " + fOpFlashesModuleLabel[s] + " (event " + std::to_string(fEventInfo.eventID) + ")";
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    art::FindManyP<recob::OpHit> flashToOpHitAssns(opflashListHandle, e, fOpFlashesModuleLabel[s]);

    if(fVerbosity>1) std::cout << "Processing OpFlashes: " << fOpFlashesModuleLabel[s] << std::endl;

    for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
      art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
      recob::OpFlash Flash = *FlashPtr;

      // Beam window filter: Reject out-of-time flashes (cosmic ray rejection)
      // MC: AbsTime in μs | DATA: AbsTime in ns → convert to μs?
      double flash_time_us = (fProcessingMode == "DATA_inference") ? Flash.AbsTime()/1000.0 : Flash.AbsTime();
      if(flash_time_us < 0.367 || flash_time_us > 1.9) continue;

      if(fVerbosity > 2) {
        std::cout << "Flash #" << fOpticalData.flash_x.size() << ": t=" << flash_time_us
                  << " μs, PE=" << Flash.TotalPE() << std::endl;
      }

      // Store flash position (x,y,z) in cm
      fOpticalData.flash_x.push_back( Flash.XCenter() );
      fOpticalData.flash_y.push_back( Flash.YCenter() );
      fOpticalData.flash_z.push_back( Flash.ZCenter() );

      if(fSaveOpHits){
        // Store OpHit data for CNN input: PE, channel, time per flash
        fOpticalData.flash_ophit_time.push_back({});
        fOpticalData.flash_ophit_pe.push_back({});
        fOpticalData.flash_ophit_ch.push_back({});

        std::vector<art::Ptr<recob::OpHit>> ophit_v = flashToOpHitAssns.at(i);
        for (auto ophit : ophit_v) {
          fOpticalData.flash_ophit_time.back().push_back(ophit->PeakTimeAbs());
          fOpticalData.flash_ophit_pe.back().push_back(ophit->PE());
          fOpticalData.flash_ophit_ch.back().push_back(ophit->OpChannel());
        }
      }
    }
  }

  LogTiming("ProcessOpFlashes", section_start);

  // --- Apply event filters: (1) Neutrino multiplicity (MC only), (2) OpHit availability
  bool passFilter1, passFilter2, passFilter;

  // Filter 1: Neutrino count (DATA always passes, MC requires exactly 1 neutrino unless skipped)
  if(fProcessingMode == "DATA_inference") {
    passFilter1 = true;
  } else {
    passFilter1 = fSkipNeutrinoFilter || (fMCData.nuvT.size() == 1);
    if(fVerbosity > 0 && !passFilter1) std::cout << "FILTER 1 FAIL: nuvT=" << fMCData.nuvT.size() << " (expected 1)" << std::endl;
  }

  // Filter 2: At least one flash with OpHits
  passFilter2 = false;
  if (!fOpticalData.flash_ophit_pe.empty()) {
    for (const auto& flash_pe : fOpticalData.flash_ophit_pe) {
      if (!flash_pe.empty()) {
        passFilter2 = true;
        break;
      }
    }
  }
  if(fVerbosity > 0 && !passFilter2) std::cout << "FILTER 2 FAIL: No flashes with OpHits" << std::endl;

  passFilter = passFilter1 && passFilter2;

  if(fVerbosity > 2) {
    std::cout << "Filter: nu=" << passFilter1 << " (nuvT=" << fMCData.nuvT.size()
              << "), optical=" << passFilter2 << " (flashes=" << fOpticalData.flash_ophit_pe.size()
              << "), overall=" << passFilter << std::endl;
  }

  LogTiming("ApplyFilters", section_start);

  // --- Process flashes: Channel classification → Energy filter → PE matrix → PE images → CNN inference
  if(passFilter) {
    if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();

    ApplyFlashSelection();  // Classify channels by TPC (even=TPC0, odd=TPC1) and PMT type
    ApplyFinalEnergyFilter();  // Apply energy/PE thresholds

    passFilter = !fProcessedData.flash_ophit_pe_final.empty();

    if(fVerbosity > 0) {
      std::cout << "Flash processing: " << (passFilter ? "PASS" : "FAIL")
                << " (" << fProcessedData.flash_ophit_pe_final.size() << " flashes retained)" << std::endl;
    }
    LogTiming("FlashSelection+EnergyFilter", section_start);

    if (passFilter) {
      if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();
      CreatePEMatrix();
      LogTiming("CreatePEMatrix", section_start);

      if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();
      CreatePEImages();
      LogTiming("CreatePEImages", section_start);
    }
  }

  // --- Run CNN inference if images were generated
  if (passFilter && !fProcessedData.pe_images.empty()) {
    if(fVerbosity > 0) section_start = std::chrono::high_resolution_clock::now();
    PixelMapVars tempPixelVars;
    RunInference(tempPixelVars);
    LogTiming("RunInference", section_start);
  }

  // --- Fill analysis TTree (always executed)
  FillInferenceTree(passFilter);

  // --- Optionally save PixelMapVars to output file (memory-intensive, disabled by default)
  if(fSavePixelMapVars) {
    auto pixelVars = std::make_unique<PixelMapVars>();
    FillPixelMapVarsConditional(*pixelVars, calculateGroundTruth, passFilter);

    if (passFilter && !fProcessedData.pe_images.empty()) {
      pixelVars->pe_images = fProcessedData.pe_images;
    }

    if(fVerbosity > 2) {
      std::cout << "PixelMapVars: " << pixelVars->flash_ophit_pe.size() << " flashes, "
                << pixelVars->pe_images.size() << " images" << std::endl;
    }

    e.put(std::move(pixelVars));
    if(fVerbosity > 1) std::cout << "PixelMapVars saved" << std::endl;
  } else {
    if(fVerbosity > 1) std::cout << "PixelMapVars NOT saved (SavePixelMapVars=false)" << std::endl;
  }

  if(fVerbosity > 0) {
    auto total_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - event_start).count();
    std::cout << "[TIMING] TOTAL: " << total_elapsed << "ms (Run " << fEventInfo.runID
              << ", Event " << fEventInfo.eventID << ")\n" << std::string(40,'=') << std::endl;
  }

}


// ============================================================================
// CreatePEMatrix - Convert flash/OpHit data to PE matrix
//   - Input: Multiple flashes with OpHits [flash][hit](PE, channel, time)
//   - Output: Single PE vector [312 channels] summing all flashes
// ============================================================================
void opdet::PosRecoCVNProducer::CreatePEMatrix() {
  fProcessedData.pe_matrix.clear();

  // Select data source: final filtered > selected > original
  const auto* pe_data = &fOpticalData.flash_ophit_pe;
  const auto* ch_data = &fOpticalData.flash_ophit_ch;

  if (!fProcessedData.flash_ophit_pe_final.empty() && !fProcessedData.flash_ophit_ch_final.empty()) {
    pe_data = &fProcessedData.flash_ophit_pe_final;
    ch_data = &fProcessedData.flash_ophit_ch_final;
    if(fVerbosity > 1) std::cout << "PE matrix: using final filtered data" << std::endl;
  } else if (!fProcessedData.flash_ophit_pe_sel.empty() && !fProcessedData.flash_ophit_ch_sel.empty()) {
    pe_data = &fProcessedData.flash_ophit_pe_sel;
    ch_data = &fProcessedData.flash_ophit_ch_sel;
    if(fVerbosity > 1) std::cout << "PE matrix: using selected data" << std::endl;
  } else {
    if(fVerbosity > 1) std::cout << "PE matrix: using original data" << std::endl;
  }

  fProcessedData.pe_matrix.resize(1, std::vector<float>(312, 0.0f));  // [1 event x 312 channels]

  if (pe_data->empty() || ch_data->empty()) {
    if(fVerbosity > 0) std::cout << "WARNING: No flash data for PE matrix" << std::endl;
    return;
  }

  // Accumulate PE from all flashes into single 312-channel vector
  for (size_t flash = 0; flash < pe_data->size(); ++flash) {
    if (flash >= ch_data->size()) {
      if(fVerbosity > 0) std::cout << "WARNING: PE/channel data size mismatch at flash " << flash << std::endl;
      break;
    }

    size_t min_size = std::min((*pe_data)[flash].size(), (*ch_data)[flash].size());
    for (size_t j = 0; j < min_size; ++j) {
      float pe = (*pe_data)[flash][j];
      int channel = (*ch_data)[flash][j];

      if (channel >= 0 && channel < 312) {
        fProcessedData.pe_matrix[0][channel] += pe;
      } else if (fVerbosity > 0) {
        std::cout << "WARNING: Invalid channel " << channel << " (flash " << flash << ", hit " << j << ")" << std::endl;
      }
    }
  }

  if(fVerbosity > 1) std::cout << "PE matrix: [1 x 312] created" << std::endl;
}


// ============================================================================
// LoadPMTMaps - Load channel-to-position mapping CSV files
//   - Coated PMTs: wavelength-shifting coating (VIS-sensitive)
//   - Uncoated PMTs: VUV-sensitive
//   - Maps define spatial layout [y][z] → channel for image generation
// ============================================================================
void opdet::PosRecoCVNProducer::LoadPMTMaps() {
  std::string coated_path = fCoatedPMTMapPath;
  std::string uncoated_path = fUncoatedPMTMapPath;

  // Search for coated PMT map in multiple fallback locations (Grid/MRB/local)
  std::ifstream test_coated(coated_path);
  if (!test_coated.is_open()) {
    size_t pos = coated_path.find_last_of("/\\");
    std::string filename = (pos != std::string::npos) ? coated_path.substr(pos + 1) : coated_path;

    std::vector<std::string> search_paths = {
      filename,
      "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename,
      "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/inference/module/" + filename,
      (getenv("MRB_INSTALL") ? std::string(getenv("MRB_INSTALL")) + "/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename : ""),
      (getenv("MRB_SOURCE") ? std::string(getenv("MRB_SOURCE")) + "/sbndcode/sbndcode/PosRecoCVN/pmt_maps/" + filename : ""),
      "../../../PosRecoCVN/pmt_maps/" + filename,
      "../../pmt_maps/" + filename,
      "sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename,
      "pmt_maps/" + filename
    };

    for (const auto& path : search_paths) {
      std::ifstream test_file(path);
      if (test_file.is_open()) {
        coated_path = path;
        test_file.close();
        break;
      }
    }
  }
  test_coated.close();

  // Search for uncoated PMT map
  std::ifstream test_uncoated(uncoated_path);
  if (!test_uncoated.is_open()) {
    size_t pos = uncoated_path.find_last_of("/\\");
    std::string filename = (pos != std::string::npos) ? uncoated_path.substr(pos + 1) : uncoated_path;

    std::vector<std::string> search_paths = {
      filename,
      "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename,
      "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/inference/module/" + filename,
      (getenv("MRB_INSTALL") ? std::string(getenv("MRB_INSTALL")) + "/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename : ""),
      (getenv("MRB_SOURCE") ? std::string(getenv("MRB_SOURCE")) + "/sbndcode/sbndcode/PosRecoCVN/pmt_maps/" + filename : ""),
      "../../../PosRecoCVN/pmt_maps/" + filename,
      "../../pmt_maps/" + filename,
      "sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename,
      "pmt_maps/" + filename
    };

    for (const auto& path : search_paths) {
      std::ifstream test_file(path);
      if (test_file.is_open()) {
        uncoated_path = path;
        test_file.close();
        break;
      }
    }
  }
  test_uncoated.close();

  std::ifstream coated_file(coated_path);
  std::ifstream uncoated_file(uncoated_path);

  if (!coated_file.is_open()) {
    throw std::runtime_error("Coated PMT map not found: " + coated_path);
  }
  if (!uncoated_file.is_open()) {
    throw std::runtime_error("Uncoated PMT map not found: " + uncoated_path);
  }

  fCoatedPMTMap.clear();
  fUncoatedPMTMap.clear();

  // Parse CSV files into 2D vectors
  std::string line;
  while (std::getline(coated_file, line)) {
    std::vector<int> row;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      row.push_back(std::stoi(cell));
    }
    fCoatedPMTMap.push_back(row);
  }
  coated_file.close();

  while (std::getline(uncoated_file, line)) {
    std::vector<int> row;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      row.push_back(std::stoi(cell));
    }
    fUncoatedPMTMap.push_back(row);
  }
  uncoated_file.close();

  if(fVerbosity > 0) {
    std::cout << "PMT maps loaded: coated[" << fCoatedPMTMap.size() << "x" << (fCoatedPMTMap.empty() ? 0 : fCoatedPMTMap[0].size())
              << "], uncoated[" << fUncoatedPMTMap.size() << "x" << (fUncoatedPMTMap.empty() ? 0 : fUncoatedPMTMap[0].size()) << "]" << std::endl;
  }
}


// ============================================================================
// SelectNonEmptyHalf - Select active TPC half based on signal distribution
//   - Methods: max (brightest pixel), sum (total PE), nonzero (hit count), mean_top (top-5 avg)
//   - Used to reduce image from full TPC to active half (top or bottom)
// ============================================================================
std::vector<std::vector<float>> opdet::PosRecoCVNProducer::SelectNonEmptyHalf(
    const std::vector<std::vector<float>>& left_half,
    const std::vector<std::vector<float>>& right_half,
    const std::string& method) {

  float left_score = 0.0f, right_score = 0.0f;

  if (method == "max") {  // Maximum pixel value
    for (const auto& row : left_half) {
      for (float val : row) left_score = std::max(left_score, val);
    }
    for (const auto& row : right_half) {
      for (float val : row) right_score = std::max(right_score, val);
    }
  }
  else if (method == "sum") {  // Total PE sum
    for (const auto& row : left_half) {
      for (float val : row) left_score += val;
    }
    for (const auto& row : right_half) {
      for (float val : row) right_score += val;
    }
  }
  else if (method == "nonzero") {  // Count non-zero pixels
    for (const auto& row : left_half) {
      for (float val : row) if (val != 0.0f) left_score += 1.0f;
    }
    for (const auto& row : right_half) {
      for (float val : row) if (val != 0.0f) right_score += 1.0f;
    }
  }
  else if (method == "mean_top") {  // Mean of top-5 pixels
    std::vector<float> left_vals, right_vals;
    for (const auto& row : left_half) {
      for (float val : row) left_vals.push_back(val);
    }
    for (const auto& row : right_half) {
      for (float val : row) right_vals.push_back(val);
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


// ============================================================================
// CreatePEImages - Generate CNN input images from PE matrix
//   1. Map 312 channels to spatial grids [ch_y, ch_z] for coated/uncoated PMTs
//   2. Normalize PE values to [0,1] range (with optional custom factor)
//   3. Select active TPC half (top/bottom based on max signal)
//   4. Output: [1, ch_y/2, ch_z, 2] images ready for CNN
// ============================================================================
void opdet::PosRecoCVNProducer::CreatePEImages() {
  fProcessedData.pe_images.clear();

  if (fProcessedData.pe_matrix.empty() || fCoatedPMTMap.empty() || fUncoatedPMTMap.empty()) {
    if(fVerbosity > 0) std::cout << "WARNING: Missing PE matrix or PMT maps" << std::endl;
    return;
  }

  int n_flashes = fProcessedData.pe_matrix.size();
  int n_events = 1;
  int ch_y = fCoatedPMTMap.size();
  int ch_z = fCoatedPMTMap[0].size();
  int map_count = 2;  // Coated + Uncoated

  // Map PE matrix to spatial images: [map_count][n_events][ch_y][ch_z]
  std::vector<std::vector<std::vector<std::vector<float>>>> pe_matrices_map(
    map_count, std::vector<std::vector<std::vector<float>>>(
      n_events, std::vector<std::vector<float>>(
        ch_y, std::vector<float>(ch_z, 0.0f))));

  std::vector<std::vector<std::vector<int>>*> maps = {&fUncoatedPMTMap, &fCoatedPMTMap};

  for (int idx = 0; idx < map_count; ++idx) {
    auto& map_ = *maps[idx];
    for (int y = 0; y < ch_y; ++y) {
      for (int z = 0; z < ch_z; ++z) {
        int channel = map_[y][z];
        if (channel >= 0 && channel < 312) {
          float total_pe = 0.0f;
          for (int flash = 0; flash < n_flashes; ++flash) {
            total_pe += fProcessedData.pe_matrix[flash][channel];
          }
          pe_matrices_map[idx][0][y][z] = total_pe;
        }
      }
    }
  }

  // Determine normalization factor
  float max_val_1 = 0.0f;
  for (int i = 0; i < n_events; ++i) {
    for (int y = 0; y < ch_y; ++y) {
      for (int z = 0; z < ch_z; ++z) {
        max_val_1 = std::max(max_val_1, std::max(pe_matrices_map[0][i][y][z], pe_matrices_map[1][i][y][z]));
      }
    }
  }

  if (fCustomNormFactor > 0) {
    max_val_1 = fCustomNormFactor;
    if(fVerbosity > 1) std::cout << "Norm factor: " << fCustomNormFactor << " (custom)" << std::endl;
  } else {
    max_val_1 = (max_val_1 > 0) ? max_val_1 : 1.0f;
    if(fVerbosity > 1) std::cout << "Norm factor: " << max_val_1 << " (auto)" << std::endl;
  }

  // Normalize and clip to [0, 1]
  int pixels_clipped = 0;
  float max_observed_value = 0.0f;

  for (int idx = 0; idx < map_count; ++idx) {
    for (int i = 0; i < n_events; ++i) {
      for (int y = 0; y < ch_y; ++y) {
        for (int z = 0; z < ch_z; ++z) {
          float original_value = pe_matrices_map[idx][i][y][z];
          pe_matrices_map[idx][i][y][z] = original_value / max_val_1;
          max_observed_value = std::max(max_observed_value, pe_matrices_map[idx][i][y][z]);

          if (pe_matrices_map[idx][i][y][z] > 1.0f) {
            pixels_clipped++;
            if(fVerbosity > 2 && pixels_clipped <= 5) {
              std::cout << "Clipping pixel (" << idx << "," << i << "," << y << "," << z
                        << "): " << pe_matrices_map[idx][i][y][z] << " → 1.0" << std::endl;
            }
            pe_matrices_map[idx][i][y][z] = 1.0f;
          }
        }
      }
    }
  }

  if (pixels_clipped > 0) {
    std::cout << "WARNING: " << pixels_clipped << " pixels clipped (max=" << max_observed_value
              << "). Use CustomNormFactor>=" << (max_val_1 * max_observed_value) << " to avoid clipping." << std::endl;
  } else if(fVerbosity > 1) {
    std::cout << "Normalization OK: max=" << max_observed_value << std::endl;
  }

  // Select active TPC half (top/bottom) and build final images: [n_events][ch_y/2][ch_z][map_count]
  fProcessedData.pe_images.resize(n_events, std::vector<std::vector<std::vector<float>>>(
    ch_y/2, std::vector<std::vector<float>>(
      ch_z, std::vector<float>(map_count, 0.0f))));

  for (int i = 0; i < n_events; ++i) {
    for (int idx = 0; idx < map_count; ++idx) {
      const auto& event_matrix = pe_matrices_map[idx][i];
      std::vector<std::vector<float>> top_half(event_matrix.begin(), event_matrix.begin() + ch_y/2);
      std::vector<std::vector<float>> bottom_half(event_matrix.begin() + ch_y/2, event_matrix.end());

      auto selected_half = SelectNonEmptyHalf(top_half, bottom_half, "max");

      for (int y = 0; y < ch_y/2; ++y) {
        for (int z = 0; z < ch_z; ++z) {
          fProcessedData.pe_images[i][y][z][idx] = selected_half[y][z];
        }
      }
    }
  }

  if(fVerbosity > 1) {
    std::cout << "PE images: [" << n_events << ", " << ch_y/2 << ", " << ch_z << ", " << map_count << "]" << std::endl;
  }
}

// ============================================================================
// RunInference - Execute TensorFlow CNN model
//   - Input: PE images [1, ch_y/2, ch_z, 2]
//   - CNN output: Normalized coordinates (x,y,z) in [0,1] or [-1,1]
//   - Apply inverse scaling to get real coordinates [cm]
//   - Calculate 3D error vs ground truth (MC_testing mode only)
// ============================================================================
void opdet::PosRecoCVNProducer::RunInference(PixelMapVars& pixelmapvars) {

  if (!fTFGraph || fProcessedData.pe_images.empty()) {
    if(fVerbosity > 1) std::cout << "Skipping inference (graph=" << (fTFGraph != nullptr)
                                 << ", images=" << !fProcessedData.pe_images.empty() << ")" << std::endl;
    return;
  }

  try {
    auto predictions = fTFGraph->run(fProcessedData.pe_images, -1);

    if(fVerbosity > 2) {
      std::cout << "Predictions shape: [" << predictions.size()
                << ", " << (predictions.empty() ? 0 : predictions[0].size())
                << ", " << (predictions.empty() || predictions[0].empty() ? 0 : predictions[0][0].size()) << "]" << std::endl;
    }

    if (predictions.empty() || predictions[0].empty()) {
      std::cout << "ERROR: Empty predictions from TensorFlow" << std::endl;
      return;
    }

    // Extract (x,y,z) predictions for each event
    size_t max_events = predictions.size();
    for (size_t i = 0; i < max_events; ++i) {
      if (predictions[i].empty()) continue;

      // Extract raw CNN outputs: predictions[event][0][x, y, z]
      std::vector<double> raw_predictions;
      if (!predictions[i][0].empty()) {
        for (size_t k = 0; k < std::min((size_t)3, predictions[i][0].size()); ++k) {
          raw_predictions.push_back(predictions[i][0][k]);
        }
      }

      if (raw_predictions.size() < 3) continue;

      // Apply inverse scaling to convert normalized outputs to real coordinates [cm]
      std::vector<double> unscaled = ApplyInverseScaling(raw_predictions);

      if(fVerbosity > 1) {
        std::cout << "Scaled [" << raw_predictions[0] << ", " << raw_predictions[1] << ", " << raw_predictions[2]
                  << "] → Unscaled [" << unscaled[0] << ", " << unscaled[1] << ", " << unscaled[2] << "]" << std::endl;
      }

      // Store predictions (X is always positive by taking absolute value)
      fTreeData.predX = std::abs(unscaled[0]);
      fTreeData.predY = unscaled[1];
      fTreeData.predZ = unscaled[2];

      // Calculate errors if ground truth available (MC_testing mode only)
      if (fProcessingMode == "MC_testing" && i < fProcessedData.dEpromx_final.size()
          && i < fProcessedData.dEpromy_final.size() && i < fProcessedData.dEpromz_final.size()) {

        double true_x = std::abs(fProcessedData.dEpromx_final[i]);
        double true_y = fProcessedData.dEpromy_final[i];
        double true_z = fProcessedData.dEpromz_final[i];

        fTreeData.diffX = fTreeData.predX - true_x;
        fTreeData.diffY = fTreeData.predY - true_y;
        fTreeData.diffZ = fTreeData.predZ - true_z;

        fTreeData.error3D = std::sqrt(fTreeData.diffX*fTreeData.diffX +
                                       fTreeData.diffY*fTreeData.diffY +
                                       fTreeData.diffZ*fTreeData.diffZ);
        pixelmapvars.error_3d.push_back(fTreeData.error3D);

        if(fVerbosity > 0) {
          std::cout << "R" << fEventInfo.runID << " E" << fEventInfo.eventID
                    << " True(" << true_x << "," << true_y << "," << true_z << ")"
                    << " Pred(" << fTreeData.predX << "," << fTreeData.predY << "," << fTreeData.predZ << ")"
                    << " Err3D=" << fTreeData.error3D << " cm" << std::endl;
        }
      } else {
        fTreeData.diffX = fTreeData.diffY = fTreeData.diffZ = fTreeData.error3D = -999.0;
        pixelmapvars.error_3d.push_back(-999.0);

        if(fVerbosity > 0) {
          std::cout << "R" << fEventInfo.runID << " E" << fEventInfo.eventID
                    << " Pred(" << fTreeData.predX << "," << fTreeData.predY << "," << fTreeData.predZ << ")" << std::endl;
        }
      }
    }

  } catch (const std::exception& e) {
    std::cerr << "ERROR: TensorFlow inference failed: " << e.what() << "\n"
              << "Check: input dimensions, model file, TF version compatibility" << std::endl;
  }
}


// ============================================================================
// ApplyInverseScaling - Convert normalized CNN outputs to detector coordinates
//   - X: [0,1] → [0,200] cm (drift distance)
//   - Y: [-1,1] → [-200,200] cm (vertical position)
//   - Z: [0,1] → [0,500] cm (beam direction)
//   - Validates predictions with configurable tolerance, clips out-of-range values
// ============================================================================
std::vector<double> opdet::PosRecoCVNProducer::ApplyInverseScaling(const std::vector<double>& scaled) {

  if (scaled.size() != 3) {
    std::cout << "WARNING: Expected 3 predictions, got " << scaled.size() << std::endl;
    return scaled;
  }

  std::vector<double> validated(3);
  const double tol = fPredictionTolerance;

  // Validate and clip X: [0,1] ± tolerance
  validated[0] = scaled[0];
  if (scaled[0] < -tol) {
    if(fVerbosity > 0) std::cout << "WARNING: X=" << scaled[0] << " clipped to 0.0" << std::endl;
    validated[0] = 0.0;
  } else if (scaled[0] > 1.0 + tol) {
    if(fVerbosity > 0) std::cout << "WARNING: X=" << scaled[0] << " clipped to 1.0" << std::endl;
    validated[0] = 1.0;
  }

  // Validate and clip Y: [-1,1] ± tolerance
  validated[1] = scaled[1];
  if (scaled[1] < -1.0 - tol) {
    if(fVerbosity > 0) std::cout << "WARNING: Y=" << scaled[1] << " clipped to -1.0" << std::endl;
    validated[1] = -1.0;
  } else if (scaled[1] > 1.0 + tol) {
    if(fVerbosity > 0) std::cout << "WARNING: Y=" << scaled[1] << " clipped to 1.0" << std::endl;
    validated[1] = 1.0;
  }

  // Validate and clip Z: [0,1] ± tolerance
  validated[2] = scaled[2];
  if (scaled[2] < -tol) {
    if(fVerbosity > 0) std::cout << "WARNING: Z=" << scaled[2] << " clipped to 0.0" << std::endl;
    validated[2] = 0.0;
  } else if (scaled[2] > 1.0 + tol) {
    if(fVerbosity > 0) std::cout << "WARNING: Z=" << scaled[2] << " clipped to 1.0" << std::endl;
    validated[2] = 1.0;
  }

  // Apply inverse scaling: X:[0,200], Y:[-200,200], Z:[0,500]
  std::vector<double> unscaled(3);
  unscaled[0] = validated[0] * 200.0;
  unscaled[1] = validated[1] * 200.0;
  unscaled[2] = validated[2] * 500.0;

  if(fVerbosity > 2) {
    std::cout << "Unscale: X=" << validated[0] << "→" << unscaled[0]
              << ", Y=" << validated[1] << "→" << unscaled[1]
              << ", Z=" << validated[2] << "→" << unscaled[2] << std::endl;
  }

  return unscaled;
}

// ============================================================================
// FillInferenceTree - Fill output TTree with event data and CNN predictions
// ============================================================================
void opdet::PosRecoCVNProducer::FillInferenceTree(bool passedFilters)
{
  // Fill basic event information for all events
  fTreeData.run = fEventInfo.runID;
  fTreeData.subrun = fEventInfo.subrunID;
  fTreeData.event = fEventInfo.eventID;
  fTreeData.passedFilters = passedFilters;

  // Initialize variables with default values
  // NOTE: fTreePred* variables are filled directly in RunInference() - do NOT reinitialize here
  fTreeData.error3D = -999.0;
  fTreeData.nuvT = -999.0;
  fTreeData.nuvZ = -999.0;
  fTreeData.dEtpc = -999.0;

  // Truth and difference variables only meaningful in MC_testing mode
  if(fProcessingMode == "MC_testing") {
    fTreeData.trueX = -999.0;
    fTreeData.trueY = -999.0;
    fTreeData.trueZ = -999.0;
    // fTreeData.diffX/Y/Z are filled directly in RunInference(), don't reinitialize
  } else {
    // In inference modes, set truth and diff to -999 (not applicable)
    fTreeData.trueX = -999.0;
    fTreeData.trueY = -999.0;
    fTreeData.trueZ = -999.0;
    fTreeData.diffX = -999.0;
    fTreeData.diffY = -999.0;
    fTreeData.diffZ = -999.0;
  }

  // Fill physics data if event passed filters
  if (passedFilters) {
    // Fill ground truth data only in MC_testing mode
    if(fProcessingMode == "MC_testing") {
      if (!fProcessedData.dEpromx_final.empty()) {
        fTreeData.trueX = fProcessedData.dEpromx_final[0];  // Note: already has abs() applied
        fTreeData.trueY = fProcessedData.dEpromy_final[0];
        fTreeData.trueZ = fProcessedData.dEpromz_final[0];
      }

      if (!fProcessedData.dEtpc_final.empty()) {
        fTreeData.dEtpc = fProcessedData.dEtpc_final[0];
      }

      if (!fProcessedData.nuvT_final.empty()) {
        fTreeData.nuvT = fProcessedData.nuvT_final[0];
      }

      if (!fProcessedData.nuvZ_final.empty()) {
        fTreeData.nuvZ = fProcessedData.nuvZ_final[0];
      }

      // Calculate 3D error if predictions are available
      if(fTreeData.predX != -999.0 && fTreeData.predY != -999.0 && fTreeData.predZ != -999.0 &&
         fTreeData.trueX != -999.0 && fTreeData.trueY != -999.0 && fTreeData.trueZ != -999.0) {
        double dx = fTreeData.predX - fTreeData.trueX;
        double dy = fTreeData.predY - fTreeData.trueY;
        double dz = fTreeData.predZ - fTreeData.trueZ;
        fTreeData.error3D = std::sqrt(dx*dx + dy*dy + dz*dz);
      }
    }

    // Prediction data (fTreeData.predX/Y/Z) is filled directly in RunInference()
  }

  // Fill the TTree entry
  fInferenceTree->Fill();

  if(fVerbosity > 1) {
    std::cout << "DEBUG: FillInferenceTree (" << fProcessingMode << " mode) - PredX=" << fTreeData.predX
              << " PredY=" << fTreeData.predY << " PredZ=" << fTreeData.predZ << std::endl;
    if(fProcessingMode == "testing") {
      std::cout << "DEBUG: FillInferenceTree - TrueX=" << fTreeData.trueX << " TrueY=" << fTreeData.trueY << " TrueZ=" << fTreeData.trueZ << std::endl;
      std::cout << "DEBUG: FillInferenceTree - DiffX=" << fTreeData.diffX << " DiffY=" << fTreeData.diffY << " DiffZ=" << fTreeData.diffZ << std::endl;
    }
  }

  if(fVerbosity > 2) {
    std::cout << "Filled TTree entry for Run=" << fTreeData.run << " Event=" << fTreeData.event
              << " PassedFilters=" << fTreeData.passedFilters << " Mode=" << fProcessingMode << std::endl;
  }
}

// ============================================================================
// FillPixelMapVarsConditional - Conditionally copy data to PixelMapVars output
// ============================================================================
void opdet::PosRecoCVNProducer::FillPixelMapVarsConditional(PixelMapVars& pixelVars, bool calculateGroundTruth, bool passFilter)
{
  // Always include the channel dictionary (independent of filter)
  pixelVars.channel_dict = fChannelDict;

  // Fill basic event information
  pixelVars.run_id = fEventInfo.runID;
  pixelVars.subrun_id = fEventInfo.subrunID;
  pixelVars.event_id = fEventInfo.eventID;
  pixelVars.passed_filters = passFilter;

  if(passFilter) {
    // Apply mask to all variables - use final filtered data after all cuts
    pixelVars.flash_ophit_pe = fProcessedData.flash_ophit_pe_final;
    pixelVars.flash_ophit_ch = fProcessedData.flash_ophit_ch_final;
    pixelVars.flash_ophit_time = fProcessedData.flash_ophit_time_final;

    if (calculateGroundTruth) {
      // Include ground truth data in testing mode
      pixelVars.nuvT = fProcessedData.nuvT_final;
      pixelVars.dEpromx = fProcessedData.dEpromx_final;
      pixelVars.dEpromy = fProcessedData.dEpromy_final;
      pixelVars.dEpromz = fProcessedData.dEpromz_final;
      pixelVars.dEtpc = fProcessedData.dEtpc_final;
      pixelVars.nuvZ = fProcessedData.nuvZ_final;

      // Extended MC Truth information
      pixelVars.nuvX = fMCData.nuvX;
      pixelVars.nuvY = fMCData.nuvY;
      pixelVars.nuvE = fMCData.nuvE;
    } else {
      // In MC_inference or DATA_inference mode, initialize ground truth vectors as empty
      pixelVars.nuvT.clear();
      pixelVars.dEpromx.clear();
      pixelVars.dEpromy.clear();
      pixelVars.dEpromz.clear();
      pixelVars.dEtpc.clear();
      pixelVars.nuvZ.clear();

      // Extended MC Truth information (empty in inference modes)
      pixelVars.nuvX.clear();
      pixelVars.nuvY.clear();
      pixelVars.nuvE.clear();
    }

    // CNN predictions and differences (filled for all modes if inference was run)
    pixelVars.dEpromx_pred.clear();
    pixelVars.dEpromy_pred.clear();
    pixelVars.dEpromz_pred.clear();
    pixelVars.dEpromx_diff.clear();
    pixelVars.dEpromy_diff.clear();
    pixelVars.dEpromz_diff.clear();

    pixelVars.dEpromx_pred.push_back(fTreeData.predX);
    pixelVars.dEpromy_pred.push_back(fTreeData.predY);
    pixelVars.dEpromz_pred.push_back(fTreeData.predZ);
    pixelVars.dEpromx_diff.push_back(fTreeData.diffX);
    pixelVars.dEpromy_diff.push_back(fTreeData.diffY);
    pixelVars.dEpromz_diff.push_back(fTreeData.diffZ);

    // PE matrix and PMT maps (for debugging and detailed analysis)
    pixelVars.pe_matrix = fProcessedData.pe_matrix;
    pixelVars.coated_pmt_map = fCoatedPMTMap;
    pixelVars.uncoated_pmt_map = fUncoatedPMTMap;

    if(fVerbosity > 1) {
      std::cout << "Event passed filter - data stored in PixelMapVars" << std::endl;
    }
  } else {
    // Initialize empty vectors for events that don't pass filter
    pixelVars.flash_ophit_pe.clear();
    pixelVars.flash_ophit_ch.clear();
    pixelVars.flash_ophit_time.clear();
    pixelVars.nuvT.clear();
    pixelVars.dEpromx.clear();
    pixelVars.dEpromy.clear();
    pixelVars.dEpromz.clear();
    pixelVars.dEtpc.clear();
    pixelVars.nuvZ.clear();

    // Extended MC Truth information (empty for failed events)
    pixelVars.nuvX.clear();
    pixelVars.nuvY.clear();
    pixelVars.nuvE.clear();

    // CNN predictions and differences (empty for failed events)
    pixelVars.dEpromx_pred.clear();
    pixelVars.dEpromy_pred.clear();
    pixelVars.dEpromz_pred.clear();
    pixelVars.dEpromx_diff.clear();
    pixelVars.dEpromy_diff.clear();
    pixelVars.dEpromz_diff.clear();

    // PE matrix and PMT maps (empty for failed events)
    pixelVars.pe_matrix.clear();
    pixelVars.coated_pmt_map.clear();
    pixelVars.uncoated_pmt_map.clear();

    if(fVerbosity > 0) {
      std::cout << "Run=" << fEventInfo.runID << " Subrun=" << fEventInfo.subrunID << " Event=" << fEventInfo.eventID << " FAILED FILTER" << std::endl;
    }
  }
}


// ============================================================================
// FillMCTruth - Extract and store MC truth data (neutrinos, particles)
// ============================================================================
void opdet::PosRecoCVNProducer::FillMCTruth(art::Event const& e){
  
  
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
                       " not found or empty in event " + std::to_string(fEventInfo.eventID);
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

          if(fVerbosity>1){
            std::cout << "    " << par.TrackId() << "  Particle PDG: " << par.PdgCode() << " E: " << par.E() << " t: " << par.T();
            std::cout << " Mother: "<<par.Mother() << " Process: "<<par.Process()<<" Status: "<<par.StatusCode()<<std::endl;
          }
          
          // Only save vertex if the PDG is specified in the fhicl list
          if( find (fMCTruthPDG.begin(), fMCTruthPDG.end(), par.PdgCode() ) != fMCTruthPDG.end() ){
       
            if(par.StatusCode()==0 && evtTruth->Origin()==1){
              nu_x=par.Vx();
              nu_y=par.Vy();
              nu_z=par.Vz();
              nu_t=par.T();
              nu_E=par.E();
            }
           
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

      fMCData.nuvT.push_back(nu_t);
      fMCData.nuvX.push_back(nu_x);
      fMCData.nuvY.push_back(nu_y);
      fMCData.nuvZ.push_back(nu_z);
      fMCData.nuvE.push_back(nu_E);

    }
  }
}


// ============================================================================
// FillAverageDepositedEnergyVariables - Compute energy-weighted centroids
// ============================================================================
void opdet::PosRecoCVNProducer::FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
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
    dEpromx_tpc0=std::abs(dEpromx_tpc0/dE_tpc0);
    dEpromy_tpc0=dEpromy_tpc0/dE_tpc0;
    dEpromz_tpc0=dEpromz_tpc0/dE_tpc0;
    spreadx_tpc0=sqrt( spreadx_tpc0/dE_tpc0-dEpromx_tpc0*dEpromx_tpc0 );
    spready_tpc0=std::sqrt(spready_tpc0/dE_tpc0-dEpromy_tpc0*dEpromy_tpc0);
    spreadz_tpc0=std::sqrt(spreadz_tpc0/dE_tpc0-dEpromz_tpc0*dEpromz_tpc0);
    dEtpc[0]=dE_tpc0;dEpromx[0]=dEpromx_tpc0;dEpromy[0]=dEpromy_tpc0;dEpromz[0]=dEpromz_tpc0;
    dEspreadx[0]=spreadx_tpc0;dEspready[0]=spready_tpc0;dEspreadz[0]=spreadz_tpc0;
  }
  if(ndeps_tpc1!=0){
    dEpromx_tpc1=std::abs(dEpromx_tpc1/dE_tpc1);
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


// ============================================================================
// InitializeChannelDict - Map OpDet channels to PMT types
// ============================================================================
void opdet::PosRecoCVNProducer::InitializeChannelDict(){
  
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


// ============================================================================
// ClassifyChannels - Separate channels by TPC (even/odd) and PMT type
// ============================================================================
void opdet::PosRecoCVNProducer::ClassifyChannels(){
  
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


// ============================================================================
// GetDetectorType - Determine PMT type from OpDet channel ID
// ============================================================================
int opdet::PosRecoCVNProducer::CategorizeFirstChannel(const std::vector<int>& channels){
  
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
std::vector<std::vector<T>> opdet::PosRecoCVNProducer::FilterByMask(
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
// ============================================================================
// ApplyFlashSelection - Filter flashes and classify OpHits by TPC/PMT type
// ============================================================================
void opdet::PosRecoCVNProducer::ApplyFlashSelection(){
  
  // Clear output vectors
  fProcessedData.flash_ophit_pe_sel.clear();
  fProcessedData.flash_ophit_ch_sel.clear(); 
  fProcessedData.flash_ophit_time_sel.clear();
  fProcessedData.categorized_flashes.clear();
  fProcessedData.dEpromx_sel.clear();
  fProcessedData.dEpromy_sel.clear();
  fProcessedData.dEpromz_sel.clear();
  fProcessedData.dEtpc_sel.clear();
  
  if(fOpticalData.flash_ophit_ch.empty()) {
    if(fVerbosity > 0) std::cout << "No flashes to process" << std::endl;
    return;
  }
  
  // 1. Categorize flashes based on first channel
  std::vector<int> categorized;
  for(const auto& flash_ch : fOpticalData.flash_ophit_ch) {
    categorized.push_back(CategorizeFirstChannel(flash_ch));
  }
  
  // 2. Calculate sum PE per flash
  std::vector<float> sum_pe;
  for(const auto& flash_pe : fOpticalData.flash_ophit_pe) {
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
  std::vector<std::vector<std::vector<float>>> pe_2d = {fOpticalData.flash_ophit_pe};
  std::vector<std::vector<std::vector<int>>> ch_2d = {fOpticalData.flash_ophit_ch};
  std::vector<std::vector<std::vector<float>>> time_2d = {fOpticalData.flash_ophit_time};
  std::vector<std::vector<bool>> mask_2d = {selected_mask};
  
  auto pe_filtered = FilterByMask(pe_2d, mask_2d);
  auto ch_filtered = FilterByMask(ch_2d, mask_2d);  
  auto time_filtered = FilterByMask(time_2d, mask_2d);
  
  if(!pe_filtered.empty()) fProcessedData.flash_ophit_pe_sel = pe_filtered[0];
  if(!ch_filtered.empty()) fProcessedData.flash_ophit_ch_sel = ch_filtered[0];
  if(!time_filtered.empty()) fProcessedData.flash_ophit_time_sel = time_filtered[0];
  
  // 6. TPC selection - only in MC_testing mode where ground truth is available
  if(fProcessingMode == "MC_testing" && !fMCData.dEpromx.empty()) {
    std::vector<std::vector<double>> selector = {{decision ? 1.0 : 0.0, decision ? 0.0 : 1.0}};
    std::vector<std::vector<double>> dEpromx_2d = {fMCData.dEpromx};
    std::vector<std::vector<double>> dEpromy_2d = {fMCData.dEpromy};
    std::vector<std::vector<double>> dEpromz_2d = {fMCData.dEpromz};
    std::vector<std::vector<double>> dEtpc_2d = {fMCData.dEtpc};

    // Select TPC values based on decision with bounds checking
    size_t max_tpc = std::min({(size_t)2, fMCData.dEpromx.size(), fMCData.dEpromy.size(), fMCData.dEpromz.size(), fMCData.dEtpc.size()});
    for(size_t i = 0; i < max_tpc; ++i) {
      if(selector[0][i] > 0.5) {
        fProcessedData.dEpromx_sel.push_back(fMCData.dEpromx[i]);
        fProcessedData.dEpromy_sel.push_back(fMCData.dEpromy[i]);
        fProcessedData.dEpromz_sel.push_back(fMCData.dEpromz[i]);
        fProcessedData.dEtpc_sel.push_back(fMCData.dEtpc[i]);
        break; // Only select one TPC
      }
    }
  } else if(fProcessingMode != "MC_testing") {
    // Inference modes: No ground truth data to select, clear the vectors
    fProcessedData.dEpromx_sel.clear();
    fProcessedData.dEpromy_sel.clear();
    fProcessedData.dEpromz_sel.clear();
    fProcessedData.dEtpc_sel.clear();
  }
  
  if(fVerbosity > 0) {
    std::cout << "Flash selection: " << n_flashes << " flashes, decision=" << decision 
              << " (even=" << sum_even << ", odd=" << sum_odd << ")" << std::endl;
    std::cout << "Selected " << fProcessedData.flash_ophit_pe_sel.size() << " flashes" << std::endl;
  }
}


// ============================================================================
// ApplyFinalEnergyFilter - Apply PE/energy thresholds to select best flash
// ============================================================================
void opdet::PosRecoCVNProducer::ApplyFinalEnergyFilter(){

  // Clear final output vectors
  fProcessedData.flash_ophit_pe_final.clear();
  fProcessedData.flash_ophit_ch_final.clear();
  fProcessedData.flash_ophit_time_final.clear();
  fProcessedData.nuvT_final.clear();
  fProcessedData.nuvZ_final.clear();
  fProcessedData.dEpromx_final.clear();
  fProcessedData.dEpromy_final.clear();
  fProcessedData.dEpromz_final.clear();
  fProcessedData.dEtpc_final.clear();

  bool passFilter = false;

  // Different filtering logic for MC_testing vs MC_inference vs DATA_inference mode
  if(fProcessingMode == "MC_testing") {
    // MC_testing mode: Use ground truth energy filter (original logic)

    // Check if we have valid selected data
    if(fProcessedData.dEpromx_sel.empty() || fProcessedData.dEpromy_sel.empty() ||
       fProcessedData.dEpromz_sel.empty() || fProcessedData.dEtpc_sel.empty()) {
      if(fVerbosity > 0) {
        std::cout << "FAILED FILTER: Final energy filter - No selected energy data available" << std::endl;
      }
      return;
    }

    // Apply energy deposition mask:
    // (dEpromx != -999) & (dEpromy != -999) & (dEpromz != -999) & (dEtpc > 50)
    // + position cuts: dEpromx in (-200,200), dEpromy in (-200,200), dEpromz in (0,500)

    // Check all selected TPC data (typically just one after flash selection)
    for(size_t i = 0; i < fProcessedData.dEpromx_sel.size(); ++i) {
      bool validX = (fProcessedData.dEpromx_sel[i] != fDefaultSimIDE);
      bool validY = (fProcessedData.dEpromy_sel[i] != fDefaultSimIDE);
      bool validZ = (fProcessedData.dEpromz_sel[i] != fDefaultSimIDE);
      bool energyCut = (fProcessedData.dEtpc_sel[i] > 50.0);

      // Position cuts - matching training data filters
      bool positionCutX = (fProcessedData.dEpromx_sel[i] >= -200.0 && fProcessedData.dEpromx_sel[i] <= 200.0);
      bool positionCutY = (fProcessedData.dEpromy_sel[i] >= -200.0 && fProcessedData.dEpromy_sel[i] <= 200.0);
      bool positionCutZ = (fProcessedData.dEpromz_sel[i] >= 0.0 && fProcessedData.dEpromz_sel[i] <= 500.0);

      if(validX && validY && validZ && energyCut && positionCutX && positionCutY && positionCutZ) {
        passFilter = true;
        // Store the passing values
        fProcessedData.dEpromx_final.push_back(fProcessedData.dEpromx_sel[i]);
        fProcessedData.dEpromy_final.push_back(fProcessedData.dEpromy_sel[i]);
        fProcessedData.dEpromz_final.push_back(fProcessedData.dEpromz_sel[i]);
        fProcessedData.dEtpc_final.push_back(fProcessedData.dEtpc_sel[i]);

        if(fVerbosity > 0) {
          std::cout << "Energy and position filters passed: dE=" << fProcessedData.dEtpc_sel[i]
                    << " MeV, pos=(" << fProcessedData.dEpromx_sel[i] << ","
                    << fProcessedData.dEpromy_sel[i] << "," << fProcessedData.dEpromz_sel[i] << ")" << std::endl;
          std::cout << "  Position cuts: X(" << fProcessedData.dEpromx_sel[i] << " in [-200,200]), "
                    << "Y(" << fProcessedData.dEpromy_sel[i] << " in [-200,200]), "
                    << "Z(" << fProcessedData.dEpromz_sel[i] << " in [0,500])" << std::endl;
        }
        break; // Only need one passing entry
      } else if(fVerbosity > 0) {
        // Log why this entry failed
        std::cout << "FAILED FILTER: Energy/position filter - dE=" << fProcessedData.dEtpc_sel[i]
                  << " MeV (need >50), pos=(" << fProcessedData.dEpromx_sel[i] << ","
                  << fProcessedData.dEpromy_sel[i] << "," << fProcessedData.dEpromz_sel[i] << ")" << std::endl;
        std::cout << "  Reasons: validX=" << validX << ", validY=" << validY << ", validZ=" << validZ
                  << ", energyCut=" << energyCut << ", posX=" << positionCutX
                  << ", posY=" << positionCutY << ", posZ=" << positionCutZ << std::endl;
      }
    }
  } else if(fProcessingMode == "MC_inference") {
    // MC_inference mode: Use PE-based quality filter (same as DATA_inference)

    if(fProcessedData.flash_ophit_pe_sel.empty()) {
      if(fVerbosity > 0) {
        std::cout << "FAILED FILTER: MC_inference PE filter - No selected flash data available" << std::endl;
      }
      return;
    }

    // Calculate total PE from selected flashes
    double total_pe = 0.0;
    int num_channels = 0;

    for(const auto& flash_pe : fProcessedData.flash_ophit_pe_sel) {
      for(float pe : flash_pe) {
        total_pe += pe;
        num_channels++;
      }
    }

    // Quality filter based on PE thresholds (same as DATA_inference)
    double min_pe_threshold = 100.0;
    int min_channels = 3;

    bool pe_cut = (total_pe > min_pe_threshold);
    bool channel_cut = (num_channels >= min_channels);

    if(fVerbosity > 0) {
      std::cout << "MC_inference PE filter evaluation: total_pe=" << total_pe << ", channels=" << num_channels << std::endl;
      std::cout << "  Thresholds: min_pe=" << min_pe_threshold << ", min_channels=" << min_channels << std::endl;
      std::cout << "  Results: pe_cut=" << pe_cut << ", channel_cut=" << channel_cut << std::endl;
    }

    if(pe_cut && channel_cut) {
      passFilter = true;

      if(fVerbosity > 0) {
        std::cout << "MC_inference PE filter PASSED: total_pe=" << total_pe
                  << " (>" << min_pe_threshold << "), channels=" << num_channels
                  << " (>=" << min_channels << ")" << std::endl;
      }
    } else {
      if(fVerbosity > 0) {
        std::cout << "FAILED FILTER: MC_inference PE filter - total_pe=" << total_pe
                  << " (need >" << min_pe_threshold << "), channels=" << num_channels
                  << " (need >=" << min_channels << ")" << std::endl;
      }
    }
  } else if(fProcessingMode == "DATA_inference") {
    // DATA_inference mode: Simplified filter with only PE and channel cuts (no beam window, no other filters)

    if(fProcessedData.flash_ophit_pe_sel.empty()) {
      if(fVerbosity > 0) {
        std::cout << "FAILED FILTER: DATA PE filter - No selected flash data available" << std::endl;
      }
      return;
    }

    // Calculate total PE from selected flashes
    double total_pe = 0.0;
    int num_channels = 0;

    for(const auto& flash_pe : fProcessedData.flash_ophit_pe_sel) {
      for(float pe : flash_pe) {
        total_pe += pe;
        num_channels++;
      }
    }

    // Simplified quality filter for data: only PE and channel thresholds
    double min_pe_threshold = 100.0;  // Same as MC_inference
    int min_channels = 3;  // Same as MC_inference

    bool pe_cut = (total_pe > min_pe_threshold);
    bool channel_cut = (num_channels >= min_channels);

    if(fVerbosity > 0) {
      std::cout << "DATA PE filter evaluation: total_pe=" << total_pe << ", channels=" << num_channels << std::endl;
      std::cout << "  Thresholds: min_pe=" << min_pe_threshold << ", min_channels=" << min_channels << std::endl;
      std::cout << "  Results: pe_cut=" << pe_cut << ", channel_cut=" << channel_cut << std::endl;
    }

    if(pe_cut && channel_cut) {
      passFilter = true;

      if(fVerbosity > 0) {
        std::cout << "DATA PE filter PASSED: total_pe=" << total_pe
                  << " (>" << min_pe_threshold << "), channels=" << num_channels
                  << " (>=" << min_channels << ")" << std::endl;
      }
    } else {
      if(fVerbosity > 0) {
        std::cout << "FAILED FILTER: DATA PE filter - total_pe=" << total_pe
                  << " (need >" << min_pe_threshold << "), channels=" << num_channels
                  << " (need >=" << min_channels << ")" << std::endl;
      }
    }
  }

  if(passFilter) {
    // Copy flash data if filter passes
    fProcessedData.flash_ophit_pe_final = fProcessedData.flash_ophit_pe_sel;
    fProcessedData.flash_ophit_ch_final = fProcessedData.flash_ophit_ch_sel;
    fProcessedData.flash_ophit_time_final = fProcessedData.flash_ophit_time_sel;

    // Copy neutrino data (may be empty in inference mode)
    fProcessedData.nuvT_final = fMCData.nuvT;
    fProcessedData.nuvZ_final = fMCData.nuvZ;

    if(fVerbosity > 0) {
      std::cout << "Final filter (" << fProcessingMode << " mode): Event PASSED - "
                << fProcessedData.flash_ophit_pe_final.size() << " flashes retained" << std::endl;
    }
  } else {
    if(fVerbosity > 0) {
      std::cout << "Final filter (" << fProcessingMode << " mode): Event FAILED" << std::endl;
    }
  }
}

// ============================================================================
// ClearEventData - Reset all event-level data structures
// ============================================================================
void opdet::PosRecoCVNProducer::ClearEventData(){
  // Optimized: Clear all data with 5 calls instead of 80+
  fEventInfo.clear();
  fMCData.clear();
  fOpticalData.clear();
  fProcessedData.clear();
  fTreeData.resetToDefaults();

  // Note: TTree variables are handled by fTreeData.resetToDefaults()
  // Channel dictionaries and TensorFlow objects are persistent across events
}

// ============================================================================
// LogTiming - Print elapsed time for operation (if verbosity enabled)
// ============================================================================
void opdet::PosRecoCVNProducer::LogTiming(const std::string& operation, const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time) {
  if(fVerbosity > 0) {
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - start_time).count();
    std::cout << "[TIMING] " << operation << ": " << elapsed << "ms" << std::endl;
  }
}