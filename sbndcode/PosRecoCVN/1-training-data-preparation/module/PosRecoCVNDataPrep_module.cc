// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// ROOT and C++ includes
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <iomanip>
#include "TTree.h"
#include "art_root_io/TFileService.h"

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Simulation includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"

// Reco includes
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

// Geometry and mapping
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// Eigen for PCA calculation
#include <Eigen/Dense>


namespace opdet {

  // ============================================================================
  // Photosensor type classification (same as PosRecoCVNProducer)
  // ============================================================================
  enum SBNDPDSDetectorType_DP {
    kPDUnknown_DP = -1,
    kPMTCoated_DP = 0,
    kPMTUncoated_DP = 1,
    kXARAPUCAVUV_DP,
    kXARAPUCAVIS_DP
  };

  // ============================================================================
  // DATA STRUCTURES (mirroring PosRecoCVNProducer structures)
  // ============================================================================

  struct EventInfo_DP {
    int eventID = 0;
    int runID = 0;
    int subrunID = 0;
    void clear() { eventID = runID = subrunID = 0; }
  };

  struct MCTruthData_DP {
    std::vector<double> nuvT, nuvX, nuvY, nuvZ, nuvE;
    std::vector<int>    nuvPDG;

    // Energy-weighted position (training targets)
    std::vector<double> dEpromx, dEpromy, dEpromz;
    std::vector<double> dEtpc;
    std::vector<double> dEspreadx, dEspready, dEspreadz;
    std::vector<double> dEdirx, dEdiry, dEdirz;
    // PCA eigenvalues (λ1 ≥ λ2 ≥ λ3, energy-weighted covariance matrix, units cm²)
    // λ1 corresponds to dEdirx/y/z (already saved). λ2,λ3 are the secondary/tertiary.
    std::vector<double> dEpcaLam1, dEpcaLam2, dEpcaLam3;
    // Secondary (v2) and tertiary (v3) PCA eigenvectors (v1 = dEdirx/y/z)
    std::vector<double> dEpcaV2x, dEpcaV2y, dEpcaV2z;
    std::vector<double> dEpcaV3x, dEpcaV3y, dEpcaV3z;
    std::vector<std::vector<double>> dElowedges, dEmaxedges;

    // Geant4 particle trajectories
    std::vector<std::vector<double>> stepX, stepY, stepZ, stepT;
    std::vector<double> dE, E;
    std::vector<int> trackID, motherID, PDGcode;
    std::vector<std::string> process;
    std::vector<double> StartPx, StartPy, StartPz, EndPx, EndPy, EndPz;
    std::vector<std::vector<double>> energydep, energydepX, energydepY, energydepZ;

    int InTimeCosmics = 0;
    std::vector<double> InTimeCosmicsTime;
    double neutrinowindow = 0.0;
    std::vector<simb::MCParticle> mcpartVec;

    void clear() {
      nuvT.clear(); nuvX.clear(); nuvY.clear(); nuvZ.clear(); nuvE.clear(); nuvPDG.clear();
      dEpromx.clear(); dEpromy.clear(); dEpromz.clear(); dEtpc.clear();
      dEspreadx.clear(); dEspready.clear(); dEspreadz.clear();
      dEdirx.clear(); dEdiry.clear(); dEdirz.clear();
      dEpcaLam1.clear(); dEpcaLam2.clear(); dEpcaLam3.clear();
      dEpcaV2x.clear(); dEpcaV2y.clear(); dEpcaV2z.clear();
      dEpcaV3x.clear(); dEpcaV3y.clear(); dEpcaV3z.clear();
      dElowedges.clear(); dEmaxedges.clear();
      stepX.clear(); stepY.clear(); stepZ.clear(); stepT.clear();
      dE.clear(); E.clear();
      trackID.clear(); motherID.clear(); PDGcode.clear(); process.clear();
      StartPx.clear(); StartPy.clear(); StartPz.clear();
      EndPx.clear(); EndPy.clear(); EndPz.clear();
      energydep.clear(); energydepX.clear(); energydepY.clear(); energydepZ.clear();
      InTimeCosmics = 0; InTimeCosmicsTime.clear(); neutrinowindow = 0.0; mcpartVec.clear();
    }
  };

  struct OpticalData_DP {
    std::vector<double> flash_x, flash_y, flash_z;
    std::vector<std::vector<float>> flash_ophit_pe;
    std::vector<std::vector<int>>   flash_ophit_ch;
    std::vector<std::vector<float>> flash_ophit_time;

    void clear() {
      flash_x.clear(); flash_y.clear(); flash_z.clear();
      flash_ophit_pe.clear(); flash_ophit_ch.clear(); flash_ophit_time.clear();
    }
  };

  struct ProcessedData_DP {
    std::vector<std::vector<float>> flash_ophit_pe_sel;
    std::vector<std::vector<int>>   flash_ophit_ch_sel;
    std::vector<std::vector<float>> flash_ophit_time_sel;
    std::vector<int> categorized_flashes;
    std::vector<double> dEpromx_sel, dEpromy_sel, dEpromz_sel, dEtpc_sel;
    int selectedTPC = -1;

    std::vector<std::vector<float>> flash_ophit_pe_final;
    std::vector<std::vector<int>>   flash_ophit_ch_final;
    std::vector<std::vector<float>> flash_ophit_time_final;
    std::vector<double> nuvT_final, nuvZ_final;
    std::vector<double> dEpromx_final, dEpromy_final, dEpromz_final, dEtpc_final;

    std::vector<std::vector<float>> pe_matrix;
    // pe_images[event][y][z][map]: 0=Uncoated, 1=Coated  — raw (not normalized)
    std::vector<std::vector<std::vector<std::vector<float>>>> pe_images;

    void clear() {
      flash_ophit_pe_sel.clear(); flash_ophit_ch_sel.clear(); flash_ophit_time_sel.clear();
      categorized_flashes.clear();
      dEpromx_sel.clear(); dEpromy_sel.clear(); dEpromz_sel.clear(); dEtpc_sel.clear();
      selectedTPC = -1;
      flash_ophit_pe_final.clear(); flash_ophit_ch_final.clear(); flash_ophit_time_final.clear();
      nuvT_final.clear(); nuvZ_final.clear();
      dEpromx_final.clear(); dEpromy_final.clear(); dEpromz_final.clear(); dEtpc_final.clear();
      pe_matrix.clear(); pe_images.clear();
    }
  };

} // namespace opdet


// ============================================================================
// PosRecoCVNDataPrep
//
// art EDAnalyzer that prepares MC data for CNN training.
//
// Mirrors the PosRecoCVNProducer processing chain up to and including
// CreatePEImages(), but saves raw (un-normalized) images and MC truth labels
// to a TTree instead of running TensorFlow inference.
//
// Eliminates the need for the preprocess_images.py step:
//   TTree output contains images + labels ready for Python training.
//
// Output TTree branches:
//   - Event ID, filter flags
//   - MC truth: nuv vertex, energy
//   - Training labels: dEpromx/y/z, dEtpc, dEdirx/y/z, dEspreadx/y/z (per TPC)
//   - pe_per_channel: flat [312] PE vector (CNN input before spatial mapping)
//   - flash_ophit_pe/ch/time: per-flash optical hit data
//   - image_uncoated / image_coated: raw 2D images (flattened), + ny, nz, max_pe
// ============================================================================

namespace opdet { class PosRecoCVNDataPrep; }

class opdet::PosRecoCVNDataPrep : public art::EDAnalyzer {
public:
  explicit PosRecoCVNDataPrep(fhicl::ParameterSet const& p);

  PosRecoCVNDataPrep(PosRecoCVNDataPrep const&) = delete;
  PosRecoCVNDataPrep(PosRecoCVNDataPrep&&) = delete;
  PosRecoCVNDataPrep& operator=(PosRecoCVNDataPrep const&) = delete;
  PosRecoCVNDataPrep& operator=(PosRecoCVNDataPrep&&) = delete;

  void analyze(art::Event const& e) override;
  void beginJob() override;

private:

  // ─── Core processing (copied from PosRecoCVNProducer) ───────────────────
  void FillMCTruth(art::Event const& e);
  void InitializeChannelDict();
  void ClassifyChannels();
  int  CategorizeFirstChannel(const std::vector<int>& channels);
  void ApplyFlashSelection();
  void ApplyFinalEnergyFilter();
  void CreatePEMatrix();
  void LoadPMTMaps();
  void ClearEventData();
  std::vector<std::vector<float>> SelectNonEmptyHalf(
    const std::vector<std::vector<float>>& left_half,
    const std::vector<std::vector<float>>& right_half,
    const std::string& method = "max");
  void FillAverageDepositedEnergyVariables(
    std::vector<std::vector<double>> fenergydep,
    std::vector<std::vector<double>> fenergydepX,
    std::vector<std::vector<double>> fenergydepY,
    std::vector<std::vector<double>> fenergydepZ,
    std::vector<std::vector<double>> fstepT,
    std::vector<double>& dEtpc,
    std::vector<double>& dEpromx, std::vector<double>& dEpromy, std::vector<double>& dEpromz,
    std::vector<double>& dEspreadx, std::vector<double>& dEspready, std::vector<double>& dEspreadz,
    std::vector<std::vector<double>>& dElowedges, std::vector<std::vector<double>>& dEmaxedges,
    std::vector<double>& dEdirx, std::vector<double>& dEdiry, std::vector<double>& dEdirz,
    std::vector<double>& dEpcaLam1, std::vector<double>& dEpcaLam2, std::vector<double>& dEpcaLam3,
    std::vector<double>& dEpcaV2x, std::vector<double>& dEpcaV2y, std::vector<double>& dEpcaV2z,
    std::vector<double>& dEpcaV3x, std::vector<double>& dEpcaV3y, std::vector<double>& dEpcaV3z);

  template<typename T>
  std::vector<std::vector<T>> FilterByMask(
    const std::vector<std::vector<T>>& array,
    const std::vector<std::vector<bool>>& mask);

  void LogTiming(const std::string& op,
                 const std::chrono::time_point<std::chrono::high_resolution_clock>& t0);

  // ─── DataPrep-specific ──────────────────────────────────────────────────
  // Creates raw (un-normalized) 2D images; computes per-map max PE
  void CreateRawPEImages();

  // Fills the TTree from the processed event data
  void FillTrainingTree(bool passedFilters);

  // ─── FCL parameters ─────────────────────────────────────────────────────
  std::vector<int>         fMCTruthOrigin;
  std::vector<int>         fMCTruthPDG;
  std::vector<std::string> fMCTruthModuleLabel;
  std::vector<std::string> fMCTruthInstanceLabel;
  std::string              fMCModuleLabel;
  std::vector<std::string> fOpHitsModuleLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  std::vector<int>         fG4BufferBoxX;
  std::vector<int>         fG4BufferBoxY;
  std::vector<int>         fG4BufferBoxZ;
  std::vector<int>         fG4BeamWindow;
  std::vector<int>         fKeepPDGCode;
  std::vector<double>      fActiveVolumeX;
  std::vector<double>      fActiveVolumeY;
  std::vector<double>      fActiveVolumeZ;
  bool                     fSaveOpHits;
  int                      fVerbosity;
  std::string              fSbndcodeVersion;   // for PMT map path search
  std::string              fCoatedPMTMapPath;
  std::string              fUncoatedPMTMapPath;

  // Always MC_testing (no FCL parameter needed)
  static constexpr const char* kProcessingMode = "MC_testing";
  static constexpr double fDefaultSimIDE = -999.0;

  // ─── Event-level data ───────────────────────────────────────────────────
  EventInfo_DP    fEventInfo;
  MCTruthData_DP  fMCData;
  OpticalData_DP  fOpticalData;
  ProcessedData_DP fProcessedData;

  // ─── Channel classification ──────────────────────────────────────────────
  std::map<int, int> fChannelDict;
  std::set<int> fPMTEven, fPMTOdd;
  std::set<int> fXASEven, fXASOdd;

  // ─── PMT spatial maps ────────────────────────────────────────────────────
  std::vector<std::vector<int>> fCoatedPMTMap;
  std::vector<std::vector<int>> fUncoatedPMTMap;

  // ─── TTree and its branch variables ─────────────────────────────────────
  TTree* fTrainingTree;

  // Event ID
  int  fRun, fSubrun, fEvent;
  bool fPassedFilters;
  int  fSelectedTPC;

  // MC truth (per neutrino interaction)
  std::vector<double> fNuvX, fNuvY, fNuvZ, fNuvT, fNuvE;
  std::vector<int>    fNuvPDG;

  // Training labels: energy-weighted centroid, energy, direction, spread (2 entries: TPC0/TPC1)
  std::vector<double> fDEpromX, fDEpromY, fDEpromZ;
  std::vector<double> fDEtpc;
  std::vector<double> fDEdirX, fDEdirY, fDEdirZ;
  std::vector<double> fDEspreadX, fDEspreadY, fDEspreadZ;
  // PCA eigenvalues and secondary/tertiary eigenvectors (2 entries: TPC0/TPC1)
  std::vector<double> fDEpcaLam1, fDEpcaLam2, fDEpcaLam3;
  std::vector<double> fDEpcaV2x, fDEpcaV2y, fDEpcaV2z;
  std::vector<double> fDEpcaV3x, fDEpcaV3y, fDEpcaV3z;

  // PE flat per channel: [312] summed over all selected flashes
  std::vector<float> fPePerChannel;

  // Per-flash optical hit data (final selection)
  std::vector<std::vector<float>> fFlashOphitPE;
  std::vector<std::vector<int>>   fFlashOphitCh;
  std::vector<std::vector<float>> fFlashOphitTime;

  // 2D CNN images (raw, un-normalized, flattened row-major [ny * nz])
  std::vector<float> fImageUncoated;
  std::vector<float> fImageCoated;
  int   fImageNy, fImageNz;
  float fMaxPeUncoated, fMaxPeCoated;
};


// ============================================================================
// Constructor
// ============================================================================
opdet::PosRecoCVNDataPrep::PosRecoCVNDataPrep(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fMCTruthOrigin       ( p.get<std::vector<int>>("MCTruthOrigin") ),
    fMCTruthPDG          ( p.get<std::vector<int>>("MCTruthPDG") ),
    fMCTruthModuleLabel  ( p.get<std::vector<std::string>>("MCTruthModuleLabel") ),
    fMCTruthInstanceLabel( p.get<std::vector<std::string>>("MCTruthInstanceLabel") ),
    fMCModuleLabel       ( p.get<std::string>("MCModuleLabel") ),
    fOpHitsModuleLabel   ( p.get<std::vector<std::string>>("OpHitsModuleLabel") ),
    fOpFlashesModuleLabel( p.get<std::vector<std::string>>("OpFlashesModuleLabel") ),
    fG4BufferBoxX        ( p.get<std::vector<int>>("G4BufferBoxX") ),
    fG4BufferBoxY        ( p.get<std::vector<int>>("G4BufferBoxY") ),
    fG4BufferBoxZ        ( p.get<std::vector<int>>("G4BufferBoxZ") ),
    fG4BeamWindow        ( p.get<std::vector<int>>("G4BeamWindow") ),
    fKeepPDGCode         ( p.get<std::vector<int>>("KeepPDGCode", {}) ),
    fActiveVolumeX       ( p.get<std::vector<double>>("ActiveVolumeX", {-200.0, 200.0}) ),
    fActiveVolumeY       ( p.get<std::vector<double>>("ActiveVolumeY", {-200.0, 200.0}) ),
    fActiveVolumeZ       ( p.get<std::vector<double>>("ActiveVolumeZ", {0.0, 500.0}) ),
    fSaveOpHits          ( p.get<bool>("SaveOpHits", true) ),
    fVerbosity           ( p.get<int>("Verbosity", 0) ),
    fSbndcodeVersion     ( p.get<std::string>("SbndcodeVersion", "v10_09_00") ),
    fCoatedPMTMapPath    ("coatedPMT_map.csv"),
    fUncoatedPMTMapPath  ("uncoatedPMT_map.csv")
{
  InitializeChannelDict();
  ClassifyChannels();
  LoadPMTMaps();
}


// ============================================================================
// beginJob — create TTree with all training data branches
// ============================================================================
void opdet::PosRecoCVNDataPrep::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTrainingTree = tfs->make<TTree>("training_tree", "PosRecoCVN Training Data");

  // ─── Event ID ─────────────────────────────────────────────────────────────
  fTrainingTree->Branch("run",             &fRun,            "run/I");
  fTrainingTree->Branch("subrun",          &fSubrun,         "subrun/I");
  fTrainingTree->Branch("event",           &fEvent,          "event/I");
  fTrainingTree->Branch("passed_filters",  &fPassedFilters,  "passed_filters/O");
  fTrainingTree->Branch("selected_tpc",    &fSelectedTPC,    "selected_tpc/I");

  // ─── MC truth (one entry per neutrino interaction) ────────────────────────
  fTrainingTree->Branch("nuv_x",   &fNuvX);
  fTrainingTree->Branch("nuv_y",   &fNuvY);
  fTrainingTree->Branch("nuv_z",   &fNuvZ);
  fTrainingTree->Branch("nuv_t",   &fNuvT);
  fTrainingTree->Branch("nuv_e",   &fNuvE);
  fTrainingTree->Branch("nuv_pdg", &fNuvPDG);

  // ─── Training labels (2 entries: TPC0, TPC1) ──────────────────────────────
  fTrainingTree->Branch("dEprom_x",    &fDEpromX);
  fTrainingTree->Branch("dEprom_y",    &fDEpromY);
  fTrainingTree->Branch("dEprom_z",    &fDEpromZ);
  fTrainingTree->Branch("dEtpc",       &fDEtpc);
  fTrainingTree->Branch("dEdir_x",     &fDEdirX);
  fTrainingTree->Branch("dEdir_y",     &fDEdirY);
  fTrainingTree->Branch("dEdir_z",     &fDEdirZ);
  fTrainingTree->Branch("dEspread_x",  &fDEspreadX);
  fTrainingTree->Branch("dEspread_y",  &fDEspreadY);
  fTrainingTree->Branch("dEspread_z",  &fDEspreadZ);
  // PCA eigenvalues (λ1≥λ2≥λ3, cm², energy-weighted covariance; λ1 eigenvector = dEdir_x/y/z)
  fTrainingTree->Branch("dEpca_lam1",  &fDEpcaLam1);
  fTrainingTree->Branch("dEpca_lam2",  &fDEpcaLam2);
  fTrainingTree->Branch("dEpca_lam3",  &fDEpcaLam3);
  // Secondary (v2) and tertiary (v3) PCA eigenvectors
  fTrainingTree->Branch("dEpca_v2x",   &fDEpcaV2x);
  fTrainingTree->Branch("dEpca_v2y",   &fDEpcaV2y);
  fTrainingTree->Branch("dEpca_v2z",   &fDEpcaV2z);
  fTrainingTree->Branch("dEpca_v3x",   &fDEpcaV3x);
  fTrainingTree->Branch("dEpca_v3y",   &fDEpcaV3y);
  fTrainingTree->Branch("dEpca_v3z",   &fDEpcaV3z);

  // ─── PE flat per channel (before image mapping) ────────────────────────────
  fTrainingTree->Branch("pe_per_channel", &fPePerChannel);

  // ─── Per-flash selected optical hit data ──────────────────────────────────
  fTrainingTree->Branch("flash_ophit_pe",   &fFlashOphitPE);
  fTrainingTree->Branch("flash_ophit_ch",   &fFlashOphitCh);
  fTrainingTree->Branch("flash_ophit_time", &fFlashOphitTime);

  // ─── 2D images (raw, un-normalized, row-major flat: ny * nz elements) ──────
  fTrainingTree->Branch("image_uncoated",   &fImageUncoated);
  fTrainingTree->Branch("image_coated",     &fImageCoated);
  fTrainingTree->Branch("image_ny",         &fImageNy,         "image_ny/I");
  fTrainingTree->Branch("image_nz",         &fImageNz,         "image_nz/I");
  fTrainingTree->Branch("max_pe_uncoated",  &fMaxPeUncoated,   "max_pe_uncoated/F");
  fTrainingTree->Branch("max_pe_coated",    &fMaxPeCoated,     "max_pe_coated/F");

  if(fVerbosity > 0)
    std::cout << "[DataPrep] TTree initialized: " << fTrainingTree->GetNbranches() << " branches" << std::endl;
}


// ============================================================================
// analyze — main per-event processing
// ============================================================================
void opdet::PosRecoCVNDataPrep::analyze(art::Event const& e)
{
  auto event_start = std::chrono::high_resolution_clock::now();
  auto section_start = event_start;

  ClearEventData();

  fEventInfo.eventID  = e.id().event();
  fEventInfo.runID    = e.id().run();
  fEventInfo.subrunID = e.id().subRun();

  if(fVerbosity > 1)
    std::cout << "[DataPrep] Run " << fEventInfo.runID
              << " | Subrun " << fEventInfo.subrunID
              << " | Event "  << fEventInfo.eventID << std::endl;

  // ─── Services ────────────────────────────────────────────────────────────
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService>       bt_serv;

  // ─── MC Truth ─────────────────────────────────────────────────────────────
  FillMCTruth(e);
  LogTiming("FillMCTruth", section_start);

  // ─── MCParticles: trajectories + SimIDEs ──────────────────────────────────
  fMCData.mcpartVec.clear();
  art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
  e.getByLabel(fMCModuleLabel, mcParticleHandle);
  if(mcParticleHandle.isValid())
    fMCData.mcpartVec = *mcParticleHandle;
  else if(fVerbosity > 0)
    std::cout << "[DataPrep] WARNING: MCParticles not found (" << fMCModuleLabel << ")" << std::endl;

  fMCData.stepX.clear(); fMCData.stepY.clear(); fMCData.stepZ.clear(); fMCData.stepT.clear();
  fMCData.dE.clear(); fMCData.E.clear();
  fMCData.trackID.clear(); fMCData.motherID.clear(); fMCData.PDGcode.clear(); fMCData.process.clear();
  fMCData.StartPx.clear(); fMCData.StartPy.clear(); fMCData.StartPz.clear();
  fMCData.EndPx.clear();   fMCData.EndPy.clear();   fMCData.EndPz.clear();
  fMCData.energydep.clear(); fMCData.energydepX.clear(); fMCData.energydepY.clear(); fMCData.energydepZ.clear();
  fMCData.InTimeCosmics = 0; fMCData.InTimeCosmicsTime.clear(); fMCData.neutrinowindow = 0.0;

  for(const auto& pPart : fMCData.mcpartVec) {
    if(!fKeepPDGCode.empty() && std::find(fKeepPDGCode.begin(), fKeepPDGCode.end(), pPart.PdgCode()) == fKeepPDGCode.end())
      continue;

    std::vector<double> xp, yp, zp, tp;
    xp.push_back(pPart.Position(0).X());
    yp.push_back(pPart.Position(0).Y());
    zp.push_back(pPart.Position(0).Z());
    tp.push_back(pPart.Position(0).T());

    for(size_t i_s = 1; i_s < pPart.NumberTrajectoryPoints(); ++i_s) {
      double t = pPart.Position(i_s).T();
      double x = pPart.Position(i_s).X();
      double y = pPart.Position(i_s).Y();
      double z = pPart.Position(i_s).Z();
      if(x < fG4BufferBoxX.at(0) || x > fG4BufferBoxX.at(1) ||
         y < fG4BufferBoxY.at(0) || y > fG4BufferBoxY.at(1) ||
         z < fG4BufferBoxZ.at(0) || z > fG4BufferBoxZ.at(1)) continue;
      if(t < fG4BeamWindow.at(0) || t > fG4BeamWindow.at(1))  continue;
      xp.push_back(x); yp.push_back(y); zp.push_back(z); tp.push_back(t);
    }

    fMCData.stepX.push_back(xp); fMCData.stepY.push_back(yp);
    fMCData.stepZ.push_back(zp); fMCData.stepT.push_back(tp);
    fMCData.E.push_back(pPart.E());
    fMCData.process.push_back(pPart.Process());
    fMCData.trackID.push_back(pPart.TrackId());
    fMCData.motherID.push_back(pPart.Mother());
    fMCData.PDGcode.push_back(pPart.PdgCode());
    fMCData.StartPx.push_back(pPart.Px(0)); fMCData.StartPy.push_back(pPart.Py(0)); fMCData.StartPz.push_back(pPart.Pz(0));
    fMCData.EndPx.push_back(pPart.EndPx()); fMCData.EndPy.push_back(pPart.EndPy()); fMCData.EndPz.push_back(pPart.EndPz());

    double endep = 0;
    std::vector<double> truedE, truedEX, truedEY, truedEZ;
    std::vector<const sim::IDE*> ides_v = bt_serv->TrackIdToSimIDEs_Ps(pPart.TrackId());
    for(auto* ide : ides_v) {
      endep += ide->energy / 3.0;
      if(pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1)) {
        truedE.push_back(ide->energy / 3.0);
        truedEX.push_back(ide->x); truedEY.push_back(ide->y); truedEZ.push_back(ide->z);
        fMCData.neutrinowindow += ide->energy / 3.0;
      }
    }
    fMCData.dE.push_back(endep);
    fMCData.energydep.push_back(truedE);
    fMCData.energydepX.push_back(truedEX); fMCData.energydepY.push_back(truedEY); fMCData.energydepZ.push_back(truedEZ);

    art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(pPart.TrackId());
    if(pPart.Position(0).T() > fG4BeamWindow.at(0) && pPart.Position(0).T() < fG4BeamWindow.at(1)
       && truth->Origin() == 2 && endep > 1.0) {
      fMCData.InTimeCosmics++;
      fMCData.InTimeCosmicsTime.push_back(pPart.Position(0).T());
    }
  }

  FillAverageDepositedEnergyVariables(
    fMCData.energydep, fMCData.energydepX, fMCData.energydepY, fMCData.energydepZ, fMCData.stepT,
    fMCData.dEtpc, fMCData.dEpromx, fMCData.dEpromy, fMCData.dEpromz,
    fMCData.dEspreadx, fMCData.dEspready, fMCData.dEspreadz,
    fMCData.dElowedges, fMCData.dEmaxedges,
    fMCData.dEdirx, fMCData.dEdiry, fMCData.dEdirz,
    fMCData.dEpcaLam1, fMCData.dEpcaLam2, fMCData.dEpcaLam3,
    fMCData.dEpcaV2x, fMCData.dEpcaV2y, fMCData.dEpcaV2z,
    fMCData.dEpcaV3x, fMCData.dEpcaV3y, fMCData.dEpcaV3z);

  LogTiming("ProcessMCParticles", section_start);

  // ─── OpHits ───────────────────────────────────────────────────────────────
  art::Handle<std::vector<recob::OpHit>> ophitListHandle;
  std::vector<art::Ptr<recob::OpHit>>    ophitlist;

  for(size_t s = 0; s < fOpHitsModuleLabel.size(); ++s) {
    e.getByLabel(fOpHitsModuleLabel[s], ophitListHandle);
    if(!ophitListHandle.isValid()) {
      std::string msg = "[DataPrep] OpHit not found: " + fOpHitsModuleLabel[s];
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    art::fill_ptr_vector(ophitlist, ophitListHandle);
  }
  LogTiming("ProcessOpHits", section_start);

  // ─── OpFlashes (beam window 0.367–1.9 µs for MC) ─────────────────────────
  art::Handle<std::vector<recob::OpFlash>> opflashListHandle;

  for(size_t s = 0; s < fOpFlashesModuleLabel.size(); ++s) {
    e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);
    if(!opflashListHandle.isValid()) {
      std::string msg = "[DataPrep] OpFlash not found: " + fOpFlashesModuleLabel[s];
      std::cout << msg << std::endl;
      throw std::runtime_error(msg);
    }
    art::FindManyP<recob::OpHit> flashToOpHitAssns(opflashListHandle, e, fOpFlashesModuleLabel[s]);

    for(unsigned int i = 0; i < opflashListHandle->size(); ++i) {
      art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
      recob::OpFlash Flash = *FlashPtr;

      double flash_time_us = Flash.AbsTime();
      if(flash_time_us < 0.367 || flash_time_us > 1.9) continue;  // MC beam window

      fOpticalData.flash_x.push_back(Flash.XCenter());
      fOpticalData.flash_y.push_back(Flash.YCenter());
      fOpticalData.flash_z.push_back(Flash.ZCenter());

      if(fSaveOpHits) {
        fOpticalData.flash_ophit_time.push_back({});
        fOpticalData.flash_ophit_pe.push_back({});
        fOpticalData.flash_ophit_ch.push_back({});

        for(auto ophit : flashToOpHitAssns.at(i)) {
          fOpticalData.flash_ophit_time.back().push_back(ophit->PeakTimeAbs());
          fOpticalData.flash_ophit_pe.back().push_back(ophit->PE());
          fOpticalData.flash_ophit_ch.back().push_back(ophit->OpChannel());
        }
      }
    }
  }
  LogTiming("ProcessOpFlashes", section_start);

  // ─── Event filters ────────────────────────────────────────────────────────
  // Filter 1: exactly 1 neutrino in active volume
  int nNeutrinosInActiveVolume = 0;
  for(size_t i = 0; i < fMCData.nuvT.size(); ++i) {
    bool inAV = (fMCData.nuvX[i] >= fActiveVolumeX[0] && fMCData.nuvX[i] <= fActiveVolumeX[1]) &&
                (fMCData.nuvY[i] >= fActiveVolumeY[0] && fMCData.nuvY[i] <= fActiveVolumeY[1]) &&
                (fMCData.nuvZ[i] >= fActiveVolumeZ[0] && fMCData.nuvZ[i] <= fActiveVolumeZ[1]);
    if(inAV) nNeutrinosInActiveVolume++;
  }
  bool passFilter1 = (nNeutrinosInActiveVolume == 1);

  // Filter 2: at least one flash with OpHits
  bool passFilter2 = false;
  for(const auto& fpe : fOpticalData.flash_ophit_pe)
    if(!fpe.empty()) { passFilter2 = true; break; }

  bool passFilter = passFilter1 && passFilter2;

  if(fVerbosity > 0 && !passFilter1)
    std::cout << "[DataPrep] FILTER 1 FAIL: " << nNeutrinosInActiveVolume << " neutrinos in AV (need 1)" << std::endl;
  if(fVerbosity > 0 && !passFilter2)
    std::cout << "[DataPrep] FILTER 2 FAIL: No flashes with OpHits" << std::endl;

  LogTiming("ApplyFilters", section_start);

  // ─── Flash selection, PE matrix, images ───────────────────────────────────
  if(passFilter) {
    ApplyFlashSelection();
    ApplyFinalEnergyFilter();
    passFilter = !fProcessedData.flash_ophit_pe_final.empty();

    if(passFilter) {
      CreatePEMatrix();
      CreateRawPEImages();
    }
  }

  LogTiming("FlashSelection+Matrix+Images", section_start);

  // ─── Fill TTree ───────────────────────────────────────────────────────────
  FillTrainingTree(passFilter);

  if(fVerbosity > 0) {
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - event_start).count();
    std::cout << "[DataPrep][TIMING] TOTAL: " << total_ms << "ms"
              << " | R" << fEventInfo.runID << " E" << fEventInfo.eventID
              << " | filter=" << (passFilter ? "PASS" : "FAIL") << std::endl;
  }
}


// ============================================================================
// CreateRawPEImages
//   - Maps 312-channel PE matrix to spatial grids [ch_y, ch_z]
//   - Selects active TPC half (top/bottom) based on max signal
//   - Saves RAW (un-normalized) images → normalization in Python
//   - Records per-map max PE separately (for later normalization)
// ============================================================================
void opdet::PosRecoCVNDataPrep::CreateRawPEImages()
{
  fProcessedData.pe_images.clear();
  fMaxPeUncoated = 0.0f;
  fMaxPeCoated   = 0.0f;
  fImageNy = 0;
  fImageNz = 0;

  if(fProcessedData.pe_matrix.empty() || fCoatedPMTMap.empty() || fUncoatedPMTMap.empty()) {
    if(fVerbosity > 0) std::cout << "[DataPrep] WARNING: Missing PE matrix or PMT maps" << std::endl;
    return;
  }

  int n_flashes = fProcessedData.pe_matrix.size();
  int ch_y      = fCoatedPMTMap.size();
  int ch_z      = fCoatedPMTMap[0].size();
  int map_count = 2;  // [0]=Uncoated, [1]=Coated

  // Build raw spatial grids: [map][event=0][y][z]
  std::vector<std::vector<std::vector<std::vector<float>>>> pe_maps(
    map_count, std::vector<std::vector<std::vector<float>>>(
      1, std::vector<std::vector<float>>(
        ch_y, std::vector<float>(ch_z, 0.0f))));

  std::vector<std::vector<std::vector<int>>*> maps = {&fUncoatedPMTMap, &fCoatedPMTMap};

  for(int idx = 0; idx < map_count; ++idx) {
    auto& map_ = *maps[idx];
    for(int y = 0; y < ch_y; ++y) {
      for(int z = 0; z < ch_z; ++z) {
        int channel = map_[y][z];
        if(channel >= 0 && channel < 312) {
          float total_pe = 0.0f;
          for(int flash = 0; flash < n_flashes; ++flash)
            total_pe += fProcessedData.pe_matrix[flash][channel];
          pe_maps[idx][0][y][z] = total_pe;
        }
      }
    }
  }

  // Per-map max PE (separately, matching preprocess_images.py: max_pe_values=[uncoated, coated])
  for(int y = 0; y < ch_y; ++y) {
    for(int z = 0; z < ch_z; ++z) {
      fMaxPeUncoated = std::max(fMaxPeUncoated, pe_maps[0][0][y][z]);
      fMaxPeCoated   = std::max(fMaxPeCoated,   pe_maps[1][0][y][z]);
    }
  }

  // Select active TPC half and build output images: [1][ch_y/2][ch_z][map_count] (RAW)
  fImageNy = ch_y / 2;
  fImageNz = ch_z;

  fProcessedData.pe_images.resize(1, std::vector<std::vector<std::vector<float>>>(
    fImageNy, std::vector<std::vector<float>>(
      fImageNz, std::vector<float>(map_count, 0.0f))));

  for(int idx = 0; idx < map_count; ++idx) {
    const auto& ev = pe_maps[idx][0];
    std::vector<std::vector<float>> top   (ev.begin(),               ev.begin() + fImageNy);
    std::vector<std::vector<float>> bottom(ev.begin() + fImageNy,    ev.end());
    auto selected = SelectNonEmptyHalf(top, bottom, "max");

    for(int y = 0; y < fImageNy; ++y)
      for(int z = 0; z < fImageNz; ++z)
        fProcessedData.pe_images[0][y][z][idx] = selected[y][z];  // RAW
  }

  if(fVerbosity > 1)
    std::cout << "[DataPrep] Raw images: [" << fImageNy << "×" << fImageNz << "×2]"
              << "  max_pe_uncoated=" << fMaxPeUncoated
              << "  max_pe_coated="   << fMaxPeCoated << std::endl;
}


// ============================================================================
// FillTrainingTree
// ============================================================================
void opdet::PosRecoCVNDataPrep::FillTrainingTree(bool passedFilters)
{
  fRun           = fEventInfo.runID;
  fSubrun        = fEventInfo.subrunID;
  fEvent         = fEventInfo.eventID;
  fPassedFilters = passedFilters;
  fSelectedTPC   = fProcessedData.selectedTPC;

  // MC truth (always filled if FillMCTruth succeeded)
  fNuvX   = fMCData.nuvX;
  fNuvY   = fMCData.nuvY;
  fNuvZ   = fMCData.nuvZ;
  fNuvT   = fMCData.nuvT;
  fNuvE   = fMCData.nuvE;
  fNuvPDG = fMCData.nuvPDG;

  // Training labels (always filled after FillAverageDepositedEnergyVariables)
  fDEpromX   = fMCData.dEpromx;
  fDEpromY   = fMCData.dEpromy;
  fDEpromZ   = fMCData.dEpromz;
  fDEtpc     = fMCData.dEtpc;
  fDEdirX    = fMCData.dEdirx;
  fDEdirY    = fMCData.dEdiry;
  fDEdirZ    = fMCData.dEdirz;
  fDEspreadX = fMCData.dEspreadx;
  fDEspreadY = fMCData.dEspready;
  fDEspreadZ = fMCData.dEspreadz;
  fDEpcaLam1 = fMCData.dEpcaLam1;
  fDEpcaLam2 = fMCData.dEpcaLam2;
  fDEpcaLam3 = fMCData.dEpcaLam3;
  fDEpcaV2x  = fMCData.dEpcaV2x;
  fDEpcaV2y  = fMCData.dEpcaV2y;
  fDEpcaV2z  = fMCData.dEpcaV2z;
  fDEpcaV3x  = fMCData.dEpcaV3x;
  fDEpcaV3y  = fMCData.dEpcaV3y;
  fDEpcaV3z  = fMCData.dEpcaV3z;

  // PE flat per channel + flash ophit data + images: only meaningful if passedFilters
  fPePerChannel.assign(312, 0.0f);
  fFlashOphitPE.clear();
  fFlashOphitCh.clear();
  fFlashOphitTime.clear();
  fImageUncoated.clear();
  fImageCoated.clear();

  if(passedFilters) {
    // pe_per_channel: flat [312] vector, already summed across all selected flashes by CreatePEMatrix
    if(!fProcessedData.pe_matrix.empty())
      fPePerChannel = fProcessedData.pe_matrix[0];

    // Per-flash selected ophit data
    fFlashOphitPE   = fProcessedData.flash_ophit_pe_final;
    fFlashOphitCh   = fProcessedData.flash_ophit_ch_final;
    fFlashOphitTime = fProcessedData.flash_ophit_time_final;

    // Flatten 2D images: row-major [ny * nz]
    if(!fProcessedData.pe_images.empty()) {
      fImageUncoated.reserve(fImageNy * fImageNz);
      fImageCoated.reserve(fImageNy * fImageNz);
      for(int y = 0; y < fImageNy; ++y)
        for(int z = 0; z < fImageNz; ++z) {
          fImageUncoated.push_back(fProcessedData.pe_images[0][y][z][0]);
          fImageCoated.push_back(fProcessedData.pe_images[0][y][z][1]);
        }
    }
  }

  fTrainingTree->Fill();

  if(fVerbosity > 1)
    std::cout << "[DataPrep] TTree filled for R" << fRun << " E" << fEvent
              << " | passed=" << passedFilters << std::endl;
}


// ============================================================================
// FillMCTruth — Extract neutrino vertex information
//   (extended to save nuvPDG)
// ============================================================================
void opdet::PosRecoCVNDataPrep::FillMCTruth(art::Event const& e)
{
  if(fMCTruthModuleLabel.size() != fMCTruthInstanceLabel.size()) {
    throw std::invalid_argument("[DataPrep] MCTruthModuleLabel and MCTruthInstanceLabel size mismatch");
  }

  art::Handle<std::vector<simb::MCTruth>> MCTruthListHandle;

  for(size_t s = 0; s < fMCTruthModuleLabel.size(); ++s) {
    e.getByLabel(fMCTruthModuleLabel[s], fMCTruthInstanceLabel[s], MCTruthListHandle);

    if(!MCTruthListHandle.isValid() || MCTruthListHandle->empty()) {
      throw std::runtime_error("[DataPrep] MCTruth not found: " + fMCTruthModuleLabel[s]);
    }

    std::vector<art::Ptr<simb::MCTruth>> mctruth_v;
    art::fill_ptr_vector(mctruth_v, MCTruthListHandle);

    for(size_t n = 0; n < mctruth_v.size(); ++n) {
      art::Ptr<simb::MCTruth> evtTruth = mctruth_v[n];
      if(std::find(fMCTruthOrigin.begin(), fMCTruthOrigin.end(), (int)evtTruth->Origin()) == fMCTruthOrigin.end())
        continue;

      double nu_x = fDefaultSimIDE, nu_y = fDefaultSimIDE, nu_z = fDefaultSimIDE;
      double nu_t = fDefaultSimIDE, nu_E = fDefaultSimIDE;
      int    nu_pdg = 0;

      for(int p = 0; p < evtTruth->NParticles(); ++p) {
        simb::MCParticle const& par = evtTruth->GetParticle(p);

        if(std::find(fMCTruthPDG.begin(), fMCTruthPDG.end(), par.PdgCode()) == fMCTruthPDG.end())
          continue;

        if(par.StatusCode() == 0 && evtTruth->Origin() == 1) {
          nu_x   = par.Vx();
          nu_y   = par.Vy();
          nu_z   = par.Vz();
          nu_t   = par.T();
          nu_E   = par.E();
          nu_pdg = par.PdgCode();
        } else if(par.StatusCode() == 1 && evtTruth->Origin() == 4) {
          nu_x   = par.Vx();
          nu_y   = par.Vy();
          nu_z   = par.Vz();
          nu_t   = par.T();
          nu_E   = par.E();
          nu_pdg = par.PdgCode();
        }
      }

      fMCData.nuvT.push_back(nu_t);
      fMCData.nuvX.push_back(nu_x);
      fMCData.nuvY.push_back(nu_y);
      fMCData.nuvZ.push_back(nu_z);
      fMCData.nuvE.push_back(nu_E);
      fMCData.nuvPDG.push_back(nu_pdg);
    }
  }
}


// ============================================================================
// FillAverageDepositedEnergyVariables — energy-weighted centroids + PCA
//   (copied verbatim from PosRecoCVNProducer)
// ============================================================================
void opdet::PosRecoCVNDataPrep::FillAverageDepositedEnergyVariables(
  std::vector<std::vector<double>> fenergydep,
  std::vector<std::vector<double>> fenergydepX,
  std::vector<std::vector<double>> fenergydepY,
  std::vector<std::vector<double>> fenergydepZ,
  std::vector<std::vector<double>> fstepT,
  std::vector<double>& dEtpc,
  std::vector<double>& dEpromx, std::vector<double>& dEpromy, std::vector<double>& dEpromz,
  std::vector<double>& dEspreadx, std::vector<double>& dEspready, std::vector<double>& dEspreadz,
  std::vector<std::vector<double>>& dElowedges, std::vector<std::vector<double>>& dEmaxedges,
  std::vector<double>& dEdirx, std::vector<double>& dEdiry, std::vector<double>& dEdirz,
  std::vector<double>& dEpcaLam1, std::vector<double>& dEpcaLam2, std::vector<double>& dEpcaLam3,
  std::vector<double>& dEpcaV2x, std::vector<double>& dEpcaV2y, std::vector<double>& dEpcaV2z,
  std::vector<double>& dEpcaV3x, std::vector<double>& dEpcaV3y, std::vector<double>& dEpcaV3z)
{
  dEtpc.clear(); dEpromx.clear(); dEpromy.clear(); dEpromz.clear();
  dEtpc.resize(2, 0);
  dEpromx.resize(2, fDefaultSimIDE); dEpromy.resize(2, fDefaultSimIDE); dEpromz.resize(2, fDefaultSimIDE);
  dEspreadx.clear(); dEspready.clear(); dEspreadz.clear();
  dEspreadx.resize(2, fDefaultSimIDE); dEspready.resize(2, fDefaultSimIDE); dEspreadz.resize(2, fDefaultSimIDE);
  dEdirx.clear(); dEdiry.clear(); dEdirz.clear();
  dEdirx.resize(2, 0); dEdiry.resize(2, 0); dEdirz.resize(2, 1);
  dEpcaLam1.clear(); dEpcaLam2.clear(); dEpcaLam3.clear();
  dEpcaLam1.resize(2, fDefaultSimIDE); dEpcaLam2.resize(2, fDefaultSimIDE); dEpcaLam3.resize(2, fDefaultSimIDE);
  dEpcaV2x.clear(); dEpcaV2y.clear(); dEpcaV2z.clear();
  dEpcaV2x.resize(2, 0); dEpcaV2y.resize(2, 0); dEpcaV2z.resize(2, 0);
  dEpcaV3x.clear(); dEpcaV3y.clear(); dEpcaV3z.clear();
  dEpcaV3x.resize(2, 0); dEpcaV3y.resize(2, 0); dEpcaV3z.resize(2, 0);

  int ndeps_tpc0 = 0, ndeps_tpc1 = 0;
  double dEpromx_tpc0=0, dEpromy_tpc0=0, dEpromz_tpc0=0;
  double dEpromx_tpc1=0, dEpromy_tpc1=0, dEpromz_tpc1=0;
  double spreadx_tpc0=0, spreadx_tpc1=0;
  double spready_tpc0=0, spready_tpc1=0;
  double spreadz_tpc0=0, spreadz_tpc1=0;
  double minz_tpc0=1e3, maxz_tpc0=-1e3, minz_tpc1=1e3, maxz_tpc1=-1e3;
  double miny_tpc0=1e3, maxy_tpc0=-1e3, miny_tpc1=1e3, maxy_tpc1=-1e3;
  double minx_tpc0=1e3, maxx_tpc0=-1e3, minx_tpc1=1e3, maxx_tpc1=-1e3;
  double dE_tpc0=0, dE_tpc1=0;

  std::vector<double> x_tpc0, y_tpc0, z_tpc0, w_tpc0;
  std::vector<double> x_tpc1, y_tpc1, z_tpc1, w_tpc1;

  for(size_t k = 0; k < fenergydep.size(); ++k) {
    for(size_t j = 0; j < fenergydep.at(k).size(); ++j) {
      if(fenergydepX.at(k).at(j) < 0) {
        ndeps_tpc0++;
        dE_tpc0     += fenergydep.at(k).at(j);
        dEpromx_tpc0 += fenergydep.at(k).at(j) * fenergydepX.at(k).at(j);
        dEpromy_tpc0 += fenergydep.at(k).at(j) * fenergydepY.at(k).at(j);
        dEpromz_tpc0 += fenergydep.at(k).at(j) * fenergydepZ.at(k).at(j);
        spreadx_tpc0 += fenergydep.at(k).at(j) * fenergydepX.at(k).at(j) * fenergydepX.at(k).at(j);
        spready_tpc0 += fenergydep.at(k).at(j) * fenergydepY.at(k).at(j) * fenergydepY.at(k).at(j);
        spreadz_tpc0 += fenergydep.at(k).at(j) * fenergydepZ.at(k).at(j) * fenergydepZ.at(k).at(j);
        x_tpc0.push_back(fenergydepX.at(k).at(j)); y_tpc0.push_back(fenergydepY.at(k).at(j));
        z_tpc0.push_back(fenergydepZ.at(k).at(j)); w_tpc0.push_back(fenergydep.at(k).at(j));
        minx_tpc0 = std::min(minx_tpc0, fenergydepX.at(k).at(j)); maxx_tpc0 = std::max(maxx_tpc0, fenergydepX.at(k).at(j));
        miny_tpc0 = std::min(miny_tpc0, fenergydepY.at(k).at(j)); maxy_tpc0 = std::max(maxy_tpc0, fenergydepY.at(k).at(j));
        minz_tpc0 = std::min(minz_tpc0, fenergydepZ.at(k).at(j)); maxz_tpc0 = std::max(maxz_tpc0, fenergydepZ.at(k).at(j));
      } else {
        ndeps_tpc1++;
        dE_tpc1     += fenergydep.at(k).at(j);
        dEpromx_tpc1 += fenergydep.at(k).at(j) * fenergydepX.at(k).at(j);
        dEpromy_tpc1 += fenergydep.at(k).at(j) * fenergydepY.at(k).at(j);
        dEpromz_tpc1 += fenergydep.at(k).at(j) * fenergydepZ.at(k).at(j);
        spreadx_tpc1 += fenergydep.at(k).at(j) * fenergydepX.at(k).at(j) * fenergydepX.at(k).at(j);
        spready_tpc1 += fenergydep.at(k).at(j) * fenergydepY.at(k).at(j) * fenergydepY.at(k).at(j);
        spreadz_tpc1 += fenergydep.at(k).at(j) * fenergydepZ.at(k).at(j) * fenergydepZ.at(k).at(j);
        x_tpc1.push_back(fenergydepX.at(k).at(j)); y_tpc1.push_back(fenergydepY.at(k).at(j));
        z_tpc1.push_back(fenergydepZ.at(k).at(j)); w_tpc1.push_back(fenergydep.at(k).at(j));
        minx_tpc1 = std::min(minx_tpc1, fenergydepX.at(k).at(j)); maxx_tpc1 = std::max(maxx_tpc1, fenergydepX.at(k).at(j));
        miny_tpc1 = std::min(miny_tpc1, fenergydepY.at(k).at(j)); maxy_tpc1 = std::max(maxy_tpc1, fenergydepY.at(k).at(j));
        minz_tpc1 = std::min(minz_tpc1, fenergydepZ.at(k).at(j)); maxz_tpc1 = std::max(maxz_tpc1, fenergydepZ.at(k).at(j));
      }
    }
  }

  if(ndeps_tpc0 != 0) {
    double cx0 = dEpromx_tpc0/dE_tpc0, cy0 = dEpromy_tpc0/dE_tpc0, cz0 = dEpromz_tpc0/dE_tpc0;
    dEpromx_tpc0 = std::abs(cx0); dEpromy_tpc0 = cy0; dEpromz_tpc0 = cz0;
    spreadx_tpc0 = std::sqrt(spreadx_tpc0/dE_tpc0 - dEpromx_tpc0*dEpromx_tpc0);
    spready_tpc0 = std::sqrt(spready_tpc0/dE_tpc0 - dEpromy_tpc0*dEpromy_tpc0);
    spreadz_tpc0 = std::sqrt(spreadz_tpc0/dE_tpc0 - dEpromz_tpc0*dEpromz_tpc0);
    dEtpc[0]=dE_tpc0; dEpromx[0]=dEpromx_tpc0; dEpromy[0]=dEpromy_tpc0; dEpromz[0]=dEpromz_tpc0;
    dEspreadx[0]=spreadx_tpc0; dEspready[0]=spready_tpc0; dEspreadz[0]=spreadz_tpc0;

    if(ndeps_tpc0 >= 2) {
      Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
      for(size_t i = 0; i < x_tpc0.size(); ++i) {
        Eigen::Vector3d pt(x_tpc0[i]-cx0, y_tpc0[i]-cy0, z_tpc0[i]-cz0);
        cov += std::abs(w_tpc0[i]) * (pt * pt.transpose());
      }
      cov /= dE_tpc0;
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
      Eigen::Vector3d eigvals = solver.eigenvalues();
      Eigen::Matrix3d eigvecs = solver.eigenvectors();
      int maxIdx; eigvals.maxCoeff(&maxIdx);
      Eigen::Vector3d dir = eigvecs.col(maxIdx);
      if(dir(2) < 0) dir = -dir;
      dEdirx[0]=dir(0); dEdiry[0]=dir(1); dEdirz[0]=dir(2);
      // Eigen sorts eigenvalues ascending: col(2)=λ1 (dominant), col(1)=λ2, col(0)=λ3
      dEpcaLam1[0]=eigvals(2); dEpcaLam2[0]=eigvals(1); dEpcaLam3[0]=eigvals(0);
      Eigen::Vector3d v2=eigvecs.col(1), v3=eigvecs.col(0);
      dEpcaV2x[0]=v2(0); dEpcaV2y[0]=v2(1); dEpcaV2z[0]=v2(2);
      dEpcaV3x[0]=v3(0); dEpcaV3y[0]=v3(1); dEpcaV3z[0]=v3(2);
    }
  }

  if(ndeps_tpc1 != 0) {
    double cx1 = dEpromx_tpc1/dE_tpc1, cy1 = dEpromy_tpc1/dE_tpc1, cz1 = dEpromz_tpc1/dE_tpc1;
    dEpromx_tpc1 = std::abs(cx1); dEpromy_tpc1 = cy1; dEpromz_tpc1 = cz1;
    spreadx_tpc1 = std::sqrt(spreadx_tpc1/dE_tpc1 - dEpromx_tpc1*dEpromx_tpc1);
    spready_tpc1 = std::sqrt(spready_tpc1/dE_tpc1 - dEpromy_tpc1*dEpromy_tpc1);
    spreadz_tpc1 = std::sqrt(spreadz_tpc1/dE_tpc1 - dEpromz_tpc1*dEpromz_tpc1);
    dEtpc[1]=dE_tpc1; dEpromx[1]=dEpromx_tpc1; dEpromy[1]=dEpromy_tpc1; dEpromz[1]=dEpromz_tpc1;
    dEspreadx[1]=spreadx_tpc1; dEspready[1]=spready_tpc1; dEspreadz[1]=spreadz_tpc1;

    if(ndeps_tpc1 >= 2) {
      Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
      for(size_t i = 0; i < x_tpc1.size(); ++i) {
        Eigen::Vector3d pt(x_tpc1[i]-cx1, y_tpc1[i]-cy1, z_tpc1[i]-cz1);
        cov += std::abs(w_tpc1[i]) * (pt * pt.transpose());
      }
      cov /= dE_tpc1;
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
      Eigen::Vector3d eigvals = solver.eigenvalues();
      Eigen::Matrix3d eigvecs = solver.eigenvectors();
      int maxIdx; eigvals.maxCoeff(&maxIdx);
      Eigen::Vector3d dir = eigvecs.col(maxIdx);
      if(dir(2) < 0) dir = -dir;
      dEdirx[1]=dir(0); dEdiry[1]=dir(1); dEdirz[1]=dir(2);
      dEpcaLam1[1]=eigvals(2); dEpcaLam2[1]=eigvals(1); dEpcaLam3[1]=eigvals(0);
      Eigen::Vector3d v2=eigvecs.col(1), v3=eigvecs.col(0);
      dEpcaV2x[1]=v2(0); dEpcaV2y[1]=v2(1); dEpcaV2z[1]=v2(2);
      dEpcaV3x[1]=v3(0); dEpcaV3y[1]=v3(1); dEpcaV3z[1]=v3(2);
    }
  }

  dElowedges.clear(); dElowedges.resize(2);
  dEmaxedges.clear(); dEmaxedges.resize(2);
  dElowedges[0] = {minx_tpc0, miny_tpc0, minz_tpc0};
  dElowedges[1] = {minx_tpc1, miny_tpc1, minz_tpc1};
  dEmaxedges[0] = {maxx_tpc0, maxy_tpc0, maxz_tpc0};
  dEmaxedges[1] = {maxx_tpc1, maxy_tpc1, maxz_tpc1};
}


// ============================================================================
// InitializeChannelDict (copied from PosRecoCVNProducer)
// ============================================================================
void opdet::PosRecoCVNDataPrep::InitializeChannelDict()
{
  opdet::sbndPDMapAlg pdsMap;
  fChannelDict.clear();
  for(size_t ch = 0; ch < pdsMap.size(); ++ch) {
    std::string pdType = pdsMap.pdType(ch);
    int typeInt = static_cast<int>(kPDUnknown_DP);
    if     (pdType == "pmt_coated")    typeInt = static_cast<int>(kPMTCoated_DP);
    else if(pdType == "pmt_uncoated")  typeInt = static_cast<int>(kPMTUncoated_DP);
    else if(pdType == "xarapuca_vuv")  typeInt = static_cast<int>(kXARAPUCAVUV_DP);
    else if(pdType == "xarapuca_vis")  typeInt = static_cast<int>(kXARAPUCAVIS_DP);
    fChannelDict[static_cast<int>(ch)] = typeInt;
  }
  if(fVerbosity > 0)
    std::cout << "[DataPrep] Channel dict: " << fChannelDict.size() << " channels" << std::endl;
}


// ============================================================================
// ClassifyChannels (copied from PosRecoCVNProducer)
// ============================================================================
void opdet::PosRecoCVNDataPrep::ClassifyChannels()
{
  fPMTEven.clear(); fPMTOdd.clear();
  fXASEven.clear(); fXASOdd.clear();
  for(const auto& pair : fChannelDict) {
    int ch = pair.first, type = pair.second;
    bool isPMT = (type == static_cast<int>(kPMTCoated_DP) || type == static_cast<int>(kPMTUncoated_DP));
    bool isXAS = (type == static_cast<int>(kXARAPUCAVUV_DP) || type == static_cast<int>(kXARAPUCAVIS_DP));
    bool isEven = (ch % 2 == 0);
    if(isPMT) { if(isEven) fPMTEven.insert(ch); else fPMTOdd.insert(ch); }
    else if(isXAS) { if(isEven) fXASEven.insert(ch); else fXASOdd.insert(ch); }
  }
}


// ============================================================================
// CategorizeFirstChannel (copied from PosRecoCVNProducer)
// ============================================================================
int opdet::PosRecoCVNDataPrep::CategorizeFirstChannel(const std::vector<int>& channels)
{
  if(channels.empty()) return -1;
  int ch = channels[0];
  if(fPMTEven.count(ch)) return 0;
  if(fPMTOdd.count(ch))  return 1;
  if(fXASEven.count(ch)) return 2;
  if(fXASOdd.count(ch))  return 3;
  return -1;
}


// ============================================================================
// LoadPMTMaps (adapted from PosRecoCVNProducer — uses fSbndcodeVersion)
// ============================================================================
void opdet::PosRecoCVNDataPrep::LoadPMTMaps()
{
  auto searchAndOpen = [&](const std::string& basePath, const std::string& label) -> std::string {
    std::string resolved = basePath;
    std::ifstream test(resolved);
    if(!test.is_open()) {
      size_t pos = basePath.find_last_of("/\\");
      std::string filename = (pos != std::string::npos) ? basePath.substr(pos+1) : basePath;
      std::vector<std::string> candidates = {
        filename,
        "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename,
        "../local/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/3-inference-larsoft-module/module/" + filename,
        (getenv("MRB_INSTALL") ? std::string(getenv("MRB_INSTALL")) + "/sbndcode/" + fSbndcodeVersion + "/scripts/PosRecoCVN/pmt_maps/" + filename : ""),
        (getenv("MRB_SOURCE")  ? std::string(getenv("MRB_SOURCE"))  + "/sbndcode/sbndcode/PosRecoCVN/pmt_maps/" + filename : ""),
        "../../../PosRecoCVN/pmt_maps/" + filename,
        "../../pmt_maps/" + filename,
        "pmt_maps/" + filename,
        filename
      };
      for(const auto& c : candidates) {
        std::ifstream f(c);
        if(f.is_open()) { resolved = c; f.close(); break; }
      }
    } else test.close();
    std::ifstream f(resolved);
    if(!f.is_open()) throw std::runtime_error("[DataPrep] " + label + " map not found: " + resolved);
    return resolved;
  };

  std::string coated_path   = searchAndOpen(fCoatedPMTMapPath,   "Coated PMT");
  std::string uncoated_path = searchAndOpen(fUncoatedPMTMapPath, "Uncoated PMT");

  auto parseCSV = [](const std::string& path) {
    std::vector<std::vector<int>> result;
    std::ifstream f(path);
    std::string line;
    while(std::getline(f, line)) {
      std::vector<int> row;
      std::stringstream ss(line);
      std::string cell;
      while(std::getline(ss, cell, ',')) row.push_back(std::stoi(cell));
      result.push_back(row);
    }
    return result;
  };

  fCoatedPMTMap   = parseCSV(coated_path);
  fUncoatedPMTMap = parseCSV(uncoated_path);

  if(fVerbosity > 0)
    std::cout << "[DataPrep] PMT maps: coated[" << fCoatedPMTMap.size() << "×"
              << (fCoatedPMTMap.empty() ? 0 : (int)fCoatedPMTMap[0].size()) << "]"
              << "  uncoated[" << fUncoatedPMTMap.size() << "×"
              << (fUncoatedPMTMap.empty() ? 0 : (int)fUncoatedPMTMap[0].size()) << "]" << std::endl;
}


// ============================================================================
// SelectNonEmptyHalf (copied verbatim from PosRecoCVNProducer)
// ============================================================================
std::vector<std::vector<float>> opdet::PosRecoCVNDataPrep::SelectNonEmptyHalf(
    const std::vector<std::vector<float>>& left_half,
    const std::vector<std::vector<float>>& right_half,
    const std::string& method)
{
  float left_score = 0.0f, right_score = 0.0f;
  if(method == "max") {
    for(const auto& row : left_half)  for(float v : row) left_score  = std::max(left_score,  v);
    for(const auto& row : right_half) for(float v : row) right_score = std::max(right_score, v);
  } else if(method == "sum") {
    for(const auto& row : left_half)  for(float v : row) left_score  += v;
    for(const auto& row : right_half) for(float v : row) right_score += v;
  } else if(method == "nonzero") {
    for(const auto& row : left_half)  for(float v : row) if(v != 0.f) left_score  += 1.f;
    for(const auto& row : right_half) for(float v : row) if(v != 0.f) right_score += 1.f;
  } else if(method == "mean_top") {
    std::vector<float> lv, rv;
    for(const auto& r : left_half)  for(float v : r) lv.push_back(v);
    for(const auto& r : right_half) for(float v : r) rv.push_back(v);
    std::sort(lv.rbegin(), lv.rend()); std::sort(rv.rbegin(), rv.rend());
    int n = std::min(5, (int)lv.size());
    for(int i = 0; i < n; ++i) { left_score += lv[i]; right_score += rv[i]; }
    left_score /= n; right_score /= n;
  }
  return left_score >= right_score ? left_half : right_half;
}


// ============================================================================
// ApplyFlashSelection (adapted from PosRecoCVNProducer — always MC_testing path)
// ============================================================================
void opdet::PosRecoCVNDataPrep::ApplyFlashSelection()
{
  fProcessedData.flash_ophit_pe_sel.clear();
  fProcessedData.flash_ophit_ch_sel.clear();
  fProcessedData.flash_ophit_time_sel.clear();
  fProcessedData.categorized_flashes.clear();
  fProcessedData.dEpromx_sel.clear(); fProcessedData.dEpromy_sel.clear();
  fProcessedData.dEpromz_sel.clear(); fProcessedData.dEtpc_sel.clear();

  if(fOpticalData.flash_ophit_ch.empty()) return;

  std::vector<int>   categorized;
  std::vector<float> sum_pe;
  for(const auto& ch  : fOpticalData.flash_ophit_ch) categorized.push_back(CategorizeFirstChannel(ch));
  for(const auto& fpe : fOpticalData.flash_ophit_pe) {
    float tot = 0; for(float p : fpe) tot += p; sum_pe.push_back(tot);
  }

  float sum_even = 0, sum_odd = 0;
  std::vector<bool> mask_even, mask_odd;
  for(size_t i = 0; i < categorized.size(); ++i) {
    bool is_even = (categorized[i] == 0 || categorized[i] == 2);
    bool is_odd  = (categorized[i] == 1 || categorized[i] == 3);
    mask_even.push_back(is_even); mask_odd.push_back(is_odd);
    if(is_even) sum_even += sum_pe[i];
    if(is_odd)  sum_odd  += sum_pe[i];
  }

  bool decision = (sum_even >= sum_odd);
  fProcessedData.selectedTPC = decision ? 0 : 1;
  size_t n = categorized.size();

  std::vector<bool> selected_mask;
  if(n <= 2) selected_mask.assign(n, true);
  else       selected_mask = decision ? mask_even : mask_odd;

  std::vector<std::vector<std::vector<float>>> pe_2d   = {fOpticalData.flash_ophit_pe};
  std::vector<std::vector<std::vector<int>>>   ch_2d   = {fOpticalData.flash_ophit_ch};
  std::vector<std::vector<std::vector<float>>> time_2d = {fOpticalData.flash_ophit_time};
  std::vector<std::vector<bool>> mask_2d = {selected_mask};

  auto pe_f   = FilterByMask(pe_2d,   mask_2d);
  auto ch_f   = FilterByMask(ch_2d,   mask_2d);
  auto time_f = FilterByMask(time_2d, mask_2d);

  if(!pe_f.empty())   fProcessedData.flash_ophit_pe_sel   = pe_f[0];
  if(!ch_f.empty())   fProcessedData.flash_ophit_ch_sel   = ch_f[0];
  if(!time_f.empty()) fProcessedData.flash_ophit_time_sel = time_f[0];

  // MC_testing: select corresponding TPC ground truth
  if(!fMCData.dEpromx.empty()) {
    std::vector<std::vector<double>> sel = {{decision ? 1.0 : 0.0, decision ? 0.0 : 1.0}};
    size_t max_tpc = std::min({(size_t)2, fMCData.dEpromx.size(), fMCData.dEpromy.size(), fMCData.dEpromz.size(), fMCData.dEtpc.size()});
    for(size_t i = 0; i < max_tpc; ++i) {
      if(sel[0][i] > 0.5) {
        fProcessedData.dEpromx_sel.push_back(fMCData.dEpromx[i]);
        fProcessedData.dEpromy_sel.push_back(fMCData.dEpromy[i]);
        fProcessedData.dEpromz_sel.push_back(fMCData.dEpromz[i]);
        fProcessedData.dEtpc_sel.push_back(fMCData.dEtpc[i]);
        break;
      }
    }
  }
}


// ============================================================================
// ApplyFinalEnergyFilter (MC_testing path only — adapted from PosRecoCVNProducer)
// ============================================================================
void opdet::PosRecoCVNDataPrep::ApplyFinalEnergyFilter()
{
  fProcessedData.flash_ophit_pe_final.clear();
  fProcessedData.flash_ophit_ch_final.clear();
  fProcessedData.flash_ophit_time_final.clear();
  fProcessedData.nuvT_final.clear();
  fProcessedData.nuvZ_final.clear();
  fProcessedData.dEpromx_final.clear();
  fProcessedData.dEpromy_final.clear();
  fProcessedData.dEpromz_final.clear();
  fProcessedData.dEtpc_final.clear();

  if(fProcessedData.dEpromx_sel.empty() || fProcessedData.dEtpc_sel.empty()) {
    if(fVerbosity > 0) std::cout << "[DataPrep] FAILED: no selected energy data" << std::endl;
    return;
  }

  bool passFilter = false;
  for(size_t i = 0; i < fProcessedData.dEpromx_sel.size(); ++i) {
    bool validX    = (fProcessedData.dEpromx_sel[i] != fDefaultSimIDE);
    bool validY    = (fProcessedData.dEpromy_sel[i] != fDefaultSimIDE);
    bool validZ    = (fProcessedData.dEpromz_sel[i] != fDefaultSimIDE);
    bool energyCut = (fProcessedData.dEtpc_sel[i] > 50.0);
    bool posX = (fProcessedData.dEpromx_sel[i] >= -200.0 && fProcessedData.dEpromx_sel[i] <= 200.0);
    bool posY = (fProcessedData.dEpromy_sel[i] >= -200.0 && fProcessedData.dEpromy_sel[i] <= 200.0);
    bool posZ = (fProcessedData.dEpromz_sel[i] >=    0.0 && fProcessedData.dEpromz_sel[i] <= 500.0);

    if(validX && validY && validZ && energyCut && posX && posY && posZ) {
      passFilter = true;
      fProcessedData.dEpromx_final.push_back(fProcessedData.dEpromx_sel[i]);
      fProcessedData.dEpromy_final.push_back(fProcessedData.dEpromy_sel[i]);
      fProcessedData.dEpromz_final.push_back(fProcessedData.dEpromz_sel[i]);
      fProcessedData.dEtpc_final.push_back(fProcessedData.dEtpc_sel[i]);
      break;
    }
  }

  if(passFilter) {
    fProcessedData.flash_ophit_pe_final   = fProcessedData.flash_ophit_pe_sel;
    fProcessedData.flash_ophit_ch_final   = fProcessedData.flash_ophit_ch_sel;
    fProcessedData.flash_ophit_time_final = fProcessedData.flash_ophit_time_sel;
    fProcessedData.nuvT_final = fMCData.nuvT;
    fProcessedData.nuvZ_final = fMCData.nuvZ;
  }
}


// ============================================================================
// CreatePEMatrix (copied verbatim from PosRecoCVNProducer)
// ============================================================================
void opdet::PosRecoCVNDataPrep::CreatePEMatrix()
{
  fProcessedData.pe_matrix.clear();
  const auto* pe_data = &fOpticalData.flash_ophit_pe;
  const auto* ch_data = &fOpticalData.flash_ophit_ch;

  if(!fProcessedData.flash_ophit_pe_final.empty() && !fProcessedData.flash_ophit_ch_final.empty()) {
    pe_data = &fProcessedData.flash_ophit_pe_final;
    ch_data = &fProcessedData.flash_ophit_ch_final;
  } else if(!fProcessedData.flash_ophit_pe_sel.empty() && !fProcessedData.flash_ophit_ch_sel.empty()) {
    pe_data = &fProcessedData.flash_ophit_pe_sel;
    ch_data = &fProcessedData.flash_ophit_ch_sel;
  }

  fProcessedData.pe_matrix.resize(1, std::vector<float>(312, 0.0f));

  if(pe_data->empty() || ch_data->empty()) return;

  for(size_t flash = 0; flash < pe_data->size(); ++flash) {
    if(flash >= ch_data->size()) break;
    size_t min_size = std::min((*pe_data)[flash].size(), (*ch_data)[flash].size());
    for(size_t j = 0; j < min_size; ++j) {
      int channel = (*ch_data)[flash][j];
      if(channel >= 0 && channel < 312)
        fProcessedData.pe_matrix[0][channel] += (*pe_data)[flash][j];
    }
  }
}


// ============================================================================
// ClearEventData
// ============================================================================
void opdet::PosRecoCVNDataPrep::ClearEventData()
{
  fEventInfo.clear();
  fMCData.clear();
  fOpticalData.clear();
  fProcessedData.clear();

  fRun = fSubrun = fEvent = 0;
  fPassedFilters = false;
  fSelectedTPC   = -1;
  fNuvX.clear(); fNuvY.clear(); fNuvZ.clear(); fNuvT.clear(); fNuvE.clear(); fNuvPDG.clear();
  fDEpromX.clear(); fDEpromY.clear(); fDEpromZ.clear();
  fDEtpc.clear();
  fDEdirX.clear(); fDEdirY.clear(); fDEdirZ.clear();
  fDEspreadX.clear(); fDEspreadY.clear(); fDEspreadZ.clear();
  fDEpcaLam1.clear(); fDEpcaLam2.clear(); fDEpcaLam3.clear();
  fDEpcaV2x.clear(); fDEpcaV2y.clear(); fDEpcaV2z.clear();
  fDEpcaV3x.clear(); fDEpcaV3y.clear(); fDEpcaV3z.clear();
  fPePerChannel.clear();
  fFlashOphitPE.clear(); fFlashOphitCh.clear(); fFlashOphitTime.clear();
  fImageUncoated.clear(); fImageCoated.clear();
  fImageNy = fImageNz = 0;
  fMaxPeUncoated = fMaxPeCoated = 0.0f;
}


// ============================================================================
// FilterByMask template (copied from PosRecoCVNProducer)
// ============================================================================
template<typename T>
std::vector<std::vector<T>> opdet::PosRecoCVNDataPrep::FilterByMask(
    const std::vector<std::vector<T>>& array,
    const std::vector<std::vector<bool>>& mask)
{
  std::vector<std::vector<T>> result;
  for(size_t i = 0; i < array.size() && i < mask.size(); ++i) {
    std::vector<T> filteredEvent;
    for(size_t j = 0; j < array[i].size() && j < mask[i].size(); ++j)
      if(mask[i][j]) filteredEvent.push_back(array[i][j]);
    result.push_back(filteredEvent);
  }
  return result;
}


// ============================================================================
// LogTiming
// ============================================================================
void opdet::PosRecoCVNDataPrep::LogTiming(
    const std::string& op,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& t0)
{
  if(fVerbosity > 0) {
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << "[DataPrep][TIMING] " << op << ": " << ms << "ms" << std::endl;
  }
}


DEFINE_ART_MODULE(opdet::PosRecoCVNDataPrep)
