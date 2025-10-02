// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// ROOT and C++ includes
#include <string.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <memory>
#include "TTree.h"
#include "art_root_io/TFileService.h"

// TensorFlow includes
#include "sbndcode/PosRecoCVN/inference/tf/tf_graph.h"
#include <sstream>

// Include the data structure definition
#include "sbndcode/PosRecoCVN/inference/module/PixelMapVars.h"

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Simulation includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// Reco includes
// PDS
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
// TPC
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"

// Cosmic rejection includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/T0.h"

// Geometry and mapping
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"


#define xdet_size 1000
#define ydet_size 1000
#define zmindet_size -500
#define zmaxdet_size 1800

namespace opdet {
  class PosRecoCVNProducer;

  // Photosensor type classification for SBND photon detection system
  enum SBNDPDSDetectorType {
    kPDUnknown = -1,      // Unknown or invalid detector type
    kPMTCoated = 0,       // PMT with wavelength-shifting coating
    kPMTUncoated = 1,     // Bare PMT sensitive to VUV scintillation light (128 nm)
    kXARAPUCAVUV,         // XARAPUCA device sensitive to VUV light
    kXARAPUCAVIS          // XARAPUCA device sensitive to visible light
  };

  // ============================================================================
  // GROUPED DATA STRUCTURES
  // ============================================================================
  // Event-level data organized into logical groups for improved code clarity
  // and memory management.

  // EventInfo: Event identification metadata (run, subrun, event numbers)
  struct EventInfo {
    int eventID = 0;      // Event number within subrun
    int runID = 0;        // Run number (detector configuration period)
    int subrunID = 0;     // Subrun number within run

    void clear() {
      eventID = runID = subrunID = 0;
    }
  };

  // Monte Carlo truth - only filled in MC modes
  struct MCTruthData {
    // Neutrino vertex
    std::vector<double> nuvT, nuvX, nuvY, nuvZ, nuvE;  // [μs], [cm], [GeV]

    // Energy-weighted position (CNN targets): dEprom = Σ(E_i * pos_i) / Σ(E_i)
    std::vector<double> dEpromx, dEpromy, dEpromz;  // [cm]
    std::vector<double> dEtpc;                      // Total energy [MeV]
    std::vector<double> dEspreadx, dEspready, dEspreadz;  // RMS spread [cm]
    std::vector<std::vector<double>> dElowedges, dEmaxedges;  // Bounding box [cm]

    // Geant4 particle trajectories
    std::vector<std::vector<double>> stepX, stepY, stepZ, stepT;  // [cm], [ns]
    std::vector<double> dE, E;  // [MeV], [GeV]
    std::vector<int> trackID, motherID, PDGcode;
    std::vector<std::string> process;
    std::vector<double> StartPx, StartPy, StartPz, EndPx, EndPy, EndPz;  // [GeV/c]

    // Energy deposits per step
    std::vector<std::vector<double>> energydep, energydepX, energydepY, energydepZ;

    // Cosmic ray info
    int InTimeCosmics = 0;
    std::vector<double> InTimeCosmicsTime;  // [μs]

    // Metadata
    double neutrinowindow = 0.0;  // [μs]
    std::vector<simb::MCParticle> mcpartVec;

    void clear() {
      nuvT.clear(); nuvX.clear(); nuvY.clear(); nuvZ.clear(); nuvE.clear();
      dEpromx.clear(); dEpromy.clear(); dEpromz.clear(); dEtpc.clear();
      dEspreadx.clear(); dEspready.clear(); dEspreadz.clear();
      dElowedges.clear(); dEmaxedges.clear();
      stepX.clear(); stepY.clear(); stepZ.clear(); stepT.clear();
      dE.clear(); E.clear();
      trackID.clear(); motherID.clear(); PDGcode.clear();
      process.clear();
      StartPx.clear(); StartPy.clear(); StartPz.clear();
      EndPx.clear(); EndPy.clear(); EndPz.clear();
      energydep.clear(); energydepX.clear(); energydepY.clear(); energydepZ.clear();
      InTimeCosmics = 0;
      InTimeCosmicsTime.clear();
      neutrinowindow = 0.0;
      mcpartVec.clear();
    }
  };

  // Optical detector data - OpFlash positions and OpHit signals
  struct OpticalData {
    std::vector<double> flash_x, flash_y, flash_z;  // OpFlash position [cm]

    // OpHits composing each flash: [flash][hit]
    std::vector<std::vector<float>> flash_ophit_pe;    // PE
    std::vector<std::vector<int>> flash_ophit_ch;      // Channel ID
    std::vector<std::vector<float>> flash_ophit_time;  // Time [μs]

    void clear() {
      flash_x.clear();
      flash_y.clear();
      flash_z.clear();

      flash_ophit_pe.clear();
      flash_ophit_ch.clear();
      flash_ophit_time.clear();
    }
  };

  // Processed data - filtered through TPC selection and quality cuts
  struct ProcessedData {
    // After TPC classification
    std::vector<std::vector<float>> flash_ophit_pe_sel;
    std::vector<std::vector<int>> flash_ophit_ch_sel;
    std::vector<std::vector<float>> flash_ophit_time_sel;
    std::vector<int> categorized_flashes;  // 0=TPC0, 1=TPC1
    std::vector<double> dEpromx_sel, dEpromy_sel, dEpromz_sel, dEtpc_sel;

    // After final quality cuts
    std::vector<std::vector<float>> flash_ophit_pe_final;
    std::vector<std::vector<int>> flash_ophit_ch_final;
    std::vector<std::vector<float>> flash_ophit_time_final;
    std::vector<double> nuvT_final, nuvZ_final;
    std::vector<double> dEpromx_final, dEpromy_final, dEpromz_final, dEtpc_final;

    // CNN inputs
    std::vector<std::vector<float>> pe_matrix;  // [flash][channel]
    std::vector<std::vector<std::vector<std::vector<float>>>> pe_images;  // [event][y][z][map]

    void clear() {
      flash_ophit_pe_sel.clear();
      flash_ophit_ch_sel.clear();
      flash_ophit_time_sel.clear();
      categorized_flashes.clear();
      dEpromx_sel.clear(); dEpromy_sel.clear(); dEpromz_sel.clear(); dEtpc_sel.clear();

      flash_ophit_pe_final.clear();
      flash_ophit_ch_final.clear();
      flash_ophit_time_final.clear();
      nuvT_final.clear(); nuvZ_final.clear();
      dEpromx_final.clear(); dEpromy_final.clear(); dEpromz_final.clear(); dEtpc_final.clear();

      pe_matrix.clear();
      pe_images.clear();
    }
  };

  // TTree output - flat structure for performance analysis
  struct TreeData {
    int run = 0, subrun = 0, event = 0;
    bool passedFilters = false;

    double trueX = -999.0, trueY = -999.0, trueZ = -999.0;  // MC truth [cm]
    double predX = -999.0, predY = -999.0, predZ = -999.0;  // CNN prediction [cm]
    double diffX = -999.0, diffY = -999.0, diffZ = -999.0;  // Residuals [cm]
    double error3D = -999.0;  // Euclidean error [cm]

    double nuvT = -999.0, nuvZ = -999.0, dEtpc = -999.0;

    void resetToDefaults() {
      run = subrun = event = 0;
      passedFilters = false;
      trueX = trueY = trueZ = -999.0;
      predX = predY = predZ = -999.0;
      diffX = diffY = diffZ = -999.0;
      error3D = -999.0;
      nuvT = nuvZ = dEtpc = -999.0;
    }
  };

}

// ============================================================================
// MAIN FUNCTION: PosRecoCVNProducer
// ============================================================================
//
// CNN-based 3D position reconstruction for neutrino interactions in SBND
//
// Art EDProducer module that reconstructs interaction positions using a
// convolutional neural network (CNN) trained on 2D PMT response images.
//
// Processing workflow:
// 1. Load OpFlashes and OpHits from photon detection system
// 2. Load MC truth information (if available)
// 3. Classify flashes by TPC volume (even/odd channel parity)
// 4. Apply quality filters (energy, PE thresholds, beam timing)
// 5. Generate 2D PMT images (coated/uncoated maps)
// 6. Execute TensorFlow CNN inference
// 7. Compute performance metrics (if ground truth available)
// 8. Store results in PixelMapVars data product and TTree
//
// Operating modes (configured via ProcessingMode FCL parameter):
// - MC_testing: Full validation with MC truth and energy filters
// - MC_inference: Inference with MC truth but minimal filtering
// - DATA_inference: Data inference without MC truth
//
// Output products:
// - PixelMapVars: Complete data product with images, predictions, metadata
// - TTree: Flat structure for inference performance analysis

class opdet::PosRecoCVNProducer : public art::EDProducer {
public:
  explicit PosRecoCVNProducer(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  PosRecoCVNProducer(PosRecoCVNProducer const&) = delete;
  PosRecoCVNProducer(PosRecoCVNProducer&&) = delete;
  PosRecoCVNProducer& operator=(PosRecoCVNProducer const&) = delete;
  PosRecoCVNProducer& operator=(PosRecoCVNProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // ========== Core processing functions ==========

  // Extract MC truth information from event (neutrino vertex, energy deposits)
  void FillMCTruth(art::Event const& e);

  // Initialize channel-to-detector-type mapping
  void InitializeChannelDict();

  // Classify channels into PMT/XARAPUCA and even/odd sets for TPC assignment
  void ClassifyChannels();

  // Determine TPC assignment (0/1)
  // Returns: 0 for even (TPC0), 1 for odd (TPC1)
  int CategorizeFirstChannel(const std::vector<int>& channels);

  // Apply flash selection based on TPC categorization
  void ApplyFlashSelection();

  // Apply final quality cuts (energy/PE thresholds depending on mode)
  void ApplyFinalEnergyFilter();

  // Create PE matrix [event][flash][channel] by summing PE across OpHits
  void CreatePEMatrix();

  // Load PMT spatial mapping files (coated/uncoated) for image generation
  void LoadPMTMaps();

  // Clear all event-level data structures (called at start of each event)
  void ClearEventData();

  // Select non-empty detector half (top or bottom) for CNN input
  // method: Selection criterion - "max" (most PE) or "sum"
  // Returns: Selected half-detector 2D array
  std::vector<std::vector<float>> SelectNonEmptyHalf(
    const std::vector<std::vector<float>>& left_half,
    const std::vector<std::vector<float>>& right_half,
    const std::string& method = "max");

  // Generate 2D CNN input images from PE matrix using PMT spatial maps
  void CreatePEImages();

  // Execute TensorFlow inference and store predictions
  // pixelmapvars: Output structure for predictions and metadata
  void RunInference(PixelMapVars& pixelmapvars);

  // Convert CNN output from normalized [0,1] to physical coordinates
  // scaled_predictions: Raw CNN outputs
  // Returns: Physical positions [cm]
  std::vector<double> ApplyInverseScaling(const std::vector<double>& scaled_predictions);

  // Fill TTree with event results (ground truth, predictions, metrics)
  // passedFilters: Whether event passed all quality cuts
  void FillInferenceTree(bool passedFilters);

  // Conditionally fill PixelMapVars data product (only if SavePixelMapVars=true)
  // pixelVars: Output data product
  // calculateGroundTruth: Include MC truth in output
  // passFilter: Event passed quality filters
  void FillPixelMapVarsConditional(PixelMapVars& pixelVars, bool calculateGroundTruth, bool passFilter);

  // Template function to filter 2D arrays using boolean mask
  template<typename T>
  std::vector<std::vector<T>> FilterByMask(
    const std::vector<std::vector<T>>& array,
    const std::vector<std::vector<bool>>& mask);

  // ========== Helper functions ==========

  // Log elapsed time for operation (if verbosity > 0)
  void LogTiming(const std::string& operation,
                 const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time);

  // Compute energy-weighted position centroids and spatial RMS from Geant4 steps
  // fenergydep: Energy deposited in each step [MeV]
  // fenergydepX,Y,Z: Spatial coordinates of deposits [cm]
  // fstepT: Time of each step [ns]
  // Outputs (by reference):
  //   dEtpc: Total deposited energy per TPC [MeV]
  //   dEpromx,y,z: Energy-weighted mean positions [cm]
  //   dEspreadx,y,z: RMS spatial spreads [cm]
  //   dElowedges, dEmaxedges: Bounding box of energy deposits [cm]
  void FillAverageDepositedEnergyVariables(
    std::vector<std::vector<double>> fenergydep,
    std::vector<std::vector<double>> fenergydepX,
    std::vector<std::vector<double>> fenergydepY,
    std::vector<std::vector<double>> fenergydepZ,
    std::vector<std::vector<double>> fstepT,
    std::vector<double> &dEtpc,
    std::vector<double> &dEpromx,
    std::vector<double> &dEpromy,
    std::vector<double> &dEpromz,
    std::vector<double> &dEspreadx,
    std::vector<double> &dEspready,
    std::vector<double> &dEspreadz,
    std::vector<std::vector<double>> &dElowedges,
    std::vector<std::vector<double>> &dEmaxedges);

  // Default value for invalid/missing MC truth data
  static constexpr double fDefaultSimIDE = -999.;

  // ========== FCL Configuration parameters ==========

  // MC truth selection
  std::vector<int> fMCTruthOrigin;              // Allowed neutrino origins (beam, cosmic)
  std::vector<int> fMCTruthPDG;                 // Allowed neutrino PDG codes
  std::vector<std::string> fMCTruthModuleLabel; // MC truth producer labels
  std::vector<std::string> fMCTruthInstanceLabel; // MC truth instance names
  std::string fMCModuleLabel;                 

  // Optical data input
  std::vector<std::string> fOpHitsModuleLabel;    // OpHit producer labels
  std::vector<std::string> fOpFlashesModuleLabel; // OpFlash producer labels

  // Detector geometry cuts
  std::vector<int> fG4BufferBoxX;               // Fiducial volume X boundaries [cm]
  std::vector<int> fG4BufferBoxY;               // Fiducial volume Y boundaries [cm]
  std::vector<int> fG4BufferBoxZ;               // Fiducial volume Z boundaries [cm]
  std::vector<int> fG4BeamWindow;               // Beam timing window [μs]
  std::vector<int> fKeepPDGCode;                // PDG codes to keep in MC particle list

  // Output configuration
  bool fSaveOpHits;                             // Include OpHit details in output
  bool fSavePixelMapVars;                       // Save full PixelMapVars data product
  int fVerbosity;                               // Verbosity level (0=quiet, 1=timing, 2+=debug)

  // PMT mapping files
  std::string fCoatedPMTMapPath;                // Path to coated PMT spatial map CSV
  std::string fUncoatedPMTMapPath;              // Path to uncoated PMT spatial map CSV

  // TensorFlow model configuration
  std::string fModelPath;                       // Path to frozen TensorFlow model (.pb)
  std::string fProcessingMode;                  // Processing mode: MC_testing, MC_inference, DATA_inference
  std::vector<std::string> fInputNames;         // TensorFlow input tensor names
  std::vector<std::string> fOutputNames;        // TensorFlow output tensor names
  double fCustomNormFactor;                     // Custom PE normalization factor (0=auto)
  double fPredictionTolerance;                  // Tolerance for out-of-range predictions (default 5%)
  bool fSkipNeutrinoFilter;                     // Skip neutrino presence filter for data
  std::string fSbndcodeVersion;                 // SBND code version string

  // CNN prediction scaling ranges: X[0,200], Y[-200,200], Z[0,500] cm

  // ========== Grouped data structures ==========
  //
  // Organizes related data into cohesive structures for easier access and manipulation
  EventInfo fEventInfo;           // Run/subrun/event identification
  MCTruthData fMCData;            // Monte Carlo truth information
  OpticalData fOpticalData;       // Raw PDS data (OpFlashes, OpHits)
  ProcessedData fProcessedData;   // Filtered data and CNN inputs
  TreeData fTreeData;             // TTree output variables

  // ========== Channel classification ==========
  std::map<int, int> fChannelDict;              // OpDetID → detector type mapping
  std::set<int> fPMTEven, fPMTOdd;              // PMT channels by parity (for TPC assignment)
  std::set<int> fXASEven, fXASOdd;              // XARAPUCA channels by parity

  // ========== PMT spatial maps ==========
  std::vector<std::vector<int>> fCoatedPMTMap;    // Coated PMT channel map [Y][Z]
  std::vector<std::vector<int>> fUncoatedPMTMap;  // Uncoated PMT channel map [Y][Z]

  // ========== TensorFlow inference ==========
  std::unique_ptr<tf::Graph> fTFGraph;          // TensorFlow graph for CNN inference

  // ========== Output TTree ==========
  TTree* fInferenceTree;                        // Simple TTree for performance analysis
};


DEFINE_ART_MODULE(opdet::PosRecoCVNProducer)
