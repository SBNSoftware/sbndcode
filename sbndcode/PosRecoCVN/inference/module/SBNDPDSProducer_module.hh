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

// TensorFlow includes
#include "sbndcode/PosRecoCVN/tf/tf_graph.h"
#include <sstream>

// Include the data structure definition
#include "sbndcode/PosRecoCVN/module/PixelMapVars.h"

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// G4 includes
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
  class SBNDPDSProducer;

  enum SBNDPDSDetectorType {
    kPDUnknown = -1,   
    kPMTCoated = 0,       
    kPMTUncoated = 1,
    kXARAPUCAVUV,
    kXARAPUCAVIS
  };

}


class opdet::SBNDPDSProducer : public art::EDProducer {
public:
  explicit SBNDPDSProducer(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  SBNDPDSProducer(SBNDPDSProducer const&) = delete;
  SBNDPDSProducer(SBNDPDSProducer&&) = delete;
  SBNDPDSProducer& operator=(SBNDPDSProducer const&) = delete;
  SBNDPDSProducer& operator=(SBNDPDSProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  // void beginJob() override;

private:

  // Functions
  void FillMCTruth(art::Event const& e);
  void InitializeChannelDict();
  void ClassifyChannels();
  int CategorizeFirstChannel(const std::vector<int>& channels);
  void ApplyFlashSelection();
  void ApplyFinalEnergyFilter();
  void CreatePEMatrix();
  void LoadPMTMaps();
  void ClearEventData();
  std::vector<std::vector<float>> SelectNonEmptyHalf(const std::vector<std::vector<float>>& left_half, 
                                                     const std::vector<std::vector<float>>& right_half,
                                                     const std::string& method = "max");
  void CreatePEImages();
  void RunInference(PixelMapVars& pixelmapvars);
  std::vector<double> ApplyInverseScaling(const std::vector<double>& scaled_predictions);
  template<typename T>
  std::vector<std::vector<T>> FilterByMask(const std::vector<std::vector<T>>& array, const std::vector<std::vector<bool>>& mask);

  void FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
  std::vector<std::vector<double>> fenergydepY, std::vector<std::vector<double>> fenergydepZ, std::vector<std::vector<double>> fstepT,
  std::vector<double> &dEtpc, std::vector<double> &dEpromx, std::vector<double> &dEpromy, std::vector<double> &dEpromz,
  std::vector<double> &dEspreadx, std::vector<double> &dEspready, std::vector<double> &dEspreadz,
  std::vector<std::vector<double>> &dElowedges, std::vector<std::vector<double>> &dEmaxedges);

  static constexpr double fDefaultSimIDE = -999.;

  // Parámetros de configuración (ordenados para coincidir con el constructor)
  std::vector <int> fMCTruthOrigin;
  std::vector <int> fMCTruthPDG;
  std::vector<std::string> fMCTruthModuleLabel;
  std::vector<std::string> fMCTruthInstanceLabel;
  std::string fMCModuleLabel;
  std::vector<std::string> fOpHitsModuleLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  std::vector<int> fG4BufferBoxX;
  std::vector<int> fG4BufferBoxY;
  std::vector<int> fG4BufferBoxZ;
  std::vector<int> fG4BeamWindow;
  std::vector<int> fKeepPDGCode;
  bool fSaveOpHits;
  int fVerbosity;
  std::string fCoatedPMTMapPath;
  std::string fUncoatedPMTMapPath;
  
  // TensorFlow model parameters
  std::string fModelPath;
  bool fRunInference;
  std::vector<std::string> fInputNames;
  std::vector<std::string> fOutputNames;
  double fCustomNormFactor;
  
  // Simple scaling ranges: X[0,200], Y[-200,200], Z[0,500]

  // Variables internas necesarias
  std::vector<double> _nuvT;
  std::vector<double> _nuvZ;
  std::vector<double> _mc_dEpromx, _mc_dEpromy, _mc_dEpromz, _mc_dEtpc;
  std::vector<std::vector<float>> _flash_ophit_pe;
  std::vector<std::vector<int>> _flash_ophit_ch;
  std::vector<std::vector<float>> _flash_ophit_time;

  // --- Variables agregadas para corregir errores de compilación ---
  // Identificadores de evento
  int _eventID;
  int _runID;
  int _subrunID;

  // Variables de neutrinos
  std::vector<double> _nuvX, _nuvY, _nuvE;

  // Variables de MCParticles
  std::vector<std::vector<double>> _mc_stepX, _mc_stepY, _mc_stepZ, _mc_stepT;
  std::vector<double> _mc_dE, _mc_E;
  std::vector<int> _mc_trackID, _mc_motherID, _mc_PDGcode;
  std::vector<std::string> _mc_process;
  std::vector<double> _mc_StartPx, _mc_StartPy, _mc_StartPz;
  std::vector<double> _mc_EndPx, _mc_EndPy, _mc_EndPz;
  std::vector<std::vector<double>> _mc_energydep, _mc_energydepX, _mc_energydepY, _mc_energydepZ;
  int _mc_InTimeCosmics;
  std::vector<double> _mc_InTimeCosmicsTime;

  // Variables para energía depositada
  std::vector<double> _mc_dEspreadx, _mc_dEspready, _mc_dEspreadz;
  std::vector<std::vector<double>> _mc_dElowedges, _mc_dEmaxedges;

  // Variables para OpHits
  int _nophits;
  std::vector<int> _ophit_opch;
  std::vector<double> _ophit_peakT, _ophit_startT, _ophit_riseT, _ophit_width, _ophit_area, _ophit_amplitude, _ophit_pe;

  // Variables para OpFlashes
  int _nopflash;
  std::vector<int> _flash_id;
  std::vector<double> _flash_time, _flash_total_pe;
  std::vector<std::vector<double>> _flash_pe_v;
  std::vector<int> _flash_tpc;
  std::vector<double> _flash_y, _flash_yerr, _flash_z, _flash_zerr, _flash_x, _flash_xerr;

  // Variables para OpHit en flashes
  std::vector<std::vector<float>> _flash_ophit_risetime, _flash_ophit_starttime, _flash_ophit_amp, _flash_ophit_area, _flash_ophit_width;

  // Variables auxiliares
  double dE_neutrinowindow;

  // Vectores auxiliares
  std::vector<simb::MCParticle> mcpartVec;
  
  // Channel dictionary mapping OpDetID to OpDetType
  std::map<int, int> fChannelDict;
  
  // Channel classification sets
  std::set<int> fPMTEven, fPMTOdd, fXASEven, fXASOdd;
  
  // Selected flash data after classification
  std::vector<std::vector<float>> _flash_ophit_pe_sel;
  std::vector<std::vector<int>> _flash_ophit_ch_sel;
  std::vector<std::vector<float>> _flash_ophit_time_sel;
  std::vector<int> _categorized_flashes;
  std::vector<double> _mc_dEpromx_sel, _mc_dEpromy_sel, _mc_dEpromz_sel, _mc_dEtpc_sel;
  
  // Final filtered data after energy deposition cuts
  std::vector<std::vector<float>> _flash_ophit_pe_final;
  std::vector<std::vector<int>> _flash_ophit_ch_final;
  std::vector<std::vector<float>> _flash_ophit_time_final;
  std::vector<double> _nuvT_final, _nuvZ_final;
  std::vector<double> _mc_dEpromx_final, _mc_dEpromy_final, _mc_dEpromz_final, _mc_dEtpc_final;
  
  // PE matrix
  std::vector<std::vector<float>> _pe_matrix;
  
  // PMT maps
  std::vector<std::vector<int>> _coated_pmt_map;
  std::vector<std::vector<int>> _uncoated_pmt_map;
  
  // Generated images
  std::vector<std::vector<std::vector<std::vector<float>>>> _pe_images;
  
  // TensorFlow inference variables
  std::unique_ptr<tf::Graph> fTFGraph;
};


DEFINE_ART_MODULE(opdet::SBNDPDSProducer)
