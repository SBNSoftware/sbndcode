////////////////////////////////////////////////////////////////////////
// Class:       SBNDPDSAnalyzer
// Plugin Type: analyzer
// File:        SBNDPDSAnalyzer_module.hh
//
// Created by Francisco Nicolas-Arnaldos using cetskelgen
////////////////////////////////////////////////////////////////////////

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
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// ROOT and C++ includes
#include <TTree.h>
#include <string.h>

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
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"


#define xdet_size 1000
#define ydet_size 1000
#define zmindet_size -500
#define zmaxdet_size 1800

namespace opdet {
  class SBNDPDSAnalyzer;

  enum SBNDPDSDetectorType {
    kPDUnknown = -1,   
    kPMTCoated = 0,       
    kPMTUncoated = 1,
    kXARAPUCAVUV,
    kXARAPUCAVIS
  };

}


class opdet::SBNDPDSAnalyzer : public art::EDAnalyzer {
public:
  explicit SBNDPDSAnalyzer(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  SBNDPDSAnalyzer(SBNDPDSAnalyzer const&) = delete;
  SBNDPDSAnalyzer(SBNDPDSAnalyzer&&) = delete;
  SBNDPDSAnalyzer& operator=(SBNDPDSAnalyzer const&) = delete;
  SBNDPDSAnalyzer& operator=(SBNDPDSAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Functions
  void FillMCTruth(art::Event const& e);

  void FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
  std::vector<std::vector<double>> fenergydepY, std::vector<std::vector<double>> fenergydepZ, std::vector<std::vector<double>> fstepT,
  std::vector<double> &dEtpc, std::vector<double> &dEpromx, std::vector<double> &dEpromy, std::vector<double> &dEpromz,
  std::vector<double> &dEspreadx, std::vector<double> &dEspready, std::vector<double> &dEspreadz,
  std::vector<std::vector<double>> &dElowedges, std::vector<std::vector<double>> &dEmaxedges);

  void ResetSimPhotons();

  void FillSimPhotons(std::vector<art::Handle<std::vector<sim::SimPhotons> >> photonHandle_list);

  void FillSimPhotonsLite(std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> photonHandle_list);

  simb::MCTruth MakeTruthMatching(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &recoHits, std::map<std::string, int>& allHitsTruthMap, unsigned int TPC, double& purity, double& completeness, simb::MCParticle &mainMCParticle);

  std::map<std::string, int> GetAllHitsTruthMatch(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &allHits);

  // TTree saving options
  bool fSaveMCTruth;
  bool fSaveMCParticles;
  bool fSaveSimPhotons;
  bool fSaveSimPhotonsArrivalTimes;
  bool fSaveRawWaveforms;
  bool fSaveDeconvolvedWaveforms;
  bool fSaveOpHits;
  bool fSaveOpFlashes;
  bool fSaveCosmicId;

  // Configuration parameters
  int fVerbosity;
  bool fMakePerTrackTree;
  bool fMakePDSGeoTree;
  bool fUseSimPhotonsLite;
  std::vector <std::string> fPDTypes;
  std::vector <int> fKeepPDGCode;
  std::vector <int> fMCTruthOrigin;
  std::vector <int> fMCTruthPDG;

  // Data product input labels
  std::vector<std::string> fMCTruthModuleLabel;
  std::vector<std::string> fMCTruthInstanceLabel;
  std::string fMCModuleLabel;
  std::vector<std::string> fSimPhotonsModuleLabel;
  std::string fRawWaveformsModuleLabel;
  std::string fDeconvolvedWaveformsModuleLabel;
  std::vector<std::string> fOpHitsModuleLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  std::string fHitsLabel;
  std::string fReco2Label;
  std::string fCosmicIdModuleLabel;
  std::string fOpT0FinderModuleLabel;
  std::string fSimpleFlashMatchModuleLabel;

  // Fiducial volume for MC Particles
  std::vector<int> fG4BufferBoxX;
  std::vector<int> fG4BufferBoxY;
  std::vector<int> fG4BufferBoxZ;
  std::vector<int> fG4BeamWindow;

  // PDS mapping and geometry
  opdet::sbndPDMapAlg fPDSMap;
  art::ServiceHandle<geo::Geometry> fGeoService;

  static constexpr double fDefaultSimIDE = -999.;

  // TTrees
  TTree *fPDSMapTree;
  TTree *fTree;
  TTree *fPerTrackTree;
  
  unsigned int _eventID, _runID, _subrunID;
  
  // Saving MCTruth
  std::vector <double> _nuvT;
  std::vector <double> _nuvX;
  std::vector <double> _nuvY;
  std::vector <double> _nuvZ;
  std::vector <double> _nuvE;
  
  // Saving MCParticles                       
  std::vector < std::vector <double> > _mc_stepX;
  std::vector < std::vector <double> > _mc_stepY;
  std::vector < std::vector <double> > _mc_stepZ;
  std::vector < std::vector <double> > _mc_stepT;
  std::vector <int> _mc_PDGcode;
  std::vector <int> _mc_trackID;
  std::vector <int> _mc_motherID;
  std::vector <std::string> _mc_process;
  std::vector <double> _mc_E;
  std::vector <double> _mc_dE;
  std::vector <double> _mc_StartPx;
  std::vector <double> _mc_StartPy;
  std::vector <double> _mc_StartPz;
  std::vector <double> _mc_EndPx;
  std::vector <double> _mc_EndPy;
  std::vector <double> _mc_EndPz;
  unsigned int _mc_InTimeCosmics;
  std::vector <double> _mc_InTimeCosmicsTime;
  std::vector < std::vector <double> > _mc_energydep;
  std::vector < std::vector <double> > _mc_energydepX;
  std::vector < std::vector <double> > _mc_energydepY;
  std::vector < std::vector <double> > _mc_energydepZ;
  std::vector<double> _mc_dEtpc, _mc_dEpromx, _mc_dEpromy, _mc_dEpromz;
  std::vector<double> _mc_dEspreadx, _mc_dEspready, _mc_dEspreadz;
  std::vector<std::vector<double>> _mc_dElowedges, _mc_dEmaxedges;

  // Saving SimPhotons (VUV and VIS)
  std::vector < std::vector <double> > _simPhotonsLiteVUV;
  std::vector < std::vector <double> > _simPhotonsLiteVIS;
  std::vector <double> _simPhotonsperOpChVUV;
  std::vector <double> _simPhotonsperOpChVIS;
  int _NPhotons;
  int _NPhotonsPMTCo;
  int _NPhotonsPMTUnco;
  int _NPhotonsPMTCoVUV;
  int _NPhotonsXARAPUCAVUV;
  int _NPhotonsXARAPUCAVIS;

  // Saving raw signals
  std::vector < std::vector <double> > _signalsDigi;
  std::vector <double> _stampTime;
  std::vector <int> _opChDigi;

  // Saving deconvolved signals
  std::vector < std::vector <double> > _signalsDeco;
  std::vector <double> _stampTimeDeco;
  std::vector <int> _opChDeco;

  // Saving all the OpHits
  int _nophits;
  std::vector<int>  _ophit_opch;           ///< OpChannel of the optical hit
  std::vector<double>  _ophit_peakT;       ///< Peak time of the optical hit
  std::vector<double>  _ophit_startT;       ///< Peak time of the optical hit
  std::vector<double>  _ophit_riseT;       ///< Peak time of the optical hit
  std::vector<double> _ophit_width;       ///< Width of the optical hit
  std::vector<double> _ophit_area;        ///< Area of the optical hit
  std::vector<double> _ophit_amplitude;   ///< Amplitude of the optical hit
  std::vector<double> _ophit_pe;          ///< PEs of the optical hit

  // Saving OpFlash
  int _nopflash;
  std::vector<int> _flash_id;
  std::vector<double> _flash_time;
  std::vector<double> _flash_total_pe;
  std::vector<std::vector<double>> _flash_pe_v;
  std::vector<double> _flash_y;
  std::vector<double> _flash_yerr ;
  std::vector<double> _flash_z;
  std::vector<double> _flash_zerr;
  std::vector<double> _flash_x;
  std::vector<double> _flash_xerr;
  std::vector<int> _flash_tpc;
  std::vector<std::vector<double>> _flash_ophit_time;
  std::vector<std::vector<double>> _flash_ophit_risetime;
  std::vector<std::vector<double>> _flash_ophit_starttime;
  std::vector<std::vector<double>>_flash_ophit_amp;
  std::vector<std::vector<double>> _flash_ophit_area;
  std::vector<std::vector<double>> _flash_ophit_width;
  std::vector<std::vector<double>> _flash_ophit_pe;
  std::vector<std::vector<int>> _flash_ophit_ch;

  // Cosmic ID
  std::vector<double> _CRUMBSScore;
  std::vector<int> _sliceOrigin;
  std::vector<double> _sliceCompleteness;
  std::vector<double> _slicePurity;

  // OpT0Finder
  std::vector<double> _opT0Chi2;
  std::vector<double> _opT0Score;
  std::vector<double> _opT0Time;
  std::vector<double> _opT0HypoPE;
  std::vector<double> _opT0MeasPE;

  // SimpleFlashMatch
  std::vector<int> _sFMSliceOrigin;
  std::vector<double> _sFMScore;
  std::vector<double> _sFMScoreY;
  std::vector<double> _sFMScoreZ;
  std::vector<double> _sFMScoreRR;
  std::vector<double> _sFMScoreRatio;
  std::vector<double> _sFMScoreSlope;
  std::vector<double> _sFMScorePEtoQ;
  std::vector<double> _sFMTime;
  std::vector<double> _sFMPE;

  // Make per track ID
  // Saving SimPhotons (VUV and VIS) per track id
  int _perTrackID;
  std::vector < std::vector <double> > _simPhotonsVUV;
  std::vector < std::vector <double> > _simPhotonsVIS;
  // maps to store simphotons per trackID
  std::map<int, std::vector < std::vector <double> > > fSimPhotonsVUVMap;
  std::map<int, std::vector < std::vector <double> > > fSimPhotonsVISMap;


  // PDS geo tree
  std::vector <int> _opDetID;
  std::vector <double> _opDetX;
  std::vector <double> _opDetY;
  std::vector <double> _opDetZ;
  std::vector <int> _opDetType;

};


DEFINE_ART_MODULE(opdet::SBNDPDSAnalyzer)
