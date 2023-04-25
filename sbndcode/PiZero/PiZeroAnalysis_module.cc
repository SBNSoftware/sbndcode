////////////////////////////////////////////////////////////////////////
// Class:       PiZeroAnalysis
// Plugin Type: analyzer
// File:        PiZeroAnalysis_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

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

#include "TTree.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include <numeric>

constexpr int def_int      = -std::numeric_limits<int>::max();
constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class PiZeroAnalysis;
}

class sbnd::PiZeroAnalysis : public art::EDAnalyzer {
public:
  explicit PiZeroAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PiZeroAnalysis(PiZeroAnalysis const&) = delete;
  PiZeroAnalysis(PiZeroAnalysis&&) = delete;
  PiZeroAnalysis& operator=(PiZeroAnalysis const&) = delete;
  PiZeroAnalysis& operator=(PiZeroAnalysis&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;

  void ClearMaps(const art::Event &e);

  void SetupMaps(const art::Event &e);

  bool SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalysePFPs(const art::Event &e, const art::Handle<std::vector<recob::PFParticle>> &PFParticleHandle);

  void AnalyseMCParticles(const art::Event &e, const std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec, const art::FindManyP<simb::MCTruth> &MCParticlesToMCTruths);

  void FillGammaTruth(const int i, const simb::MCParticle* gamma, const bool lead);

  void FillGammaTruthDefault(const int i);

  void FillGammaReco(const int i, const art::Event &e, const int trueID, const bool lead);

  void FillGammaRecoDefault(const int i);

  void FillBestReco(const int i);

  void ResizeVectors(const int size);

  void ResizeRecoSubVectors(const int i, const int size, const bool lead);

  double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  std::string fMCParticleModuleLabel, fPFParticleModuleLabel, fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fTrackCalorimetryModuleLabel;
  bool fDebug;

  std::map<int,int> fHitsMap;
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;

  TTree* fPiZeroTree;

  // Tree variables

  int  _run;
  int  _subrun;
  int  _event;
  bool _signal;

  std::vector<int>    _mct_origin;
  std::vector<bool>   _mc_primary;
  std::vector<double> _mc_vx;
  std::vector<double> _mc_vy;
  std::vector<double> _mc_vz;
  std::vector<double> _mc_vt;
  std::vector<bool>   _mc_av;
  std::vector<bool>   _mc_fv;
  std::vector<bool>   _mc_uboone_fv;
  std::vector<double> _mc_endx;
  std::vector<double> _mc_endy;
  std::vector<double> _mc_endz;
  std::vector<double> _mc_endt;
  std::vector<double> _mc_vpx;
  std::vector<double> _mc_vpy;
  std::vector<double> _mc_vpz;
  std::vector<double> _mc_vtheta;
  std::vector<double> _mc_vphi;
  std::vector<double> _mc_ve;
  std::vector<double> _mc_endpx;
  std::vector<double> _mc_endpy;
  std::vector<double> _mc_endpz;
  std::vector<double> _mc_ende;

  std::vector<bool>   _mc_digamma;

  std::vector<int>    _mc_gamma1_pdg;
  std::vector<int>    _mc_gamma1_status;
  std::vector<double> _mc_gamma1_vx;
  std::vector<double> _mc_gamma1_vy;
  std::vector<double> _mc_gamma1_vz;
  std::vector<double> _mc_gamma1_vt;
  std::vector<double> _mc_gamma1_endx;
  std::vector<double> _mc_gamma1_endy;
  std::vector<double> _mc_gamma1_endz;
  std::vector<double> _mc_gamma1_endt;
  std::vector<double> _mc_gamma1_vpx;
  std::vector<double> _mc_gamma1_vpy;
  std::vector<double> _mc_gamma1_vpz;
  std::vector<double> _mc_gamma1_vtheta;
  std::vector<double> _mc_gamma1_vphi;
  std::vector<double> _mc_gamma1_ve;
  std::vector<double> _mc_gamma1_endpx;
  std::vector<double> _mc_gamma1_endpy;
  std::vector<double> _mc_gamma1_endpz;
  std::vector<double> _mc_gamma1_ende;

  std::vector<int>    _mc_gamma2_pdg;
  std::vector<int>    _mc_gamma2_status;
  std::vector<double> _mc_gamma2_vx;
  std::vector<double> _mc_gamma2_vy;
  std::vector<double> _mc_gamma2_vz;
  std::vector<double> _mc_gamma2_vt;
  std::vector<double> _mc_gamma2_endx;
  std::vector<double> _mc_gamma2_endy;
  std::vector<double> _mc_gamma2_endz;
  std::vector<double> _mc_gamma2_endt;
  std::vector<double> _mc_gamma2_vpx;
  std::vector<double> _mc_gamma2_vpy;
  std::vector<double> _mc_gamma2_vpz;
  std::vector<double> _mc_gamma2_vtheta;
  std::vector<double> _mc_gamma2_vphi;
  std::vector<double> _mc_gamma2_ve;
  std::vector<double> _mc_gamma2_endpx;
  std::vector<double> _mc_gamma2_endpy;
  std::vector<double> _mc_gamma2_endpz;
  std::vector<double> _mc_gamma2_ende;

  std::vector<double> _mc_open_angle;

  std::vector<int> _reco_gamma1_nhits;
  std::vector<int> _reco_gamma1_npfps;

  std::vector<std::vector<double>> _reco_gamma1_pfp_track_score;
  std::vector<std::vector<double>> _reco_gamma1_comp;
  std::vector<std::vector<double>> _reco_gamma1_pur;

  std::vector<std::vector<bool>>   _reco_gamma1_has_good_track;
  std::vector<std::vector<double>> _reco_gamma1_track_start_x;
  std::vector<std::vector<double>> _reco_gamma1_track_start_y;
  std::vector<std::vector<double>> _reco_gamma1_track_start_z;
  std::vector<std::vector<double>> _reco_gamma1_track_end_x;
  std::vector<std::vector<double>> _reco_gamma1_track_end_y;
  std::vector<std::vector<double>> _reco_gamma1_track_end_z;
  std::vector<std::vector<double>> _reco_gamma1_track_dir_x;
  std::vector<std::vector<double>> _reco_gamma1_track_dir_y;
  std::vector<std::vector<double>> _reco_gamma1_track_dir_z;
  std::vector<std::vector<double>> _reco_gamma1_track_start_theta;
  std::vector<std::vector<double>> _reco_gamma1_track_start_phi;
  std::vector<std::vector<double>> _reco_gamma1_track_length;
  std::vector<std::vector<double>> _reco_gamma1_track_best_plane_energy;

  std::vector<std::vector<bool>>   _reco_gamma1_has_good_shower;
  std::vector<std::vector<double>> _reco_gamma1_shower_start_x;
  std::vector<std::vector<double>> _reco_gamma1_shower_start_y;
  std::vector<std::vector<double>> _reco_gamma1_shower_start_z;
  std::vector<std::vector<double>> _reco_gamma1_shower_dir_x;
  std::vector<std::vector<double>> _reco_gamma1_shower_dir_y;
  std::vector<std::vector<double>> _reco_gamma1_shower_dir_z;
  std::vector<std::vector<double>> _reco_gamma1_shower_start_theta;
  std::vector<std::vector<double>> _reco_gamma1_shower_start_phi;
  std::vector<std::vector<double>> _reco_gamma1_shower_length;
  std::vector<std::vector<double>> _reco_gamma1_shower_open_angle;
  std::vector<std::vector<double>> _reco_gamma1_shower_best_plane_energy;
  std::vector<std::vector<double>> _reco_gamma1_shower_best_plane_dedx;

  std::vector<double> _reco_gamma1_bestpfp_pfp_track_score;
  std::vector<double> _reco_gamma1_bestpfp_comp;
  std::vector<double> _reco_gamma1_bestpfp_pur;

  std::vector<bool>   _reco_gamma1_bestpfp_has_good_track;
  std::vector<double> _reco_gamma1_bestpfp_track_start_x;
  std::vector<double> _reco_gamma1_bestpfp_track_start_y;
  std::vector<double> _reco_gamma1_bestpfp_track_start_z;
  std::vector<double> _reco_gamma1_bestpfp_track_end_x;
  std::vector<double> _reco_gamma1_bestpfp_track_end_y;
  std::vector<double> _reco_gamma1_bestpfp_track_end_z;
  std::vector<double> _reco_gamma1_bestpfp_track_dir_x;
  std::vector<double> _reco_gamma1_bestpfp_track_dir_y;
  std::vector<double> _reco_gamma1_bestpfp_track_dir_z;
  std::vector<double> _reco_gamma1_bestpfp_track_start_theta;
  std::vector<double> _reco_gamma1_bestpfp_track_start_phi;
  std::vector<double> _reco_gamma1_bestpfp_track_length;
  std::vector<double> _reco_gamma1_bestpfp_track_best_plane_energy;

  std::vector<bool>   _reco_gamma1_bestpfp_has_good_shower;
  std::vector<double> _reco_gamma1_bestpfp_shower_start_x;
  std::vector<double> _reco_gamma1_bestpfp_shower_start_y;
  std::vector<double> _reco_gamma1_bestpfp_shower_start_z;
  std::vector<double> _reco_gamma1_bestpfp_shower_dir_x;
  std::vector<double> _reco_gamma1_bestpfp_shower_dir_y;
  std::vector<double> _reco_gamma1_bestpfp_shower_dir_z;
  std::vector<double> _reco_gamma1_bestpfp_shower_start_theta;
  std::vector<double> _reco_gamma1_bestpfp_shower_start_phi;
  std::vector<double> _reco_gamma1_bestpfp_shower_length;
  std::vector<double> _reco_gamma1_bestpfp_shower_open_angle;
  std::vector<double> _reco_gamma1_bestpfp_shower_best_plane_energy;
  std::vector<double> _reco_gamma1_bestpfp_shower_best_plane_dedx;

  std::vector<int> _reco_gamma2_nhits;
  std::vector<int> _reco_gamma2_npfps;

  std::vector<std::vector<double>> _reco_gamma2_pfp_track_score;
  std::vector<std::vector<double>> _reco_gamma2_comp;
  std::vector<std::vector<double>> _reco_gamma2_pur;

  std::vector<std::vector<bool>>   _reco_gamma2_has_good_track;
  std::vector<std::vector<double>> _reco_gamma2_track_start_x;
  std::vector<std::vector<double>> _reco_gamma2_track_start_y;
  std::vector<std::vector<double>> _reco_gamma2_track_start_z;
  std::vector<std::vector<double>> _reco_gamma2_track_end_x;
  std::vector<std::vector<double>> _reco_gamma2_track_end_y;
  std::vector<std::vector<double>> _reco_gamma2_track_end_z;
  std::vector<std::vector<double>> _reco_gamma2_track_dir_x;
  std::vector<std::vector<double>> _reco_gamma2_track_dir_y;
  std::vector<std::vector<double>> _reco_gamma2_track_dir_z;
  std::vector<std::vector<double>> _reco_gamma2_track_start_theta;
  std::vector<std::vector<double>> _reco_gamma2_track_start_phi;
  std::vector<std::vector<double>> _reco_gamma2_track_length;
  std::vector<std::vector<double>> _reco_gamma2_track_best_plane_energy;

  std::vector<std::vector<bool>>   _reco_gamma2_has_good_shower;
  std::vector<std::vector<double>> _reco_gamma2_shower_start_x;
  std::vector<std::vector<double>> _reco_gamma2_shower_start_y;
  std::vector<std::vector<double>> _reco_gamma2_shower_start_z;
  std::vector<std::vector<double>> _reco_gamma2_shower_dir_x;
  std::vector<std::vector<double>> _reco_gamma2_shower_dir_y;
  std::vector<std::vector<double>> _reco_gamma2_shower_dir_z;
  std::vector<std::vector<double>> _reco_gamma2_shower_start_theta;
  std::vector<std::vector<double>> _reco_gamma2_shower_start_phi;
  std::vector<std::vector<double>> _reco_gamma2_shower_length;
  std::vector<std::vector<double>> _reco_gamma2_shower_open_angle;
  std::vector<std::vector<double>> _reco_gamma2_shower_best_plane_energy;
  std::vector<std::vector<double>> _reco_gamma2_shower_best_plane_dedx;

  std::vector<double> _reco_gamma2_bestpfp_pfp_track_score;
  std::vector<double> _reco_gamma2_bestpfp_comp;
  std::vector<double> _reco_gamma2_bestpfp_pur;

  std::vector<bool>   _reco_gamma2_bestpfp_has_good_track;
  std::vector<double> _reco_gamma2_bestpfp_track_start_x;
  std::vector<double> _reco_gamma2_bestpfp_track_start_y;
  std::vector<double> _reco_gamma2_bestpfp_track_start_z;
  std::vector<double> _reco_gamma2_bestpfp_track_end_x;
  std::vector<double> _reco_gamma2_bestpfp_track_end_y;
  std::vector<double> _reco_gamma2_bestpfp_track_end_z;
  std::vector<double> _reco_gamma2_bestpfp_track_dir_x;
  std::vector<double> _reco_gamma2_bestpfp_track_dir_y;
  std::vector<double> _reco_gamma2_bestpfp_track_dir_z;
  std::vector<double> _reco_gamma2_bestpfp_track_start_theta;
  std::vector<double> _reco_gamma2_bestpfp_track_start_phi;
  std::vector<double> _reco_gamma2_bestpfp_track_length;
  std::vector<double> _reco_gamma2_bestpfp_track_best_plane_energy;

  std::vector<bool>   _reco_gamma2_bestpfp_has_good_shower;
  std::vector<double> _reco_gamma2_bestpfp_shower_start_x;
  std::vector<double> _reco_gamma2_bestpfp_shower_start_y;
  std::vector<double> _reco_gamma2_bestpfp_shower_start_z;
  std::vector<double> _reco_gamma2_bestpfp_shower_dir_x;
  std::vector<double> _reco_gamma2_bestpfp_shower_dir_y;
  std::vector<double> _reco_gamma2_bestpfp_shower_dir_z;
  std::vector<double> _reco_gamma2_bestpfp_shower_start_theta;
  std::vector<double> _reco_gamma2_bestpfp_shower_start_phi;
  std::vector<double> _reco_gamma2_bestpfp_shower_length;
  std::vector<double> _reco_gamma2_bestpfp_shower_open_angle;
  std::vector<double> _reco_gamma2_bestpfp_shower_best_plane_energy;
  std::vector<double> _reco_gamma2_bestpfp_shower_best_plane_dedx;
};

sbnd::PiZeroAnalysis::PiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCParticleModuleLabel       = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fPFParticleModuleLabel       = p.get<std::string>("PFParticleModuleLabel", "pandoraSCE");
    fHitModuleLabel              = p.get<std::string>("HitModuleLabel", "gaushit");
    fTrackModuleLabel            = p.get<std::string>("TrackModuleLabel", "pandoraSCETrack");
    fShowerModuleLabel           = p.get<std::string>("ShowerModuleLabel", "pandoraSCEShower");
    fTrackCalorimetryModuleLabel = p.get<std::string>("TrackCalorimetryModuleLabel", "pandoraSCECalo");
    fDebug                       = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fPiZeroTree = fs->make<TTree>("pizeros","");
    fPiZeroTree->Branch("run", &_run);
    fPiZeroTree->Branch("subrun", &_subrun);
    fPiZeroTree->Branch("event", &_event);
    fPiZeroTree->Branch("signal", &_signal);

    fPiZeroTree->Branch("mct_origin", "std::vector<int>", &_mct_origin);
    fPiZeroTree->Branch("mc_primary", "std::vector<bool>", &_mc_primary);
    fPiZeroTree->Branch("mc_vx", "std::vector<double>", &_mc_vx);
    fPiZeroTree->Branch("mc_vy", "std::vector<double>", &_mc_vy);
    fPiZeroTree->Branch("mc_vz", "std::vector<double>", &_mc_vz);
    fPiZeroTree->Branch("mc_vt", "std::vector<double>", &_mc_vt);
    fPiZeroTree->Branch("mc_av", "std::vector<bool>", &_mc_av);
    fPiZeroTree->Branch("mc_fv", "std::vector<bool>", &_mc_fv);
    fPiZeroTree->Branch("mc_uboone_fv", "std::vector<bool>", &_mc_uboone_fv);
    fPiZeroTree->Branch("mc_endx", "std::vector<double>", &_mc_endx);
    fPiZeroTree->Branch("mc_endy", "std::vector<double>", &_mc_endy);
    fPiZeroTree->Branch("mc_endz", "std::vector<double>", &_mc_endz);
    fPiZeroTree->Branch("mc_endt", "std::vector<double>", &_mc_endt);
    fPiZeroTree->Branch("mc_vpx", "std::vector<double>", &_mc_vpx);
    fPiZeroTree->Branch("mc_vpy", "std::vector<double>", &_mc_vpy);
    fPiZeroTree->Branch("mc_vpz", "std::vector<double>", &_mc_vpz);
    fPiZeroTree->Branch("mc_vtheta", "std::vector<double>", &_mc_vtheta);
    fPiZeroTree->Branch("mc_vphi", "std::vector<double>", &_mc_vphi);
    fPiZeroTree->Branch("mc_ve", "std::vector<double>", &_mc_ve);
    fPiZeroTree->Branch("mc_endpx", "std::vector<double>", &_mc_endpx);
    fPiZeroTree->Branch("mc_endpy", "std::vector<double>", &_mc_endpy);
    fPiZeroTree->Branch("mc_endpz", "std::vector<double>", &_mc_endpz);
    fPiZeroTree->Branch("mc_ende", "std::vector<double>", &_mc_ende);

    fPiZeroTree->Branch("mc_digamma", "std::vector<bool>", &_mc_digamma);

    fPiZeroTree->Branch("mc_gamma1_pdg", "std::vector<int>", &_mc_gamma1_pdg);
    fPiZeroTree->Branch("mc_gamma1_status", "std::vector<int>", &_mc_gamma1_status);
    fPiZeroTree->Branch("mc_gamma1_vx", "std::vector<double>", &_mc_gamma1_vx);
    fPiZeroTree->Branch("mc_gamma1_vy", "std::vector<double>", &_mc_gamma1_vy);
    fPiZeroTree->Branch("mc_gamma1_vz", "std::vector<double>", &_mc_gamma1_vz);
    fPiZeroTree->Branch("mc_gamma1_vt", "std::vector<double>", &_mc_gamma1_vt);
    fPiZeroTree->Branch("mc_gamma1_endx", "std::vector<double>", &_mc_gamma1_endx);
    fPiZeroTree->Branch("mc_gamma1_endy", "std::vector<double>", &_mc_gamma1_endy);
    fPiZeroTree->Branch("mc_gamma1_endz", "std::vector<double>", &_mc_gamma1_endz);
    fPiZeroTree->Branch("mc_gamma1_endt", "std::vector<double>", &_mc_gamma1_endt);
    fPiZeroTree->Branch("mc_gamma1_vpx", "std::vector<double>", &_mc_gamma1_vpx);
    fPiZeroTree->Branch("mc_gamma1_vpy", "std::vector<double>", &_mc_gamma1_vpy);
    fPiZeroTree->Branch("mc_gamma1_vpz", "std::vector<double>", &_mc_gamma1_vpz);
    fPiZeroTree->Branch("mc_gamma1_vtheta", "std::vector<double>", &_mc_gamma1_vtheta);
    fPiZeroTree->Branch("mc_gamma1_vphi", "std::vector<double>", &_mc_gamma1_vphi);
    fPiZeroTree->Branch("mc_gamma1_ve", "std::vector<double>", &_mc_gamma1_ve);
    fPiZeroTree->Branch("mc_gamma1_endpx", "std::vector<double>", &_mc_gamma1_endpx);
    fPiZeroTree->Branch("mc_gamma1_endpy", "std::vector<double>", &_mc_gamma1_endpy);
    fPiZeroTree->Branch("mc_gamma1_endpz", "std::vector<double>", &_mc_gamma1_endpz);
    fPiZeroTree->Branch("mc_gamma1_ende", "std::vector<double>", &_mc_gamma1_ende);

    fPiZeroTree->Branch("mc_gamma2_pdg", "std::vector<int>", &_mc_gamma2_pdg);
    fPiZeroTree->Branch("mc_gamma2_status", "std::vector<int>", &_mc_gamma2_status);
    fPiZeroTree->Branch("mc_gamma2_vx", "std::vector<double>", &_mc_gamma2_vx);
    fPiZeroTree->Branch("mc_gamma2_vy", "std::vector<double>", &_mc_gamma2_vy);
    fPiZeroTree->Branch("mc_gamma2_vz", "std::vector<double>", &_mc_gamma2_vz);
    fPiZeroTree->Branch("mc_gamma2_vt", "std::vector<double>", &_mc_gamma2_vt);
    fPiZeroTree->Branch("mc_gamma2_endx", "std::vector<double>", &_mc_gamma2_endx);
    fPiZeroTree->Branch("mc_gamma2_endy", "std::vector<double>", &_mc_gamma2_endy);
    fPiZeroTree->Branch("mc_gamma2_endz", "std::vector<double>", &_mc_gamma2_endz);
    fPiZeroTree->Branch("mc_gamma2_endt", "std::vector<double>", &_mc_gamma2_endt);
    fPiZeroTree->Branch("mc_gamma2_vpx", "std::vector<double>", &_mc_gamma2_vpx);
    fPiZeroTree->Branch("mc_gamma2_vpy", "std::vector<double>", &_mc_gamma2_vpy);
    fPiZeroTree->Branch("mc_gamma2_vpz", "std::vector<double>", &_mc_gamma2_vpz);
    fPiZeroTree->Branch("mc_gamma2_vtheta", "std::vector<double>", &_mc_gamma2_vtheta);
    fPiZeroTree->Branch("mc_gamma2_vphi", "std::vector<double>", &_mc_gamma2_vphi);
    fPiZeroTree->Branch("mc_gamma2_ve", "std::vector<double>", &_mc_gamma2_ve);
    fPiZeroTree->Branch("mc_gamma2_endpx", "std::vector<double>", &_mc_gamma2_endpx);
    fPiZeroTree->Branch("mc_gamma2_endpy", "std::vector<double>", &_mc_gamma2_endpy);
    fPiZeroTree->Branch("mc_gamma2_endpz", "std::vector<double>", &_mc_gamma2_endpz);
    fPiZeroTree->Branch("mc_gamma2_ende", "std::vector<double>", &_mc_gamma2_ende);

    fPiZeroTree->Branch("mc_open_angle", "std::vector<double>", &_mc_open_angle);

    fPiZeroTree->Branch("reco_gamma1_nhits", "std::vector<int>", &_reco_gamma1_nhits);
    fPiZeroTree->Branch("reco_gamma1_npfps", "std::vector<int>", &_reco_gamma1_npfps);

    fPiZeroTree->Branch("reco_gamma1_pfp_track_score", "std::vector<std::vector<double>>", &_reco_gamma1_pfp_track_score);
    fPiZeroTree->Branch("reco_gamma1_comp", "std::vector<std::vector<double>>", &_reco_gamma1_comp);
    fPiZeroTree->Branch("reco_gamma1_pur", "std::vector<std::vector<double>>", &_reco_gamma1_pur);

    fPiZeroTree->Branch("reco_gamma1_has_good_track", "std::vector<std::vector<bool>>", &_reco_gamma1_has_good_track);
    fPiZeroTree->Branch("reco_gamma1_track_start_x", "std::vector<std::vector<double>>", &_reco_gamma1_track_start_x);
    fPiZeroTree->Branch("reco_gamma1_track_start_y", "std::vector<std::vector<double>>", &_reco_gamma1_track_start_y);
    fPiZeroTree->Branch("reco_gamma1_track_start_z", "std::vector<std::vector<double>>", &_reco_gamma1_track_start_z);
    fPiZeroTree->Branch("reco_gamma1_track_end_x", "std::vector<std::vector<double>>", &_reco_gamma1_track_end_x);
    fPiZeroTree->Branch("reco_gamma1_track_end_y", "std::vector<std::vector<double>>", &_reco_gamma1_track_end_y);
    fPiZeroTree->Branch("reco_gamma1_track_end_z", "std::vector<std::vector<double>>", &_reco_gamma1_track_end_z);
    fPiZeroTree->Branch("reco_gamma1_track_dir_x", "std::vector<std::vector<double>>", &_reco_gamma1_track_dir_x);
    fPiZeroTree->Branch("reco_gamma1_track_dir_y", "std::vector<std::vector<double>>", &_reco_gamma1_track_dir_y);
    fPiZeroTree->Branch("reco_gamma1_track_dir_z", "std::vector<std::vector<double>>", &_reco_gamma1_track_dir_z);
    fPiZeroTree->Branch("reco_gamma1_track_start_theta", "std::vector<std::vector<double>>", &_reco_gamma1_track_start_theta);
    fPiZeroTree->Branch("reco_gamma1_track_start_phi", "std::vector<std::vector<double>>", &_reco_gamma1_track_start_phi);
    fPiZeroTree->Branch("reco_gamma1_track_length", "std::vector<std::vector<double>>", &_reco_gamma1_track_length);
    fPiZeroTree->Branch("reco_gamma1_track_best_plane_energy", "std::vector<std::vector<double>>", &_reco_gamma1_track_best_plane_energy);

    fPiZeroTree->Branch("reco_gamma1_has_good_shower", "std::vector<std::vector<bool>>", &_reco_gamma1_has_good_shower);
    fPiZeroTree->Branch("reco_gamma1_shower_start_x", "std::vector<std::vector<double>>", &_reco_gamma1_shower_start_x);
    fPiZeroTree->Branch("reco_gamma1_shower_start_y", "std::vector<std::vector<double>>", &_reco_gamma1_shower_start_y);
    fPiZeroTree->Branch("reco_gamma1_shower_start_z", "std::vector<std::vector<double>>", &_reco_gamma1_shower_start_z);
    fPiZeroTree->Branch("reco_gamma1_shower_dir_x", "std::vector<std::vector<double>>", &_reco_gamma1_shower_dir_x);
    fPiZeroTree->Branch("reco_gamma1_shower_dir_y", "std::vector<std::vector<double>>", &_reco_gamma1_shower_dir_y);
    fPiZeroTree->Branch("reco_gamma1_shower_dir_z", "std::vector<std::vector<double>>", &_reco_gamma1_shower_dir_z);
    fPiZeroTree->Branch("reco_gamma1_shower_start_theta", "std::vector<std::vector<double>>", &_reco_gamma1_shower_start_theta);
    fPiZeroTree->Branch("reco_gamma1_shower_start_phi", "std::vector<std::vector<double>>", &_reco_gamma1_shower_start_phi);
    fPiZeroTree->Branch("reco_gamma1_shower_length", "std::vector<std::vector<double>>", &_reco_gamma1_shower_length);
    fPiZeroTree->Branch("reco_gamma1_shower_open_angle", "std::vector<std::vector<double>>", &_reco_gamma1_shower_open_angle);
    fPiZeroTree->Branch("reco_gamma1_shower_best_plane_energy", "std::vector<std::vector<double>>", &_reco_gamma1_shower_best_plane_energy);
    fPiZeroTree->Branch("reco_gamma1_shower_best_plane_dedx", "std::vector<std::vector<double>>", &_reco_gamma1_shower_best_plane_dedx);

    fPiZeroTree->Branch("reco_gamma1_bestpfp_pfp_track_score", "std::vector<double>", &_reco_gamma1_bestpfp_pfp_track_score);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_comp", "std::vector<double>", &_reco_gamma1_bestpfp_comp);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_pur", "std::vector<double>", &_reco_gamma1_bestpfp_pur);

    fPiZeroTree->Branch("reco_gamma1_bestpfp_has_good_track", "std::vector<bool>", &_reco_gamma1_bestpfp_has_good_track);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_start_x", "std::vector<double>", &_reco_gamma1_bestpfp_track_start_x);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_start_y", "std::vector<double>", &_reco_gamma1_bestpfp_track_start_y);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_start_z", "std::vector<double>", &_reco_gamma1_bestpfp_track_start_z);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_end_x", "std::vector<double>", &_reco_gamma1_bestpfp_track_end_x);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_end_y", "std::vector<double>", &_reco_gamma1_bestpfp_track_end_y);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_end_z", "std::vector<double>", &_reco_gamma1_bestpfp_track_end_z);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_dir_x", "std::vector<double>", &_reco_gamma1_bestpfp_track_dir_x);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_dir_y", "std::vector<double>", &_reco_gamma1_bestpfp_track_dir_y);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_dir_z", "std::vector<double>", &_reco_gamma1_bestpfp_track_dir_z);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_start_theta", "std::vector<double>", &_reco_gamma1_bestpfp_track_start_theta);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_start_phi", "std::vector<double>", &_reco_gamma1_bestpfp_track_start_phi);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_length", "std::vector<double>", &_reco_gamma1_bestpfp_track_length);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_track_best_plane_energy", "std::vector<double>", &_reco_gamma1_bestpfp_track_best_plane_energy);

    fPiZeroTree->Branch("reco_gamma1_bestpfp_has_good_shower", "std::vector<bool>", &_reco_gamma1_bestpfp_has_good_shower);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_start_x", "std::vector<double>", &_reco_gamma1_bestpfp_shower_start_x);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_start_y", "std::vector<double>", &_reco_gamma1_bestpfp_shower_start_y);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_start_z", "std::vector<double>", &_reco_gamma1_bestpfp_shower_start_z);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_dir_x", "std::vector<double>", &_reco_gamma1_bestpfp_shower_dir_x);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_dir_y", "std::vector<double>", &_reco_gamma1_bestpfp_shower_dir_y);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_dir_z", "std::vector<double>", &_reco_gamma1_bestpfp_shower_dir_z);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_start_theta", "std::vector<double>", &_reco_gamma1_bestpfp_shower_start_theta);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_start_phi", "std::vector<double>", &_reco_gamma1_bestpfp_shower_start_phi);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_length", "std::vector<double>", &_reco_gamma1_bestpfp_shower_length);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_open_angle", "std::vector<double>", &_reco_gamma1_bestpfp_shower_open_angle);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_best_plane_energy", "std::vector<double>", &_reco_gamma1_bestpfp_shower_best_plane_energy);
    fPiZeroTree->Branch("reco_gamma1_bestpfp_shower_best_plane_dedx", "std::vector<double>", &_reco_gamma1_bestpfp_shower_best_plane_dedx);

    fPiZeroTree->Branch("reco_gamma2_nhits", "std::vector<int>", &_reco_gamma2_nhits);
    fPiZeroTree->Branch("reco_gamma2_npfps", "std::vector<int>", &_reco_gamma2_npfps);

    fPiZeroTree->Branch("reco_gamma2_pfp_track_score", "std::vector<std::vector<double>>", &_reco_gamma2_pfp_track_score);
    fPiZeroTree->Branch("reco_gamma2_comp", "std::vector<std::vector<double>>", &_reco_gamma2_comp);
    fPiZeroTree->Branch("reco_gamma2_pur", "std::vector<std::vector<double>>", &_reco_gamma2_pur);

    fPiZeroTree->Branch("reco_gamma2_has_good_track", "std::vector<std::vector<bool>>", &_reco_gamma2_has_good_track);
    fPiZeroTree->Branch("reco_gamma2_track_start_x", "std::vector<std::vector<double>>", &_reco_gamma2_track_start_x);
    fPiZeroTree->Branch("reco_gamma2_track_start_y", "std::vector<std::vector<double>>", &_reco_gamma2_track_start_y);
    fPiZeroTree->Branch("reco_gamma2_track_start_z", "std::vector<std::vector<double>>", &_reco_gamma2_track_start_z);
    fPiZeroTree->Branch("reco_gamma2_track_end_x", "std::vector<std::vector<double>>", &_reco_gamma2_track_end_x);
    fPiZeroTree->Branch("reco_gamma2_track_end_y", "std::vector<std::vector<double>>", &_reco_gamma2_track_end_y);
    fPiZeroTree->Branch("reco_gamma2_track_end_z", "std::vector<std::vector<double>>", &_reco_gamma2_track_end_z);
    fPiZeroTree->Branch("reco_gamma2_track_dir_x", "std::vector<std::vector<double>>", &_reco_gamma2_track_dir_x);
    fPiZeroTree->Branch("reco_gamma2_track_dir_y", "std::vector<std::vector<double>>", &_reco_gamma2_track_dir_y);
    fPiZeroTree->Branch("reco_gamma2_track_dir_z", "std::vector<std::vector<double>>", &_reco_gamma2_track_dir_z);
    fPiZeroTree->Branch("reco_gamma2_track_start_theta", "std::vector<std::vector<double>>", &_reco_gamma2_track_start_theta);
    fPiZeroTree->Branch("reco_gamma2_track_start_phi", "std::vector<std::vector<double>>", &_reco_gamma2_track_start_phi);
    fPiZeroTree->Branch("reco_gamma2_track_length", "std::vector<std::vector<double>>", &_reco_gamma2_track_length);
    fPiZeroTree->Branch("reco_gamma2_track_best_plane_energy", "std::vector<std::vector<double>>", &_reco_gamma2_track_best_plane_energy);

    fPiZeroTree->Branch("reco_gamma2_has_good_shower", "std::vector<std::vector<bool>>", &_reco_gamma2_has_good_shower);
    fPiZeroTree->Branch("reco_gamma2_shower_start_x", "std::vector<std::vector<double>>", &_reco_gamma2_shower_start_x);
    fPiZeroTree->Branch("reco_gamma2_shower_start_y", "std::vector<std::vector<double>>", &_reco_gamma2_shower_start_y);
    fPiZeroTree->Branch("reco_gamma2_shower_start_z", "std::vector<std::vector<double>>", &_reco_gamma2_shower_start_z);
    fPiZeroTree->Branch("reco_gamma2_shower_dir_x", "std::vector<std::vector<double>>", &_reco_gamma2_shower_dir_x);
    fPiZeroTree->Branch("reco_gamma2_shower_dir_y", "std::vector<std::vector<double>>", &_reco_gamma2_shower_dir_y);
    fPiZeroTree->Branch("reco_gamma2_shower_dir_z", "std::vector<std::vector<double>>", &_reco_gamma2_shower_dir_z);
    fPiZeroTree->Branch("reco_gamma2_shower_start_theta", "std::vector<std::vector<double>>", &_reco_gamma2_shower_start_theta);
    fPiZeroTree->Branch("reco_gamma2_shower_start_phi", "std::vector<std::vector<double>>", &_reco_gamma2_shower_start_phi);
    fPiZeroTree->Branch("reco_gamma2_shower_length", "std::vector<std::vector<double>>", &_reco_gamma2_shower_length);
    fPiZeroTree->Branch("reco_gamma2_shower_open_angle", "std::vector<std::vector<double>>", &_reco_gamma2_shower_open_angle);
    fPiZeroTree->Branch("reco_gamma2_shower_best_plane_energy", "std::vector<std::vector<double>>", &_reco_gamma2_shower_best_plane_energy);
    fPiZeroTree->Branch("reco_gamma2_shower_best_plane_dedx", "std::vector<std::vector<double>>", &_reco_gamma2_shower_best_plane_dedx);

    fPiZeroTree->Branch("reco_gamma2_bestpfp_pfp_track_score", "std::vector<double>", &_reco_gamma2_bestpfp_pfp_track_score);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_comp", "std::vector<double>", &_reco_gamma2_bestpfp_comp);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_pur", "std::vector<double>", &_reco_gamma2_bestpfp_pur);

    fPiZeroTree->Branch("reco_gamma2_bestpfp_has_good_track", "std::vector<bool>", &_reco_gamma2_bestpfp_has_good_track);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_start_x", "std::vector<double>", &_reco_gamma2_bestpfp_track_start_x);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_start_y", "std::vector<double>", &_reco_gamma2_bestpfp_track_start_y);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_start_z", "std::vector<double>", &_reco_gamma2_bestpfp_track_start_z);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_end_x", "std::vector<double>", &_reco_gamma2_bestpfp_track_end_x);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_end_y", "std::vector<double>", &_reco_gamma2_bestpfp_track_end_y);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_end_z", "std::vector<double>", &_reco_gamma2_bestpfp_track_end_z);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_dir_x", "std::vector<double>", &_reco_gamma2_bestpfp_track_dir_x);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_dir_y", "std::vector<double>", &_reco_gamma2_bestpfp_track_dir_y);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_dir_z", "std::vector<double>", &_reco_gamma2_bestpfp_track_dir_z);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_start_theta", "std::vector<double>", &_reco_gamma2_bestpfp_track_start_theta);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_start_phi", "std::vector<double>", &_reco_gamma2_bestpfp_track_start_phi);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_length", "std::vector<double>", &_reco_gamma2_bestpfp_track_length);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_track_best_plane_energy", "std::vector<double>", &_reco_gamma2_bestpfp_track_best_plane_energy);

    fPiZeroTree->Branch("reco_gamma2_bestpfp_has_good_shower", "std::vector<bool>", &_reco_gamma2_bestpfp_has_good_shower);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_start_x", "std::vector<double>", &_reco_gamma2_bestpfp_shower_start_x);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_start_y", "std::vector<double>", &_reco_gamma2_bestpfp_shower_start_y);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_start_z", "std::vector<double>", &_reco_gamma2_bestpfp_shower_start_z);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_dir_x", "std::vector<double>", &_reco_gamma2_bestpfp_shower_dir_x);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_dir_y", "std::vector<double>", &_reco_gamma2_bestpfp_shower_dir_y);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_dir_z", "std::vector<double>", &_reco_gamma2_bestpfp_shower_dir_z);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_start_theta", "std::vector<double>", &_reco_gamma2_bestpfp_shower_start_theta);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_start_phi", "std::vector<double>", &_reco_gamma2_bestpfp_shower_start_phi);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_length", "std::vector<double>", &_reco_gamma2_bestpfp_shower_length);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_open_angle", "std::vector<double>", &_reco_gamma2_bestpfp_shower_open_angle);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_best_plane_energy", "std::vector<double>", &_reco_gamma2_bestpfp_shower_best_plane_energy);
    fPiZeroTree->Branch("reco_gamma2_bestpfp_shower_best_plane_dedx", "std::vector<double>", &_reco_gamma2_bestpfp_shower_best_plane_dedx);
  }

void sbnd::PiZeroAnalysis::analyze(const art::Event &e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if(fDebug)
    std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  ClearMaps(e);
  SetupMaps(e);

  // Get PFParticles
  art::Handle<std::vector<recob::PFParticle>> PFParticleHandle;
  e.getByLabel(fPFParticleModuleLabel, PFParticleHandle);
  if(!PFParticleHandle.isValid()){
    std::cout << "PFParticle product " << fPFParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  AnalysePFPs(e, PFParticleHandle);

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();

  _signal = SignalEvent(e, MCTruthHandles);

  // Get MCParticles
  art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
  e.getByLabel(fMCParticleModuleLabel, MCParticleHandle);
  if(!MCParticleHandle.isValid()){
    std::cout << "MCParticle product " << fMCParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCParticle>> MCParticleVec;
  art::fill_ptr_vector(MCParticleVec, MCParticleHandle);

  // Get MCTruth to MCParticles Assns
  art::FindManyP<simb::MCTruth> MCParticlesToMCTruths(MCParticleHandle, e, fMCParticleModuleLabel);

  // Fill MCParticle variables
  AnalyseMCParticles(e, MCParticleVec, MCParticlesToMCTruths);

  // Fill the Tree
  fPiZeroTree->Fill();
}

void sbnd::PiZeroAnalysis::ClearMaps(const art::Event &e)
{
  fHitsMap.clear();
  fRecoPFPMap.clear();
}

void sbnd::PiZeroAnalysis::SetupMaps(const art::Event &e)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  art::Handle<std::vector<recob::Hit> > hitsHandle;
  e.getByLabel(fHitModuleLabel,hitsHandle);

  for(unsigned hit_i = 0; hit_i < hitsHandle->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(hitsHandle,hit_i);
    fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
  }
}

bool sbnd::PiZeroAnalysis::SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  std::vector<std::pair<bool, double>> events;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
	{
	  if(mct->Origin() != 1)
	    continue;

	  const simb::MCNeutrino mcn = mct->GetNeutrino();
	  const simb::MCParticle nu  = mcn.Nu();

	  const bool nc = mcn.CCNC() == 1;
	  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

	  unsigned pizeros = 0;

	  for(int i = 0; i < mct->NParticles(); ++i)
	    {
	      const auto mcp = mct->GetParticle(i);

	      if(mcp.PdgCode() == 111 && mcp.StatusCode() != 1)
		++pizeros;
	    }

	  const bool pizero = pizeros > 0;

	  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
	  double total_en = 0.;

	  for(auto const& mcp : MCParticleVec)
	    {
	      std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

	      for(auto const& ide : ides)
		total_en += ide->energy / 1000.;
	    }

	  events.push_back({nc && fv && pizero, total_en});
	}
    }

  if(events.size() == 0)
    return false;

  std::sort(events.begin(), events.end(),
	    [](const auto &a, const auto &b)
	    { return a.second > b.second; });

  return events.at(0).first;
}

void sbnd::PiZeroAnalysis::AnalysePFPs(const art::Event &e, const art::Handle<std::vector<recob::PFParticle>> &PFParticleHandle)
{
  std::vector<art::Ptr<recob::PFParticle>> PFParticleVec;
  art::fill_ptr_vector(PFParticleVec, PFParticleHandle);

  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackModuleLabel, trackHandle);
  if(!trackHandle.isValid()){
    std::cout << "Track product " << fTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<recob::Shower>> showerHandle;
  e.getByLabel(fShowerModuleLabel, showerHandle);
  if(!showerHandle.isValid()){
    std::cout << "Shower product " << fShowerModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::FindOneP<recob::Track>  pfpTrackAssn(PFParticleHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower> pfpShowerAssn(PFParticleHandle, e, fShowerModuleLabel);
  art::FindManyP<recob::Hit>   trackHitAssn(trackHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::Hit>   showerHitAssn(showerHandle, e, fShowerModuleLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  for(auto const& pfp : PFParticleVec)
    {
      const art::Ptr<recob::Track>  track  = pfpTrackAssn.at(pfp.key());
      const art::Ptr<recob::Shower> shower = pfpShowerAssn.at(pfp.key());

      if(track.isNonnull() && shower.isNonnull())
	{
	  const std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssn.at(track.key());
	  const int trackID      = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, trackHits, true);
	  const double trackComp = Completeness(e, trackHits, trackID);
	  const double trackPur  = Purity(e, trackHits, trackID);

	  const std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssn.at(shower.key());
	  const int showerID      = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, trackHits, true);
	  const double showerComp = Completeness(e, showerHits, showerID);
	  const double showerPur  = Purity(e, showerHits, showerID);

	  if(trackID != showerID ||
	     std::abs(trackComp - showerComp) > std::numeric_limits<double>::max() ||
	     std::abs(trackPur - showerPur) > std::numeric_limits<double>::max())
	    {
	      std::cout << "\nSomething is broken\n"
			<< "TrackID: " << trackID << "\tShowerID: " << showerID << '\n'
			<< "TrackComp: " << trackComp << "\tShowerComp: " << showerComp << '\n'
			<< "TrackPur: " << trackPur << "\tShowerPur: " << showerPur << '\n'
			<< std::endl;

	      throw std::exception();
	    }

	  fRecoPFPMap[trackID].insert(pfp);
	}
      else if(track.isNonnull())
	{
	  const std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssn.at(track.key());
          const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, trackHits, true);
	  fRecoPFPMap[trackID].insert(pfp);
	}
      else if(shower.isNonnull())
	{
	  const std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssn.at(shower.key());
          const int showerID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, true);
	  fRecoPFPMap[showerID].insert(pfp);
	}
    }
}

void sbnd::PiZeroAnalysis::AnalyseMCParticles(const art::Event &e, const std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec, const art::FindManyP<simb::MCTruth> &MCParticlesToMCTruths)
{
  const unsigned nMCParticles = MCParticleVec.size();
  unsigned nPiZeros           = 0;

  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      if(mcp->PdgCode() != 111 || mcp->StatusCode() != 1)
	continue;

      ++nPiZeros;
    }

  ResizeVectors(nPiZeros);

  unsigned pizero = 0;
  
  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      if(mcp->PdgCode() != 111 || mcp->StatusCode() != 1)
	continue;

      const auto mcts = MCParticlesToMCTruths.at(mcp.key());
      if(mcts.size() != 1)
	std::cout << "There are multiple truths" << std::endl;
      const auto mct  = mcts[0];

      _mct_origin[pizero] = mct->Origin();

      _mc_primary[pizero]   = mcp->Mother() == mct->GetNeutrino().Nu().TrackId();
      _mc_vx[pizero]        = mcp->Vx();
      _mc_vy[pizero]        = mcp->Vy();
      _mc_vz[pizero]        = mcp->Vz();
      _mc_vt[pizero]        = mcp->T();
      _mc_av[pizero]        = VolumeCheck(mcp->Position().Vect());
      _mc_fv[pizero]        = VolumeCheck(mcp->Position().Vect(), 20., 5., 10., 50.);
      _mc_uboone_fv[pizero] = VolumeCheck(mcp->Position().Vect(), 5., 5., 5., 5.);
      _mc_endx[pizero]      = mcp->EndX();
      _mc_endy[pizero]      = mcp->EndY();
      _mc_endz[pizero]      = mcp->EndZ();
      _mc_endt[pizero]      = mcp->EndT();
      _mc_vpx[pizero]       = mcp->Px();
      _mc_vpy[pizero]       = mcp->Py();
      _mc_vpz[pizero]       = mcp->Pz();
      _mc_vtheta[pizero]    = TMath::RadToDeg() * mcp->Momentum().Theta();
      _mc_vphi[pizero]      = TMath::RadToDeg() * mcp->Momentum().Phi();
      _mc_ve[pizero]        = mcp->E();
      _mc_endpx[pizero]     = mcp->EndPx();
      _mc_endpy[pizero]     = mcp->EndPy();
      _mc_endpz[pizero]     = mcp->EndPz();
      _mc_ende[pizero]      = mcp->EndE();

      if(mcp->NumberDaughters() == 2)
	{
	  const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
	  const simb::MCParticle* gamma2 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

	  if(gamma1->PdgCode() != 22 || gamma2->PdgCode() != 22)
	    {
	      _mc_digamma[pizero] = false;
	      FillGammaTruthDefault(pizero);
	      FillGammaRecoDefault(pizero);
	    }
	  else
	    {
	      _mc_digamma[pizero] = true;

	      const simb::MCParticle* gammalead    = gamma2->E() > gamma1->E() ? gamma2 : gamma1;
	      const simb::MCParticle* gammasublead = gamma2->E() > gamma1->E() ? gamma1 : gamma2;

	      FillGammaTruth(pizero, gammalead, true);
	      FillGammaTruth(pizero, gammasublead, false);

	      const TVector3 gamma1_mom = gammalead->Momentum().Vect();
	      const TVector3 gamma2_mom = gammasublead->Momentum().Vect();
	      _mc_open_angle[pizero]    = TMath::RadToDeg() * gamma1_mom.Angle(gamma2_mom);

	      FillGammaReco(pizero, e, gammalead->TrackId(), true);
	      FillGammaReco(pizero, e, gammasublead->TrackId(), false);
	    }
	}
      else
	{
	  _mc_digamma[pizero] = false;
	  FillGammaTruthDefault(pizero);
	  FillGammaRecoDefault(pizero);
	}

      FillBestReco(pizero);

      ++pizero;
    }
}

void sbnd::PiZeroAnalysis::FillGammaTruth(const int i, const simb::MCParticle* gamma, const bool lead)
{
  if(lead)
    {
      _mc_gamma1_pdg[i]    = gamma->PdgCode();
      _mc_gamma1_status[i] = gamma->StatusCode();
      _mc_gamma1_vx[i]     = gamma->Vx();
      _mc_gamma1_vy[i]     = gamma->Vy();
      _mc_gamma1_vz[i]     = gamma->Vz();
      _mc_gamma1_vt[i]     = gamma->T();
      _mc_gamma1_endx[i]   = gamma->EndX();
      _mc_gamma1_endy[i]   = gamma->EndY();
      _mc_gamma1_endz[i]   = gamma->EndZ();
      _mc_gamma1_endt[i]   = gamma->EndT();
      _mc_gamma1_vpx[i]    = gamma->Px();
      _mc_gamma1_vpy[i]    = gamma->Py();
      _mc_gamma1_vpz[i]    = gamma->Pz();
      _mc_gamma1_vtheta[i] = TMath::RadToDeg() * gamma->Momentum().Theta();
      _mc_gamma1_vphi[i]   = TMath::RadToDeg() * gamma->Momentum().Phi();
      _mc_gamma1_ve[i]     = gamma->E();
      _mc_gamma1_endpx[i]  = gamma->EndPx();
      _mc_gamma1_endpy[i]  = gamma->EndPy();
      _mc_gamma1_endpz[i]  = gamma->EndPz();
      _mc_gamma1_ende[i]   = gamma->EndE();
    }
  else
    {
      _mc_gamma2_pdg[i]    = gamma->PdgCode();
      _mc_gamma2_status[i] = gamma->StatusCode();
      _mc_gamma2_vx[i]     = gamma->Vx();
      _mc_gamma2_vy[i]     = gamma->Vy();
      _mc_gamma2_vz[i]     = gamma->Vz();
      _mc_gamma2_vt[i]     = gamma->T();
      _mc_gamma2_endx[i]   = gamma->EndX();
      _mc_gamma2_endy[i]   = gamma->EndY();
      _mc_gamma2_endz[i]   = gamma->EndZ();
      _mc_gamma2_endt[i]   = gamma->EndT();
      _mc_gamma2_vpx[i]    = gamma->Px();
      _mc_gamma2_vpy[i]    = gamma->Py();
      _mc_gamma2_vpz[i]    = gamma->Pz();
      _mc_gamma2_vtheta[i] = TMath::RadToDeg() * gamma->Momentum().Theta();
      _mc_gamma2_vphi[i]   = TMath::RadToDeg() * gamma->Momentum().Phi();
      _mc_gamma2_ve[i]     = gamma->E();
      _mc_gamma2_endpx[i]  = gamma->EndPx();
      _mc_gamma2_endpy[i]  = gamma->EndPy();
      _mc_gamma2_endpz[i]  = gamma->EndPz();
      _mc_gamma2_ende[i]   = gamma->EndE();
    }
}

void sbnd::PiZeroAnalysis::FillGammaTruthDefault(const int i)
{
  _mc_gamma1_pdg[i]    = def_int;
  _mc_gamma1_status[i] = def_int;
  _mc_gamma1_vx[i]     = def_double;
  _mc_gamma1_vy[i]     = def_double;
  _mc_gamma1_vz[i]     = def_double;
  _mc_gamma1_vt[i]     = def_double;
  _mc_gamma1_endx[i]   = def_double;
  _mc_gamma1_endy[i]   = def_double;
  _mc_gamma1_endz[i]   = def_double;
  _mc_gamma1_endt[i]   = def_double;
  _mc_gamma1_vpx[i]    = def_double;
  _mc_gamma1_vpy[i]    = def_double;
  _mc_gamma1_vpz[i]    = def_double;
  _mc_gamma1_ve[i]     = def_double;
  _mc_gamma1_vtheta[i] = def_double;
  _mc_gamma1_vphi[i]   = def_double;
  _mc_gamma1_endpx[i]  = def_double;
  _mc_gamma1_endpy[i]  = def_double;
  _mc_gamma1_endpz[i]  = def_double;
  _mc_gamma1_ende[i]   = def_double;

  _mc_gamma2_pdg[i]    = def_int;
  _mc_gamma2_status[i] = def_int;
  _mc_gamma2_vx[i]     = def_double;
  _mc_gamma2_vy[i]     = def_double;
  _mc_gamma2_vz[i]     = def_double;
  _mc_gamma2_vt[i]     = def_double;
  _mc_gamma2_endx[i]   = def_double;
  _mc_gamma2_endy[i]   = def_double;
  _mc_gamma2_endz[i]   = def_double;
  _mc_gamma2_endt[i]   = def_double;
  _mc_gamma2_vpx[i]    = def_double;
  _mc_gamma2_vpy[i]    = def_double;
  _mc_gamma2_vpz[i]    = def_double;
  _mc_gamma2_ve[i]     = def_double;
  _mc_gamma2_vtheta[i] = def_double;
  _mc_gamma2_vphi[i]   = def_double;
  _mc_gamma2_endpx[i]  = def_double;
  _mc_gamma2_endpy[i]  = def_double;
  _mc_gamma2_endpz[i]  = def_double;
  _mc_gamma2_ende[i]   = def_double;
}

void sbnd::PiZeroAnalysis::FillGammaReco(const int i, const art::Event &e, const int trueID, const bool lead)
{
  art::Handle<std::vector<recob::PFParticle>> PFParticleHandle;
  e.getByLabel(fPFParticleModuleLabel, PFParticleHandle);
  if(!PFParticleHandle.isValid()){
    std::cout << "PFParticle product " << fPFParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackModuleLabel, trackHandle);
  if(!trackHandle.isValid()){
    std::cout << "Track product " << fTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<recob::Shower>> showerHandle;
  e.getByLabel(fShowerModuleLabel, showerHandle);
  if(!showerHandle.isValid()){
    std::cout << "Shower product " << fShowerModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  const std::set<art::Ptr<recob::PFParticle>> pfps = fRecoPFPMap.count(trueID) == 0 ? std::set<art::Ptr<recob::PFParticle>>() : fRecoPFPMap.at(trueID);
  const int nPFPs = fRecoPFPMap.count(trueID) == 0 ? 0 : pfps.size();

  art::FindOneP<recob::Track>                      pfpTrackAssn(PFParticleHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpShowerAssn(PFParticleHandle, e, fShowerModuleLabel);
  art::FindManyP<recob::Hit>                       trackHitAssn(trackHandle, e, fTrackModuleLabel);
  art::FindManyP<anab::Calorimetry>                trackCaloAssn(trackHandle, e, fTrackCalorimetryModuleLabel);
  art::FindManyP<recob::Hit>                       showerHitAssn(showerHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpsToMetadata(PFParticleHandle, e, fPFParticleModuleLabel);

  ResizeRecoSubVectors(i, nPFPs, lead);

  if(lead)
    {
      _reco_gamma1_nhits[i] = fHitsMap[trueID];
      _reco_gamma1_npfps[i] = nPFPs;

      int pfp_i = 0;

      for(auto const& pfp : pfps)
	{
	  const auto meta  = pfpsToMetadata.at(pfp.key());
	  const auto props = meta->GetPropertiesMap();
	  const auto trackscore = props.find("TrackScore");

	  if(trackscore != props.end())
	    _reco_gamma1_pfp_track_score[i][pfp_i] = trackscore->second;
	  else
	    _reco_gamma1_pfp_track_score[i][pfp_i] = def_double;

	  const art::Ptr<recob::Track> track = pfpTrackAssn.at(pfp.key());

	  if(track.isNonnull())
	    {
	      const std::vector<art::Ptr<recob::Hit>> trackHits         = trackHitAssn.at(track.key());
	      const std::vector<art::Ptr<anab::Calorimetry>> trackCalos = trackCaloAssn.at(track.key());

	      _reco_gamma1_comp[i][pfp_i] = Completeness(e, trackHits, trueID);
	      _reco_gamma1_pur[i][pfp_i]  = Purity(e, trackHits, trueID);

	      _reco_gamma1_has_good_track[i][pfp_i]    = true;
	      _reco_gamma1_track_start_x[i][pfp_i]     = track->Start().X();
	      _reco_gamma1_track_start_y[i][pfp_i]     = track->Start().Y();
	      _reco_gamma1_track_start_z[i][pfp_i]     = track->Start().Z();
	      _reco_gamma1_track_end_x[i][pfp_i]       = track->End().X();
	      _reco_gamma1_track_end_y[i][pfp_i]       = track->End().Y();
	      _reco_gamma1_track_end_z[i][pfp_i]       = track->End().Z();
	      _reco_gamma1_track_dir_x[i][pfp_i]       = track->StartDirection().X();
	      _reco_gamma1_track_dir_y[i][pfp_i]       = track->StartDirection().Y();
	      _reco_gamma1_track_dir_z[i][pfp_i]       = track->StartDirection().Z();
	      _reco_gamma1_track_start_theta[i][pfp_i] = TMath::RadToDeg() * track->StartDirection().Theta();
	      _reco_gamma1_track_start_phi[i][pfp_i]   = TMath::RadToDeg() * track->StartDirection().Phi();
	      _reco_gamma1_track_length[i][pfp_i]      = track->Length();

	      if(trackCalos.size() == 3)
		{
		  std::vector<std::pair<int, double>> pitches;
		  int calo_i = 0;

		  for(auto const& calo : trackCalos)
		    {
		      const std::vector<float> pitchVec = calo->TrkPitchVec();
		      const double avePitch             = std::reduce(pitchVec.begin(), pitchVec.end()) / pitchVec.size();
		      pitches.push_back({calo_i, avePitch});
		      ++calo_i;
		    }

		  std::sort(pitches.begin(), pitches.end(),
			    [](const auto &a, const auto &b)
			    { return a.second < b.second; });

		  const int bestplane = pitches.size() == 0 ? 0 : pitches[0].first;

		  _reco_gamma1_track_best_plane_energy[i][pfp_i] = trackCalos[bestplane]->KineticEnergy();
		}
	    }
	  else
	    {
	      _reco_gamma1_comp[i][pfp_i] = def_double;
	      _reco_gamma1_pur[i][pfp_i]  = def_double;

	      _reco_gamma1_has_good_track[i][pfp_i]          = false;
	      _reco_gamma1_track_start_x[i][pfp_i]           = def_double;
	      _reco_gamma1_track_start_y[i][pfp_i]           = def_double;
	      _reco_gamma1_track_start_z[i][pfp_i]           = def_double;
	      _reco_gamma1_track_end_x[i][pfp_i]             = def_double;
	      _reco_gamma1_track_end_y[i][pfp_i]             = def_double;
	      _reco_gamma1_track_end_z[i][pfp_i]             = def_double;
	      _reco_gamma1_track_dir_x[i][pfp_i]             = def_double;
	      _reco_gamma1_track_dir_y[i][pfp_i]             = def_double;
	      _reco_gamma1_track_dir_z[i][pfp_i]             = def_double;
	      _reco_gamma1_track_start_theta[i][pfp_i]       = def_double;
	      _reco_gamma1_track_start_phi[i][pfp_i]         = def_double;
	      _reco_gamma1_track_best_plane_energy[i][pfp_i] = def_double;
	      _reco_gamma1_track_length[i][pfp_i]            = def_double;
	    }

	  const art::Ptr<recob::Shower> shower = pfpShowerAssn.at(pfp.key());

	  if(shower.isNonnull())
	    {
	      const std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssn.at(shower.key());

	      _reco_gamma1_comp[i][pfp_i] = Completeness(e, showerHits, trueID);
	      _reco_gamma1_pur[i][pfp_i]  = Purity(e, showerHits, trueID);

	      _reco_gamma1_has_good_shower[i][pfp_i]          = true;
	      _reco_gamma1_shower_start_x[i][pfp_i]           = shower->ShowerStart().X();
	      _reco_gamma1_shower_start_y[i][pfp_i]           = shower->ShowerStart().Y();
	      _reco_gamma1_shower_start_z[i][pfp_i]           = shower->ShowerStart().Z();
	      _reco_gamma1_shower_dir_x[i][pfp_i]             = shower->Direction().X();
	      _reco_gamma1_shower_dir_y[i][pfp_i]             = shower->Direction().Y();
	      _reco_gamma1_shower_dir_z[i][pfp_i]             = shower->Direction().Z();
	      _reco_gamma1_shower_start_theta[i][pfp_i]       = TMath::RadToDeg() * shower->Direction().Theta();
	      _reco_gamma1_shower_start_phi[i][pfp_i]         = TMath::RadToDeg() * shower->Direction().Phi();
	      _reco_gamma1_shower_length[i][pfp_i]            = shower->Length();
	      _reco_gamma1_shower_open_angle[i][pfp_i]        = shower->OpenAngle();
	      _reco_gamma1_shower_best_plane_energy[i][pfp_i] = shower->Energy()[shower->best_plane()];
	      _reco_gamma1_shower_best_plane_dedx[i][pfp_i]   = shower->dEdx()[shower->best_plane()];
	    }
	  else
	    {
	      if(track.isNull())
		{
		  _reco_gamma1_comp[i][pfp_i] = def_double;
		  _reco_gamma1_pur[i][pfp_i]  = def_double;
		}

	      _reco_gamma1_has_good_shower[i][pfp_i]          = false;
	      _reco_gamma1_shower_start_x[i][pfp_i]           = def_double;
	      _reco_gamma1_shower_start_y[i][pfp_i]           = def_double;
	      _reco_gamma1_shower_start_z[i][pfp_i]           = def_double;
	      _reco_gamma1_shower_dir_x[i][pfp_i]             = def_double;
	      _reco_gamma1_shower_dir_y[i][pfp_i]             = def_double;
	      _reco_gamma1_shower_dir_z[i][pfp_i]             = def_double;
	      _reco_gamma1_shower_start_theta[i][pfp_i]       = def_double;
	      _reco_gamma1_shower_start_phi[i][pfp_i]         = def_double;
	      _reco_gamma1_shower_length[i][pfp_i]            = def_double;
	      _reco_gamma1_shower_open_angle[i][pfp_i]        = def_double;
	      _reco_gamma1_shower_best_plane_energy[i][pfp_i] = def_double;
	      _reco_gamma1_shower_best_plane_dedx[i][pfp_i]   = def_double;
	    }
	  ++pfp_i;
	}
    }
  else
    {
      _reco_gamma2_nhits[i] = fHitsMap[trueID];
      _reco_gamma2_npfps[i] = nPFPs;

      int pfp_i = 0;

      for(auto const& pfp : pfps)
	{
	  const auto meta  = pfpsToMetadata.at(pfp.key());
	  const auto props = meta->GetPropertiesMap();
	  const auto trackscore = props.find("TrackScore");

	  if(trackscore != props.end())
	    _reco_gamma2_pfp_track_score[i][pfp_i] = trackscore->second;
	  else
	    _reco_gamma2_pfp_track_score[i][pfp_i] = def_double;

	  const art::Ptr<recob::Track> track  = pfpTrackAssn.at(pfp.key());

	  if(track.isNonnull())
	    {
	      const std::vector<art::Ptr<recob::Hit>> trackHits         = trackHitAssn.at(track.key());
	      const std::vector<art::Ptr<anab::Calorimetry>> trackCalos = trackCaloAssn.at(track.key());

	      _reco_gamma2_comp[i][pfp_i] = Completeness(e, trackHits, trueID);
	      _reco_gamma2_pur[i][pfp_i]  = Purity(e, trackHits, trueID);

	      _reco_gamma2_has_good_track[i][pfp_i]    = true;
	      _reco_gamma2_track_start_x[i][pfp_i]     = track->Start().X();
	      _reco_gamma2_track_start_y[i][pfp_i]     = track->Start().Y();
	      _reco_gamma2_track_start_z[i][pfp_i]     = track->Start().Z();
	      _reco_gamma2_track_end_x[i][pfp_i]       = track->End().X();
	      _reco_gamma2_track_end_y[i][pfp_i]       = track->End().Y();
	      _reco_gamma2_track_end_z[i][pfp_i]       = track->End().Z();
	      _reco_gamma2_track_dir_x[i][pfp_i]       = track->StartDirection().X();
	      _reco_gamma2_track_dir_y[i][pfp_i]       = track->StartDirection().Y();
	      _reco_gamma2_track_dir_z[i][pfp_i]       = track->StartDirection().Z();
	      _reco_gamma2_track_start_theta[i][pfp_i] = TMath::RadToDeg() * track->StartDirection().Theta();
	      _reco_gamma2_track_start_phi[i][pfp_i]   = TMath::RadToDeg() * track->StartDirection().Phi();
	      _reco_gamma2_track_length[i][pfp_i]      = track->Length();

	      if(trackCalos.size() == 3)
		{
		  std::vector<std::pair<int, double>> pitches;
		  int calo_i = 0;

		  for(auto const& calo : trackCalos)
		    {
		      const std::vector<float> pitchVec = calo->TrkPitchVec();
		      const double avePitch             = std::reduce(pitchVec.begin(), pitchVec.end()) / pitchVec.size();
		      pitches.push_back({calo_i, avePitch});
		      ++calo_i;
		    }

		  std::sort(pitches.begin(), pitches.end(),
			    [](const auto &a, const auto &b)
			    { return a.second < b.second; });

		  const int bestplane = pitches.size() == 0 ? 0 : pitches[0].first;

		  _reco_gamma2_track_best_plane_energy[i][pfp_i] = trackCalos[bestplane]->KineticEnergy();
		}
	    }
	  else
	    {
	      _reco_gamma2_comp[i][pfp_i] = def_double;
	      _reco_gamma2_pur[i][pfp_i]  = def_double;

	      _reco_gamma2_has_good_track[i][pfp_i]          = false;
	      _reco_gamma2_track_start_x[i][pfp_i]           = def_double;
	      _reco_gamma2_track_start_y[i][pfp_i]           = def_double;
	      _reco_gamma2_track_start_z[i][pfp_i]           = def_double;
	      _reco_gamma2_track_end_x[i][pfp_i]             = def_double;
	      _reco_gamma2_track_end_y[i][pfp_i]             = def_double;
	      _reco_gamma2_track_end_z[i][pfp_i]             = def_double;
	      _reco_gamma2_track_dir_x[i][pfp_i]             = def_double;
	      _reco_gamma2_track_dir_y[i][pfp_i]             = def_double;
	      _reco_gamma2_track_dir_z[i][pfp_i]             = def_double;
	      _reco_gamma2_track_start_theta[i][pfp_i]       = def_double;
	      _reco_gamma2_track_start_phi[i][pfp_i]         = def_double;
	      _reco_gamma2_track_best_plane_energy[i][pfp_i] = def_double;
	      _reco_gamma2_track_length[i][pfp_i]            = def_double;
	    }

	  const art::Ptr<recob::Shower> shower = pfpShowerAssn.at(pfp.key());

	  if(shower.isNonnull())
	    {
	      const std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssn.at(shower.key());

	      _reco_gamma2_comp[i][pfp_i] = Completeness(e, showerHits, trueID);
	      _reco_gamma2_pur[i][pfp_i]  = Purity(e, showerHits, trueID);

	      _reco_gamma2_has_good_shower[i][pfp_i]          = true;
	      _reco_gamma2_shower_start_x[i][pfp_i]           = shower->ShowerStart().X();
	      _reco_gamma2_shower_start_y[i][pfp_i]           = shower->ShowerStart().Y();
	      _reco_gamma2_shower_start_z[i][pfp_i]           = shower->ShowerStart().Z();
	      _reco_gamma2_shower_dir_x[i][pfp_i]             = shower->Direction().X();
	      _reco_gamma2_shower_dir_y[i][pfp_i]             = shower->Direction().Y();
	      _reco_gamma2_shower_dir_z[i][pfp_i]             = shower->Direction().Z();
	      _reco_gamma2_shower_start_theta[i][pfp_i]       = TMath::RadToDeg() * shower->Direction().Theta();
	      _reco_gamma2_shower_start_phi[i][pfp_i]         = TMath::RadToDeg() * shower->Direction().Phi();
	      _reco_gamma2_shower_length[i][pfp_i]            = shower->Length();
	      _reco_gamma2_shower_open_angle[i][pfp_i]        = shower->OpenAngle();
	      _reco_gamma2_shower_best_plane_energy[i][pfp_i] = shower->Energy()[shower->best_plane()];
	      _reco_gamma2_shower_best_plane_dedx[i][pfp_i]   = shower->dEdx()[shower->best_plane()];
	    }
	  else
	    {
	      if(track.isNull())
		{
		  _reco_gamma2_comp[i][pfp_i] = def_double;
		  _reco_gamma2_pur[i][pfp_i]  = def_double;
		}

	      _reco_gamma2_has_good_shower[i][pfp_i]          = false;
	      _reco_gamma2_shower_start_x[i][pfp_i]           = def_double;
	      _reco_gamma2_shower_start_y[i][pfp_i]           = def_double;
	      _reco_gamma2_shower_start_z[i][pfp_i]           = def_double;
	      _reco_gamma2_shower_dir_x[i][pfp_i]             = def_double;
	      _reco_gamma2_shower_dir_y[i][pfp_i]             = def_double;
	      _reco_gamma2_shower_dir_z[i][pfp_i]             = def_double;
	      _reco_gamma2_shower_start_theta[i][pfp_i]       = def_double;
	      _reco_gamma2_shower_start_phi[i][pfp_i]         = def_double;
	      _reco_gamma2_shower_length[i][pfp_i]            = def_double;
	      _reco_gamma2_shower_open_angle[i][pfp_i]        = def_double;
	      _reco_gamma2_shower_best_plane_energy[i][pfp_i] = def_double;
	      _reco_gamma2_shower_best_plane_dedx[i][pfp_i]   = def_double;
	    }
	  ++pfp_i;
	}
    }
}

void sbnd::PiZeroAnalysis::FillGammaRecoDefault(const int i)
{
  _reco_gamma1_nhits[i] = 0;
  _reco_gamma1_npfps[i] = 0;
  _reco_gamma2_nhits[i] = 0;
  _reco_gamma2_npfps[i] = 0;

  _reco_gamma1_bestpfp_pfp_track_score[i] = def_double;
  _reco_gamma1_bestpfp_comp[i]            = def_double;
  _reco_gamma1_bestpfp_pur[i]             = def_double;

  _reco_gamma1_bestpfp_has_good_track[i]          = false;
  _reco_gamma1_bestpfp_track_start_x[i]           = def_double;
  _reco_gamma1_bestpfp_track_start_y[i]           = def_double;
  _reco_gamma1_bestpfp_track_start_z[i]           = def_double;
  _reco_gamma1_bestpfp_track_end_x[i]             = def_double;
  _reco_gamma1_bestpfp_track_end_y[i]             = def_double;
  _reco_gamma1_bestpfp_track_end_z[i]             = def_double;
  _reco_gamma1_bestpfp_track_dir_x[i]             = def_double;
  _reco_gamma1_bestpfp_track_dir_y[i]             = def_double;
  _reco_gamma1_bestpfp_track_dir_z[i]             = def_double;
  _reco_gamma1_bestpfp_track_start_theta[i]       = def_double;
  _reco_gamma1_bestpfp_track_start_phi[i]         = def_double;
  _reco_gamma1_bestpfp_track_length[i]            = def_double;
  _reco_gamma1_bestpfp_track_best_plane_energy[i] = def_double;

  _reco_gamma1_bestpfp_has_good_shower[i]          = false;
  _reco_gamma1_bestpfp_shower_start_x[i]           = def_double;
  _reco_gamma1_bestpfp_shower_start_y[i]           = def_double;
  _reco_gamma1_bestpfp_shower_start_z[i]           = def_double;
  _reco_gamma1_bestpfp_shower_dir_x[i]             = def_double;
  _reco_gamma1_bestpfp_shower_dir_y[i]             = def_double;
  _reco_gamma1_bestpfp_shower_dir_z[i]             = def_double;
  _reco_gamma1_bestpfp_shower_start_theta[i]       = def_double;
  _reco_gamma1_bestpfp_shower_start_phi[i]         = def_double;
  _reco_gamma1_bestpfp_shower_length[i]            = def_double;
  _reco_gamma1_bestpfp_shower_open_angle[i]        = def_double;
  _reco_gamma1_bestpfp_shower_best_plane_energy[i] = def_double;
  _reco_gamma1_bestpfp_shower_best_plane_dedx[i]   = def_double;

  _reco_gamma2_bestpfp_pfp_track_score[i] = def_double;
  _reco_gamma2_bestpfp_comp[i]            = def_double;
  _reco_gamma2_bestpfp_pur[i]             = def_double;

  _reco_gamma2_bestpfp_has_good_track[i]          = false;
  _reco_gamma2_bestpfp_track_start_x[i]           = def_double;
  _reco_gamma2_bestpfp_track_start_y[i]           = def_double;
  _reco_gamma2_bestpfp_track_start_z[i]           = def_double;
  _reco_gamma2_bestpfp_track_end_x[i]             = def_double;
  _reco_gamma2_bestpfp_track_end_y[i]             = def_double;
  _reco_gamma2_bestpfp_track_end_z[i]             = def_double;
  _reco_gamma2_bestpfp_track_dir_x[i]             = def_double;
  _reco_gamma2_bestpfp_track_dir_y[i]             = def_double;
  _reco_gamma2_bestpfp_track_dir_z[i]             = def_double;
  _reco_gamma2_bestpfp_track_start_theta[i]       = def_double;
  _reco_gamma2_bestpfp_track_start_phi[i]         = def_double;
  _reco_gamma2_bestpfp_track_length[i]            = def_double;
  _reco_gamma2_bestpfp_track_best_plane_energy[i] = def_double;

  _reco_gamma2_bestpfp_has_good_shower[i]          = false;
  _reco_gamma2_bestpfp_shower_start_x[i]           = def_double;
  _reco_gamma2_bestpfp_shower_start_y[i]           = def_double;
  _reco_gamma2_bestpfp_shower_start_z[i]           = def_double;
  _reco_gamma2_bestpfp_shower_dir_x[i]             = def_double;
  _reco_gamma2_bestpfp_shower_dir_y[i]             = def_double;
  _reco_gamma2_bestpfp_shower_dir_z[i]             = def_double;
  _reco_gamma2_bestpfp_shower_start_theta[i]       = def_double;
  _reco_gamma2_bestpfp_shower_start_phi[i]         = def_double;
  _reco_gamma2_bestpfp_shower_length[i]            = def_double;
  _reco_gamma2_bestpfp_shower_open_angle[i]        = def_double;
  _reco_gamma2_bestpfp_shower_best_plane_energy[i] = def_double;
  _reco_gamma2_bestpfp_shower_best_plane_dedx[i]   = def_double;

  ResizeRecoSubVectors(i, 0, true);
  ResizeRecoSubVectors(i, 0, false);
}

void sbnd::PiZeroAnalysis::FillBestReco(const int i)
{
  if(_reco_gamma1_comp[i].size() > 0)
    {
      int gamma1_best = def_int;
      double best_comp = def_double;

      for(unsigned ii = 0; ii < _reco_gamma1_comp[i].size(); ++ii)
	{
	  if(_reco_gamma1_comp[i][ii] > best_comp)
	    {
	      best_comp   = _reco_gamma1_comp[i][ii];
	      gamma1_best = ii;
	    }
	}

      _reco_gamma1_bestpfp_pfp_track_score[i] = _reco_gamma1_pfp_track_score[i][gamma1_best];
      _reco_gamma1_bestpfp_comp[i]            = _reco_gamma1_comp[i][gamma1_best];
      _reco_gamma1_bestpfp_pur[i]             = _reco_gamma1_pur[i][gamma1_best];

      _reco_gamma1_bestpfp_has_good_track[i]          = _reco_gamma1_has_good_track[i][gamma1_best];
      _reco_gamma1_bestpfp_track_start_x[i]           = _reco_gamma1_track_start_x[i][gamma1_best];
      _reco_gamma1_bestpfp_track_start_y[i]           = _reco_gamma1_track_start_y[i][gamma1_best];
      _reco_gamma1_bestpfp_track_start_z[i]           = _reco_gamma1_track_start_z[i][gamma1_best];
      _reco_gamma1_bestpfp_track_end_x[i]             = _reco_gamma1_track_end_x[i][gamma1_best];
      _reco_gamma1_bestpfp_track_end_y[i]             = _reco_gamma1_track_end_y[i][gamma1_best];
      _reco_gamma1_bestpfp_track_end_z[i]             = _reco_gamma1_track_end_z[i][gamma1_best];
      _reco_gamma1_bestpfp_track_dir_x[i]             = _reco_gamma1_track_dir_x[i][gamma1_best];
      _reco_gamma1_bestpfp_track_dir_y[i]             = _reco_gamma1_track_dir_y[i][gamma1_best];
      _reco_gamma1_bestpfp_track_dir_z[i]             = _reco_gamma1_track_dir_z[i][gamma1_best];
      _reco_gamma1_bestpfp_track_start_theta[i]       = _reco_gamma1_track_start_theta[i][gamma1_best];
      _reco_gamma1_bestpfp_track_start_phi[i]         = _reco_gamma1_track_start_phi[i][gamma1_best];
      _reco_gamma1_bestpfp_track_length[i]            = _reco_gamma1_track_length[i][gamma1_best];
      _reco_gamma1_bestpfp_track_best_plane_energy[i] = _reco_gamma1_track_best_plane_energy[i][gamma1_best];

      _reco_gamma1_bestpfp_has_good_shower[i]          = _reco_gamma1_has_good_shower[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_start_x[i]           = _reco_gamma1_shower_start_x[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_start_y[i]           = _reco_gamma1_shower_start_y[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_start_z[i]           = _reco_gamma1_shower_start_z[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_dir_x[i]             = _reco_gamma1_shower_dir_x[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_dir_y[i]             = _reco_gamma1_shower_dir_y[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_dir_z[i]             = _reco_gamma1_shower_dir_z[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_start_theta[i]       = _reco_gamma1_shower_start_theta[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_start_phi[i]         = _reco_gamma1_shower_start_phi[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_length[i]            = _reco_gamma1_shower_length[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_open_angle[i]        = _reco_gamma1_shower_open_angle[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_best_plane_energy[i] = _reco_gamma1_shower_best_plane_energy[i][gamma1_best];
      _reco_gamma1_bestpfp_shower_best_plane_dedx[i]   = _reco_gamma1_shower_best_plane_dedx[i][gamma1_best];
    }

  if(_reco_gamma2_comp[i].size() > 0)
    {
      int gamma2_best = def_int;
      double best_comp = def_double;

      for(unsigned ii = 0; ii < _reco_gamma2_comp[i].size(); ++ii)
	{
	  if(_reco_gamma2_comp[i][ii] > best_comp)
	    {
	      best_comp   = _reco_gamma2_comp[i][ii];
	      gamma2_best = ii;
	    }
	}

      _reco_gamma2_bestpfp_pfp_track_score[i] = _reco_gamma2_pfp_track_score[i][gamma2_best];
      _reco_gamma2_bestpfp_comp[i]            = _reco_gamma2_comp[i][gamma2_best];
      _reco_gamma2_bestpfp_pur[i]             = _reco_gamma2_pur[i][gamma2_best];

      _reco_gamma2_bestpfp_has_good_track[i]          = _reco_gamma2_has_good_track[i][gamma2_best];
      _reco_gamma2_bestpfp_track_start_x[i]           = _reco_gamma2_track_start_x[i][gamma2_best];
      _reco_gamma2_bestpfp_track_start_y[i]           = _reco_gamma2_track_start_y[i][gamma2_best];
      _reco_gamma2_bestpfp_track_start_z[i]           = _reco_gamma2_track_start_z[i][gamma2_best];
      _reco_gamma2_bestpfp_track_end_x[i]             = _reco_gamma2_track_end_x[i][gamma2_best];
      _reco_gamma2_bestpfp_track_end_y[i]             = _reco_gamma2_track_end_y[i][gamma2_best];
      _reco_gamma2_bestpfp_track_end_z[i]             = _reco_gamma2_track_end_z[i][gamma2_best];
      _reco_gamma2_bestpfp_track_dir_x[i]             = _reco_gamma2_track_dir_x[i][gamma2_best];
      _reco_gamma2_bestpfp_track_dir_y[i]             = _reco_gamma2_track_dir_y[i][gamma2_best];
      _reco_gamma2_bestpfp_track_dir_z[i]             = _reco_gamma2_track_dir_z[i][gamma2_best];
      _reco_gamma2_bestpfp_track_start_theta[i]       = _reco_gamma2_track_start_theta[i][gamma2_best];
      _reco_gamma2_bestpfp_track_start_phi[i]         = _reco_gamma2_track_start_phi[i][gamma2_best];
      _reco_gamma2_bestpfp_track_length[i]            = _reco_gamma2_track_length[i][gamma2_best];
      _reco_gamma2_bestpfp_track_best_plane_energy[i] = _reco_gamma2_track_best_plane_energy[i][gamma2_best];

      _reco_gamma2_bestpfp_has_good_shower[i]          = _reco_gamma2_has_good_shower[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_start_x[i]           = _reco_gamma2_shower_start_x[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_start_y[i]           = _reco_gamma2_shower_start_y[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_start_z[i]           = _reco_gamma2_shower_start_z[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_dir_x[i]             = _reco_gamma2_shower_dir_x[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_dir_y[i]             = _reco_gamma2_shower_dir_y[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_dir_z[i]             = _reco_gamma2_shower_dir_z[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_start_theta[i]       = _reco_gamma2_shower_start_theta[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_start_phi[i]         = _reco_gamma2_shower_start_phi[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_length[i]            = _reco_gamma2_shower_length[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_open_angle[i]        = _reco_gamma2_shower_open_angle[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_best_plane_energy[i] = _reco_gamma2_shower_best_plane_energy[i][gamma2_best];
      _reco_gamma2_bestpfp_shower_best_plane_dedx[i]   = _reco_gamma2_shower_best_plane_dedx[i][gamma2_best];
    }
}

void sbnd::PiZeroAnalysis::ResizeVectors(const int size)
{
  _mct_origin.resize(size);
  _mc_primary.resize(size);
  _mc_vx.resize(size);
  _mc_vy.resize(size);
  _mc_vz.resize(size);
  _mc_vt.resize(size);
  _mc_av.resize(size);
  _mc_fv.resize(size);
  _mc_uboone_fv.resize(size);
  _mc_endx.resize(size);
  _mc_endy.resize(size);
  _mc_endz.resize(size);
  _mc_endt.resize(size);
  _mc_vpx.resize(size);
  _mc_vpy.resize(size);
  _mc_vpz.resize(size);
  _mc_vtheta.resize(size);
  _mc_vphi.resize(size);
  _mc_ve.resize(size);
  _mc_endpx.resize(size);
  _mc_endpy.resize(size);
  _mc_endpz.resize(size);
  _mc_ende.resize(size);

  _mc_digamma.resize(size);

  _mc_gamma1_pdg.resize(size);
  _mc_gamma1_status.resize(size);
  _mc_gamma1_vx.resize(size);
  _mc_gamma1_vy.resize(size);
  _mc_gamma1_vz.resize(size);
  _mc_gamma1_vt.resize(size);
  _mc_gamma1_endx.resize(size);
  _mc_gamma1_endy.resize(size);
  _mc_gamma1_endz.resize(size);
  _mc_gamma1_endt.resize(size);
  _mc_gamma1_vpx.resize(size);
  _mc_gamma1_vpy.resize(size);
  _mc_gamma1_vpz.resize(size);
  _mc_gamma1_vtheta.resize(size);
  _mc_gamma1_vphi.resize(size);
  _mc_gamma1_ve.resize(size);
  _mc_gamma1_endpx.resize(size);
  _mc_gamma1_endpy.resize(size);
  _mc_gamma1_endpz.resize(size);
  _mc_gamma1_ende.resize(size);

  _mc_gamma2_pdg.resize(size);
  _mc_gamma2_status.resize(size);
  _mc_gamma2_vx.resize(size);
  _mc_gamma2_vy.resize(size);
  _mc_gamma2_vz.resize(size);
  _mc_gamma2_vt.resize(size);
  _mc_gamma2_endx.resize(size);
  _mc_gamma2_endy.resize(size);
  _mc_gamma2_endz.resize(size);
  _mc_gamma2_endt.resize(size);
  _mc_gamma2_vpx.resize(size);
  _mc_gamma2_vpy.resize(size);
  _mc_gamma2_vpz.resize(size);
  _mc_gamma2_vtheta.resize(size);
  _mc_gamma2_vphi.resize(size);
  _mc_gamma2_ve.resize(size);
  _mc_gamma2_endpx.resize(size);
  _mc_gamma2_endpy.resize(size);
  _mc_gamma2_endpz.resize(size);
  _mc_gamma2_ende.resize(size);

  _mc_open_angle.resize(size);

  _reco_gamma1_nhits.resize(size);
  _reco_gamma1_npfps.resize(size);

  _reco_gamma1_pfp_track_score.resize(size);
  _reco_gamma1_comp.resize(size);
  _reco_gamma1_pur.resize(size);

  _reco_gamma1_has_good_track.resize(size);
  _reco_gamma1_track_start_x.resize(size);
  _reco_gamma1_track_start_y.resize(size);
  _reco_gamma1_track_start_z.resize(size);
  _reco_gamma1_track_end_x.resize(size);
  _reco_gamma1_track_end_y.resize(size);
  _reco_gamma1_track_end_z.resize(size);
  _reco_gamma1_track_dir_x.resize(size);
  _reco_gamma1_track_dir_y.resize(size);
  _reco_gamma1_track_dir_z.resize(size);
  _reco_gamma1_track_start_theta.resize(size);
  _reco_gamma1_track_start_phi.resize(size);
  _reco_gamma1_track_length.resize(size);
  _reco_gamma1_track_best_plane_energy.resize(size);

  _reco_gamma1_has_good_shower.resize(size);
  _reco_gamma1_shower_start_x.resize(size);
  _reco_gamma1_shower_start_y.resize(size);
  _reco_gamma1_shower_start_z.resize(size);
  _reco_gamma1_shower_dir_x.resize(size);
  _reco_gamma1_shower_dir_y.resize(size);
  _reco_gamma1_shower_dir_z.resize(size);
  _reco_gamma1_shower_start_theta.resize(size);
  _reco_gamma1_shower_start_phi.resize(size);
  _reco_gamma1_shower_length.resize(size);
  _reco_gamma1_shower_open_angle.resize(size);
  _reco_gamma1_shower_best_plane_energy.resize(size);
  _reco_gamma1_shower_best_plane_dedx.resize(size);

  _reco_gamma1_bestpfp_pfp_track_score.resize(size);
  _reco_gamma1_bestpfp_comp.resize(size);
  _reco_gamma1_bestpfp_pur.resize(size);

  _reco_gamma1_bestpfp_has_good_track.resize(size);
  _reco_gamma1_bestpfp_track_start_x.resize(size);
  _reco_gamma1_bestpfp_track_start_y.resize(size);
  _reco_gamma1_bestpfp_track_start_z.resize(size);
  _reco_gamma1_bestpfp_track_end_x.resize(size);
  _reco_gamma1_bestpfp_track_end_y.resize(size);
  _reco_gamma1_bestpfp_track_end_z.resize(size);
  _reco_gamma1_bestpfp_track_dir_x.resize(size);
  _reco_gamma1_bestpfp_track_dir_y.resize(size);
  _reco_gamma1_bestpfp_track_dir_z.resize(size);
  _reco_gamma1_bestpfp_track_start_theta.resize(size);
  _reco_gamma1_bestpfp_track_start_phi.resize(size);
  _reco_gamma1_bestpfp_track_length.resize(size);
  _reco_gamma1_bestpfp_track_best_plane_energy.resize(size);

  _reco_gamma1_bestpfp_has_good_shower.resize(size);
  _reco_gamma1_bestpfp_shower_start_x.resize(size);
  _reco_gamma1_bestpfp_shower_start_y.resize(size);
  _reco_gamma1_bestpfp_shower_start_z.resize(size);
  _reco_gamma1_bestpfp_shower_dir_x.resize(size);
  _reco_gamma1_bestpfp_shower_dir_y.resize(size);
  _reco_gamma1_bestpfp_shower_dir_z.resize(size);
  _reco_gamma1_bestpfp_shower_start_theta.resize(size);
  _reco_gamma1_bestpfp_shower_start_phi.resize(size);
  _reco_gamma1_bestpfp_shower_length.resize(size);
  _reco_gamma1_bestpfp_shower_open_angle.resize(size);
  _reco_gamma1_bestpfp_shower_best_plane_energy.resize(size);
  _reco_gamma1_bestpfp_shower_best_plane_dedx.resize(size);

  _reco_gamma2_nhits.resize(size);
  _reco_gamma2_npfps.resize(size);

  _reco_gamma2_pfp_track_score.resize(size);
  _reco_gamma2_comp.resize(size);
  _reco_gamma2_pur.resize(size);

  _reco_gamma2_has_good_track.resize(size);
  _reco_gamma2_track_start_x.resize(size);
  _reco_gamma2_track_start_y.resize(size);
  _reco_gamma2_track_start_z.resize(size);
  _reco_gamma2_track_end_x.resize(size);
  _reco_gamma2_track_end_y.resize(size);
  _reco_gamma2_track_end_z.resize(size);
  _reco_gamma2_track_dir_x.resize(size);
  _reco_gamma2_track_dir_y.resize(size);
  _reco_gamma2_track_dir_z.resize(size);
  _reco_gamma2_track_start_theta.resize(size);
  _reco_gamma2_track_start_phi.resize(size);
  _reco_gamma2_track_length.resize(size);
  _reco_gamma2_track_best_plane_energy.resize(size);

  _reco_gamma2_has_good_shower.resize(size);
  _reco_gamma2_shower_start_x.resize(size);
  _reco_gamma2_shower_start_y.resize(size);
  _reco_gamma2_shower_start_z.resize(size);
  _reco_gamma2_shower_dir_x.resize(size);
  _reco_gamma2_shower_dir_y.resize(size);
  _reco_gamma2_shower_dir_z.resize(size);
  _reco_gamma2_shower_start_theta.resize(size);
  _reco_gamma2_shower_start_phi.resize(size);
  _reco_gamma2_shower_length.resize(size);
  _reco_gamma2_shower_open_angle.resize(size);
  _reco_gamma2_shower_best_plane_energy.resize(size);
  _reco_gamma2_shower_best_plane_dedx.resize(size);

  _reco_gamma2_bestpfp_pfp_track_score.resize(size);
  _reco_gamma2_bestpfp_comp.resize(size);
  _reco_gamma2_bestpfp_pur.resize(size);

  _reco_gamma2_bestpfp_has_good_track.resize(size);
  _reco_gamma2_bestpfp_track_start_x.resize(size);
  _reco_gamma2_bestpfp_track_start_y.resize(size);
  _reco_gamma2_bestpfp_track_start_z.resize(size);
  _reco_gamma2_bestpfp_track_end_x.resize(size);
  _reco_gamma2_bestpfp_track_end_y.resize(size);
  _reco_gamma2_bestpfp_track_end_z.resize(size);
  _reco_gamma2_bestpfp_track_dir_x.resize(size);
  _reco_gamma2_bestpfp_track_dir_y.resize(size);
  _reco_gamma2_bestpfp_track_dir_z.resize(size);
  _reco_gamma2_bestpfp_track_start_theta.resize(size);
  _reco_gamma2_bestpfp_track_start_phi.resize(size);
  _reco_gamma2_bestpfp_track_length.resize(size);
  _reco_gamma2_bestpfp_track_best_plane_energy.resize(size);

  _reco_gamma2_bestpfp_has_good_shower.resize(size);
  _reco_gamma2_bestpfp_shower_start_x.resize(size);
  _reco_gamma2_bestpfp_shower_start_y.resize(size);
  _reco_gamma2_bestpfp_shower_start_z.resize(size);
  _reco_gamma2_bestpfp_shower_dir_x.resize(size);
  _reco_gamma2_bestpfp_shower_dir_y.resize(size);
  _reco_gamma2_bestpfp_shower_dir_z.resize(size);
  _reco_gamma2_bestpfp_shower_start_theta.resize(size);
  _reco_gamma2_bestpfp_shower_start_phi.resize(size);
  _reco_gamma2_bestpfp_shower_length.resize(size);
  _reco_gamma2_bestpfp_shower_open_angle.resize(size);
  _reco_gamma2_bestpfp_shower_best_plane_energy.resize(size);
  _reco_gamma2_bestpfp_shower_best_plane_dedx.resize(size);
}

void sbnd::PiZeroAnalysis::ResizeRecoSubVectors(const int i, const int size, const bool lead)
{
  if(lead)
    {
      _reco_gamma1_pfp_track_score[i].resize(size);
      _reco_gamma1_comp[i].resize(size);
      _reco_gamma1_pur[i].resize(size);

      _reco_gamma1_has_good_track[i].resize(size);
      _reco_gamma1_track_start_x[i].resize(size);
      _reco_gamma1_track_start_y[i].resize(size);
      _reco_gamma1_track_start_z[i].resize(size);
      _reco_gamma1_track_end_x[i].resize(size);
      _reco_gamma1_track_end_y[i].resize(size);
      _reco_gamma1_track_end_z[i].resize(size);
      _reco_gamma1_track_dir_x[i].resize(size);
      _reco_gamma1_track_dir_y[i].resize(size);
      _reco_gamma1_track_dir_z[i].resize(size);
      _reco_gamma1_track_start_theta[i].resize(size);
      _reco_gamma1_track_start_phi[i].resize(size);
      _reco_gamma1_track_length[i].resize(size);
      _reco_gamma1_track_best_plane_energy[i].resize(size);

      _reco_gamma1_has_good_shower[i].resize(size);
      _reco_gamma1_shower_start_x[i].resize(size);
      _reco_gamma1_shower_start_y[i].resize(size);
      _reco_gamma1_shower_start_z[i].resize(size);
      _reco_gamma1_shower_dir_x[i].resize(size);
      _reco_gamma1_shower_dir_y[i].resize(size);
      _reco_gamma1_shower_dir_z[i].resize(size);
      _reco_gamma1_shower_start_theta[i].resize(size);
      _reco_gamma1_shower_start_phi[i].resize(size);
      _reco_gamma1_shower_length[i].resize(size);
      _reco_gamma1_shower_open_angle[i].resize(size);
      _reco_gamma1_shower_best_plane_energy[i].resize(size);
      _reco_gamma1_shower_best_plane_dedx[i].resize(size);
    }
  else
    {
      _reco_gamma2_pfp_track_score[i].resize(size);
      _reco_gamma2_comp[i].resize(size);
      _reco_gamma2_pur[i].resize(size);

      _reco_gamma2_has_good_track[i].resize(size);
      _reco_gamma2_track_start_x[i].resize(size);
      _reco_gamma2_track_start_y[i].resize(size);
      _reco_gamma2_track_start_z[i].resize(size);
      _reco_gamma2_track_end_x[i].resize(size);
      _reco_gamma2_track_end_y[i].resize(size);
      _reco_gamma2_track_end_z[i].resize(size);
      _reco_gamma2_track_dir_x[i].resize(size);
      _reco_gamma2_track_dir_y[i].resize(size);
      _reco_gamma2_track_dir_z[i].resize(size);
      _reco_gamma2_track_start_theta[i].resize(size);
      _reco_gamma2_track_start_phi[i].resize(size);
      _reco_gamma2_track_length[i].resize(size);
      _reco_gamma2_track_best_plane_energy[i].resize(size);

      _reco_gamma2_has_good_shower[i].resize(size);
      _reco_gamma2_shower_start_x[i].resize(size);
      _reco_gamma2_shower_start_y[i].resize(size);
      _reco_gamma2_shower_start_z[i].resize(size);
      _reco_gamma2_shower_dir_x[i].resize(size);
      _reco_gamma2_shower_dir_y[i].resize(size);
      _reco_gamma2_shower_dir_z[i].resize(size);
      _reco_gamma2_shower_start_theta[i].resize(size);
      _reco_gamma2_shower_start_phi[i].resize(size);
      _reco_gamma2_shower_length[i].resize(size);
      _reco_gamma2_shower_open_angle[i].resize(size);
      _reco_gamma2_shower_best_plane_energy[i].resize(size);
      _reco_gamma2_shower_best_plane_dedx[i].resize(size);
    }
}

double sbnd::PiZeroAnalysis::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (objectHits.size() == 0) ? def_double : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}

double sbnd::PiZeroAnalysis::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (fHitsMap[trackID] == 0) ? def_double : objectHitsMap[trackID]/static_cast<double>(fHitsMap[trackID]);
}

bool sbnd::PiZeroAnalysis::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(sbnd::PiZeroAnalysis)
