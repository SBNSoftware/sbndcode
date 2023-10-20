////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroAnalysis
// Plugin Type: analyzer
// File:        NCPiZeroAnalysis_module.cc
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
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

#include "TTree.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include "NCPiZeroStructs.h"

#include <numeric>

constexpr int def_int       = std::numeric_limits<int>::min();
constexpr size_t def_size   = std::numeric_limits<size_t>::max();
constexpr float def_float   = -std::numeric_limits<float>::max();
constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NCPiZeroAnalysis;
}

class sbnd::NCPiZeroAnalysis : public art::EDAnalyzer {
public:
  explicit NCPiZeroAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCPiZeroAnalysis(NCPiZeroAnalysis const&) = delete;
  NCPiZeroAnalysis(NCPiZeroAnalysis&&) = delete;
  NCPiZeroAnalysis& operator=(NCPiZeroAnalysis const&) = delete;
  NCPiZeroAnalysis& operator=(NCPiZeroAnalysis&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;
  virtual void beginSubRun(const art::SubRun& sr);
  virtual void endSubRun(const art::SubRun& sr);

  void ResetSubRunVars();
  void ResetEventVars();
  void ClearMaps();

  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                 const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  int GetTotalGenEvents(const art::Event &e);

  void SetupBranches(VecVarMap &map);

  void ResizeVectors(VecVarMap &map, const int size);
  void ResizeSubVectors(VecVarMap &map, const std::string &subname, const int pos, const int size);

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalyseMCTruth(const art::Event &e, VecVarMap &vars, const art::Ptr<simb::MCTruth> &mct, const int counter, const std::string prefix);

  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                     const art::Handle<std::vector<recob::Track>> &trackHandle,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle);

  void AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const art::Ptr<recob::Vertex> &vtx, const int slcCounter,
                   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                   const art::Handle<std::vector<recob::Shower>> &showerHandle);
  void AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, const int slcCounter, const int pfpCounter,
                    const art::Handle<std::vector<recob::Track>> &trackHandle);
  void AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle, const art::Ptr<recob::Vertex> &vtx,
                     const std::vector<art::Ptr<recob::Hit>> &hits);
  void AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, const int slcCounter,
                         const art::Handle<std::vector<recob::Slice>> &sliceHandle);

  void ExtractDazzle(const art::Ptr<sbn::MVAPID> &dazzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<anab::Calorimetry> &calo, const int slcCounter, const int pfpCounter);
  void ExtractChi2PID(const art::Ptr<anab::ParticleID> &chi2pid, const int slcCounter, const int pfpCounter);
  void ExtractMCS(const art::Ptr<recob::MCSFitResult> &mcs, const int slcCounter, const int pfpCounter);
  void ExtractStoppingChi2(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2, const int slcCounter, const int pfpCounter);

  void ExtractRazzle(const art::Ptr<sbn::MVAPID> &razzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                   const std::vector<art::Ptr<recob::Hit>> &hits);

  void ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter);

  void SelectSlice(const int counter);

  void ProduceMultiSliceCandidates(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  void ProducePiZeroCandidates(VecVarMap &vars, const std::string &prefix,
                               const int counter, const std::vector<int> slc_ids);

  void ProducePiZeroCandidate(VecVarMap &vars, const std::string &prefix, const int counter, const int pzcCounter,
                              const TVector3 &dir0, const TVector3 &dir1, const double &en0, const double &en1);

  float Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  float Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  double DCA(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
             const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const size_t key0, const size_t key1);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  template<typename T>
  void FillElement(VecVar *vec, const int pos, const T value);
  template<typename T>
  void FillElement(VecVar *vec, const int posA, const int posB, const T value);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

  template<typename T>
  void AccessElement(VecVar *vec, const int pos, T &value);
  template<typename T>
  void AccessElement(VecVar *vec, const int posA, const int posB, T &value);
  template<typename T>
  void GetVar(VecVar *vec, std::vector<T> &var);
  template<typename T>
  void GetVar(VecVar *vec, std::vector<std::vector<T>> &var);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  art::InputTag fMCParticleModuleLabel, fSliceModuleLabel, fPFParticleModuleLabel, fVertexModuleLabel,
    fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fTrackCalorimetryModuleLabel,
    fCRUMBSModuleLabel, fDazzleModuleLabel, fCaloModuleLabel, fMCSModuleLabel, fChi2ModuleLabel, fRangeModuleLabel,
    fClosestApproachModuleLabel, fStoppingChi2ModuleLabel, fRazzleModuleLabel, fCosmicDistModuleLabel,
    fShowerTrackFitModuleLabel, fShowerDensityFitModuleLabel, fPOTModuleLabel, fOpT0ModuleLabel, fRazzledModuleLabel,
    fSpacePointModuleLabel;
  bool fDebug, fBeamOff;

  std::map<int, int> fHitsMap;
  std::map<const art::Ptr<simb::MCTruth>, int> fNuHitsMap;
  std::map<int, art::Ptr<recob::PFParticle>> fPFPMap;
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;

  TTree* fSubRunTree;

  double _pot;
  int    _spills, _ngenevts;

  TTree* fEventTree;

  int  _run;
  int  _subrun;
  int  _event;

  int _n_nu;

  VecVarMap nuVars = {
    { "nu_mctruth_id", new InhVecVar<size_t>("nu_mctruth_id") },
    { "nu_event_type", new InhVecVar<int>("nu_event_type") },
    { "nu_signal", new InhVecVar<bool>("nu_signal") },
    { "nu_en_dep", new InhVecVar<float>("nu_en_dep") },
    { "nu_pdg", new InhVecVar<int>("nu_pdg") },
    { "nu_ccnc", new InhVecVar<int>("nu_ccnc") },
    { "nu_av", new InhVecVar<bool>("nu_av") },
    { "nu_fv", new InhVecVar<bool>("nu_fv") },
    { "nu_mode", new InhVecVar<int>("nu_mode") },
    { "nu_int_type", new InhVecVar<int>("nu_int_type") },
    { "nu_n_protons", new InhVecVar<int>("nu_n_protons") },
    { "nu_n_neutrons", new InhVecVar<int>("nu_n_neutrons") },
    { "nu_n_charged_pions", new InhVecVar<int>("nu_n_charged_pions") },
    { "nu_n_neutral_pions", new InhVecVar<int>("nu_n_neutral_pions") },
    { "nu_n_photons", new InhVecVar<int>("nu_n_photons") },
    { "nu_n_other", new InhVecVar<int>("nu_n_other") },
    { "nu_w", new InhVecVar<double>("nu_w") },
    { "nu_x", new InhVecVar<double>("nu_x") },
    { "nu_y", new InhVecVar<double>("nu_y") },
    { "nu_q_sqr", new InhVecVar<double>("nu_q_sqr") },
    { "nu_pt", new InhVecVar<double>("nu_pt") },
    { "nu_theta", new InhVecVar<double>("nu_theta") },
    { "nu_e", new InhVecVar<double>("nu_e") },
    { "nu_vtx_x", new InhVecVar<double>("nu_vtx_x") },
    { "nu_vtx_y", new InhVecVar<double>("nu_vtx_y") },
    { "nu_vtx_z", new InhVecVar<double>("nu_vtx_z") },
    { "nu_n_pzs", new InhVecVar<size_t>("nu_n_pzs") },
    { "nu_pz_invariant_mass", new InhVecVecVar<double>("nu_pz_invariant_mass") },
    { "nu_pz_pizero_mom", new InhVecVecVar<double>("nu_pz_pizero_mom") },
    { "nu_pz_cos_theta_pizero", new InhVecVecVar<double>("nu_pz_cos_theta_pizero") },
    { "nu_pz_cos_com", new InhVecVecVar<double>("nu_pz_cos_com") },
    { "nu_pz_decay_asymmetry", new InhVecVecVar<double>("nu_pz_decay_asymmetry") },
    { "nu_pz_two_gamma_decay", new InhVecVecVar<bool>("nu_pz_two_gamma_decay") },
  };

  int _n_slc;

  VecVarMap slcVars = {
    { "slc_key", new InhVecVar<size_t>("slc_key") },
    { "slc_n_pfps", new InhVecVar<size_t>("slc_n_pfps") },
    { "slc_primary_pfp_id", new InhVecVar<size_t>("slc_primary_pfp_id") },
    { "slc_primary_pfp_pdg", new InhVecVar<int>("slc_primary_pfp_pdg") },
    { "slc_is_clear_cosmic", new InhVecVar<bool>("slc_is_clear_cosmic") },
    { "slc_n_primary_daughters", new InhVecVar<int>("slc_n_primary_daughters") },
    { "slc_n_trks", new InhVecVar<int>("slc_n_trks") },
    { "slc_n_shws", new InhVecVar<int>("slc_n_shws") },
    { "slc_n_dazzle_muons", new InhVecVar<int>("slc_n_dazzle_muons") },
    { "slc_n_dazzle_muons_cut_based", new InhVecVar<int>("slc_n_dazzle_muons_cut_based") },
    { "slc_n_dazzle_pions", new InhVecVar<int>("slc_n_dazzle_pions") },
    { "slc_n_dazzle_protons", new InhVecVar<int>("slc_n_dazzle_protons") },
    { "slc_n_dazzle_other", new InhVecVar<int>("slc_n_dazzle_other") },
    { "slc_n_razzle_electrons", new InhVecVar<int>("slc_n_razzle_electrons") },
    { "slc_n_razzle_photons", new InhVecVar<int>("slc_n_razzle_photons") },
    { "slc_n_razzle_photons_cut_based", new InhVecVar<int>("slc_n_razzle_photons_cut_based") },
    { "slc_n_razzle_other", new InhVecVar<int>("slc_n_razzle_other") },
    { "slc_n_razzled_electrons", new InhVecVar<int>("slc_n_razzled_electrons") },
    { "slc_n_razzled_muons", new InhVecVar<int>("slc_n_razzled_muons") },
    { "slc_n_razzled_photons", new InhVecVar<int>("slc_n_razzled_photons") },
    { "slc_n_razzled_pions", new InhVecVar<int>("slc_n_razzled_pions") },
    { "slc_n_razzled_protons", new InhVecVar<int>("slc_n_razzled_protons") },
    { "slc_true_mctruth_id", new InhVecVar<size_t>("slc_true_mctruth_id") },
    { "slc_true_event_type", new InhVecVar<int>("slc_true_event_type") },
    { "slc_true_signal", new InhVecVar<bool>("slc_true_signal") },
    { "slc_comp", new InhVecVar<float>("slc_comp") },
    { "slc_pur", new InhVecVar<float>("slc_pur") },
    { "slc_true_en_dep", new InhVecVar<float>("slc_true_en_dep") },
    { "slc_true_pdg", new InhVecVar<int>("slc_true_pdg") },
    { "slc_true_ccnc", new InhVecVar<int>("slc_true_ccnc") },
    { "slc_true_av", new InhVecVar<bool>("slc_true_av") },
    { "slc_true_fv", new InhVecVar<bool>("slc_true_fv") },
    { "slc_true_mode", new InhVecVar<int>("slc_true_mode") },
    { "slc_true_int_type", new InhVecVar<int>("slc_true_int_type") },
    { "slc_true_n_protons", new InhVecVar<int>("slc_true_n_protons") },
    { "slc_true_n_neutrons", new InhVecVar<int>("slc_true_n_neutrons") },
    { "slc_true_n_charged_pions", new InhVecVar<int>("slc_true_n_charged_pions") },
    { "slc_true_n_neutral_pions", new InhVecVar<int>("slc_true_n_neutral_pions") },
    { "slc_true_n_photons", new InhVecVar<int>("slc_true_n_photons") },
    { "slc_true_n_other", new InhVecVar<int>("slc_true_n_other") },
    { "slc_true_w", new InhVecVar<double>("slc_true_w") },
    { "slc_true_x", new InhVecVar<double>("slc_true_x") },
    { "slc_true_y", new InhVecVar<double>("slc_true_y") },
    { "slc_true_q_sqr", new InhVecVar<double>("slc_true_q_sqr") },
    { "slc_true_pt", new InhVecVar<double>("slc_true_pt") },
    { "slc_true_theta", new InhVecVar<double>("slc_true_theta") },
    { "slc_true_e", new InhVecVar<double>("slc_true_e") },
    { "slc_true_vtx_x", new InhVecVar<double>("slc_true_vtx_x") },
    { "slc_true_vtx_y", new InhVecVar<double>("slc_true_vtx_y") },
    { "slc_true_vtx_z", new InhVecVar<double>("slc_true_vtx_z") },
    { "slc_true_n_pzs", new InhVecVar<size_t>("slc_true_n_pzs") },
    { "slc_true_pz_invariant_mass", new InhVecVecVar<double>("slc_true_pz_invariant_mass") },
    { "slc_true_pz_pizero_mom", new InhVecVecVar<double>("slc_true_pz_pizero_mom") },
    { "slc_true_pz_cos_theta_pizero", new InhVecVecVar<double>("slc_true_pz_cos_theta_pizero") },
    { "slc_true_pz_cos_com", new InhVecVecVar<double>("slc_true_pz_cos_com") },
    { "slc_true_pz_decay_asymmetry", new InhVecVecVar<double>("slc_true_pz_decay_asymmetry") },
    { "slc_true_pz_two_gamma_decay", new InhVecVecVar<bool>("slc_true_pz_two_gamma_decay") },
    { "slc_vtx_x", new InhVecVar<double>("slc_vtx_x") },
    { "slc_vtx_y", new InhVecVar<double>("slc_vtx_y") },
    { "slc_vtx_z", new InhVecVar<double>("slc_vtx_z") },
    { "slc_is_fv", new InhVecVar<bool>("slc_is_fv") },
    { "slc_crumbs_score", new InhVecVar<float>("slc_crumbs_score") },
    { "slc_crumbs_nc_score", new InhVecVar<float>("slc_crumbs_nc_score") },
    { "slc_crumbs_ccnue_score", new InhVecVar<float>("slc_crumbs_ccnue_score") },
    { "slc_crumbs_ccnumu_score", new InhVecVar<float>("slc_crumbs_ccnumu_score") },
    { "slc_opt0_time", new InhVecVar<double>("slc_opt0_time") },
    { "slc_opt0_score", new InhVecVar<double>("slc_opt0_score") },
    { "slc_opt0_measPE", new InhVecVar<double>("slc_opt0_measPE") },
    { "slc_opt0_hypPE", new InhVecVar<double>("slc_opt0_hypPE") },
    { "slc_pfp_id", new InhVecVecVar<size_t>("slc_pfp_id") },
    { "slc_pfp_pdg", new InhVecVecVar<int>("slc_pfp_pdg") },
    { "slc_pfp_track_score", new InhVecVecVar<float>("slc_pfp_track_score") },
    { "slc_pfp_good_track", new InhVecVecVar<bool>("slc_pfp_good_track") },
    { "slc_pfp_good_shower", new InhVecVecVar<bool>("slc_pfp_good_shower") },
    { "slc_pfp_true_trackid", new InhVecVecVar<int>("slc_pfp_true_trackid") },
    { "slc_pfp_true_pdg", new InhVecVecVar<int>("slc_pfp_true_pdg") },
    { "slc_pfp_comp", new InhVecVecVar<float>("slc_pfp_comp") },
    { "slc_pfp_pur", new InhVecVecVar<float>("slc_pfp_pur") },
    { "slc_pfp_track_start_x", new InhVecVecVar<double>("slc_pfp_track_start_x") },
    { "slc_pfp_track_start_y", new InhVecVecVar<double>("slc_pfp_track_start_y") },
    { "slc_pfp_track_start_z", new InhVecVecVar<double>("slc_pfp_track_start_z") },
    { "slc_pfp_track_dir_x", new InhVecVecVar<double>("slc_pfp_track_dir_x") },
    { "slc_pfp_track_dir_y", new InhVecVecVar<double>("slc_pfp_track_dir_y") },
    { "slc_pfp_track_dir_z", new InhVecVecVar<double>("slc_pfp_track_dir_z") },
    { "slc_pfp_track_length", new InhVecVecVar<double>("slc_pfp_track_length") },
    { "slc_pfp_track_dazzle_muon_score", new InhVecVecVar<float>("slc_pfp_track_dazzle_muon_score") },
    { "slc_pfp_track_dazzle_pion_score", new InhVecVecVar<float>("slc_pfp_track_dazzle_pion_score") },
    { "slc_pfp_track_dazzle_proton_score", new InhVecVecVar<float>("slc_pfp_track_dazzle_proton_score") },
    { "slc_pfp_track_dazzle_other_score", new InhVecVecVar<float>("slc_pfp_track_dazzle_other_score") },
    { "slc_pfp_track_dazzle_pdg", new InhVecVecVar<int>("slc_pfp_track_dazzle_pdg") },
    { "slc_pfp_track_ke", new InhVecVecVar<float>("slc_pfp_track_ke") },
    { "slc_pfp_track_charge", new InhVecVecVar<float>("slc_pfp_track_charge") },
    { "slc_pfp_track_chi2_muon", new InhVecVecVar<float>("slc_pfp_track_chi2_muon") },
    { "slc_pfp_track_chi2_pion", new InhVecVecVar<float>("slc_pfp_track_chi2_pion") },
    { "slc_pfp_track_chi2_kaon", new InhVecVecVar<float>("slc_pfp_track_chi2_kaon") },
    { "slc_pfp_track_chi2_proton", new InhVecVecVar<float>("slc_pfp_track_chi2_proton") },
    { "slc_pfp_track_chi2_pdg", new InhVecVecVar<int>("slc_pfp_track_chi2_pdg") },
    { "slc_pfp_track_mcs_mom", new InhVecVecVar<float>("slc_pfp_track_mcs_mom") },
    { "slc_pfp_track_mcs_mean_scatter", new InhVecVecVar<float>("slc_pfp_track_mcs_mean_scatter") },
    { "slc_pfp_track_mcs_max_scatter_ratio", new InhVecVecVar<float>("slc_pfp_track_mcs_max_scatter_ratio") },
    { "slc_pfp_track_range_p", new InhVecVecVar<float>("slc_pfp_track_range_p") },
    { "slc_pfp_track_closest_approach_mean_dca", new InhVecVecVar<float>("slc_pfp_track_closest_approach_mean_dca") },
    { "slc_pfp_track_stopping_dedx_chi2_ratio", new InhVecVecVar<float>("slc_pfp_track_stopping_dedx_chi2_ratio") },
    { "slc_pfp_track_stopping_dedx_pol0_fit", new InhVecVecVar<float>("slc_pfp_track_stopping_dedx_pol0_fit") },
    { "slc_pfp_shower_start_x", new InhVecVecVar<double>("slc_pfp_shower_start_x") },
    { "slc_pfp_shower_start_y", new InhVecVecVar<double>("slc_pfp_shower_start_y") },
    { "slc_pfp_shower_start_z", new InhVecVecVar<double>("slc_pfp_shower_start_z") },
    { "slc_pfp_shower_conv_gap", new InhVecVecVar<double>("slc_pfp_shower_conv_gap") },
    { "slc_pfp_shower_dir_x", new InhVecVecVar<double>("slc_pfp_shower_dir_x") },
    { "slc_pfp_shower_dir_y", new InhVecVecVar<double>("slc_pfp_shower_dir_y") },
    { "slc_pfp_shower_dir_z", new InhVecVecVar<double>("slc_pfp_shower_dir_z") },
    { "slc_pfp_shower_length", new InhVecVecVar<double>("slc_pfp_shower_length") },
    { "slc_pfp_shower_open_angle", new InhVecVecVar<double>("slc_pfp_shower_open_angle") },
    { "slc_pfp_shower_energy", new InhVecVecVar<double>("slc_pfp_shower_energy") },
    { "slc_pfp_shower_dedx", new InhVecVecVar<double>("slc_pfp_shower_dedx") },
    { "slc_pfp_shower_sqrt_energy_density", new InhVecVecVar<double>("slc_pfp_shower_sqrt_energy_density") },
    { "slc_pfp_shower_modified_hit_density", new InhVecVecVar<double>("slc_pfp_shower_modified_hit_density") },
    { "slc_pfp_shower_razzle_electron_score", new InhVecVecVar<float>("slc_pfp_shower_razzle_electron_score") },
    { "slc_pfp_shower_razzle_photon_score", new InhVecVecVar<float>("slc_pfp_shower_razzle_photon_score") },
    { "slc_pfp_shower_razzle_other_score", new InhVecVecVar<float>("slc_pfp_shower_razzle_other_score") },
    { "slc_pfp_shower_razzle_pdg", new InhVecVecVar<int>("slc_pfp_shower_razzle_pdg") },
    { "slc_pfp_shower_cosmic_dist", new InhVecVecVar<float>("slc_pfp_shower_cosmic_dist") },
    { "slc_pfp_shower_track_length", new InhVecVecVar<double>("slc_pfp_shower_track_length") },
    { "slc_pfp_shower_track_width", new InhVecVecVar<double>("slc_pfp_shower_track_width") },
    { "slc_pfp_shower_density_grad", new InhVecVecVar<double>("slc_pfp_shower_density_grad") },
    { "slc_pfp_shower_density_pow", new InhVecVecVar<double>("slc_pfp_shower_density_pow") },
    { "slc_pfp_razzled_electron_score", new InhVecVecVar<float>("slc_pfp_razzled_electron_score") },
    { "slc_pfp_razzled_muon_score", new InhVecVecVar<float>("slc_pfp_razzled_muon_score") },
    { "slc_pfp_razzled_photon_score", new InhVecVecVar<float>("slc_pfp_razzled_photon_score") },
    { "slc_pfp_razzled_pion_score", new InhVecVecVar<float>("slc_pfp_razzled_pion_score") },
    { "slc_pfp_razzled_proton_score", new InhVecVecVar<float>("slc_pfp_razzled_proton_score") },
    { "slc_pfp_razzled_pdg", new InhVecVecVar<int>("slc_pfp_razzled_pdg") },
    { "slc_n_pzcs", new InhVecVar<size_t>("slc_n_pzcs") },
    { "slc_pzc_photon_0_id", new InhVecVecVar<int>("slc_pzc_photon_0_id") },
    { "slc_pzc_photon_1_id", new InhVecVecVar<int>("slc_pzc_photon_1_id") },
    { "slc_pzc_good_kinematics", new InhVecVecVar<bool>("slc_pzc_good_kinematics") },
    { "slc_pzc_invariant_mass", new InhVecVecVar<double>("slc_pzc_invariant_mass") },
    { "slc_pzc_pizero_mom", new InhVecVecVar<double>("slc_pzc_pizero_mom") },
    { "slc_pzc_cos_theta_pizero", new InhVecVecVar<double>("slc_pzc_cos_theta_pizero") },
    { "slc_pzc_cos_com", new InhVecVecVar<double>("slc_pzc_cos_com") },
    { "slc_pzc_decay_asymmetry", new InhVecVecVar<double>("slc_pzc_decay_asymmetry") },
    { "slc_sel1", new InhVecVar<bool>("slc_sel1") },
    { "slc_sel2", new InhVecVar<bool>("slc_sel2") },
  };

  int _n_multi_slice_candidates;

  VecVarMap mscVars = {
    { "msc_0_slc_id", new InhVecVar<int>("msc_0_slc_id") },
    { "msc_0_true_mctruth_id", new InhVecVar<size_t>("msc_0_true_mctruth_id") },
    { "msc_0_true_event_type", new InhVecVar<int>("msc_0_true_event_type") },
    { "msc_0_comp", new InhVecVar<float>("msc_0_comp") },
    { "msc_0_pur", new InhVecVar<float>("msc_0_pur") },
    { "msc_0_vtx_x", new InhVecVar<double>("msc_0_vtx_x") },
    { "msc_0_vtx_y", new InhVecVar<double>("msc_0_vtx_y") },
    { "msc_0_vtx_z", new InhVecVar<double>("msc_0_vtx_z") },
    { "msc_0_opt0_measPE", new InhVecVar<double>("msc_0_opt0_measPE") },
    { "msc_0_opt0_hypPE", new InhVecVar<double>("msc_0_opt0_hypPE") },
    { "msc_0_crumbs_score", new InhVecVar<float>("msc_0_crumbs_score") },
    { "msc_0_good_opt0", new InhVecVar<bool>("msc_0_good_opt0") },
    { "msc_0_opt0_frac", new InhVecVar<double>("msc_0_opt0_frac") },
    { "msc_0_n_pfps", new InhVecVar<size_t>("msc_0_n_pfps") },
    { "msc_0_n_razzle_electrons", new InhVecVar<int>("msc_0_n_razzle_electrons") },
    { "msc_0_n_razzle_photons", new InhVecVar<int>("msc_0_n_razzle_photons") },
    { "msc_0_n_razzle_other", new InhVecVar<int>("msc_0_n_razzle_other") },
    { "msc_0_n_dazzle_muons", new InhVecVar<int>("msc_0_n_dazzle_muons") },
    { "msc_0_n_dazzle_pions", new InhVecVar<int>("msc_0_n_dazzle_pions") },
    { "msc_0_n_dazzle_protons", new InhVecVar<int>("msc_0_n_dazzle_protons") },
    { "msc_0_n_dazzle_other", new InhVecVar<int>("msc_0_n_dazzle_other") },
    { "msc_0_n_razzled_electrons", new InhVecVar<int>("msc_0_n_razzled_electrons") },
    { "msc_0_n_razzled_muons", new InhVecVar<int>("msc_0_n_razzled_muons") },
    { "msc_0_n_razzled_photons", new InhVecVar<int>("msc_0_n_razzled_photons") },
    { "msc_0_n_razzled_pions", new InhVecVar<int>("msc_0_n_razzled_pions") },
    { "msc_0_n_razzled_protons", new InhVecVar<int>("msc_0_n_razzled_protons") },
    { "msc_1_slc_id", new InhVecVar<int>("msc_1_slc_id") },
    { "msc_1_true_mctruth_id", new InhVecVar<size_t>("msc_1_true_mctruth_id") },
    { "msc_1_true_event_type", new InhVecVar<int>("msc_1_true_event_type") },
    { "msc_1_comp", new InhVecVar<float>("msc_1_comp") },
    { "msc_1_pur", new InhVecVar<float>("msc_1_pur") },
    { "msc_1_vtx_x", new InhVecVar<double>("msc_1_vtx_x") },
    { "msc_1_vtx_y", new InhVecVar<double>("msc_1_vtx_y") },
    { "msc_1_vtx_z", new InhVecVar<double>("msc_1_vtx_z") },
    { "msc_1_opt0_measPE", new InhVecVar<double>("msc_1_opt0_measPE") },
    { "msc_1_opt0_hypPE", new InhVecVar<double>("msc_1_opt0_hypPE") },
    { "msc_1_crumbs_score", new InhVecVar<float>("msc_1_crumbs_score") },
    { "msc_1_good_opt0", new InhVecVar<bool>("msc_1_good_opt0") },
    { "msc_1_opt0_frac", new InhVecVar<double>("msc_1_opt0_frac") },
    { "msc_1_n_pfps", new InhVecVar<size_t>("msc_1_n_pfps") },
    { "msc_1_n_razzle_electrons", new InhVecVar<int>("msc_1_n_razzle_electrons") },
    { "msc_1_n_razzle_photons", new InhVecVar<int>("msc_1_n_razzle_photons") },
    { "msc_1_n_razzle_other", new InhVecVar<int>("msc_1_n_razzle_other") },
    { "msc_1_n_dazzle_muons", new InhVecVar<int>("msc_1_n_dazzle_muons") },
    { "msc_1_n_dazzle_pions", new InhVecVar<int>("msc_1_n_dazzle_pions") },
    { "msc_1_n_dazzle_protons", new InhVecVar<int>("msc_1_n_dazzle_protons") },
    { "msc_1_n_dazzle_other", new InhVecVar<int>("msc_1_n_dazzle_other") },
    { "msc_1_n_razzled_electrons", new InhVecVar<int>("msc_1_n_razzled_electrons") },
    { "msc_1_n_razzled_muons", new InhVecVar<int>("msc_1_n_razzled_muons") },
    { "msc_1_n_razzled_photons", new InhVecVar<int>("msc_1_n_razzled_photons") },
    { "msc_1_n_razzled_pions", new InhVecVar<int>("msc_1_n_razzled_pions") },
    { "msc_1_n_razzled_protons", new InhVecVar<int>("msc_1_n_razzled_protons") },
    { "msc_signal", new InhVecVar<bool>("msc_signal") },
    { "msc_sep", new InhVecVar<double>("msc_sep") },
    { "msc_matching_flash_pe", new InhVecVar<bool>("msc_matching_flash_pe") },
    { "msc_same_tpc", new InhVecVar<bool>("msc_same_tpc") },
    { "msc_sum_opt0_frac", new InhVecVar<double>("msc_sum_opt0_frac") },
    { "msc_dca", new InhVecVar<double>("msc_dca") },
    { "msc_n_pzcs", new InhVecVar<size_t>("msc_n_pzcs") },
    { "msc_pzc_photon_0_id", new InhVecVecVar<int>("msc_pzc_photon_0_id") },
    { "msc_pzc_photon_1_id", new InhVecVecVar<int>("msc_pzc_photon_1_id") },
    { "msc_pzc_photon_0_slc_id", new InhVecVecVar<int>("msc_pzc_photon_0_slc_id") },
    { "msc_pzc_photon_1_slc_id", new InhVecVecVar<int>("msc_pzc_photon_1_slc_id") },
    { "msc_pzc_good_kinematics", new InhVecVecVar<bool>("msc_pzc_good_kinematics") },
    { "msc_pzc_invariant_mass", new InhVecVecVar<double>("msc_pzc_invariant_mass") },
    { "msc_pzc_pizero_mom", new InhVecVecVar<double>("msc_pzc_pizero_mom") },
    { "msc_pzc_cos_theta_pizero", new InhVecVecVar<double>("msc_pzc_cos_theta_pizero") },
    { "msc_pzc_cos_com", new InhVecVecVar<double>("msc_pzc_cos_com") },
    { "msc_pzc_decay_asymmetry", new InhVecVecVar<double>("msc_pzc_decay_asymmetry") },
  };
};

sbnd::NCPiZeroAnalysis::NCPiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCParticleModuleLabel       = p.get<art::InputTag>("MCParticleModuleLabel", "largeant");
    fSliceModuleLabel            = p.get<art::InputTag>("SliceModuleLabel", "pandoraSCE");
    fPFParticleModuleLabel       = p.get<art::InputTag>("PFParticleModuleLabel", "pandoraSCE");
    fVertexModuleLabel           = p.get<art::InputTag>("VertexModuleLabel", "pandoraSCE");
    fHitModuleLabel              = p.get<art::InputTag>("HitModuleLabel", "gaushit");
    fTrackModuleLabel            = p.get<art::InputTag>("TrackModuleLabel", "pandoraTrack");
    fShowerModuleLabel           = p.get<art::InputTag>("ShowerModuleLabel", "pandoraShowerSBN");
    fTrackCalorimetryModuleLabel = p.get<art::InputTag>("TrackCalorimetryModuleLabel", "pandoraCalo");
    fCRUMBSModuleLabel           = p.get<art::InputTag>("CRUMBSModuleLabel", "crumbs");
    fDazzleModuleLabel           = p.get<art::InputTag>("DazzleModuleLabel", "dazzle");
    fCaloModuleLabel             = p.get<art::InputTag>("CaloModuleLabel", "pandoraCalo");
    fMCSModuleLabel              = p.get<art::InputTag>("MCSModuleLabel", "pandoraTrackMCS:muon");
    fChi2ModuleLabel             = p.get<art::InputTag>("Chi2ModuleLabel", "pandoraPid");
    fRangeModuleLabel            = p.get<art::InputTag>("RangeModuleLabel", "pandoraTrackRange:muon");
    fClosestApproachModuleLabel  = p.get<art::InputTag>("ClosestApproachModuleLabel", "pandoraTrackClosestApproach");
    fStoppingChi2ModuleLabel     = p.get<art::InputTag>("StoppingChi2ModuleLabel", "pandoraTrackStoppingChi2");
    fRazzleModuleLabel           = p.get<art::InputTag>("RazzleModuleLabel", "razzle");
    fCosmicDistModuleLabel       = p.get<art::InputTag>("CosmicDistModuleLabel", "pandoraShowerCosmicDist");
    fShowerTrackFitModuleLabel   = p.get<art::InputTag>("ShowerTrackFitModuleLabel", "pandoraShowerSelectionVars");
    fShowerDensityFitModuleLabel = p.get<art::InputTag>("ShowerDensityFitModuleLabel", "pandoraShowerSelectionVars");
    fPOTModuleLabel              = p.get<art::InputTag>("POTModuleLabel", "generator");
    fOpT0ModuleLabel             = p.get<art::InputTag>("OpT0ModuleLabel", "opt0finder");
    fRazzledModuleLabel          = p.get<art::InputTag>("RazzledModuleLabel", "razzled");
    fSpacePointModuleLabel       = p.get<art::InputTag>("SpacePointModuleLabel", "pandoraSCE");
    fDebug                       = p.get<bool>("Debug", false);
    fBeamOff                     = p.get<bool>("BeamOff", false);

    art::ServiceHandle<art::TFileService> fs;

    fSubRunTree = fs->make<TTree>("subruns","");

    fSubRunTree->Branch("pot", &_pot);
    fSubRunTree->Branch("spills", &_spills);
    fSubRunTree->Branch("ngenevts", &_ngenevts);

    fEventTree = fs->make<TTree>("events","");
    fEventTree->Branch("run", &_run);
    fEventTree->Branch("subrun", &_subrun);
    fEventTree->Branch("event", &_event);

    fEventTree->Branch("n_nu", &_n_nu);
    SetupBranches(nuVars);

    fEventTree->Branch("n_slc", &_n_slc);
    SetupBranches(slcVars);

    fEventTree->Branch("n_multi_slice_candidates", &_n_multi_slice_candidates);
    SetupBranches(mscVars);
  }

void sbnd::NCPiZeroAnalysis::SetupBranches(VecVarMap &map)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kOneD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<bool>*>(var)->Var()));
              break;
            case kInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<int>*>(var)->Var()));
              break;
            case kUInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<size_t>*>(var)->Var()));
              break;
            case kFloat:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<float>*>(var)->Var()));
              break;
            case kDouble:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<double>*>(var)->Var()));
              break;
            case kUnknownVar:
              break;
            }
        }
      else if(var->IdentifyVec() == kTwoD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<bool>*>(var)->Var()));
              break;
            case kInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<int>*>(var)->Var()));
              break;
            case kUInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<size_t>*>(var)->Var()));
              break;
            case kFloat:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<float>*>(var)->Var()));
              break;
            case kDouble:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<double>*>(var)->Var()));
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

void sbnd::NCPiZeroAnalysis::beginSubRun(const art::SubRun &sr)
{
  ResetSubRunVars();

  if(fBeamOff)
    return;

  // Get POT
  art::Handle<sumdata::POTSummary> potHandle;
  sr.getByLabel(fPOTModuleLabel, potHandle);
  if(!potHandle.isValid()){
    std::cout << "POT product " << fPOTModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  _pot    = potHandle->totpot;
  _spills = potHandle->totspills;
}

void sbnd::NCPiZeroAnalysis::endSubRun(const art::SubRun &sr)
{
  fSubRunTree->Fill();
}

void sbnd::NCPiZeroAnalysis::analyze(const art::Event &e)
{
  ResetEventVars();
  ClearMaps();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  // Note this can only be accessed from the event object but is a subrun level quantity.
  // Hence, we override it every event but it is only filled in the subrun tree.
  _ngenevts = GetTotalGenEvents(e);

  if(fDebug)
    std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();


  // Get Hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitModuleLabel, hitHandle);
  if(!hitHandle.isValid()){
    std::cout << "Hit product " << fHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Slices
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceModuleLabel, sliceHandle);
  if(!sliceHandle.isValid()){
    std::cout << "Slice product " << fSliceModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get PFParticles
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPFParticleModuleLabel, pfpHandle);
  if(!pfpHandle.isValid()){
    std::cout << "PFParticle product " << fPFParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Tracks
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackModuleLabel, trackHandle);
  if(!trackHandle.isValid()){
    std::cout << "Track product " << fTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Showers
  art::Handle<std::vector<recob::Shower>> showerHandle;
  e.getByLabel(fShowerModuleLabel, showerHandle);
  if(!showerHandle.isValid()){
    std::cout << "Shower product " << fShowerModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  SetupMaps(e, hitHandle, pfpHandle);
  AnalyseNeutrinos(e, MCTruthHandles);
  AnalyseSlices(e, sliceHandle, pfpHandle, trackHandle, showerHandle);
  ProduceMultiSliceCandidates(e, sliceHandle, pfpHandle);

  // Fill the Tree
  fEventTree->Fill();
}

void sbnd::NCPiZeroAnalysis::ClearMaps()
{
  fHitsMap.clear();
  fNuHitsMap.clear();
  fPFPMap.clear();
  fRecoPFPMap.clear();
}

void sbnd::NCPiZeroAnalysis::SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                                       const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  for(auto const& hit : hitVec)
    {
      const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      fHitsMap[trackID]++;
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      fNuHitsMap[mct]++;
    }

  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  art::fill_ptr_vector(pfpVec, pfpHandle);

  for(auto const& pfp : pfpVec)
    fPFPMap[pfp->Self()] = pfp;
}

void sbnd::NCPiZeroAnalysis::AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          ++_n_nu;
        }
    }

  ResizeVectors(nuVars, _n_nu);

  int nuCounter = 0;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          AnalyseMCTruth(e, nuVars, mct, nuCounter, "nu");

          ++nuCounter;
        }
    }
}

void sbnd::NCPiZeroAnalysis::AnalyseMCTruth(const art::Event &e, VecVarMap &vars, const art::Ptr<simb::MCTruth> &mct, const int counter, const std::string prefix)
{
  if(mct->Origin() == 2)
    {
      FillElement(vars[prefix + "_event_type"], counter, (int) kCosmic);
      return;
    }
  else if(mct->Origin() != 1)
    {
      FillElement(vars[prefix + "_event_type"], counter, (int) kUnknownEv);
      return;
    }

  FillElement(vars[prefix + "_mctruth_id"], counter, mct.key());

  const simb::MCNeutrino mcn = mct->GetNeutrino();
  const simb::MCParticle nu  = mcn.Nu();

  const bool nc = mcn.CCNC() == 1;
  const bool av = VolumeCheck(nu.Position().Vect());
  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

  art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);

  int protons = 0, neutrons = 0, charged_pions = 0, neutral_pions = 0, photons = 0, other = 0;
  float trueEnDep = 0.;

  for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
        {
          switch(abs(mcp->PdgCode()))
            {
            case 2212:
              if(mcp->P() > .25)
                ++protons;
              break;
            case 2112:
              ++neutrons;
              break;
            case 211:
              if(mcp->P() > .1)
                ++charged_pions;
              break;
            case 111:
              ++neutral_pions;
              break;
            case 22:
              ++photons;
              break;
            default:
              ++other;
              break;
            }
        }

      std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

      for(auto const& ide : ides)
        trueEnDep += ide->energy / 1000.;
    }

  const bool pizero = neutral_pions > 0;

  if(nc && fv && pizero)
    {
      FillElement(vars[prefix + "_event_type"], counter, (int) kNCPiZero);
      FillElement(vars[prefix + "_signal"], counter, true);
    }
  else
    {
      FillElement(vars[prefix + "_signal"], counter, false);

      if(nc && fv)
        FillElement(vars[prefix + "_event_type"], counter, (int) kOtherNC);
      else if(abs(nu.PdgCode()) == 14 && !nc && fv)
        FillElement(vars[prefix + "_event_type"], counter, (int) kCCNuMu);
      else if(abs(nu.PdgCode()) == 12 && !nc && fv)
        FillElement(vars[prefix + "_event_type"], counter, (int) kCCNuE);
      else if(!fv && av)
        FillElement(vars[prefix + "_event_type"], counter, (int) kNonFV);
      else if(!av)
        FillElement(vars[prefix + "_event_type"], counter, (int) kDirt);
      else
        FillElement(vars[prefix + "_event_type"], counter, (int) kUnknownEv);
    }

  FillElement(vars[prefix + "_en_dep"], counter, trueEnDep);
  FillElement(vars[prefix + "_pdg"], counter, nu.PdgCode());
  FillElement(vars[prefix + "_ccnc"], counter, mcn.CCNC());
  FillElement(vars[prefix + "_av"], counter, av);
  FillElement(vars[prefix + "_fv"], counter, fv);
  FillElement(vars[prefix + "_mode"], counter, mcn.Mode());
  FillElement(vars[prefix + "_int_type"], counter, mcn.InteractionType());
  FillElement(vars[prefix + "_n_protons"], counter, protons);
  FillElement(vars[prefix + "_n_neutrons"], counter, neutrons);
  FillElement(vars[prefix + "_n_charged_pions"], counter, charged_pions);
  FillElement(vars[prefix + "_n_neutral_pions"], counter, neutral_pions);
  FillElement(vars[prefix + "_n_photons"], counter, photons);
  FillElement(vars[prefix + "_n_other"], counter, other);
  FillElement(vars[prefix + "_w"], counter, mcn.W());
  FillElement(vars[prefix + "_x"], counter, mcn.X());
  FillElement(vars[prefix + "_y"], counter, mcn.Y());
  FillElement(vars[prefix + "_q_sqr"], counter, mcn.QSqr());
  FillElement(vars[prefix + "_pt"], counter, mcn.Pt());
  FillElement(vars[prefix + "_theta"], counter, mcn.Theta());
  FillElement(vars[prefix + "_e"], counter, nu.E());
  FillElement(vars[prefix + "_vtx_x"], counter, nu.Vx());
  FillElement(vars[prefix + "_vtx_y"], counter, nu.Vy());
  FillElement(vars[prefix + "_vtx_z"], counter, nu.Vz());

  FillElement(vars[prefix + "_n_pzs"], counter, (size_t) neutral_pions);
  ResizeSubVectors(vars, prefix + "_pz", counter, neutral_pions);

  int pzCounter = 0;

  for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
        {
          if(abs(mcp->PdgCode()) == 111)
            {
	      FillElement(vars[prefix + "_pz_invariant_mass"], counter, pzCounter, mcp->Mass());
	      FillElement(vars[prefix + "_pz_pizero_mom"], counter, pzCounter, mcp->P());
	      FillElement(vars[prefix + "_pz_cos_theta_pizero"], counter, pzCounter, mcp->Pz() / mcp->P());

	      bool two_gamma_decay = mcp->NumberDaughters() == 2;
	      if(!two_gamma_decay)
		{
		  FillElement(vars[prefix + "_pz_two_gamma_decay"], counter, pzCounter, two_gamma_decay);
		  ++pzCounter;
		  continue;
		}

	      const simb::MCParticle* gamma0 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
	      const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

	      two_gamma_decay &= (gamma0->PdgCode() == 22 && gamma1->PdgCode() == 22);
	      if(!two_gamma_decay)
		{
		  FillElement(vars[prefix + "_pz_two_gamma_decay"], counter, pzCounter, two_gamma_decay);
		  ++pzCounter;
		  continue;
		}

	      const double en0 = gamma0->E();
	      const double en1 = gamma1->E();

	      FillElement(vars[prefix + "_pz_cos_com"], counter, pzCounter, std::abs(en0 - en1) / mcp->P());
	      FillElement(vars[prefix + "_pz_decay_asymmetry"], counter, pzCounter, std::abs(en0 - en1) / (en0 + en1));
	      FillElement(vars[prefix + "_pz_two_gamma_decay"], counter, pzCounter, two_gamma_decay);

	      ++pzCounter;
	    }
	}
    }
}

void sbnd::NCPiZeroAnalysis::AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                           const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                                           const art::Handle<std::vector<recob::Track>> &trackHandle,
                                           const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  _n_slc = sliceVec.size();
  ResizeVectors(slcVars, _n_slc);

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>      pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  art::FindOneP<sbn::CRUMBSResult>  slicesToCRUMBS(sliceHandle, e, fCRUMBSModuleLabel);
  art::FindManyP<sbn::OpT0Finder>   slicesToOpT0(sliceHandle, e, fOpT0ModuleLabel);

  for (auto&& [slcCounter, slc] : enumerate(sliceVec))
    {
      FillElement(slcVars["slc_key"], slcCounter, slc.key());

      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slc.key());
      FillElement(slcVars["slc_n_pfps"], slcCounter, pfps.size());

      if(pfps.size() == 0)
        {
          FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, true);
          continue;
        }

      const art::Ptr<recob::PFParticle> prim = GetPrimaryPFP(pfps);
      if(prim.isNull())
        continue;

      FillElement(slcVars["slc_primary_pfp_id"], slcCounter, prim->Self());
      FillElement(slcVars["slc_primary_pfp_pdg"], slcCounter, prim->PdgCode());
      FillElement(slcVars["slc_n_primary_daughters"], slcCounter, prim->NumDaughters());

      if(abs(prim->PdgCode()) == 13 || abs(prim->PdgCode()) == 11)
        FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, true);
      else
        {
          FillElement(slcVars["slc_n_pfps"], slcCounter, pfps.size() - 1);
          FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, false);
        }

      const art::Ptr<recob::Vertex> vtx = pfpToVertices.at(prim.key());
      geo::Point_t vtxPos = vtx.isNonnull() ? vtx->position() : geo::Point_t(def_double, def_double, def_double);
      FillElement(slcVars["slc_vtx_x"], slcCounter, vtxPos.X());
      FillElement(slcVars["slc_vtx_y"], slcCounter, vtxPos.Y());
      FillElement(slcVars["slc_vtx_z"], slcCounter, vtxPos.Z());
      FillElement(slcVars["slc_is_fv"], slcCounter, VolumeCheck(vtxPos, 20., 5., 10., 50.));

      const art::Ptr<sbn::CRUMBSResult> crumbs = slicesToCRUMBS.at(slc.key());
      if(crumbs.isNonnull())
        {
          FillElement(slcVars["slc_crumbs_score"], slcCounter, crumbs->score);
          FillElement(slcVars["slc_crumbs_nc_score"], slcCounter, crumbs->ncscore);
          FillElement(slcVars["slc_crumbs_ccnue_score"], slcCounter, crumbs->ccnuescore);
          FillElement(slcVars["slc_crumbs_ccnumu_score"], slcCounter, crumbs->ccnumuscore);
        }

      std::vector<art::Ptr<sbn::OpT0Finder>> opT0Vec = slicesToOpT0.at(slc.key());
      if(opT0Vec.size() > 0)
        {
          std::sort(opT0Vec.begin(), opT0Vec.end(),
                    [](auto const& a, auto const& b)
                    { return a->score > b->score; });

          FillElement(slcVars["slc_opt0_time"], slcCounter, opT0Vec[0]->time);
          FillElement(slcVars["slc_opt0_score"], slcCounter, opT0Vec[0]->score);
          FillElement(slcVars["slc_opt0_measPE"], slcCounter, opT0Vec[0]->measPE);
          FillElement(slcVars["slc_opt0_hypPE"], slcCounter, opT0Vec[0]->hypoPE);
        }

      ResizeSubVectors(slcVars, "slc_pfp", slcCounter, prim->NumDaughters());

      AnalysePFPs(e, prim, vtx, slcCounter, pfpHandle, trackHandle, showerHandle);

      if(abs(prim->PdgCode()) == 13 || abs(prim->PdgCode()) == 11)
        ResizeSubVectors(slcVars, "slc_pzc", slcCounter, 0);
      else
        ProducePiZeroCandidates(slcVars, "slc", slcCounter, { (int) slcCounter });

      SelectSlice(slcCounter);

      AnalyseSliceTruth(e, slc, slcCounter, sliceHandle);
    }
}

void sbnd::NCPiZeroAnalysis::AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const art::Ptr<recob::Vertex> &vtx, const int slcCounter,
                                         const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                                         const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  art::FindOneP<recob::Track>                      pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpToMeta(pfpHandle, e, fPFParticleModuleLabel);
  art::FindManyP<recob::Hit>                       showersToHits(showerHandle, e, fShowerModuleLabel);
  art::FindOneP<sbn::MVAPID>                       pfpsToRazzled(pfpHandle, e, fRazzledModuleLabel);

  int ntrks = 0, nshws = 0, ndazzlemuons = 0, ndazzlepions = 0, ndazzleprotons = 0, ndazzleother = 0,
    nrazzleelectrons = 0, nrazzlephotons = 0, nrazzleother = 0, ndazzlemuonscut = 0, nrazzlephotonscut = 0,
    nrazzledelectrons = 0, nrazzledmuons = 0, nrazzledphotons = 0, nrazzledpions = 0, nrazzledprotons = 0;

  for(auto&& [pfpCounter, id] : enumerate(prim->Daughters()))
    {
      const art::Ptr<recob::PFParticle> pfp = fPFPMap[id];
      FillElement(slcVars["slc_pfp_id"], slcCounter, pfpCounter, id);
      FillElement(slcVars["slc_pfp_pdg"], slcCounter, pfpCounter, pfp->PdgCode());

      if(abs(pfp->PdgCode()) == 11)
        ++nshws;
      else if(abs(pfp->PdgCode()) == 13)
        ++ntrks;
      else
        {
          std::cout << "PFP with PDG: " << pfp->PdgCode() << std::endl;
          throw std::exception();
        }

      const art::Ptr<larpandoraobj::PFParticleMetadata> meta         = pfpToMeta.at(pfp.key());
      const larpandoraobj::PFParticleMetadata::PropertiesMap metaMap = meta->GetPropertiesMap();
      const std::map<std::string, float>::const_iterator scoreIter   = metaMap.find("TrackScore");

      if(scoreIter != metaMap.end())
        FillElement(slcVars["slc_pfp_track_score"], slcCounter, pfpCounter, scoreIter->second);

      const art::Ptr<recob::Track> track   = pfpToTrack.at(pfp.key());
      const art::Ptr<recob::Shower> shower = pfpToShower.at(pfp.key());

      FillElement(slcVars["slc_pfp_good_track"], slcCounter, pfpCounter, track.isNonnull());
      FillElement(slcVars["slc_pfp_good_shower"], slcCounter, pfpCounter, shower.isNonnull());

      const art::Ptr<sbn::MVAPID> razzled = pfpsToRazzled.at(pfp.key());
      if(razzled.isNonnull())
        ExtractRazzled(razzled, slcCounter, pfpCounter);

      const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
      const std::vector<art::Ptr<recob::Hit>> hits = showersToHits.at(shower.key());
      const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, true);
      FillElement(slcVars["slc_pfp_true_trackid"], slcCounter, pfpCounter, trackID);
      FillElement(slcVars["slc_pfp_comp"], slcCounter, pfpCounter, Completeness(e, hits, trackID));
      FillElement(slcVars["slc_pfp_pur"], slcCounter, pfpCounter, Purity(e, hits, trackID));

      if(trackID != def_int)
        {
          const simb::MCParticle* mcp = particleInv->TrackIdToParticle_P(trackID);
          if(mcp != NULL)
            FillElement(slcVars["slc_pfp_true_pdg"], slcCounter, pfpCounter, mcp->PdgCode());
        }

      if(track.isNonnull())
        AnalyseTrack(e, track, slcCounter, pfpCounter, trackHandle);

      if(pfp->PdgCode() == 13)
        {
          int dazzlepdg;
          AccessElement(slcVars["slc_pfp_track_dazzle_pdg"], slcCounter, pfpCounter, dazzlepdg);

          if(dazzlepdg == 13)
            ++ndazzlemuons;
          else if(dazzlepdg == 211)
            ++ndazzlepions;
          else if(dazzlepdg == 2212)
            ++ndazzleprotons;
          else if(dazzlepdg == 0)
            ++ndazzleother;

          float dazzlemuonscore;
          AccessElement(slcVars["slc_pfp_track_dazzle_muon_score"], slcCounter, pfpCounter, dazzlemuonscore);
          if(dazzlemuonscore > 0.0125)
            ++ndazzlemuonscut;
        }

      if(shower.isNonnull())
        AnalyseShower(e, shower, slcCounter, pfpCounter, showerHandle, vtx, hits);

      if(pfp->PdgCode() == 11)
        {
          int razzlepdg;
          AccessElement(slcVars["slc_pfp_shower_razzle_pdg"], slcCounter, pfpCounter, razzlepdg);

          if(razzlepdg == 11)
            ++nrazzleelectrons;
          else if(razzlepdg == 22)
            ++nrazzlephotons;
          else if(razzlepdg == 0)
            ++nrazzleother;

          float razzlephotonscore;
          AccessElement(slcVars["slc_pfp_shower_razzle_photon_score"], slcCounter, pfpCounter, razzlephotonscore);
          if(razzlephotonscore > 0.05)
            ++nrazzlephotonscut;
        }

      int razzledpdg;
      AccessElement(slcVars["slc_pfp_razzled_pdg"], slcCounter, pfpCounter, razzledpdg);

      if(razzledpdg == 11)
        ++nrazzledelectrons;
      else if(razzledpdg == 13)
        ++nrazzledmuons;
      else if(razzledpdg == 22)
        ++nrazzledphotons;
      else if(razzledpdg == 211)
        ++nrazzledpions;
      else if(razzledpdg == 2212)
        ++nrazzledprotons;
    }

  FillElement(slcVars["slc_n_trks"], slcCounter, ntrks);
  FillElement(slcVars["slc_n_shws"], slcCounter, nshws);
  FillElement(slcVars["slc_n_dazzle_muons"], slcCounter, ndazzlemuons);
  FillElement(slcVars["slc_n_dazzle_muons_cut_based"], slcCounter, ndazzlemuonscut);
  FillElement(slcVars["slc_n_dazzle_pions"], slcCounter, ndazzlepions);
  FillElement(slcVars["slc_n_dazzle_protons"], slcCounter, ndazzleprotons);
  FillElement(slcVars["slc_n_dazzle_other"], slcCounter, ndazzleother);
  FillElement(slcVars["slc_n_razzle_electrons"], slcCounter, nrazzleelectrons);
  FillElement(slcVars["slc_n_razzle_photons"], slcCounter, nrazzlephotons);
  FillElement(slcVars["slc_n_razzle_photons_cut_based"], slcCounter, nrazzlephotonscut);
  FillElement(slcVars["slc_n_razzle_other"], slcCounter, nrazzleother);
  FillElement(slcVars["slc_n_razzled_electrons"], slcCounter, nrazzledelectrons);
  FillElement(slcVars["slc_n_razzled_muons"], slcCounter, nrazzledmuons);
  FillElement(slcVars["slc_n_razzled_photons"], slcCounter, nrazzledphotons);
  FillElement(slcVars["slc_n_razzled_pions"], slcCounter, nrazzledpions);
  FillElement(slcVars["slc_n_razzled_protons"], slcCounter, nrazzledprotons);
}

void sbnd::NCPiZeroAnalysis::AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, const int slcCounter, const int pfpCounter,
                                          const art::Handle<std::vector<recob::Track>> &trackHandle)
{
  art::FindOneP<sbn::MVAPID> tracksToDazzle(trackHandle, e, fDazzleModuleLabel);
  art::FindManyP<anab::Calorimetry> tracksToCalos(trackHandle, e, fCaloModuleLabel);
  art::FindOneP<recob::MCSFitResult> tracksToMCSs(trackHandle, e, fMCSModuleLabel);
  art::FindManyP<anab::ParticleID> tracksToChi2s(trackHandle, e, fChi2ModuleLabel);
  art::FindOneP<sbn::RangeP> tracksToRangePs(trackHandle, e, fRangeModuleLabel);
  art::FindOneP<sbn::ScatterClosestApproach> tracksToClosestApproaches(trackHandle, e, fClosestApproachModuleLabel);
  art::FindOneP<sbn::StoppingChi2Fit> tracksToStoppingChi2s(trackHandle, e, fStoppingChi2ModuleLabel);

  geo::Point_t start = track->Start();
  FillElement(slcVars["slc_pfp_track_start_x"], slcCounter, pfpCounter, start.X());
  FillElement(slcVars["slc_pfp_track_start_y"], slcCounter, pfpCounter, start.Y());
  FillElement(slcVars["slc_pfp_track_start_z"], slcCounter, pfpCounter, start.Z());

  geo::Vector_t dir = track->StartDirection();
  FillElement(slcVars["slc_pfp_track_dir_x"], slcCounter, pfpCounter, dir.X());
  FillElement(slcVars["slc_pfp_track_dir_y"], slcCounter, pfpCounter, dir.Y());
  FillElement(slcVars["slc_pfp_track_dir_z"], slcCounter, pfpCounter, dir.Z());

  FillElement(slcVars["slc_pfp_track_length"], slcCounter, pfpCounter, track->Length());

  const art::Ptr<sbn::MVAPID> dazzle = tracksToDazzle.at(track.key());
  if(dazzle.isNonnull())
    ExtractDazzle(dazzle, slcCounter, pfpCounter);

  const std::vector<art::Ptr<anab::Calorimetry>> calos = tracksToCalos.at(track.key());
  const size_t maxHits = calos.size() != 3 ? -1 : std::max({calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size()});
  const int bestPlane  = calos.size() != 3 ? -1 : (calos[2]->dEdx().size() == maxHits) ? 2 : (calos[0]->dEdx().size() == maxHits) ? 0 : 
    (calos[1]->dEdx().size() == maxHits) ? 1 : -1;

  if(calos.size() == 3)
    ExtractCalo(calos[bestPlane], slcCounter, pfpCounter);

  const std::vector<art::Ptr<anab::ParticleID>> chi2s = tracksToChi2s.at(track.key());
  if(chi2s.size() == 3)
    ExtractChi2PID(chi2s[bestPlane], slcCounter, pfpCounter);

  const art::Ptr<recob::MCSFitResult> mcs = tracksToMCSs.at(track.key());
  if(mcs.isNonnull())
    ExtractMCS(mcs, slcCounter, pfpCounter);

  const art::Ptr<sbn::RangeP> rangeP = tracksToRangePs.at(track.key());
  if(rangeP.isNonnull())
    FillElement(slcVars["slc_pfp_track_range_p"], slcCounter, pfpCounter, rangeP->range_p);

  const art::Ptr<sbn::ScatterClosestApproach> closestApproach = tracksToClosestApproaches.at(track.key());
  if(closestApproach.isNonnull())
    FillElement(slcVars["slc_pfp_track_closest_approach_mean_dca"], slcCounter, pfpCounter, closestApproach->mean);

  const art::Ptr<sbn::StoppingChi2Fit> stoppingChi2 = tracksToStoppingChi2s.at(track.key());
  if(stoppingChi2.isNonnull())
    ExtractStoppingChi2(stoppingChi2, slcCounter, pfpCounter);
}

void sbnd::NCPiZeroAnalysis::ExtractDazzle(const art::Ptr<sbn::MVAPID> &dazzle, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = dazzle->mvaScoreMap;

  FillElement(slcVars["slc_pfp_track_dazzle_muon_score"], slcCounter, pfpCounter, map.at(13));
  FillElement(slcVars["slc_pfp_track_dazzle_pion_score"], slcCounter, pfpCounter, map.at(211));
  FillElement(slcVars["slc_pfp_track_dazzle_proton_score"], slcCounter, pfpCounter, map.at(2212));
  FillElement(slcVars["slc_pfp_track_dazzle_other_score"], slcCounter, pfpCounter, map.at(0));
  FillElement(slcVars["slc_pfp_track_dazzle_pdg"], slcCounter, pfpCounter, dazzle->BestPDG());
}

void sbnd::NCPiZeroAnalysis::ExtractCalo(const art::Ptr<anab::Calorimetry> &calo, const int slcCounter, const int pfpCounter)
{
  const std::vector<float> &dQdx = calo->dQdx();
  const std::vector<float> &dEdx = calo->dEdx();
  const std::vector<float> &pitch = calo->TrkPitchVec();

  float ke = 0., charge = 0.;

  for(size_t i = 0; i < dQdx.size(); ++i)
    {
      ke     += dEdx[i] * pitch[i];
      charge += dQdx[i] * pitch[i];
    }

  FillElement(slcVars["slc_pfp_track_ke"], slcCounter, pfpCounter, ke);
  FillElement(slcVars["slc_pfp_track_charge"], slcCounter, pfpCounter, charge);
}

void sbnd::NCPiZeroAnalysis::ExtractChi2PID(const art::Ptr<anab::ParticleID> &chi2pid, const int slcCounter, const int pfpCounter)
{
  const std::vector<anab::sParticleIDAlgScores> AlgScoresVec = chi2pid->ParticleIDAlgScores();
  std::vector<std::pair<int, float>> chi2s;

  for(size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
    {
      const anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

      if(AlgScore.fAlgName == "Chi2")
        {
          chi2s.push_back({AlgScore.fAssumedPdg, AlgScore.fValue});

          switch(AlgScore.fAssumedPdg)
            {
            case 13:
              FillElement(slcVars["slc_pfp_track_chi2_muon"], slcCounter, pfpCounter, AlgScore.fValue);
              break;
            case 211:
              FillElement(slcVars["slc_pfp_track_chi2_pion"], slcCounter, pfpCounter, AlgScore.fValue);
              break;
            case 321:
              FillElement(slcVars["slc_pfp_track_chi2_kaon"], slcCounter, pfpCounter, AlgScore.fValue);
              break;
            case 2212:
              FillElement(slcVars["slc_pfp_track_chi2_proton"], slcCounter, pfpCounter, AlgScore.fValue);
              break;
            }
        }
    }

  if(chi2s.size() > 0)
    {
      std::sort(chi2s.begin(), chi2s.end(),
                [](auto const& a, auto const& b)
                { return a.second < b.second; });

      FillElement(slcVars["slc_pfp_track_chi2_pdg"], slcCounter, pfpCounter, chi2s[0].first);
    }
}

void sbnd::NCPiZeroAnalysis::ExtractMCS(const art::Ptr<recob::MCSFitResult> &mcs, const int slcCounter, const int pfpCounter)
{
  if(mcs->scatterAngles().empty())
    return;

  unsigned int counter = 0;
  float maxScatter = 0.f, sumScatter = 0.f;

  for(auto const& angle : mcs->scatterAngles())
    {
      if(angle < 0)
        continue;

      maxScatter = std::max(maxScatter, angle);
      sumScatter += angle;
      counter++;
    }

  if(!counter)
    return;

  FillElement(slcVars["slc_pfp_track_mcs_mom"], slcCounter, pfpCounter, mcs->fwdMomentum());
  FillElement(slcVars["slc_pfp_track_mcs_mean_scatter"], slcCounter, pfpCounter, sumScatter / counter);
  FillElement(slcVars["slc_pfp_track_mcs_max_scatter_ratio"], slcCounter, pfpCounter, maxScatter / sumScatter);
}

void sbnd::NCPiZeroAnalysis::ExtractStoppingChi2(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2, const int slcCounter, const int pfpCounter)
{
  const float pol0Chi2 = stoppingChi2->pol0Chi2;
  const float expChi2  = stoppingChi2->expChi2;
  const float ratio    = (pol0Chi2 > 0.f && expChi2 > 0.f) ? pol0Chi2 / expChi2 : -5.f;

  FillElement(slcVars["slc_pfp_track_stopping_dedx_chi2_ratio"], slcCounter, pfpCounter, ratio);
  FillElement(slcVars["slc_pfp_track_stopping_dedx_pol0_fit"], slcCounter, pfpCounter, stoppingChi2->pol0Fit);
}

void sbnd::NCPiZeroAnalysis::AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                                           const art::Handle<std::vector<recob::Shower>> &showerHandle, const art::Ptr<recob::Vertex> &vtx,
                                           const std::vector<art::Ptr<recob::Hit>> &hits)
{
  art::FindOneP<sbn::MVAPID> showersToRazzle(showerHandle, e, fRazzleModuleLabel);
  art::FindOneP<float> showersToCosmicDist(showerHandle, e, fCosmicDistModuleLabel);
  art::FindOneP<sbn::ShowerTrackFit> showersToTrackFit(showerHandle, e, fShowerTrackFitModuleLabel);
  art::FindOneP<sbn::ShowerDensityFit> showersToDensityFit(showerHandle, e, fShowerDensityFitModuleLabel);

  geo::Point_t start(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
  FillElement(slcVars["slc_pfp_shower_start_x"], slcCounter, pfpCounter, start.X());
  FillElement(slcVars["slc_pfp_shower_start_y"], slcCounter, pfpCounter, start.Y());
  FillElement(slcVars["slc_pfp_shower_start_z"], slcCounter, pfpCounter, start.Z());

  double convGap = vtx.isNonnull() ? (start - vtx->position()).R() : def_double;
  FillElement(slcVars["slc_pfp_shower_conv_gap"], slcCounter, pfpCounter, convGap);

  geo::Vector_t dir(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
  FillElement(slcVars["slc_pfp_shower_dir_x"], slcCounter, pfpCounter, dir.X());
  FillElement(slcVars["slc_pfp_shower_dir_y"], slcCounter, pfpCounter, dir.Y());
  FillElement(slcVars["slc_pfp_shower_dir_z"], slcCounter, pfpCounter, dir.Z());

  FillElement(slcVars["slc_pfp_shower_length"], slcCounter, pfpCounter, shower->Length());
  FillElement(slcVars["slc_pfp_shower_open_angle"], slcCounter, pfpCounter, shower->OpenAngle());

  ExtractCalo(shower, slcCounter, pfpCounter, hits);

  const art::Ptr<sbn::MVAPID> razzle = showersToRazzle.at(shower.key());
  if(razzle.isNonnull())
    ExtractRazzle(razzle, slcCounter, pfpCounter);

  const art::Ptr<float> cosmicDist = showersToCosmicDist.at(shower.key());
  if(cosmicDist.isNonnull())
    FillElement(slcVars["slc_pfp_shower_cosmic_dist"], slcCounter, pfpCounter, *cosmicDist);

  const art::Ptr<sbn::ShowerTrackFit> trackFit = showersToTrackFit.at(shower.key());
  if(trackFit.isNonnull())
    {
      FillElement(slcVars["slc_pfp_shower_track_length"], slcCounter, pfpCounter, trackFit->mTrackLength);
      FillElement(slcVars["slc_pfp_shower_track_width"], slcCounter, pfpCounter, trackFit->mTrackWidth);
    }

  const art::Ptr<sbn::ShowerDensityFit> densityFit = showersToDensityFit.at(shower.key());
  if(densityFit.isNonnull())
    {
      FillElement(slcVars["slc_pfp_shower_density_grad"], slcCounter, pfpCounter, densityFit->mDensityGrad);
      FillElement(slcVars["slc_pfp_shower_density_pow"], slcCounter, pfpCounter, densityFit->mDensityPow);
    }
}

void sbnd::NCPiZeroAnalysis::ExtractRazzle(const art::Ptr<sbn::MVAPID> &razzle, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = razzle->mvaScoreMap;

  FillElement(slcVars["slc_pfp_shower_razzle_electron_score"], slcCounter, pfpCounter, map.at(11));
  FillElement(slcVars["slc_pfp_shower_razzle_photon_score"], slcCounter, pfpCounter, map.at(22));
  FillElement(slcVars["slc_pfp_shower_razzle_other_score"], slcCounter, pfpCounter, map.at(0));
  FillElement(slcVars["slc_pfp_shower_razzle_pdg"], slcCounter, pfpCounter, razzle->BestPDG());
}

void sbnd::NCPiZeroAnalysis::ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                                         const std::vector<art::Ptr<recob::Hit>> &hits)
{
  const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();

  std::array<int, 3> showerPlaneHits      = { 0, 0, 0 };
  std::array<double, 3> showerPlanePitches = { -1., -1., -1. };

  for(auto const& hit : hits)
    showerPlaneHits[hit->WireID().Plane]++;

  for(geo::PlaneGeo const& plane : geom->Iterate<geo::PlaneGeo>())
    {
      const double angleToVert = geom->WireAngleToVertical(plane.View(), plane.ID()) - 0.5 * M_PI;
      const double cosgamma    = std::abs(std::sin(angleToVert) * shower->Direction().Y() + std::cos(angleToVert) * shower->Direction().Z());

      showerPlanePitches[plane.ID().Plane] = plane.WirePitch() / cosgamma;
    }

  int bestPlane = shower->best_plane();

  FillElement(slcVars["slc_pfp_shower_energy"], slcCounter, pfpCounter, shower->Energy()[bestPlane]);
  FillElement(slcVars["slc_pfp_shower_dedx"], slcCounter, pfpCounter, shower->dEdx()[bestPlane]);

  const double length      = shower->Length();
  const double bestEnergy  = shower->Energy()[bestPlane];
  const int bestPlaneHits  = showerPlaneHits[bestPlane];
  const double bestPitch   = showerPlanePitches[bestPlane];
  const double wiresHit    = bestPitch > std::numeric_limits<double>::epsilon() ? length / bestPitch : -5.;

  FillElement(slcVars["slc_pfp_shower_sqrt_energy_density"], slcCounter, pfpCounter,
              (length > 0 && bestEnergy > 0) ? std::sqrt(bestEnergy) / length : -5.);
  FillElement(slcVars["slc_pfp_shower_modified_hit_density"], slcCounter, pfpCounter,
              wiresHit > 1. ? bestPlaneHits / wiresHit : -5.);
}

void sbnd::NCPiZeroAnalysis::ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = razzled->mvaScoreMap;

  FillElement(slcVars["slc_pfp_razzled_electron_score"], slcCounter, pfpCounter, map.at(11));
  FillElement(slcVars["slc_pfp_razzled_muon_score"], slcCounter, pfpCounter, map.at(13));
  FillElement(slcVars["slc_pfp_razzled_photon_score"], slcCounter, pfpCounter, map.at(22));
  FillElement(slcVars["slc_pfp_razzled_pion_score"], slcCounter, pfpCounter, map.at(211));
  FillElement(slcVars["slc_pfp_razzled_proton_score"], slcCounter, pfpCounter, map.at(2212));
  FillElement(slcVars["slc_pfp_razzled_pdg"], slcCounter, pfpCounter, razzled->BestPDG());
}

void sbnd::NCPiZeroAnalysis::AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, const int slcCounter,
                                               const art::Handle<std::vector<recob::Slice>> &sliceHandle)
{
  art::FindManyP<recob::Hit> slicesToHits(sliceHandle, e, fSliceModuleLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const std::vector<art::Ptr<recob::Hit>> sliceHits = slicesToHits.at(slc.key());

  std::map<int, int> objectHitMap;
  for(auto const &hit : sliceHits)
    objectHitMap[TruthMatchUtils::TrueParticleID(clockData, hit, true)]++;

  std::map<const art::Ptr<simb::MCTruth>, int> mcTruthHitMap;
  for(auto const& [trackID, nhits] : objectHitMap)
    {
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      mcTruthHitMap[mct] += nhits;
    }

  int maxHits = def_int;
  art::Ptr<simb::MCTruth> bestMCT = art::Ptr<simb::MCTruth>();

  for(auto const& [mct, nhits] : mcTruthHitMap)
    {
      if(nhits > maxHits)
        {
          maxHits = nhits;
          bestMCT  = mct;
        }
    }

  const float comp = fNuHitsMap[bestMCT] == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(fNuHitsMap[bestMCT]);
  const float pur  = sliceHits.size() == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(sliceHits.size());

  FillElement(slcVars["slc_comp"], slcCounter, comp);
  FillElement(slcVars["slc_pur"], slcCounter, pur);

  if(fBeamOff)
    FillElement(slcVars["slc_true_event_type"], slcCounter, (int) kCosmic);
  else if(bestMCT.isNonnull())
    AnalyseMCTruth(e, slcVars, bestMCT, slcCounter, "slc_true");
  else
    FillElement(slcVars["slc_true_event_type"], slcCounter, (int) kFailedTruthMatch);
}

void sbnd::NCPiZeroAnalysis::SelectSlice(const int counter)
{
  bool is_clear_cosmic;
  AccessElement(slcVars["slc_is_clear_cosmic"], counter, is_clear_cosmic);

  bool is_fv;
  AccessElement(slcVars["slc_is_fv"], counter, is_fv);

  float crumbs;
  AccessElement(slcVars["slc_crumbs_score"], counter, crumbs);
  const bool passes_crumbs = crumbs > -0.025;

  int ndazzlemuonscut;
  AccessElement(slcVars["slc_n_dazzle_muons_cut_based"], counter, ndazzlemuonscut);
  const bool passes_dazzle = ndazzlemuonscut == 0;

  int nshws;
  AccessElement(slcVars["slc_n_shws"], counter, nshws);
  const bool passes_shws = nshws > 1;

  int nrazzlephotons;
  AccessElement(slcVars["slc_n_razzle_photons"], counter, nrazzlephotons);
  const bool passes_razzle = nrazzlephotons > 1;

  int nrazzlephotonscut;
  AccessElement(slcVars["slc_n_razzle_photons_cut_based"], counter, nrazzlephotonscut);
  const bool passes_razzle_cut = nrazzlephotonscut > 0;

  const bool sel1 = is_clear_cosmic && is_fv && passes_crumbs && passes_dazzle && passes_shws && passes_razzle;
  FillElement(slcVars["slc_sel1"], counter, sel1);

  const bool sel2 = is_clear_cosmic && is_fv && passes_crumbs && passes_dazzle && passes_shws && passes_razzle_cut;
  FillElement(slcVars["slc_sel2"], counter, sel2);
}

void sbnd::NCPiZeroAnalysis::ProduceMultiSliceCandidates(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                                         const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  std::vector<bool> slc_is_clear_cosmic;
  GetVar(slcVars["slc_is_clear_cosmic"], slc_is_clear_cosmic);

  _n_multi_slice_candidates = 0;

  for(int i = 0; i < _n_slc; ++i)
    {
      if(slc_is_clear_cosmic.at(i))
        continue;

      for(int j = i + 1; j < _n_slc; ++j)
        {
          if(slc_is_clear_cosmic.at(j))
            continue;

          ++_n_multi_slice_candidates;
        }
    }

  ResizeVectors(mscVars, _n_multi_slice_candidates);

  int mscCounter = 0;

  std::vector<size_t> slc_key, slc_true_mctruth_id, slc_n_pfps;
  std::vector<int> slc_true_event_type, slc_n_razzle_electrons, slc_n_razzle_photons,
    slc_n_razzle_other, slc_n_dazzle_muons, slc_n_dazzle_pions, slc_n_dazzle_protons,
    slc_n_dazzle_other, slc_n_razzled_electrons, slc_n_razzled_muons,
    slc_n_razzled_photons, slc_n_razzled_pions, slc_n_razzled_protons;
  std::vector<float> slc_comp, slc_pur, slc_crumbs_score;
  std::vector<double> slc_vtx_x, slc_vtx_y, slc_vtx_z,
    slc_opt0_measPE, slc_opt0_hypPE;

  GetVar(slcVars["slc_key"], slc_key);
  GetVar(slcVars["slc_true_mctruth_id"], slc_true_mctruth_id);
  GetVar(slcVars["slc_true_event_type"], slc_true_event_type);
  GetVar(slcVars["slc_comp"], slc_comp);
  GetVar(slcVars["slc_pur"], slc_pur);
  GetVar(slcVars["slc_vtx_x"], slc_vtx_x);
  GetVar(slcVars["slc_vtx_y"], slc_vtx_y);
  GetVar(slcVars["slc_vtx_z"], slc_vtx_z);
  GetVar(slcVars["slc_opt0_measPE"], slc_opt0_measPE);
  GetVar(slcVars["slc_opt0_hypPE"], slc_opt0_hypPE);
  GetVar(slcVars["slc_crumbs_score"], slc_crumbs_score);
  GetVar(slcVars["slc_n_pfps"], slc_n_pfps);
  GetVar(slcVars["slc_n_razzle_electrons"], slc_n_razzle_electrons);
  GetVar(slcVars["slc_n_razzle_photons"], slc_n_razzle_photons);
  GetVar(slcVars["slc_n_razzle_other"], slc_n_razzle_other);
  GetVar(slcVars["slc_n_dazzle_muons"], slc_n_dazzle_muons);
  GetVar(slcVars["slc_n_dazzle_pions"], slc_n_dazzle_pions);
  GetVar(slcVars["slc_n_dazzle_protons"], slc_n_dazzle_protons);
  GetVar(slcVars["slc_n_dazzle_other"], slc_n_dazzle_other);
  GetVar(slcVars["slc_n_razzled_electrons"], slc_n_razzled_electrons);
  GetVar(slcVars["slc_n_razzled_muons"], slc_n_razzled_muons);
  GetVar(slcVars["slc_n_razzled_photons"], slc_n_razzled_photons);
  GetVar(slcVars["slc_n_razzled_pions"], slc_n_razzled_pions);
  GetVar(slcVars["slc_n_razzled_protons"], slc_n_razzled_protons);

  for(int i = 0; i < _n_slc; ++i)
    {
      if(slc_is_clear_cosmic.at(i))
        continue;

      for(int j = i + 1; j < _n_slc; ++j)
        {
          if(slc_is_clear_cosmic.at(j))
            continue;

          auto FillHalfCandidate = [&] (VecVarMap &vars, const int counter, const int index, const std::string prefix)
            {
              FillElement(vars[prefix + "_slc_id"], counter, index);
              FillElement(vars[prefix + "_true_mctruth_id"], counter, slc_true_mctruth_id.at(index));
              FillElement(vars[prefix + "_true_event_type"], counter, slc_true_event_type.at(index));
              FillElement(vars[prefix + "_comp"], counter, slc_comp.at(index));
              FillElement(vars[prefix + "_pur"], counter, slc_pur.at(index));
              FillElement(vars[prefix + "_vtx_x"], counter, slc_vtx_x.at(index));
              FillElement(vars[prefix + "_vtx_y"], counter, slc_vtx_y.at(index));
              FillElement(vars[prefix + "_vtx_z"], counter, slc_vtx_z.at(index));
              FillElement(vars[prefix + "_opt0_measPE"], counter, slc_opt0_measPE.at(index));
              FillElement(vars[prefix + "_opt0_hypPE"], counter, slc_opt0_hypPE.at(index));
              FillElement(vars[prefix + "_crumbs_score"], counter, slc_crumbs_score.at(index));
              FillElement(vars[prefix + "_n_pfps"], counter, slc_n_pfps.at(index));
              FillElement(vars[prefix + "_n_razzle_electrons"], counter, slc_n_razzle_electrons.at(index));
              FillElement(vars[prefix + "_n_razzle_photons"], counter, slc_n_razzle_photons.at(index));
              FillElement(vars[prefix + "_n_razzle_other"], counter, slc_n_razzle_other.at(index));
              FillElement(vars[prefix + "_n_dazzle_muons"], counter, slc_n_dazzle_muons.at(index));
              FillElement(vars[prefix + "_n_dazzle_pions"], counter, slc_n_dazzle_pions.at(index));
              FillElement(vars[prefix + "_n_dazzle_protons"], counter, slc_n_dazzle_protons.at(index));
              FillElement(vars[prefix + "_n_dazzle_other"], counter, slc_n_dazzle_other.at(index));
              FillElement(vars[prefix + "_n_razzled_electrons"], counter, slc_n_razzled_electrons.at(index));
              FillElement(vars[prefix + "_n_razzled_muons"], counter, slc_n_razzled_muons.at(index));
              FillElement(vars[prefix + "_n_razzled_photons"], counter, slc_n_razzled_photons.at(index));
              FillElement(vars[prefix + "_n_razzled_pions"], counter, slc_n_razzled_pions.at(index));
              FillElement(vars[prefix + "_n_razzled_protons"], counter, slc_n_razzled_protons.at(index));

              const bool good_opt0 = slc_opt0_measPE.at(index) > 0 && slc_opt0_hypPE.at(index) > 0;
              FillElement(vars[prefix + "_good_opt0"], counter, good_opt0);

              const double opt0_frac = (slc_opt0_hypPE.at(index) - slc_opt0_measPE.at(index)) / slc_opt0_measPE.at(index);
              FillElement(vars[prefix + "_opt0_frac"], counter, opt0_frac);
            };

          FillHalfCandidate(mscVars, mscCounter, i, "msc_0");
          FillHalfCandidate(mscVars, mscCounter, j, "msc_1");

          const bool signal = (slc_comp.at(i) + slc_comp.at(j)) > .8
            && (slc_comp.at(i) + slc_comp.at(j)) <= 1.0001
            && slc_comp.at(i) > .15 && slc_comp.at(j) > .15
            && (slc_pur.at(i) + slc_pur.at(j)) / 2. > .8
            && slc_true_event_type.at(i) == 0 && slc_true_event_type.at(j) == 0
            && slc_true_mctruth_id.at(i) == slc_true_mctruth_id.at(j);

          const double sep = TMath::Sqrt(TMath::Power(slc_vtx_x.at(i) - slc_vtx_x.at(j), 2) +
                                         TMath::Power(slc_vtx_y.at(i) - slc_vtx_y.at(j), 2) +
                                         TMath::Power(slc_vtx_z.at(i) - slc_vtx_z.at(j), 2));

          const bool matching_flash_pe = slc_opt0_measPE.at(i) == slc_opt0_measPE.at(j);

          const bool same_tpc = (slc_vtx_x.at(i) > 0 && slc_vtx_x.at(j) > 0) + (slc_vtx_x.at(i) < 0 && slc_vtx_x.at(j) < 0);

          const double sum_opt0_frac = (slc_opt0_hypPE.at(i) + slc_opt0_hypPE.at(j) - slc_opt0_measPE.at(i)) / slc_opt0_measPE.at(i);

          FillElement(mscVars["msc_signal"], mscCounter, signal);
          FillElement(mscVars["msc_sep"], mscCounter, sep);
          FillElement(mscVars["msc_matching_flash_pe"], mscCounter, matching_flash_pe);
          FillElement(mscVars["msc_same_tpc"], mscCounter, same_tpc);
          FillElement(mscVars["msc_sum_opt0_frac"], mscCounter, sum_opt0_frac);

          const double dca = DCA(e, sliceHandle, pfpHandle, slc_key.at(i), slc_key.at(j));
          FillElement(mscVars["msc_dca"], mscCounter, dca);

          ProducePiZeroCandidates(mscVars, "msc", mscCounter, { i, j });

          ++mscCounter;
        }
    }
}

void sbnd::NCPiZeroAnalysis::ProducePiZeroCandidates(VecVarMap &vars, const std::string &prefix,
                                                     const int counter, const std::vector<int> slc_ids)
{
  std::vector<int> n_razzled_photons(slc_ids.size(), 0);
  for(auto&& [i, slc_id] : enumerate(slc_ids))
    AccessElement(slcVars["slc_n_razzled_photons"], slc_id, n_razzled_photons[i]);

  const int n_photons = std::reduce(n_razzled_photons.begin(), n_razzled_photons.end());

  const size_t n_pzcs = (n_photons * (n_photons - 1)) / 2;

  FillElement(vars[prefix + "_n_pzcs"], counter, n_pzcs);

  ResizeSubVectors(vars, prefix + "_pzc", counter, n_pzcs);

  if(n_pzcs == 0)
    return;

  std::vector<int> slc_n_primary_daughters;
  std::vector<std::vector<int>> slc_pfp_razzled_pdg;
  std::vector<std::vector<double>> slc_pfp_shower_dir_x, slc_pfp_shower_dir_y, slc_pfp_shower_dir_z,
    slc_pfp_shower_energy;

  GetVar(slcVars["slc_n_primary_daughters"], slc_n_primary_daughters);
  GetVar(slcVars["slc_pfp_razzled_pdg"], slc_pfp_razzled_pdg);
  GetVar(slcVars["slc_pfp_shower_dir_x"], slc_pfp_shower_dir_x);
  GetVar(slcVars["slc_pfp_shower_dir_y"], slc_pfp_shower_dir_y);
  GetVar(slcVars["slc_pfp_shower_dir_z"], slc_pfp_shower_dir_z);
  GetVar(slcVars["slc_pfp_shower_energy"], slc_pfp_shower_energy);

  int pzcCounter = 0;

  for(auto&& [i, slc_id_a] : enumerate(slc_ids))
    {
      const int n_primary_daughters_a = slc_n_primary_daughters.at(slc_id_a);

      for(int ii = 0; ii < n_primary_daughters_a; ++ii)
        {
          if(slc_pfp_razzled_pdg.at(slc_id_a).at(ii) != 22)
            continue;

          for(auto&& [j, slc_id_b] : enumerate(slc_ids))
            {
              const int n_primary_daughters_b = slc_n_primary_daughters.at(slc_id_b);

              for(int jj = 0; jj < n_primary_daughters_b; ++jj)
                {
                  if(slc_pfp_razzled_pdg.at(slc_id_b).at(jj) != 22)
                    continue;

                  if((slc_id_a == slc_id_b && ii == jj)
                     || j < i
                     || (i == j && jj < ii))
                    continue;

                  if(slc_ids.size() > 1)
                    {
                      FillElement(vars[prefix + "_pzc_photon_0_slc_id"], counter, pzcCounter, slc_id_a);
                      FillElement(vars[prefix + "_pzc_photon_1_slc_id"], counter, pzcCounter, slc_id_b);
                    }

                  FillElement(vars[prefix + "_pzc_photon_0_id"], counter, pzcCounter, ii);
                  FillElement(vars[prefix + "_pzc_photon_1_id"], counter, pzcCounter, jj);

                  const TVector3 dir0(slc_pfp_shower_dir_x.at(slc_id_a).at(ii), slc_pfp_shower_dir_y.at(slc_id_a).at(ii), slc_pfp_shower_dir_z.at(slc_id_a).at(ii));
                  const TVector3 dir1(slc_pfp_shower_dir_x.at(slc_id_b).at(jj), slc_pfp_shower_dir_y.at(slc_id_b).at(jj), slc_pfp_shower_dir_z.at(slc_id_b).at(jj));

                  const double en0 = slc_pfp_shower_energy.at(slc_id_a).at(ii);
                  const double en1 = slc_pfp_shower_energy.at(slc_id_b).at(jj);

                  ProducePiZeroCandidate(vars, prefix, counter, pzcCounter, dir0, dir1, en0, en1);

                  ++pzcCounter;
                }
            }
        }
    }
}

void sbnd::NCPiZeroAnalysis::ProducePiZeroCandidate(VecVarMap &vars, const std::string &prefix, const int counter, const int pzcCounter,
                                                    const TVector3 &dir0, const TVector3 &dir1, const double &en0, const double &en1)
{
  const bool goodKinematics = !(dir0.X() == -999 || dir1.X() == -999 || en0 < 0 || en1 < 0);

  const double cosineThetaGammaGamma = dir0.Dot(dir1) / (dir0.Mag() * dir1.Mag());
  const TVector3 pizeroDir           = (en0 * dir0) + (en1 * dir1);

  const double invariantMass  = sqrt(2 * en0 * en1 * (1 - cosineThetaGammaGamma));
  const double pizeroMom      = pizeroDir.Mag();
  const double pizeroCosTheta = pizeroDir.Z() / pizeroMom;
  const double cosCOM         = std::abs(en0 - en1) / pizeroMom;
  const double decayAsym      = std::abs(en0 - en1) / (en0 + en1);

  FillElement(vars[prefix + "_pzc_good_kinematics"], counter, pzcCounter, goodKinematics);
  FillElement(vars[prefix + "_pzc_invariant_mass"], counter, pzcCounter, invariantMass);
  FillElement(vars[prefix + "_pzc_pizero_mom"], counter, pzcCounter, pizeroMom);
  FillElement(vars[prefix + "_pzc_cos_theta_pizero"], counter, pzcCounter, pizeroCosTheta);
  FillElement(vars[prefix + "_pzc_cos_com"], counter, pzcCounter, cosCOM);
  FillElement(vars[prefix + "_pzc_decay_asymmetry"], counter, pzcCounter, decayAsym);
}

float sbnd::NCPiZeroAnalysis::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

float sbnd::NCPiZeroAnalysis::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (fHitsMap[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(fHitsMap[trackID]);
}

double sbnd::NCPiZeroAnalysis::DCA(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const size_t key0, const size_t key1)
{
  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindManyP<recob::SpacePoint> pfpsToSpacePoints(pfpHandle, e, fSpacePointModuleLabel);

  const std::vector<art::Ptr<recob::PFParticle>> pfps0 = slicesToPFPs.at(key0);
  const std::vector<art::Ptr<recob::PFParticle>> pfps1 = slicesToPFPs.at(key1);

  std::vector<art::Ptr<recob::SpacePoint>> spacepoints0;
  std::vector<art::Ptr<recob::SpacePoint>> spacepoints1;

  for(auto const& pfp : pfps0)
    {
      const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpsToSpacePoints.at(pfp.key());
      spacepoints0.insert(spacepoints0.end(), spacepoints.begin(), spacepoints.end());
    }

  for(auto const& pfp : pfps1)
    {
      const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpsToSpacePoints.at(pfp.key());
      spacepoints1.insert(spacepoints1.end(), spacepoints.begin(), spacepoints.end());
    }

  double dca = std::numeric_limits<double>::max();

  for(auto const& sp0 : spacepoints0)
    {
      const geo::Point_t pos0 = sp0->position();

      for(auto const& sp1 : spacepoints1)
        {
          const geo::Point_t pos1 = sp1->position();

          const double dist = TMath::Sqrt(TMath::Power(pos1.X() - pos0.X(), 2) +
                                          TMath::Power(pos1.Y() - pos0.Y(), 2) +
                                          TMath::Power(pos1.Z() - pos0.Z(), 2));

          if(dist < dca)
            dca = dist;
        }
    }

  return dca;
}


bool sbnd::NCPiZeroAnalysis::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
}

bool sbnd::NCPiZeroAnalysis::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

void sbnd::NCPiZeroAnalysis::ResetSubRunVars()
{
  _pot = 0.; _spills = 0; _ngenevts = 0;
}

void sbnd::NCPiZeroAnalysis::ResetEventVars()
{
  _run = -1; _subrun = -1; _event  = -1;

  _n_nu = 0; _n_slc  = 0;
}

int sbnd::NCPiZeroAnalysis::GetTotalGenEvents(const art::Event &e)
{
  int nGenEvt = 0;
  for (const art::ProcessConfiguration &process: e.processHistory()) {
    std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
    if (genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") && genConfig->has_key("source.module_type") ) {
      int maxEvents = genConfig->get<int>("source.maxEvents");
      std::string moduleType = genConfig->get<std::string>("source.module_type");
      if (moduleType == "EmptyEvent") {
        nGenEvt += maxEvents;
      }
    }
  }

  return nGenEvt;
}

void sbnd::NCPiZeroAnalysis::ResizeVectors(VecVarMap &map, const int size)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kOneD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVar<bool>*>(var)->Assign(size, false);
              break;
            case kInt:
              dynamic_cast<InhVecVar<int>*>(var)->Assign(size, def_int);
              break;
            case kUInt:
              dynamic_cast<InhVecVar<size_t>*>(var)->Assign(size, def_size);
              break;
            case kFloat:
              dynamic_cast<InhVecVar<float>*>(var)->Assign(size, def_float);
              break;
            case kDouble:
              dynamic_cast<InhVecVar<double>*>(var)->Assign(size, def_double);
              break;
            case kUnknownVar:
              break;
            }
        }
      else if(var->IdentifyVec() == kTwoD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVecVar<bool>*>(var)->Resize(size);
              break;
            case kInt:
              dynamic_cast<InhVecVecVar<int>*>(var)->Resize(size);
              break;
            case kUInt:
              dynamic_cast<InhVecVecVar<size_t>*>(var)->Resize(size);
              break;
            case kFloat:
              dynamic_cast<InhVecVecVar<float>*>(var)->Resize(size);
              break;
            case kDouble:
              dynamic_cast<InhVecVecVar<double>*>(var)->Resize(size);
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

void sbnd::NCPiZeroAnalysis::ResizeSubVectors(VecVarMap &map, const std::string &subname, const int pos, const int size)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kTwoD && var->Name().find(subname) != std::string::npos)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVecVar<bool>*>(var)->Assign(pos, size, false);
              break;
            case kInt:
              dynamic_cast<InhVecVecVar<int>*>(var)->Assign(pos, size, def_int);
              break;
            case kUInt:
              dynamic_cast<InhVecVecVar<size_t>*>(var)->Assign(pos, size, def_size);
              break;
            case kFloat:
              dynamic_cast<InhVecVecVar<float>*>(var)->Assign(pos, size, def_float);
              break;
            case kDouble:
              dynamic_cast<InhVecVecVar<double>*>(var)->Assign(pos, size, def_double);
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

template<typename T>
void sbnd::NCPiZeroAnalysis::FillElement(VecVar *vec, const int pos, const T value)
{
  dynamic_cast<InhVecVar<T>*>(vec)->SetVal(pos, value);
}

template<typename T>
void sbnd::NCPiZeroAnalysis::FillElement(VecVar *vec, const int posA, const int posB, const T value)
{
  dynamic_cast<InhVecVecVar<T>*>(vec)->SetVal(posA, posB, value);
}

art::Ptr<recob::PFParticle> sbnd::NCPiZeroAnalysis::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
{
  for(auto pfp : pfps)
    if(pfp->IsPrimary())
      return pfp;

  return art::Ptr<recob::PFParticle>();
}

template<typename T>
void sbnd::NCPiZeroAnalysis::AccessElement(VecVar *vec, const int pos, T &value)
{
  value = dynamic_cast<InhVecVar<T>*>(vec)->GetVal(pos);
}

template<typename T>
void sbnd::NCPiZeroAnalysis::AccessElement(VecVar *vec, const int posA, const int posB, T &value)
{
  value = dynamic_cast<InhVecVecVar<T>*>(vec)->GetVal(posA, posB);
}

template<typename T>
void sbnd::NCPiZeroAnalysis::GetVar(VecVar *vec, std::vector<T> &var)
{
  var = dynamic_cast<InhVecVar<T>*>(vec)->Var();
}

template<typename T>
void sbnd::NCPiZeroAnalysis::GetVar(VecVar *vec, std::vector<std::vector<T>> &var)
{
  var = dynamic_cast<InhVecVecVar<T>*>(vec)->Var();
}

DEFINE_ART_MODULE(sbnd::NCPiZeroAnalysis)
