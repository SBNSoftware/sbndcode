////////////////////////////////////////////////////////////////////////
// Class:       DPreco
// Module Type: analyzer
// File:        DPreco_module.cc
//
// Generated at Wed Jul 17 13:53:28 2024 by Rohan Rajagopalan using artmod
// from cetpkgsupport v1_14_01.
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
#include "TProfile.h"
#include "TFile.h"
#include "TSpline.h"

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
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

#include "DPrecoStructs.h"

#include <numeric>



constexpr int def_int       = std::numeric_limits<int>::min();
constexpr size_t def_size   = std::numeric_limits<size_t>::max();
constexpr float def_float   = -std::numeric_limits<float>::max();
constexpr double def_double = -std::numeric_limits<double>::max();



namespace sbnd {
  class DPreco;
}

class sbnd::DPreco : public art::EDAnalyzer {
public:
  explicit DPreco(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DPreco(DPreco const&) = delete;
  DPreco(DPreco&&) = delete;
  DPreco& operator=(DPreco const&) = delete;
  DPreco& operator=(DPreco&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;
  virtual void beginSubRun(const art::SubRun& sr);
  virtual void endSubRun(const art::SubRun& sr);

  void ResetSubRunVars();
  void ResetEventVars();
  void ClearMaps();

  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                 const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  void SetupBranches(VecVarMap &map);

  void ResizeVectors(VecVarMap &map, const int size);
  void ResizeSubVectors(VecVarMap &map, const std::string &subname, const int pos, const int size);

  

  void AnalyseMCTruth(const art::Event &e, VecVarMap &vars, const art::Ptr<simb::MCTruth> &mct, const int counter, const std::string prefix);



  void AnalysePhotonReco(const std::string name, const int trackid, const int nuCounter, const int pzCounter);

  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                     const art::Handle<std::vector<recob::Track>> &trackHandle,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle,
                     const art::Handle<std::vector<recob::Hit>> &hitHandle);

  void AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const std::vector<art::Ptr<recob::PFParticle>> &pfps,
                   const std::vector<art::Ptr<recob::Hit>> &sliceHits, const art::Ptr<recob::Vertex> &vtx, const int slcCounter,
                   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                   const art::Handle<std::vector<recob::Shower>> &showerHandle);
  void AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, const int slcCounter, const int pfpCounter,
                    const art::Handle<std::vector<recob::Track>> &trackHandle);
  void AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle, const art::Ptr<recob::Vertex> &vtx,
                     const std::vector<art::Ptr<recob::Hit>> &hits);
  void AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, const int slcCounter,
                         const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                         const art::Handle<std::vector<recob::Hit>> &hitHandle);

  
  void ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                   const std::vector<art::Ptr<recob::Hit>> &hits);
  void ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter);
  
  

  void SelectSlice(const int counter);

  
  float Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  float Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  int GetTotalGenEvents(const art::Event &e);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

  double CorrectEnergy(const double &energy);

  template<typename T>
  void FillElement(VecVar *vec, const int pos, const T value);
  template<typename T>
  void FillElement(VecVar *vec, const int posA, const int posB, const T value);
  template<typename T>
  void FillElement(VecVar *vec, const int pos, const std::vector<T> &value);

  void TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB, const int posA, const int posB);
  void TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB,
                       const int posA0, const int posA1, const int posB0, const int posB1);
  void TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB,
                       const int posA, const int posB0, const int posB1);

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
    fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fCRUMBSModuleLabel, fDazzleModuleLabel,
    fCaloModuleLabel, fMCSModuleLabel, fChi2ModuleLabel, fRangeModuleLabel, fClosestApproachModuleLabel,
    fStoppingChi2ModuleLabel, fRazzleModuleLabel, fCosmicDistModuleLabel, fShowerTrackFitModuleLabel,
    fShowerDensityFitModuleLabel, fPOTModuleLabel, fOpT0ModuleLabel, fRazzledModuleLabel, fSpacePointModuleLabel;
  std::vector<art::InputTag> fEventWeightModuleLabels;
  bool fDebug, fBeamOff;
  

  

  std::string fShowerEnergyCorrectionFileName;
  TProfile* fShowerEnergyCorrectionHist;

  std::map<int, int> fHitsMap;
  std::map<const art::Ptr<simb::MCTruth>, int> fNuHitsMap, fNuHitsMapSP;
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
    { "nu_event_type_incl", new InhVecVar<int>("nu_event_type_incl") },
    { "nu_event_type_0p0pi", new InhVecVar<int>("nu_event_type_0p0pi") },
    { "nu_event_type_1p0pi", new InhVecVar<int>("nu_event_type_1p0pi") },
    { "nu_event_type_Np0pi", new InhVecVar<int>("nu_event_type_Np0pi") },
    { "nu_event_type_Xp0pi", new InhVecVar<int>("nu_event_type_Xp0pi") },
    { "nu_event_type_cc", new InhVecVar<int>("nu_event_type_cc") },
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
    { "nu_n_dalitz_neutral_pions", new InhVecVar<int>("nu_n_dalitz_neutral_pions") },
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
    { "nu_pz_gamma0_trackid", new InhVecVecVar<int>("nu_pz_gamma0_trackid") },
    { "nu_pz_gamma0_n_hits", new InhVecVecVar<int>("nu_pz_gamma0_n_hits") },
    { "nu_pz_gamma0_energy", new InhVecVecVar<double>("nu_pz_gamma0_energy") },
    { "nu_pz_gamma0_dir_x", new InhVecVecVar<double>("nu_pz_gamma0_dir_x") },
    { "nu_pz_gamma0_dir_y", new InhVecVecVar<double>("nu_pz_gamma0_dir_y") },
    { "nu_pz_gamma0_dir_z", new InhVecVecVar<double>("nu_pz_gamma0_dir_z") },
    { "nu_pz_gamma0_best_pfp_comp", new InhVecVecVar<float>("nu_pz_gamma0_best_pfp_comp") },
    { "nu_pz_gamma0_best_pfp_pur", new InhVecVecVar<float>("nu_pz_gamma0_best_pfp_pur") },
    { "nu_pz_gamma0_best_pfp_pdg", new InhVecVecVar<int>("nu_pz_gamma0_best_pfp_pdg") },
    { "nu_pz_gamma0_best_pfp_shower_razzle_pdg", new InhVecVecVar<int>("nu_pz_gamma0_best_pfp_shower_razzle_pdg") },
    { "nu_pz_gamma0_best_pfp_razzled_pdg", new InhVecVecVar<int>("nu_pz_gamma0_best_pfp_razzled_pdg") },
    { "nu_pz_gamma1_trackid", new InhVecVecVar<int>("nu_pz_gamma1_trackid") },
    { "nu_pz_gamma1_n_hits", new InhVecVecVar<int>("nu_pz_gamma1_n_hits") },
    { "nu_pz_gamma1_energy", new InhVecVecVar<double>("nu_pz_gamma1_energy") },
    { "nu_pz_gamma1_dir_x", new InhVecVecVar<double>("nu_pz_gamma1_dir_x") },
    { "nu_pz_gamma1_dir_y", new InhVecVecVar<double>("nu_pz_gamma1_dir_y") },
    { "nu_pz_gamma1_dir_z", new InhVecVecVar<double>("nu_pz_gamma1_dir_z") },
    { "nu_pz_gamma1_best_pfp_comp", new InhVecVecVar<float>("nu_pz_gamma1_best_pfp_comp") },
    { "nu_pz_gamma1_best_pfp_pur", new InhVecVecVar<float>("nu_pz_gamma1_best_pfp_pur") },
    { "nu_pz_gamma1_best_pfp_pdg", new InhVecVecVar<int>("nu_pz_gamma1_best_pfp_pdg") },
    { "nu_pz_gamma1_best_pfp_shower_razzle_pdg", new InhVecVecVar<int>("nu_pz_gamma1_best_pfp_shower_razzle_pdg") },
    { "nu_pz_gamma1_best_pfp_razzled_pdg", new InhVecVecVar<int>("nu_pz_gamma1_best_pfp_razzled_pdg") },
    { "nu_pz_open_angle", new InhVecVecVar<double>("nu_pz_open_angle")},
    { "nu_best_slc_comp", new InhVecVar<float>("nu_best_slc_comp") },
    { "nu_best_slc_pur", new InhVecVar<float>("nu_best_slc_pur") },
    { "nu_best_slc_is_clear_cosmic", new InhVecVar<bool>("nu_best_slc_is_clear_cosmic") },
    { "nu_best_slc_n_pfps", new InhVecVar<size_t>("nu_best_slc_n_pfps") },
    { "nu_best_slc_n_dazzle_muons", new InhVecVar<int>("nu_best_slc_n_dazzle_muons") },
    { "nu_best_slc_n_dazzle_pions", new InhVecVar<int>("nu_best_slc_n_dazzle_pions") },
    { "nu_best_slc_n_dazzle_pions_thresh", new InhVecVar<int>("nu_best_slc_n_dazzle_pions_thresh") },
    { "nu_best_slc_n_dazzle_protons", new InhVecVar<int>("nu_best_slc_n_dazzle_protons") },
    { "nu_best_slc_n_dazzle_protons_thresh", new InhVecVar<int>("nu_best_slc_n_dazzle_protons_thresh") },
    { "nu_best_slc_n_dazzle_other", new InhVecVar<int>("nu_best_slc_n_dazzle_other") },
    { "nu_best_slc_n_razzle_electrons", new InhVecVar<int>("nu_best_slc_n_razzle_electrons") },
    { "nu_best_slc_n_razzle_photons", new InhVecVar<int>("nu_best_slc_n_razzle_photons") },
    { "nu_best_slc_n_razzle_other", new InhVecVar<int>("nu_best_slc_n_razzle_other") },
    { "nu_best_slc_n_razzled_electrons", new InhVecVar<int>("nu_best_slc_n_razzled_electrons") },
    { "nu_best_slc_n_razzled_muons", new InhVecVar<int>("nu_best_slc_n_razzled_muons") },
    { "nu_best_slc_n_razzled_photons", new InhVecVar<int>("nu_best_slc_n_razzled_photons") },
    { "nu_best_slc_n_razzled_pions", new InhVecVar<int>("nu_best_slc_n_razzled_pions") },
    { "nu_best_slc_n_razzled_pions_thresh", new InhVecVar<int>("nu_best_slc_n_razzled_pions_thresh") },
    { "nu_best_slc_n_razzled_protons", new InhVecVar<int>("nu_best_slc_n_razzled_protons") },
    { "nu_best_slc_n_razzled_protons_thresh", new InhVecVar<int>("nu_best_slc_n_razzled_protons_thresh") },
    { "nu_best_slc_is_fv", new InhVecVar<bool>("nu_best_slc_is_fv") },
    { "nu_best_slc_crumbs_score", new InhVecVar<float>("nu_best_slc_crumbs_score") },
    { "nu_best_slc_crumbs_nc_score", new InhVecVar<float>("nu_best_slc_crumbs_nc_score") },
    { "nu_best_slc_crumbs_ccnue_score", new InhVecVar<float>("nu_best_slc_crumbs_ccnue_score") },
    { "nu_best_slc_crumbs_ccnumu_score", new InhVecVar<float>("nu_best_slc_crumbs_ccnumu_score") },
    { "nu_best_slc_best_pzc_invariant_mass", new InhVecVar<double>("nu_best_slc_best_pzc_invariant_mass") },
    { "nu_best_slc_best_pzc_pizero_mom", new InhVecVar<double>("nu_best_slc_best_pzc_pizero_mom") },
    { "nu_best_slc_best_pzc_cos_theta_pizero", new InhVecVar<double>("nu_best_slc_best_pzc_cos_theta_pizero") },
    { "nu_best_slc_best_pzc_cos_com", new InhVecVar<double>("nu_best_slc_best_pzc_cos_com") },
    { "nu_best_slc_best_pzc_decay_asymmetry", new InhVecVar<double>("nu_best_slc_best_pzc_decay_asymmetry") },
    { "nu_best_slc_best_pzc_invariant_mass_corr", new InhVecVar<double>("nu_best_slc_best_pzc_invariant_mass_corr") },
    { "nu_best_slc_best_pzc_pizero_mom_corr", new InhVecVar<double>("nu_best_slc_best_pzc_pizero_mom_corr") },
    { "nu_best_slc_best_pzc_cos_theta_pizero_corr", new InhVecVar<double>("nu_best_slc_best_pzc_cos_theta_pizero_corr") },
    { "nu_best_slc_best_pzc_cos_com_corr", new InhVecVar<double>("nu_best_slc_best_pzc_cos_com_corr") },
    { "nu_best_slc_best_pzc_decay_asymmetry_corr", new InhVecVar<double>("nu_best_slc_best_pzc_decay_asymmetry_corr") },
    { "nu_best_slc_best_corr_pzc_invariant_mass", new InhVecVar<double>("nu_best_slc_best_corr_pzc_invariant_mass") },
    { "nu_best_slc_best_corr_pzc_pizero_mom", new InhVecVar<double>("nu_best_slc_best_corr_pzc_pizero_mom") },
    { "nu_best_slc_best_corr_pzc_cos_theta_pizero", new InhVecVar<double>("nu_best_slc_best_corr_pzc_cos_theta_pizero") },
    { "nu_best_slc_best_corr_pzc_cos_com", new InhVecVar<double>("nu_best_slc_best_corr_pzc_cos_com") },
    { "nu_best_slc_best_corr_pzc_decay_asymmetry", new InhVecVar<double>("nu_best_slc_best_corr_pzc_decay_asymmetry") },
    { "nu_best_slc_best_corr_pzc_invariant_mass_corr", new InhVecVar<double>("nu_best_slc_best_corr_pzc_invariant_mass_corr") },
    { "nu_best_slc_best_corr_pzc_pizero_mom_corr", new InhVecVar<double>("nu_best_slc_best_corr_pzc_pizero_mom_corr") },
    { "nu_best_slc_best_corr_pzc_cos_theta_pizero_corr", new InhVecVar<double>("nu_best_slc_best_corr_pzc_cos_theta_pizero_corr") },
    { "nu_best_slc_best_corr_pzc_cos_com_corr", new InhVecVar<double>("nu_best_slc_best_corr_pzc_cos_com_corr") },
    { "nu_best_slc_best_corr_pzc_decay_asymmetry_corr", new InhVecVar<double>("nu_best_slc_best_corr_pzc_decay_asymmetry_corr") },
  };

  int _n_slc;

  VecVarMap slcVars = {
    { "slc_key", new InhVecVar<size_t>("slc_key") },
    { "slc_n_hits", new InhVecVar<size_t>("slc_n_hits") },
    { "slc_n_used_hits", new InhVecVar<size_t>("slc_n_used_hits") },
    { "slc_n_pfps", new InhVecVar<size_t>("slc_n_pfps") },
    { "slc_primary_pfp_id", new InhVecVar<size_t>("slc_primary_pfp_id") },
    { "slc_primary_pfp_pdg", new InhVecVar<int>("slc_primary_pfp_pdg") },
    { "slc_is_clear_cosmic", new InhVecVar<bool>("slc_is_clear_cosmic") },
    { "slc_n_primary_children", new InhVecVar<int>("slc_n_primary_children") },
    { "slc_n_trks", new InhVecVar<int>("slc_n_trks") },
    { "slc_n_shws", new InhVecVar<int>("slc_n_shws") },
    { "slc_all_trks_contained", new InhVecVar<bool>("slc_all_trks_contained") },
    { "slc_all_shws_contained", new InhVecVar<bool>("slc_all_shws_contained") },
    { "slc_all_other_trks_contained", new InhVecVar<bool>("slc_all_other_trks_contained") },
    { "slc_all_other_shws_contained", new InhVecVar<bool>("slc_all_other_shws_contained") },
    { "slc_n_dazzle_muons", new InhVecVar<int>("slc_n_dazzle_muons") },
    { "slc_n_dazzle_pions", new InhVecVar<int>("slc_n_dazzle_pions") },
    { "slc_n_dazzle_pions_thresh", new InhVecVar<int>("slc_n_dazzle_pions_thresh") },
    { "slc_n_dazzle_protons", new InhVecVar<int>("slc_n_dazzle_protons") },
    { "slc_n_dazzle_protons_thresh", new InhVecVar<int>("slc_n_dazzle_protons_thresh") },
    { "slc_n_dazzle_other", new InhVecVar<int>("slc_n_dazzle_other") },
    { "slc_n_razzle_electrons", new InhVecVar<int>("slc_n_razzle_electrons") },
    { "slc_n_razzle_photons", new InhVecVar<int>("slc_n_razzle_photons") },
    { "slc_n_razzle_other", new InhVecVar<int>("slc_n_razzle_other") },
    { "slc_n_razzled_electrons", new InhVecVar<int>("slc_n_razzled_electrons") },
    { "slc_n_razzled_muons", new InhVecVar<int>("slc_n_razzled_muons") },
    { "slc_n_razzled_photons", new InhVecVar<int>("slc_n_razzled_photons") },
    { "slc_n_razzled_pions", new InhVecVar<int>("slc_n_razzled_pions") },
    { "slc_n_razzled_pions_thresh", new InhVecVar<int>("slc_n_razzled_pions_thresh") },
    { "slc_n_razzled_protons", new InhVecVar<int>("slc_n_razzled_protons") },
    { "slc_n_razzled_protons_thresh", new InhVecVar<int>("slc_n_razzled_protons_thresh") },
    { "slc_n_primary_trks", new InhVecVar<int>("slc_n_primary_trks") },
    { "slc_n_primary_shws", new InhVecVar<int>("slc_n_primary_shws") },
    { "slc_all_primary_trks_contained", new InhVecVar<bool>("slc_all_primary_trks_contained") },
    { "slc_all_primary_shws_contained", new InhVecVar<bool>("slc_all_primary_shws_contained") },
    { "slc_all_other_primary_trks_contained", new InhVecVar<bool>("slc_all_other_primary_trks_contained") },
    { "slc_all_other_primary_shws_contained", new InhVecVar<bool>("slc_all_other_primary_shws_contained") },
    { "slc_n_primary_dazzle_muons", new InhVecVar<int>("slc_n_primary_dazzle_muons") },
    { "slc_n_primary_dazzle_pions", new InhVecVar<int>("slc_n_primary_dazzle_pions") },
    { "slc_n_primary_dazzle_pions_thresh", new InhVecVar<int>("slc_n_primary_dazzle_pions_thresh") },
    { "slc_n_primary_dazzle_protons", new InhVecVar<int>("slc_n_primary_dazzle_protons") },
    { "slc_n_primary_dazzle_protons_thresh", new InhVecVar<int>("slc_n_primary_dazzle_protons_thresh") },
    { "slc_n_primary_dazzle_other", new InhVecVar<int>("slc_n_primary_dazzle_other") },
    { "slc_n_primary_razzle_electrons", new InhVecVar<int>("slc_n_primary_razzle_electrons") },
    { "slc_n_primary_razzle_photons", new InhVecVar<int>("slc_n_primary_razzle_photons") },
    { "slc_n_primary_razzle_other", new InhVecVar<int>("slc_n_primary_razzle_other") },
    { "slc_n_primary_razzled_electrons", new InhVecVar<int>("slc_n_primary_razzled_electrons") },
    { "slc_n_primary_razzled_muons", new InhVecVar<int>("slc_n_primary_razzled_muons") },
    { "slc_n_primary_razzled_photons", new InhVecVar<int>("slc_n_primary_razzled_photons") },
    { "slc_n_primary_razzled_pions", new InhVecVar<int>("slc_n_primary_razzled_pions") },
    { "slc_n_primary_razzled_pions_thresh", new InhVecVar<int>("slc_n_primary_razzled_pions_thresh") },
    { "slc_n_primary_razzled_protons", new InhVecVar<int>("slc_n_primary_razzled_protons") },
    { "slc_n_primary_razzled_protons_thresh", new InhVecVar<int>("slc_n_primary_razzled_protons_thresh") },
    { "slc_true_mctruth_id", new InhVecVar<size_t>("slc_true_mctruth_id") },
    { "slc_true_event_type_incl", new InhVecVar<int>("slc_true_event_type_incl") },
    { "slc_true_event_type_0p0pi", new InhVecVar<int>("slc_true_event_type_0p0pi") },
    { "slc_true_event_type_1p0pi", new InhVecVar<int>("slc_true_event_type_1p0pi") },
    { "slc_true_event_type_Np0pi", new InhVecVar<int>("slc_true_event_type_Np0pi") },
    { "slc_true_event_type_Xp0pi", new InhVecVar<int>("slc_true_event_type_Xp0pi") },
    { "slc_true_event_type_cc", new InhVecVar<int>("slc_true_event_type_cc") },
    { "slc_true_signal", new InhVecVar<bool>("slc_true_signal") },
    { "slc_comp", new InhVecVar<float>("slc_comp") },
    { "slc_comp_sp_only", new InhVecVar<float>("slc_comp_sp_only") },
    { "slc_pur", new InhVecVar<float>("slc_pur") },
    { "slc_pur_sp_only", new InhVecVar<float>("slc_pur_sp_only") },
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
    { "slc_true_n_dalitz_neutral_pions", new InhVecVar<int>("slc_true_n_dalitz_neutral_pions") },
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
    { "slc_true_pz_gamma0_trackid", new InhVecVecVar<int>("slc_true_pz_gamma0_trackid") },
    { "slc_true_pz_gamma0_n_hits", new InhVecVecVar<int>("slc_true_pz_gamma0_n_hits") },
    { "slc_true_pz_gamma0_energy", new InhVecVecVar<double>("slc_true_pz_gamma0_energy") },
    { "slc_true_pz_gamma0_dir_x", new InhVecVecVar<double>("slc_true_pz_gamma0_dir_x") },
    { "slc_true_pz_gamma0_dir_y", new InhVecVecVar<double>("slc_true_pz_gamma0_dir_y") },
    { "slc_true_pz_gamma0_dir_z", new InhVecVecVar<double>("slc_true_pz_gamma0_dir_z") },
    { "slc_true_pz_gamma1_trackid", new InhVecVecVar<int>("slc_true_pz_gamma1_trackid") },
    { "slc_true_pz_gamma1_n_hits", new InhVecVecVar<int>("slc_true_pz_gamma1_n_hits") },
    { "slc_true_pz_gamma1_energy", new InhVecVecVar<double>("slc_true_pz_gamma1_energy") },
    { "slc_true_pz_gamma1_dir_x", new InhVecVecVar<double>("slc_true_pz_gamma1_dir_x") },
    { "slc_true_pz_gamma1_dir_y", new InhVecVecVar<double>("slc_true_pz_gamma1_dir_y") },
    { "slc_true_pz_gamma1_dir_z", new InhVecVecVar<double>("slc_true_pz_gamma1_dir_z") },
    { "slc_true_pz_open_angle", new InhVecVecVar<double>("slc_true_pz_open_angle")},
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
    { "slc_opt0_fracPE", new InhVecVar<double>("slc_opt0_fracPE") },
    { "slc_pfp_id", new InhVecVecVar<size_t>("slc_pfp_id") },
    { "slc_pfp_primary", new InhVecVecVar<bool>("slc_pfp_primary") },
    { "slc_pfp_primary_child", new InhVecVecVar<bool>("slc_pfp_primary_child") },
    { "slc_pfp_pdg", new InhVecVecVar<int>("slc_pfp_pdg") },
    { "slc_pfp_track_score", new InhVecVecVar<float>("slc_pfp_track_score") },
    { "slc_pfp_n_children", new InhVecVecVar<int>("slc_pfp_n_children") },
    { "slc_pfp_good_track", new InhVecVecVar<bool>("slc_pfp_good_track") },
    { "slc_pfp_good_shower", new InhVecVecVar<bool>("slc_pfp_good_shower") },
    { "slc_pfp_true_trackid", new InhVecVecVar<int>("slc_pfp_true_trackid") },
    { "slc_pfp_true_pdg", new InhVecVecVar<int>("slc_pfp_true_pdg") },
    { "slc_pfp_true_energy", new InhVecVecVar<double>("slc_pfp_true_energy") },
    { "slc_pfp_true_ke", new InhVecVecVar<double>("slc_pfp_true_ke") },
    { "slc_pfp_true_p_x", new InhVecVecVar<double>("slc_pfp_true_p_x") },
    { "slc_pfp_true_p_y", new InhVecVecVar<double>("slc_pfp_true_p_y") },
    { "slc_pfp_true_p_z", new InhVecVecVar<double>("slc_pfp_true_p_z") },
    { "slc_pfp_comp", new InhVecVecVar<float>("slc_pfp_comp") },
    { "slc_pfp_pur", new InhVecVecVar<float>("slc_pfp_pur") },
    { "slc_pfp_n_sps", new InhVecVecVar<size_t>("slc_pfp_n_sps") },
    { "slc_pfp_n_hits", new InhVecVecVar<size_t>("slc_pfp_n_hits") },
    { "slc_pfp_track_start_x", new InhVecVecVar<double>("slc_pfp_track_start_x") },
    { "slc_pfp_track_start_y", new InhVecVecVar<double>("slc_pfp_track_start_y") },
    { "slc_pfp_track_start_z", new InhVecVecVar<double>("slc_pfp_track_start_z") },
    { "slc_pfp_track_dir_x", new InhVecVecVar<double>("slc_pfp_track_dir_x") },
    { "slc_pfp_track_dir_y", new InhVecVecVar<double>("slc_pfp_track_dir_y") },
    { "slc_pfp_track_dir_z", new InhVecVecVar<double>("slc_pfp_track_dir_z") },
    { "slc_pfp_track_contained", new InhVecVecVar<bool>("slc_pfp_track_contained") },
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
    { "slc_pfp_shower_contained", new InhVecVecVar<bool>("slc_pfp_shower_contained") },
    { "slc_pfp_shower_length", new InhVecVecVar<double>("slc_pfp_shower_length") },
    { "slc_pfp_shower_open_angle", new InhVecVecVar<double>("slc_pfp_shower_open_angle") },
    { "slc_pfp_shower_energy", new InhVecVecVar<double>("slc_pfp_shower_energy") },
    { "slc_pfp_shower_energy_coll", new InhVecVecVar<double>("slc_pfp_shower_energy_coll") },
    { "slc_pfp_shower_energy_corr", new InhVecVecVar<double>("slc_pfp_shower_energy_corr") },
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
    { "slc_pzc_photon_0_id", new InhVecVecVar<size_t>("slc_pzc_photon_0_id") },
    { "slc_pzc_photon_1_id", new InhVecVecVar<size_t>("slc_pzc_photon_1_id") },
    { "slc_pzc_good_kinematics", new InhVecVecVar<bool>("slc_pzc_good_kinematics") },
    { "slc_pzc_invariant_mass", new InhVecVecVar<double>("slc_pzc_invariant_mass") },
    { "slc_pzc_pizero_mom", new InhVecVecVar<double>("slc_pzc_pizero_mom") },
    { "slc_pzc_cos_theta_pizero", new InhVecVecVar<double>("slc_pzc_cos_theta_pizero") },
    { "slc_pzc_cos_com", new InhVecVecVar<double>("slc_pzc_cos_com") },
    { "slc_pzc_decay_asymmetry", new InhVecVecVar<double>("slc_pzc_decay_asymmetry") },
    { "slc_pzc_good_kinematics_corr", new InhVecVecVar<bool>("slc_pzc_good_kinematics_corr") },
    { "slc_pzc_invariant_mass_corr", new InhVecVecVar<double>("slc_pzc_invariant_mass_corr") },
    { "slc_pzc_pizero_mom_corr", new InhVecVecVar<double>("slc_pzc_pizero_mom_corr") },
    { "slc_pzc_cos_theta_pizero_corr", new InhVecVecVar<double>("slc_pzc_cos_theta_pizero_corr") },
    { "slc_pzc_cos_com_corr", new InhVecVecVar<double>("slc_pzc_cos_com_corr") },
    { "slc_pzc_decay_asymmetry_corr", new InhVecVecVar<double>("slc_pzc_decay_asymmetry_corr") },
    { "slc_pzc_photon_0_true_trackid", new InhVecVecVar<int>("slc_pzc_photon_0_true_trackid") },
    { "slc_pzc_photon_0_true_pdg", new InhVecVecVar<int>("slc_pzc_photon_0_true_pdg") },
    { "slc_pzc_photon_0_comp", new InhVecVecVar<float>("slc_pzc_photon_0_comp") },
    { "slc_pzc_photon_0_pur", new InhVecVecVar<float>("slc_pzc_photon_0_pur") },
    { "slc_pzc_photon_1_true_trackid", new InhVecVecVar<int>("slc_pzc_photon_1_true_trackid") },
    { "slc_pzc_photon_1_true_pdg", new InhVecVecVar<int>("slc_pzc_photon_1_true_pdg") },
    { "slc_pzc_photon_1_comp", new InhVecVecVar<float>("slc_pzc_photon_1_comp") },
    { "slc_pzc_photon_1_pur", new InhVecVecVar<float>("slc_pzc_photon_1_pur") },
    { "slc_pzc_invariant_mass_fit", new InhVecVecVar<double>("slc_pzc_invariant_mass_fit") },
    { "slc_pzc_pizero_mom_fit", new InhVecVecVar<double>("slc_pzc_pizero_mom_fit") },
    { "slc_pzc_gamma0_energy_fit", new InhVecVecVar<double>("slc_pzc_gamma0_energy_fit") },
    { "slc_pzc_gamma1_energy_fit", new InhVecVecVar<double>("slc_pzc_gamma1_energy_fit") },
    { "slc_pzc_open_angle_fit", new InhVecVecVar<double>("slc_pzc_open_angle_fit") },
    { "slc_best_pzc_photon_0_id", new InhVecVar<size_t>("slc_best_pzc_photon_0_id") },
    { "slc_best_pzc_photon_1_id", new InhVecVar<size_t>("slc_best_pzc_photon_1_id") },
    { "slc_best_pzc_good_kinematics", new InhVecVar<bool>("slc_best_pzc_good_kinematics") },
    { "slc_best_pzc_invariant_mass", new InhVecVar<double>("slc_best_pzc_invariant_mass") },
    { "slc_best_pzc_pizero_mom", new InhVecVar<double>("slc_best_pzc_pizero_mom") },
    { "slc_best_pzc_cos_theta_pizero", new InhVecVar<double>("slc_best_pzc_cos_theta_pizero") },
    { "slc_best_pzc_cos_com", new InhVecVar<double>("slc_best_pzc_cos_com") },
    { "slc_best_pzc_decay_asymmetry", new InhVecVar<double>("slc_best_pzc_decay_asymmetry") },
    { "slc_best_pzc_good_kinematics_corr", new InhVecVar<bool>("slc_best_pzc_good_kinematics_corr") },
    { "slc_best_pzc_invariant_mass_corr", new InhVecVar<double>("slc_best_pzc_invariant_mass_corr") },
    { "slc_best_pzc_pizero_mom_corr", new InhVecVar<double>("slc_best_pzc_pizero_mom_corr") },
    { "slc_best_pzc_cos_theta_pizero_corr", new InhVecVar<double>("slc_best_pzc_cos_theta_pizero_corr") },
    { "slc_best_pzc_cos_com_corr", new InhVecVar<double>("slc_best_pzc_cos_com_corr") },
    { "slc_best_pzc_decay_asymmetry_corr", new InhVecVar<double>("slc_best_pzc_decay_asymmetry_corr") },
    { "slc_best_pzc_invariant_mass_fit", new InhVecVar<double>("slc_best_pzc_invariant_mass_fit") },
    { "slc_best_pzc_pizero_mom_fit", new InhVecVar<double>("slc_best_pzc_pizero_mom_fit") },
    { "slc_best_pzc_gamma0_energy_fit", new InhVecVar<double>("slc_best_pzc_gamma0_energy_fit") },
    { "slc_best_pzc_gamma1_energy_fit", new InhVecVar<double>("slc_best_pzc_gamma1_energy_fit") },
    { "slc_best_pzc_open_angle_fit", new InhVecVar<double>("slc_best_pzc_open_angle_fit") },
    { "slc_best_corr_pzc_photon_0_id", new InhVecVar<size_t>("slc_best_corr_pzc_photon_0_id") },
    { "slc_best_corr_pzc_photon_1_id", new InhVecVar<size_t>("slc_best_corr_pzc_photon_1_id") },
    { "slc_best_corr_pzc_good_kinematics", new InhVecVar<bool>("slc_best_corr_pzc_good_kinematics") },
    { "slc_best_corr_pzc_invariant_mass", new InhVecVar<double>("slc_best_corr_pzc_invariant_mass") },
    { "slc_best_corr_pzc_pizero_mom", new InhVecVar<double>("slc_best_corr_pzc_pizero_mom") },
    { "slc_best_corr_pzc_cos_theta_pizero", new InhVecVar<double>("slc_best_corr_pzc_cos_theta_pizero") },
    { "slc_best_corr_pzc_cos_com", new InhVecVar<double>("slc_best_corr_pzc_cos_com") },
    { "slc_best_corr_pzc_decay_asymmetry", new InhVecVar<double>("slc_best_corr_pzc_decay_asymmetry") },
    { "slc_best_corr_pzc_good_kinematics_corr", new InhVecVar<bool>("slc_best_corr_pzc_good_kinematics_corr") },
    { "slc_best_corr_pzc_invariant_mass_corr", new InhVecVar<double>("slc_best_corr_pzc_invariant_mass_corr") },
    { "slc_best_corr_pzc_pizero_mom_corr", new InhVecVar<double>("slc_best_corr_pzc_pizero_mom_corr") },
    { "slc_best_corr_pzc_cos_theta_pizero_corr", new InhVecVar<double>("slc_best_corr_pzc_cos_theta_pizero_corr") },
    { "slc_best_corr_pzc_cos_com_corr", new InhVecVar<double>("slc_best_corr_pzc_cos_com_corr") },
    { "slc_best_corr_pzc_decay_asymmetry_corr", new InhVecVar<double>("slc_best_corr_pzc_decay_asymmetry_corr") },
    { "slc_best_corr_pzc_invariant_mass_fit", new InhVecVar<double>("slc_best_corr_pzc_invariant_mass_fit") },
    { "slc_best_corr_pzc_pizero_mom_fit", new InhVecVar<double>("slc_best_corr_pzc_pizero_mom_fit") },
    { "slc_best_corr_pzc_gamma0_energy_fit", new InhVecVar<double>("slc_best_corr_pzc_gamma0_energy_fit") },
    { "slc_best_corr_pzc_gamma1_energy_fit", new InhVecVar<double>("slc_best_corr_pzc_gamma1_energy_fit") },
    { "slc_best_corr_pzc_open_angle_fit", new InhVecVar<double>("slc_best_corr_pzc_open_angle_fit") },
    { "slc_best_pzc_photon_0_true_trackid", new InhVecVar<int>("slc_best_pzc_photon_0_true_trackid") },
    { "slc_best_pzc_photon_0_true_pdg", new InhVecVar<int>("slc_best_pzc_photon_0_true_pdg") },
    { "slc_best_pzc_photon_0_comp", new InhVecVar<float>("slc_best_pzc_photon_0_comp") },
    { "slc_best_pzc_photon_0_pur", new InhVecVar<float>("slc_best_pzc_photon_0_pur") },
    { "slc_best_pzc_photon_1_true_trackid", new InhVecVar<int>("slc_best_pzc_photon_1_true_trackid") },
    { "slc_best_pzc_photon_1_true_pdg", new InhVecVar<int>("slc_best_pzc_photon_1_true_pdg") },
    { "slc_best_pzc_photon_1_comp", new InhVecVar<float>("slc_best_pzc_photon_1_comp") },
    { "slc_best_pzc_photon_1_pur", new InhVecVar<float>("slc_best_pzc_photon_1_pur") },
    { "slc_ssss_n_u_clusters", new InhVecVar<size_t>("slc_ssss_n_u_clusters") },
    { "slc_ssss_u_cluster_n_hits", new InhVecVecVar<size_t>("slc_ssss_u_cluster_n_hits") },
    { "slc_ssss_n_v_clusters", new InhVecVar<size_t>("slc_ssss_n_v_clusters") },
    { "slc_ssss_v_cluster_n_hits", new InhVecVecVar<size_t>("slc_ssss_v_cluster_n_hits") },
    { "slc_ssss_n_w_clusters", new InhVecVar<size_t>("slc_ssss_n_w_clusters") },
    { "slc_ssss_w_cluster_n_hits", new InhVecVecVar<size_t>("slc_ssss_w_cluster_n_hits") },
    { "slc_sel_incl", new InhVecVar<bool>("slc_sel_incl") },
    { "slc_sel_0p0pi", new InhVecVar<bool>("slc_sel_0p0pi") },
    { "slc_sel_1p0pi", new InhVecVar<bool>("slc_sel_1p0pi") },
    { "slc_sel_Np0pi", new InhVecVar<bool>("slc_sel_Np0pi") },
    { "slc_sel_Xp0pi", new InhVecVar<bool>("slc_sel_Xp0pi") },
    { "slc_sel_cc", new InhVecVar<bool>("slc_sel_cc") },
  };

  std::map<std::string, std::map<int, double>> genie_multisigma_universe_weights;
};

sbnd::DPreco::DPreco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCParticleModuleLabel          (p.get<art::InputTag>("MCParticleModuleLabel"))
  , fSliceModuleLabel               (p.get<art::InputTag>("SliceModuleLabel"))
  , fPFParticleModuleLabel          (p.get<art::InputTag>("PFParticleModuleLabel"))
  , fVertexModuleLabel              (p.get<art::InputTag>("VertexModuleLabel"))
  , fHitModuleLabel                 (p.get<art::InputTag>("HitModuleLabel"))
  , fTrackModuleLabel               (p.get<art::InputTag>("TrackModuleLabel"))
  , fShowerModuleLabel              (p.get<art::InputTag>("ShowerModuleLabel"))
  , fCRUMBSModuleLabel              (p.get<art::InputTag>("CRUMBSModuleLabel"))
  , fDazzleModuleLabel              (p.get<art::InputTag>("DazzleModuleLabel"))
  , fCaloModuleLabel                (p.get<art::InputTag>("CaloModuleLabel"))
  , fMCSModuleLabel                 (p.get<art::InputTag>("MCSModuleLabel"))
  , fChi2ModuleLabel                (p.get<art::InputTag>("Chi2ModuleLabel"))
  , fRangeModuleLabel               (p.get<art::InputTag>("RangeModuleLabel"))
  , fClosestApproachModuleLabel     (p.get<art::InputTag>("ClosestApproachModuleLabel"))
  , fStoppingChi2ModuleLabel        (p.get<art::InputTag>("StoppingChi2ModuleLabel"))
  , fRazzleModuleLabel              (p.get<art::InputTag>("RazzleModuleLabel"))
  , fCosmicDistModuleLabel          (p.get<art::InputTag>("CosmicDistModuleLabel"))
  , fShowerTrackFitModuleLabel      (p.get<art::InputTag>("ShowerTrackFitModuleLabel"))
  , fShowerDensityFitModuleLabel    (p.get<art::InputTag>("ShowerDensityFitModuleLabel"))
  , fPOTModuleLabel                 (p.get<art::InputTag>("POTModuleLabel"))
  , fOpT0ModuleLabel                (p.get<art::InputTag>("OpT0ModuleLabel"))
  , fRazzledModuleLabel             (p.get<art::InputTag>("RazzledModuleLabel"))
  , fSpacePointModuleLabel          (p.get<art::InputTag>("SpacePointModuleLabel"))
  , fEventWeightModuleLabels        (p.get<std::vector<art::InputTag>>("EventWeightModuleLabels"))
  , fDebug                          (p.get<bool>("Debug", false))
  , fBeamOff                        (p.get<bool>("BeamOff", false))
  {
    for(auto const& name : flux_weight_names)
      {
        nuVars["nu_weight_" + name ] = new InhVecVecVar<float>("nu_weight_" + name);
        slcVars["slc_true_weight_" + name ] = new InhVecVecVar<float>("slc_true_weight_" + name);
      }

    nuVars["nu_weight_all_flux"] = new InhVecVecVar<float>("nu_weight_all_flux");
    slcVars["slc_true_weight_all_flux"] = new InhVecVecVar<float>("slc_true_weight_all_flux");

    for(auto const& name : genie_weight_names)
      {
        nuVars["nu_weight_" + name ] = new InhVecVecVar<float>("nu_weight_" + name);
        slcVars["slc_true_weight_" + name ] = new InhVecVecVar<float>("slc_true_weight_" + name);

        if(name.find("multisigma") != std::string::npos)
          {
            nuVars["nu_weight_" + name + "_multisim"] = new InhVecVecVar<float>("nu_weight_" + name + "_multisim");
            slcVars["slc_true_weight_" + name + "_multisim" ] = new InhVecVecVar<float>("slc_true_weight_" + name + "_multisim");
          }
      }

    for(auto const& name : geant4_weight_names)
      {
        nuVars["nu_weight_" + name ] = new InhVecVecVar<float>("nu_weight_" + name);
        slcVars["slc_true_weight_" + name ] = new InhVecVecVar<float>("slc_true_weight_" + name);
      }

    nuVars["nu_weight_all_genie_multisim"] = new InhVecVecVar<float>("nu_weight_all_genie_multisim");
    slcVars["slc_true_weight_all_genie_multisim"] = new InhVecVecVar<float>("slc_true_weight_all_genie_multisim");

    cet::search_path sp("FW_SEARCH_PATH");
    std::string showerEnergyCorrectionFileFullPath;
    //if (!sp.find_file(fShowerEnergyCorrectionFileName, showerEnergyCorrectionFileFullPath))
    //throw cet::exception("DPreco") << "Could not find shower energy correction file: \n"
    //                                               << fShowerEnergyCorrectionFileName;

    //TFile* file = TFile::Open(showerEnergyCorrectionFileFullPath.c_str());
    //fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

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

    

    for(auto const& name : genie_weight_names)
      {
        if(name.find("multisigma") != std::string::npos)
          {
            genie_multisigma_universe_weights[name] = std::map<int, double>();
	    
	    
          }
      }
  }

void sbnd::DPreco::SetupBranches(VecVarMap &map)
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

void sbnd::DPreco::beginSubRun(const art::SubRun &sr)
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

void sbnd::DPreco::endSubRun(const art::SubRun &sr)
{
  fSubRunTree->Fill();
}

void sbnd::DPreco::analyze(const art::Event &e)
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
  AnalyseSlices(e, sliceHandle, pfpHandle, trackHandle, showerHandle, hitHandle);
  

  // Fill the Tree
  fEventTree->Fill();
}

void sbnd::DPreco::ClearMaps()
{
  fHitsMap.clear();
  fNuHitsMap.clear();
  fNuHitsMapSP.clear();
  fPFPMap.clear();
  fRecoPFPMap.clear();
}

void sbnd::DPreco::SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                                       const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  art::FindOneP<recob::SpacePoint> hitsToSPs(hitHandle, e, fSpacePointModuleLabel);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  for(auto const& hit : hitVec)
    {
      const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      fHitsMap[trackID]++;
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      fNuHitsMap[mct]++;

      const art::Ptr<recob::SpacePoint> spacepoint = hitsToSPs.at(hit.key());
      if(spacepoint.isNonnull())
        fNuHitsMapSP[mct]++;
    }

  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  art::fill_ptr_vector(pfpVec, pfpHandle);

  for(auto const& pfp : pfpVec)
    fPFPMap[pfp->Self()] = pfp;
}



void sbnd::DPreco::AnalyseMCTruth(const art::Event &e, VecVarMap &vars, const art::Ptr<simb::MCTruth> &mct, const int counter, const std::string prefix)
{
  if(mct->Origin() == 2)
    {
      FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kCosmic);
      FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kCosmic);
      FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kCosmic);
      FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kCosmic);
      FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kCosmic);
      FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kCosmic);
      return;
    }
  else if(mct->Origin() != 1)
    {
      FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kUnknownEv);
      FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kUnknownEv);
      FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kUnknownEv);
      FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kUnknownEv);
      FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kUnknownEv);
      FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kUnknownEv);
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

  for(auto const& weightModuleLabel : fEventWeightModuleLabels)
    {
      art::FindManyP<sbn::evwgh::EventWeightMap> MCTruthToWeights( { mct }, e, weightModuleLabel);
      const std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> ewms = MCTruthToWeights.at(0);

      int n_univs = 0;
      if(weightModuleLabel == "fluxweight")
        n_univs = n_fluxweight_univs;
      else if(weightModuleLabel == "systtools")
        n_univs = n_genieweight_univs;
      else if(weightModuleLabel == "geant4weight")
        n_univs = n_geant4weight_univs;

      std::vector<float> all(n_univs, 1.);

      for(auto const& ewm : ewms)
        {
          for(auto const& [ name, weights ] : *ewm)
            {
              if((weightModuleLabel == "fluxweight") || (weightModuleLabel == "geant4weight") || (weightModuleLabel == "systtools" && name.find("multisim") != std::string::npos))
                {
                  FillElement(vars[prefix + "_weight_" + name], counter, weights);

                  for(int univ = 0; univ < n_univs; ++univ)
                    all[univ] *= weights[univ];
                }
              else if(weightModuleLabel == "systtools" && name.find("multisigma") != std::string::npos)
                {
                  FillElement(vars[prefix + "_weight_" + name], counter, weights);

                  std::vector<float> thrown_weights(n_univs, 1.);

                  if(weights.size() == 6)
                    {
                      double multisigma_sigmas[6] = { -1, 1, -2, 2, -3, 3 };
                      double multisigma_vals[6]   = { weights[0], weights[1], weights[2], weights[3], weights[4], weights[5] };

                      TSpline3 *spline = new TSpline3(Form("%s_spline", name.c_str()), multisigma_sigmas, multisigma_vals, 6);

                      for(int univ = 0; univ < n_univs; ++univ)
                        {
                          thrown_weights[univ] = spline->Eval(genie_multisigma_universe_weights[name][univ]);
                          all[univ] *= thrown_weights[univ];
                        }
                    }
                  else if(weights.size() == 1)
                    {
                      for(int univ = 0; univ < n_univs; ++univ)
                        {
                          thrown_weights[univ] = 1 + (weights[0] - 1) * 2 * genie_multisigma_universe_weights[name][univ];
                          all[univ] *= thrown_weights[univ];
                        }
                    }
                  else
                    std::cout << "Whoaaaaaaaaaa, multisigma of size " << weights.size() << std::endl;

                  FillElement(vars[prefix + "_weight_" + name + "_multisim"], counter, thrown_weights);
                }
            }
        }

      if(weightModuleLabel == "fluxweight")
        FillElement(vars[prefix + "_weight_all_flux" ], counter, all);
      else if(weightModuleLabel == "systtools")
        FillElement(vars[prefix + "_weight_all_genie_multisim" ], counter, all);
    }

  int protons = 0, neutrons = 0, charged_pions = 0, neutral_pions = 0, dalitz_neutral_pions = 0, photons = 0, other = 0;
  float trueEnDep = 0.;

  for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
        {
          switch(abs(mcp->PdgCode()))
            {
            case 2212:
              if(mcp->P() > .4)
                ++protons;
              break;
            case 2112:
              ++neutrons;
              break;
            case 211:
              if(mcp->P() > .15)
                ++charged_pions;
              break;
            case 111:
              if(mcp->NumberDaughters() == 2)
                ++neutral_pions;
              else
                ++dalitz_neutral_pions;
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

  const bool pizero = neutral_pions == 1;

  if(nc && fv && pizero)
    {
      FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kSignalNCPiZero);
      FillElement(vars[prefix + "_signal"], counter, true);

      if(charged_pions == 0)
        {
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kSignalNCPiZero);

          if(protons == 0)
            FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kSignalNCPiZero);
          else
            FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kOtherNCPiZero);

          if(protons == 1)
            FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kSignalNCPiZero);
          else
            FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kOtherNCPiZero);

          if(protons > 0)
            FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kSignalNCPiZero);
          else
            FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kOtherNCPiZero);
        }
      else
        {
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kOtherNCPiZero);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kOtherNCPiZero);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kOtherNCPiZero);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kOtherNCPiZero);
        }

      FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kNC);
    }
  else
    {
      FillElement(vars[prefix + "_signal"], counter, false);

      if(nc && fv)
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kOtherNC);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kOtherNC);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kOtherNC);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kOtherNC);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kOtherNC);
          FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kNC);
        }
      else if(abs(nu.PdgCode()) == 14 && !nc && fv)
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kCCNuMu);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kCCNuMu);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kCCNuMu);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kCCNuMu);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kCCNuMu);
          if(pizero)
            FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kSignalCCPiZero);
          else
            FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kOtherCCNuMu);

        }
      else if(abs(nu.PdgCode()) == 12 && !nc && fv)
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kCCNuE);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kCCNuE);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kCCNuE);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kCCNuE);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kCCNuE);
          FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kCCNuE);
        }
      else if(!fv && av)
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kNonFV);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kNonFV);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kNonFV);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kNonFV);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kNonFV);
          FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kNonFV);
        }
      else if(!av)
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kDirt);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kDirt);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kDirt);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kDirt);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kDirt);
          FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kDirt);
        }
      else
        {
          FillElement(vars[prefix + "_event_type_incl"], counter, (int) NC::kUnknownEv);
          FillElement(vars[prefix + "_event_type_0p0pi"], counter, (int) NC::kUnknownEv);
          FillElement(vars[prefix + "_event_type_1p0pi"], counter, (int) NC::kUnknownEv);
          FillElement(vars[prefix + "_event_type_Np0pi"], counter, (int) NC::kUnknownEv);
          FillElement(vars[prefix + "_event_type_Xp0pi"], counter, (int) NC::kUnknownEv);
          FillElement(vars[prefix + "_event_type_cc"], counter, (int) CC::kUnknownEv);
        }
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
  FillElement(vars[prefix + "_n_dalitz_neutral_pions"], counter, dalitz_neutral_pions);
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
          if(abs(mcp->PdgCode()) == 111 && mcp->NumberDaughters() == 2)
            {
              FillElement(vars[prefix + "_pz_invariant_mass"], counter, pzCounter, mcp->Mass());
              FillElement(vars[prefix + "_pz_pizero_mom"], counter, pzCounter, mcp->P());
              FillElement(vars[prefix + "_pz_cos_theta_pizero"], counter, pzCounter, mcp->Pz() / mcp->P());

              const simb::MCParticle* gamma0 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
              const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

              const bool two_gamma_decay = (gamma0->PdgCode() == 22 && gamma1->PdgCode() == 22);
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
              FillElement(vars[prefix + "_pz_gamma0_trackid"], counter, pzCounter, gamma0->TrackId());
              FillElement(vars[prefix + "_pz_gamma0_n_hits"], counter, pzCounter, fHitsMap[gamma0->TrackId()]);
              FillElement(vars[prefix + "_pz_gamma0_energy"], counter, pzCounter, en0);
              FillElement(vars[prefix + "_pz_gamma0_dir_x"], counter, pzCounter, gamma0->Px() / gamma0->P());
              FillElement(vars[prefix + "_pz_gamma0_dir_y"], counter, pzCounter, gamma0->Py() / gamma0->P());
              FillElement(vars[prefix + "_pz_gamma0_dir_z"], counter, pzCounter, gamma0->Pz() / gamma0->P());
              FillElement(vars[prefix + "_pz_gamma1_trackid"], counter, pzCounter, gamma1->TrackId());
              FillElement(vars[prefix + "_pz_gamma1_n_hits"], counter, pzCounter, fHitsMap[gamma1->TrackId()]);
              FillElement(vars[prefix + "_pz_gamma1_energy"], counter, pzCounter, en1);
              FillElement(vars[prefix + "_pz_gamma1_dir_x"], counter, pzCounter, gamma1->Px() / gamma1->P());
              FillElement(vars[prefix + "_pz_gamma1_dir_y"], counter, pzCounter, gamma1->Py() / gamma1->P());
              FillElement(vars[prefix + "_pz_gamma1_dir_z"], counter, pzCounter, gamma1->Pz() / gamma1->P());
              FillElement(vars[prefix + "_pz_open_angle"], counter, pzCounter, TMath::RadToDeg() * gamma0->Momentum().Vect().Angle(gamma1->Momentum().Vect()));

              if(prefix.find("nu") != std::string::npos)
                {
                  AnalysePhotonReco("gamma0", gamma0->TrackId(), counter, pzCounter);
                  AnalysePhotonReco("gamma1", gamma1->TrackId(), counter, pzCounter);
                }

              ++pzCounter;
            }
        }
    }
}



void sbnd::DPreco::AnalysePhotonReco(const std::string name, const int trackid, const int nuCounter, const int pzCounter)
{
  std::vector<std::vector<float>> slc_pfp_comp;
  std::vector<std::vector<int>> slc_pfp_true_trackid;

  GetVar(slcVars["slc_pfp_true_trackid"], slc_pfp_true_trackid);
  GetVar(slcVars["slc_pfp_comp"], slc_pfp_comp);

  float bestComp = std::numeric_limits<float>::lowest();
  std::pair<int, int> bestPFP = { -1, -1 };

  for(size_t i = 0; i < slc_pfp_comp.size(); ++i)
    {
      for(size_t j = 0; j < slc_pfp_comp.at(i).size(); ++j)
        {
          if(slc_pfp_true_trackid.at(i).at(j) == trackid && slc_pfp_comp.at(i).at(j) > bestComp)
            {
              bestComp = slc_pfp_comp.at(i).at(j);
              bestPFP = { i, j };
            }
        }
    }

  if(bestPFP.first ==  -1 || bestPFP.second == -1 )
    return;

  for(auto& [varName, var] : nuVars)
    {
      if(varName.find("nu_pz_" + name + "_best_pfp") != std::string::npos)
        TransferElement(var, slcVars, "nu_pz_" + name + "_best_pfp", "slc_pfp", nuCounter, pzCounter, bestPFP.first, bestPFP.second);
    }
}

void sbnd::DPreco::AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                           const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                                           const art::Handle<std::vector<recob::Track>> &trackHandle,
                                           const art::Handle<std::vector<recob::Shower>> &showerHandle,
                                           const art::Handle<std::vector<recob::Hit>> &hitHandle)
{
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  _n_slc = sliceVec.size();
  ResizeVectors(slcVars, _n_slc);

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>      pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  
  art::FindManyP<sbn::OpT0Finder>   slicesToOpT0(sliceHandle, e, fOpT0ModuleLabel);
  art::FindManyP<recob::Hit>        slicesToHits(sliceHandle, e, fSliceModuleLabel);

  for (auto&& [slcCounter, slc] : enumerate(sliceVec))
    {
      FillElement(slcVars["slc_key"], slcCounter, slc.key());

      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slc.key());
      FillElement(slcVars["slc_n_pfps"], slcCounter, pfps.size());

      const std::vector<art::Ptr<recob::Hit>> hits = slicesToHits.at(slc.key());
      FillElement(slcVars["slc_n_hits"], slcCounter, hits.size());

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
      FillElement(slcVars["slc_n_primary_children"], slcCounter, prim->NumDaughters());

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
          FillElement(slcVars["slc_opt0_fracPE"], slcCounter, (opT0Vec[0]->hypoPE - opT0Vec[0]->measPE) / opT0Vec[0]->measPE);
        }

      ResizeSubVectors(slcVars, "slc_pfp", slcCounter, pfps.size());

      AnalysePFPs(e, prim, pfps, hits, vtx, slcCounter, pfpHandle, trackHandle, showerHandle);

      

      SelectSlice(slcCounter);

      AnalyseSliceTruth(e, slc, slcCounter, sliceHandle, hitHandle);
    }
}

void sbnd::DPreco::AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const std::vector<art::Ptr<recob::PFParticle>> &pfps,
                                         const std::vector<art::Ptr<recob::Hit>> &sliceHits, const art::Ptr<recob::Vertex> &vtx, const int slcCounter,
                                         const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                                         const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  art::FindOneP<recob::Track>                      pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpToMeta(pfpHandle, e, fPFParticleModuleLabel);
  art::FindManyP<recob::Hit>                       showersToHits(showerHandle, e, fShowerModuleLabel);
  art::FindOneP<sbn::MVAPID>                       pfpToRazzled(pfpHandle, e, fRazzledModuleLabel);
  art::FindManyP<recob::SpacePoint>                pfpToSpacePoints(pfpHandle, e, fSpacePointModuleLabel);

  int ntrks = 0, nshws = 0, ndazzlemuons = 0, ndazzlepions = 0, ndazzlepionsthresh = 0, ndazzleprotons = 0,
    ndazzleprotonsthresh = 0, ndazzleother = 0, nrazzleelectrons = 0, nrazzlephotons = 0, nrazzleother = 0,
    nrazzledelectrons = 0, nrazzledmuons = 0, nrazzledphotons = 0, nrazzledpions = 0, nrazzledpionsthresh = 0,
    nrazzledprotons = 0, nrazzledprotonsthresh = 0;

  int nprimtrks = 0, nprimshws = 0, nprimdazzlemuons = 0, nprimdazzlepions = 0, nprimdazzlepionsthresh = 0, nprimdazzleprotons = 0,
    nprimdazzleprotonsthresh = 0, nprimdazzleother = 0, nprimrazzleelectrons = 0, nprimrazzlephotons = 0, nprimrazzleother = 0,
    nprimrazzledelectrons = 0, nprimrazzledmuons = 0, nprimrazzledphotons = 0, nprimrazzledpions = 0, nprimrazzledpionsthresh = 0,
    nprimrazzledprotons = 0, nprimrazzledprotonsthresh = 0;

  std::vector<art::Ptr<recob::Hit>> used_hits;

  for(auto&& [pfpCounter, pfp] : enumerate(pfps))
    {
      FillElement(slcVars["slc_pfp_id"], slcCounter, pfpCounter, pfp->Self());
      FillElement(slcVars["slc_pfp_pdg"], slcCounter, pfpCounter, pfp->PdgCode());
      FillElement(slcVars["slc_pfp_n_children"], slcCounter, pfpCounter, pfp->NumDaughters());

      const bool primary_child = prim->Self() == pfp->Parent();
      FillElement(slcVars["slc_pfp_primary"], slcCounter, pfpCounter, pfp->IsPrimary());
      FillElement(slcVars["slc_pfp_primary_child"], slcCounter, pfpCounter, primary_child);

      if(abs(pfp->PdgCode()) == 11)
        {
          ++nshws;
          if(primary_child)
            ++nprimshws;
        }
      else if(abs(pfp->PdgCode()) == 13)
        {
          ++ntrks;
          if(primary_child)
            ++nprimtrks;
        }
      else
        {
          continue;
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

      const art::Ptr<sbn::MVAPID> razzled = pfpToRazzled.at(pfp.key());
      if(razzled.isNonnull())
        ExtractRazzled(razzled, slcCounter, pfpCounter);

      const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
      const std::vector<art::Ptr<recob::Hit>> hits = showersToHits.at(shower.key());

      FillElement(slcVars["slc_pfp_n_hits"], slcCounter, pfpCounter, hits.size());
      used_hits.insert(used_hits.end(), hits.begin(), hits.end());

      const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, true);
      FillElement(slcVars["slc_pfp_true_trackid"], slcCounter, pfpCounter, trackID);
      FillElement(slcVars["slc_pfp_comp"], slcCounter, pfpCounter, Completeness(e, hits, trackID));
      FillElement(slcVars["slc_pfp_pur"], slcCounter, pfpCounter, Purity(e, hits, trackID));

      const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpToSpacePoints.at(pfp.key());
      FillElement(slcVars["slc_pfp_n_sps"], slcCounter, pfpCounter, spacepoints.size());

      
      
      const simb::MCParticle* mcp = particleInv->TrackIdToParticle_P(trackID);
      if(mcp != NULL){
	  FillElement(slcVars["slc_pfp_true_pdg"], slcCounter, pfpCounter, mcp->PdgCode());
	  FillElement(slcVars["slc_pfp_true_energy"], slcCounter, pfpCounter, mcp->E());
	  FillElement(slcVars["slc_pfp_true_ke"], slcCounter, pfpCounter, mcp->E() - mcp->Mass());
	  FillElement(slcVars["slc_pfp_true_p_x"], slcCounter, pfpCounter, mcp->Px());
	  FillElement(slcVars["slc_pfp_true_p_y"], slcCounter, pfpCounter, mcp->Py());
	  FillElement(slcVars["slc_pfp_true_p_z"], slcCounter, pfpCounter, mcp->Pz());
      }
      
      
      if(track.isNonnull())
        AnalyseTrack(e, track, slcCounter, pfpCounter, trackHandle);

      double pfpenergy = def_double;

      if(pfp->PdgCode() == 13)
        {
          int dazzlepdg;
          AccessElement(slcVars["slc_pfp_track_dazzle_pdg"], slcCounter, pfpCounter, dazzlepdg);

          float trkenergy;
          AccessElement(slcVars["slc_pfp_track_ke"], slcCounter, pfpCounter, trkenergy);
          pfpenergy = trkenergy;

          if(dazzlepdg == 13)
            ++ndazzlemuons;
          if(dazzlepdg == 211)
            ++ndazzlepions;
          if(dazzlepdg == 211 && pfpenergy > 65.3)
            ++ndazzlepionsthresh;
          if(dazzlepdg == 2212)
            ++ndazzleprotons;
          if(dazzlepdg == 2212 && pfpenergy > 81.7)
            ++ndazzleprotonsthresh;
          if(dazzlepdg == 0)
            ++ndazzleother;

          if(primary_child)
            {
              if(dazzlepdg == 13)
                ++nprimdazzlemuons;
              if(dazzlepdg == 211)
                ++nprimdazzlepions;
              if(dazzlepdg == 211 && pfpenergy > 32.1)
                ++nprimdazzlepionsthresh;
              if(dazzlepdg == 2212)
                ++nprimdazzleprotons;
              if(dazzlepdg == 2212 && pfpenergy > 32.7)
                ++nprimdazzleprotonsthresh;
              if(dazzlepdg == 0)
                ++nprimdazzleother;
            }
        }

      if(shower.isNonnull())
        AnalyseShower(e, shower, slcCounter, pfpCounter, showerHandle, vtx, hits);

      if(pfp->PdgCode() == 11)
        {
          int razzlepdg;
          AccessElement(slcVars["slc_pfp_shower_razzle_pdg"], slcCounter, pfpCounter, razzlepdg);

          AccessElement(slcVars["slc_pfp_shower_energy"], slcCounter, pfpCounter, pfpenergy);

          if(razzlepdg == 11)
            ++nrazzleelectrons;
          if(razzlepdg == 22)
            ++nrazzlephotons;
          if(razzlepdg == 0)
            ++nrazzleother;

          if(primary_child)
            {
              if(razzlepdg == 11)
                ++nprimrazzleelectrons;
              if(razzlepdg == 22)
                ++nprimrazzlephotons;
              if(razzlepdg == 0)
                ++nprimrazzleother;
            }
        }

      int razzledpdg;
      AccessElement(slcVars["slc_pfp_razzled_pdg"], slcCounter, pfpCounter, razzledpdg);

      if(razzledpdg == 11)
        ++nrazzledelectrons;
      if(razzledpdg == 13)
        ++nrazzledmuons;
      if(razzledpdg == 22)
        ++nrazzledphotons;
      if(razzledpdg == 211)
        ++nrazzledpions;
      if(razzledpdg == 211 && pfpenergy > 32.1)
        ++nrazzledpionsthresh;
      if(razzledpdg == 2212)
        ++nrazzledprotons;
      if(razzledpdg == 2212 && pfpenergy > 32.7)
        ++nrazzledprotonsthresh;

      if(primary_child)
        {
          if(razzledpdg == 11)
            ++nprimrazzledelectrons;
          if(razzledpdg == 13)
            ++nprimrazzledmuons;
          if(razzledpdg == 22)
            ++nprimrazzledphotons;
          if(razzledpdg == 211)
            ++nprimrazzledpions;
          if(razzledpdg == 211 && pfpenergy > 32.1)
            ++nprimrazzledpionsthresh;
          if(razzledpdg == 2212)
            ++nprimrazzledprotons;
          if(razzledpdg == 2212 && pfpenergy > 32.7)
            ++nprimrazzledprotonsthresh;
        }
    }

  FillElement(slcVars["slc_n_trks"], slcCounter, ntrks);
  FillElement(slcVars["slc_n_shws"], slcCounter, nshws);
  FillElement(slcVars["slc_n_dazzle_muons"], slcCounter, ndazzlemuons);
  FillElement(slcVars["slc_n_dazzle_pions"], slcCounter, ndazzlepions);
  FillElement(slcVars["slc_n_dazzle_pions_thresh"], slcCounter, ndazzlepionsthresh);
  FillElement(slcVars["slc_n_dazzle_protons"], slcCounter, ndazzleprotons);
  FillElement(slcVars["slc_n_dazzle_protons_thresh"], slcCounter, ndazzleprotonsthresh);
  FillElement(slcVars["slc_n_dazzle_other"], slcCounter, ndazzleother);
  FillElement(slcVars["slc_n_razzle_electrons"], slcCounter, nrazzleelectrons);
  FillElement(slcVars["slc_n_razzle_photons"], slcCounter, nrazzlephotons);
  FillElement(slcVars["slc_n_razzle_other"], slcCounter, nrazzleother);
  FillElement(slcVars["slc_n_razzled_electrons"], slcCounter, nrazzledelectrons);
  FillElement(slcVars["slc_n_razzled_muons"], slcCounter, nrazzledmuons);
  FillElement(slcVars["slc_n_razzled_photons"], slcCounter, nrazzledphotons);
  FillElement(slcVars["slc_n_razzled_pions"], slcCounter, nrazzledpions);
  FillElement(slcVars["slc_n_razzled_pions_thresh"], slcCounter, nrazzledpionsthresh);
  FillElement(slcVars["slc_n_razzled_protons"], slcCounter, nrazzledprotons);
  FillElement(slcVars["slc_n_razzled_protons_thresh"], slcCounter, nrazzledprotonsthresh);

  FillElement(slcVars["slc_n_primary_trks"], slcCounter, nprimtrks);
  FillElement(slcVars["slc_n_primary_shws"], slcCounter, nprimshws);
  FillElement(slcVars["slc_n_primary_dazzle_muons"], slcCounter, nprimdazzlemuons);
  FillElement(slcVars["slc_n_primary_dazzle_pions"], slcCounter, nprimdazzlepions);
  FillElement(slcVars["slc_n_primary_dazzle_pions_thresh"], slcCounter, nprimdazzlepionsthresh);
  FillElement(slcVars["slc_n_primary_dazzle_protons"], slcCounter, nprimdazzleprotons);
  FillElement(slcVars["slc_n_primary_dazzle_protons_thresh"], slcCounter, nprimdazzleprotonsthresh);
  FillElement(slcVars["slc_n_primary_dazzle_other"], slcCounter, nprimdazzleother);
  FillElement(slcVars["slc_n_primary_razzle_electrons"], slcCounter, nprimrazzleelectrons);
  FillElement(slcVars["slc_n_primary_razzle_photons"], slcCounter, nprimrazzlephotons);
  FillElement(slcVars["slc_n_primary_razzle_other"], slcCounter, nprimrazzleother);
  FillElement(slcVars["slc_n_primary_razzled_electrons"], slcCounter, nprimrazzledelectrons);
  FillElement(slcVars["slc_n_primary_razzled_muons"], slcCounter, nprimrazzledmuons);
  FillElement(slcVars["slc_n_primary_razzled_photons"], slcCounter, nprimrazzledphotons);
  FillElement(slcVars["slc_n_primary_razzled_pions"], slcCounter, nprimrazzledpions);
  FillElement(slcVars["slc_n_primary_razzled_pions_thresh"], slcCounter, nprimrazzledpionsthresh);
  FillElement(slcVars["slc_n_primary_razzled_protons"], slcCounter, nprimrazzledprotons);
  FillElement(slcVars["slc_n_primary_razzled_protons_thresh"], slcCounter, nprimrazzledprotonsthresh);

  FillElement(slcVars["slc_n_used_hits"], slcCounter, used_hits.size());

  if((prim->PdgCode() == 12 || prim->PdgCode() == 14) && nprimrazzledmuons == 0 && (nprimrazzledphotons==1 || nprimrazzledelectrons==1))
    {
      std::vector<art::Ptr<recob::Hit>> unused_hits;

      for(auto const& hit : sliceHits)
        {
          if(std::find(used_hits.begin(), used_hits.end(), hit) == used_hits.end())
            unused_hits.push_back(hit);
        }

      

      
    }
}

bool sbnd::DPreco::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
}

bool sbnd::DPreco::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}


void sbnd::DPreco::AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, const int slcCounter, const int pfpCounter,
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

  geo::Point_t end = track->End();
  const bool contained = VolumeCheck(start, 10, 0, 10, 10) && VolumeCheck(end, 10, 0, 10, 10);
  FillElement(slcVars["slc_pfp_track_contained"], slcCounter, pfpCounter, contained);

  FillElement(slcVars["slc_pfp_track_length"], slcCounter, pfpCounter, track->Length());

  
  

  const std::vector<art::Ptr<anab::Calorimetry>> calos = tracksToCalos.at(track.key());
  
  

  
  

  const art::Ptr<sbn::RangeP> rangeP = tracksToRangePs.at(track.key());
  if(rangeP.isNonnull())
    FillElement(slcVars["slc_pfp_track_range_p"], slcCounter, pfpCounter, rangeP->range_p);

  const art::Ptr<sbn::ScatterClosestApproach> closestApproach = tracksToClosestApproaches.at(track.key());
  if(closestApproach.isNonnull())
    FillElement(slcVars["slc_pfp_track_closest_approach_mean_dca"], slcCounter, pfpCounter, closestApproach->mean);

  
}



void sbnd::DPreco::AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
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

  geo::Point_t end = start + shower->Length() * dir;
  const bool contained = VolumeCheck(start, 10, 0, 10, 10) && VolumeCheck(end, 10, 0, 10, 10);
  FillElement(slcVars["slc_pfp_shower_contained"], slcCounter, pfpCounter, contained);

  ExtractCalo(shower, slcCounter, pfpCounter, hits);


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

void sbnd::DPreco::ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,const std::vector<art::Ptr<recob::Hit>> &hits)
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

  int bestPlane = -1;

  if(showerPlaneHits[2] >= showerPlaneHits[1] && showerPlaneHits[2] >= showerPlaneHits[0])
    bestPlane = 2;
  else if(showerPlaneHits[0] >= showerPlaneHits[1])
    bestPlane = 0;
  else
    bestPlane = 1;

  FillElement(slcVars["slc_pfp_shower_energy"], slcCounter, pfpCounter, shower->Energy()[bestPlane]);
  FillElement(slcVars["slc_pfp_shower_energy_coll"], slcCounter, pfpCounter, shower->Energy()[2]);
  FillElement(slcVars["slc_pfp_shower_energy_corr"], slcCounter, pfpCounter, CorrectEnergy(shower->Energy()[bestPlane]));
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

void sbnd::DPreco::ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = razzled->mvaScoreMap;

  FillElement(slcVars["slc_pfp_razzled_electron_score"], slcCounter, pfpCounter, map.at(11));
  FillElement(slcVars["slc_pfp_razzled_muon_score"], slcCounter, pfpCounter, map.at(13));
  FillElement(slcVars["slc_pfp_razzled_photon_score"], slcCounter, pfpCounter, map.at(22));
  FillElement(slcVars["slc_pfp_razzled_pion_score"], slcCounter, pfpCounter, map.at(211));
  FillElement(slcVars["slc_pfp_razzled_proton_score"], slcCounter, pfpCounter, map.at(2212));
  FillElement(slcVars["slc_pfp_razzled_pdg"], slcCounter, pfpCounter, razzled->BestPDG());
}

void sbnd::DPreco::AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, const int slcCounter,
                                               const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                               const art::Handle<std::vector<recob::Hit>> &hitHandle)
{
  art::FindManyP<recob::Hit> slicesToHits(sliceHandle, e, fSliceModuleLabel);
  art::FindOneP<recob::SpacePoint> hitsToSPs(hitHandle, e, fSpacePointModuleLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const std::vector<art::Ptr<recob::Hit>> sliceHits = slicesToHits.at(slc.key());

  int nSPHits = 0;
  std::map<int, int> objectHitMap, objectHitMapSP;
  for(auto const &hit : sliceHits)
    {
      const art::Ptr<recob::SpacePoint> spacepoint = hitsToSPs.at(hit.key());
      if(spacepoint.isNonnull())
        {
          objectHitMapSP[TruthMatchUtils::TrueParticleID(clockData, hit, true)]++;
          ++nSPHits;
        }

      objectHitMap[TruthMatchUtils::TrueParticleID(clockData, hit, true)]++;
    }

  std::map<const art::Ptr<simb::MCTruth>, int> mcTruthHitMap, mcTruthHitMapSP;
  for(auto const& [trackID, nhits] : objectHitMap)
    {
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      mcTruthHitMap[mct] += nhits;
    }
  for(auto const& [trackID, nhits] : objectHitMapSP)
    {
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      mcTruthHitMapSP[mct] += nhits;
    }

  int maxHits = def_int, maxHitsSP = def_int;
  art::Ptr<simb::MCTruth> bestMCT = art::Ptr<simb::MCTruth>(), bestMCTSP = art::Ptr<simb::MCTruth>();

  for(auto const& [mct, nhits] : mcTruthHitMap)
    {
      if(nhits > maxHits)
        {
          maxHits = nhits;
          bestMCT = mct;
        }
    }

  for(auto const& [mct, nhits] : mcTruthHitMapSP)
    {
      if(nhits > maxHitsSP)
        {
          maxHitsSP = nhits;
          bestMCTSP = mct;
        }
    }

  const float comp = fNuHitsMap[bestMCT] == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(fNuHitsMap[bestMCT]);
  const float pur  = sliceHits.size() == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(sliceHits.size());

  const float compSP = fNuHitsMapSP[bestMCTSP] == 0 ? def_float : mcTruthHitMapSP[bestMCTSP] / static_cast<float>(fNuHitsMapSP[bestMCTSP]);
  const float purSP  = nSPHits == 0 ? def_float : mcTruthHitMapSP[bestMCTSP] / static_cast<float>(nSPHits);

  FillElement(slcVars["slc_comp"], slcCounter, comp);
  FillElement(slcVars["slc_pur"], slcCounter, pur);
  FillElement(slcVars["slc_comp_sp_only"], slcCounter, compSP);
  FillElement(slcVars["slc_pur_sp_only"], slcCounter, purSP);

  if(fBeamOff)
    {
      FillElement(slcVars["slc_true_event_type_incl"], slcCounter, (int) NC::kCosmic);
      FillElement(slcVars["slc_true_event_type_0p0pi"], slcCounter, (int) NC::kCosmic);
      FillElement(slcVars["slc_true_event_type_1p0pi"], slcCounter, (int) NC::kCosmic);
      FillElement(slcVars["slc_true_event_type_Np0pi"], slcCounter, (int) NC::kCosmic);
      FillElement(slcVars["slc_true_event_type_Xp0pi"], slcCounter, (int) NC::kCosmic);
      FillElement(slcVars["slc_true_event_type_cc"], slcCounter, (int) CC::kCosmic);
    }
  else if(bestMCT.isNonnull())
    AnalyseMCTruth(e, slcVars, bestMCT, slcCounter, "slc_true");
  else
    {
      FillElement(slcVars["slc_true_event_type_incl"], slcCounter, (int) NC::kFailedTruthMatch);
      FillElement(slcVars["slc_true_event_type_0p0pi"], slcCounter, (int) NC::kFailedTruthMatch);
      FillElement(slcVars["slc_true_event_type_1p0pi"], slcCounter, (int) NC::kFailedTruthMatch);
      FillElement(slcVars["slc_true_event_type_Np0pi"], slcCounter, (int) NC::kFailedTruthMatch);
      FillElement(slcVars["slc_true_event_type_Xp0pi"], slcCounter, (int) NC::kFailedTruthMatch);
      FillElement(slcVars["slc_true_event_type_cc"], slcCounter, (int) CC::kFailedTruthMatch);
    }
}

float sbnd::DPreco::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

float sbnd::DPreco::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (fHitsMap[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(fHitsMap[trackID]);
}

void sbnd::DPreco::SelectSlice(const int counter)
{
  bool is_clear_cosmic;
  AccessElement(slcVars["slc_is_clear_cosmic"], counter, is_clear_cosmic);

  bool is_fv;
  AccessElement(slcVars["slc_is_fv"], counter, is_fv);

  float crumbs;
  AccessElement(slcVars["slc_crumbs_nc_score"], counter, crumbs);
  const bool passes_crumbs_incl  = crumbs > -0.005;
  const bool passes_crumbs_0p0pi = crumbs > -0.005;
  const bool passes_crumbs_Np0pi = crumbs > 0.235;

  float crumbs_ccnumu;
  AccessElement(slcVars["slc_crumbs_ccnumu_score"], counter, crumbs_ccnumu);
  const bool passes_crumbs_cc = crumbs_ccnumu > 0;

  int nrazzledmuons;
  AccessElement(slcVars["slc_n_primary_razzled_muons"], counter, nrazzledmuons);
  const bool passes_razzled_muons    = nrazzledmuons == 0;
  const bool passes_razzled_muons_cc = nrazzledmuons == 1;

  size_t npfps;
  AccessElement(slcVars["slc_n_pfps"], counter, npfps);
  const bool passes_pfps = npfps > 1;

  int nrazzledphotons;
  AccessElement(slcVars["slc_n_primary_razzled_photons"], counter, nrazzledphotons);
  const bool passes_razzled_photons = nrazzledphotons > 1;

  bool bestpzcgoodkinematics;
  AccessElement(slcVars["slc_best_corr_pzc_good_kinematics"], counter, bestpzcgoodkinematics);

  double opt0frac;
  AccessElement(slcVars["slc_opt0_fracPE"], counter, opt0frac);
  const bool passes_opt0frac_incl  = opt0frac < 0.512 && opt0frac > -0.676;
  const bool passes_opt0frac_0p0pi = opt0frac < 0.4 && opt0frac > -0.676;
  const bool passes_opt0frac_Np0pi = opt0frac < 0.512 && opt0frac > -0.44;
  const bool passes_opt0frac_cc    = opt0frac < 0.836 && opt0frac > -0.376;

  double opt0score;
  AccessElement(slcVars["slc_opt0_score"], counter, opt0score);
  const bool passes_opt0score_incl  = opt0score > 45;
  const bool passes_opt0score_0p0pi = opt0score > 110;
  const bool passes_opt0score_Np0pi = opt0score > 80;
  const bool passes_opt0score_cc    = opt0score > 210;

  int nrazzledpions;
  AccessElement(slcVars["slc_n_primary_razzled_pions_thresh"], counter, nrazzledpions);
  const bool passes_razzled_pions = nrazzledpions == 0;

  int nrazzledprotons;
  AccessElement(slcVars["slc_n_primary_razzled_protons_thresh"], counter, nrazzledprotons);

  bool allothertrkscontained;
  AccessElement(slcVars["slc_all_other_trks_contained"], counter, allothertrkscontained);

  const bool sel_incl = !is_clear_cosmic && is_fv && passes_crumbs_incl && passes_razzled_muons && passes_pfps
    && passes_razzled_photons && bestpzcgoodkinematics && passes_opt0frac_incl && passes_opt0score_incl
    && allothertrkscontained;
  FillElement(slcVars["slc_sel_incl"], counter, sel_incl);

  const bool sel_0p0pi = !is_clear_cosmic && is_fv && passes_crumbs_0p0pi && passes_razzled_muons && passes_pfps
    && passes_razzled_photons && bestpzcgoodkinematics && passes_opt0frac_0p0pi && passes_opt0score_0p0pi
    && passes_razzled_pions && nrazzledprotons == 0;
  FillElement(slcVars["slc_sel_0p0pi"], counter, sel_0p0pi);

  const bool sel_Np0pi = !is_clear_cosmic && is_fv && passes_crumbs_Np0pi && passes_razzled_muons && passes_pfps
    && passes_razzled_photons && bestpzcgoodkinematics && passes_opt0frac_Np0pi && passes_opt0score_Np0pi
    && passes_razzled_pions && nrazzledprotons > 0;
  FillElement(slcVars["slc_sel_Np0pi"], counter, sel_Np0pi);

  const bool sel_1p0pi = sel_incl && passes_razzled_pions && nrazzledprotons == 1;
  FillElement(slcVars["slc_sel_1p0pi"], counter, sel_1p0pi);

  const bool sel_Xp0pi = sel_incl && passes_razzled_pions;
  FillElement(slcVars["slc_sel_Xp0pi"], counter, sel_Xp0pi);

  const bool sel_cc = !is_clear_cosmic && is_fv && passes_crumbs_cc && passes_razzled_muons_cc && passes_pfps
    && passes_razzled_photons && bestpzcgoodkinematics && passes_opt0frac_cc && passes_opt0score_cc;
  FillElement(slcVars["slc_sel_cc"], counter, sel_cc);
}



void sbnd::DPreco::ResetSubRunVars()
{
  _pot = 0.; _spills = 0; _ngenevts = 0;
}

void sbnd::DPreco::ResetEventVars()
{
  _run = -1; _subrun = -1; _event  = -1;

  _n_nu = 0; _n_slc  = 0;
}

int sbnd::DPreco::GetTotalGenEvents(const art::Event &e)
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

art::Ptr<recob::PFParticle> sbnd::DPreco::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
{
  for(auto pfp : pfps)
    if(pfp->IsPrimary())
      return pfp;

  return art::Ptr<recob::PFParticle>();
}

double sbnd::DPreco::CorrectEnergy(const double &energy)
{
  //const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return 1;//energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}

void sbnd::DPreco::ResizeVectors(VecVarMap &map, const int size)
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

void sbnd::DPreco::ResizeSubVectors(VecVarMap &map, const std::string &subname, const int pos, const int size)
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

void sbnd::DPreco::TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB, const int posA, const int posB)
{
  std::string name = var->Name();
  name.erase(name.find(prefixA), prefixA.size());

  switch(var->IdentifyVar())
    {
    case kBool:
      bool bool_val;
      AccessElement(vars[prefixB + name], posB, bool_val);
      FillElement(var, posA, bool_val);
      break;
    case kInt:
      int int_val;
      AccessElement(vars[prefixB + name], posB, int_val);
      FillElement(var, posA, int_val);
      break;
    case kUInt:
      size_t uint_val;
      AccessElement(vars[prefixB + name], posB, uint_val);
      FillElement(var, posA, uint_val);
      break;
    case kFloat:
      float float_val;
      AccessElement(vars[prefixB + name], posB, float_val);
      FillElement(var, posA, float_val);
      break;
    case kDouble:
      double double_val;
      AccessElement(vars[prefixB + name], posB, double_val);
      FillElement(var, posA, double_val);
      break;
    case kUnknownVar:
      break;
    }
}

void sbnd::DPreco::TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB,
                                             const int posA0, const int posA1, const int posB0, const int posB1)
{
  std::string name = var->Name();
  name.erase(name.find(prefixA), prefixA.size());

  switch(var->IdentifyVar())
    {
    case kBool:
      bool bool_val;
      AccessElement(vars[prefixB + name], posB0, posB1, bool_val);
      FillElement(var, posA0, posA1, bool_val);
      break;
    case kInt:
      int int_val;
      AccessElement(vars[prefixB + name], posB0, posB1, int_val);
      FillElement(var, posA0, posA1, int_val);
      break;
    case kUInt:
      size_t uint_val;
      AccessElement(vars[prefixB + name], posB0, posB1, uint_val);
      FillElement(var, posA0, posA1, uint_val);
      break;
    case kFloat:
      float float_val;
      AccessElement(vars[prefixB + name], posB0, posB1, float_val);
      FillElement(var, posA0, posA1, float_val);
      break;
    case kDouble:
      double double_val;
      AccessElement(vars[prefixB + name], posB0, posB1, double_val);
      FillElement(var, posA0, posA1, double_val);
      break;
    case kUnknownVar:
      break;
    }
}

void sbnd::DPreco::TransferElement(VecVar *var, VecVarMap &vars, const std::string prefixA, const std::string prefixB,
                                             const int posA, const int posB0, const int posB1)
{
  std::string name = var->Name();
  name.erase(name.find(prefixA), prefixA.size());

  switch(var->IdentifyVar())
    {
    case kBool:
      bool bool_val;
      AccessElement(vars[prefixB + name], posB0, posB1, bool_val);
      FillElement(var, posA, bool_val);
      break;
    case kInt:
      int int_val;
      AccessElement(vars[prefixB + name], posB0, posB1, int_val);
      FillElement(var, posA, int_val);
      break;
    case kUInt:
      size_t uint_val;
      AccessElement(vars[prefixB + name], posB0, posB1, uint_val);
      FillElement(var, posA, uint_val);
      break;
    case kFloat:
      float float_val;
      AccessElement(vars[prefixB + name], posB0, posB1, float_val);
      FillElement(var, posA, float_val);
      break;
    case kDouble:
      double double_val;
      AccessElement(vars[prefixB + name], posB0, posB1, double_val);
      FillElement(var, posA, double_val);
      break;
    case kUnknownVar:
      break;
    }
}

template<typename T>
void sbnd::DPreco::FillElement(VecVar *vec, const int pos, const T value)
{
  dynamic_cast<InhVecVar<T>*>(vec)->SetVal(pos, value);
}

template<typename T>
void sbnd::DPreco::FillElement(VecVar *vec, const int pos, const std::vector<T> &value)
{
  dynamic_cast<InhVecVecVar<T>*>(vec)->SetVal(pos, value);
}

template<typename T>
void sbnd::DPreco::FillElement(VecVar *vec, const int posA, const int posB, const T value)
{
  dynamic_cast<InhVecVecVar<T>*>(vec)->SetVal(posA, posB, value);
}

template<typename T>
void sbnd::DPreco::AccessElement(VecVar *vec, const int pos, T &value)
{
  value = dynamic_cast<InhVecVar<T>*>(vec)->GetVal(pos);
}

template<typename T>
void sbnd::DPreco::AccessElement(VecVar *vec, const int posA, const int posB, T &value)
{
  value = dynamic_cast<InhVecVecVar<T>*>(vec)->GetVal(posA, posB);
}

template<typename T>
void sbnd::DPreco::GetVar(VecVar *vec, std::vector<T> &var)
{
  var = dynamic_cast<InhVecVar<T>*>(vec)->Var();
}

template<typename T>
void sbnd::DPreco::GetVar(VecVar *vec, std::vector<std::vector<T>> &var)
{
  var = dynamic_cast<InhVecVecVar<T>*>(vec)->Var();
}

DEFINE_ART_MODULE(sbnd::DPreco)
