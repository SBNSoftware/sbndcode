////////////////////////////////////////////////////////////////////////
// Class:       HNLPiZeroAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        HNLPiZeroAnalysis_module.cc
//
// Generated at Thu Sep 28 09:58:07 2023 by Vu Chi Lan Nguyen using cetskelgen
// from  version .
// This code was heavily stolen from Henry Lay ;)
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
#include "art_root_io/TFileService.h"

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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CNNScore.h"
#include "sbnobj/Common/Reco/Stub.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

#include <numeric>
#include <fstream>
#include <iostream>

constexpr int def_int       = std::numeric_limits<int>::min();
constexpr size_t def_size   = std::numeric_limits<size_t>::max();
constexpr float def_float   = -std::numeric_limits<float>::max();
constexpr double def_double = -std::numeric_limits<double>::max();
constexpr double light_speed_cm_ns = 29.9792458; // cm/ns
constexpr double light_speed_cm_us = 29979.2458; // cm/ns
constexpr double pi0_mass = 134.9770; // MeV/c^2
					   
template <typename T,
  typename TIter = decltype(std::begin(std::declval<T>())),
  typename = decltype(std::end(std::declval<T>()))>
  constexpr auto enumerate(T && iterable)
{
  struct iterator
  {
    size_t i;
    TIter iter;
    bool operator != (const iterator & other) const { return iter != other.iter; }
    void operator ++ () { ++i; ++iter; }
    auto operator * () const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper
  {
    T iterable;
    auto begin() { return iterator{ 0, std::begin(iterable) }; }
    auto end() { return iterator{ 0, std::end(iterable) }; }
  };
  return iterable_wrapper{ std::forward<T>(iterable) };
}

enum EventType
{
  kHNL,
  kHNLNonFV,
  kHNLDirt,
  kNCPiZero,
  kOtherNC,
  kCCNuMu,
  kCCNuE,
  kNuNonFV,
  kNuDirt,
  kCosmic,
  kUnknownEv = -1
};

namespace sbnd {
  class HNLPiZeroAnalysis;
}

class sbnd::HNLPiZeroAnalysis : public art::EDAnalyzer {

public:
  explicit HNLPiZeroAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HNLPiZeroAnalysis(HNLPiZeroAnalysis const&) = delete;
  HNLPiZeroAnalysis(HNLPiZeroAnalysis&&) = delete;
  HNLPiZeroAnalysis& operator=(HNLPiZeroAnalysis const&) = delete;
  HNLPiZeroAnalysis& operator=(HNLPiZeroAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(const art::SubRun& sr);
  virtual void endSubRun(const art::SubRun& sr);

  int GetTotalGenEvents(const art::Event &e);
  void ResetSubRunVars();
  void ResetEventVars();
  
  void ResizeMCTruth1DVector(const int col);
  void ResizeSlice1DVector(const int col);
  void ResizeSlice2DVectorRow(const int row);
  void ResizeSlice2DVectorCol(const int row, const int col);
  void ResizeWeight2DVectorRow(const int row);
  void ResizeWeight2DVectorCol(const int row);

  void ClearMaps();
  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                 const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);
  void AnalyseMCTruthHandle(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);
  void AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct, const int mctCounter);
  void AnalyseMeVPrtlTruth(const art::Event &e, const art::Handle<std::vector<evgen::ldm::MeVPrtlTruth>> mevptHandle);
  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                     const art::Handle<std::vector<recob::Track>> &trackHandle,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle);
  
  void AnalysePFPs(const art::Event &e, const int slcCounter,
		   const art::Ptr<recob::PFParticle> &prim, 
		   const art::Ptr<recob::Vertex> &vtx, 
		   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, 
		   const art::Handle<std::vector<recob::Track>> &trackHandle,
		   const art::Handle<std::vector<recob::Shower>> &showerHandle);
  
  void AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, 
                    const int slcCounter, const int pfpCounter,
	            const art::Handle<std::vector<recob::Track>> &trackHandle);
  
  void AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower,
		     const int slcCounter, const int pfpCounter,
		     const art::Handle<std::vector<recob::Shower>> &showerHandle, 
		     const art::Ptr<recob::Vertex> &vtx,
                     const std::vector<art::Ptr<recob::Hit>> &hits);
  
  void AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, const int slcCounter,
                         const art::Handle<std::vector<recob::Slice>> &sliceHandle);
  
  void AnalyseSliceMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct, const int slcCounter);

  void AnalysePiZero(const int slcCounter);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

  float Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  float Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  void ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter); 
  void ExtractDazzle(const art::Ptr<sbn::MVAPID> &dazzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<anab::Calorimetry> &calo, const int slcCounter, const int pfpCounter);
  void ExtractChi2PID(const art::Ptr<anab::ParticleID> &chi2pid, const int slcCounter, const int pfpCounter);
  void ExtractMCS(const art::Ptr<recob::MCSFitResult> &mcs, const int slcCounter, const int pfpCounter);
  void ExtractStoppingChi2(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2, const int slcCounter, const int pfpCounter);
  void ExtractRazzle(const art::Ptr<sbn::MVAPID> &razzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                   const std::vector<art::Ptr<recob::Hit>> &hits);  
  void ExtractCNNScores(const sbn::PFPCNNScore *cnnscore, const int slcCounter, const int pfpCounter);

private:
  
  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  // Module Labels
  art::InputTag fPOTModuleLabel, fMeVPrtlTruthModuleLabel, fMCParticleModuleLabel,
	  fHitModuleLabel, fPFParticleModuleLabel, fSliceModuleLabel, 
	  fVertexModuleLabel, fTrackModuleLabel, fShowerModuleLabel,
	  fCRUMBSModuleLabel, fOpT0ModuleLabel,fDazzleModuleLabel, 
	  fCaloModuleLabel, fMCSModuleLabel, fChi2ModuleLabel, 
	  fRangeModuleLabel, fClosestApproachModuleLabel, fStoppingChi2ModuleLabel,
	  fRazzleModuleLabel, fCosmicDistModuleLabel, fShowerDensityFitModuleLabel,
	  fShowerTrackFitModuleLabel, fCNNScoreModuleLabel, fRazzledModuleLabel,
	  fSpacePointModuleLabel, fStubModuleLabel;
  
  std::vector<art::InputTag> fEventWeightModuleLabels;
  std::vector<std::string> fOpFlashesModuleLabel;
 
  bool fBeamOff; 
  bool fMeVPrtl;
  bool fGenie;
  bool fDebug;
  bool fPiZeroAna;

  // Map
  std::map<int, int> fHitsMap; //track ID , nHits
  std::map<const art::Ptr<simb::MCTruth>, int> fMCTruthHitsMap; //mctruth id, nHits 
  std::map<int, art::Ptr<recob::PFParticle>> fPFPMap; //pfp ID, pfp
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap; 

  // Event Tree
  TTree *fEventTree;
  int _run, _subrun, _event;

  // Event Tree: MeVPrtl Truth
  int _n_hnl;
  std::vector<double> mevprtl_decay_pos_x, mevprtl_decay_pos_y, mevprtl_decay_pos_z, mevprtl_decay_pos_t;
  std::vector<double> mevprtl_mom_x, mevprtl_mom_y, mevprtl_mom_z, mevprtl_mom_e;
  std::vector<double> mevprtl_mass;
  std::vector<double> mevprtl_flux_weight, mevprtl_ray_weight, mevprtl_decay_weight;
  std::vector<double> mevprtl_total_decay_width, mevprtl_total_mean_lifetime, mevprtl_total_mean_distance, mevprtl_allowed_decay_fraction;
  std::vector<double> mevprtl_C1, mevprtl_C2, mevprtl_C3, mevprtl_C4, mevprtl_C5;

  // Event Tree: MCTruth Neutrino
  int _n_mctruth, _n_nu;
  std::vector<size_t> 	nu_mctruth_id;
  std::vector<int> 	nu_event_type;
  std::vector<bool> 	nu_signal;
  std::vector<float> 	nu_en_dep;
  std::vector<int> 	nu_ccnc;
  std::vector<int> 	nu_mode;
  std::vector<int> 	nu_int_type;
  std::vector<double> 	nu_w;
  std::vector<double> 	nu_x;
  std::vector<double> 	nu_y;
  std::vector<double> 	nu_q_sqr;
  std::vector<double> 	nu_pt;
  std::vector<double> 	nu_theta;
  std::vector<double> 	nu_e;
  std::vector<double> 	nu_vtx_x, nu_vtx_y, nu_vtx_z, nu_vtx_t;
  std::vector<int>	nu_n_protons, nu_n_neutrons, nu_n_charged_pions, nu_n_neutral_pions, nu_n_photons, nu_n_others;

  // Event Tree: Slice -- 1D vector
  int _n_slc;
  std::vector<int> 	slc_id;
  std::vector<size_t> 	slc_n_pfps;
  std::vector<bool> 	slc_is_clear_cosmics;
  std::vector<size_t>  	slc_primary_pfp_id;
  std::vector<int> 	slc_primary_pfp_pdg;
  std::vector<int>	slc_n_primary_daughters;
  std::vector<double>	slc_vtx_x, slc_vtx_y, slc_vtx_z;
  std::vector<bool>	slc_is_fv;
  std::vector<float> 	slc_crumbs_score, slc_crumbs_nc_score, slc_crumbs_ccnue_score;
  std::vector<double>	slc_opt0_time, slc_opt0_score, slc_opt0_measPE, slc_opt0_hypoPE, slc_opt0_frac, slc_opt0_time_corrected_Z_pandora;
  std::vector<int>	slc_n_trks, slc_n_shws;
  std::vector<int> 	slc_n_primary_trks, slc_n_primary_shws;
  std::vector<int>	slc_n_dazzle_muons, slc_n_dazzle_pions, slc_n_dazzle_pions_thresh, slc_n_dazzle_protons, slc_n_dazzle_protons_thresh, slc_n_dazzle_other;
  std::vector<int>	slc_n_primary_dazzle_muons, slc_n_primary_dazzle_pions, slc_n_primary_dazzle_pions_thresh, slc_n_primary_dazzle_protons, slc_n_primary_dazzle_protons_thresh, slc_n_primary_dazzle_other;
  std::vector<int>	slc_n_razzle_electrons, slc_n_razzle_photons, slc_n_razzle_other;
  std::vector<int>	slc_n_razzled_electrons, slc_n_razzled_muons, slc_n_razzled_photons, slc_n_razzled_pions, slc_n_razzled_pions_thresh, slc_n_razzled_protons, slc_n_razzled_protons_thresh;
  std::vector<int>	slc_n_primary_razzle_electrons, slc_n_primary_razzle_photons, slc_n_primary_razzle_other;
  std::vector<int>	slc_n_primary_razzled_electrons, slc_n_primary_razzled_muons, slc_n_primary_razzled_photons, slc_n_primary_razzled_pions, slc_n_primary_razzled_pions_thresh, slc_n_primary_razzled_protons, slc_n_primary_razzled_protons_thresh;
  std::vector<int>	slc_n_stub;
  std::vector<float>	slc_total_shower_E, slc_total_track_E;
 
  // Event Tree: Slice Truth
  std::vector<float> 	slc_comp, slc_pur;
  std::vector<size_t> 	slc_true_mctruth_id;
  std::vector<int> 	slc_true_event_type;
  std::vector<float> 	slc_true_en_dep;
  std::vector<float> 	slc_true_vtx_x, slc_true_vtx_y, slc_true_vtx_z, slc_true_vtx_t, slc_true_vtx_t_corrected_Z;
 
  // Event Tree: Slice Truth Pi0
  std::vector<float>	slc_true_pizero_mom, slc_true_pizero_cosThetaZ;

  // Event Tree: Slice -> PFP -- 2D vector
  std::vector<std::vector<size_t>>	slc_pfp_id;
  std::vector<std::vector<int>> 	slc_pfp_pdg;
  std::vector<std::vector<int>>  	slc_pfp_n_children;
  std::vector<std::vector<bool>>	slc_pfp_primary;
  std::vector<std::vector<bool>>	slc_pfp_primary_child;
  std::vector<std::vector<size_t>>	slc_pfp_n_hits;
  std::vector<std::vector<size_t>>	slc_pfp_n_sps;
  std::vector<std::vector<float>> 	slc_pfp_track_score;
  std::vector<std::vector<bool>> 	slc_pfp_good_track;
  std::vector<std::vector<bool>> 	slc_pfp_good_shower;
  std::vector<std::vector<float>>	slc_pfp_comp;
  std::vector<std::vector<float>>	slc_pfp_pur;
  std::vector<std::vector<float>> 	slc_pfp_cnnscore_track, slc_pfp_cnnscore_shower, slc_pfp_cnnscore_noise, slc_pfp_cnnscore_michel, slc_pfp_cnnscore_endmichel;
  std::vector<std::vector<int>>		slc_pfp_cnnscore_nclusters;
  std::vector<std::vector<float>>	slc_pfp_razzled_electron_score, slc_pfp_razzled_muon_score, slc_pfp_razzled_photon_score, slc_pfp_razzled_pion_score, slc_pfp_razzled_proton_score;
  std::vector<std::vector<int>> 	slc_pfp_razzled_pdg;

  // Event Tree: Slice -> PFP -> True MCParticle -- 2D vector
  std::vector<std::vector<int>>		slc_pfp_true_trackid;
  std::vector<std::vector<int>>        	slc_pfp_true_pdg;
  std::vector<std::vector<float>>      	slc_pfp_true_energy;
  std::vector<std::vector<float>>      	slc_pfp_true_p_x;
  std::vector<std::vector<float>>      	slc_pfp_true_p_y;
  std::vector<std::vector<float>>      	slc_pfp_true_p_z;

  // Event Tree: Slice -> PFP -> Track -- 2D vector
  std::vector<std::vector<double>>	slc_pfp_track_start_x, slc_pfp_track_start_y, slc_pfp_track_start_z;
  std::vector<std::vector<double>>	slc_pfp_track_end_x, slc_pfp_track_end_y, slc_pfp_track_end_z;
  std::vector<std::vector<double>> 	slc_pfp_track_dir_x, slc_pfp_track_dir_y, slc_pfp_track_dir_z;
  std::vector<std::vector<double>> 	slc_pfp_track_length;
  std::vector<std::vector<bool>> 	slc_pfp_track_contained;
  std::vector<std::vector<float>>   	slc_pfp_track_dazzle_muon_score, slc_pfp_track_dazzle_pion_score, slc_pfp_track_dazzle_proton_score, slc_pfp_track_dazzle_other_score;
  std::vector<std::vector<int>>		slc_pfp_track_dazzle_pdg;
  std::vector<std::vector<float>>	slc_pfp_track_ke, slc_pfp_track_charge;
  std::vector<std::vector<float>>	slc_pfp_track_chi2_muon, slc_pfp_track_chi2_pion, slc_pfp_track_chi2_kaon, slc_pfp_track_chi2_proton;	
  std::vector<std::vector<int>>		slc_pfp_track_chi2_pdg;
  std::vector<std::vector<float>> 	slc_pfp_track_mcs_mean_scatter, slc_pfp_track_mcs_max_scatter_ratio;
  std::vector<std::vector<float>> 	slc_pfp_track_range_p;
  std::vector<std::vector<float>>	slc_pfp_track_closest_approach_mean_dca;
  std::vector<std::vector<float>>	slc_pfp_track_stopping_dedx_chi2_ratio, slc_pfp_track_stopping_dedx_pol0_fit;
  std::vector<std::vector<float>> 	slc_pfp_track_theta, slc_pfp_track_phi;

  // Event Tree: Slice -> PFP -> Shower -- 2D vector
  std::vector<std::vector<double>>	slc_pfp_shower_start_x, slc_pfp_shower_start_y, slc_pfp_shower_start_z;
  std::vector<std::vector<double>>	slc_pfp_shower_end_x, slc_pfp_shower_end_y, slc_pfp_shower_end_z;
  std::vector<std::vector<double>>	slc_pfp_shower_conv_gap;
  std::vector<std::vector<double>>	slc_pfp_shower_dir_x, slc_pfp_shower_dir_y, slc_pfp_shower_dir_z;
  std::vector<std::vector<double>>	slc_pfp_shower_length;
  std::vector<std::vector<double>>	slc_pfp_shower_open_angle;
  std::vector<std::vector<double>>	slc_pfp_shower_energy, slc_pfp_shower_dedx;
  std::vector<std::vector<double>>	slc_pfp_shower_sqrt_energy_density, slc_pfp_shower_modified_hit_density;
  std::vector<std::vector<float>>	slc_pfp_shower_razzle_electron_score, slc_pfp_shower_razzle_photon_score, slc_pfp_shower_razzle_other_score;
  std::vector<std::vector<int>>		slc_pfp_shower_razzle_pdg;
  std::vector<std::vector<float>>	slc_pfp_shower_cosmic_dist;
  std::vector<std::vector<double>>	slc_pfp_shower_track_length, slc_pfp_shower_track_width;
  std::vector<std::vector<double>>	slc_pfp_shower_density_grad, slc_pfp_shower_density_pow;
  std::vector<std::vector<bool>> 	slc_pfp_shower_contained;
  std::vector<std::vector<float>> 	slc_pfp_shower_theta, slc_pfp_shower_phi;

  // Event Tree: OpFlashes
  std::vector<int> _flash_id;
  std::vector<double> _flash_time;
  std::vector<double> _flash_total_pe;
  std::vector<std::vector<double>> _flash_pe_v;
  std::vector<double> _flash_y;
  std::vector<double> _flash_yerr ;
  std::vector<double> _flash_z;
  std::vector<double> _flash_zerr;
  std::vector<int> _flash_tpc;

  // Event Tree: Pizero Variables
  std::vector<bool> slc_pizero_good;
  std::vector<int> slc_photon0_pfp_id;
  std::vector<int> slc_photon1_pfp_id;
  std::vector<int> slc_photon2_pfp_id;
  std::vector<float> slc_photon0_energy;
  std::vector<float> slc_photon1_energy;
  std::vector<float> slc_photon2_energy;
  std::vector<float> slc_pizero_open_angle;
  std::vector<float> slc_pizero_mass;
  std::vector<float> slc_pizero_mom;
  std::vector<float> slc_pizero_CosThetaZ;

  //Event Tree: Flux Weight
  std::vector<std::vector<float>> slc_flux_weight_expskin, slc_flux_weight_horncurrent, slc_flux_weight_kminus,
  	  			  slc_flux_weight_kplus, slc_flux_weight_kzero, slc_flux_weight_nucleoninexsec,
  				  slc_flux_weight_nucleonqexsec, slc_flux_weight_nucleontotxsec, slc_flux_weight_piminus,
  				  slc_flux_weight_pioninexsex, slc_flux_weight_pionqexsec, slc_flux_weight_piontotxsec,
  				  slc_flux_weight_piplus;
  std::vector<std::vector<float>> slc_flux_weight_total; //store multiplicative of all the flux weights

  //Event Tree: xsection
  std::vector<std::vector<float>> slc_xsec_unisim_DecayAngMEC, 
	  			  slc_xsec_unisim_ThetaDelta2NRad, 
	  			  slc_xsec_unisim_Theta_Delta2Npi, 
				  slc_xsec_unisim_VecFFCCQEshape,
                                  slc_xsec_multisigma_CoulombCCQE,
                                  slc_xsec_multisigma_NonRESBGvbarnCC1pi,
                                  slc_xsec_multisigma_NonRESBGvbarnCC2pi,
                                  slc_xsec_multisigma_NonRESBGvbarnNC1pi,
                                  slc_xsec_multisigma_NonRESBGvbarnNC2pi,
                                  slc_xsec_multisigma_NonRESBGvbarpCC1pi,
                                  slc_xsec_multisigma_NonRESBGvbarpCC2pi,
                                  slc_xsec_multisigma_NonRESBGvbarpNC1pi,
                                  slc_xsec_multisigma_NonRESBGvbarpNC2pi,
                                  slc_xsec_multisigma_NonRESBGvnCC1pi,
                                  slc_xsec_multisigma_NonRESBGvnCC2pi,
                                  slc_xsec_multisigma_NonRESBGvnNC1pi,
                                  slc_xsec_multisigma_NonRESBGvnNC2pi,
                                  slc_xsec_multisigma_NonRESBGvpCC1pi,
                                  slc_xsec_multisigma_NonRESBGvpCC2pi,
                                  slc_xsec_multisigma_NonRESBGvpNC1pi,
                                  slc_xsec_multisigma_NonRESBGvpNC2pi,
                                  slc_xsec_multisigma_NormCCMEC,
                                  slc_xsec_multisigma_NormNCMEC,
                                  slc_xsec_multisigma_RDecBR1eta,
                                  slc_xsec_multisigma_RDecBR1gamma,
                                  slc_xsec_multisigma_RPA_CCQE,
				  slc_xsec_multisigma_NormNCCOH,
				  slc_xsec_multisigma_NormCCCOH;
  std::vector<std::vector<float>> slc_xsec_multisim_ZExpA,
	     			  slc_xsec_multisim_NCEL,
	      			  slc_xsec_multisim_CCRES,
	      			  slc_xsec_multisim_NCRES,
	      			  slc_xsec_multisim_DISBY,
	      			  slc_xsec_multisim_FSI_pi, 
	      			  slc_xsec_multisim_FSI_N; 
  std::vector<std::vector<float>> slc_xsec_multisim_total;
  
  std::vector<std::vector<float>> slc_geant4_multisim_reinteractions;

  //Sub Run Tree
  TTree *fSubRunTree;
  double _pot;
  int _spills, _ngenevts;

  //CSV File
  std::ofstream outputFile;
};

sbnd::HNLPiZeroAnalysis::HNLPiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Module Labels
  fPOTModuleLabel              	= p.get<art::InputTag>("POTModuleLabel", "generator");
  fMeVPrtlTruthModuleLabel      = p.get<art::InputTag>("MeVPrtlTruthModuleLabel", "generator");
  fMCParticleModuleLabel        = p.get<art::InputTag>("MCParticleModuleLabel", "largeant");
  fHitModuleLabel 		= p.get<art::InputTag>("HitModuleLabel", "gaushit");
  fPFParticleModuleLabel 	= p.get<art::InputTag>("PFParticleModuleLabel", "pandoraSCE");
  fSliceModuleLabel 		= p.get<art::InputTag>("SlideModuleLabel", "pandoraSCE");
  fTrackModuleLabel		= p.get<art::InputTag>("TrackModuleLabel", "pandoraTrack");
  fShowerModuleLabel 		= p.get<art::InputTag>("ShowerModuleLabel", "pandoraShower");
  fOpFlashesModuleLabel = p.get<std::vector<std::string>>("OpFlashesModuleLabel",   {"opflashtpc0", "opflashtpc1"});
  fVertexModuleLabel 		= p.get<art::InputTag>("VertexModuleLabel", "pandoraSCE");
  fCRUMBSModuleLabel 		= p.get<art::InputTag>("CRUMBSModuleLabelr", "crumbs");
  fOpT0ModuleLabel 		= p.get<art::InputTag>("OpT0ModuleLabel", "opt0finder");
  fDazzleModuleLabel 		= p.get<art::InputTag>("DazzleModuleLabel", "dazzle");
  fCaloModuleLabel 		= p.get<art::InputTag>("CaloModuleLabel", "pandoraCalo");
  fMCSModuleLabel 		= p.get<art::InputTag>("MCSModuleLabel", "pandoraTrackMCS:muon");
  fChi2ModuleLabel 		= p.get<art::InputTag>("Chi2ModuleLabel", "pandoraPid");
  fRangeModuleLabel 		= p.get<art::InputTag>("RangeModuleLabel", "pandoraTrackRange:muon");
  fClosestApproachModuleLabel 	= p.get<art::InputTag>("ClosestApproachModuleLabel", "pandoraTrackClosestApproach");
  fStoppingChi2ModuleLabel 	= p.get<art::InputTag>("StoppingChi2ModuleLabel", "pandoraTrackStoppingChi2");
  fRazzleModuleLabel 		= p.get<art::InputTag>("RazzleModuleLabel", "razzle");
  fCosmicDistModuleLabel 	= p.get<art::InputTag>("CosmicDistModuleLabel", "pandoraShowerCosmicDist");
  fShowerTrackFitModuleLabel 	= p.get<art::InputTag>("ShowerTrackFitModuleLabel", "pandoraShowerSelectionVars");
  fShowerDensityFitModuleLabel 	= p.get<art::InputTag>("ShowerDensityFitModuleLabel", "pandoraShowerSelectionVars");
  fCNNScoreModuleLabel        	= p.get<art::InputTag>("CNNScoreModuleLabel", "cnnid");
  fStubModuleLabel 		= p.get<art::InputTag>("StubModuleLabel", "vertexStub");
  fRazzledModuleLabel           = p.get<art::InputTag>("RazzledModuleLabel", "razzled");
  fSpacePointModuleLabel        = p.get<art::InputTag>("SpacePointModuleLabel", "pandoraSCE");
  fEventWeightModuleLabels 	= p.get<std::vector<art::InputTag>>("EventWeightModuleLabels");
  fGenie 			= p.get<bool>("Genie", true);
  fBeamOff			= p.get<bool>("BeamOff", false);
  fMeVPrtl			= p.get<bool>("MeVPrtl", false);
  fDebug			= p.get<bool>("Debug", false);
  fPiZeroAna			= p.get<bool>("PiZeroAna", false);

  art::ServiceHandle<art::TFileService> tfs;

  //Sub Run Tree
  fSubRunTree = tfs->make<TTree>("subruns","");

  fSubRunTree->Branch("run", &_run);
  fSubRunTree->Branch("subrun", &_subrun);
  fSubRunTree->Branch("pot", &_pot);
  fSubRunTree->Branch("spills", &_spills);
  fSubRunTree->Branch("ngenevts", &_ngenevts);

  //Event Tree
  fEventTree = tfs->make<TTree>("events", "");

  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);

  //Event Tree: MeVPrtl Truth
  fEventTree->Branch("n_hnl",&_n_hnl);
  fEventTree->Branch("mevprtl_decay_pos_x", &mevprtl_decay_pos_x);
  fEventTree->Branch("mevprtl_decay_pos_y", &mevprtl_decay_pos_y);
  fEventTree->Branch("mevprtl_decay_pos_z", &mevprtl_decay_pos_z);
  fEventTree->Branch("mevprtl_decay_pos_t", &mevprtl_decay_pos_t);
  fEventTree->Branch("mevprtl_mom_x", &mevprtl_mom_x); 
  fEventTree->Branch("mevprtl_mom_y", &mevprtl_mom_y); 
  fEventTree->Branch("mevprtl_mom_z", &mevprtl_mom_z); 
  fEventTree->Branch("mevprtl_mom_e", &mevprtl_mom_e); 
  fEventTree->Branch("mevprtl_mass", &mevprtl_mass);
  fEventTree->Branch("mevprtl_flux_weight", &mevprtl_flux_weight);
  fEventTree->Branch("mevprtl_ray_weight", &mevprtl_ray_weight);
  fEventTree->Branch("mevprtl_decay_weight", &mevprtl_decay_weight);
  fEventTree->Branch("mevprtl_total_decay_width", &mevprtl_total_decay_width);
  fEventTree->Branch("mevprtl_total_mean_lifetime", &mevprtl_total_mean_lifetime);
  fEventTree->Branch("mevprtl_total_mean_distance", &mevprtl_total_mean_distance);
  fEventTree->Branch("mevprtl_allowed_decay_fraction", &mevprtl_allowed_decay_fraction);
  fEventTree->Branch("mevprtl_C1", &mevprtl_C1);
  fEventTree->Branch("mevprtl_C2", &mevprtl_C2);
  fEventTree->Branch("mevprtl_C3", &mevprtl_C3);
  fEventTree->Branch("mevprtl_C4", &mevprtl_C4);
  fEventTree->Branch("mevprtl_C5", &mevprtl_C5);

  //Event Tree: MCTruth 
  fEventTree->Branch("n_nu",&_n_nu);
  fEventTree->Branch("nu_mctruth_id", &nu_mctruth_id); 
  fEventTree->Branch("nu_event_type", &nu_event_type); 
  fEventTree->Branch("nu_signal", &nu_signal); 
  fEventTree->Branch("nu_en_dep", &nu_en_dep); 
  fEventTree->Branch("nu_ccnc", &nu_ccnc); 
  fEventTree->Branch("nu_mode", &nu_mode); 
  fEventTree->Branch("nu_int_type", &nu_int_type); 
  fEventTree->Branch("nu_w", &nu_w); 
  fEventTree->Branch("nu_x", &nu_x); 
  fEventTree->Branch("nu_y", &nu_y); 
  fEventTree->Branch("nu_q_sqr", &nu_q_sqr); 
  fEventTree->Branch("nu_pt", &nu_pt); 
  fEventTree->Branch("nu_theta", &nu_theta); 
  fEventTree->Branch("nu_e", &nu_e); 
  fEventTree->Branch("nu_vtx_x", &nu_vtx_x); 
  fEventTree->Branch("nu_vtx_y", &nu_vtx_y); 
  fEventTree->Branch("nu_vtx_z", &nu_vtx_z); 
  fEventTree->Branch("nu_vtx_t", &nu_vtx_t); 
  fEventTree->Branch("nu_n_protons", &nu_n_protons); 
  fEventTree->Branch("nu_n_neutrons", &nu_n_neutrons); 
  fEventTree->Branch("nu_n_charged_pions", &nu_n_charged_pions); 
  fEventTree->Branch("nu_n_neutral_pions", &nu_n_neutral_pions); 
  fEventTree->Branch("nu_n_photons", &nu_n_photons); 
  fEventTree->Branch("nu_n_others", &nu_n_others); 
  
  //Event Tree: Slice
  fEventTree->Branch("slc_id", &slc_id);
  fEventTree->Branch("n_slc", &_n_slc);
  fEventTree->Branch("slc_n_pfps", &slc_n_pfps);
  fEventTree->Branch("slc_is_clear_cosmics", &slc_is_clear_cosmics);
  fEventTree->Branch("slc_primary_pfp_id", &slc_primary_pfp_id);
  fEventTree->Branch("slc_primary_pfp_pdg", &slc_primary_pfp_pdg);
  fEventTree->Branch("slc_n_primary_daughters", &slc_n_primary_daughters);
  fEventTree->Branch("slc_vtx_x", &slc_vtx_x);
  fEventTree->Branch("slc_vtx_y", &slc_vtx_y);
  fEventTree->Branch("slc_vtx_z", &slc_vtx_z);
  fEventTree->Branch("slc_is_fv", &slc_is_fv);
  fEventTree->Branch("slc_crumbs_score", &slc_crumbs_score);
  fEventTree->Branch("slc_crumbs_nc_score", &slc_crumbs_nc_score);
  fEventTree->Branch("slc_crumbs_ccnue_score", &slc_crumbs_ccnue_score);
  fEventTree->Branch("slc_opt0_time", &slc_opt0_time);
  fEventTree->Branch("slc_opt0_time_corrected_Z_pandora", &slc_opt0_time_corrected_Z_pandora);
  fEventTree->Branch("slc_opt0_score", &slc_opt0_score);
  fEventTree->Branch("slc_opt0_measPE", &slc_opt0_measPE);
  fEventTree->Branch("slc_opt0_hypoPE", &slc_opt0_hypoPE);
  fEventTree->Branch("slc_opt0_frac", &slc_opt0_frac);
  fEventTree->Branch("slc_n_trks", &slc_n_trks);
  fEventTree->Branch("slc_n_shws", &slc_n_shws);
  fEventTree->Branch("slc_n_primary_trks", &slc_n_primary_trks);
  fEventTree->Branch("slc_n_primary_shws", &slc_n_primary_shws);
  fEventTree->Branch("slc_n_dazzle_muons", &slc_n_dazzle_muons);
  fEventTree->Branch("slc_n_dazzle_pions", &slc_n_dazzle_pions);
  fEventTree->Branch("slc_n_dazzle_pions_thresh", &slc_n_dazzle_pions_thresh);
  fEventTree->Branch("slc_n_dazzle_protons", &slc_n_dazzle_protons);
  fEventTree->Branch("slc_n_dazzle_protons_thresh", &slc_n_dazzle_protons_thresh);
  fEventTree->Branch("slc_n_dazzle_other", &slc_n_dazzle_other);
  fEventTree->Branch("slc_n_primary_dazzle_muons", &slc_n_primary_dazzle_muons);
  fEventTree->Branch("slc_n_primary_dazzle_pions", &slc_n_primary_dazzle_pions);
  fEventTree->Branch("slc_n_primary_dazzle_pions_thresh", &slc_n_primary_dazzle_pions_thresh);
  fEventTree->Branch("slc_n_primary_dazzle_protons", &slc_n_primary_dazzle_protons);
  fEventTree->Branch("slc_n_primary_dazzle_protons_thresh", &slc_n_primary_dazzle_protons_thresh);
  fEventTree->Branch("slc_n_primary_dazzle_other", &slc_n_primary_dazzle_other);
  fEventTree->Branch("slc_n_razzle_electrons", &slc_n_razzle_electrons);
  fEventTree->Branch("slc_n_razzle_photons", &slc_n_razzle_photons);
  fEventTree->Branch("slc_n_razzle_other", &slc_n_razzle_other);
  fEventTree->Branch("slc_n_primary_razzle_electrons", &slc_n_primary_razzle_electrons);
  fEventTree->Branch("slc_n_primary_razzle_photons", &slc_n_primary_razzle_photons);
  fEventTree->Branch("slc_n_primary_razzle_other", &slc_n_primary_razzle_other);
  fEventTree->Branch("slc_n_razzled_electrons", &slc_n_razzled_electrons);
  fEventTree->Branch("slc_n_razzled_muons", &slc_n_razzled_muons);
  fEventTree->Branch("slc_n_razzled_photons", &slc_n_razzled_photons);
  fEventTree->Branch("slc_n_razzled_pions", &slc_n_razzled_pions);
  fEventTree->Branch("slc_n_razzled_pions_thresh", &slc_n_razzled_pions_thresh);
  fEventTree->Branch("slc_n_razzled_protons", &slc_n_razzled_protons);
  fEventTree->Branch("slc_n_razzled_protons_thresh", &slc_n_razzled_protons_thresh);
  fEventTree->Branch("slc_n_primary_razzled_electrons", &slc_n_primary_razzled_electrons);
  fEventTree->Branch("slc_n_primary_razzled_muons", &slc_n_primary_razzled_muons);
  fEventTree->Branch("slc_n_primary_razzled_photons", &slc_n_primary_razzled_photons);
  fEventTree->Branch("slc_n_primary_razzled_pions", &slc_n_primary_razzled_pions);
  fEventTree->Branch("slc_n_primary_razzled_pions_thresh", &slc_n_primary_razzled_pions_thresh);
  fEventTree->Branch("slc_n_primary_razzled_protons", &slc_n_primary_razzled_protons);
  fEventTree->Branch("slc_n_primary_razzled_protons_thresh", &slc_n_primary_razzled_protons_thresh);
  fEventTree->Branch("slc_n_stub", &slc_n_stub);
  fEventTree->Branch("slc_total_track_E", &slc_total_track_E); 
  fEventTree->Branch("slc_total_shower_E", &slc_total_shower_E); 
  
  //Event Tree: Slice Truth
  fEventTree->Branch("slc_comp" , &slc_comp);
  fEventTree->Branch("slc_pur" , &slc_pur);
  fEventTree->Branch("slc_true_mctruth_id", &slc_true_mctruth_id);
  fEventTree->Branch("slc_true_event_type", &slc_true_event_type);
  fEventTree->Branch("slc_true_en_dep", &slc_true_en_dep);
  fEventTree->Branch("slc_true_vtx_x", &slc_true_vtx_x);
  fEventTree->Branch("slc_true_vtx_y", &slc_true_vtx_y);
  fEventTree->Branch("slc_true_vtx_z", &slc_true_vtx_z);
  fEventTree->Branch("slc_true_vtx_t", &slc_true_vtx_t);
  fEventTree->Branch("slc_true_vtx_t_corrected_Z", &slc_true_vtx_t_corrected_Z);

  // Event Tree: Slice -> PFP
  fEventTree->Branch("slc_pfp_id", &slc_pfp_id);
  fEventTree->Branch("slc_pfp_pdg", &slc_pfp_pdg);
  fEventTree->Branch("slc_pfp_n_children", &slc_pfp_n_children);
  fEventTree->Branch("slc_pfp_primary", &slc_pfp_primary);
  fEventTree->Branch("slc_pfp_primary_child", &slc_pfp_primary_child);
  fEventTree->Branch("slc_pfp_n_hits", &slc_pfp_n_hits);
  fEventTree->Branch("slc_pfp_n_sps", &slc_pfp_n_sps);
  fEventTree->Branch("slc_pfp_track_score", &slc_pfp_track_score);
  fEventTree->Branch("slc_pfp_good_track", &slc_pfp_good_track);
  fEventTree->Branch("slc_pfp_good_shower", &slc_pfp_good_shower);
  fEventTree->Branch("slc_pfp_true_trackid", &slc_pfp_true_trackid);
  fEventTree->Branch("slc_pfp_comp", &slc_pfp_comp);
  fEventTree->Branch("slc_pfp_pur", &slc_pfp_pur);
//  fEventTree->Branch("slc_pfp_cnnscore_track", &slc_pfp_cnnscore_track);
//  fEventTree->Branch("slc_pfp_cnnscore_shower", &slc_pfp_cnnscore_shower);
//  fEventTree->Branch("slc_pfp_cnnscore_noise", &slc_pfp_cnnscore_noise);
//  fEventTree->Branch("slc_pfp_cnnscore_michel", &slc_pfp_cnnscore_michel);
//  fEventTree->Branch("slc_pfp_cnnscore_endmichel", &slc_pfp_cnnscore_endmichel);
//  fEventTree->Branch("slc_pfp_cnnscore_nclusters", &slc_pfp_cnnscore_nclusters);
  fEventTree->Branch("slc_pfp_razzled_electron_score", &slc_pfp_razzled_electron_score);
  fEventTree->Branch("slc_pfp_razzled_muon_score", &slc_pfp_razzled_muon_score);
  fEventTree->Branch("slc_pfp_razzled_photon_score", &slc_pfp_razzled_photon_score);
  fEventTree->Branch("slc_pfp_razzled_pion_score", &slc_pfp_razzled_pion_score);
  fEventTree->Branch("slc_pfp_razzled_proton_score", &slc_pfp_razzled_proton_score);
  fEventTree->Branch("slc_pfp_razzled_pdg", &slc_pfp_razzled_pdg);
  
  // Event Tree: Slice -> PFP -> MCParticle
  fEventTree->Branch("slc_pfp_true_trackid", &slc_pfp_true_trackid);
  fEventTree->Branch("slc_pfp_true_pdg", &slc_pfp_true_pdg);
  fEventTree->Branch("slc_pfp_true_energy", &slc_pfp_true_energy);
  fEventTree->Branch("slc_pfp_true_p_x", &slc_pfp_true_p_x);
  fEventTree->Branch("slc_pfp_true_p_y", &slc_pfp_true_p_y);
  fEventTree->Branch("slc_pfp_true_p_z", &slc_pfp_true_p_z);

  // Event Tree: Slice -> PFP -> Track
  fEventTree->Branch("slc_pfp_track_start_x", &slc_pfp_track_start_x);
  fEventTree->Branch("slc_pfp_track_start_y", &slc_pfp_track_start_y);
  fEventTree->Branch("slc_pfp_track_start_z", &slc_pfp_track_start_z);
  fEventTree->Branch("slc_pfp_track_end_x", &slc_pfp_track_end_x);
  fEventTree->Branch("slc_pfp_track_end_y", &slc_pfp_track_end_y);
  fEventTree->Branch("slc_pfp_track_end_z", &slc_pfp_track_end_z);
  fEventTree->Branch("slc_pfp_track_dir_x", &slc_pfp_track_dir_x);
  fEventTree->Branch("slc_pfp_track_dir_y", &slc_pfp_track_dir_y);
  fEventTree->Branch("slc_pfp_track_dir_z", &slc_pfp_track_dir_z);
  fEventTree->Branch("slc_pfp_track_length", &slc_pfp_track_length);
  fEventTree->Branch("slc_pfp_track_contained", &slc_pfp_track_contained);
  fEventTree->Branch("slc_pfp_track_dazzle_muon_score", &slc_pfp_track_dazzle_muon_score);
  fEventTree->Branch("slc_pfp_track_dazzle_pion_score", &slc_pfp_track_dazzle_pion_score);
  fEventTree->Branch("slc_pfp_track_dazzle_proton_score", &slc_pfp_track_dazzle_proton_score);
  fEventTree->Branch("slc_pfp_track_dazzle_other_score", &slc_pfp_track_dazzle_other_score);
  fEventTree->Branch("slc_pfp_track_dazzle_pdg", &slc_pfp_track_dazzle_pdg);
  fEventTree->Branch("slc_pfp_track_ke", &slc_pfp_track_ke);
  fEventTree->Branch("slc_pfp_track_charge", &slc_pfp_track_charge);
  fEventTree->Branch("slc_pfp_track_chi2_muon", &slc_pfp_track_chi2_muon);
  fEventTree->Branch("slc_pfp_track_chi2_pion", &slc_pfp_track_chi2_pion);
  fEventTree->Branch("slc_pfp_track_chi2_kaon", &slc_pfp_track_chi2_kaon);
  fEventTree->Branch("slc_pfp_track_chi2_proton", &slc_pfp_track_chi2_proton);
  fEventTree->Branch("slc_pfp_track_chi2_pdg", &slc_pfp_track_chi2_pdg);
  fEventTree->Branch("slc_pfp_track_mcs_mean_scatter", &slc_pfp_track_mcs_mean_scatter);
  fEventTree->Branch("slc_pfp_track_mcs_max_scatter_ratio", &slc_pfp_track_mcs_max_scatter_ratio);
  fEventTree->Branch("slc_pfp_track_range_p", &slc_pfp_track_range_p);
  fEventTree->Branch("slc_pfp_track_closest_approach_mean_dca", &slc_pfp_track_closest_approach_mean_dca);
  fEventTree->Branch("slc_pfp_track_stopping_dedx_chi2_ratio", &slc_pfp_track_stopping_dedx_chi2_ratio);
  fEventTree->Branch("slc_pfp_track_stopping_dedx_pol0_fit", &slc_pfp_track_stopping_dedx_pol0_fit);
  fEventTree->Branch("slc_pfp_track_theta", &slc_pfp_track_theta); 
  fEventTree->Branch("slc_pfp_track_phi", &slc_pfp_track_phi); 

  // Event Tree: Slice -> PFP -> Shower -- 2D vector
  fEventTree->Branch("slc_pfp_shower_start_x", &slc_pfp_shower_start_x);
  fEventTree->Branch("slc_pfp_shower_start_y", &slc_pfp_shower_start_y);
  fEventTree->Branch("slc_pfp_shower_start_z", &slc_pfp_shower_start_z);
  fEventTree->Branch("slc_pfp_shower_end_x", &slc_pfp_shower_end_x);
  fEventTree->Branch("slc_pfp_shower_end_y", &slc_pfp_shower_end_y);
  fEventTree->Branch("slc_pfp_shower_end_z", &slc_pfp_shower_end_z);
  fEventTree->Branch("slc_pfp_shower_conv_gap", &slc_pfp_shower_conv_gap);
  fEventTree->Branch("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
  fEventTree->Branch("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
  fEventTree->Branch("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
  fEventTree->Branch("slc_pfp_shower_length", &slc_pfp_shower_length);
  fEventTree->Branch("slc_pfp_shower_open_angle", &slc_pfp_shower_open_angle);
  fEventTree->Branch("slc_pfp_shower_contained", &slc_pfp_shower_contained);
  fEventTree->Branch("slc_pfp_shower_energy", &slc_pfp_shower_energy);
  fEventTree->Branch("slc_pfp_shower_dedx", &slc_pfp_shower_dedx);
  fEventTree->Branch("slc_pfp_shower_sqrt_energy_density", &slc_pfp_shower_sqrt_energy_density);
  fEventTree->Branch("slc_pfp_shower_modified_hit_density", &slc_pfp_shower_modified_hit_density);
  fEventTree->Branch("slc_pfp_shower_razzle_electron_score", &slc_pfp_shower_razzle_electron_score);
  fEventTree->Branch("slc_pfp_shower_razzle_photon_score", &slc_pfp_shower_razzle_photon_score);
  fEventTree->Branch("slc_pfp_shower_razzle_other_score", &slc_pfp_shower_razzle_other_score);
  fEventTree->Branch("slc_pfp_shower_razzle_pdg", &slc_pfp_shower_razzle_pdg);
  fEventTree->Branch("slc_pfp_shower_cosmic_dist", &slc_pfp_shower_cosmic_dist);
  fEventTree->Branch("slc_pfp_shower_track_length", &slc_pfp_shower_track_length);
  fEventTree->Branch("slc_pfp_shower_track_width", &slc_pfp_shower_track_width);
  fEventTree->Branch("slc_pfp_shower_density_grad", &slc_pfp_shower_density_grad);
  fEventTree->Branch("slc_pfp_shower_density_pow", &slc_pfp_shower_density_pow);
  fEventTree->Branch("slc_pfp_shower_theta", &slc_pfp_shower_theta); 
  fEventTree->Branch("slc_pfp_shower_phi", &slc_pfp_shower_phi); 

  // Event Tree: OpFlashes
  fEventTree->Branch("flash_time","std::vector<double>", &_flash_time);
  fEventTree->Branch("flash_total_pe", "std::vector<double>", &_flash_total_pe);
  fEventTree->Branch("flash_pe_v","std::vector<std::vector<double>>", &_flash_pe_v);
  fEventTree->Branch("flash_tpc", "std::vector<int>", &_flash_tpc);
  fEventTree->Branch("flash_y","std::vector<double>", &_flash_y);
  fEventTree->Branch("flash_yerr", "std::vector<double>", &_flash_yerr);
  fEventTree->Branch("flash_z","std::vector<double>", &_flash_z);
  fEventTree->Branch("flash_zerr", "std::vector<double>", &_flash_zerr);

  // Event Tree: Flux Weight
  fEventTree->Branch("slc_flux_weight_expskin", &slc_flux_weight_expskin);
  fEventTree->Branch("slc_flux_weight_horncurrent", &slc_flux_weight_horncurrent);
  fEventTree->Branch("slc_flux_weight_kminus", &slc_flux_weight_kminus);
  fEventTree->Branch("slc_flux_weight_kplus", &slc_flux_weight_kplus);
  fEventTree->Branch("slc_flux_weight_kzero", &slc_flux_weight_kzero);
  fEventTree->Branch("slc_flux_weight_nucleoninexsec", &slc_flux_weight_nucleoninexsec);
  fEventTree->Branch("slc_flux_weight_nucleonqexsec", &slc_flux_weight_nucleonqexsec);
  fEventTree->Branch("slc_flux_weight_nucleontotxsec", &slc_flux_weight_nucleontotxsec);
  fEventTree->Branch("slc_flux_weight_piminus", &slc_flux_weight_piminus);
  fEventTree->Branch("slc_flux_weight_pioninexsex", &slc_flux_weight_pioninexsex);
  fEventTree->Branch("slc_flux_weight_pionqexsec", &slc_flux_weight_pionqexsec);
  fEventTree->Branch("slc_flux_weight_piontotxsec", &slc_flux_weight_piontotxsec);
  fEventTree->Branch("slc_flux_weight_piplus", &slc_flux_weight_piplus);
  fEventTree->Branch("slc_flux_weight_total", &slc_flux_weight_total);

  // Event Tree: XSection Weight
  fEventTree->Branch("slc_xsec_unisim_DecayAngMEC", &slc_xsec_unisim_DecayAngMEC);
  fEventTree->Branch("slc_xsec_unisim_ThetaDelta2NRad", &slc_xsec_unisim_ThetaDelta2NRad); 
  fEventTree->Branch("slc_xsec_unisim_Theta_Delta2Npi", &slc_xsec_unisim_Theta_Delta2Npi);
  fEventTree->Branch("slc_xsec_unisim_VecFFCCQEshape", &slc_xsec_unisim_VecFFCCQEshape);

  fEventTree->Branch("slc_xsec_multisigma_CoulombCCQE", &slc_xsec_multisigma_CoulombCCQE);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarnCC1pi", &slc_xsec_multisigma_NonRESBGvbarnCC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarnCC2pi", &slc_xsec_multisigma_NonRESBGvbarnCC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarnNC1pi", &slc_xsec_multisigma_NonRESBGvbarnNC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarnNC2pi", &slc_xsec_multisigma_NonRESBGvbarnNC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarpCC1pi", &slc_xsec_multisigma_NonRESBGvbarpCC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarpCC2pi", &slc_xsec_multisigma_NonRESBGvbarpCC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarpNC1pi", &slc_xsec_multisigma_NonRESBGvbarpNC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvbarpNC2pi", &slc_xsec_multisigma_NonRESBGvbarpNC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvnCC1pi", &slc_xsec_multisigma_NonRESBGvnCC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvnCC2pi", &slc_xsec_multisigma_NonRESBGvnCC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvnNC1pi", &slc_xsec_multisigma_NonRESBGvnNC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvnNC2pi", &slc_xsec_multisigma_NonRESBGvnNC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvpCC1pi", &slc_xsec_multisigma_NonRESBGvpCC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvpCC2pi", &slc_xsec_multisigma_NonRESBGvpCC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvpNC1pi", &slc_xsec_multisigma_NonRESBGvpNC1pi);
  fEventTree->Branch("slc_xsec_multisigma_NonRESBGvpNC2pi", &slc_xsec_multisigma_NonRESBGvpNC2pi);
  fEventTree->Branch("slc_xsec_multisigma_NormCCMEC", &slc_xsec_multisigma_NormCCMEC);
  fEventTree->Branch("slc_xsec_multisigma_NormNCMEC", &slc_xsec_multisigma_NormNCMEC);
  fEventTree->Branch("slc_xsec_multisigma_RDecBR1eta", &slc_xsec_multisigma_RDecBR1eta);
  fEventTree->Branch("slc_xsec_multisigma_RDecBR1gamma", &slc_xsec_multisigma_RDecBR1gamma);
  fEventTree->Branch("slc_xsec_multisigma_RPA_CCQE", &slc_xsec_multisigma_RPA_CCQE);
  fEventTree->Branch("slc_xsec_multisigma_NormNCCOH", &slc_xsec_multisigma_NormNCCOH);
  fEventTree->Branch("slc_xsec_multisigma_NormCCCOH", &slc_xsec_multisigma_NormCCCOH);
  fEventTree->Branch("slc_xsec_multisim_ZExpA", &slc_xsec_multisim_ZExpA);
  fEventTree->Branch("slc_xsec_multisim_NCEL", &slc_xsec_multisim_NCEL);
  fEventTree->Branch("slc_xsec_multisim_CCRES", &slc_xsec_multisim_CCRES);
  fEventTree->Branch("slc_xsec_multisim_NCRES", &slc_xsec_multisim_NCRES);
  fEventTree->Branch("slc_xsec_multisim_DISBY", &slc_xsec_multisim_DISBY);
  fEventTree->Branch("slc_xsec_multisim_FSI_pi", &slc_xsec_multisim_FSI_pi);
  fEventTree->Branch("slc_xsec_multisim_FSI_N", &slc_xsec_multisim_FSI_N);
  fEventTree->Branch("slc_xsec_multisim_total", &slc_xsec_multisim_total);
  
  fEventTree->Branch("slc_geant4_multisim_reinteractions", &slc_geant4_multisim_reinteractions);

  if (fPiZeroAna){
    fEventTree->Branch("slc_pizero_good", &slc_pizero_good);
    fEventTree->Branch("slc_photon0_pfp_id", &slc_photon0_pfp_id);
    fEventTree->Branch("slc_photon1_pfp_id", &slc_photon1_pfp_id);
    fEventTree->Branch("slc_photon2_pfp_id", &slc_photon2_pfp_id);
    fEventTree->Branch("slc_photon0_energy", &slc_photon0_energy);
    fEventTree->Branch("slc_photon1_energy", &slc_photon1_energy);
    fEventTree->Branch("slc_photon2_energy", &slc_photon2_energy);
    fEventTree->Branch("slc_pizero_open_angle", &slc_pizero_open_angle);
    fEventTree->Branch("slc_pizero_mass", &slc_pizero_mass);
    fEventTree->Branch("slc_pizero_mom", &slc_pizero_mom);
    fEventTree->Branch("slc_pizero_CosThetaZ", &slc_pizero_CosThetaZ);
   
    fEventTree->Branch("slc_true_pizero_mom", &slc_true_pizero_mom);
    fEventTree->Branch("slc_true_pizero_cosThetaZ", &slc_true_pizero_cosThetaZ);

  }

}

void sbnd::HNLPiZeroAnalysis::analyze(art::Event const& e)
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
    std::cout << std::endl << ">>>> This is event " << _run << "-" << _subrun << "-" << _event << std::endl;
  
  // Get MeVPrtl
  art::Handle<std::vector<evgen::ldm::MeVPrtlTruth>> mevptHandle;
  e.getByLabel(fMeVPrtlTruthModuleLabel, mevptHandle);

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

  // Get OpFlashes
  art::Handle<std::vector<recob::OpFlash>> opflashListHandle; // 2 opflashes objects, one for each TPC

  //if (fDebug) std::cout<<"Saving OpFlashes..."<<std::endl;
  _flash_id.clear();
  _flash_time.clear();
  _flash_total_pe.clear();
  _flash_pe_v.clear();
  _flash_tpc.clear();
  _flash_y.clear();
  _flash_yerr.clear();
  _flash_z.clear();
  _flash_zerr.clear();

  for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) 
  {
    //if (fDebug) std::cout<<"  --OpFlash Module Label:"<<fOpFlashesModuleLabel[s]<<std::endl;
    e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);
    if(!opflashListHandle.isValid())
    {
      //if (fDebug) std::cout<<"Did not find any OpFlash for label "<<fOpFlashesModuleLabel[s]<<std::endl;
    }
    else
    {
      for (unsigned int i = 0; i < opflashListHandle->size(); ++i) 
      {
        // Get OpFlash
        art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
        recob::OpFlash Flash = *FlashPtr;
        _flash_time.push_back( Flash.AbsTime() );
        _flash_total_pe.push_back( Flash.TotalPE() );
        _flash_pe_v.push_back( Flash.PEs() );
        _flash_tpc.push_back( s );
        _flash_y.push_back( Flash.YCenter() );
        _flash_yerr.push_back( Flash.YWidth() );
        _flash_z.push_back( Flash.ZCenter() );
        _flash_zerr.push_back( Flash.ZWidth() );
      }
    }
  }

  SetupMaps(e, hitHandle, pfpHandle);
  AnalyseMeVPrtlTruth(e, mevptHandle);
  AnalyseMCTruthHandle(e, MCTruthHandles);
  AnalyseSlices(e, sliceHandle, pfpHandle, trackHandle, showerHandle);

  // Fill Tree
  fEventTree->Fill();

}

void sbnd::HNLPiZeroAnalysis::beginSubRun(const art::SubRun &sr)
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

  if(fDebug){
    std::cout << "---POT this subrun is " << _pot << std::endl;
    std::cout << "---SPILL this subrun is " << _spills << std::endl;
  }
}

void sbnd::HNLPiZeroAnalysis::endSubRun(const art::SubRun &sr)
{
  fSubRunTree->Fill();
}

void sbnd::HNLPiZeroAnalysis::ResetSubRunVars()
{
  _pot = 0.; _spills = 0; _ngenevts = 0;
}

void sbnd::HNLPiZeroAnalysis::ResetEventVars()
{
  _run = -1; _subrun = -1; _event = -1;

  _n_hnl = 0; _n_mctruth = 0; _n_nu = 0; _n_slc = 0;

  mevprtl_decay_pos_x.clear(); mevprtl_decay_pos_y.clear(); mevprtl_decay_pos_z.clear(); mevprtl_decay_pos_t.clear();
  mevprtl_mom_x.clear(); mevprtl_mom_y.clear(); mevprtl_mom_z.clear(); mevprtl_mom_e.clear();
  mevprtl_mass.clear();
  mevprtl_flux_weight.clear(); mevprtl_ray_weight.clear(); mevprtl_decay_weight.clear();
  mevprtl_total_decay_width.clear(); mevprtl_total_mean_lifetime.clear(); mevprtl_total_mean_distance.clear(); mevprtl_allowed_decay_fraction.clear();
  mevprtl_C1.clear(); mevprtl_C2.clear(); mevprtl_C3.clear(); mevprtl_C4.clear(); mevprtl_C5.clear();
  
  nu_mctruth_id.clear();
  nu_event_type.clear(); 
  nu_signal.clear(); 
  nu_en_dep.clear();
  nu_ccnc.clear();
  nu_mode.clear();
  nu_int_type.clear();
  nu_w.clear();
  nu_x.clear();
  nu_y.clear();
  nu_q_sqr.clear();
  nu_pt.clear();
  nu_theta.clear();
  nu_e.clear();
  nu_vtx_x.clear(); nu_vtx_y.clear(); nu_vtx_z.clear(); nu_vtx_t.clear();
  nu_n_protons.clear();
  nu_n_neutrons.clear();
  nu_n_charged_pions.clear();
  nu_n_neutral_pions.clear();
  nu_n_photons.clear();
  nu_n_others.clear();

  slc_id.clear();
  slc_n_pfps.clear();
  slc_is_clear_cosmics.clear();
  slc_primary_pfp_id.clear();
  slc_primary_pfp_pdg.clear();
  slc_n_primary_daughters.clear();
  slc_vtx_x.clear();
  slc_vtx_y.clear();
  slc_vtx_z.clear();
  slc_is_fv.clear();
  slc_crumbs_score.clear();
  slc_crumbs_nc_score.clear();
  slc_crumbs_ccnue_score.clear();
  slc_opt0_time.clear();
  slc_opt0_time_corrected_Z_pandora.clear();
  slc_opt0_score.clear();
  slc_opt0_measPE.clear();
  slc_opt0_hypoPE.clear();
  slc_opt0_frac.clear();
  slc_n_trks.clear();
  slc_n_shws.clear();
  slc_n_primary_trks.clear(); 
  slc_n_primary_shws.clear(); 
  slc_n_dazzle_muons.clear(); 
  slc_n_dazzle_pions.clear(); 
  slc_n_dazzle_pions_thresh.clear();
  slc_n_dazzle_protons.clear();
  slc_n_dazzle_protons_thresh.clear();
  slc_n_dazzle_other.clear();
  slc_n_primary_dazzle_muons.clear();
  slc_n_primary_dazzle_pions.clear();
  slc_n_primary_dazzle_pions_thresh.clear(); 
  slc_n_primary_dazzle_protons.clear();
  slc_n_primary_dazzle_protons_thresh.clear(); 
  slc_n_primary_dazzle_other.clear();
  slc_n_razzle_electrons.clear();
  slc_n_razzle_photons.clear();
  slc_n_razzle_other.clear();
  slc_n_primary_razzle_electrons.clear(); 
  slc_n_primary_razzle_photons.clear();
  slc_n_primary_razzle_other.clear();
  slc_n_razzled_electrons.clear();
  slc_n_razzled_muons.clear();
  slc_n_razzled_photons.clear();
  slc_n_razzled_pions.clear();
  slc_n_razzled_pions_thresh.clear();
  slc_n_razzled_protons.clear();
  slc_n_razzled_protons_thresh.clear();
  slc_n_primary_razzled_electrons.clear(); 
  slc_n_primary_razzled_muons.clear();
  slc_n_primary_razzled_photons.clear();
  slc_n_primary_razzled_pions.clear();
  slc_n_primary_razzled_pions_thresh.clear(); 
  slc_n_primary_razzled_protons.clear();
  slc_n_primary_razzled_protons_thresh.clear(); 
  slc_n_stub.clear();
  slc_total_track_E.clear();
  slc_total_shower_E.clear();

  slc_comp.clear();
  slc_pur.clear();
  slc_true_mctruth_id.clear();
  slc_true_event_type.clear();
  slc_true_en_dep.clear();
  slc_true_vtx_x.clear();
  slc_true_vtx_y.clear();
  slc_true_vtx_z.clear();
  slc_true_vtx_t.clear();
  slc_true_vtx_t_corrected_Z.clear();

  slc_pfp_id.clear();
  slc_pfp_pdg.clear();
  slc_pfp_n_children.clear();
  slc_pfp_primary.clear();
  slc_pfp_primary_child.clear();
  slc_pfp_n_hits.clear();
  slc_pfp_n_sps.clear();
  slc_pfp_track_score.clear();
  slc_pfp_good_track.clear();
  slc_pfp_good_shower.clear();
  slc_pfp_comp.clear();
  slc_pfp_pur.clear();
  slc_pfp_cnnscore_track.clear();
  slc_pfp_cnnscore_shower.clear();
  slc_pfp_cnnscore_noise.clear();
  slc_pfp_cnnscore_michel.clear();
  slc_pfp_cnnscore_endmichel.clear();
  slc_pfp_cnnscore_nclusters.clear();
  slc_pfp_razzled_electron_score.clear();
  slc_pfp_razzled_muon_score.clear();
  slc_pfp_razzled_photon_score.clear();
  slc_pfp_razzled_pion_score.clear();
  slc_pfp_razzled_proton_score.clear();
  slc_pfp_razzled_pdg.clear();

  slc_pfp_true_trackid.clear();
  slc_pfp_true_pdg.clear();
  slc_pfp_true_energy.clear();
  slc_pfp_true_p_x.clear();
  slc_pfp_true_p_y.clear();
  slc_pfp_true_p_z.clear();
  
  slc_pfp_track_start_x.clear();
  slc_pfp_track_start_y.clear();
  slc_pfp_track_start_z.clear();
  slc_pfp_track_end_x.clear();
  slc_pfp_track_end_y.clear();
  slc_pfp_track_end_z.clear();
  slc_pfp_track_dir_x.clear();
  slc_pfp_track_dir_y.clear();
  slc_pfp_track_dir_z.clear();
  slc_pfp_track_length.clear();
  slc_pfp_track_contained.clear();
  slc_pfp_track_dazzle_muon_score.clear();
  slc_pfp_track_dazzle_pion_score.clear();
  slc_pfp_track_dazzle_proton_score.clear();
  slc_pfp_track_dazzle_other_score.clear();
  slc_pfp_track_dazzle_pdg.clear();
  slc_pfp_track_ke.clear();
  slc_pfp_track_charge.clear();
  slc_pfp_track_chi2_muon.clear();
  slc_pfp_track_chi2_pion.clear();
  slc_pfp_track_chi2_kaon.clear();
  slc_pfp_track_chi2_muon.clear();
  slc_pfp_track_chi2_pdg.clear();
  slc_pfp_track_mcs_mean_scatter.clear();
  slc_pfp_track_mcs_max_scatter_ratio.clear();
  slc_pfp_track_range_p.clear();
  slc_pfp_track_closest_approach_mean_dca.clear();
  slc_pfp_track_stopping_dedx_chi2_ratio.clear();
  slc_pfp_track_stopping_dedx_pol0_fit.clear();
  slc_pfp_track_theta.clear();
  slc_pfp_track_phi.clear();

  slc_pfp_shower_start_x.clear();
  slc_pfp_shower_start_y.clear();
  slc_pfp_shower_start_z.clear();
  slc_pfp_shower_end_x.clear();
  slc_pfp_shower_end_y.clear();
  slc_pfp_shower_end_z.clear();
  slc_pfp_shower_conv_gap.clear();
  slc_pfp_shower_dir_x.clear();
  slc_pfp_shower_dir_y.clear();
  slc_pfp_shower_dir_z.clear();
  slc_pfp_shower_length.clear();
  slc_pfp_shower_open_angle.clear();
  slc_pfp_shower_energy.clear();
  slc_pfp_shower_dedx.clear();
  slc_pfp_shower_sqrt_energy_density.clear();
  slc_pfp_shower_modified_hit_density.clear();
  slc_pfp_shower_razzle_electron_score.clear();
  slc_pfp_shower_razzle_photon_score.clear();
  slc_pfp_shower_razzle_other_score.clear();
  slc_pfp_shower_razzle_pdg.clear();
  slc_pfp_shower_cosmic_dist.clear();
  slc_pfp_shower_track_length.clear();
  slc_pfp_shower_track_width.clear();
  slc_pfp_shower_density_grad.clear();
  slc_pfp_shower_density_pow.clear();
  slc_pfp_shower_contained.clear();
  slc_pfp_shower_theta.clear();
  slc_pfp_shower_phi.clear();

  slc_flux_weight_expskin.clear();
  slc_flux_weight_horncurrent.clear();
  slc_flux_weight_kminus.clear();
  slc_flux_weight_kplus.clear();
  slc_flux_weight_kzero.clear();
  slc_flux_weight_nucleoninexsec.clear();
  slc_flux_weight_nucleonqexsec.clear();
  slc_flux_weight_nucleontotxsec.clear();
  slc_flux_weight_piminus.clear();
  slc_flux_weight_pioninexsex.clear();
  slc_flux_weight_pionqexsec.clear();
  slc_flux_weight_piontotxsec.clear();
  slc_flux_weight_piplus.clear();
  slc_flux_weight_total.clear();

  slc_xsec_unisim_DecayAngMEC.clear();
  slc_xsec_unisim_ThetaDelta2NRad.clear(); 
  slc_xsec_unisim_Theta_Delta2Npi.clear();
  slc_xsec_unisim_VecFFCCQEshape.clear();
  
  slc_xsec_multisigma_CoulombCCQE.clear();
  slc_xsec_multisigma_NonRESBGvbarnCC1pi.clear();
  slc_xsec_multisigma_NonRESBGvbarnCC2pi.clear();
  slc_xsec_multisigma_NonRESBGvbarnNC1pi.clear();
  slc_xsec_multisigma_NonRESBGvbarnNC2pi.clear();
  slc_xsec_multisigma_NonRESBGvbarpCC1pi.clear();
  slc_xsec_multisigma_NonRESBGvbarpCC2pi.clear();
  slc_xsec_multisigma_NonRESBGvbarpNC1pi.clear();
  slc_xsec_multisigma_NonRESBGvbarpNC2pi.clear();
  slc_xsec_multisigma_NonRESBGvnCC1pi.clear();
  slc_xsec_multisigma_NonRESBGvnCC2pi.clear();
  slc_xsec_multisigma_NonRESBGvnNC1pi.clear();
  slc_xsec_multisigma_NonRESBGvnNC2pi.clear();
  slc_xsec_multisigma_NonRESBGvpCC1pi.clear();
  slc_xsec_multisigma_NonRESBGvpCC2pi.clear();
  slc_xsec_multisigma_NonRESBGvpNC1pi.clear();
  slc_xsec_multisigma_NonRESBGvpNC2pi.clear();
  slc_xsec_multisigma_NormCCMEC.clear();
  slc_xsec_multisigma_NormNCMEC.clear();
  slc_xsec_multisigma_RDecBR1eta.clear();
  slc_xsec_multisigma_RDecBR1gamma.clear();
  slc_xsec_multisigma_RPA_CCQE.clear();
  slc_xsec_multisigma_NormNCCOH.clear();
  slc_xsec_multisigma_NormCCCOH.clear();

  slc_xsec_multisim_ZExpA.clear();
  slc_xsec_multisim_NCEL.clear();
  slc_xsec_multisim_CCRES.clear();
  slc_xsec_multisim_NCRES.clear();
  slc_xsec_multisim_DISBY.clear();
  slc_xsec_multisim_FSI_pi.clear();
  slc_xsec_multisim_FSI_N.clear();
  slc_xsec_multisim_total.clear();
  
  slc_geant4_multisim_reinteractions.clear();

  if (fPiZeroAna){
    slc_pizero_good.clear();
    slc_photon0_pfp_id.clear();
    slc_photon1_pfp_id.clear();
    slc_photon2_pfp_id.clear();
    slc_photon0_energy.clear();
    slc_photon1_energy.clear();
    slc_photon2_energy.clear();
    slc_pizero_open_angle.clear();
    slc_pizero_mass.clear();
    slc_pizero_mom.clear();
    slc_pizero_CosThetaZ.clear();

    slc_true_pizero_mom.clear();
    slc_true_pizero_cosThetaZ.clear();
  }
}

void sbnd::HNLPiZeroAnalysis::ClearMaps()
{
  fHitsMap.clear();
  fMCTruthHitsMap.clear();
  fPFPMap.clear();
  fRecoPFPMap.clear();
}

int sbnd::HNLPiZeroAnalysis::GetTotalGenEvents(const art::Event &e)
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

bool sbnd::HNLPiZeroAnalysis::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
}

bool sbnd::HNLPiZeroAnalysis::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

void sbnd::HNLPiZeroAnalysis::SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
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
    fMCTruthHitsMap[mct]++;
  }

  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  art::fill_ptr_vector(pfpVec, pfpHandle);

  for(auto const& pfp : pfpVec)
    fPFPMap[pfp->Self()] = pfp;

}

void sbnd::HNLPiZeroAnalysis::AnalyseMeVPrtlTruth(const art::Event &e, const art::Handle<std::vector<evgen::ldm::MeVPrtlTruth>> mevptHandle)
{
  //if(fDebug)
  //  std::cout << "Analysing MeVPrtlTruth...";

  std::vector<art::Ptr<evgen::ldm::MeVPrtlTruth>> mevptVec;
  if (mevptHandle.isValid()) art::fill_ptr_vector(mevptVec, mevptHandle);

  //if(fDebug)
  //  std::cout << "Found " << mevptVec.size() << " hnl(s)..." << std::endl;

  for (auto const &mevpt: mevptVec){

    mevprtl_decay_pos_x.push_back(mevpt->decay_pos.X()); 
    mevprtl_decay_pos_y.push_back(mevpt->decay_pos.Y()); 
    mevprtl_decay_pos_z.push_back(mevpt->decay_pos.Z()); 
    mevprtl_decay_pos_t.push_back(mevpt->decay_pos.T());

    mevprtl_mom_x.push_back(mevpt->mevprtl_mom.X());
    mevprtl_mom_y.push_back(mevpt->mevprtl_mom.Y());
    mevprtl_mom_z.push_back(mevpt->mevprtl_mom.Z());
    mevprtl_mom_e.push_back(mevpt->mevprtl_mom.E());

    mevprtl_mass.push_back(mevpt->mass);

    mevprtl_flux_weight.push_back(mevpt->flux_weight);
    mevprtl_ray_weight.push_back(mevpt->ray_weight);
    mevprtl_decay_weight.push_back(mevpt->decay_weight);

    mevprtl_total_decay_width.push_back(mevpt->total_decay_width);
    mevprtl_total_mean_distance.push_back(mevpt->total_mean_distance);
    mevprtl_total_mean_lifetime.push_back(mevpt->total_mean_lifetime);
    mevprtl_allowed_decay_fraction.push_back(mevpt->allowed_decay_fraction);

    mevprtl_C1.push_back(mevpt->C1);
    mevprtl_C2.push_back(mevpt->C2);
    mevprtl_C3.push_back(mevpt->C3);
    mevprtl_C4.push_back(mevpt->C4);
    mevprtl_C5.push_back(mevpt->C5);
    _n_hnl++;
  }
}

void sbnd::HNLPiZeroAnalysis::ResizeMCTruth1DVector(const int col){
  nu_mctruth_id.resize(col);
  nu_event_type.resize(col);
  nu_signal.resize(col);
  nu_en_dep.resize(col);
  nu_ccnc.resize(col);
  nu_mode.resize(col);
  nu_int_type.resize(col);
  nu_w.resize(col);
  nu_x.resize(col);
  nu_y.resize(col);
  nu_q_sqr.resize(col);
  nu_pt.resize(col);
  nu_theta.resize(col);
  nu_e.resize(col);
  nu_vtx_x.resize(col);
  nu_vtx_y.resize(col);
  nu_vtx_z.resize(col);
  nu_vtx_t.resize(col);
  nu_n_protons.resize(col);
  nu_n_neutrons.resize(col);
  nu_n_charged_pions.resize(col);
  nu_n_neutral_pions.resize(col);
  nu_n_photons.resize(col);
  nu_n_others.resize(col);
}

void sbnd::HNLPiZeroAnalysis::AnalyseMCTruthHandle(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  //if(fDebug){
  //  std::cout << "Analysing MCTruth Handle...";
  //  std::cout << "Found " << MCTruthHandles.size() << " MCTruth Handle(s)..." << std::endl;
  //}

  for(auto const& MCTruthHandle : MCTruthHandles)
  {
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    art::fill_ptr_vector(MCTruthVec, MCTruthHandle);
    
    for(auto const& mct : MCTruthVec)
    {
      if (fMeVPrtl){
        if(mct->Origin() == 2) continue;
        ++_n_mctruth;
      }
      if (fBeamOff | fGenie){
	if(mct->Origin() != 1) continue;
        ++_n_mctruth;
      }
	
    }
  }

  ResizeMCTruth1DVector(_n_mctruth); 
 
  int mctCounter = 0;

  for(auto const& MCTruthHandle : MCTruthHandles)
  {
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    art::fill_ptr_vector(MCTruthVec, MCTruthHandle);
    
    for(auto const& mct : MCTruthVec)
    {
      if (fMeVPrtl){
        if(mct->Origin() == 2) continue;
      	AnalyseMCTruth(e, mct, mctCounter);
      	++mctCounter;
      }
      if (fBeamOff | fGenie){
	if(mct->Origin() != 1) continue;
      	AnalyseMCTruth(e, mct, mctCounter);
      }
    }
  }
}

void sbnd::HNLPiZeroAnalysis::AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct, const int mctCounter)
{
  //if(fDebug){
  //  std::cout << "Analysing MCTruth ...";
  //  std::cout << "This MCTruth Origin is " << mct->Origin() << "..." << std::endl;
  //}

  if (mct->Origin() == 1) _n_nu++; //beam neutrino counter
  
  const simb::MCNeutrino mcn = mct->GetNeutrino();
  const simb::MCParticle nu  = mcn.Nu();
  
  const bool nc = mcn.CCNC() == 1;
  const bool av = VolumeCheck(nu.Position().Vect());
  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);
  
  art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);
  
  //if(fDebug)
  //  std::cout << "Found " << MCParticleVec.size() << " MCParticle(s)..." << std::endl;
  
  int protons = 0, neutrons = 0, charged_pions = 0, neutral_pions = 0, photons = 0, others = 0;
  float trueEnDep = 0.;
  
  for(auto const& mcp : MCParticleVec)
  {
    if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
    {
      switch(abs(mcp->PdgCode()))
      { 
        case 2212:
          ++protons;
          break;
        case 2112:
          ++neutrons;
          break;
        case 211:
          ++charged_pions;
          break;
        case 111:
          ++neutral_pions;
          break;
        case 22:
          ++photons;
          break;
        default:
          ++others;
          break;
      }
    }
    std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());
    for(auto const& ide : ides)
      trueEnDep += ide->energy / 1000.;
  }
  
  const bool pizero = neutral_pions > 0;
  const bool _is_hnl = mct->Origin() == 0 && _n_hnl > 0; //MeVPrtl check always run before truth matching so n_hnl == 1 if running on HNL sample
  const bool _is_nu = mct->Origin() == 1? 1:0;
  const bool _is_cosmic = mct->Origin() == 2? 1:0;
  
  // Event Classification
  if (_is_hnl && fv)
  {
    nu_event_type[mctCounter] = (int) kHNL;
    nu_signal[mctCounter] = true;
  }else
  {
    nu_signal[mctCounter] = false;

    if(_is_hnl && !fv && av)
      nu_event_type[mctCounter] = (int) kHNLNonFV;
    else if(_is_hnl && !av)
      nu_event_type[mctCounter] = (int) kHNLDirt;
    else if(_is_nu && nc && pizero && fv)
      nu_event_type[mctCounter] = (int) kNCPiZero;
    else if(_is_nu && nc && !pizero && fv)
      nu_event_type[mctCounter] = (int) kOtherNC;
    else if(_is_nu && abs(nu.PdgCode()) == 14 && fv)
      nu_event_type[mctCounter] = (int) kCCNuMu;
    else if(_is_nu && abs(nu.PdgCode()) == 12 && fv)
      nu_event_type[mctCounter] = (int) kCCNuE;
    else if(_is_nu && !fv && av)
      nu_event_type[mctCounter] = (int) kNuNonFV;
    else if(_is_nu && !av)
      nu_event_type[mctCounter] = (int) kNuDirt;
    else if(_is_cosmic)
      nu_event_type[mctCounter] = (int) kCosmic;
    else
      nu_event_type[mctCounter] = (int) kUnknownEv;
  }

  nu_mctruth_id[mctCounter] = mct.key();
  nu_en_dep[mctCounter] = trueEnDep;
  nu_ccnc[mctCounter] = mcn.CCNC();
  nu_mode[mctCounter] = mcn.Mode();
  nu_int_type[mctCounter] = mcn.InteractionType();
  nu_w[mctCounter] = mcn.W();
  nu_x[mctCounter] = mcn.X();
  nu_y[mctCounter] = mcn.Y();
  nu_q_sqr[mctCounter] = mcn.QSqr();
  nu_pt[mctCounter] = mcn.Pt();
  nu_theta[mctCounter] = mcn.Theta();
  nu_e[mctCounter] = nu.E();
  nu_vtx_x[mctCounter] = nu.Vx();
  nu_vtx_y[mctCounter] = nu.Vy();
  nu_vtx_z[mctCounter] = nu.Vz();
  nu_vtx_t[mctCounter] = nu.T();
  nu_n_protons[mctCounter] = protons;  
  nu_n_neutrons[mctCounter] = neutrons;  
  nu_n_charged_pions[mctCounter] = charged_pions;  
  nu_n_neutral_pions[mctCounter] = neutral_pions;  
  nu_n_photons[mctCounter] = photons;  
  nu_n_others[mctCounter] = others;  
}

art::Ptr<recob::PFParticle> sbnd::HNLPiZeroAnalysis::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
{
  for(auto pfp : pfps)
    if(pfp->IsPrimary())
      return pfp;

  return art::Ptr<recob::PFParticle>();
}

void sbnd::HNLPiZeroAnalysis::ResizeSlice1DVector(const int col){

  slc_id.resize(col, -99999);
  slc_n_pfps.resize(col, -99999);
  slc_is_clear_cosmics.resize(col, -99999);
  slc_primary_pfp_id.resize(col, -99999);
  slc_primary_pfp_pdg.resize(col, -99999);
  slc_n_primary_daughters.resize(col, -99999);
  slc_vtx_x.resize(col, -99999);
  slc_vtx_y.resize(col, -99999);
  slc_vtx_z.resize(col, -99999);
  slc_is_fv.resize(col, -99999);
  slc_crumbs_score.resize(col, -99999);
  slc_crumbs_nc_score.resize(col, -99999);
  slc_crumbs_ccnue_score.resize(col, -99999);
  slc_opt0_time.resize(col, -99999);
  slc_opt0_time_corrected_Z_pandora.resize(col, -99999);
  slc_opt0_score.resize(col, -99999);
  slc_opt0_measPE.resize(col, -99999);
  slc_opt0_hypoPE.resize(col, -99999);
  slc_opt0_frac.resize(col, -99999);

  slc_n_trks.resize(col, -99999);
  slc_n_shws.resize(col, -99999);
  slc_n_primary_trks.resize(col, -99999);
  slc_n_primary_shws.resize(col, -99999);
  slc_n_dazzle_muons.resize(col, -99999);
  slc_n_dazzle_pions.resize(col, -99999);
  slc_n_dazzle_pions_thresh.resize(col, -99999);
  slc_n_dazzle_protons.resize(col, -99999);
  slc_n_dazzle_protons_thresh.resize(col, -99999);
  slc_n_dazzle_other.resize(col, -99999);
  slc_n_primary_dazzle_muons.resize(col, -99999);
  slc_n_primary_dazzle_pions.resize(col, -99999);
  slc_n_primary_dazzle_pions_thresh.resize(col, -99999);
  slc_n_primary_dazzle_protons.resize(col, -99999);
  slc_n_primary_dazzle_protons_thresh.resize(col, -99999);
  slc_n_primary_dazzle_other.resize(col, -99999);
  slc_n_razzle_electrons.resize(col, -99999);
  slc_n_razzle_photons.resize(col, -99999);
  slc_n_razzle_other.resize(col, -99999);
  slc_n_razzled_electrons.resize(col, -99999);
  slc_n_razzled_muons.resize(col, -99999);
  slc_n_razzled_photons.resize(col, -99999);
  slc_n_razzled_pions.resize(col, -99999);
  slc_n_razzled_pions_thresh.resize(col, -99999);
  slc_n_razzled_protons.resize(col, -99999);
  slc_n_razzled_protons_thresh.resize(col, -99999);
  slc_n_primary_razzle_electrons.resize(col, -99999);
  slc_n_primary_razzle_photons.resize(col, -99999);
  slc_n_primary_razzle_other.resize(col, -99999);
  slc_n_primary_razzled_electrons.resize(col, -99999);
  slc_n_primary_razzled_muons.resize(col, -99999);
  slc_n_primary_razzled_photons.resize(col, -99999);
  slc_n_primary_razzled_pions.resize(col, -99999);
  slc_n_primary_razzled_pions_thresh.resize(col, -99999);
  slc_n_primary_razzled_protons.resize(col, -99999);
  slc_n_primary_razzled_protons_thresh.resize(col, -99999);
  slc_n_stub.resize(col, -99999);
  slc_total_track_E.resize(col, -99999);
  slc_total_shower_E.resize(col, -99999);

  slc_comp.resize(col, -99999);
  slc_pur.resize(col, -99999);
  slc_true_mctruth_id.resize(col, 99999);
  slc_true_event_type.resize(col, -1);
  slc_true_en_dep.resize(col, -99999);
  slc_true_vtx_x.resize(col, -99999);
  slc_true_vtx_y.resize(col, -99999);
  slc_true_vtx_z.resize(col, -99999);
  slc_true_vtx_t.resize(col, -99999);
  slc_true_vtx_t_corrected_Z.resize(col, -99999);

  if (fPiZeroAna){
    slc_pizero_good.resize(col, 0);
    slc_photon0_pfp_id.resize(col, -99999);
    slc_photon1_pfp_id.resize(col, -99999);
    slc_photon2_pfp_id.resize(col, -99999);
    slc_photon0_energy.resize(col, -99999);
    slc_photon1_energy.resize(col, -99999);
    slc_photon2_energy.resize(col, -99999);
    slc_pizero_open_angle.resize(col, -99999);
    slc_pizero_mass.resize(col, -99999);
    slc_pizero_mom.resize(col, -99999);
    slc_pizero_CosThetaZ.resize(col, -99999);

    slc_true_pizero_mom.resize(col, -99999);
    slc_true_pizero_cosThetaZ.resize(col, -99999);

  }
}

void sbnd::HNLPiZeroAnalysis::ResizeSlice2DVectorRow(const int row){

   slc_pfp_id.resize(row); 
   slc_pfp_pdg.resize(row);
   slc_pfp_n_children.resize(row);
   slc_pfp_primary.resize(row);
   slc_pfp_primary_child.resize(row);
   slc_pfp_n_hits.resize(row);
   slc_pfp_n_sps.resize(row);
   slc_pfp_track_score.resize(row); 
   slc_pfp_good_track.resize(row); 
   slc_pfp_good_shower.resize(row);
   slc_pfp_comp.resize(row);
   slc_pfp_pur.resize(row);
   slc_pfp_cnnscore_track.resize(row);
   slc_pfp_cnnscore_shower.resize(row);
   slc_pfp_cnnscore_noise.resize(row);
   slc_pfp_cnnscore_michel.resize(row);
   slc_pfp_cnnscore_endmichel.resize(row);
   slc_pfp_cnnscore_nclusters.resize(row);
   slc_pfp_razzled_electron_score.resize(row);
   slc_pfp_razzled_muon_score.resize(row);
   slc_pfp_razzled_photon_score.resize(row);
   slc_pfp_razzled_pion_score.resize(row);
   slc_pfp_razzled_proton_score.resize(row);
   slc_pfp_razzled_pdg.resize(row);

   slc_pfp_true_trackid.resize(row);
   slc_pfp_true_pdg.resize(row);
   slc_pfp_true_energy.resize(row);
   slc_pfp_true_p_x.resize(row);
   slc_pfp_true_p_y.resize(row);
   slc_pfp_true_p_z.resize(row);

   slc_pfp_track_start_x.resize(row); 
   slc_pfp_track_start_y.resize(row); 
   slc_pfp_track_start_z.resize(row); 
   slc_pfp_track_end_x.resize(row); 
   slc_pfp_track_end_y.resize(row); 
   slc_pfp_track_end_z.resize(row); 
   slc_pfp_track_dir_x.resize(row); 
   slc_pfp_track_dir_y.resize(row); 
   slc_pfp_track_dir_z.resize(row); 

   slc_pfp_track_length.resize(row);
   slc_pfp_track_contained.resize(row);
   slc_pfp_track_dazzle_muon_score.resize(row);
   slc_pfp_track_dazzle_pion_score.resize(row);
   slc_pfp_track_dazzle_proton_score.resize(row);
   slc_pfp_track_dazzle_other_score.resize(row);
   slc_pfp_track_dazzle_pdg.resize(row);
   slc_pfp_track_ke.resize(row);
   slc_pfp_track_charge.resize(row);
   slc_pfp_track_chi2_muon.resize(row);
   slc_pfp_track_chi2_pion.resize(row);
   slc_pfp_track_chi2_kaon.resize(row);
   slc_pfp_track_chi2_proton.resize(row);
   slc_pfp_track_chi2_pdg.resize(row);
   slc_pfp_track_mcs_mean_scatter.resize(row);
   slc_pfp_track_mcs_max_scatter_ratio.resize(row);
   slc_pfp_track_range_p.resize(row);
   slc_pfp_track_closest_approach_mean_dca.resize(row);
   slc_pfp_track_stopping_dedx_chi2_ratio.resize(row);
   slc_pfp_track_stopping_dedx_pol0_fit.resize(row);
   slc_pfp_track_theta.resize(row);
   slc_pfp_track_phi.resize(row);

   slc_pfp_shower_start_x.resize(row);
   slc_pfp_shower_start_y.resize(row);
   slc_pfp_shower_start_z.resize(row);
   slc_pfp_shower_end_x.resize(row);
   slc_pfp_shower_end_y.resize(row);
   slc_pfp_shower_end_z.resize(row);
   slc_pfp_shower_conv_gap.resize(row);
   slc_pfp_shower_dir_x.resize(row);
   slc_pfp_shower_dir_y.resize(row);
   slc_pfp_shower_dir_z.resize(row);
   slc_pfp_shower_length.resize(row);
   slc_pfp_shower_open_angle.resize(row);
   slc_pfp_shower_energy.resize(row);
   slc_pfp_shower_dedx.resize(row);
   slc_pfp_shower_sqrt_energy_density.resize(row);
   slc_pfp_shower_modified_hit_density.resize(row);
   slc_pfp_shower_razzle_electron_score.resize(row);
   slc_pfp_shower_razzle_photon_score.resize(row);
   slc_pfp_shower_razzle_other_score.resize(row);
   slc_pfp_shower_razzle_pdg.resize(row);
   slc_pfp_shower_cosmic_dist.resize(row);
   slc_pfp_shower_track_length.resize(row);
   slc_pfp_shower_track_width.resize(row);
   slc_pfp_shower_density_grad.resize(row);
   slc_pfp_shower_density_pow.resize(row);
   slc_pfp_shower_contained.resize(row);
   slc_pfp_shower_theta.resize(row);
   slc_pfp_shower_phi.resize(row);
}

void sbnd::HNLPiZeroAnalysis::ResizeSlice2DVectorCol(const int row, const int col){

   slc_pfp_id[row].resize(col, -99999); 
   slc_pfp_pdg[row].resize(col, -99999);
   slc_pfp_n_children[row].resize(col, -99999);
   slc_pfp_primary[row].resize(col, -99999);
   slc_pfp_primary_child[row].resize(col, -99999);
   slc_pfp_n_hits[row].resize(col, 0);
   slc_pfp_n_sps[row].resize(col, -99999);
   slc_pfp_track_score[row].resize(col, -99999); 
   slc_pfp_good_track[row].resize(col, -99999); 
   slc_pfp_good_shower[row].resize(col, -99999);
   slc_pfp_comp[row].resize(col, -99999);
   slc_pfp_pur[row].resize(col, -99999);
   slc_pfp_cnnscore_track[row].resize(col, -99999);
   slc_pfp_cnnscore_shower[row].resize(col, -99999);
   slc_pfp_cnnscore_noise[row].resize(col, -99999);
   slc_pfp_cnnscore_michel[row].resize(col, -99999);
   slc_pfp_cnnscore_endmichel[row].resize(col, -99999);
   slc_pfp_cnnscore_nclusters[row].resize(col, -99999);
   slc_pfp_razzled_electron_score[row].resize(col, -99999);
   slc_pfp_razzled_muon_score[row].resize(col, -99999);
   slc_pfp_razzled_photon_score[row].resize(col, -99999);
   slc_pfp_razzled_pion_score[row].resize(col, -99999);
   slc_pfp_razzled_proton_score[row].resize(col, -99999);
   slc_pfp_razzled_pdg[row].resize(col, -99999);

   slc_pfp_true_trackid[row].resize(col, -99999);
   slc_pfp_true_pdg[row].resize(col, -99999);
   slc_pfp_true_energy[row].resize(col, -99999);
   slc_pfp_true_p_x[row].resize(col, -99999);
   slc_pfp_true_p_y[row].resize(col, -99999);
   slc_pfp_true_p_z[row].resize(col, -99999);

   slc_pfp_track_start_x[row].resize(col, -99999); 
   slc_pfp_track_start_y[row].resize(col, -99999); 
   slc_pfp_track_start_z[row].resize(col, -99999); 
   slc_pfp_track_end_x[row].resize(col, -99999); 
   slc_pfp_track_end_y[row].resize(col, -99999); 
   slc_pfp_track_end_z[row].resize(col, -99999); 
   slc_pfp_track_dir_x[row].resize(col, -99999); 
   slc_pfp_track_dir_y[row].resize(col, -99999); 
   slc_pfp_track_dir_z[row].resize(col, -99999); 

   slc_pfp_track_length[row].resize(col, -99999);
   slc_pfp_track_contained[row].resize(col, 0); 
   slc_pfp_track_dazzle_muon_score[row].resize(col, -99999);
   slc_pfp_track_dazzle_pion_score[row].resize(col, -99999);
   slc_pfp_track_dazzle_proton_score[row].resize(col, -99999);
   slc_pfp_track_dazzle_other_score[row].resize(col, -99999);
   slc_pfp_track_dazzle_pdg[row].resize(col, -99999);
   slc_pfp_track_ke[row].resize(col, -99999);
   slc_pfp_track_charge[row].resize(col, -99999);
   slc_pfp_track_chi2_muon[row].resize(col, -99999);
   slc_pfp_track_chi2_pion[row].resize(col, -99999);
   slc_pfp_track_chi2_kaon[row].resize(col, -99999);
   slc_pfp_track_chi2_proton[row].resize(col, -99999);
   slc_pfp_track_chi2_pdg[row].resize(col, -99999);
   slc_pfp_track_mcs_mean_scatter[row].resize(col, -99999);
   slc_pfp_track_mcs_max_scatter_ratio[row].resize(col, -99999);
   slc_pfp_track_range_p[row].resize(col, -99999);
   slc_pfp_track_closest_approach_mean_dca[row].resize(col, -99999);
   slc_pfp_track_stopping_dedx_chi2_ratio[row].resize(col, -99999);
   slc_pfp_track_stopping_dedx_pol0_fit[row].resize(col, -99999);
   slc_pfp_track_theta[row].resize(col, -99999);
   slc_pfp_track_phi[row].resize(col, -99999);

   slc_pfp_shower_start_x[row].resize(col, -99999);
   slc_pfp_shower_start_y[row].resize(col, -99999);
   slc_pfp_shower_start_z[row].resize(col, -99999);
   slc_pfp_shower_end_x[row].resize(col, -99999);
   slc_pfp_shower_end_y[row].resize(col, -99999);
   slc_pfp_shower_end_z[row].resize(col, -99999);
   slc_pfp_shower_conv_gap[row].resize(col, -99999);
   slc_pfp_shower_dir_x[row].resize(col, -99999);
   slc_pfp_shower_dir_y[row].resize(col, -99999);
   slc_pfp_shower_dir_z[row].resize(col, -99999);
   slc_pfp_shower_length[row].resize(col, -99999);
   slc_pfp_shower_open_angle[row].resize(col, -99999);
   slc_pfp_shower_energy[row].resize(col, -99999);
   slc_pfp_shower_dedx[row].resize(col, -99999);
   slc_pfp_shower_sqrt_energy_density[row].resize(col, -99999);
   slc_pfp_shower_modified_hit_density[row].resize(col, -99999);
   slc_pfp_shower_razzle_electron_score[row].resize(col, -99999);
   slc_pfp_shower_razzle_photon_score[row].resize(col, -99999);
   slc_pfp_shower_razzle_other_score[row].resize(col, -99999);
   slc_pfp_shower_razzle_pdg[row].resize(col, -99999);
   slc_pfp_shower_cosmic_dist[row].resize(col, -99999);
   slc_pfp_shower_track_length[row].resize(col, -99999);
   slc_pfp_shower_track_width[row].resize(col, -99999);
   slc_pfp_shower_density_grad[row].resize(col, -99999);
   slc_pfp_shower_density_pow[row].resize(col, -99999);
   slc_pfp_shower_contained[row].resize(col, 0);
   slc_pfp_shower_theta[row].resize(col, -99999);
   slc_pfp_shower_phi[row].resize(col, -99999);
}

void sbnd::HNLPiZeroAnalysis::ResizeWeight2DVectorRow(const int row){
  slc_flux_weight_expskin.resize(row);
  slc_flux_weight_horncurrent.resize(row);
  slc_flux_weight_kminus.resize(row);
  slc_flux_weight_kplus.resize(row);
  slc_flux_weight_kzero.resize(row);
  slc_flux_weight_nucleoninexsec.resize(row);
  slc_flux_weight_nucleonqexsec.resize(row);
  slc_flux_weight_nucleontotxsec.resize(row);
  slc_flux_weight_piminus.resize(row);
  slc_flux_weight_pioninexsex.resize(row);
  slc_flux_weight_pionqexsec.resize(row);
  slc_flux_weight_piontotxsec.resize(row);
  slc_flux_weight_piplus.resize(row);
  slc_flux_weight_total.resize(row);

  slc_xsec_unisim_DecayAngMEC.resize(row);
  slc_xsec_unisim_ThetaDelta2NRad.resize(row); 
  slc_xsec_unisim_Theta_Delta2Npi.resize(row);
  slc_xsec_unisim_VecFFCCQEshape.resize(row);

  slc_xsec_multisigma_CoulombCCQE.resize(row);
  slc_xsec_multisigma_NonRESBGvbarnCC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarnCC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarnNC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarnNC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarpCC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarpCC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarpNC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvbarpNC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvnCC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvnCC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvnNC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvnNC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvpCC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvpCC2pi.resize(row);
  slc_xsec_multisigma_NonRESBGvpNC1pi.resize(row);
  slc_xsec_multisigma_NonRESBGvpNC2pi.resize(row);
  slc_xsec_multisigma_NormCCMEC.resize(row);
  slc_xsec_multisigma_NormNCMEC.resize(row);
  slc_xsec_multisigma_RDecBR1eta.resize(row);
  slc_xsec_multisigma_RDecBR1gamma.resize(row);
  slc_xsec_multisigma_RPA_CCQE.resize(row);
  slc_xsec_multisigma_NormCCCOH.resize(row);        
  slc_xsec_multisigma_NormNCCOH.resize(row);        

  slc_xsec_multisim_ZExpA.resize(row);
  slc_xsec_multisim_NCEL.resize(row);
  slc_xsec_multisim_CCRES.resize(row);
  slc_xsec_multisim_NCRES.resize(row);
  slc_xsec_multisim_DISBY.resize(row);
  slc_xsec_multisim_FSI_pi.resize(row);
  slc_xsec_multisim_FSI_N.resize(row);
  slc_xsec_multisim_total.resize(row);
  
  slc_geant4_multisim_reinteractions.resize(row);
}

void sbnd::HNLPiZeroAnalysis::ResizeWeight2DVectorCol(const int row){
  //flux weight multi-sim 1000 universes
  slc_flux_weight_expskin[row].resize(1000,-9999);
  slc_flux_weight_horncurrent[row].resize(1000,-9999);
  slc_flux_weight_kminus[row].resize(1000,-9999);
  slc_flux_weight_kplus[row].resize(1000,-9999);
  slc_flux_weight_kzero[row].resize(1000,-9999);
  slc_flux_weight_nucleoninexsec[row].resize(1000,-9999);
  slc_flux_weight_nucleonqexsec[row].resize(1000,-9999);
  slc_flux_weight_nucleontotxsec[row].resize(1000,-9999);
  slc_flux_weight_piminus[row].resize(1000,-9999);
  slc_flux_weight_pioninexsex[row].resize(1000,-9999);
  slc_flux_weight_pionqexsec[row].resize(1000,-9999);
  slc_flux_weight_piontotxsec[row].resize(1000,-9999);
  slc_flux_weight_piplus[row].resize(1000,-9999);

  slc_flux_weight_total[row].resize(1000, 1);

  //xsection unisim 1 direction variation
  slc_xsec_unisim_DecayAngMEC[row].resize(1, -9999);
  slc_xsec_unisim_ThetaDelta2NRad[row].resize(1, -9999); 
  slc_xsec_unisim_Theta_Delta2Npi[row].resize(1, -9999);
  slc_xsec_unisim_VecFFCCQEshape[row].resize(1, -9999);

  //xsection multisigmal +/- 1 sigma
  slc_xsec_multisigma_CoulombCCQE[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarnCC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarnCC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarnNC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarnNC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarpCC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarpCC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarpNC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvbarpNC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvnCC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvnCC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvnNC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvnNC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvpCC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvpCC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvpNC1pi[row].resize(6, -9999);
  slc_xsec_multisigma_NonRESBGvpNC2pi[row].resize(6, -9999);
  slc_xsec_multisigma_NormCCMEC[row].resize(6, -9999);
  slc_xsec_multisigma_NormNCMEC[row].resize(6, -9999);
  slc_xsec_multisigma_RDecBR1eta[row].resize(6, -9999);
  slc_xsec_multisigma_RDecBR1gamma[row].resize(6, -9999);
  slc_xsec_multisigma_RPA_CCQE[row].resize(6, -9999);
  slc_xsec_multisigma_NormCCCOH[row].resize(6, -9999);        
  slc_xsec_multisigma_NormNCCOH[row].resize(6, -9999);        

  //xsec 500 universe
  slc_xsec_multisim_ZExpA[row].resize(500,-9999);
  slc_xsec_multisim_NCEL[row].resize(500,-9999);
  slc_xsec_multisim_CCRES[row].resize(500,-9999);
  slc_xsec_multisim_NCRES[row].resize(500,-9999);
  slc_xsec_multisim_DISBY[row].resize(500,-9999);
  slc_xsec_multisim_FSI_pi[row].resize(500,-9999);
  slc_xsec_multisim_FSI_N[row].resize(500,-9999);

  slc_xsec_multisim_total[row].resize(500, 1);
            
  slc_geant4_multisim_reinteractions[row].resize(1000, -9999);
}

void sbnd::HNLPiZeroAnalysis::AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                           const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                                           const art::Handle<std::vector<recob::Track>> &trackHandle,
                                           const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  if(fDebug)
    std::cout << "Analysing slices...";

  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);
    
  if(fDebug)
      std::cout << "Found " << sliceVec.size() << " slice(s)..." << std::endl;

  _n_slc = sliceVec.size();

  ResizeSlice1DVector(_n_slc);
  ResizeSlice2DVectorRow(_n_slc);
  ResizeWeight2DVectorRow(_n_slc); //Flux weight default has 1000 universes

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>      pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  art::FindOneP<sbn::CRUMBSResult>  slicesToCRUMBS(sliceHandle, e, fCRUMBSModuleLabel);
  art::FindManyP<sbn::OpT0Finder>   slicesToOpT0(sliceHandle, e, fOpT0ModuleLabel);
  art::FindManyP<sbn::Stub> 	    slicesToStubs(sliceHandle, e, fStubModuleLabel);
  
  for (auto&& [slcCounter, slc] : enumerate(sliceVec))
  {
    const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slc.key());
    slc_n_pfps[slcCounter] = pfps.size();

    if(fDebug){
      std::cout << std::endl << "Slice #" << slcCounter << std::endl;
      std::cout << "Found " << pfps.size() << " pfp(s) with this slice..." << std::endl;
    }

    slc_id[slcCounter] = slcCounter;

    if(pfps.size() == 0)
    {
      slc_is_clear_cosmics[slcCounter] = true;
      continue;
    }

    const art::Ptr<recob::PFParticle> prim = GetPrimaryPFP(pfps);

    if(prim.isNull())
      continue;

    slc_primary_pfp_id[slcCounter] = prim->Self();
    slc_primary_pfp_pdg[slcCounter] = prim->PdgCode();
    slc_n_primary_daughters[slcCounter] = prim->NumDaughters();
    
    if(abs(prim->PdgCode()) == 13 || abs(prim->PdgCode()) == 11)
	slc_is_clear_cosmics[slcCounter] = true;
    else
	slc_is_clear_cosmics[slcCounter] = false;

    const art::Ptr<recob::Vertex> vtx = pfpToVertices.at(prim.key());
    geo::Point_t vtxPos = vtx.isNonnull() ? vtx->position() : geo::Point_t(def_double, def_double, def_double);
    slc_vtx_x[slcCounter] = vtxPos.X();
    slc_vtx_y[slcCounter] = vtxPos.Y();
    slc_vtx_z[slcCounter] = vtxPos.Z();
    slc_is_fv[slcCounter] = VolumeCheck(vtxPos, 20., 5., 10., 50.);

    const art::Ptr<sbn::CRUMBSResult> crumbs = slicesToCRUMBS.at(slc.key());
    if(crumbs.isNonnull())
    {
      slc_crumbs_score[slcCounter] = crumbs->score;  
      slc_crumbs_nc_score[slcCounter] = crumbs->ncscore;  
      slc_crumbs_ccnue_score[slcCounter] = crumbs->ccnuescore;  
    }

    std::vector<art::Ptr<sbn::Stub>> stubsVec;
    if (slicesToStubs.isValid()) {
      stubsVec = slicesToStubs.at(0);
    } 

    // Lookup stubs to overlaid PFP
    art::FindManyP<recob::PFParticle> stubsToPFPs(stubsVec, e, fStubModuleLabel);

    int nstub = 0;

    for (size_t iStub = 0; iStub < stubsVec.size(); iStub++) {

      art::Ptr<recob::PFParticle> thisStubPFP;

      if (!stubsToPFPs.at(iStub).empty()){

        thisStubPFP = stubsToPFPs.at(iStub).at(0);

        for(auto&& [pfpCounter, id]: enumerate(prim->Daughters())){
          if (thisStubPFP->Self() == id) {
            nstub++;
            break;
          }
        }
      } 

      //std::cout << "primary pfp id = " << prim->Self() << std::endl;
      //
    }

    slc_n_stub[slcCounter] = nstub;

    std::vector<art::Ptr<sbn::OpT0Finder>> opT0Vec = slicesToOpT0.at(slc.key());
    if(opT0Vec.size() > 0)
    {
      std::sort(opT0Vec.begin(), opT0Vec.end(),
                [](auto const& a, auto const& b)
                { return a->score > b->score; });
      slc_opt0_time[slcCounter] = opT0Vec[0]->time;
      slc_opt0_score[slcCounter] = opT0Vec[0]->score;
      slc_opt0_measPE[slcCounter] = opT0Vec[0]->measPE;
      slc_opt0_hypoPE[slcCounter] = opT0Vec[0]->hypoPE;

      slc_opt0_frac[slcCounter] = (opT0Vec[0]->hypoPE - opT0Vec[0]->measPE) / opT0Vec[0]->measPE;

      // Add corrected Z vertex variable
      slc_opt0_time_corrected_Z_pandora [slcCounter] = slc_opt0_time[slcCounter] - slc_vtx_z[slcCounter]/light_speed_cm_us;

    }

    ResizeSlice2DVectorCol(slcCounter, prim->Daughters().size());
    AnalysePFPs(e, slcCounter, prim, vtx, pfpHandle, trackHandle, showerHandle);
    AnalyseSliceTruth(e, slc, slcCounter, sliceHandle);

    //Analyse Pi0 kinematics on slice basis once all the shower variables are filled
    if(fPiZeroAna) AnalysePiZero(slcCounter);

  } // End of loop sliceVec
}

void sbnd::HNLPiZeroAnalysis::AnalysePFPs(const art::Event &e, const int slcCounter,
					const art::Ptr<recob::PFParticle> &prim, 
					const art::Ptr<recob::Vertex> &vtx, 
					const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, 
					const art::Handle<std::vector<recob::Track>> &trackHandle,
                                        const art::Handle<std::vector<recob::Shower>> &showerHandle
					)
{

  if(fDebug){
    std::cout << std::endl << "Analysing Primary PFP..." << std::endl;
    std::cout << "Found " << prim->Daughters().size() << " pfp daughters..."<< std::endl;
  }
  
  art::FindOneP<recob::Track>                      pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpToMeta(pfpHandle, e, fPFParticleModuleLabel);
  art::FindManyP<recob::Hit>                       showersToHits(showerHandle, e, fShowerModuleLabel);
  art::FindOneP<sbn::MVAPID>                       pfpToRazzled(pfpHandle, e, fRazzledModuleLabel);
  art::FindManyP<recob::SpacePoint>                pfpToSpacePoints(pfpHandle, e, fSpacePointModuleLabel);
  art::FindOneP<sbn::PFPCNNScore> 		   pfpToCNNScore(pfpHandle, e, fCNNScoreModuleLabel);
  art::FindManyP<sbn::Stub> 			   pfpToStub(pfpHandle, e, fStubModuleLabel);

  int ntrks = 0, nshws = 0, ndazzlemuons = 0, ndazzlepions = 0, ndazzlepionsthresh = 0, ndazzleprotons = 0,
    ndazzleprotonsthresh = 0, ndazzleother = 0, nrazzleelectrons = 0, nrazzlephotons = 0, nrazzleother = 0,
    nrazzledelectrons = 0, nrazzledmuons = 0, nrazzledphotons = 0, nrazzledpions = 0, nrazzledpionsthresh = 0,
    nrazzledprotons = 0, nrazzledprotonsthresh = 0;

  int nprimtrks = 0, nprimshws = 0, nprimdazzlemuons = 0, nprimdazzlepions = 0, nprimdazzlepionsthresh = 0, nprimdazzleprotons = 0,
    nprimdazzleprotonsthresh = 0, nprimdazzleother = 0, nprimrazzleelectrons = 0, nprimrazzlephotons = 0, nprimrazzleother = 0,
    nprimrazzledelectrons = 0, nprimrazzledmuons = 0, nprimrazzledphotons = 0, nprimrazzledpions = 0, nprimrazzledpionsthresh = 0,
    nprimrazzledprotons = 0, nprimrazzledprotonsthresh = 0;

  for(auto&& [pfpCounter, id]: enumerate(prim->Daughters()))
  {
    //if(fDebug)
    //  std::cout << "Looking at slc#" << slcCounter << " & pfp#" << pfpCounter << std::endl;

    const art::Ptr<recob::PFParticle> pfp = fPFPMap[id];
    slc_pfp_id[slcCounter][pfpCounter] = id;
    slc_pfp_pdg[slcCounter][pfpCounter] = pfp->PdgCode(); 
    slc_pfp_n_children[slcCounter][pfpCounter] = pfp->NumDaughters();

    const bool primary_child = prim->Self() == pfp->Parent();
    slc_pfp_primary[slcCounter][pfpCounter] = pfp->IsPrimary();
    slc_pfp_primary_child[slcCounter][pfpCounter] = primary_child;

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
        slc_pfp_track_score[slcCounter][pfpCounter] = scoreIter->second;

    const art::Ptr<recob::Track> track   = pfpToTrack.at(pfp.key());
    const art::Ptr<recob::Shower> shower = pfpToShower.at(pfp.key());

    slc_pfp_good_track[slcCounter][pfpCounter] = track.isNonnull();
    slc_pfp_good_shower[slcCounter][pfpCounter] = shower.isNonnull();

    const art::Ptr<sbn::MVAPID> razzled = pfpToRazzled.at(pfp.key());
    if(razzled.isNonnull())
      ExtractRazzled(razzled, slcCounter, pfpCounter);

    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    const std::vector<art::Ptr<recob::Hit>> hits = showersToHits.at(shower.key());

    slc_pfp_n_hits[slcCounter][pfpCounter] = hits.size();

    const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, true);

    slc_pfp_true_trackid[slcCounter][pfpCounter] = trackID;
    slc_pfp_comp[slcCounter][pfpCounter] = Completeness(e, hits, trackID);
    slc_pfp_pur[slcCounter][pfpCounter] = Purity(e, hits, trackID);

    // TODO: Why Seg Fault?
    //std::cout << "pfp2cnn valid? = " << pfpToCNNScore.isValid() << std::endl;
    //if(pfpToCNNScore.isValid() && !slc_is_clear_cosmics[slcCounter]){
    //  const sbn::PFPCNNScore *cnnscore = pfpToCNNScore.at(id).get();
    //  ExtractCNNScores(cnnscore, slcCounter, pfpCounter);
    //}

    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpToSpacePoints.at(pfp.key());
    slc_pfp_n_sps[slcCounter][pfpCounter] = spacepoints.size();
    
    if(trackID != def_int)
    {
      const simb::MCParticle* mcp = particleInv->TrackIdToParticle_P(trackID);
      if(mcp != NULL)
        {
          slc_pfp_true_pdg[slcCounter][pfpCounter] = mcp->PdgCode();
          slc_pfp_true_energy[slcCounter][pfpCounter] = mcp->E();
          slc_pfp_true_p_x[slcCounter][pfpCounter] = mcp->Px();
          slc_pfp_true_p_y[slcCounter][pfpCounter] = mcp->Py();
          slc_pfp_true_p_z[slcCounter][pfpCounter] = mcp->Pz();
        }
    }

    if(track.isNonnull())
      AnalyseTrack(e, track, slcCounter, pfpCounter, trackHandle);

    double pfpenergy = def_double;
    
    if(pfp->PdgCode() == 13)
    {
      int dazzlepdg = slc_pfp_track_dazzle_pdg[slcCounter][pfpCounter];
      float trkenergy = slc_pfp_track_ke[slcCounter][pfpCounter];
      pfpenergy = trkenergy;

      if(dazzlepdg == 13)
        ++ndazzlemuons;
      if(dazzlepdg == 211)
        ++ndazzlepions;
      if(dazzlepdg == 211 && pfpenergy > 32.1)
        ++ndazzlepionsthresh;
      if(dazzlepdg == 2212)
        ++ndazzleprotons;
      if(dazzlepdg == 2212 && pfpenergy > 32.7)
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
        int razzlepdg = slc_pfp_shower_razzle_pdg[slcCounter][pfpCounter];
        pfpenergy = slc_pfp_shower_energy[slcCounter][pfpCounter];

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

     int razzledpdg = slc_pfp_razzled_pdg[slcCounter][pfpCounter];

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
  } // End of loop pfp daughters
    
  slc_n_trks[slcCounter] = ntrks;
  slc_n_shws[slcCounter] = nshws;
  slc_n_dazzle_muons[slcCounter] = ndazzlemuons;
  slc_n_dazzle_pions[slcCounter] = ndazzlepions;
  slc_n_dazzle_pions_thresh[slcCounter] = ndazzlepionsthresh;
  slc_n_dazzle_protons[slcCounter] = ndazzleprotons;
  slc_n_dazzle_protons_thresh[slcCounter] = ndazzleprotonsthresh;
  slc_n_dazzle_other[slcCounter] = ndazzleother;
  slc_n_razzle_electrons[slcCounter] = nrazzleelectrons;
  slc_n_razzle_photons[slcCounter] = nrazzlephotons;
  slc_n_razzle_other[slcCounter] = nrazzleother;
  slc_n_razzled_electrons[slcCounter] = nrazzledelectrons;
  slc_n_razzled_muons[slcCounter] = nrazzledmuons;
  slc_n_razzled_photons[slcCounter] = nrazzledphotons;
  slc_n_razzled_pions[slcCounter] = nrazzledpions;
  slc_n_razzled_pions_thresh[slcCounter] = nrazzledpionsthresh;
  slc_n_razzled_protons[slcCounter] = nrazzledprotons;
  slc_n_razzled_protons_thresh[slcCounter] = nrazzledprotonsthresh;

  slc_n_primary_trks[slcCounter] = nprimtrks;
  slc_n_primary_shws[slcCounter] = nprimshws;
  slc_n_primary_dazzle_muons[slcCounter] = nprimdazzlemuons;
  slc_n_primary_dazzle_pions[slcCounter] = nprimdazzlepions;
  slc_n_primary_dazzle_pions_thresh[slcCounter] = nprimdazzlepionsthresh;
  slc_n_primary_dazzle_protons[slcCounter] = nprimdazzleprotons;
  slc_n_primary_dazzle_protons_thresh[slcCounter] = nprimdazzleprotonsthresh;
  slc_n_primary_dazzle_other[slcCounter] = nprimdazzleother;
  slc_n_primary_razzle_electrons[slcCounter] = nprimrazzleelectrons;
  slc_n_primary_razzle_photons[slcCounter] = nprimrazzlephotons;
  slc_n_primary_razzle_other[slcCounter] = nprimrazzleother;
  slc_n_primary_razzled_electrons[slcCounter] = nprimrazzledelectrons;
  slc_n_primary_razzled_muons[slcCounter] = nprimrazzledmuons;
  slc_n_primary_razzled_photons[slcCounter] = nprimrazzledphotons;
  slc_n_primary_razzled_pions[slcCounter] = nprimrazzledpions;
  slc_n_primary_razzled_pions_thresh[slcCounter] = nprimrazzledpionsthresh;
  slc_n_primary_razzled_protons[slcCounter] = nprimrazzledprotons;
  slc_n_primary_razzled_protons_thresh[slcCounter] = nprimrazzledprotonsthresh;

  float totalShowerE = 0;
  float totalTrackE = 0;

  if(slc_n_shws[slcCounter] > 0){
    for (int pfpCounter=0; pfpCounter < (int)slc_n_pfps[slcCounter]; pfpCounter++){
      if((slc_pfp_shower_energy[slcCounter][pfpCounter] > 0) & (slc_pfp_shower_energy[slcCounter][pfpCounter] < 1e9)) totalShowerE += slc_pfp_shower_energy[slcCounter][pfpCounter];
    } 
  }

  if(slc_n_trks[slcCounter] > 0){
    for (int pfpCounter=0; pfpCounter < (int)slc_n_pfps[slcCounter]; pfpCounter++){
      if((slc_pfp_track_ke[slcCounter][pfpCounter] > 0) & (slc_pfp_track_ke[slcCounter][pfpCounter] < 1e9)) totalTrackE += slc_pfp_track_ke[slcCounter][pfpCounter];
    }
  }

  slc_total_shower_E[slcCounter] = totalShowerE;
  slc_total_track_E[slcCounter] = totalTrackE;


}

void sbnd::HNLPiZeroAnalysis::ExtractRazzled(const art::Ptr<sbn::MVAPID> &razzled, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = razzled->mvaScoreMap;

  slc_pfp_razzled_electron_score[slcCounter][pfpCounter] = map.at(11);
  slc_pfp_razzled_muon_score[slcCounter][pfpCounter] = map.at(13);
  slc_pfp_razzled_photon_score[slcCounter][pfpCounter] = map.at(22);
  slc_pfp_razzled_pion_score[slcCounter][pfpCounter] = map.at(211);
  slc_pfp_razzled_proton_score[slcCounter][pfpCounter] = map.at(2212);
  slc_pfp_razzled_pdg[slcCounter][pfpCounter] = razzled->BestPDG();
}


float sbnd::HNLPiZeroAnalysis::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID){

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (fHitsMap[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(fHitsMap[trackID]);
}

float sbnd::HNLPiZeroAnalysis::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

void sbnd::HNLPiZeroAnalysis::AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, 
					  const int slcCounter, const int pfpCounter,
					  const art::Handle<std::vector<recob::Track>> &trackHandle)
{
  //if(fDebug)
  //  std::cout << "Analysing track..." << std::endl;

  art::FindOneP<sbn::MVAPID> tracksToDazzle(trackHandle, e, fDazzleModuleLabel);
  art::FindManyP<anab::Calorimetry> tracksToCalos(trackHandle, e, fCaloModuleLabel);
  art::FindOneP<recob::MCSFitResult> tracksToMCSs(trackHandle, e, fMCSModuleLabel);
  art::FindManyP<anab::ParticleID> tracksToChi2s(trackHandle, e, fChi2ModuleLabel);
  art::FindOneP<sbn::RangeP> tracksToRangePs(trackHandle, e, fRangeModuleLabel);
  art::FindOneP<sbn::ScatterClosestApproach> tracksToClosestApproaches(trackHandle, e, fClosestApproachModuleLabel);
  art::FindOneP<sbn::StoppingChi2Fit> tracksToStoppingChi2s(trackHandle, e, fStoppingChi2ModuleLabel);

  geo::Point_t start = track->Start();
  slc_pfp_track_start_x[slcCounter][pfpCounter] = start.X();
  slc_pfp_track_start_y[slcCounter][pfpCounter] = start.Y();
  slc_pfp_track_start_z[slcCounter][pfpCounter] = start.Z();

  geo::Point_t end = track->End();
  slc_pfp_track_end_x[slcCounter][pfpCounter] = end.X();
  slc_pfp_track_end_y[slcCounter][pfpCounter] = end.Y();
  slc_pfp_track_end_z[slcCounter][pfpCounter] = end.Z();
 
  geo::Vector_t dir = track->StartDirection();
  slc_pfp_track_dir_x[slcCounter][pfpCounter] = dir.X();
  slc_pfp_track_dir_y[slcCounter][pfpCounter] = dir.Y();
  slc_pfp_track_dir_z[slcCounter][pfpCounter] = dir.Z();

  bool isContained = false;
  if( VolumeCheck(start, 5.0, 5.0, 5.0, 5.0) & VolumeCheck(end, 5.0, 5.0, 5.0, 5.0) ) isContained = true;
  slc_pfp_track_contained[slcCounter][pfpCounter] = isContained;

  slc_pfp_track_theta[slcCounter][pfpCounter] = track->Theta() * 180 / TMath::Pi();
  slc_pfp_track_phi[slcCounter][pfpCounter] = track->Phi() * 180 / TMath::Pi();
  slc_pfp_track_length[slcCounter][pfpCounter] = track->Length();

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
    slc_pfp_track_range_p[slcCounter][pfpCounter] = rangeP->range_p;

  const art::Ptr<sbn::ScatterClosestApproach> closestApproach = tracksToClosestApproaches.at(track.key());
  if(closestApproach.isNonnull())
    slc_pfp_track_closest_approach_mean_dca[slcCounter][pfpCounter] = closestApproach->mean;

  const art::Ptr<sbn::StoppingChi2Fit> stoppingChi2 = tracksToStoppingChi2s.at(track.key());
  if(stoppingChi2.isNonnull())
    ExtractStoppingChi2(stoppingChi2, slcCounter, pfpCounter);
}

void sbnd::HNLPiZeroAnalysis::ExtractDazzle(const art::Ptr<sbn::MVAPID> &dazzle, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = dazzle->mvaScoreMap;

  slc_pfp_track_dazzle_muon_score[slcCounter][pfpCounter] = map.at(13);
  slc_pfp_track_dazzle_pion_score[slcCounter][pfpCounter] = map.at(211);
  slc_pfp_track_dazzle_proton_score[slcCounter][pfpCounter] = map.at(211);
  slc_pfp_track_dazzle_other_score[slcCounter][pfpCounter] = map.at(0);
  slc_pfp_track_dazzle_pdg[slcCounter][pfpCounter] = dazzle->BestPDG();
}

void sbnd::HNLPiZeroAnalysis::ExtractCalo(const art::Ptr<anab::Calorimetry> &calo, const int slcCounter, const int pfpCounter)
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

  slc_pfp_track_ke[slcCounter][pfpCounter] = ke;
  slc_pfp_track_charge[slcCounter][pfpCounter] = charge;
}

void sbnd::HNLPiZeroAnalysis::ExtractChi2PID(const art::Ptr<anab::ParticleID> &chi2pid, const int slcCounter, const int pfpCounter)
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
            slc_pfp_track_chi2_muon[slcCounter][pfpCounter] = AlgScore.fValue;
            break;
          case 211:
            slc_pfp_track_chi2_pion[slcCounter][pfpCounter] = AlgScore.fValue;
            break;
          case 321:
            slc_pfp_track_chi2_kaon[slcCounter][pfpCounter] = AlgScore.fValue;
            break;
          case 2212:
            slc_pfp_track_chi2_proton[slcCounter][pfpCounter] = AlgScore.fValue;
            break;
        }
      }
  }

  if(chi2s.size() > 0)
  {
    std::sort(chi2s.begin(), chi2s.end(),
              [](auto const& a, auto const& b)
              { return a.second < b.second; });

    slc_pfp_track_chi2_pdg[slcCounter][pfpCounter] = chi2s[0].first;
  }
}

void sbnd::HNLPiZeroAnalysis::ExtractMCS(const art::Ptr<recob::MCSFitResult> &mcs, const int slcCounter, const int pfpCounter)
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

  slc_pfp_track_mcs_mean_scatter[slcCounter][pfpCounter] = sumScatter / counter;
  slc_pfp_track_mcs_max_scatter_ratio[slcCounter][pfpCounter] = maxScatter / counter;
}

void sbnd::HNLPiZeroAnalysis::ExtractStoppingChi2(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2, const int slcCounter, const int pfpCounter)
{
  const float pol0Chi2 = stoppingChi2->pol0Chi2;
  const float expChi2  = stoppingChi2->expChi2;
  const float ratio    = (pol0Chi2 > 0.f && expChi2 > 0.f) ? pol0Chi2 / expChi2 : -5.f;

  slc_pfp_track_stopping_dedx_chi2_ratio[slcCounter][pfpCounter] = ratio;
  slc_pfp_track_stopping_dedx_pol0_fit[slcCounter][pfpCounter] = stoppingChi2->pol0Fit;
}

void sbnd::HNLPiZeroAnalysis::AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                                           const art::Handle<std::vector<recob::Shower>> &showerHandle, const art::Ptr<recob::Vertex> &vtx,
                                           const std::vector<art::Ptr<recob::Hit>> &hits)
{
  //if(fDebug)
  //  std::cout << "Analysing shower..." << std::endl;
  
  art::FindOneP<sbn::MVAPID> showersToRazzle(showerHandle, e, fRazzleModuleLabel);
  art::FindOneP<float> showersToCosmicDist(showerHandle, e, fCosmicDistModuleLabel);
  art::FindOneP<sbn::ShowerTrackFit> showersToTrackFit(showerHandle, e, fShowerTrackFitModuleLabel);
  art::FindOneP<sbn::ShowerDensityFit> showersToDensityFit(showerHandle, e, fShowerDensityFitModuleLabel);

  geo::Point_t start(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
  slc_pfp_shower_start_x[slcCounter][pfpCounter] = start.X();
  slc_pfp_shower_start_y[slcCounter][pfpCounter] = start.Y();
  slc_pfp_shower_start_z[slcCounter][pfpCounter] = start.Z();
  
  bool isContained = false;

  if (shower->Direction().Z()>-990 && shower->ShowerStart().Z()>-990 && shower->Length()>0) {
    auto end = shower->ShowerStart() + (shower->Length() * shower->Direction());
    slc_pfp_shower_end_x[slcCounter][pfpCounter] = end.X();
    slc_pfp_shower_end_y[slcCounter][pfpCounter] = end.Y();
    slc_pfp_shower_end_z[slcCounter][pfpCounter] = end.Z();
    
    geo::Point_t start(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    if( VolumeCheck(start, 5.0, 5.0, 5.0, 5.0) & VolumeCheck(end, 5.0, 5.0, 5.0, 5.0) ) isContained = true;
  }

  slc_pfp_shower_contained[slcCounter][pfpCounter] = isContained;

  if (shower->Direction().Z()>-990 && shower->ShowerStart().Z()>-990 && shower->Length()>0) {
    const TVector3 dir(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());

    slc_pfp_shower_theta[slcCounter][pfpCounter] = dir.Theta() * 180 / TMath::Pi();
    slc_pfp_shower_phi[slcCounter][pfpCounter] = dir.Phi() * 180 / TMath::Pi();

  }

  double convGap = vtx.isNonnull() ? (start - vtx->position()).R() : def_double;
  slc_pfp_shower_conv_gap[slcCounter][pfpCounter] = convGap;

  geo::Vector_t dir(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
  slc_pfp_shower_dir_x[slcCounter][pfpCounter] = dir.X();
  slc_pfp_shower_dir_y[slcCounter][pfpCounter] = dir.Y();
  slc_pfp_shower_dir_z[slcCounter][pfpCounter] = dir.Z();

  slc_pfp_shower_length[slcCounter][pfpCounter] = shower->Length();
  slc_pfp_shower_open_angle[slcCounter][pfpCounter] = shower->OpenAngle();

  ExtractCalo(shower, slcCounter, pfpCounter, hits);

  const art::Ptr<sbn::MVAPID> razzle = showersToRazzle.at(shower.key());
  if(razzle.isNonnull())
    ExtractRazzle(razzle, slcCounter, pfpCounter);

  const art::Ptr<float> cosmicDist = showersToCosmicDist.at(shower.key());
  if(cosmicDist.isNonnull())
    slc_pfp_shower_cosmic_dist[slcCounter][pfpCounter] = *cosmicDist;

  const art::Ptr<sbn::ShowerTrackFit> trackFit = showersToTrackFit.at(shower.key());
  if(trackFit.isNonnull())
  {
    slc_pfp_shower_track_length[slcCounter][pfpCounter] = trackFit->mTrackLength;
    slc_pfp_shower_track_width[slcCounter][pfpCounter] = trackFit->mTrackWidth;
  }

  const art::Ptr<sbn::ShowerDensityFit> densityFit = showersToDensityFit.at(shower.key());
  if(densityFit.isNonnull())
  {
    slc_pfp_shower_density_grad[slcCounter][pfpCounter] = densityFit->mDensityGrad;
    slc_pfp_shower_density_pow[slcCounter][pfpCounter] = densityFit->mDensityPow;
  }
}

void sbnd::HNLPiZeroAnalysis::ExtractRazzle(const art::Ptr<sbn::MVAPID> &razzle, const int slcCounter, const int pfpCounter)
{
  const std::map<int, float> map = razzle->mvaScoreMap;

  slc_pfp_shower_razzle_electron_score[slcCounter][pfpCounter] = map.at(11);
  slc_pfp_shower_razzle_photon_score[slcCounter][pfpCounter] = map.at(22);
  slc_pfp_shower_razzle_other_score[slcCounter][pfpCounter] = map.at(0);
  slc_pfp_shower_razzle_pdg[slcCounter][pfpCounter] = razzle->BestPDG();
}

void sbnd::HNLPiZeroAnalysis::ExtractCalo(const art::Ptr<recob::Shower> &shower, 
					const int slcCounter, const int pfpCounter,
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

  slc_pfp_shower_energy[slcCounter][pfpCounter] = shower->Energy()[bestPlane];
  slc_pfp_shower_dedx[slcCounter][pfpCounter] = shower->dEdx()[bestPlane];

  const double length      = shower->Length();
  const double bestEnergy  = shower->Energy()[bestPlane];
  const int bestPlaneHits  = showerPlaneHits[bestPlane];
  const double bestPitch   = showerPlanePitches[bestPlane];
  const double wiresHit    = bestPitch > std::numeric_limits<double>::epsilon() ? length / bestPitch : -5.;

  slc_pfp_shower_sqrt_energy_density[slcCounter][pfpCounter] = (length > 0 && bestEnergy > 0) ? std::sqrt(bestEnergy) / length : -5.;
  slc_pfp_shower_modified_hit_density[slcCounter][pfpCounter] = wiresHit > 1. ? bestPlaneHits / wiresHit : -5.;
}

void sbnd::HNLPiZeroAnalysis::ExtractCNNScores(const sbn::PFPCNNScore *cnnscore, const int slcCounter, const int pfpCounter)
{
  slc_pfp_cnnscore_track[slcCounter][pfpCounter] = cnnscore->pfpTrackScore;
  slc_pfp_cnnscore_shower[slcCounter][pfpCounter] = cnnscore->pfpShowerScore;
  slc_pfp_cnnscore_noise[slcCounter][pfpCounter] = cnnscore->pfpNoiseScore;
  slc_pfp_cnnscore_michel[slcCounter][pfpCounter] = cnnscore->pfpMichelScore;
  slc_pfp_cnnscore_endmichel[slcCounter][pfpCounter] = cnnscore->pfpEndMichelScore;
  slc_pfp_cnnscore_nclusters[slcCounter][pfpCounter] = cnnscore->nClusters;
}

void sbnd::HNLPiZeroAnalysis::AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc, 
						const int slcCounter,
						const art::Handle<std::vector<recob::Slice>> &sliceHandle)
{
  //if(fDebug)
  //  std::cout << "Analysing Slice Truth..." << std::endl;

  art::FindManyP<recob::Hit> slicesToHits(sliceHandle, e, fSliceModuleLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const std::vector<art::Ptr<recob::Hit>> sliceHits = slicesToHits.at(slc.key());

  //Get map of Track Id - nhits matched to slice
  std::map<int, int> objectHitMap;
  for(auto const &hit : sliceHits)
    objectHitMap[TruthMatchUtils::TrueParticleID(clockData, hit, true)]++;

  int maxHits  = def_int;
  int bestTrackID = def_int;

  //Get map of MCTruth - nhits matched to slice
  std::map<const art::Ptr<simb::MCTruth>, int> mcTruthHitMap;
  for(auto const& [trackID, nhits] : objectHitMap)
  {
    const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
    mcTruthHitMap[mct] += nhits;

    if(nhits > maxHits){
      bestTrackID = trackID; 
    }
  }
  
  maxHits  = def_int;
  art::Ptr<simb::MCTruth> bestMCT = art::Ptr<simb::MCTruth>();

  //Get MCTruth that deposits the most hits
  for(auto const& [mct, nhits] : mcTruthHitMap)
  {
    if(nhits > maxHits)
    {
      maxHits = nhits;
      bestMCT  = mct;
    }
  }
  
  const float comp = fMCTruthHitsMap[bestMCT] == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(fMCTruthHitsMap[bestMCT]);
  const float pur  = sliceHits.size() == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(sliceHits.size());

  slc_comp[slcCounter] = comp;
  slc_pur[slcCounter] = pur;

  if(bestMCT.isNonnull())
  {
    if (fGenie){
      if(bestMCT->Origin() == 2){
        auto const mcp = particleInv->TrackIdToParticle_P(bestTrackID);

        slc_true_event_type[slcCounter] = (int) kCosmic;
        slc_true_vtx_x[slcCounter] = mcp->Vx();
        slc_true_vtx_y[slcCounter] = mcp->Vy();
        slc_true_vtx_z[slcCounter] = mcp->Vz();
        slc_true_vtx_t[slcCounter] = mcp->T();
  	slc_true_vtx_t_corrected_Z[slcCounter] = slc_true_vtx_t[slcCounter] - slc_true_vtx_z[slcCounter]/light_speed_cm_ns;
      }else if(bestMCT->Origin() == 0){
        slc_true_event_type[slcCounter] = (int) kUnknownEv;
      }else{
        AnalyseSliceMCTruth(e, bestMCT, slcCounter);
      }
    }
    else if (fMeVPrtl){
      if(bestMCT->Origin() == 2){
	//std::cout << "This slice cosmics " << std::endl;
        auto const mcp = particleInv->TrackIdToParticle_P(bestTrackID);

        slc_true_event_type[slcCounter] = (int) kCosmic;
        slc_true_vtx_x[slcCounter] = mcp->Vx();
        slc_true_vtx_y[slcCounter] = mcp->Vy();
        slc_true_vtx_z[slcCounter] = mcp->Vz();
        slc_true_vtx_t[slcCounter] = mcp->T();
  	slc_true_vtx_t_corrected_Z[slcCounter] = slc_true_vtx_t[slcCounter] - slc_true_vtx_z[slcCounter]/light_speed_cm_ns;
      }else{
	//std::cout << "This slice hnl " << std::endl;
        AnalyseSliceMCTruth(e, bestMCT, slcCounter);
      }

    }
    else if (fBeamOff){
      if(bestMCT->Origin() == 2){
        auto const mcp = particleInv->TrackIdToParticle_P(bestTrackID);

        slc_true_event_type[slcCounter] = (int) kCosmic;
        slc_true_vtx_x[slcCounter] = mcp->Vx();
        slc_true_vtx_y[slcCounter] = mcp->Vy();
        slc_true_vtx_z[slcCounter] = mcp->Vz();
        slc_true_vtx_t[slcCounter] = mcp->T();
  	slc_true_vtx_t_corrected_Z[slcCounter] = slc_true_vtx_t[slcCounter] - slc_true_vtx_z[slcCounter]/light_speed_cm_ns;
      }else{
        slc_true_event_type[slcCounter] = (int) kUnknownEv;
      }
    } 
  } //End of best MCT
}

void sbnd::HNLPiZeroAnalysis::AnalyseSliceMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct, const int slcCounter)
{
  if(fDebug){
    std::cout << "Analysing Slice MCTruth ..." << std::endl;
  }

  const simb::MCNeutrino mcn = mct->GetNeutrino();
  const simb::MCParticle nu  = mcn.Nu();
  
  const bool nc = mcn.CCNC() == 1;
  const bool av = VolumeCheck(nu.Position().Vect());
  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);
 
  // Only save flux weights for now
  ResizeWeight2DVectorCol(slcCounter); //Vector is initialised for 1000 columns

  for(auto const& weightModuleLabel : fEventWeightModuleLabels)
  {
    art::FindManyP<sbn::evwgh::EventWeightMap> MCTruthToWeights( { mct }, e, weightModuleLabel);
    const std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> ewms = MCTruthToWeights.at(0);

    for(auto const& ewm : ewms)
    {
      
      for(auto const& [ name, weights ] : *ewm)
      {
	// Multi-Sigma Cross Section
	if (name.find("multisigma") != std::string::npos){
	  //uni-sim only has one varible
	  if (weights.size() == 1){
	    for(auto const w: weights){
	      if(name.find("DecayAngMEC") != std::string::npos )
	  	    slc_xsec_unisim_DecayAngMEC[slcCounter][0] = w;
	      if(name.find("ThetaDelta2NRad") != std::string::npos )
	  	    slc_xsec_unisim_ThetaDelta2NRad[slcCounter][0] = w;
              if(name.find("Theta_Delta2Npi") != std::string::npos )
	  	    slc_xsec_unisim_Theta_Delta2Npi[slcCounter][0] = w;
              if(name.find("VecFFCCQEshape") != std::string::npos )
	  	    slc_xsec_unisim_VecFFCCQEshape[slcCounter][0] = w;
	    }
	  }
	  //multi-sigma has 6 variables  -1, +1, -2, +2, -3, +3 sigma
	  else if (weights.size() == 6){
	    unsigned int wCounter = 0;  //only need to save index 0 and 1 for -/+1 sigma    
	    for(auto const w: weights){
 	      if(name.find("CoulombCCQE") != std::string::npos )
		      slc_xsec_multisigma_CoulombCCQE[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarnCC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarnCC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarnCC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarnCC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarnNC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarnNC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarnNC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarnNC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarpCC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarpCC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarpCC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarpCC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarpNC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarpNC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvbarpNC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvbarpNC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvnCC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvnCC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvnCC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvnCC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvnNC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvnNC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvnNC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvnNC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvpCC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvpCC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvpCC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvpCC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvpNC1pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvpNC1pi[slcCounter][wCounter] = w;
 	      if(name.find("NonRESBGvpNC2pi") != std::string::npos )
		      slc_xsec_multisigma_NonRESBGvpNC2pi[slcCounter][wCounter] = w;
 	      if(name.find("NormCCMEC") != std::string::npos )
		      slc_xsec_multisigma_NormCCMEC[slcCounter][wCounter] = w;
 	      if(name.find("NormNCMEC") != std::string::npos )
		      slc_xsec_multisigma_NormNCMEC[slcCounter][wCounter] = w;
 	      if(name.find("RDecBR1eta") != std::string::npos )
		      slc_xsec_multisigma_RDecBR1eta[slcCounter][wCounter] = w;
 	      if(name.find("RDecBR1gamma") != std::string::npos )
		      slc_xsec_multisigma_RDecBR1gamma[slcCounter][wCounter] = w;
 	      if(name.find("RPA_CCQE") != std::string::npos )
		      slc_xsec_multisigma_RPA_CCQE[slcCounter][wCounter] = w;
	      if(name.find("NormCCCOH") != std::string::npos )
		      slc_xsec_multisigma_NormCCCOH[slcCounter][wCounter] = w;
	      if(name.find("NormNCCOH") != std::string::npos )
		      slc_xsec_multisigma_NormNCCOH[slcCounter][wCounter] = w;

	      wCounter++;
	    }
	  }
	  
	}
	//Multi-Sim xsection weight
	else if (name.find("multisim") != std::string::npos){
	  unsigned int wCounter = 0; //Multisim xsection size = 500
	  for (auto const w: weights){
	    if (name.find("ZExpA") != std::string::npos )  slc_xsec_multisim_ZExpA[slcCounter][wCounter] = w;
	    if (name.find("NCEL") != std::string::npos )  slc_xsec_multisim_NCEL[slcCounter][wCounter] = w;
	    if (name.find("CCRES") != std::string::npos )  slc_xsec_multisim_CCRES[slcCounter][wCounter] = w;
	    if (name.find("NCRES") != std::string::npos )  slc_xsec_multisim_NCRES[slcCounter][wCounter] = w;
	    if (name.find("DISBY") != std::string::npos )  slc_xsec_multisim_DISBY[slcCounter][wCounter] = w;
	    if (name.find("FSI_pi") != std::string::npos )  slc_xsec_multisim_FSI_pi[slcCounter][wCounter] = w;
	    if (name.find("FSI_N") != std::string::npos )  slc_xsec_multisim_FSI_N[slcCounter][wCounter] = w;

            slc_xsec_multisim_total[slcCounter][wCounter] *= w;
	    wCounter++;
	  }
	}

	// Multi-Sim Flux Weights
	else if (name.find("_Flux") != std::string::npos){
	  unsigned int wCounter = 0;  //Flux weight size = 1000    
	  for(auto const w: weights){
	    if (name == "expskin_Flux") slc_flux_weight_expskin[slcCounter][wCounter] = w;
	    if (name == "horncurrent_Flux") slc_flux_weight_horncurrent[slcCounter][wCounter] = w;
	    if (name == "kminus_Flux") slc_flux_weight_kminus[slcCounter][wCounter] = w;
	    if (name == "kplus_Flux") slc_flux_weight_kplus[slcCounter][wCounter] = w;
	    if (name == "kzero_Flux") slc_flux_weight_kzero[slcCounter][wCounter] = w;
	    if (name == "nucleoninexsec_Flux") slc_flux_weight_nucleoninexsec[slcCounter][wCounter] = w;
	    if (name == "nucleonqexsec_Flux") slc_flux_weight_nucleonqexsec[slcCounter][wCounter] = w;
	    if (name == "nucleontotxsec_Flux") slc_flux_weight_nucleontotxsec[slcCounter][wCounter] = w;
	    if (name == "piminus_Flux") slc_flux_weight_piminus[slcCounter][wCounter] = w;
	    if (name == "pioninexsec_Flux") slc_flux_weight_pioninexsex[slcCounter][wCounter] = w;
	    if (name == "pionqexsec_Flux") slc_flux_weight_pionqexsec[slcCounter][wCounter] = w;
	    if (name == "piontotxsec_Flux") slc_flux_weight_piontotxsec[slcCounter][wCounter] = w;
	    if (name == "piplus_Flux") slc_flux_weight_piplus[slcCounter][wCounter] = w;
	    
	    slc_flux_weight_total[slcCounter][wCounter] *= w;
	    wCounter++;
          }
	}
	// Geant4 weight
	else if (name.find("Geant4") != std::string::npos){
	  unsigned int wCounter = 0; //Multisim geant4 size = 1000
	  for (auto const w: weights){
            slc_geant4_multisim_reinteractions[slcCounter][wCounter] = w;
	    wCounter++;
	  }
	}
      } //Loop over weight names
    } 
  } //Loop over weight labels

  art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);
  
  int neutral_pions = 0;
  float trueEnDep = 0.;
  
  for(auto const& mcp : MCParticleVec)
  {
    if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
    {
      switch(abs(mcp->PdgCode()))
      { 
        case 111:
          ++neutral_pions;
          break;
      }
    }
    std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());
    for(auto const& ide : ides)
      trueEnDep += ide->energy / 1000.;
  }
  
  const bool pizero = neutral_pions > 0;
  const bool _is_hnl = mct->Origin() == 0 && _n_hnl > 0; //MeVPrtl check always run before truth matching so n_hnl == 1 if running on HNL sample
  const bool _is_nu = mct->Origin() == 1? 1:0;
  
  // Event Classification
  if (_is_hnl && fv)
  {
    slc_true_event_type[slcCounter] = (int) kHNL;
  }else
  {
    if(_is_hnl && !fv && av)
      slc_true_event_type[slcCounter] = (int) kHNLNonFV;
    else if(_is_hnl && !av)
      slc_true_event_type[slcCounter] = (int) kHNLDirt;
    else if(_is_nu && nc && pizero && fv)
      slc_true_event_type[slcCounter] = (int) kNCPiZero;
    else if(_is_nu && nc && !pizero && fv)
      slc_true_event_type[slcCounter] = (int) kOtherNC;
    else if(_is_nu && abs(nu.PdgCode()) == 14 && fv)
      slc_true_event_type[slcCounter] = (int) kCCNuMu;
    else if(_is_nu && abs(nu.PdgCode()) == 12 && fv)
      slc_true_event_type[slcCounter] = (int) kCCNuE;
    else if(_is_nu && !fv && av)
      slc_true_event_type[slcCounter] = (int) kNuNonFV;
    else if(_is_nu && !av)
      slc_true_event_type[slcCounter] = (int) kNuDirt;
    else
      slc_true_event_type[slcCounter] = (int) kUnknownEv;
  }

  if (fDebug) std::cout << "This is event type = " << slc_true_event_type[slcCounter] << std::endl;

  slc_true_mctruth_id[slcCounter] = mct.key();
  slc_true_en_dep[slcCounter] = trueEnDep;
  slc_true_vtx_x[slcCounter] = nu.Vx(); //cm
  slc_true_vtx_y[slcCounter] = nu.Vy(); //cm
  slc_true_vtx_z[slcCounter] = nu.Vz(); //cm 
  slc_true_vtx_t[slcCounter] = nu.T(); //ns ?
  slc_true_vtx_t_corrected_Z[slcCounter] = slc_true_vtx_t[slcCounter] - slc_true_vtx_z[slcCounter]/light_speed_cm_ns;
  
  if (fPiZeroAna){
    
    float pizeroMom = -99999;
    float pizeroCosTheta = -99999;

    for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
      {
        if(abs(mcp->PdgCode()) == 111)
        { 
           TVector3 pizeroDir(mcp->Px(), mcp->Py(), mcp->Pz()); 
           pizeroMom = pizeroDir.Mag() * 1000;
           pizeroCosTheta = pizeroDir.Z() / pizeroMom;
        }
      }

    }

    slc_true_pizero_mom[slcCounter] = pizeroMom;
    slc_true_pizero_cosThetaZ[slcCounter] = pizeroCosTheta;

    //if(fDebug){
    //  std::cout << "====================================" << std::endl;
    //  std::cout << "Event " << _event << " Slice " << slcCounter << " has " << slc_n_shws[slcCounter] << " showers" << std::endl;
    //  std::cout << "mom = " << pizeroMom << std::endl;
    //  std::cout << "CosThetaZ = " << pizeroCosTheta << std::endl;
    //  std::cout << "====================================" << std::endl;
    //} 
  }
}

void sbnd::HNLPiZeroAnalysis::AnalysePiZero(const int slcCounter){
  
  if(slc_is_clear_cosmics[slcCounter] ==1) return;

  //Need >= 1 showers
  if (slc_n_shws[slcCounter] < 0 ) return;

  float E_best = -999;
  int pfp0 = -999;

  //Chose best photon candidate for the slice: most deposited energy for now
  for (unsigned int pfpCounter = 0; pfpCounter < slc_pfp_shower_energy[slcCounter].size(); pfpCounter++){
    
    //Check if shower 1 is good
    if ( (slc_pfp_shower_energy[slcCounter][pfpCounter] < 0 )
	 || (slc_pfp_shower_start_z[slcCounter][pfpCounter] < -990)
	 || (slc_pfp_shower_dir_z[slcCounter][pfpCounter] < -990)
	 || (slc_pfp_shower_length[slcCounter][pfpCounter] < 0)
	 || (slc_pfp_shower_dedx[slcCounter][pfpCounter] < -990) ){
      continue;
    }
    
    if(slc_pfp_shower_energy[slcCounter][pfpCounter] > E_best){
      E_best = slc_pfp_shower_energy[slcCounter][pfpCounter];
      pfp0 = pfpCounter;
    }
  }

  //Need >= 2 showers
  if (slc_n_shws[slcCounter] < 1 ) return;

  float diff = 99999;
  float m_best = -999;
  int pfp1 = -999;
  int pfp2 = -999;

  //Chose best 2 photons candidate for the slice
  for (unsigned int pfpCounter1 = 0; pfpCounter1 < slc_pfp_shower_energy[slcCounter].size(); pfpCounter1++){
    
    //Check if shower 1 is good
    if ( (slc_pfp_shower_energy[slcCounter][pfpCounter1] < 0 )
	 || (slc_pfp_shower_start_z[slcCounter][pfpCounter1] < -990)
	 || (slc_pfp_shower_dir_z[slcCounter][pfpCounter1] < -990)
	 || (slc_pfp_shower_length[slcCounter][pfpCounter1] < 0)
	 || (slc_pfp_shower_dedx[slcCounter][pfpCounter1] < -990)){

      continue;
    }

    for (unsigned int pfpCounter2 = 0; pfpCounter2 < slc_pfp_shower_energy[slcCounter].size(); pfpCounter2++){

      //Need to look at a different pfp
      if(pfpCounter2 == pfpCounter1) continue; 				

      //Check if shower 2 is good
      if ( (slc_pfp_shower_energy[slcCounter][pfpCounter2] < 0 )
           || (slc_pfp_shower_start_z[slcCounter][pfpCounter2] < -990)
           || (slc_pfp_shower_dir_z[slcCounter][pfpCounter2] < -990)
           || (slc_pfp_shower_length[slcCounter][pfpCounter2] < 0)
	   || (slc_pfp_shower_dedx[slcCounter][pfpCounter2] < -990)){

        continue;
      }

      float E1 = slc_pfp_shower_energy[slcCounter][pfpCounter1];
      float E2 = slc_pfp_shower_energy[slcCounter][pfpCounter2];

      const TVector3 dir1(slc_pfp_shower_dir_x[slcCounter][pfpCounter1], slc_pfp_shower_dir_y[slcCounter][pfpCounter1], slc_pfp_shower_dir_z[slcCounter][pfpCounter1]);
      const TVector3 dir2(slc_pfp_shower_dir_x[slcCounter][pfpCounter2], slc_pfp_shower_dir_y[slcCounter][pfpCounter2], slc_pfp_shower_dir_z[slcCounter][pfpCounter2]);

      float cosTheta = dir1.Dot(dir2) / ( dir1.Mag() * dir2.Mag() );

      float m_inv = sqrt (2 * E1 * E2 * (1 - cosTheta ));
    
      if (TMath::Abs( m_inv - pi0_mass) < diff){
        diff = TMath::Abs( m_inv - pi0_mass);
	m_best = m_inv;
	pfp1 = pfpCounter1;
	pfp2 = pfpCounter2;
      }
    }
  }
 
  bool goodPiZero = false;

  if (E_best > -999){
    
    slc_pizero_good[slcCounter] = goodPiZero;
    slc_photon0_pfp_id[slcCounter] = pfp0;
    slc_photon0_energy[slcCounter] = E_best;
    
    //if (fDebug){
    //  std::cout << "====================================" << std::endl;
    //  std::cout << "Event " << _event << " Slice " << slcCounter << "  One Shower Candidate" << std::endl;
    //  std::cout << "------------------------------------" << std::endl;
    //  std::cout << "pfp id " << slc_pfp_id[slcCounter][pfp0] << std::endl;
    //  std::cout << "razzled photon score " << slc_pfp_razzled_photon_score[slcCounter][pfp0] << std::endl;
    //  std::cout << "shower energy " << slc_pfp_shower_energy[slcCounter][pfp0] << std::endl;
    //  std::cout << "dEdx " << slc_pfp_shower_dedx[slcCounter][pfp0] << std::endl;
    //  std::cout << "====================================" << std::endl;
    //}
  }
  else if ( m_best > 0 ) {
 
    goodPiZero = true;

    const float E1 = slc_pfp_shower_energy[slcCounter][pfp1];
    const float E2 = slc_pfp_shower_energy[slcCounter][pfp2];
  
    const TVector3 dir1(slc_pfp_shower_dir_x[slcCounter][pfp1], slc_pfp_shower_dir_y[slcCounter][pfp1], slc_pfp_shower_dir_z[slcCounter][pfp1]);
    const TVector3 dir2(slc_pfp_shower_dir_x[slcCounter][pfp2], slc_pfp_shower_dir_y[slcCounter][pfp2], slc_pfp_shower_dir_z[slcCounter][pfp2]);
  
    const float cosTheta = dir1.Dot(dir2) / ( dir1.Mag() * dir2.Mag() );
  
    const float pizeroMass = sqrt (2 * E1 * E2 * (1 - cosTheta ));
  
    const TVector3 pizeroDir = (E1 * dir1) + (E2 * dir2);

    const float pizeroMom = pizeroDir.Mag();
    const float pizeroCosTheta = pizeroDir.Z() / pizeroMom;
 
    if (fDebug){
      std::cout << "====================================" << std::endl;
      std::cout << "Event " << _event << " Slice " << slcCounter << "  Two Shower Candidate" << std::endl;
      std::cout << "------------------------------------" << std::endl;
      std::cout << "pfp id " << slc_pfp_id[slcCounter][pfp1] << std::endl;
      std::cout << "razzled photon score " << slc_pfp_razzled_photon_score[slcCounter][pfp1] << std::endl;
      std::cout << "shower energy " << slc_pfp_shower_energy[slcCounter][pfp1] << std::endl;
      std::cout << "dEdx " << slc_pfp_shower_dedx[slcCounter][pfp1] << std::endl;
      std::cout << "------------------------------------" << std::endl;
      std::cout << "pfp id " << slc_pfp_id[slcCounter][pfp2] << std::endl;
      std::cout << "razzled photon score " << slc_pfp_razzled_photon_score[slcCounter][pfp2] << std::endl;
      std::cout << "shower energy " << slc_pfp_shower_energy[slcCounter][pfp2] << std::endl;
      std::cout << "dEdx " << slc_pfp_shower_dedx[slcCounter][pfp2] << std::endl;
      std::cout << "------------------------------------" << std::endl;
      std::cout << "OpenAngle = " << TMath::ACos(cosTheta) * 180 / TMath::Pi() << std::endl;
      std::cout << "mass = " << pizeroMass << std::endl;
      std::cout << "mom = " << pizeroMom << std::endl;
      std::cout << "CosThetaZ = " << pizeroCosTheta << std::endl;
      std::cout << "GoodPizero = " << goodPiZero << std::endl;
      std::cout << "====================================" << std::endl;
      std::cout << std::endl;
    }

    slc_pizero_good[slcCounter] = goodPiZero;

    slc_photon1_pfp_id[slcCounter] = pfp1;
    slc_photon2_pfp_id[slcCounter] = pfp2;
  
    slc_photon1_energy[slcCounter] = E1;
    slc_photon2_energy[slcCounter] = E2;
   
    slc_pizero_open_angle[slcCounter] = TMath::ACos(cosTheta) * 180 / TMath::Pi();
    slc_pizero_mass[slcCounter] = pizeroMass;
    slc_pizero_mom[slcCounter] = pizeroMom;
    slc_pizero_CosThetaZ[slcCounter] = pizeroCosTheta;
  }
  else{
    slc_pizero_good[slcCounter] = goodPiZero;
  }
}

DEFINE_ART_MODULE(sbnd::HNLPiZeroAnalysis)
