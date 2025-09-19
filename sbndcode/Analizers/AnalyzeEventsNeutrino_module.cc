////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEventsNeutrino
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEventsNeutrino_module.cc
//
// Generated at Thu Jan 26 04:44:06 2023 by Luis Pelegrina gutierrez using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////
// Additional Framework includes
#include "art_root_io/TFileService.h"

// Additional LArSoft includes
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/MCBase/MCTrack.h"

// Additional LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"


#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/Common/Reco/Stub.h"


// ROOT includes                                                                                                                                                                        
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <larcorealg/Geometry/Exceptions.h>

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

//Truth Matching
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

constexpr int def_int = std::numeric_limits<int>::min();

namespace test {
  class AnalyzeEventsNeutrino;
}


class test::AnalyzeEventsNeutrino : public art::EDAnalyzer {
public:
  explicit AnalyzeEventsNeutrino(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEventsNeutrino(AnalyzeEventsNeutrino const &) = delete;

  AnalyzeEventsNeutrino(AnalyzeEventsNeutrino &&) = delete;

  AnalyzeEventsNeutrino &operator=(AnalyzeEventsNeutrino const &) = delete;

  AnalyzeEventsNeutrino &operator=(AnalyzeEventsNeutrino &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void respondToOpenInputFile(art::FileBlock const &);
  void beginSubRun(const art::SubRun &sr);
  void endSubRun(const art::SubRun &sr);


  void reset_g4_vars();
  void reset_gen_vars();
  void reset_slc_vars();

  int get_total_gen_evts(const art::Event &e);
  std::vector<int> get_point_in_wires(double x, double y, double z);
  int vertex_to_drift_tick(double vt, double vx);

  void do_hit_matching(std::vector<art::Ptr<recob::Hit>> hit_vector, int pfp_ID);
  //void do_sp_matching(std::vector<art::Ptr<recob::SpacePoint>> sp_vector, int pfp_ID);

  void analize_truth_gen(const art::Ptr< simb::MCTruth> MCT, bool is_non_null, int best_track_ID);
  void analize_truth_g4_particles(std::vector<art::Ptr < simb::MCParticle>> MCP_vec);
  void analize_slc_pfps(art::Event const &e, art::Handle< std::vector<recob::PFParticle> > pfp_handle, std::vector< art::Ptr<recob::PFParticle> > pfp_vec, art::Handle< std::vector<recob::Track> > track_handle, art::Handle< std::vector<recob::Shower>> shower_handle, art::FindManyP<recob::Hit> sp_hit_assns);

  bool slc_check_clear_cosmic(std::vector<art::Ptr < recob::PFParticle>> pfp_vec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns);
  double get_slice_nu_score(std::vector<art::Ptr<recob::PFParticle>> pfp_vec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns);
  bool get_slice_nu_vertex(art::Event const &e,std::vector<art::Ptr<recob::PFParticle>> slice_pfp_vec, art::FindManyP<recob::Vertex> pfp_vertex_assns);

  void do_truth_matching(const detinfo::DetectorClocksData& clock_data_evt, const std::vector<art::Ptr<recob::Hit> >& track_hit_vec);
  int NumberofHitsFromTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits);

  void save_slc_trk_calo(std::vector<art::Ptr<anab::Calorimetry>> track_calo);
  void save_slc_pfp_vtx(std::vector<art::Ptr < recob::Vertex>> pfp_vertex_vec);
  void save_slc_pfp_razzled(std::vector<art::Ptr<sbn::MVAPID> >  pfp_razzled_vec);

  //int TrueParticleID(const detinfo::DetectorClocksData& clockData, const art::Ptr<recob::Hit>& hit);
  //int TrueParticleIDFromTotalRecoCharge(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits);
  std::map<int,std::map<geo::PlaneID,int> > NumberofPlaneHitsPerTrack(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

  void save_track_info(std::vector<art::Ptr<recob::Track>> pfp_track_vec, std::vector<art::Ptr<recob::Shower>> pfp_shower_vec, art::FindManyP<anab::Calorimetry> track_calo_assns, art::FindManyP<anab::ParticleID> track_chi2_assn, art::FindManyP<recob::MCSFitResult> tracksToMCSs_muon, art::FindManyP<recob::MCSFitResult> tracksToMCSs_pion , art::FindManyP<sbn::RangeP> tracksToRangePs_muon, art::FindManyP<sbn::RangeP> tracksToRangePs_pion, art::FindManyP<sbn::ScatterClosestApproach> tracksToClosestApproaches, art::FindManyP<sbn::StoppingChi2Fit> tracksToStoppingChi2s);
  void save_shower_info(std::vector<art::Ptr<recob::Track>> pfp_track_vec, std::vector<art::Ptr<recob::Shower>> pfp_shower_vec);
  void save_pfp_track_start_end(TVector3 start, TVector3 end);

  void FillMomDiff(const art::Ptr<sbn::RangeP> &rangeP_muon, const art::Ptr<sbn::RangeP> &rangeP_pion, const art::Ptr<recob::MCSFitResult> &mcs_muon, const art::Ptr<recob::MCSFitResult> &mcs_pion);
  void FillStoppingChi2Metrics(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2);
  void FillMCSMetrics(const art::Ptr<recob::MCSFitResult> &mcs);

  void analyze_slice_truth(art::Event const& e, std::vector<art::Ptr<recob::Hit>> slice_hit_vec,   std::map<const art::Ptr<simb::MCTruth>, int> MCTruth_hit_map);
  void FillPandoraTrackIDScoreVars(const std::map<std::string, float> &propertiesMap);
  void FillPandoraTrackIDScoreVar(const std::map<std::string, float> &propertiesMap, std::vector<double> &var, std::string name);

  void fill_empty_slice();
  void fill_empty_truth();
  void fill_empty_slice_pfp_track_info();
  void fill_empty_slice_pfp_shower_info();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv = art::ServiceHandle<cheat::ParticleInventoryService>();
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();


  std::map<int,std::map<geo::PlaneID,int> > slice_true_hit_map;

  TTree *fSubRunTree;
  double POT;
  double spill;

  // Create out output tree.
  TTree *fTree;
  std::string file_name;

  //Event variables
  unsigned int event_ID;
  unsigned int run_ID;
  unsigned int subrun_ID;
  int num_gen_evts;

  //MCTruth variables
  //Consult https://internal.dunescience.org/doxygen/classsimb_1_1MCNeutrino.html#a8e213917487f1f5cd52ca5280d33a777 for information on the numbers meaning
  int nu_PDG;
  double nu_E0;
  double nu_weight;
  int gen_index;
  unsigned int nu_interaction_mode;
  unsigned int nu_interaction_type;
  int nu_CC_NC;
  int nu_target;
  int nu_HitNuc;
  int nu_HitQuark;
  double nu_W;
  double nu_X;
  double nu_Y;
  double nu_QSqr;

  //MCTruth Particles variables
  std::vector<int> gen_part_trackID;
  std::vector<int> gen_part_statusCode;
  std::vector<int> gen_part_mother;
  std::vector<int> gen_part_PDGcode;
  std::vector<double> gen_part_mass;
  std::vector<double> gen_part_E0;
  std::vector<double> gen_part_start_pos_X, gen_part_start_pos_Y, gen_part_start_pos_Z;
  std::vector<double> gen_part_P0_Y, gen_part_P0_Z, gen_part_P0_X;

  //Geant4 Particles variables
  std::vector<int> g4_part_trackID;
  std::vector<int> g4_part_mother;
  std::vector<int> g4_part_PDGcode;
  std::vector<double> g4_part_mass;
  std::vector<double> g4_part_TL;
  std::vector<double> g4_part_E0;
  std::vector<double> g4_part_Ef;

  std::vector<double> g4_part_start_pos_X, g4_part_start_pos_Y, g4_part_start_pos_Z;
  std::vector<double> g4_part_start_T;
  std::vector<double> g4_part_end_pos_X, g4_part_end_pos_Y, g4_part_end_pos_Z;
  std::vector<double> g4_part_end_T;

  std::vector<double> g4_part_P0_X, g4_part_P0_Y, g4_part_P0_Z;
  std::vector<double> g4_part_Pf_X, g4_part_Pf_Y, g4_part_Pf_Z;

  std::vector<std::string> g4_part_process;
  std::vector<std::string> g4_part_end_process;

  //Slices related parameters
  bool slc_is_reconstructed;
  bool slc_is_obvious_cosmic;
  int slc_ID;

  float slc_true_prim_vtx_x, slc_true_prim_vtx_y, slc_true_prim_vtx_z;
  int slc_true_prim_vtx_wires_U, slc_true_prim_vtx_wires_V, slc_true_prim_vtx_wires_C;
  double slc_true_prim_vtx_time_tick;
  int slc_true_prim_vtx_channel_U, slc_true_prim_vtx_channel_V, slc_true_prim_vtx_channel_C;


  float slc_prim_vtx_x, slc_prim_vtx_y, slc_prim_vtx_z;
  int slc_prim_vtx_wires_U, slc_prim_vtx_wires_V, slc_prim_vtx_wires_C;
  double slc_prim_vtx_time_tick;
  int slc_prim_vtx_channel_U, slc_prim_vtx_channel_V, slc_prim_vtx_channel_C;

  double slc_truth_matching_completeness;
  double slc_truth_matching_purity;
  double slc_crumbs_score;
  double slc_nu_score;

  double slc_t0;
  double slc_t0_score;
  double slc_OpT0Fraction;

  double slc_bcmatch_flash_time;
  double slc_bcmatch_first_hit;
  double slc_bcmatch_flash_PE;
  double slc_bcmatch_flash_asymmetry;


  std::vector<double> slc_stub_vtx_x;
  std::vector<double> slc_stub_vtx_y;
  std::vector<double> slc_stub_vtx_z;
  std::vector<double> slc_stub_end_x;
  std::vector<double> slc_stub_end_y;
  std::vector<double> slc_stub_end_z;
  std::vector<std::vector<double>> slc_stub_hit_charge;
  std::vector<std::vector<double>> slc_stub_hit_wire;
  std::vector<std::vector<int>> slc_stub_hit_ontrack;

  std::vector<double> slc_sp_x;
  std::vector<double> slc_sp_y;
  std::vector<double> slc_sp_z;
  std::vector<double> slc_sp_chi2;
  std::vector<double> slc_sp_integral;
  std::vector<double> slc_sp_sigma_integral;
  std::vector<int> slc_sp_associated_pfp_ID;

  std::vector<long unsigned int> slc_hits_wire_ID;
  std::vector<long unsigned int> slc_hits_plane_ID;
  std::vector<long unsigned int> slc_hits_TPC_ID;
  std::vector<long unsigned int> slc_hits_channel_ID;
  std::vector<double> slc_hits_integral;
  std::vector<double> slc_hits_sigma_integral;
  std::vector<double> slc_hits_peaktime;
  std::vector<double> slc_hits_RMS;
  std::vector<int> slc_hits_associated_pfp_ID;

  std::vector<std::vector<double>> slc_pfp_calo_dEdx;
  std::vector<std::vector<double>> slc_pfp_calo_dQdx;
  std::vector<std::vector<double>> slc_pfp_calo_pitch;
  std::vector<std::vector<double>> slc_pfp_calo_residual_range;
  std::vector<std::vector<double>> slc_pfp_calo_x;
  std::vector<std::vector<double>> slc_pfp_calo_y;
  std::vector<std::vector<double>> slc_pfp_calo_z;

  std::vector<int> slc_pfp_ID;
  std::vector<int> slc_pfp_parent_ID;
  std::vector<bool> slc_pfp_mother_is_primary;
  std::vector<double> slc_pfp_daughter_max_hits;
  std::vector<double> slc_pfp_track_score;
  std::vector<double> slc_pfp_charge_end_frac;
  std::vector<double> slc_pfp_charge_frac_spread;
  std::vector<double> slc_pfp_linear_fit_diff;
  std::vector<double> slc_pfp_linear_fit_length;
  std::vector<double> slc_pfp_linear_fit_gap_length;
  std::vector<double> slc_pfp_linear_fit_RMS;
  std::vector<double> slc_pfp_open_angle_diff;
  std::vector<double> slc_pfp_secondary_PCA_ratio;
  std::vector<double> slc_pfp_tertiary_PCA_ratio;
  std::vector<double> slc_pfp_vertex_dist;
  std::vector<bool> slc_pfp_has_track;
  std::vector<int> slc_pfp_pandora_PDG;

  std::vector<double> slc_pfp_razzled_score_e;
  std::vector<double> slc_pfp_razzled_score_mu;
  std::vector<double> slc_pfp_razzled_score_pho;
  std::vector<double> slc_pfp_razzled_score_pi;
  std::vector<double> slc_pfp_razzled_score_p;
  std::vector<int> slc_pfp_razzled_pdg;

  std::vector<double> slc_pfp_vtx_x, slc_pfp_vtx_y, slc_pfp_vtx_z;
  std::vector<int> slc_pfp_vtx_wires_U, slc_pfp_vtx_wires_V, slc_pfp_vtx_wires_C;
  std::vector<double> slc_pfp_vtx_time_tick;
  std::vector<int> slc_pfp_vtx_channel_U, slc_pfp_vtx_channel_V, slc_pfp_vtx_channel_C;

  std::vector<int> slc_pfp_true_trackid;
  std::vector<int> slc_pfp_true_pdg;
  std::vector<double> slc_pfp_true_energy;
  std::vector<double> slc_pfp_true_TL;
  std::vector<std::string> slc_pfp_true_end_process;
  std::vector<double> slc_pfp_true_p_x, slc_pfp_true_p_y, slc_pfp_true_p_z;
  std::vector<double> slc_pfp_completeness;
  std::vector<double> slc_pfp_purity;

  std::vector<double> slc_pfp_trk_start_x, slc_pfp_trk_start_y, slc_pfp_trk_start_z;
  std::vector<int> slc_pfp_trk_start_wires_U, slc_pfp_trk_start_wires_V, slc_pfp_trk_start_wires_C;
  std::vector<double> slc_pfp_trk_start_time_tick;
  std::vector<int> slc_pfp_trk_start_channel_U, slc_pfp_trk_start_channel_V, slc_pfp_trk_start_channel_C;

  std::vector<double> slc_pfp_trk_end_x, slc_pfp_trk_end_y, slc_pfp_trk_end_z;
  std::vector<int> slc_pfp_trk_end_wires_U, slc_pfp_trk_end_wires_V, slc_pfp_trk_end_wires_C;
  std::vector<double> slc_pfp_trk_end_time_tick;
  std::vector<int> slc_pfp_trk_end_channel_U, slc_pfp_trk_end_channel_V, slc_pfp_trk_end_channel_C;

  std::vector<double> slc_pfp_track_dir_x, slc_pfp_track_dir_y, slc_pfp_track_dir_z;
  std::vector<std::vector<double>> slc_pfp_track_dir_vec_x, slc_pfp_track_dir_vec_y, slc_pfp_track_dir_vec_z;
  std::vector<double>	slc_pfp_track_length;

  std::vector<double>  slc_pfp_track_kinetic_energy, slc_pfp_track_visible_energy, slc_pfp_track_charge;

  std::vector<float>	slc_pfp_track_chi2_muon, slc_pfp_track_chi2_pion, slc_pfp_track_chi2_kaon, slc_pfp_track_chi2_proton;
  std::vector<float> 	slc_pfp_track_theta, slc_pfp_track_phi;



  std::vector<double> slc_pfp_trk_mcs_mom_muon;
  std::vector<double> slc_pfp_trk_mcs_mom_pion;
  std::vector<double> slc_pfp_trk_range_mom_muon;
  std::vector<double> slc_pfp_trk_range_mom_pion;

  std::vector<double> slc_pfp_trk_meanDCA;  
  std::vector<double> slc_pfp_trk_range_mom;
  std::vector<double> slc_pfp_trk_mcs_mom;
  std::vector<double> slc_pfp_trk_mom_diff;
  std::vector<double> slc_pfp_trk_stopping_dEdx_chi2_ratio;
  std::vector<double> slc_pfp_trk_chi2_pol0_dEdx_fit;
  std::vector<double> slc_pfp_trk_mcs_scatter_mean;
  std::vector<double> slc_pfp_trk_mcs_scatter_max_ratio;

  std::vector<double> slc_pfp_shower_dir_x, slc_pfp_shower_dir_y, slc_pfp_shower_dir_z;
  std::vector<double>	slc_pfp_shower_length;
  std::vector<double>	slc_pfp_shower_opening_angle;
  std::vector<double>	slc_pfp_shower_energy;
  std::vector<double>	slc_pfp_shower_dEdx;


  //Declaration of fhicl parameters
  const std::string file_type;

  const std::string MCTruth_label;
  const std::string MCParticle_label;
  const std::string POTModule_label;

  const std::string slice_label;
  const std::string pfp_label;
  const std::string track_label;
  const std::string shower_label;
  const std::string hit_label;
  const std::string space_point_label;

  const std::string pfp_to_razzled_label;

  const std::string slice_to_hit_label;
  const std::string slice_to_pfp_label;
  const std::string slice_to_crumbs_label;
  const std::string slice_to_opt0_label;
  const std::string slice_to_stub_label;
  const std::string slice_to_barycenter_matching_label;
  const std::string flash_to_barycenter_matching_label;

  const std::string pfp_to_pfpmetadata_label;
  const std::string pfp_to_vertex_label;
  const std::string pfp_to_track_label;
  const std::string pfp_to_shower_label;
  const std::string pfp_to_space_point_label;
  const std::string pfp_to_cnn_id_label;

  const std::string track_to_hits_label;
  const std::string shower_to_hits_label;
  const std::string track_to_calo_label;
  const std::string track_to_chi2_label;
  const std::string shower_to_calo_label;

  art::InputTag fMCSLabel_muon;
  art::InputTag fMCSLabel_pion;
  art::InputTag fRangeLabel_muon;
  art::InputTag fRangeLabel_pion;
  art::InputTag fClosestApproachLabel;
  art::InputTag fStoppingChi2Label;


  const std::string space_point_to_hit_label;

  bool verbose_general = false;
  bool verbose_g4 = false;
  bool verbose_gen = false;
  bool verbose_reco = false;
  bool verbose_reco_hits = false;
  bool verbose_reco_pfps = false;
  bool verbose_reco_calo_deposition = false;
  bool verbose_reco_space_points = false;

  bool save_truth_gen_particles = false;
  bool save_truth_g4 = false;
  bool save_reco_hits = false;
  bool save_reco = false;
  bool save_track_calo_deposition = false;
  bool save_reco_space_points = true;
  bool save_reco_stubs = true;
  bool save_crumbs = true;
  bool save_razzled = true;
  bool save_3d_barycenter_matching = true;

  bool do_precuts = true;

  std::vector<double> PosNull = {0., 0., 0.};
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;
  double fOpticalClock;
};

//Initialize some variables usingn the fhicl file "analysisConfigTruth"
test::AnalyzeEventsNeutrino::AnalyzeEventsNeutrino(fhicl::ParameterSet const &p)
  : EDAnalyzer{p},
  // More initializers here.
    file_type(p.get<std::string>("file_type")),
    MCTruth_label(p.get<std::string>("MCTruth_label")),
    MCParticle_label(p.get<std::string>("MCParticle_label")),
    POTModule_label(p.get<std::string>("POTModule_label")),

    slice_label(p.get<std::string>("slice_label")),
    pfp_label(p.get<std::string>("pfp_label")),
    track_label(p.get<std::string>("track_label")),
    shower_label(p.get<std::string>("shower_label")),
    hit_label(p.get<std::string>("hit_label")),
    space_point_label(p.get<std::string>("space_point_label")),

    pfp_to_razzled_label(p.get<std::string>("pfp_to_razzled_label")),

    slice_to_hit_label(p.get<std::string>("slice_to_hit_label")),
    slice_to_pfp_label(p.get<std::string>("slice_to_pfp_label")),
    slice_to_crumbs_label(p.get<std::string>("slice_to_crumbs_label")),
    slice_to_opt0_label(p.get<std::string>("slice_to_opt0_label")),
    slice_to_stub_label(p.get<std::string>("slice_to_stub_label")),
    slice_to_barycenter_matching_label(p.get<std::string>("slice_to_barycenter_matching_label")),
    flash_to_barycenter_matching_label(p.get<std::string>("flash_to_barycenter_matching_label")),

    pfp_to_pfpmetadata_label(p.get<std::string>("pfp_to_pfpmetadata_label")),
    pfp_to_vertex_label(p.get<std::string>("pfp_to_vertex_label")),
    pfp_to_track_label(p.get<std::string>("pfp_to_track_label")),
    pfp_to_shower_label(p.get<std::string>("pfp_to_shower_label")),
    pfp_to_space_point_label(p.get<std::string>("pfp_to_space_point_label")),
    pfp_to_cnn_id_label(p.get<std::string>("pfp_to_cnn_id_label")),

    track_to_hits_label(p.get<std::string>("track_to_hits_label")),
    shower_to_hits_label(p.get<std::string>("shower_to_hits_label")),
    track_to_calo_label(p.get<std::string>("track_to_calo_label")),
    track_to_chi2_label(p.get<std::string>("track_to_chi2_label")),
    shower_to_calo_label(p.get<std::string>("shower_to_calo_label")),

    fMCSLabel_muon(p.get<std::string>("fMCSLabel"), std::string("muon")),
    fMCSLabel_pion(p.get<std::string>("fMCSLabel"), std::string("pion")),
    fRangeLabel_muon(p.get<std::string>("fRangeLabel"), std::string("muon")),
    fRangeLabel_pion(p.get<std::string>("fRangeLabel"), std::string("pion")),
    fClosestApproachLabel(p.get<std::string>("fClosestApproachLabel")),
    fStoppingChi2Label(p.get<std::string>("fStoppingChi2Label")),

    space_point_to_hit_label(p.get<std::string>("space_point_to_hit_label")),

    verbose_general(p.get<bool>("Verbose_general")),
    verbose_g4(p.get<bool>("Verbose_g4")),
    verbose_gen(p.get<bool>("Verbose_gen")),
    verbose_reco(p.get<bool>("Verbose_reco")),
    verbose_reco_hits(p.get<bool>("Verbose_reco_hits")),
    verbose_reco_pfps(p.get<bool>("Verbose_reco_pfps")),
    verbose_reco_calo_deposition(p.get<bool>("Verbose_reco_calo_deposition")),
    verbose_reco_space_points(p.get<bool>("Verbose_reco_space_points")),

    save_truth_gen_particles(p.get<bool>("save_truth_gen_particles")),
    save_truth_g4(p.get<bool>("save_truth_g4")),
    save_reco_hits(p.get<bool>("save_reco_hits")),
    save_reco(p.get<bool>("save_reco")),
    save_track_calo_deposition(p.get<bool>("save_track_calo_deposition")),
    save_reco_space_points(p.get<bool>("save_reco_space_points")),
    save_reco_stubs(p.get<bool>("save_reco_stubs")),
    save_crumbs(p.get<bool>("save_crumbs")),
    save_razzled(p.get<bool>("save_razzled")),
    save_3d_barycenter_matching(p.get<bool>("save_3d_barycenter_matching")),

    do_precuts(p.get<bool>("do_precuts"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  detinfo::DetectorClocksData const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fOpticalClock =  clockData.TriggerTime();
  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  geo::PlaneID PID(0,1,0);
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  fWirePlanePosition = std::abs( wireReadout.Plane(PID).GetCenter().X() );
}

//Gets the number of events generated in each subrun
int test::AnalyzeEventsNeutrino::get_total_gen_evts(const art::Event &e) {
  int nGenEvt = 0;
  for (const art::ProcessConfiguration &process: e.processHistory()) {
    std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
    if (genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") &&
    genConfig->has_key("source.module_type")) {
      int maxEvents = genConfig->get<int>("source.maxEvents");
      std::string moduleType = genConfig->get<std::string>("source.module_type");
      if (moduleType == "EmptyEvent") {
        nGenEvt += maxEvents;
      }
    }
  }
  return nGenEvt;
}

//Code ran at the start of each subrun
void test::AnalyzeEventsNeutrino::beginSubRun(const art::SubRun &sr) {
  //Load the POT
  art::Handle <sumdata::POTSummary> POT_handle;
  sr.getByLabel(POTModule_label, POT_handle);
  if (!POT_handle.isValid()) {
    std::cout << "POT product " << POTModule_label << " not found..." << std::endl;
    //throw std::exception();
  }

  //Save the POT information of the subrun on the specified variables, values are given for events with no valid POT
  POT = 0;
  spill = 0;
  if (!POT_handle.isValid()) {
    POT = 0;
    spill = 1;
  } else {
    POT = POT_handle->totpot;
    spill = POT_handle->totspills;
  }
}

//Code ran at the end of each subrun
void test::AnalyzeEventsNeutrino::endSubRun(const art::SubRun &sr) {
  fSubRunTree->Fill();
}

void test::AnalyzeEventsNeutrino::reset_g4_vars() {
  //Reset the vectors with geant4 particles information
  g4_part_trackID.clear();
  g4_part_mother.clear();
  g4_part_PDGcode.clear();
  g4_part_mass.clear();
  g4_part_E0.clear();
  g4_part_Ef.clear();
  g4_part_TL.clear();

  g4_part_start_pos_X.clear();
  g4_part_start_pos_Y.clear();
  g4_part_start_pos_Z.clear();
  g4_part_start_T.clear();
  g4_part_end_pos_X.clear();
  g4_part_end_pos_Y.clear();
  g4_part_end_pos_Z.clear();
  g4_part_end_T.clear();

  g4_part_P0_X.clear();
  g4_part_P0_Y.clear();
  g4_part_P0_Z.clear();
  g4_part_Pf_X.clear();
  g4_part_Pf_Y.clear();
  g4_part_Pf_Z.clear();

  g4_part_process.clear();
  g4_part_end_process.clear();
}

void test::AnalyzeEventsNeutrino::reset_gen_vars() {
  //Reset all variables at the start of the analysis of an event
  nu_PDG = 0;
  nu_E0 = 0;
  nu_weight = 0;
  nu_interaction_mode = 0;
  nu_interaction_type = 0;
  nu_CC_NC = 0;
  nu_target = 0;
  nu_HitNuc = 0;
  nu_HitQuark = 0;
  nu_W = 0;
  nu_X = 0;
  nu_Y = 0;
  nu_QSqr = 0;
  gen_index = 0;


  slc_true_prim_vtx_x = 0;
  slc_true_prim_vtx_y = 0;
  slc_true_prim_vtx_z = 0;
  slc_true_prim_vtx_wires_U = 0;
  slc_true_prim_vtx_wires_V = 0;
  slc_true_prim_vtx_wires_C = 0;
  slc_true_prim_vtx_time_tick = 0;
  slc_true_prim_vtx_channel_U = 0;
  slc_true_prim_vtx_channel_V = 0;
  slc_true_prim_vtx_channel_C = 0;


  //Reset the vectors with MCTruth particles information
  gen_part_trackID.clear();
  gen_part_statusCode.clear();
  gen_part_mother.clear();
  gen_part_PDGcode.clear();
  gen_part_mass.clear();
  gen_part_E0.clear();
  gen_part_start_pos_X.clear();
  gen_part_start_pos_Y.clear();
  gen_part_start_pos_Z.clear();
  gen_part_P0_X.clear();
  gen_part_P0_Y.clear();
  gen_part_P0_Z.clear();
}

void test::AnalyzeEventsNeutrino::fill_empty_slice_pfp_track_info() {
  slc_pfp_true_trackid.push_back(-1);
  slc_pfp_true_pdg.push_back(-1);
  slc_pfp_true_TL.push_back(-1);
  slc_pfp_true_energy.push_back(-1);
  slc_pfp_true_end_process.push_back("null");
  slc_pfp_true_p_x.push_back(-1);
  slc_pfp_true_p_y.push_back(-1);
  slc_pfp_true_p_z.push_back(-1);
  slc_pfp_completeness.push_back(-1);
  slc_pfp_purity.push_back(-1);

  slc_pfp_has_track.push_back(false);

  slc_pfp_trk_start_x.push_back(-1);
  slc_pfp_trk_start_y.push_back(-1);
  slc_pfp_trk_start_z.push_back(-1);
  slc_pfp_trk_start_wires_U.push_back(-1);
  slc_pfp_trk_start_wires_V.push_back(-1);
  slc_pfp_trk_start_wires_C.push_back(-1);
  slc_pfp_trk_start_time_tick.push_back(-1);
  slc_pfp_trk_start_channel_U.push_back(-1);
  slc_pfp_trk_start_channel_V.push_back(-1);
  slc_pfp_trk_start_channel_C.push_back(-1);

  slc_pfp_trk_end_x.push_back(-1);
  slc_pfp_trk_end_y.push_back(-1);
  slc_pfp_trk_end_z.push_back(-1);
  slc_pfp_trk_end_wires_U.push_back(-1);
  slc_pfp_trk_end_wires_V.push_back(-1);
  slc_pfp_trk_end_wires_C.push_back(-1);
  slc_pfp_trk_end_time_tick.push_back(-1);
  slc_pfp_trk_end_channel_U.push_back(-1);
  slc_pfp_trk_end_channel_V.push_back(-1);
  slc_pfp_trk_end_channel_C.push_back(-1);

  slc_pfp_track_dir_x.push_back(-1);
  slc_pfp_track_dir_y.push_back(-1);
  slc_pfp_track_dir_z.push_back(-1);
  slc_pfp_track_dir_vec_x.push_back({-1});
  slc_pfp_track_dir_vec_y.push_back({-1});
  slc_pfp_track_dir_vec_z.push_back({-1});

  slc_pfp_track_length.push_back(-1);
  slc_pfp_track_kinetic_energy.push_back(-1);
  slc_pfp_track_visible_energy.push_back(-1);
  slc_pfp_track_charge.push_back(-1);

  slc_pfp_track_chi2_muon.push_back(-1);
  slc_pfp_track_chi2_pion.push_back(-1);
  slc_pfp_track_chi2_kaon.push_back(-1);
  slc_pfp_track_chi2_proton.push_back(-1);
  slc_pfp_track_theta.push_back(-1);
  slc_pfp_track_phi.push_back(-1);

  slc_pfp_trk_meanDCA.push_back(-1);
  slc_pfp_trk_range_mom.push_back(-1);
  slc_pfp_trk_mcs_mom.push_back(-1);
  slc_pfp_trk_mcs_mom_muon.push_back(-1);
  slc_pfp_trk_mcs_mom_pion.push_back(-1);
  slc_pfp_trk_range_mom_muon.push_back(-1);
  slc_pfp_trk_range_mom_pion.push_back(-1);
  slc_pfp_trk_mom_diff.push_back(-1);
  slc_pfp_trk_stopping_dEdx_chi2_ratio.push_back(-1);
  slc_pfp_trk_chi2_pol0_dEdx_fit.push_back(-1);
  slc_pfp_trk_mcs_scatter_mean.push_back(-1);
  slc_pfp_trk_mcs_scatter_max_ratio.push_back(-1);

  slc_pfp_calo_dEdx.push_back({-1});
  slc_pfp_calo_dQdx.push_back({-1});
  slc_pfp_calo_pitch.push_back({-1});
  slc_pfp_calo_residual_range.push_back({-1});
  slc_pfp_calo_x.push_back({-1});
  slc_pfp_calo_y.push_back({-1});
  slc_pfp_calo_z.push_back({-1});
}

void test::AnalyzeEventsNeutrino::FillPandoraTrackIDScoreVar(const std::map<std::string, float> &propertiesMap, std::vector<double> &var, std::string name)
{
  auto propertiesMapIter = propertiesMap.find(name);
  if(propertiesMapIter == propertiesMap.end()) {
    var.push_back(-1);
  } else var.push_back(propertiesMapIter->second);
}


void test::AnalyzeEventsNeutrino::FillPandoraTrackIDScoreVars(const std::map<std::string, float> &propertiesMap)
{
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_track_score, "TrackScore");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_charge_end_frac, "LArThreeDChargeFeatureTool_EndFraction");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_charge_frac_spread, "LArThreeDChargeFeatureTool_FractionalSpread");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_linear_fit_diff, "LArThreeDLinearFitFeatureTool_DiffStraightLineMean");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_linear_fit_length, "LArThreeDLinearFitFeatureTool_Length");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_linear_fit_gap_length, "LArThreeDLinearFitFeatureTool_MaxFitGapLength");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_linear_fit_RMS, "LArThreeDLinearFitFeatureTool_SlidingLinearFitRMS");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_open_angle_diff, "LArThreeDOpeningAngleFeatureTool_AngleDiff");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_secondary_PCA_ratio, "LArThreeDPCAFeatureTool_SecondaryPCARatio");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_tertiary_PCA_ratio, "LArThreeDPCAFeatureTool_TertiaryPCARatio");
  this->FillPandoraTrackIDScoreVar(propertiesMap, slc_pfp_vertex_dist, "LArThreeDVertexDistanceFeatureTool_VertexDistance");
}


void test::AnalyzeEventsNeutrino::fill_empty_slice_pfp_shower_info() {
  slc_pfp_shower_dir_x.push_back(-1);
  slc_pfp_shower_dir_y.push_back(-1);
  slc_pfp_shower_dir_z.push_back(-1);
  slc_pfp_shower_length.push_back(-1);
  slc_pfp_shower_opening_angle.push_back(-1);
  slc_pfp_shower_energy.push_back(-1);
  slc_pfp_shower_dEdx.push_back(-1);
}

void test::AnalyzeEventsNeutrino::fill_empty_slice() {
  slc_is_reconstructed = false;

  slc_prim_vtx_x = -1;
  slc_prim_vtx_y = -1;
  slc_prim_vtx_z = -1;
  slc_prim_vtx_wires_U = -1;
  slc_prim_vtx_wires_V = -1;
  slc_prim_vtx_wires_C = -1;
  slc_prim_vtx_time_tick = -1;
  slc_prim_vtx_channel_U = -1;
  slc_prim_vtx_channel_V = -1;
  slc_prim_vtx_channel_C = -1;

  //Fill with -1 all slice parameters
  slc_nu_score = -1;
  slc_crumbs_score = -1;

  slc_t0 = -1;
  slc_t0_score = -1;
  slc_OpT0Fraction = -1;

  slc_bcmatch_flash_time = -10000;
  slc_bcmatch_first_hit = -10000;
  slc_bcmatch_flash_PE = -10000;
  slc_bcmatch_flash_asymmetry = -10000;

  slc_sp_x.push_back(-1);
  slc_sp_y.push_back(-1);
  slc_sp_z.push_back(-1);
  slc_sp_chi2.push_back(-1);
  slc_sp_integral.push_back(-1);
  slc_sp_sigma_integral.push_back(-1);
  slc_sp_associated_pfp_ID.push_back(-1);

  slc_stub_vtx_x.push_back(-1);
  slc_stub_vtx_y.push_back(-1);
  slc_stub_vtx_z.push_back(-1);
  slc_stub_end_x.push_back(-1);
  slc_stub_end_y.push_back(-1);
  slc_stub_end_z.push_back(-1);
  slc_stub_hit_charge.push_back({-1});
  slc_stub_hit_wire.push_back({-1});
  slc_stub_hit_ontrack.push_back({0});

  slc_hits_wire_ID.push_back(-1);
  slc_hits_plane_ID.push_back(-1);
  slc_hits_TPC_ID.push_back(-1);
  slc_hits_channel_ID.push_back(-1);
  slc_hits_integral.push_back(-1);
  slc_hits_sigma_integral.push_back(-1);
  slc_hits_peaktime.push_back(-1);
  slc_hits_RMS.push_back(-1);
  slc_hits_associated_pfp_ID.push_back(-1);

  slc_pfp_ID.push_back(-1);
  slc_pfp_daughter_max_hits.push_back(-1);
  slc_pfp_pandora_PDG.push_back(-1);
  slc_pfp_parent_ID.push_back(-1);
  slc_pfp_mother_is_primary.push_back(-1);
  slc_pfp_track_score.push_back(-1);

  slc_pfp_charge_end_frac.push_back(-1);
  slc_pfp_charge_frac_spread.push_back(-1);
  slc_pfp_linear_fit_diff.push_back(-1);
  slc_pfp_linear_fit_length.push_back(-1);
  slc_pfp_linear_fit_gap_length.push_back(-1);
  slc_pfp_linear_fit_RMS.push_back(-1);
  slc_pfp_open_angle_diff.push_back(-1);
  slc_pfp_secondary_PCA_ratio.push_back(-1);
  slc_pfp_tertiary_PCA_ratio.push_back(-1);
  slc_pfp_vertex_dist.push_back(-1);

  slc_pfp_razzled_score_e.push_back(-1);
  slc_pfp_razzled_score_mu.push_back(-1);
  slc_pfp_razzled_score_pho.push_back(-1);
  slc_pfp_razzled_score_pi.push_back(-1);
  slc_pfp_razzled_score_p.push_back(-1);
  slc_pfp_razzled_pdg.push_back(-1);

  slc_pfp_vtx_x.push_back(-1);
  slc_pfp_vtx_y.push_back(-1);
  slc_pfp_vtx_z.push_back(-1);
  slc_pfp_vtx_wires_U.push_back(-1);
  slc_pfp_vtx_wires_V.push_back(-1);
  slc_pfp_vtx_wires_C.push_back(-1);
  slc_pfp_vtx_time_tick.push_back(-1);
  slc_pfp_vtx_channel_U.push_back(-1);
  slc_pfp_vtx_channel_V.push_back(-1);
  slc_pfp_vtx_channel_C.push_back(-1);

  fill_empty_slice_pfp_track_info();
  fill_empty_slice_pfp_shower_info();
}

void test::AnalyzeEventsNeutrino::reset_slc_vars() {
  //Reset all variables concerning genreator stage
  reset_gen_vars();

  slc_ID = -1;
  slc_is_reconstructed = true;
  slc_is_obvious_cosmic = false;

  slc_prim_vtx_x = 0;
  slc_prim_vtx_y = 0;
  slc_prim_vtx_z = 0;
  slc_prim_vtx_wires_U = 0;
  slc_prim_vtx_wires_V = 0;
  slc_prim_vtx_wires_C = 0;
  slc_prim_vtx_time_tick = 0;
  slc_prim_vtx_channel_U = 0;
  slc_prim_vtx_channel_V = 0;
  slc_prim_vtx_channel_C = 0;

  slc_truth_matching_completeness = 0;
  slc_truth_matching_purity = 0;
  slc_nu_score = 0;
  slc_crumbs_score = 0;

  slc_t0 = 0;
  slc_t0_score = 0;
  slc_OpT0Fraction = 0;

  slc_bcmatch_flash_time = 0;
  slc_bcmatch_first_hit = 0;
  slc_bcmatch_flash_PE = 0;
  slc_bcmatch_flash_asymmetry = 0;

  slc_sp_x.clear();
  slc_sp_y.clear();
  slc_sp_z.clear();
  slc_sp_chi2.clear();
  slc_sp_integral.clear();
  slc_sp_sigma_integral.clear();
  slc_sp_associated_pfp_ID.clear();

  slc_stub_vtx_x.clear();
  slc_stub_vtx_y.clear();
  slc_stub_vtx_z.clear();
  slc_stub_end_x.clear();
  slc_stub_end_y.clear();
  slc_stub_end_z.clear();
  slc_stub_hit_charge.clear();
  slc_stub_hit_wire.clear();
  slc_stub_hit_ontrack.clear();

  slc_hits_wire_ID.clear();
  slc_hits_plane_ID.clear();
  slc_hits_TPC_ID.clear();
  slc_hits_channel_ID.clear();
  slc_hits_integral.clear();
  slc_hits_sigma_integral.clear();
  slc_hits_peaktime.clear();
  slc_hits_RMS.clear();
  slc_hits_associated_pfp_ID.clear();

  slc_pfp_ID.clear();
  slc_pfp_daughter_max_hits.clear();
  slc_pfp_has_track.clear();
  slc_pfp_parent_ID.clear();
  slc_pfp_mother_is_primary.clear();
  slc_pfp_track_score.clear();
  slc_pfp_pandora_PDG.clear();

  slc_pfp_charge_end_frac.clear();
  slc_pfp_charge_frac_spread.clear();
  slc_pfp_linear_fit_diff.clear();
  slc_pfp_linear_fit_length.clear();
  slc_pfp_linear_fit_gap_length.clear();
  slc_pfp_linear_fit_RMS.clear();
  slc_pfp_open_angle_diff.clear();
  slc_pfp_secondary_PCA_ratio.clear();
  slc_pfp_tertiary_PCA_ratio.clear();
  slc_pfp_vertex_dist.clear();

  slc_pfp_razzled_score_e.clear();
  slc_pfp_razzled_score_mu.clear();
  slc_pfp_razzled_score_pho.clear();
  slc_pfp_razzled_score_pi.clear();
  slc_pfp_razzled_score_p.clear();
  slc_pfp_razzled_pdg.clear();

  slc_pfp_vtx_x.clear();
  slc_pfp_vtx_y.clear();
  slc_pfp_vtx_z.clear();

  slc_pfp_vtx_wires_U.clear();
  slc_pfp_vtx_wires_V.clear();
  slc_pfp_vtx_wires_C.clear();
  slc_pfp_vtx_time_tick.clear();
  slc_pfp_vtx_channel_U.clear();
  slc_pfp_vtx_channel_V.clear();
  slc_pfp_vtx_channel_C.clear();

  slc_pfp_true_trackid.clear();
  slc_pfp_true_pdg.clear();
  slc_pfp_true_TL.clear();
  slc_pfp_true_energy.clear();
  slc_pfp_true_end_process.clear();
  slc_pfp_true_p_x.clear();
  slc_pfp_true_p_y.clear();
  slc_pfp_true_p_z.clear();
  slc_pfp_completeness.clear();
  slc_pfp_purity.clear();

  slc_pfp_calo_dEdx.clear();
  slc_pfp_calo_dQdx.clear();
  slc_pfp_calo_pitch.clear();
  slc_pfp_calo_residual_range.clear();
  slc_pfp_calo_x.clear();
  slc_pfp_calo_y.clear();
  slc_pfp_calo_z.clear();

  slc_pfp_trk_start_x.clear();
  slc_pfp_trk_start_y.clear();
  slc_pfp_trk_start_z.clear();

  slc_pfp_trk_start_wires_U.clear();
  slc_pfp_trk_start_wires_V.clear();
  slc_pfp_trk_start_wires_C.clear();
  slc_pfp_trk_start_time_tick.clear();
  slc_pfp_trk_start_channel_U.clear();
  slc_pfp_trk_start_channel_V.clear();
  slc_pfp_trk_start_channel_C.clear();

  slc_pfp_trk_end_x.clear();
  slc_pfp_trk_end_y.clear();
  slc_pfp_trk_end_z.clear();

  slc_pfp_trk_end_wires_U.clear();
  slc_pfp_trk_end_wires_V.clear();
  slc_pfp_trk_end_wires_C.clear();
  slc_pfp_trk_end_time_tick.clear();
  slc_pfp_trk_end_channel_U.clear();
  slc_pfp_trk_end_channel_V.clear();
  slc_pfp_trk_end_channel_C.clear();

  slc_pfp_track_dir_x.clear();
  slc_pfp_track_dir_y.clear();
  slc_pfp_track_dir_z.clear();
  slc_pfp_track_dir_vec_x.clear();
  slc_pfp_track_dir_vec_y.clear();
  slc_pfp_track_dir_vec_z.clear();

  slc_pfp_track_length.clear();
  slc_pfp_track_kinetic_energy.clear();
  slc_pfp_track_visible_energy.clear();
  slc_pfp_track_charge.clear();

  slc_pfp_track_chi2_muon.clear();
  slc_pfp_track_chi2_pion.clear();
  slc_pfp_track_chi2_kaon.clear();
  slc_pfp_track_chi2_proton.clear();
  slc_pfp_track_theta.clear();
  slc_pfp_track_phi.clear();

  slc_pfp_trk_meanDCA.clear();
  slc_pfp_trk_range_mom.clear();
  slc_pfp_trk_mcs_mom.clear();
  slc_pfp_trk_mcs_mom_muon.clear();
  slc_pfp_trk_mcs_mom_pion.clear();
  slc_pfp_trk_range_mom_muon.clear();
  slc_pfp_trk_range_mom_pion.clear();
  slc_pfp_trk_mom_diff.clear();
  slc_pfp_trk_stopping_dEdx_chi2_ratio.clear();
  slc_pfp_trk_chi2_pol0_dEdx_fit.clear();
  slc_pfp_trk_mcs_scatter_mean.clear();
  slc_pfp_trk_mcs_scatter_max_ratio.clear();

  slc_pfp_shower_dir_x.clear();
  slc_pfp_shower_dir_y.clear();
  slc_pfp_shower_dir_z.clear();
  slc_pfp_shower_length.clear();
  slc_pfp_shower_opening_angle.clear();
  slc_pfp_shower_energy.clear();
  slc_pfp_shower_dEdx.clear();
}


bool test::AnalyzeEventsNeutrino::slc_check_clear_cosmic(std::vector<art::Ptr<recob::PFParticle>> pfp_vec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns) {
  bool is_cosmic = false;
  //std::cout << "   Cosmic start " << pfp_vec.size() <<  std::endl;
  for(const art::Ptr <recob::PFParticle> &pfp: pfp_vec) {
    bool pfp_is_cosmic = false;
    const std::vector<art::Ptr < larpandoraobj::PFParticleMetadata>> pfp_meta_vec = pfp_pfpmeta_assns.at(pfp.key());
    art::Ptr <larpandoraobj::PFParticleMetadata> pfp_meta = pfp_meta_vec.front();
    std::map<std::string, float> pfp_prop_map = pfp_meta->GetPropertiesMap();

    if(pfp_prop_map.find("IsClearCosmic") != pfp_prop_map.end()) {
      pfp_is_cosmic = pfp_prop_map["IsClearCosmic"];
    }

    if (pfp_is_cosmic) is_cosmic = true;
  }
  //std::cout << "   Cosmic end" << std::endl;
  return is_cosmic;
}

double test::AnalyzeEventsNeutrino::get_slice_nu_score(std::vector<art::Ptr<recob::PFParticle>> pfp_vec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns) {
  double nu_score = -1.;

  for(const art::Ptr<recob::PFParticle> &pfp : pfp_vec){
    if(pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 )) {
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfp_meta_vec = pfp_pfpmeta_assns.at(pfp.key());

      art::Ptr<larpandoraobj::PFParticleMetadata> nu_meta = pfp_meta_vec[0];
      std::map<std::string, float> nu_prop_map = nu_meta->GetPropertiesMap();

      if ( nu_prop_map.find("NuScore") != nu_prop_map.end() ) {
        nu_score = nu_prop_map["NuScore"];
      }
    }
  }
  std::cout << "NuScore: " << nu_score << std::endl;

  return nu_score;
}

int test::AnalyzeEventsNeutrino::vertex_to_drift_tick(double vt, double vx) {
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}

std::vector<int> test::AnalyzeEventsNeutrino::get_point_in_wires(double x, double y, double z) {
  //Get Vertex in wires
  int tpc = 0;
  if(x > 0)  tpc = 1;

  geo::PlaneID UPlane_ID(0,tpc,0);
  geo::PlaneID VPlane_ID(0,tpc,1);
  geo::PlaneID CPlane_ID(0,tpc,2);


  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  geo::PlaneGeo const& plane =  wireReadout.Plane(UPlane_ID);
  geo::WireID wire_ID;
  try {
    wire_ID = plane.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wire_ID = plane.ClosestWireID(e.suggestedWireID());
  }
  int reco_VU = wire_ID.deepestIndex();
  int channel_VU = wireReadout.PlaneWireToChannel(wire_ID);

  geo::PlaneGeo const& plane1 =  wireReadout.Plane(VPlane_ID);
  try {
    wire_ID = plane1.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wire_ID = plane1.ClosestWireID(e.suggestedWireID());
  }
  int reco_VV = wire_ID.deepestIndex();
  int channel_VV = wireReadout.PlaneWireToChannel(wire_ID);

  geo::PlaneGeo const& plane2 =  wireReadout.Plane(CPlane_ID);
  try {
    wire_ID = plane2.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wire_ID = plane2.ClosestWireID(e.suggestedWireID());
  }
  int reco_VC = wire_ID.deepestIndex();
  int channel_VC = wireReadout.PlaneWireToChannel(wire_ID);

  int reco_time_tick = test::AnalyzeEventsNeutrino::vertex_to_drift_tick(0, x);

  return {reco_VU, reco_VV, reco_VC, reco_time_tick, channel_VU, channel_VV, channel_VC};
}

bool test::AnalyzeEventsNeutrino::get_slice_nu_vertex(art::Event const &e,std::vector<art::Ptr<recob::PFParticle>> slice_pfp_vec, art::FindManyP<recob::Vertex> pfp_vertex_assns) {
  bool vertex_found = false;
  //Loop to find the RecoVertex
/*
  art::FindManyP<recob::Vertex> fmVertex = art::FindManyP<recob::Vertex>(slice_pfp_vec, e,pfp_to_vertex_label);
  // get the primary particle
  size_t iPart;
  for (iPart = 0; iPart < slice_pfp_vec.size(); ++iPart ) {
    const recob::PFParticle &thisParticle = *slice_pfp_vec[iPart];
    if (thisParticle.IsPrimary()) break;
  }
  // get the primary vertex
  const recob::Vertex *vertex = (iPart == slice_pfp_vec.size() || !fmVertex.at(iPart).size()) ? NULL : fmVertex.at(iPart).at(0).get();
  if (vertex != NULL) {
    slc_prim_vtx_x = vertex->position().X();
    slc_prim_vtx_y = vertex->position().Y();
    slc_prim_vtx_z = vertex->position().Z();
    return true;
  } else {
    return false;
  }
*/



  for(const art::Ptr<recob::PFParticle> &pfp : slice_pfp_vec) {
    //Only continue if it is a neutrino and a primary particle
    //if (!pfp->IsPrimary()) continue;
    std::vector<art::Ptr < recob::Vertex>> pfp_vertex_vec = pfp_vertex_assns.at(pfp.key());

    for (const art::Ptr <recob::Vertex> &vertex: pfp_vertex_vec) {
      std::cout << "vtx: " << vertex->position().X() << " " << vertex->position().Y() << " "<<  vertex->position().Z() << std::endl;
      if(vertex_found) continue;
      vertex_found = true;
      //Get the vertex information
      slc_prim_vtx_x = vertex->position().X();
      slc_prim_vtx_y = vertex->position().Y();
      slc_prim_vtx_z = vertex->position().Z();
      double xyz_vertex[3];
      vertex->XYZ(xyz_vertex);

      if (!fGeom->HasTPC(fGeom->FindTPCAtPosition(vertex->position())))  {
        slc_prim_vtx_wires_U = -1;
        slc_prim_vtx_wires_V = -1;
        slc_prim_vtx_wires_C = -1;
        slc_prim_vtx_time_tick = -1;
        slc_prim_vtx_channel_U = -1;
        slc_prim_vtx_channel_V = -1;
        slc_prim_vtx_channel_C = -1;
      } else {
        std::vector<int> reco_vertex_nu_wire = get_point_in_wires(xyz_vertex[0], xyz_vertex[1], xyz_vertex[2]);
        slc_prim_vtx_wires_U = reco_vertex_nu_wire[0];
        slc_prim_vtx_wires_V = reco_vertex_nu_wire[1];
        slc_prim_vtx_wires_C = reco_vertex_nu_wire[2];
        slc_prim_vtx_time_tick = reco_vertex_nu_wire[3];
        slc_prim_vtx_channel_U = reco_vertex_nu_wire[4];
        slc_prim_vtx_channel_V = reco_vertex_nu_wire[5];
        slc_prim_vtx_channel_C = reco_vertex_nu_wire[6];
      }
    }
    if(vertex_found) break;

  }

  return vertex_found;
}


void test::AnalyzeEventsNeutrino::fill_empty_truth() {
  slc_truth_matching_completeness = -1;
  slc_truth_matching_purity = -1;

  nu_PDG = -1;
  nu_E0 = -1;
  nu_weight = -1;

  nu_interaction_mode = -1;
  nu_interaction_type = -1;
  nu_CC_NC = -1;

  nu_HitNuc = -1;
  nu_target = -1;
  nu_HitQuark = -1;
  nu_W = -1;
  nu_X = -1;
  nu_Y = -1;
  nu_QSqr = -1;

  slc_true_prim_vtx_x = -1;
  slc_true_prim_vtx_y = -1;
  slc_true_prim_vtx_z = -1;
  slc_true_prim_vtx_wires_U = -1;
  slc_true_prim_vtx_wires_V = -1;
  slc_true_prim_vtx_wires_C = -1;
  slc_true_prim_vtx_time_tick = -1;
  slc_true_prim_vtx_channel_U = -1;
  slc_true_prim_vtx_channel_V = -1;
  slc_true_prim_vtx_channel_C = -1;

  gen_part_trackID.push_back(-1);
  gen_part_statusCode.push_back(-1);
  gen_part_mother.push_back(-1);
  gen_part_PDGcode.push_back(-1);
  gen_part_mass.push_back(-1);
  gen_part_E0.push_back(-1);
  gen_part_start_pos_X.push_back(-1);
  gen_part_start_pos_Y.push_back(-1);
  gen_part_start_pos_Z.push_back(-1);
  gen_part_P0_X.push_back(-1);
  gen_part_P0_Y.push_back(-1);
  gen_part_P0_Z.push_back(-1);
}


void test::AnalyzeEventsNeutrino::analize_truth_gen(const art::Ptr<simb::MCTruth> MCT, bool is_non_null, int best_track_ID)
{
  if(is_non_null && (MCT->Origin() != 0)) {
    if(MCT->Origin() == 1) {
      gen_index = 1;

      nu_PDG = MCT->GetNeutrino().Nu().PdgCode();
      nu_E0 = MCT->GetNeutrino().Nu().E(0);
      nu_weight = MCT->GetNeutrino().Nu().Weight();

      nu_interaction_mode = MCT->GetNeutrino().Mode();
      nu_interaction_type = MCT->GetNeutrino().InteractionType();
      nu_CC_NC = MCT->GetNeutrino().CCNC();

      nu_HitNuc = MCT->GetNeutrino().HitNuc();
      nu_target = MCT->GetNeutrino().Target();
      nu_HitQuark = MCT->GetNeutrino().HitQuark();
      nu_W = MCT->GetNeutrino().W();
      nu_X = MCT->GetNeutrino().X();
      nu_Y = MCT->GetNeutrino().Y();
      nu_QSqr = MCT->GetNeutrino().QSqr();

      double xyz_true_vertex[3] = {MCT->GetParticle(0).Vx(0), MCT->GetParticle(0).Vy(0), MCT->GetParticle(0).Vz(0)};
      slc_true_prim_vtx_x = xyz_true_vertex[0];
      slc_true_prim_vtx_y = xyz_true_vertex[1];
      slc_true_prim_vtx_z = xyz_true_vertex[2];

      geo::Point_t xyz_true_vertex_point_t = {xyz_true_vertex[0], xyz_true_vertex[1], xyz_true_vertex[2]};

      if (!fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_true_vertex_point_t))) {
        slc_true_prim_vtx_wires_U = -1;
        slc_true_prim_vtx_wires_V = -1;
        slc_true_prim_vtx_wires_C = -1;
        slc_true_prim_vtx_time_tick = -1;
        slc_true_prim_vtx_channel_U = -1;
        slc_true_prim_vtx_channel_V = -1;
        slc_true_prim_vtx_channel_C = -1;
      } else{
        std::vector<int> true_vertex_nu_wire = get_point_in_wires(xyz_true_vertex[0], xyz_true_vertex[1], xyz_true_vertex[2]);

        slc_true_prim_vtx_wires_U = true_vertex_nu_wire[0];
        slc_true_prim_vtx_wires_V = true_vertex_nu_wire[1];
        slc_true_prim_vtx_wires_C = true_vertex_nu_wire[2];
        slc_true_prim_vtx_time_tick = true_vertex_nu_wire[3];
        slc_true_prim_vtx_channel_U = true_vertex_nu_wire[4];
        slc_true_prim_vtx_channel_V = true_vertex_nu_wire[5];
        slc_true_prim_vtx_channel_C = true_vertex_nu_wire[6];
      }

      if (verbose_gen) {
        std::cout << "---- event " << event_ID << " ----" << std::endl;
        std::cout << "Neutrino Information: " << std::endl;
        std::cout << "Neutrino properties: " << std::endl;
        std::cout << "PDG: " << nu_PDG << " Energy: " << nu_E0 << " weight: " << nu_weight << std::endl;
        std::cout << "Interaction mode: " << nu_interaction_mode << " Interaction type: " << nu_interaction_type << " CCorNC: " << nu_CC_NC << std::endl;
        std::cout << "Target: " << nu_target << " HitNuc: " << nu_HitNuc << " HitQuark: " << nu_HitQuark << std::endl;
        std::cout << " W: " << nu_W << " X: " << nu_X << " Y: " << nu_Y << " Qsqr: " << nu_QSqr << std::endl;
        std::cout << std::endl;

        std::cout << "True vertex information: " << std::endl;
        std::cout << "   Position in cm (x,y,z): " << slc_true_prim_vtx_x << " " << slc_true_prim_vtx_y << " " << slc_true_prim_vtx_z << " "<< std::endl;
        std::cout << "   Position in wires (U,V,C, drift ticks): " << slc_true_prim_vtx_wires_U << " " << slc_true_prim_vtx_wires_V << " " << slc_true_prim_vtx_wires_C << " " << slc_true_prim_vtx_time_tick << " " << std::endl;
        std::cout << "   Position in Channels (U,V,C): " << slc_true_prim_vtx_channel_U << " " << slc_true_prim_vtx_channel_V << " " << slc_true_prim_vtx_channel_C  << std::endl;
      }

      if (save_truth_gen_particles) {
        for (int i = 0;i<MCT->NParticles();i++) {
          simb::MCParticle MCP = MCT->GetParticle(i);
          gen_part_trackID.push_back(MCP.TrackId());
          gen_part_statusCode.push_back(MCP.StatusCode());
          gen_part_mother.push_back(MCP.Mother());
          gen_part_PDGcode.push_back(MCP.PdgCode());
          gen_part_mass.push_back(MCP.Mass());
          gen_part_E0.push_back(MCP.E(0));
          gen_part_start_pos_X.push_back(MCP.Vx(0));
          gen_part_start_pos_Y.push_back(MCP.Vy(0));
          gen_part_start_pos_Z.push_back(MCP.Vz(0));
          gen_part_P0_X.push_back(MCP.Px(0));
          gen_part_P0_Y.push_back(MCP.Py(0));
          gen_part_P0_Z.push_back(MCP.Pz(0));

          if (verbose_gen) {
            std::cout << "ID: " << MCP.TrackId() << " Status code: "<< MCP.StatusCode() <<  " Mother: " << MCP.Mother()<< std::endl;
            std::cout << "   PDGCode: " << MCP.PdgCode() << " Mass: " << MCP.Mass() << " Energy: " << MCP.E(0) << std::endl;
            std::cout << "   Start Pos (x,y,z) : (" << MCP.Vx(0) << ", " << MCP.Vy(0) << ", " << MCP.Vz(0)  << ")" << std::endl;
            std::cout << "   Start Momentum (x,y,z) : (" << MCP.Px(0) << ", " << MCP.Py(0) << ", " << MCP.Pz(0)  << ")" <<std::endl;
          }
        }
      }
      if (verbose_gen) std::cout <<std::endl;
    } else if(MCT->Origin() == 2) {
      gen_index = 2;
      fill_empty_truth();

      auto const mcp = pi_serv->TrackIdToParticle_P(best_track_ID);
      double xyz_true_vertex[3] = {mcp->Vx(0), mcp->Vy(0), mcp->Vz(0)};
      slc_true_prim_vtx_x = xyz_true_vertex[0];
      slc_true_prim_vtx_y = xyz_true_vertex[1];
      slc_true_prim_vtx_z = xyz_true_vertex[2];
    }

  } else {
    gen_index = -1;
    fill_empty_truth();
  }
}


void test::AnalyzeEventsNeutrino::analize_truth_g4_particles(std::vector<art::Ptr < simb::MCParticle>> MCP_vec)
{
  //Save all the information of the g4 particles
  for (const art::Ptr <simb::MCParticle> &MCP: MCP_vec) {
    const size_t numberTrajectoryPoints = MCP->NumberTrajectoryPoints();
    const int last = numberTrajectoryPoints - 1;

    g4_part_trackID.push_back(MCP->TrackId());
    g4_part_mother.push_back(MCP->Mother());
    g4_part_PDGcode.push_back(MCP->PdgCode());
    g4_part_mass.push_back(MCP->Mass());
    g4_part_E0.push_back(MCP->E(0));
    g4_part_Ef.push_back(MCP->E(last));

    double track_lenght = 0;
    for(int i_t = 1;i_t <=last;i_t++) {
      TVector3 last_position(MCP->Vx(i_t - 1), MCP->Vy(i_t - 1), MCP->Vz(i_t - 1));
      TVector3 current_position(MCP->Vx(i_t), MCP->Vy(i_t), MCP->Vz(i_t));
      track_lenght += (last_position - current_position).Mag();
    }

    g4_part_TL.push_back(track_lenght);

    g4_part_start_pos_X.push_back(MCP->Vx(0));
    g4_part_start_pos_Y.push_back(MCP->Vy(0));
    g4_part_start_pos_Z.push_back(MCP->Vz(0));
    g4_part_start_T.push_back(MCP->T(0));
    g4_part_end_pos_X.push_back(MCP->Vx(last));
    g4_part_end_pos_Y.push_back(MCP->Vy(last));
    g4_part_end_pos_Z.push_back(MCP->Vz(last));
    g4_part_end_T.push_back(MCP->T(last));

    g4_part_P0_X.push_back(MCP->Px(0));
    g4_part_P0_Y.push_back(MCP->Py(0));
    g4_part_P0_Z.push_back(MCP->Pz(0));
    g4_part_Pf_X.push_back(MCP->Px(last));
    g4_part_Pf_Y.push_back(MCP->Py(last));
    g4_part_Pf_Z.push_back(MCP->Pz(last));

    g4_part_process.push_back(MCP->Process());
    g4_part_end_process.push_back(MCP->EndProcess());

    if (verbose_g4) {
      std::cout << "ID: " << MCP->TrackId()<< " Mother: " << MCP->Mother()<<std::endl;
      std::cout << "   PDGCode: " << MCP->PdgCode()<< " Mass: " << MCP->Mass() << " Track Lenght: " << track_lenght
                << " Initial Energy: " << MCP->E(0) << " Final Energy: " << MCP->E(last) << std::endl;
      std::cout << "   Start Pos (x,y,z,t) : (" <<  MCP->Vx(0) << ", " << MCP->Vy(0) << ", " << MCP->Vz(0) << ", " << MCP->T(0) << ")" <<std::endl;
      std::cout << "   End Pos (x,y,z,t) : (" << MCP->Vx(last)<< ", " << MCP->Vy(last)<< ", "<< MCP->Vz(last)<< ", " << MCP->T(last)<< ")" <<std::endl;
      std::cout << "   Start Momentum (x,y,z) : (" << MCP->Px(0) << ", " << MCP->Py(0) << ", "<< MCP->Pz(0) << ")" <<std::endl;
      std::cout << "   End Momentum (x,y,z) : (" << MCP->Px(last)<< ", " << MCP->Py(last)<< ", "<< MCP->Pz(last)<< ")" <<std::endl;
      std::cout << "   Process: " << MCP->Process()<<std::endl;
      std::cout << "   End_process: " << MCP->EndProcess()<< std::endl;
    }
  }

  if (verbose_g4) std::cout <<std::endl;
}

void test::AnalyzeEventsNeutrino::do_hit_matching(std::vector<art::Ptr<recob::Hit>> hit_vector, int pfp_ID) {
  for(const art::Ptr<recob::Hit> &hit : hit_vector) {
    double wire_ID = hit->WireID().deepestIndex();
    double plane_ID = hit->WireID().planeID().deepestIndex();
    double TPC_ID = hit->WireID().planeID().parentID().deepestIndex();
    double channel_ID = hit->Channel();
    double integral = hit->Integral();
    double sigma_integral = hit->SigmaIntegral();
    double peak_time = hit->PeakTime();
    double RMS = hit->RMS();

    bool hit_found = false;
    for(long unsigned int i_h = 0; i_h < slc_hits_wire_ID.size();i_h++) {
      if(hit_found) continue;
      if((wire_ID == slc_hits_wire_ID.at(i_h)) && (plane_ID==slc_hits_plane_ID.at(i_h))
      && (TPC_ID==slc_hits_TPC_ID.at(i_h)) && (channel_ID==slc_hits_channel_ID.at(i_h))
      && (integral==slc_hits_integral.at(i_h)) && (sigma_integral==slc_hits_sigma_integral.at(i_h))
      && (peak_time==slc_hits_peaktime.at(i_h)) && (RMS==slc_hits_RMS.at(i_h))) {

        hit_found = true;
        slc_hits_associated_pfp_ID.at(i_h) = pfp_ID;
      }
    }
  }
}


//Save calorimetry of the track
void test::AnalyzeEventsNeutrino::save_slc_trk_calo(std::vector<art::Ptr<anab::Calorimetry>> track_calo) {
  if(track_calo.size() == 3) {



    int num_hits[3] = {0,0,0};
    for (unsigned i = 0; i < track_calo.size(); i++) {
      const anab::Calorimetry &calo = *track_calo[i];
      if (calo.PlaneID()) {
        unsigned plane_id = calo.PlaneID().Plane;
        assert(plane_id < 3);
        for(size_t i_h = 0; i_h < calo.dEdx().size(); ++i_h) {
          if (calo.dEdx()[i_h] > 1000.) continue;
          num_hits[plane_id]++;
        }
      }
    }
    int max_hits = 0;
    for(int plane: {2, 0, 1}){
      if(num_hits[plane] > max_hits) {
        max_hits = num_hits[plane];
      }
    }

    //SAVING
    std::vector<double> dEdx_vec;
    std::vector<double> dQdx_vec;
    std::vector<double> pitch_vec;
    std::vector<double> residual_range_vec;
    std::vector<double> x_vec;
    std::vector<double> y_vec;
    std::vector<double> z_vec;


    bool first = true;
    for (unsigned i = 0; i < track_calo.size(); i++) {
      const anab::Calorimetry &calo = *track_calo[i];
      if (calo.PlaneID()) {

        int num_hits_local = 0;
        for(size_t i_h = 0; i_h < calo.dEdx().size(); ++i_h) {
          if (calo.dEdx()[i_h] > 1000.) continue;
          num_hits_local++;
        }

        if(num_hits_local != max_hits) continue;
        if(!first) continue;
        first = false;

        const std::vector<float> &dQdx = calo.dQdx();
        const std::vector<float> &dEdx = calo.dEdx();
        const std::vector<float> &pitch = calo.TrkPitchVec();
        //Visible and kinetic energy are the same depending on the calo module used
        double visible_energy = 0;
        double charge = 0.;

        for(size_t i = 0; i < dEdx.size(); ++i) {
          if (dEdx[i] > 1000.) continue;
          visible_energy     += dEdx[i] * pitch[i];
          charge += dQdx[i] * pitch[i];
        }

        slc_pfp_track_visible_energy.push_back(visible_energy);
        slc_pfp_track_kinetic_energy.push_back(calo.KineticEnergy());
        slc_pfp_track_charge.push_back(charge);

        for(long unsigned int i_d = 0; i_d < dQdx.size();i_d++) {
          dEdx_vec.push_back(calo.dEdx().at(i_d));
          dQdx_vec.push_back(calo.dQdx().at(i_d));
          pitch_vec.push_back(calo.TrkPitchVec().at(i_d));
          residual_range_vec.push_back(calo.ResidualRange().at(i_d));
          x_vec.push_back(calo.XYZ().at(i_d).X());
          y_vec.push_back(calo.XYZ().at(i_d).Y());
          z_vec.push_back(calo.XYZ().at(i_d).Z());
        }
      }
    }

    slc_pfp_calo_dEdx.push_back(dEdx_vec);
    slc_pfp_calo_dQdx.push_back(dQdx_vec);
    slc_pfp_calo_pitch.push_back(pitch_vec);
    slc_pfp_calo_residual_range.push_back(residual_range_vec);
    slc_pfp_calo_x.push_back(x_vec);
    slc_pfp_calo_y.push_back(y_vec);
    slc_pfp_calo_z.push_back(z_vec);

  } else {
    slc_pfp_track_kinetic_energy.push_back(-1);
    slc_pfp_track_visible_energy.push_back(-1);
    slc_pfp_track_charge.push_back(-1);
  }



}

std::map<int,std::map<geo::PlaneID,int> > test::AnalyzeEventsNeutrino::NumberofPlaneHitsPerTrack(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit> >& hits){

  std::map<int, std::map<geo::PlaneID, int> > HitNum;

  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    geo::PlaneID  PlaneID = hit->WireID().planeID();
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);

    //Find which track deposited the most energy.
    int likelytrack = TruthMatchUtils::TrueParticleID(clockData,hit,true);
    ++HitNum[likelytrack][PlaneID];
  }
  return HitNum;
}



void test::AnalyzeEventsNeutrino::save_slc_pfp_razzled(std::vector<art::Ptr<sbn::MVAPID> >  pfp_razzled_vec) {
  if (pfp_razzled_vec.size() != 0) {
    std::map<int, float> razzled_map = pfp_razzled_vec.at(0)->mvaScoreMap;
    slc_pfp_razzled_score_e.push_back(razzled_map[11]);
    slc_pfp_razzled_score_mu.push_back(razzled_map[13]);
    slc_pfp_razzled_score_pho.push_back(razzled_map[22]);
    slc_pfp_razzled_score_pi.push_back(razzled_map[211]);
    slc_pfp_razzled_score_p.push_back(razzled_map[2212]);

    slc_pfp_razzled_pdg.push_back(pfp_razzled_vec.at(0)->BestPDG());
  } else {
    slc_pfp_razzled_score_e.push_back(-1);
    slc_pfp_razzled_score_mu.push_back(-1);
    slc_pfp_razzled_score_pho.push_back(-1);
    slc_pfp_razzled_score_pi.push_back(-1);
    slc_pfp_razzled_score_p.push_back(-1);
    slc_pfp_razzled_pdg.push_back(-1);
  }
}

int test::AnalyzeEventsNeutrino::NumberofHitsFromTrack(const detinfo::DetectorClocksData& clockData, int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits){
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  int HitNum = 0;
  //Loop over the hits and find the IDE
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int likelytrack = TruthMatchUtils::TrueParticleID(clockData,hit,true);
    if(likelytrack == TrackID){++HitNum;}
  }
  return HitNum;
}



void test::AnalyzeEventsNeutrino::save_slc_pfp_vtx(std::vector<art::Ptr < recob::Vertex>> pfp_vertex_vec) {

  if(pfp_vertex_vec.size() != 0) {
    for (const art::Ptr <recob::Vertex> &vertex: pfp_vertex_vec) {
      double xyz_vertex[3];
      vertex->XYZ(xyz_vertex);

      slc_pfp_vtx_x.push_back(xyz_vertex[0]);
      slc_pfp_vtx_y.push_back(xyz_vertex[1]);
      slc_pfp_vtx_z.push_back(xyz_vertex[2]);

      if (!fGeom->HasTPC(fGeom->FindTPCAtPosition(vertex->position()))) {
        slc_pfp_vtx_wires_U.push_back(-1);
        slc_pfp_vtx_wires_V.push_back(-1);
        slc_pfp_vtx_wires_C.push_back(-1);
        slc_pfp_vtx_time_tick.push_back(-1);
        slc_pfp_vtx_channel_U.push_back(-1);
        slc_pfp_vtx_channel_V.push_back(-1);
        slc_pfp_vtx_channel_C.push_back(-1);
      } else {
        std::vector<int> reco_vertex_nu_wire = get_point_in_wires(xyz_vertex[0], xyz_vertex[1], xyz_vertex[2]);

        slc_pfp_vtx_wires_U.push_back(reco_vertex_nu_wire[0]);
        slc_pfp_vtx_wires_V.push_back(reco_vertex_nu_wire[1]);
        slc_pfp_vtx_wires_C.push_back(reco_vertex_nu_wire[2]);
        slc_pfp_vtx_time_tick.push_back(reco_vertex_nu_wire[3]);
        slc_pfp_vtx_channel_U.push_back(reco_vertex_nu_wire[4]);
        slc_pfp_vtx_channel_V.push_back(reco_vertex_nu_wire[5]);
        slc_pfp_vtx_channel_C.push_back(reco_vertex_nu_wire[6]);
     }


    }
  } else {
    slc_pfp_vtx_x.push_back(-1);
    slc_pfp_vtx_y.push_back(-1);
    slc_pfp_vtx_z.push_back(-1);

    slc_pfp_vtx_wires_U.push_back(-1);
    slc_pfp_vtx_wires_V.push_back(-1);
    slc_pfp_vtx_wires_C.push_back(-1);
    slc_pfp_vtx_time_tick.push_back(-1);
    slc_pfp_vtx_channel_U.push_back(-1);
    slc_pfp_vtx_channel_V.push_back(-1);
    slc_pfp_vtx_channel_C.push_back(-1);

  }
}


void test::AnalyzeEventsNeutrino::do_truth_matching(const detinfo::DetectorClocksData& clock_data_evt, const std::vector<art::Ptr<recob::Hit> >& track_hit_vec){
  int track_id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clock_data_evt, track_hit_vec, true);
  if(track_id != def_int) {
  //if(false) {
    
    const simb::MCParticle* MC_part = pi_serv->TrackIdToParticle_P(track_id);
    if(MC_part != NULL) {
      slc_pfp_true_trackid.push_back(MC_part->TrackId());
      slc_pfp_true_pdg.push_back(MC_part->PdgCode());
      double track_lenght = 0;

      const size_t numberTrajectoryPoints = MC_part->NumberTrajectoryPoints();
      const int last = numberTrajectoryPoints - 1;

      for(int i_t = 1;i_t <=last;i_t++) {
        TVector3 last_position(MC_part->Vx(i_t - 1), MC_part->Vy(i_t - 1), MC_part->Vz(i_t - 1));
        TVector3 current_position(MC_part->Vx(i_t), MC_part->Vy(i_t), MC_part->Vz(i_t));
        track_lenght += (last_position - current_position).Mag();
      }
      slc_pfp_true_TL.push_back(track_lenght);
      slc_pfp_true_energy.push_back(MC_part->E());
      slc_pfp_true_end_process.push_back(MC_part->EndProcess());
      slc_pfp_true_p_x.push_back(MC_part->Px());
      slc_pfp_true_p_y.push_back(MC_part->Py());
      slc_pfp_true_p_z.push_back(MC_part->Pz());
    } else {
      slc_pfp_true_trackid.push_back(-1);
      slc_pfp_true_pdg.push_back(-1);
      slc_pfp_true_TL.push_back(-1);
      slc_pfp_true_energy.push_back(-1);
      slc_pfp_true_end_process.push_back("null");
      slc_pfp_true_p_x.push_back(-1);
      slc_pfp_true_p_y.push_back(-1);
      slc_pfp_true_p_z.push_back(-1);
    }

    int num_hits = NumberofHitsFromTrack(clock_data_evt, track_id, track_hit_vec);
    int true_hits = 0;
    for (int i_tpc = 0; i_tpc < 2; i_tpc++) {
      for (int i_plane = 0; i_plane < 3; i_plane++) {
        geo::PlaneID plane_ID(0,i_tpc,i_plane);
        true_hits += slice_true_hit_map[track_id][plane_ID];
      }
    }
    slc_pfp_purity.push_back(num_hits*1.0/track_hit_vec.size());
    slc_pfp_completeness.push_back(num_hits*1.0/true_hits);
    
  } else {
    slc_pfp_completeness.push_back(-1);
    slc_pfp_purity.push_back(-1);
    slc_pfp_true_trackid.push_back(-1);
    slc_pfp_true_pdg.push_back(-1);
    slc_pfp_true_TL.push_back(-1);
    slc_pfp_true_energy.push_back(-1);
    slc_pfp_true_end_process.push_back("null");
    slc_pfp_true_p_x.push_back(-1);
    slc_pfp_true_p_y.push_back(-1);
    slc_pfp_true_p_z.push_back(-1);
  }

}


void test::AnalyzeEventsNeutrino::save_pfp_track_start_end(TVector3 start, TVector3 end) {
//Get track start and end
  double xyz_start[3] = {start.X(), start.Y(), start.Z()};
  geo::Point_t start_pointT = {start.X(), start.Y(), start.Z()};

  slc_pfp_trk_start_x.push_back(xyz_start[0]);
  slc_pfp_trk_start_y.push_back(xyz_start[1]);
  slc_pfp_trk_start_z.push_back(xyz_start[2]);

  double xyz_end[3] = {end.X(), end.Y(), end.Z()};
  geo::Point_t end_pointT = {end.X(), end.Y(), end.Z()};

  slc_pfp_trk_end_x.push_back(xyz_end[0]);
  slc_pfp_trk_end_y.push_back(xyz_end[1]);
  slc_pfp_trk_end_z.push_back(xyz_end[2]);
  if (fGeom->HasTPC(fGeom->FindTPCAtPosition(start_pointT))) {
    std::vector<int> reco_wire = get_point_in_wires(xyz_start[0], xyz_start[1], xyz_start[2]);
    slc_pfp_trk_start_wires_U.push_back(reco_wire[0]);
    slc_pfp_trk_start_wires_V.push_back(reco_wire[1]);
    slc_pfp_trk_start_wires_C.push_back(reco_wire[2]);
    slc_pfp_trk_start_time_tick.push_back(reco_wire[3]);
    slc_pfp_trk_start_channel_U.push_back(reco_wire[4]);
    slc_pfp_trk_start_channel_V.push_back(reco_wire[5]);
    slc_pfp_trk_start_channel_C.push_back(reco_wire[6]);
  } else {
    slc_pfp_trk_start_wires_U.push_back(-1);
    slc_pfp_trk_start_wires_V.push_back(-1);
    slc_pfp_trk_start_wires_C.push_back(-1);
    slc_pfp_trk_start_time_tick.push_back(-1);
    slc_pfp_trk_start_channel_U.push_back(-1);
    slc_pfp_trk_start_channel_V.push_back(-1);
    slc_pfp_trk_start_channel_C.push_back(-1);
  }

  if (fGeom->HasTPC(fGeom->FindTPCAtPosition(end_pointT))) {
    std::vector<int> reco_wire = get_point_in_wires(xyz_end[0], xyz_end[1], xyz_end[2]);
    slc_pfp_trk_end_wires_U.push_back(reco_wire[0]);
    slc_pfp_trk_end_wires_V.push_back(reco_wire[1]);
    slc_pfp_trk_end_wires_C.push_back(reco_wire[2]);
    slc_pfp_trk_end_time_tick.push_back(reco_wire[3]);
    slc_pfp_trk_end_channel_U.push_back(reco_wire[4]);
    slc_pfp_trk_end_channel_V.push_back(reco_wire[5]);
    slc_pfp_trk_end_channel_C.push_back(reco_wire[6]);
  } else {
    slc_pfp_trk_end_wires_U.push_back(-1);
    slc_pfp_trk_end_wires_V.push_back(-1);
    slc_pfp_trk_end_wires_C.push_back(-1);
    slc_pfp_trk_end_time_tick.push_back(-1);
    slc_pfp_trk_end_channel_U.push_back(-1);
    slc_pfp_trk_end_channel_V.push_back(-1);
    slc_pfp_trk_end_channel_C.push_back(-1);
  }
}

void  test::AnalyzeEventsNeutrino::FillMomDiff(const art::Ptr<sbn::RangeP> &rangeP_muon, const art::Ptr<sbn::RangeP> &rangeP_pion, const art::Ptr<recob::MCSFitResult> &mcs_muon, const art::Ptr<recob::MCSFitResult> &mcs_pion)
{
  const float rangeMom = rangeP_muon->range_p;
  const float mcsMom   = mcs_muon->fwdMomentum();


  slc_pfp_trk_range_mom.push_back(rangeP_muon->range_p);
  slc_pfp_trk_mcs_mom.push_back(mcs_muon->fwdMomentum());

  slc_pfp_trk_range_mom_muon.push_back(rangeP_muon->range_p);
  slc_pfp_trk_range_mom_pion.push_back(rangeP_pion->range_p);
  slc_pfp_trk_mcs_mom_muon.push_back(mcs_muon->fwdMomentum());
  slc_pfp_trk_mcs_mom_pion.push_back(mcs_pion->fwdMomentum());

  slc_pfp_trk_mom_diff.push_back( (rangeMom > 0 && mcsMom > 0) ? (mcsMom - rangeMom) / rangeMom : -10.f);
}

void test::AnalyzeEventsNeutrino::FillStoppingChi2Metrics(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2)
{
  const float pol0Chi2 = stoppingChi2->pol0Chi2;
  const float expChi2  = stoppingChi2->expChi2;

  slc_pfp_trk_stopping_dEdx_chi2_ratio.push_back((pol0Chi2 > 0.f && expChi2 > 0.f) ? pol0Chi2 / expChi2 : -5.f);
  slc_pfp_trk_chi2_pol0_dEdx_fit.push_back(stoppingChi2->pol0Fit);
}
void  test::AnalyzeEventsNeutrino::FillMCSMetrics(const art::Ptr<recob::MCSFitResult> &mcs)
{
  if(mcs->scatterAngles().empty()) {
    slc_pfp_trk_mcs_scatter_mean.push_back(-1);
    slc_pfp_trk_mcs_scatter_max_ratio.push_back(-1);
    return;
  }

  unsigned int counter = 0;
  float maxScatter = 0.f, meanScatter = 0.f;

  for(auto const& angle : mcs->scatterAngles())
  {
    if(angle < 0)
      continue;

    maxScatter = std::max(maxScatter, angle);
    meanScatter += angle;
    counter++;
  }

  if(!counter) {
    slc_pfp_trk_mcs_scatter_mean.push_back(-1);
    slc_pfp_trk_mcs_scatter_max_ratio.push_back(-1);
    return;
  }

  slc_pfp_trk_mcs_scatter_mean.push_back(meanScatter / counter);
  slc_pfp_trk_mcs_scatter_max_ratio.push_back(maxScatter / meanScatter);
}


void test::AnalyzeEventsNeutrino::save_track_info(std::vector<art::Ptr<recob::Track>> pfp_track_vec, std::vector<art::Ptr<recob::Shower>> pfp_shower_vec, art::FindManyP<anab::Calorimetry> track_calo_assns, art::FindManyP<anab::ParticleID> track_chi2_assn, art::FindManyP<recob::MCSFitResult> tracksToMCSs_muon, art::FindManyP<recob::MCSFitResult> tracksToMCSs_pion , art::FindManyP<sbn::RangeP> tracksToRangePs_muon, art::FindManyP<sbn::RangeP> tracksToRangePs_pion, art::FindManyP<sbn::ScatterClosestApproach> tracksToClosestApproaches, art::FindManyP<sbn::StoppingChi2Fit> tracksToStoppingChi2s) {
    if(pfp_track_vec.size() != 0) {
      slc_pfp_has_track.push_back(true);

      std::vector<art::Ptr<anab::Calorimetry>> track_calo = track_calo_assns.at(pfp_track_vec.front().key());

      //Some basic track variables
      art::Ptr<recob::Track> track = pfp_track_vec[0];
      geo::Vector_t start_dir = track->StartDirection();

      slc_pfp_track_dir_x.push_back(start_dir.X());
      slc_pfp_track_dir_y.push_back(start_dir.Y());
      slc_pfp_track_dir_z.push_back(start_dir.Z());


      slc_pfp_track_theta.push_back(track->Theta() * 180 / TMath::Pi());
      slc_pfp_track_phi.push_back(track->Phi() * 180 / TMath::Pi());
      slc_pfp_track_length.push_back(track->Length());

      TVector3 start = {track->Start().X(), track->Start().Y(), track->Start().Z()};
      TVector3 end = {track->End().X(), track->End().Y(), track->End().Z()};
      save_pfp_track_start_end(start, end);

      //Save calorimetry of the track
      save_slc_trk_calo(track_calo);

      int best_plane = 0;
      int num_hits[3] = {0,0,0};
      if(track_calo.size() == 3) {
        for (unsigned i = 0; i < track_calo.size(); i++) {
          const anab::Calorimetry &calo = *track_calo[i];
          if (calo.PlaneID()) {
            unsigned plane_id = calo.PlaneID().Plane;
            assert(plane_id < 3);
            for (size_t i_h = 0; i_h < calo.dEdx().size(); ++i_h) {
              if (calo.dEdx()[i_h] > 1000.) continue;
              num_hits[plane_id]++;
            }
         //   std::cout << calo.PlaneID() << " " << num_hits[plane_id] << std::endl;
          }
        }

        int bestnhit = -1;
        for(int plane: {2, 0, 1}){
          if(num_hits[plane] > bestnhit){
            best_plane = plane;
            bestnhit = num_hits[plane];
          }
        }
      }

      //std::cout << best_plane << " " << num_hits[0] << "  " <<  num_hits[1] << "  " <<  num_hits[2] << std::endl;




    const std::vector<art::Ptr<recob::MCSFitResult>> mcs_muon = tracksToMCSs_muon.at(track.key());
    const std::vector<art::Ptr<recob::MCSFitResult>> mcs_pion = tracksToMCSs_pion.at(track.key());
    //std::cout << "msc size:" << mcs.size() << std::endl;
    if(mcs_muon.size() != 0) {
      this->FillMCSMetrics(mcs_muon[0]);
    } else {
      slc_pfp_trk_mcs_scatter_mean.push_back(-1);
      slc_pfp_trk_mcs_scatter_max_ratio.push_back(-1);
    }

    const std::vector<art::Ptr<sbn::RangeP>> rangeP_muon = tracksToRangePs_muon.at(track.key());
    const std::vector<art::Ptr<sbn::RangeP>> rangeP_pion = tracksToRangePs_pion.at(track.key());
    //std::cout << "rangeP size:" << rangeP.size() << std::endl;
    if(rangeP_muon.size() != 0 && rangeP_pion.size() != 0 && mcs_muon.size() != 0 && mcs_pion.size() != 0) {
      this->FillMomDiff(rangeP_muon[0], rangeP_pion[0], mcs_muon[0], mcs_pion[0]);
    } else {
      slc_pfp_trk_mom_diff.push_back(-1);
      slc_pfp_trk_range_mom.push_back(-1);
      slc_pfp_trk_mcs_mom.push_back(-1);
      slc_pfp_trk_range_mom_muon.push_back(-1);
      slc_pfp_trk_range_mom_pion.push_back(-1);
      slc_pfp_trk_mcs_mom_muon.push_back(-1);
      slc_pfp_trk_mcs_mom_pion.push_back(-1);
    }

     std::cout << "TL: " << track->Length() << " TS: " << "not_included" << std::endl;


    const std::vector<art::Ptr<sbn::ScatterClosestApproach>> closestApproach = tracksToClosestApproaches.at(track.key());
//std::cout << "closestApproach size:" << closestApproach.size() << std::endl;
    if(closestApproach.size() != 0) {
      slc_pfp_trk_meanDCA.push_back(closestApproach[0]->mean);
    } else {
      slc_pfp_trk_meanDCA.push_back(-1);
    }

    const std::vector<art::Ptr<sbn::StoppingChi2Fit>> stoppingChi2 = tracksToStoppingChi2s.at(track.key());
//std::cout << "StoppingChi2Fit size:" << mcs.size() << std::endl;
    if(stoppingChi2.size() != 0) {
      this->FillStoppingChi2Metrics(stoppingChi2[0]);
    } else {
      slc_pfp_trk_stopping_dEdx_chi2_ratio.push_back(-1);
      slc_pfp_trk_chi2_pol0_dEdx_fit.push_back(-1);
    }

    //save Chi2
    bool chi2_selected = false;

    //Get Chi2
    const std::vector<art::Ptr<anab::ParticleID>> chi2_vec = track_chi2_assn.at(track.key());

     best_plane = 1;
      if(track_calo.size() == 3) {
        for (unsigned i = 0; i < chi2_vec.size(); i++) {
          const anab::ParticleID &chi2_pid = *chi2_vec[i];
          if(chi2_pid.PlaneID()) {
            int plane_id  = chi2_pid.PlaneID().Plane;
            assert(plane_id < 3);
            if(plane_id == best_plane) {

              std::vector<anab::sParticleIDAlgScores> alg_score_vector = chi2_pid.ParticleIDAlgScores();
              for (size_t i_as = 0; i_as < alg_score_vector.size(); i_as++) {
                const anab::sParticleIDAlgScores alg_score = alg_score_vector.at(i_as);
                if (alg_score.fAlgName == "Chi2") {
                  chi2_selected = true;
                  if (TMath::Abs(alg_score.fAssumedPdg) == 13) { // chi2mu
                    std::cout<< " chi2mu: " <<alg_score.fValue;
                    slc_pfp_track_chi2_muon.push_back(alg_score.fValue);
                  } else if (TMath::Abs(alg_score.fAssumedPdg) == 211) { // chi2pi
                    slc_pfp_track_chi2_pion.push_back(alg_score.fValue);
                  } else if (TMath::Abs(alg_score.fAssumedPdg) == 321) { // chi2ka
                    slc_pfp_track_chi2_kaon.push_back(alg_score.fValue);
                  } else if (TMath::Abs(alg_score.fAssumedPdg) == 2212) { // chi2pr
                    std::cout<< " chi2p: " <<alg_score.fValue;
                    slc_pfp_track_chi2_proton.push_back(alg_score.fValue);
                  }
                }
              }

            }
          }
        }
      }


    if(!chi2_selected) {
      slc_pfp_track_chi2_muon.push_back(-1);
      slc_pfp_track_chi2_pion.push_back(-1);
      slc_pfp_track_chi2_kaon.push_back(-1);
      slc_pfp_track_chi2_proton.push_back(-1);
    }
  } else if(pfp_shower_vec.size() != 0) {
    slc_pfp_has_track.push_back(false);

    slc_pfp_track_theta.push_back(-1);
    slc_pfp_track_phi.push_back(-1);

    TVector3 start = pfp_shower_vec.front()->ShowerStart();
    TVector3 end = TVector3(-1, -1, -1);
    save_pfp_track_start_end(start, end);

    slc_pfp_track_dir_x.push_back(pfp_shower_vec.front()->Direction().x());
    slc_pfp_track_dir_y.push_back(pfp_shower_vec.front()->Direction().y());
    slc_pfp_track_dir_z.push_back(pfp_shower_vec.front()->Direction().z());

    slc_pfp_track_length.push_back(pfp_shower_vec.front()->Length());

    slc_pfp_track_kinetic_energy.push_back(pfp_shower_vec.front()->Energy().at(pfp_shower_vec.front()->best_plane()));
    slc_pfp_track_visible_energy.push_back(-1);
    slc_pfp_track_charge.push_back(-1);

    slc_pfp_calo_dEdx.push_back({-1});
    slc_pfp_calo_dQdx.push_back({-1});
    slc_pfp_calo_pitch.push_back({-1});
    slc_pfp_calo_residual_range.push_back({-1});
    slc_pfp_calo_x.push_back({-1});
    slc_pfp_calo_y.push_back({-1});
    slc_pfp_calo_z.push_back({-1});

    slc_pfp_track_chi2_muon.push_back(-1);
    slc_pfp_track_chi2_pion.push_back(-1);
    slc_pfp_track_chi2_kaon.push_back(-1);
    slc_pfp_track_chi2_proton.push_back(-1);

    slc_pfp_trk_meanDCA.push_back(-1);
    slc_pfp_trk_range_mom.push_back(-1);
    slc_pfp_trk_mcs_mom.push_back(-1);
    slc_pfp_trk_range_mom_muon.push_back(-1);
    slc_pfp_trk_range_mom_pion.push_back(-1);
    slc_pfp_trk_mcs_mom_muon.push_back(-1);
    slc_pfp_trk_mcs_mom_pion.push_back(-1);
    slc_pfp_trk_mom_diff.push_back(-1);
    slc_pfp_trk_stopping_dEdx_chi2_ratio.push_back(-1);
    slc_pfp_trk_chi2_pol0_dEdx_fit.push_back(-1);
    slc_pfp_trk_mcs_scatter_mean.push_back(-1);
    slc_pfp_trk_mcs_scatter_max_ratio.push_back(-1);

  } else {
    fill_empty_slice_pfp_track_info();
  }
}






void test::AnalyzeEventsNeutrino::save_shower_info(std::vector<art::Ptr<recob::Track>> pfp_track_vec, std::vector<art::Ptr<recob::Shower>> pfp_shower_vec) {

  if(pfp_shower_vec.size() != 0) {

    art::Ptr<recob::Shower> shower = pfp_shower_vec.front();
    TVector3 start_dir = shower->Direction();

    slc_pfp_shower_dir_y.push_back(start_dir.Y());
    slc_pfp_shower_dir_x.push_back(start_dir.X());
    slc_pfp_shower_dir_z.push_back(start_dir.Z());

    slc_pfp_shower_length.push_back(shower->Length());
    slc_pfp_shower_opening_angle.push_back(shower->OpenAngle());

    int best_plane = shower->best_plane();
    slc_pfp_shower_energy.push_back(shower->Energy()[best_plane]);
    slc_pfp_shower_dEdx.push_back(shower->dEdx()[best_plane]);
  } else {
    fill_empty_slice_pfp_shower_info();
  }
}

void test::AnalyzeEventsNeutrino::analize_slc_pfps(art::Event const &e, art::Handle< std::vector<recob::PFParticle> > pfp_handle, std::vector< art::Ptr<recob::PFParticle> > pfp_vec, art::Handle< std::vector<recob::Track> > track_handle, art::Handle< std::vector<recob::Shower> > shower_handle, art::FindManyP<recob::Hit> sp_hit_assns)
{
  auto const clock_data_evt = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  //Loop trhough all the particles
  std::cout << "pfps:" << std::endl;
  for(const art::Ptr <recob::PFParticle> &pfp: pfp_vec) {
    //Load associations
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfp_handle, e, pfp_to_pfpmetadata_label); //pfp asociation to MetaData
    art::FindManyP<sbn::MVAPID> pfp_razzled_assns (pfp_handle, e, pfp_to_razzled_label); // pfp association to razzled
    art::FindManyP<recob::Track> pfp_track_assns (pfp_handle, e, pfp_to_track_label); //pfp association to tracks
    art::FindManyP<recob::Shower> pfp_shower_assns (pfp_handle, e, pfp_to_shower_label); //pfp association to showers
    art::FindManyP<recob::Vertex> pfp_vertex_assns (pfp_handle, e, pfp_to_vertex_label); //pfparticle asociation to vertex
    art::FindManyP<recob::SpacePoint> pfp_sp_assns (pfp_handle, e, pfp_to_space_point_label); //pfparticle asociation to space points

    art::FindManyP<anab::Calorimetry> track_calo_assns (track_handle, e, track_to_calo_label); //track association to calo
    art::FindManyP<recob::Hit> track_hit_assns(track_handle, e, track_to_hits_label); // track association to hits
    art::FindManyP<recob::Hit> shower_hit_assns(shower_handle, e, shower_to_hits_label); //shower association to hits
    art::FindManyP<anab::ParticleID> track_chi2_assn(track_handle, e, track_to_chi2_label); // Track association to Chi2

    art::FindManyP<recob::MCSFitResult> tracksToMCSs_muon(track_handle, e, fMCSLabel_muon);
    art::FindManyP<recob::MCSFitResult> tracksToMCSs_pion(track_handle, e, fMCSLabel_pion);

    art::FindManyP<sbn::RangeP> tracksToRangePs_muon(track_handle, e, fRangeLabel_muon);
    art::FindManyP<sbn::RangeP> tracksToRangePs_pion(track_handle, e, fRangeLabel_pion);
    art::FindManyP<sbn::ScatterClosestApproach> tracksToClosestApproaches(track_handle, e, fClosestApproachLabel);
    art::FindManyP<sbn::StoppingChi2Fit> tracksToStoppingChi2s(track_handle, e, fStoppingChi2Label);

    //Get the track score
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta_vec = pfp_pfpmeta_assns.at(pfp.key());

    //If track score is not found skip the particle
    std::map<std::string, float> pfp_prop_map = pfp_meta_vec.at(0)->GetPropertiesMap();
    if((pfp_prop_map.find("TrackScore") == pfp_prop_map.end()) || (pfp_prop_map["TrackScore"] < 0) ) continue;

    //slc_pfp_track_score.push_back( pfp_prop_map["TrackScore"]);
    FillPandoraTrackIDScoreVars(pfp_prop_map);


    //Get the pfp ID
    slc_pfp_ID.push_back(pfp->Self());

    //Get the pfp parent
    slc_pfp_parent_ID.push_back(pfp->Parent());

    //Get the pfp pandora PDG
    slc_pfp_pandora_PDG.push_back(pfp->PdgCode());

    //Check if it is daughter of primary particle
    bool is_primary = false;
    for(const art::Ptr<recob::PFParticle> &pfp_2 : pfp_vec) {
        if(is_primary) continue;
        if((pfp_2->Self() == pfp->Parent()) && (pfp_2->IsPrimary())) {
          is_primary = true;
        }
    }
    slc_pfp_mother_is_primary.push_back(is_primary);

    //Save daughter max hits
    double max_hits = 0;
    for(auto const& daughterID : pfp->Daughters()) {
      int hits = 0;

      for(const art::Ptr <recob::PFParticle> &daughter_pfp: pfp_vec) {
        if(daughter_pfp->Self() != daughterID) continue;

        std::vector<art::Ptr<recob::Track>> daughter_tracks = pfp_track_assns.at(daughter_pfp.key());
        std::vector<art::Ptr<recob::Hit>> track_hit_vec = track_hit_assns.at(daughter_tracks.front().key());

        hits += track_hit_vec.size();
      }
    if(hits > max_hits) max_hits = hits;
    }
    slc_pfp_daughter_max_hits.push_back(max_hits);

    //Save razzled information
    if(save_razzled) {
      art::FindManyP<sbn::MVAPID> pfp_razzled_assns (pfp_handle, e, pfp_to_razzled_label); // pfp association to razzled
      std::vector<art::Ptr<sbn::MVAPID> > pfp_razzled_vec = pfp_razzled_assns.at(pfp.key());
      save_slc_pfp_razzled(pfp_razzled_vec);
    }

    //NOT ADDED
    //Get vertex of the particle
    std::vector<art::Ptr < recob::Vertex>> pfp_vertex_vec = pfp_vertex_assns.at(pfp.key());
    save_slc_pfp_vtx(pfp_vertex_vec);

    //START the track and shower analysis
    std::vector<art::Ptr<recob::Track>> pfp_track_vec = pfp_track_assns.at(pfp.key());
    std::vector<art::Ptr<recob::Shower>> pfp_shower_vec = pfp_shower_assns.at(pfp.key());

    //Do Space point matching of the pfp
    if(save_reco_space_points) {
      std::vector<art::Ptr<recob::SpacePoint> > pfp_sp_vec = pfp_sp_assns.at(pfp.key());

      for(art::Ptr<recob::SpacePoint> pfp_sp: pfp_sp_vec) {
        slc_sp_x.push_back(pfp_sp->XYZ()[0]);
        slc_sp_y.push_back(pfp_sp->XYZ()[1]);
        slc_sp_z.push_back(pfp_sp->XYZ()[2]);
        slc_sp_chi2.push_back(pfp_sp->Chisq());

        std::vector<art::Ptr < recob::Hit>> sp_hit_vec = sp_hit_assns.at(pfp_sp.key());
        double integral = 0;
        double sigma_integral = 0;
        for (const art::Ptr <recob::Hit> &sp_hit: sp_hit_vec) {
          integral += sp_hit->Integral();
          sigma_integral += pow(sp_hit->SigmaIntegral(), 2);
        }
        slc_sp_integral.push_back(integral);
        slc_sp_sigma_integral.push_back(sqrt(sigma_integral));
        slc_sp_associated_pfp_ID.push_back(pfp->Self());
        //do_sp_matching(pfp_sp_vec, pfp->Self());
      }
    }
    std::cout << "pfpID: " << pfp->Self();
    //save track and shower info
    save_track_info(pfp_track_vec, pfp_shower_vec, track_calo_assns, track_chi2_assn, tracksToMCSs_muon, tracksToMCSs_pion , tracksToRangePs_muon, tracksToRangePs_pion, tracksToClosestApproaches, tracksToStoppingChi2s);
    save_shower_info(pfp_track_vec, pfp_shower_vec);

    //Do hit matchign and truth matching
    if(pfp_track_vec.size() != 0) {
      std::vector<art::Ptr<recob::Hit>> track_hit_vec = track_hit_assns.at(pfp_track_vec.front().key());

      //Save hits of the track
      if(save_reco_hits) {
        do_hit_matching(track_hit_vec, pfp->Self());
      }
      //Do the truth matching
      do_truth_matching(clock_data_evt, track_hit_vec);

    } else if(pfp_shower_vec.size() != 0) {
      std::vector<art::Ptr<recob::Hit>> shower_hit_vec = shower_hit_assns.at(pfp_shower_vec.front().key());

      //Save hits of the shower
      if(save_reco_hits) {
        do_hit_matching(shower_hit_vec, pfp->Self());
      }

      //Do the truth matching
      do_truth_matching(clock_data_evt, shower_hit_vec);
    } else {
      slc_pfp_true_trackid.push_back(-1);
      slc_pfp_true_pdg.push_back(-1);
      slc_pfp_true_TL.push_back(-1);
      slc_pfp_true_energy.push_back(-1);
      slc_pfp_true_end_process.push_back("null");
      slc_pfp_true_p_x.push_back(-1);
      slc_pfp_true_p_y.push_back(-1);
      slc_pfp_true_p_z.push_back(-1);
      slc_pfp_completeness.push_back(-1);
      slc_pfp_purity.push_back(-1);
    }
  }
}


  //Save the generator information (MCTruth) using truth matxhing to get the MCTruth for the slice
void test::AnalyzeEventsNeutrino::analyze_slice_truth(art::Event const& e, std::vector<art::Ptr<recob::Hit>> slice_hit_vec,   std::map<const art::Ptr<simb::MCTruth>, int> MCTruth_hit_map)
{
  auto const clock_data_evt = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);


  //Do the truth matching of the slice (to know if it is cosmic or other things)
  //Get a map of what true tracks the hits correspond to

  std::map<int, int> map_trackId_hits;
  for(auto const &hit : slice_hit_vec)
    map_trackId_hits[TruthMatchUtils::TrueParticleID(clock_data_evt, hit, true)]++;
  
  //Create a map of the number of hits contributed to each MCTruth
  std::map<const art::Ptr<simb::MCTruth>, int> map_mc_hits;
  int max_hits = def_int;
  int best_track_ID = def_int;

  for(auto const& [track_ID, nhits] : map_trackId_hits) {
    art::Ptr<simb::MCTruth> mct;
    try {
      mct = track_ID == def_int ? art::Ptr<simb::MCTruth>() : pi_serv->TrackIdToMCTruth_P(track_ID);
    }
    catch(...) {
      continue;
    }

    map_mc_hits[mct] += nhits;
    if(nhits > max_hits){
      best_track_ID = track_ID;
    }
  }

  max_hits = def_int;
  art::Ptr<simb::MCTruth> best_MCT = art::Ptr<simb::MCTruth>();
  for(auto const& [mct, nhits] : map_mc_hits)
  {
    //if(mct.isNonnull()) std::cout << int(mct->Origin()) << " " << nhits << std::endl;
    if(nhits > max_hits) {
      max_hits = nhits;
      best_MCT  = mct;
    }
  }
   
  slc_truth_matching_completeness = MCTruth_hit_map[best_MCT] == 0.0 ? -1.0 : map_mc_hits[best_MCT] / static_cast<float>(MCTruth_hit_map[best_MCT]);
  slc_truth_matching_purity  = slice_hit_vec.size() == 0.0 ? -1.0 : map_mc_hits[best_MCT] / static_cast<float>(slice_hit_vec.size());

  analize_truth_gen(best_MCT, best_MCT.isNonnull(), best_track_ID);
}


//Bulk of the analizer
void test::AnalyzeEventsNeutrino::analyze(art::Event const &e) {

  //Fill the event ID related information
  event_ID = e.id().event();
  run_ID = e.id().run();
  subrun_ID = e.id().subRun();
  //Overload the subrun tree information
  num_gen_evts = get_total_gen_evts(e);
  if (verbose_general) {
    std::cout << std::endl;
    std::cout << "Run number: " << run_ID << std::endl;
    std::cout << "Subrun number: " << subrun_ID << std::endl;
    std::cout << "Event number: " << event_ID << std::endl;
  }
  std::cout << "EvtID: " << event_ID << " Run: " << run_ID << " Subrun: "<< subrun_ID << std::endl;


  if(file_type == "single") {
    reset_gen_vars();

    art::Handle <std::vector<simb::MCTruth>> MCT_handle;
    std::vector<art::Ptr < simb::MCTruth>> MCT_vec;
    if (e.getByLabel(MCTruth_label, MCT_handle))
      art::fill_ptr_vector(MCT_vec, MCT_handle);

    analize_truth_gen(MCT_vec.front(), true, 0);
  }

  //Save the geant4 information (MCParticles)
  if (save_truth_g4) {
    //Reset all variables concerning the g4 stage
    reset_g4_vars();

    art::Handle <std::vector<simb::MCParticle>> MCP_handle;
    std::vector<art::Ptr < simb::MCParticle>> MCP_vec;
    if (e.getByLabel(MCParticle_label, MCP_handle))
      art::fill_ptr_vector(MCP_vec, MCP_handle);

    analize_truth_g4_particles(MCP_vec);
  }




  if (save_reco) {
    //Access to Slice
    art::Handle <std::vector<recob::Slice>> slice_handle;
    std::vector<art::Ptr < recob::Slice> > slice_vec;
    if (e.getByLabel(slice_label, slice_handle))
      art::fill_ptr_vector(slice_vec, slice_handle);

    //Access to pfpParticles
    art::Handle< std::vector<recob::PFParticle> > pfp_handle;
    std::vector< art::Ptr<recob::PFParticle> > pfp_vec;
    if(e.getByLabel(pfp_label, pfp_handle))
      art::fill_ptr_vector(pfp_vec, pfp_handle);

    //Access to tracks
    art::Handle< std::vector<recob::Track> > track_handle;
    std::vector< art::Ptr<recob::Track> > track_vec;
    if(e.getByLabel(track_label, track_handle))
      art::fill_ptr_vector(track_vec, track_handle);

    //Access to showers
    art::Handle< std::vector<recob::Shower> > shower_handle;
    std::vector< art::Ptr<recob::Shower> > shower_vec;
    if(e.getByLabel(shower_label, shower_handle))
      art::fill_ptr_vector(shower_vec, shower_handle);

    //Access to hits
    art::Handle< std::vector<recob::Hit> > hit_handle;
    std::vector< art::Ptr<recob::Hit> > hit_vec;
    if(e.getByLabel(hit_label, hit_handle))
      art::fill_ptr_vector(hit_vec, hit_handle);

    //save the Truth Hit map for slice-matching
    auto const clock_data_evt = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    std::map<const art::Ptr<simb::MCTruth>, int> MCTruth_hit_map; //mctruth id, nHits


    for(auto const& hit : hit_vec) {
      const int track_ID = TruthMatchUtils::TrueParticleID(clock_data_evt, hit, true);
      art::Ptr<simb::MCTruth> mct;
      try {
        mct = track_ID == def_int ? art::Ptr<simb::MCTruth>() : pi_serv->TrackIdToMCTruth_P(track_ID);
      }
      catch(...) {
        continue;
      }
      MCTruth_hit_map[mct]++;
    }




    //Access to space points
    art::Handle< std::vector<recob::SpacePoint> > sp_handle;
    std::vector< art::Ptr<recob::SpacePoint> > sp_vec;
    if(e.getByLabel(space_point_label, sp_handle))
      art::fill_ptr_vector(sp_vec, sp_handle);

    //Get all the associations
    art::FindManyP<recob::Hit> slice_hit_assns(slice_handle, e, slice_to_hit_label); //Slice associations to hits particles
    art::FindManyP<recob::PFParticle> slice_pfp_assns(slice_handle, e, slice_to_pfp_label); //Slice associations to pfp particles
    art::FindManyP<sbn::OpT0Finder> slice_T0_assns(slice_handle, e, slice_to_opt0_label); //Slice associations to T0

    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfp_handle, e, pfp_to_pfpmetadata_label); //pfp asociation to MetaData
    art::FindManyP<recob::Vertex> pfp_vertex_assns (pfp_handle, e, pfp_to_vertex_label); //pfparticle asociation to vertex

    art::FindManyP<recob::Hit> sp_hit_assns(sp_handle, e, space_point_to_hit_label); //Space Point associations to hit



    //Reset all variables concerning the generator stage
    if (slice_vec.size() == 0) {
      reset_slc_vars();

      slc_is_obvious_cosmic = false;
      gen_index = -1;
      fill_empty_truth();
      fill_empty_slice();
      if(!do_precuts) fTree->Fill();
    }

    for (const art::Ptr <recob::Slice> &slice: slice_vec) {
      std::cout << std::endl << "SliceID: " <<  slice->ID() << std::endl;
      //Reset slice variables
      reset_slc_vars();

      //Get Slice ID
      slc_ID = slice->ID();

      //Get pfparticles and hits
      std::vector<art::Ptr<recob::PFParticle>> slice_pfp_vec = slice_pfp_assns.at(slice.key());
      std::vector<art::Ptr<recob::Hit>> slice_hit_vec = slice_hit_assns.at(slice.key());

      //Check if it is an obvious cosmic
      bool slc_has_neutrino_vertex = true;
      //Find neutrino vertex, and save it (if we dont find it we stop the analysis here)
      slc_is_obvious_cosmic = slc_check_clear_cosmic(slice_pfp_vec, pfp_pfpmeta_assns);
      slc_has_neutrino_vertex = get_slice_nu_vertex(e, slice_pfp_vec, pfp_vertex_assns);

      if(verbose_reco) {
        std::cout << "Slice: " << slc_ID << std::endl;
        std::cout << "Is Obvious Cosmic: " << slc_is_obvious_cosmic << std::endl;
        std::cout << "Has Neutrino vertex: " << slc_has_neutrino_vertex << std::endl;
        if(slc_has_neutrino_vertex) {
          std::cout << "Neutrino Vertex Info: " << std::endl;
          std::cout << "   Position in cm (x,y,z): " << slc_prim_vtx_x << " " << slc_prim_vtx_y << " " << slc_prim_vtx_z << " "<< std::endl;
          std::cout << "   Position in wires (U,V,C, drift ticks): " << slc_prim_vtx_wires_U << " " << slc_prim_vtx_wires_V << " " << slc_prim_vtx_wires_C << " " << slc_prim_vtx_time_tick << " " << std::endl;
          std::cout << "   Position in Channels (U,V,C): " << slc_prim_vtx_channel_U << " " << slc_prim_vtx_channel_V << " " << slc_prim_vtx_channel_C  << std::endl;
        }
      }

      //Save the truth of the interaction
      if(file_type == "spill") {
        if(slice_hit_vec.size() != 0) {
           analyze_slice_truth(e, slice_hit_vec, MCTruth_hit_map);
        } else {
          gen_index = -1;
          fill_empty_truth();
        }
      }

    

      if((slc_is_obvious_cosmic) || ((slice_hit_vec.size() == 0)) || (!slc_has_neutrino_vertex)) {
        fill_empty_slice();
      } else {
        //Get the general Slice characteristics: completeness, purity, nu_score, crumvs_score and if it has a neutrino vertex (get the vertex information if so)
        //Fill the truth matching of the slice (Single interaction has no truth matching)
        slice_true_hit_map = NumberofPlaneHitsPerTrack(clock_data_evt, slice_hit_vec); //Returns a map of all the number of hits and the respetive track id they are associated to.

        //Fill crumbs score of the slice
        if(save_crumbs) {
          art::FindManyP<sbn::CRUMBSResult> slice_crumbs_assns(slice_handle, e, slice_to_crumbs_label); // Slice association with crumbs
          std::vector<art::Ptr<sbn::CRUMBSResult>> crumbs_result_vec = slice_crumbs_assns.at(slice.key());
          if(crumbs_result_vec.size() == 0) {
            slc_crumbs_score = -1;
          } else {
            slc_crumbs_score = crumbs_result_vec.at(0)->score;
          }
        }

        //Fill the nu_score of the slice
        slc_nu_score = get_slice_nu_score(slice_pfp_vec, pfp_pfpmeta_assns);
        if(verbose_reco) {
          std::cout << "Truth matching completeness: " << slc_truth_matching_completeness << " Truth matching purity: " << slc_truth_matching_purity << std::endl;
          std::cout << "Crumbs score: " << slc_crumbs_score << " Nu score: " << slc_nu_score<< std::endl;
        }

        //Get the timing information
        std::vector<art::Ptr<sbn::OpT0Finder>> slice_op_T0 = slice_T0_assns.at(slice.key());
        std::sort(slice_op_T0.begin(), slice_op_T0.end(), [](auto const &a, auto const &b){return a->score > b->score;});
        if(slice_op_T0.size() == 0) {
          slc_t0 = -1;
          slc_t0_score = -1;
          slc_OpT0Fraction = -1;
        } else {
          slc_t0 = slice_op_T0[0]->time;
          slc_t0_score = slice_op_T0[0]->score;
          slc_OpT0Fraction = (slice_op_T0[0]->measPE - slice_op_T0[0]->hypoPE)/slice_op_T0[0]->measPE;
        }

        if(save_3d_barycenter_matching) {
          //Access to BM
          art::FindManyP<sbn::TPCPMTBarycenterMatch> slice_BC_assns(slice_handle, e, slice_to_barycenter_matching_label); //Slice associations to BM
          std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> slice_BM = slice_BC_assns.at(slice.key());

          if(slice_BM.size() != 0) {
            slc_bcmatch_flash_time = slice_BM[0]->flashTime;
            slc_bcmatch_first_hit = slice_BM[0]->flashFirstHit ;
            slc_bcmatch_flash_PE = slice_BM[0]->flashPEs;
            slc_bcmatch_flash_asymmetry =slice_BM[0]->flashAsymmetry;
          } else {
            slc_bcmatch_flash_time = -10000;
            slc_bcmatch_first_hit = -10000;
            slc_bcmatch_flash_PE = -10000;
            slc_bcmatch_flash_asymmetry = -10000;
          }
        }

        //Fill hits of the interaction
        if(save_reco_hits) {
          for (const art::Ptr <recob::Hit> &slice_hit: slice_hit_vec) {
            slc_hits_wire_ID.push_back(slice_hit->WireID().deepestIndex());
            slc_hits_plane_ID.push_back(slice_hit->WireID().planeID().deepestIndex());
            slc_hits_TPC_ID.push_back(slice_hit->WireID().planeID().parentID().deepestIndex());
            slc_hits_channel_ID.push_back(slice_hit->Channel());
            slc_hits_integral.push_back(slice_hit->Integral());
            slc_hits_sigma_integral.push_back(slice_hit->SigmaIntegral());
            slc_hits_peaktime.push_back(slice_hit->PeakTime());
            slc_hits_RMS.push_back(slice_hit->RMS());
            slc_hits_associated_pfp_ID.push_back(-1);
          }
        }

        if(save_reco_stubs) {
          art::FindManyP<sbn::Stub> slice_stub_assns(slice_handle, e, slice_to_stub_label); //Slice associations to stubs
          std::vector<art::Ptr<sbn::Stub>> slice_stub_vec = slice_stub_assns.at(slice.key());
          for (const art::Ptr <sbn::Stub> &slice_stub: slice_stub_vec) {
            slc_stub_vtx_x.push_back(slice_stub->vtx.x());
            slc_stub_vtx_y.push_back(slice_stub->vtx.y());
            slc_stub_vtx_z.push_back(slice_stub->vtx.z());
            slc_stub_end_x.push_back(slice_stub->end.x());
            slc_stub_end_y.push_back(slice_stub->end.y());
            slc_stub_end_z.push_back(slice_stub->end.z());
            std::vector<double> hit_charge_vec;
            std::vector<double> hit_wire_vec;
            std::vector<int> ontrack_vec;

            for (unsigned ih = 0; ih < slice_stub->hits[0].size(); ih++) {
              hit_charge_vec.push_back(slice_stub->hits[0][ih].charge);
              hit_wire_vec.push_back(slice_stub->hits[0][ih].wire);
              ontrack_vec.push_back(slice_stub->hits[0][ih].ontrack);
            }

            slc_stub_hit_charge.push_back(hit_charge_vec);
            slc_stub_hit_wire.push_back(hit_wire_vec);
            slc_stub_hit_ontrack.push_back(ontrack_vec);
          }
        }

        //Analize PFPs
        analize_slc_pfps(e, pfp_handle, slice_pfp_vec, track_handle, shower_handle, sp_hit_assns);

        if(verbose_reco_hits) {
          std::cout << "Hit information: "<< std::endl;
          for(long unsigned int i_h = 0; i_h < slc_hits_wire_ID.size(); i_h++) {
            std::cout <<"Hit number: " << i_h << " Wire ID: " << slc_hits_wire_ID.at(i_h) << " Plane ID: " << slc_hits_plane_ID.at(i_h) << " TPC ID: " << slc_hits_TPC_ID.at(i_h) <<  " Channel ID: " << slc_hits_channel_ID.at(i_h) << std::endl;
            std::cout << "Integral: " << slc_hits_integral.at(i_h) << " Sigma Integral: " << slc_hits_sigma_integral.at(i_h) << std::endl;
            std::cout << "Peak Time: " << slc_hits_peaktime.at(i_h) << " Time RMS: " << slc_hits_RMS.at(i_h) << std::endl;
            std::cout << "Associated pfp: " << slc_hits_associated_pfp_ID.at(i_h) << std::endl << std::endl;
          }
        }

        if(verbose_reco_space_points) {
          std::cout << "Space points information: "<< std::endl;
          for(long unsigned int i_sp = 0; i_sp < slc_sp_x.size(); i_sp++) {
            std::cout <<"Space point number: " << i_sp << std::endl;
            std::cout << "Position (X,Y,Z): " << slc_sp_x.at(i_sp) << " " << slc_sp_y.at(i_sp) << " " << slc_sp_z.at(i_sp) <<  std::endl;
            std::cout << "Integral: " << slc_sp_integral.at(i_sp) << " Sigma Integral: " << slc_sp_sigma_integral.at(i_sp) << std::endl;
            std::cout << "Associated pfp: " << slc_sp_associated_pfp_ID.at(i_sp) << " Chi2: " << slc_sp_chi2.at(i_sp) <<  std::endl << std::endl;
          }
        }

        if(verbose_reco_pfps) {
          for(long unsigned int i_p = 0; i_p < slc_pfp_ID.size(); i_p++) {
            std::cout <<"PFP ID: " << slc_pfp_ID.at(i_p) << " PFP parent ID: " << slc_pfp_parent_ID.at(i_p) << " Is mother primary: " << slc_pfp_mother_is_primary.at(i_p) << std::endl;
            std::cout << "Track Score: " << slc_pfp_track_score.at(i_p) << " razzled PDG: " << slc_pfp_razzled_pdg.at(i_p) << std::endl;
            std::cout << "razzled e score: " << slc_pfp_razzled_score_e.at(i_p) << " razzled pho score: " << slc_pfp_razzled_score_pho.at(i_p) << std::endl;
            std::cout << " razzled mu score: " << slc_pfp_razzled_score_mu.at(i_p) << " razzled pi score: " << slc_pfp_razzled_score_pi.at(i_p) << " razzled proton score: " << slc_pfp_razzled_score_p.at(i_p) << std::endl ;
            std::cout << "Has track: " << slc_pfp_has_track.at(i_p) << " Pandora PDG: " << slc_pfp_pandora_PDG.at(i_p) <<std::endl;
            std::cout << "pfp Vertex Info: " << std::endl;
            std::cout << "   Position in cm (x,y,z): " << slc_pfp_vtx_x.at(i_p) << " " << slc_pfp_vtx_y.at(i_p) << " " << slc_pfp_vtx_z.at(i_p) << " "<< std::endl;
            std::cout << "   Position in wires (U,V,C, drift ticks): " << slc_pfp_vtx_wires_U.at(i_p) << " " << slc_pfp_vtx_wires_V.at(i_p) << " " << slc_pfp_vtx_wires_C.at(i_p) << " " << slc_pfp_vtx_time_tick.at(i_p) << " " << std::endl;
            std::cout << "   Position in Channels (U,V,C): " << slc_pfp_vtx_channel_U.at(i_p) << " " << slc_pfp_vtx_channel_V.at(i_p) << " " << slc_pfp_vtx_channel_C.at(i_p)  << std::endl;

            std::cout << "Truth matching info: " << std::endl;
            std::cout <<"True ID: " << slc_pfp_true_trackid.at(i_p) << " True PDG: " << slc_pfp_true_pdg.at(i_p) << std::endl;
            std::cout <<"Purity: " << slc_pfp_purity.at(i_p) << " Completeness: " << slc_pfp_completeness.at(i_p) << std::endl;
            std::cout <<"True Energy: " << slc_pfp_purity.at(i_p) << std::endl;
            std::cout << "True momentum (x,y,z): " << slc_pfp_true_p_x.at(i_p) << " " << slc_pfp_true_p_y.at(i_p) << " " << slc_pfp_true_p_z.at(i_p)  << std::endl;

            fTree->Branch("slc_pfp_track_length", &slc_pfp_track_length);
            std::cout << "Calorimetry info: " << std::endl;
            std::cout << "Track Lenght: " << slc_pfp_track_length.at(i_p) << std::endl;
            std::cout << "Kinetic Energy: " << slc_pfp_track_kinetic_energy.at(i_p) <<" Visible Energy: " << slc_pfp_track_visible_energy.at(i_p) <<" Charge: " << slc_pfp_track_charge.at(i_p) << std::endl;
            if(verbose_reco_calo_deposition) {
              for(long unsigned int i_d = 0; i_d < slc_pfp_calo_dEdx.at(i_p).size(); i_d++) {
                std::cout <<"Hit number: " << i_d << " dEdx: " << slc_pfp_calo_dEdx.at(i_p).at(i_d) << " dQdx: " << slc_pfp_calo_dQdx.at(i_p).at(i_d) << " pitch: " << slc_pfp_calo_pitch.at(i_p).at(i_d) << std::endl;
                std::cout << "Residual range: " << slc_pfp_calo_residual_range.at(i_p).at(i_d) << std::endl;
                std::cout << "Position (x,y,z):" << slc_pfp_calo_x.at(i_p).at(i_d) << " " << slc_pfp_calo_y.at(i_p).at(i_d) << " "<< slc_pfp_calo_z.at(i_p).at(i_d) << std::endl;
              }
            }
            std::cout << "Reco direction (x,y,z): " << slc_pfp_track_dir_x.at(i_p) << " " << slc_pfp_track_dir_y.at(i_p) << " " << slc_pfp_track_dir_z.at(i_p)  << std::endl;
            std::cout << "Reco direction (theta,phi): " << slc_pfp_track_theta.at(i_p) << " " << slc_pfp_track_phi.at(i_p) << std::endl;

            std::cout << "pfp Start Info: " << std::endl;
            std::cout << "   Position in cm (x,y,z): " << slc_pfp_trk_start_x.at(i_p) << " " << slc_pfp_trk_start_y.at(i_p) << " " << slc_pfp_trk_start_z.at(i_p) << " "<< std::endl;
            std::cout << "   Position in wires (U,V,C, drift ticks): " << slc_pfp_trk_start_wires_U.at(i_p) << " " << slc_pfp_trk_start_wires_V.at(i_p) << " " << slc_pfp_trk_start_wires_C.at(i_p) << " " << slc_pfp_trk_start_time_tick.at(i_p) << " " << std::endl;
            std::cout << "   Position in Channels (U,V,C): " << slc_pfp_trk_start_channel_U.at(i_p) << " " << slc_pfp_trk_start_channel_V.at(i_p) << " " << slc_pfp_trk_start_channel_C.at(i_p)  << std::endl;

            std::cout << "pfp End Info: " << std::endl;
            std::cout << "   Position in cm (x,y,z): " << slc_pfp_trk_end_x.at(i_p) << " " << slc_pfp_trk_end_y.at(i_p) << " " << slc_pfp_trk_end_z.at(i_p) << " "<< std::endl;
            std::cout << "   Position in wires (U,V,C, drift ticks): " << slc_pfp_trk_end_wires_U.at(i_p) << " " << slc_pfp_trk_end_wires_V.at(i_p) << " " << slc_pfp_trk_end_wires_C.at(i_p) << " " << slc_pfp_trk_end_time_tick.at(i_p) << " " << std::endl;
            std::cout << "   Position in Channels (U,V,C): " << slc_pfp_trk_end_channel_U.at(i_p) << " " << slc_pfp_trk_end_channel_V.at(i_p) << " " << slc_pfp_trk_end_channel_C.at(i_p)  << std::endl;

            std::cout << "Chi2 Info: " << std::endl;
            std::cout <<"Muon: " << slc_pfp_track_chi2_muon.at(i_p) << " Pion: " << slc_pfp_track_chi2_pion.at(i_p) << std::endl;
            std::cout <<"Kaon: " << slc_pfp_track_chi2_kaon.at(i_p) << " Proton: " << slc_pfp_track_chi2_proton.at(i_p) << std::endl;

            std::cout << "Shower Info: " << std::endl;
            std::cout << " Direction (x,y,z): " << slc_pfp_shower_dir_x.at(i_p) << " " << slc_pfp_shower_dir_y.at(i_p) << " " << slc_pfp_shower_dir_z.at(i_p) << std::endl;
            std::cout << " Length: " << slc_pfp_shower_length.at(i_p) << " Opening angle: " << slc_pfp_shower_opening_angle.at(i_p) << std::endl;
            std::cout << " Energy: " << slc_pfp_shower_energy.at(i_p) << " dEdx: " << slc_pfp_shower_dEdx.at(i_p) << std::endl;

            std::cout << std::endl;

          }
        }

      }

      if(verbose_reco) {
        std::cout << "Is reconstructed: " << slc_is_reconstructed << std::endl;
        std::cout << std::endl;
      }




      if(do_precuts) {
        if((slc_is_reconstructed) && (!slc_is_obvious_cosmic)) fTree->Fill();
      } else if(!slc_is_obvious_cosmic){
         fTree->Fill();
      }
    }
  } else {
    if(do_precuts) {
      if((slc_is_reconstructed) && (!slc_is_obvious_cosmic)) fTree->Fill();
    } else if(!slc_is_obvious_cosmic) {
      fTree->Fill();
    }
  }

}


//Change the name of the file saved on the tree when a new file is opened
void test::AnalyzeEventsNeutrino::respondToOpenInputFile(art::FileBlock const &fb) {
  file_name = fb.fileName();
}

//Create a tree when the code is initialized
void test::AnalyzeEventsNeutrino::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle <art::TFileService> tfs;

  //Create a tree to save the subrun information
  fSubRunTree = tfs->make<TTree>("subrun_tree", "SubRun Tree");
  fSubRunTree->Branch("run_ID", &run_ID);
  fSubRunTree->Branch("subrun_ID", &subrun_ID);
  fSubRunTree->Branch("POT", &POT);
  fSubRunTree->Branch("spill", &spill);
  fSubRunTree->Branch("num_gen_evts", &num_gen_evts);

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Add branches to the TTree
  fTree->Branch("file_name", &file_name);
  fTree->Branch("event_ID", &event_ID);
  fTree->Branch("run_ID", &run_ID);
  fTree->Branch("subrun_ID", &subrun_ID);

  fTree->Branch("gen_index", &gen_index);
  fTree->Branch("nu_PDG", &nu_PDG);
  fTree->Branch("nu_E0", &nu_E0);
  fTree->Branch("nu_weight", &nu_weight);
  fTree->Branch("nu_interaction_mode", &nu_interaction_mode);
  fTree->Branch("nu_interaction_type", &nu_interaction_type);
  fTree->Branch("nu_CC_NC", &nu_CC_NC);
  fTree->Branch("nu_target", &nu_target);
  fTree->Branch("nu_HitNuc", &nu_HitNuc);
  fTree->Branch("nu_HitQuark", &nu_HitQuark);
  fTree->Branch("nu_W", &nu_W);
  fTree->Branch("nu_X", &nu_X);
  fTree->Branch("nu_Y", &nu_Y);
  fTree->Branch("nu_QSqr", &nu_QSqr);

  if (save_truth_gen_particles) {
    //MCTruth Particles variables
    fTree->Branch("gen_part_trackID", &gen_part_trackID);
    fTree->Branch("gen_part_statusCode", &gen_part_statusCode);
    fTree->Branch("gen_part_mother", &gen_part_mother);
    fTree->Branch("gen_part_PDGcode", &gen_part_PDGcode);
    fTree->Branch("gen_part_mass", &gen_part_mass);
    fTree->Branch("gen_part_E0", &gen_part_E0);
    fTree->Branch("gen_part_start_pos_X", &gen_part_start_pos_X);
    fTree->Branch("gen_part_start_pos_Y", &gen_part_start_pos_Y);
    fTree->Branch("gen_part_start_pos_Z", &gen_part_start_pos_Z);
    fTree->Branch("gen_part_P0_X", &gen_part_P0_X);
    fTree->Branch("gen_part_P0_Y", &gen_part_P0_Y);
    fTree->Branch("gen_part_P0_Z", &gen_part_P0_Z);
  }

  if (save_truth_g4) {
    //Geant4 Particles variables
    fTree->Branch("g4_part_trackID", &g4_part_trackID);
    fTree->Branch("g4_part_mother", &g4_part_mother);
    fTree->Branch("g4_part_PDGcode", &g4_part_PDGcode);
    fTree->Branch("g4_part_mass", &g4_part_mass);
    fTree->Branch("g4_part_E0", &g4_part_E0);
    fTree->Branch("g4_part_Ef", &g4_part_Ef);
    fTree->Branch("g4_part_TL", &g4_part_TL);

    fTree->Branch("g4_part_start_pos_X", &g4_part_start_pos_X);
    fTree->Branch("g4_part_start_pos_Y", &g4_part_start_pos_Y);
    fTree->Branch("g4_part_start_pos_Z", &g4_part_start_pos_Z);
    fTree->Branch("g4_part_start_T", &g4_part_start_T);
    fTree->Branch("g4_part_end_pos_X", &g4_part_end_pos_X);
    fTree->Branch("g4_part_end_pos_Y", &g4_part_end_pos_Y);
    fTree->Branch("g4_part_end_pos_Z", &g4_part_end_pos_Z);
    fTree->Branch("g4_part_end_T", &g4_part_end_T);

    fTree->Branch("g4_part_P0_X", &g4_part_P0_X);
    fTree->Branch("g4_part_P0_Y", &g4_part_P0_Y);
    fTree->Branch("g4_part_P0_Z", &g4_part_P0_Z);
    fTree->Branch("g4_part_Pf_X", &g4_part_Pf_X);
    fTree->Branch("g4_part_Pf_Y", &g4_part_Pf_Y);
    fTree->Branch("g4_part_Pf_Z", &g4_part_Pf_Z);

    fTree->Branch("g4_part_process", &g4_part_process);
    fTree->Branch("g4_part_end_process", &g4_part_end_process);
  }

  //Slice realted branches
  if (save_reco) {
    fTree->Branch("slc_ID", &slc_ID);
    fTree->Branch("slc_is_reconstructed", &slc_is_reconstructed);
    fTree->Branch("slc_is_obvious_cosmic", &slc_is_obvious_cosmic);

    fTree->Branch("slc_prim_vtx_x", &slc_prim_vtx_x);
    fTree->Branch("slc_prim_vtx_y", &slc_prim_vtx_y);
    fTree->Branch("slc_prim_vtx_z", &slc_prim_vtx_z);
    if (save_reco_hits) {
      fTree->Branch("slc_prim_vtx_wires_U", &slc_prim_vtx_wires_U);
      fTree->Branch("slc_prim_vtx_wires_V", &slc_prim_vtx_wires_V);
      fTree->Branch("slc_prim_vtx_wires_C", &slc_prim_vtx_wires_C);
      fTree->Branch("slc_prim_vtx_time_tick", &slc_prim_vtx_time_tick);
      fTree->Branch("slc_prim_vtx_channel_U", &slc_prim_vtx_channel_U);
      fTree->Branch("slc_prim_vtx_channel_V", &slc_prim_vtx_channel_V);
      fTree->Branch("slc_prim_vtx_channel_C", &slc_prim_vtx_channel_C);
    }

    fTree->Branch("slc_true_prim_vtx_x", &slc_true_prim_vtx_x);
    fTree->Branch("slc_true_prim_vtx_y", &slc_true_prim_vtx_y);
    fTree->Branch("slc_true_prim_vtx_z", &slc_true_prim_vtx_z);
    if (save_reco_hits) {
      fTree->Branch("slc_true_prim_vtx_wires_U", &slc_true_prim_vtx_wires_U);
      fTree->Branch("slc_true_prim_vtx_wires_V", &slc_true_prim_vtx_wires_V);
      fTree->Branch("slc_true_prim_vtx_wires_C", &slc_true_prim_vtx_wires_C);
      fTree->Branch("slc_true_prim_vtx_time_tick", &slc_true_prim_vtx_time_tick);
      fTree->Branch("slc_true_prim_vtx_channel_U", &slc_true_prim_vtx_channel_U);
      fTree->Branch("slc_true_prim_vtx_channel_V", &slc_true_prim_vtx_channel_V);
      fTree->Branch("slc_true_prim_vtx_channel_C", &slc_true_prim_vtx_channel_C);
    }

    fTree->Branch("slc_truth_matching_completeness", &slc_truth_matching_completeness);
    fTree->Branch("slc_truth_matching_purity", &slc_truth_matching_purity);
    fTree->Branch("slc_nu_score", &slc_nu_score);
    if(save_crumbs) fTree->Branch("slc_crumbs_score", &slc_crumbs_score);

    fTree->Branch("slc_t0", &slc_t0);
    fTree->Branch("slc_t0_score", &slc_t0_score);
    fTree->Branch("slc_OpT0Fraction", &slc_OpT0Fraction);

    if(save_3d_barycenter_matching) {
      fTree->Branch("slc_bcmatch_flash_time", &slc_bcmatch_flash_time);
      fTree->Branch("slc_bcmatch_first_hit", &slc_bcmatch_first_hit);
      fTree->Branch("slc_bcmatch_flash_PE", &slc_bcmatch_flash_PE);
      fTree->Branch("slc_bcmatch_flash_asymmetry", &slc_bcmatch_flash_asymmetry);
    }

    if(save_reco_space_points) {
      fTree->Branch("slc_sp_x", &slc_sp_x);
      fTree->Branch("slc_sp_y", &slc_sp_y);
      fTree->Branch("slc_sp_z", &slc_sp_z);
      fTree->Branch("slc_sp_chi2", &slc_sp_chi2);
      fTree->Branch("slc_sp_integral", &slc_sp_integral);
      fTree->Branch("slc_sp_sigma_integral", &slc_sp_sigma_integral);
      fTree->Branch("slc_sp_associated_pfp_ID", &slc_sp_associated_pfp_ID);
    }

    if(save_reco_stubs) {
      fTree->Branch("slc_stub_vtx_x", &slc_stub_vtx_x);
      fTree->Branch("slc_stub_vtx_y", &slc_stub_vtx_y);
      fTree->Branch("slc_stub_vtx_z", &slc_stub_vtx_z);
      fTree->Branch("slc_stub_end_x", &slc_stub_end_x);
      fTree->Branch("slc_stub_end_y", &slc_stub_end_y);
      fTree->Branch("slc_stub_end_z", &slc_stub_end_z);
      fTree->Branch("slc_stub_hit_charge", &slc_stub_hit_charge);
      fTree->Branch("slc_stub_hit_wire", &slc_stub_hit_wire);
      fTree->Branch("slc_stub_hit_ontrack", &slc_stub_hit_ontrack);
    }

    if (save_reco_hits) {
      fTree->Branch("slc_hits_wire_ID", &slc_hits_wire_ID);
      fTree->Branch("slc_hits_plane_ID", &slc_hits_plane_ID);
      fTree->Branch("slc_hits_TPC_ID", &slc_hits_TPC_ID);
      fTree->Branch("slc_hits_channel_ID", &slc_hits_channel_ID);
      fTree->Branch("slc_hits_integral", &slc_hits_integral);
      fTree->Branch("slc_hits_sigma_integral", &slc_hits_sigma_integral);
      fTree->Branch("slc_hits_peaktime", &slc_hits_peaktime);
      fTree->Branch("slc_hits_RMS", &slc_hits_RMS);
      fTree->Branch("slc_hits_associated_pfp_ID", &slc_hits_associated_pfp_ID);
    }

    fTree->Branch("slc_pfp_ID", &slc_pfp_ID);
    fTree->Branch("slc_pfp_daughter_max_hits", &slc_pfp_daughter_max_hits);
    fTree->Branch("slc_pfp_has_track", &slc_pfp_has_track);
    fTree->Branch("slc_pfp_pandora_PDG", &slc_pfp_pandora_PDG);
    fTree->Branch("slc_pfp_parent_ID", &slc_pfp_parent_ID);
    fTree->Branch("slc_pfp_mother_is_primary", &slc_pfp_mother_is_primary);
    fTree->Branch("slc_pfp_track_score", &slc_pfp_track_score);

    fTree->Branch("slc_pfp_charge_end_frac", &slc_pfp_charge_end_frac);
    fTree->Branch("slc_pfp_charge_frac_spread", &slc_pfp_charge_frac_spread);
    fTree->Branch("slc_pfp_linear_fit_diff", &slc_pfp_linear_fit_diff);
    fTree->Branch("slc_pfp_linear_fit_length", &slc_pfp_linear_fit_length);
    fTree->Branch("slc_pfp_linear_fit_gap_length", &slc_pfp_linear_fit_gap_length);
    fTree->Branch("slc_pfp_linear_fit_RMS", &slc_pfp_linear_fit_RMS);
    fTree->Branch("slc_pfp_open_angle_diff", &slc_pfp_open_angle_diff);
    fTree->Branch("slc_pfp_secondary_PCA_ratio", &slc_pfp_secondary_PCA_ratio);
    fTree->Branch("slc_pfp_tertiary_PCA_ratio", &slc_pfp_tertiary_PCA_ratio);
    fTree->Branch("slc_pfp_vertex_dist", &slc_pfp_vertex_dist);

    if(save_razzled) {
      fTree->Branch("slc_pfp_razzled_score_e", &slc_pfp_razzled_score_e);
      fTree->Branch("slc_pfp_razzled_score_mu", &slc_pfp_razzled_score_mu);
      fTree->Branch("slc_pfp_razzled_score_pho", &slc_pfp_razzled_score_pho);
      fTree->Branch("slc_pfp_razzled_score_pi", &slc_pfp_razzled_score_pi);
      fTree->Branch("slc_pfp_razzled_score_p", &slc_pfp_razzled_score_p);
      fTree->Branch("slc_pfp_razzled_pdg", &slc_pfp_razzled_pdg);
    }

    fTree->Branch("slc_pfp_vtx_x", &slc_pfp_vtx_x);
    fTree->Branch("slc_pfp_vtx_y", &slc_pfp_vtx_y);
    fTree->Branch("slc_pfp_vtx_z", &slc_pfp_vtx_z);
    if (save_reco_hits) {
      fTree->Branch("slc_pfp_vtx_wires_U", &slc_pfp_vtx_wires_U);
      fTree->Branch("slc_pfp_vtx_wires_V", &slc_pfp_vtx_wires_V);
      fTree->Branch("slc_pfp_vtx_wires_C", &slc_pfp_vtx_wires_C);
      fTree->Branch("slc_pfp_vtx_time_tick", &slc_pfp_vtx_time_tick);
      fTree->Branch("slc_pfp_vtx_channel_U", &slc_pfp_vtx_channel_U);
      fTree->Branch("slc_pfp_vtx_channel_V", &slc_pfp_vtx_channel_V);
      fTree->Branch("slc_pfp_vtx_channel_C", &slc_pfp_vtx_channel_C);
    }

    fTree->Branch("slc_pfp_true_trackid", &slc_pfp_true_trackid);
    fTree->Branch("slc_pfp_true_pdg", &slc_pfp_true_pdg);
    fTree->Branch("slc_pfp_true_TL", &slc_pfp_true_TL);
    fTree->Branch("slc_pfp_true_energy", &slc_pfp_true_energy);
    fTree->Branch("slc_pfp_true_end_process", &slc_pfp_true_end_process);
    fTree->Branch("slc_pfp_true_p_x", &slc_pfp_true_p_x);
    fTree->Branch("slc_pfp_true_p_y", &slc_pfp_true_p_y);
    fTree->Branch("slc_pfp_true_p_z", &slc_pfp_true_p_z);
    fTree->Branch("slc_pfp_completeness", &slc_pfp_completeness);
    fTree->Branch("slc_pfp_purity", &slc_pfp_purity);

    fTree->Branch("slc_pfp_trk_start_x", &slc_pfp_trk_start_x);
    fTree->Branch("slc_pfp_trk_start_y", &slc_pfp_trk_start_y);
    fTree->Branch("slc_pfp_trk_start_z", &slc_pfp_trk_start_z);
    if (save_reco_hits) {
      fTree->Branch("slc_pfp_trk_start_wires_U", &slc_pfp_trk_start_wires_U);
      fTree->Branch("slc_pfp_trk_start_wires_V", &slc_pfp_trk_start_wires_V);
      fTree->Branch("slc_pfp_trk_start_wires_C", &slc_pfp_trk_start_wires_C);
      fTree->Branch("slc_pfp_trk_start_time_tick", &slc_pfp_trk_start_time_tick);
      fTree->Branch("slc_pfp_trk_start_channel_U", &slc_pfp_trk_start_channel_U);
      fTree->Branch("slc_pfp_trk_start_channel_V", &slc_pfp_trk_start_channel_V);
      fTree->Branch("slc_pfp_trk_start_channel_C", &slc_pfp_trk_start_channel_C);
    }

    fTree->Branch("slc_pfp_trk_end_x", &slc_pfp_trk_end_x);
    fTree->Branch("slc_pfp_trk_end_y", &slc_pfp_trk_end_y);
    fTree->Branch("slc_pfp_trk_end_z", &slc_pfp_trk_end_z);
    if (save_reco_hits) {
      fTree->Branch("slc_pfp_trk_end_wires_U", &slc_pfp_trk_end_wires_U);
      fTree->Branch("slc_pfp_trk_end_wires_V", &slc_pfp_trk_end_wires_V);
      fTree->Branch("slc_pfp_trk_end_wires_C", &slc_pfp_trk_end_wires_C);
      fTree->Branch("slc_pfp_trk_end_time_tick", &slc_pfp_trk_end_time_tick);
      fTree->Branch("slc_pfp_trk_end_channel_U", &slc_pfp_trk_end_channel_U);
      fTree->Branch("slc_pfp_trk_end_channel_V", &slc_pfp_trk_end_channel_V);
      fTree->Branch("slc_pfp_trk_end_channel_C", &slc_pfp_trk_end_channel_C);
    }

    fTree->Branch("slc_pfp_trk_mcs_scatter_mean", &slc_pfp_trk_mcs_scatter_mean);
    fTree->Branch("slc_pfp_trk_mcs_scatter_max_ratio", &slc_pfp_trk_mcs_scatter_max_ratio);

    fTree->Branch("slc_pfp_trk_meanDCA", &slc_pfp_trk_meanDCA);

    fTree->Branch("slc_pfp_trk_stopping_dEdx_chi2_ratio", &slc_pfp_trk_stopping_dEdx_chi2_ratio);
    fTree->Branch("slc_pfp_trk_chi2_pol0_dEdx_fit", &slc_pfp_trk_chi2_pol0_dEdx_fit);

    fTree->Branch("slc_pfp_trk_range_mom", &slc_pfp_trk_range_mom);
    fTree->Branch("slc_pfp_trk_mcs_mom", &slc_pfp_trk_mcs_mom);

    fTree->Branch("slc_pfp_trk_range_mom_muon", &slc_pfp_trk_range_mom_muon);
    fTree->Branch("slc_pfp_trk_range_mom_pion", &slc_pfp_trk_range_mom_pion);
    fTree->Branch("slc_pfp_trk_mcs_mom_muon", &slc_pfp_trk_mcs_mom_muon);
    fTree->Branch("slc_pfp_trk_mcs_mom_pion", &slc_pfp_trk_mcs_mom_pion);

    fTree->Branch("slc_pfp_trk_mom_diff", &slc_pfp_trk_mom_diff);

    fTree->Branch("slc_pfp_track_dir_x", &slc_pfp_track_dir_x);
    fTree->Branch("slc_pfp_track_dir_y", &slc_pfp_track_dir_y);
    fTree->Branch("slc_pfp_track_dir_z", &slc_pfp_track_dir_z);
    fTree->Branch("slc_pfp_track_dir_vec_x", &slc_pfp_track_dir_vec_x);
    fTree->Branch("slc_pfp_track_dir_vec_y", &slc_pfp_track_dir_vec_y);
    fTree->Branch("slc_pfp_track_dir_vec_z", &slc_pfp_track_dir_vec_z);
    fTree->Branch("slc_pfp_track_length", &slc_pfp_track_length);
    fTree->Branch("slc_pfp_track_kinetic_energy", &slc_pfp_track_kinetic_energy);
    fTree->Branch("slc_pfp_track_visible_energy", &slc_pfp_track_visible_energy);
    fTree->Branch("slc_pfp_track_charge", &slc_pfp_track_charge);

    fTree->Branch("slc_pfp_track_chi2_muon", &slc_pfp_track_chi2_muon);
    fTree->Branch("slc_pfp_track_chi2_pion", &slc_pfp_track_chi2_pion);
    fTree->Branch("slc_pfp_track_chi2_kaon", &slc_pfp_track_chi2_kaon);
    fTree->Branch("slc_pfp_track_chi2_proton", &slc_pfp_track_chi2_proton);
    fTree->Branch("slc_pfp_track_theta", &slc_pfp_track_theta);
    fTree->Branch("slc_pfp_track_phi", &slc_pfp_track_phi);

    fTree->Branch("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
    fTree->Branch("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
    fTree->Branch("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
    fTree->Branch("slc_pfp_shower_length", &slc_pfp_shower_length);
    fTree->Branch("slc_pfp_shower_opening_angle", &slc_pfp_shower_opening_angle);
    fTree->Branch("slc_pfp_shower_energy", &slc_pfp_shower_energy);
    fTree->Branch("slc_pfp_shower_dEdx", &slc_pfp_shower_dEdx);

    if(save_track_calo_deposition) {
      fTree->Branch("slc_pfp_calo_dEdx", &slc_pfp_calo_dEdx);
      fTree->Branch("slc_pfp_calo_dQdx", &slc_pfp_calo_dQdx);
      fTree->Branch("slc_pfp_calo_pitch", &slc_pfp_calo_pitch);
      fTree->Branch("slc_pfp_calo_residual_range", &slc_pfp_calo_residual_range);
      fTree->Branch("slc_pfp_calo_x", &slc_pfp_calo_x);
      fTree->Branch("slc_pfp_calo_y", &slc_pfp_calo_y);
      fTree->Branch("slc_pfp_calo_z", &slc_pfp_calo_z);
    }
  }
}


void test::AnalyzeEventsNeutrino::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEventsNeutrino)