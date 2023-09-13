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

#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"

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

  void ResetVars();
  void ClearMaps();
  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                 const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  void SetupBranches(VecVarMap &map);

  void ResizeVectors(VecVarMap &map, const int size);
  void ResizeSubVectors(VecVarMap &map, const int pos, const int size);

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
                     const art::Handle<std::vector<recob::Track>> &trackHandle,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle);

  void AnalyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track, const int slcCounter, const int pfpCounter,
                    const art::Handle<std::vector<recob::Track>> &trackHandle);
  void AnalyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle, const art::Ptr<recob::Vertex> &vtx,
                     const std::vector<art::Ptr<recob::Hit>> &hits);

  void ExtractDazzle(const art::Ptr<sbn::MVAPID> &dazzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<anab::Calorimetry> &calo, const int slcCounter, const int pfpCounter);
  void ExtractChi2PID(const art::Ptr<anab::ParticleID> &chi2pid, const int slcCounter, const int pfpCounter);
  void ExtractMCS(const art::Ptr<recob::MCSFitResult> &mcs, const int slcCounter, const int pfpCounter);
  void ExtractStoppingChi2(const art::Ptr<sbn::StoppingChi2Fit> &stoppingChi2, const int slcCounter, const int pfpCounter);

  void ExtractRazzle(const art::Ptr<sbn::MVAPID> &razzle, const int slcCounter, const int pfpCounter);
  void ExtractCalo(const art::Ptr<recob::Shower> &shower, const int slcCounter, const int pfpCounter,
                   const std::vector<art::Ptr<recob::Hit>> &hits);

  float Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  float Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  template<typename T>
  void FillElement(VecVar *vec, const int pos, const T value);
  template<typename T>
  void FillElement(VecVar *vec, const int posA, const int posB, const T value);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  art::InputTag fMCParticleModuleLabel, fSliceModuleLabel, fPFParticleModuleLabel, fVertexModuleLabel,
    fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fTrackCalorimetryModuleLabel,
    fCRUMBSModuleLabel, fDazzleModuleLabel, fCaloModuleLabel, fMCSModuleLabel, fChi2ModuleLabel, fRangeModuleLabel,
    fClosestApproachModuleLabel, fStoppingChi2ModuleLabel, fRazzleModuleLabel, fCosmicDistModuleLabel,
    fShowerTrackFitModuleLabel, fShowerDensityFitModuleLabel;
  bool fDebug;

  std::map<int,int> fHitsMap;
  std::map<int, art::Ptr<recob::PFParticle>> fPFPMap;
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;

  TTree* fEventTree;

  // Tree variables

  int  _run;
  int  _subrun;
  int  _event;

  int _n_nu;

  VecVarMap nuVars = {
    { "nu_event_type", new InhVecVar<int>("nu_event_type") },
    { "nu_signal", new InhVecVar<bool>("nu_signal") },
    { "nu_true_en_dep", new InhVecVar<float>("nu_true_en_dep") },
    { "nu_ccnc", new InhVecVar<int>("nu_ccnc") },
    { "nu_mode", new InhVecVar<int>("nu_mode") },
    { "nu_int_type", new InhVecVar<int>("nu_int_type") },
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
  };

  int _n_slc;

  VecVarMap slcVars = {
    { "slc_n_pfps", new InhVecVar<size_t>("slc_n_pfps") },
    { "slc_primary_pfp_id", new InhVecVar<size_t>("slc_primary_pfp_id") },
    { "slc_primary_pfp_pdg", new InhVecVar<int>("slc_primary_pfp_pdg") },
    { "slc_is_clear_cosmic", new InhVecVar<bool>("slc_is_clear_cosmic") },
    { "slc_n_primary_daughters", new InhVecVar<int>("slc_n_primary_daughters") },
    { "slc_n_trks", new InhVecVar<int>("slc_n_trks") },
    { "slc_n_shws", new InhVecVar<int>("slc_n_shws") },
    { "slc_n_dazzle_muons", new InhVecVar<int>("slc_n_dazzle_muons") },
    { "slc_n_dazzle_pions", new InhVecVar<int>("slc_n_dazzle_pions") },
    { "slc_n_dazzle_protons", new InhVecVar<int>("slc_n_dazzle_protons") },
    { "slc_n_dazzle_other", new InhVecVar<int>("slc_n_dazzle_other") },
    { "slc_vtx_x", new InhVecVar<double>("slc_vtx_x") },
    { "slc_vtx_y", new InhVecVar<double>("slc_vtx_y") },
    { "slc_vtx_z", new InhVecVar<double>("slc_vtx_z") },
    { "slc_is_fv", new InhVecVar<bool>("slc_is_fv") },
    { "slc_crumbs_score", new InhVecVar<float>("slc_crumbs_score") },
    { "slc_crumbs_nc_score", new InhVecVar<float>("slc_crumbs_nc_score") },
    { "slc_crumbs_ccnue_score", new InhVecVar<float>("slc_crumbs_ccnue_score") },
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
    fTrackModuleLabel            = p.get<art::InputTag>("TrackModuleLabel", "pandoraSCETrack");
    fShowerModuleLabel           = p.get<art::InputTag>("ShowerModuleLabel", "pandoraSCEShower");
    fTrackCalorimetryModuleLabel = p.get<art::InputTag>("TrackCalorimetryModuleLabel", "pandoraSCECalo");
    fCRUMBSModuleLabel           = p.get<art::InputTag>("CRUMBSModuleLabel", "crumbs");
    fDazzleModuleLabel           = p.get<art::InputTag>("DazzleModuleLabel", "dazzle");
    fCaloModuleLabel             = p.get<art::InputTag>("CaloModuleLabel", "pandoraSCECalo");
    fMCSModuleLabel              = p.get<art::InputTag>("MCSModuleLabel", "pandoraTrackMCS:muon");
    fChi2ModuleLabel             = p.get<art::InputTag>("Chi2ModuleLabel", "pandoraSCEPid");
    fRangeModuleLabel            = p.get<art::InputTag>("RangeModuleLabel", "pandoraTrackRange:muon");
    fClosestApproachModuleLabel  = p.get<art::InputTag>("ClosestApproachModuleLabel", "pandoraTrackClosestApproach");
    fStoppingChi2ModuleLabel     = p.get<art::InputTag>("StoppingChi2ModuleLabel", "pandoraTrackStoppingChi2");
    fRazzleModuleLabel           = p.get<art::InputTag>("RazzleModuleLabel", "razzle");
    fCosmicDistModuleLabel       = p.get<art::InputTag>("CosmicDistModuleLabel", "pandoraShowerCosmicDist");
    fShowerTrackFitModuleLabel   = p.get<art::InputTag>("ShowerTrackFitModuleLabel", "pandoraShowerSelectionVars");
    fShowerDensityFitModuleLabel = p.get<art::InputTag>("ShowerDensityFitModuleLabel", "pandoraShowerSelectionVars");
    fDebug                       = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fEventTree = fs->make<TTree>("events","");
    fEventTree->Branch("run", &_run);
    fEventTree->Branch("subrun", &_subrun);
    fEventTree->Branch("event", &_event);

    fEventTree->Branch("n_nu", &_n_nu);
    SetupBranches(nuVars);

    fEventTree->Branch("n_slc", &_n_slc);
    SetupBranches(slcVars);
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

void sbnd::NCPiZeroAnalysis::analyze(const art::Event &e)
{
  ResetVars();
  ClearMaps();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

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

  // Fill the Tree
  fEventTree->Fill();
}

void sbnd::NCPiZeroAnalysis::ClearMaps()
{
  fHitsMap.clear();
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
    fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;

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

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

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

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          const simb::MCNeutrino mcn = mct->GetNeutrino();
          const simb::MCParticle nu  = mcn.Nu();

          const bool nc = mcn.CCNC() == 1;
          const bool av = VolumeCheck(nu.Position().Vect());
          const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

          unsigned pizeros = 0;

          for(int i = 0; i < mct->NParticles(); ++i)
            {
              const auto mcp = mct->GetParticle(i);

              if(mcp.PdgCode() == 111 && mcp.StatusCode() != 1)
                ++pizeros;
            }

          const bool pizero = pizeros > 0;

          if(nc && fv && pizero)
            {
              FillElement(nuVars["nu_event_type"], nuCounter, (int) kNCPiZero);
              FillElement(nuVars["nu_signal"], nuCounter, true);
            }
          else
            {
              FillElement(nuVars["nu_signal"], nuCounter, false);

              if(nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kOtherNC);
              else if(abs(nu.PdgCode()) == 14 && !nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kCCNuMu);
              else if(abs(nu.PdgCode()) == 12 && !nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kCCNuE);
              else if(!fv && av)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kNonFV);
              else if(!av)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kDirt);
              else
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kUnknownEv);
            }

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
          float trueEnDep = 0.;

          for(auto const& mcp : MCParticleVec)
            {
              std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

              for(auto const& ide : ides)
                trueEnDep += ide->energy / 1000.;
            }

          FillElement(nuVars["nu_true_en_dep"], nuCounter, trueEnDep);
          FillElement(nuVars["nu_ccnc"], nuCounter, mcn.CCNC());
          FillElement(nuVars["nu_mode"], nuCounter, mcn.Mode());
          FillElement(nuVars["nu_int_type"], nuCounter, mcn.InteractionType());
          FillElement(nuVars["nu_w"], nuCounter, mcn.W());
          FillElement(nuVars["nu_x"], nuCounter, mcn.X());
          FillElement(nuVars["nu_y"], nuCounter, mcn.Y());
          FillElement(nuVars["nu_q_sqr"], nuCounter, mcn.QSqr());
          FillElement(nuVars["nu_pt"], nuCounter, mcn.Pt());
          FillElement(nuVars["nu_theta"], nuCounter, mcn.Theta());
          FillElement(nuVars["nu_e"], nuCounter, nu.E());
          FillElement(nuVars["nu_vtx_x"], nuCounter, nu.Vx());
          FillElement(nuVars["nu_vtx_y"], nuCounter, nu.Vy());
          FillElement(nuVars["nu_vtx_z"], nuCounter, nu.Vz());

          ++nuCounter;
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

  art::FindManyP<recob::PFParticle>                slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>                     pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  art::FindOneP<sbn::CRUMBSResult>                 sliceToCRUMBS(sliceHandle, e, fCRUMBSModuleLabel);
  art::FindOneP<recob::Track>                      pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpToMeta(pfpHandle, e, fPFParticleModuleLabel);
  art::FindManyP<recob::Hit>                       tracksToHits(trackHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::Hit>                       showersToHits(showerHandle, e, fShowerModuleLabel);

  for (auto&& [slcCounter, slc] : enumerate(sliceVec))
    {
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
        FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, false);

      const art::Ptr<recob::Vertex> vtx = pfpToVertices.at(prim.key());
      geo::Point_t vtxPos = vtx.isNonnull() ? vtx->position() : geo::Point_t(def_double, def_double, def_double);
      FillElement(slcVars["slc_vtx_x"], slcCounter, vtxPos.X());
      FillElement(slcVars["slc_vtx_y"], slcCounter, vtxPos.Y());
      FillElement(slcVars["slc_vtx_z"], slcCounter, vtxPos.Z());
      FillElement(slcVars["slc_is_fv"], slcCounter, VolumeCheck(vtxPos, 20., 5., 10., 50.));

      const art::Ptr<sbn::CRUMBSResult> crumbs = sliceToCRUMBS.at(slc.key());
      if(crumbs.isNonnull())
        {
          FillElement(slcVars["slc_crumbs_score"], slcCounter, crumbs->score);
          FillElement(slcVars["slc_crumbs_nc_score"], slcCounter, crumbs->ncscore);
          FillElement(slcVars["slc_crumbs_ccnue_score"], slcCounter, crumbs->ccnuescore);
        }

      int ntrks = 0, nshws = 0, ndazzlemuons = 0, ndazzlepions = 0, ndazzleprotons = 0, ndazzleother = 0;

      ResizeSubVectors(slcVars, slcCounter, prim->NumDaughters());

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

          const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
          const std::vector<art::Ptr<recob::Hit>> showerHits = showersToHits.at(shower.key());
          const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, true);
          FillElement(slcVars["slc_pfp_true_trackid"], slcCounter, pfpCounter, trackID);
          FillElement(slcVars["slc_pfp_comp"], slcCounter, pfpCounter, Completeness(e, showerHits, trackID));
          FillElement(slcVars["slc_pfp_pur"], slcCounter, pfpCounter, Purity(e, showerHits, trackID));

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
              int dazzlepdg = dynamic_cast<InhVecVecVar<int>*>(slcVars["slc_pfp_track_dazzle_pdg"])->GetVal(slcCounter, pfpCounter);

              if(dazzlepdg == 13)
                ++ndazzlemuons;
              else if(dazzlepdg == 211)
                ++ndazzlepions;
              else if(dazzlepdg == 2212)
                ++ndazzleprotons;
              else if(dazzlepdg == 0)
                ++ndazzleother;
            }

          if(shower.isNonnull())
            AnalyseShower(e, shower, slcCounter, pfpCounter, showerHandle, vtx, showerHits);
        }

      FillElement(slcVars["slc_n_trks"], slcCounter, ntrks);
      FillElement(slcVars["slc_n_shws"], slcCounter, nshws);
      FillElement(slcVars["slc_n_dazzle_muons"], slcCounter, ndazzlemuons);
      FillElement(slcVars["slc_n_dazzle_pions"], slcCounter, ndazzlepions);
      FillElement(slcVars["slc_n_dazzle_protons"], slcCounter, ndazzleprotons);
      FillElement(slcVars["slc_n_dazzle_other"], slcCounter, ndazzleother);
    }
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

  double convGap = (start - vtx->position()).R();
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

void sbnd::NCPiZeroAnalysis::ResetVars()
{
  _run = -1; _subrun = -1; _event  = -1;

  _n_nu = 0; _n_slc  = 0;
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

void sbnd::NCPiZeroAnalysis::ResizeSubVectors(VecVarMap &map, const int pos, const int size)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kTwoD)
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

DEFINE_ART_MODULE(sbnd::NCPiZeroAnalysis)
