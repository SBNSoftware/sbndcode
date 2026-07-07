////////////////////////////////////////////////////////////////////////
// Class:       LightCaloProducer
// Plugin Type: producer (Unknown Unknown)
// File:        LightCaloProducer_module.cc
//
// Generated at Wed Apr 19 16:49:20 2023 by Lynn Tung using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft includes 
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Calibration database includes for electron lifetime
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "wda.h"

// SBND includes
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/Common/Reco/LightCalo.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

// Calibration database includes
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <memory>
#include <algorithm> // sort 
#include <cmath>
#include <map>
#include <string> 

namespace sbnd {
  class LightCaloProducer;
}


class sbnd::LightCaloProducer : public art::EDProducer {
public:
  explicit LightCaloProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LightCaloProducer(LightCaloProducer const&) = delete;
  LightCaloProducer(LightCaloProducer&&) = delete;
  LightCaloProducer& operator=(LightCaloProducer const&) = delete;
  LightCaloProducer& operator=(LightCaloProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void CalculateCalorimetry(art::Event& e, 
                            std::unique_ptr< std::vector<sbn::LightCalo> >& lightcalo_v, 
                            std::unique_ptr< art::Assns<recob::Slice, sbn::LightCalo> >& slice_assn_v,
                            std::unique_ptr< art::Assns<recob::OpFlash, sbn::LightCalo> >& flash_assn_v);

  // Templated member function to collect matched slices and opflashes
  // The CheckFunc is a callable that takes a Ptr<MatchT> and returns bool
  // to determine if the match passes selection criteria
  template <typename MatchT, typename CheckFunc>
  void CollectMatches(const art::Handle<std::vector<MatchT>> &handle,
                      const std::vector<art::Ptr<MatchT>> &fm_v,
                      const std::string &label,
                      art::Event &e,
                      std::vector<art::Ptr<recob::Slice>> &match_slices_v,
                      std::vector<art::Ptr<recob::OpFlash>> &match_op0,
                      std::vector<art::Ptr<recob::OpFlash>> &match_op1,
                      CheckFunc check);

  // Finds the opflash in flash_v that is closest in time to ref_time                      
  art::Ptr<recob::OpFlash> FindMatchingFlash(const std::vector<art::Ptr<recob::OpFlash>> &flash_v,
                                             double ref_time);                      
  // Returns visibility vector for all opdets given charge/position information 
  std::vector<std::vector<double>> CalcVisibility(std::vector<geo::Point_t> xyz_v, 
                                                  std::vector<double> charge_v);

  // Fills reconstructed photon count vector (total_gamma_v) for all opdets given charge/position information 
  void CalcLight(std::vector<double>  flash_pe_v, 
                 std::vector<double>  dir_visibility,
                 std::vector<double>  ref_visibility, 
                 std::vector<double> &total_gamma_v);

  // Returns the weighted mean of the light vector 
  double CalcMean(std::vector<double> total_light, std::vector<double> total_err);

  // Struct to hold electron lifetime data from database
  struct ELifetimeInfo {
    double tau_tpc0;
    double tau_tpc1;
  };

  // Helper to get electron lifetime from database
  ELifetimeInfo GetELifetimeFromDB(uint64_t run);

  // fcl parameters 
  std::vector<std::string> fopflash_producer_v;
  std::vector<std::string> fopflash_ara_producer_v;
  std::vector<std::string> fpd_types; 
  std::string fslice_producer;
  std::string fopt0_producer;
  std::string fbcfm_producer;

  bool fuse_bcfm;
  bool fuse_opt0;
  bool fverbose;
  bool ffill_tree;

  float fbcfmscore_cut; 
  float fopt0score_cut;

  float fopflash_min;
  float fopflash_max;
  float fopt0_frac_diff_cut;

  float fflash_offset; 
  std::vector<float> fnoise_thresh;
  std::vector<float> fupper_thresh; 

  std::vector<float> fcal_area_const; 

  double fpmtcoated_vuveff_tpc0; 
  double fpmtcoated_viseff_tpc0;
  double fpmtuncoated_eff_tpc0;

  double fpmtcoated_vuveff_tpc1; 
  double fpmtcoated_viseff_tpc1;
  double fpmtuncoated_eff_tpc1;

  double fxarapucavuv_vuveff;
  double fxarapucavuv_viseff;
  double fxarapucavis_eff; 

  bool fuse_elifetime_db;
  std::string felifetime_db_file;
  std::string felifetime_db_tag;
  std::unique_ptr<lariov::DBFolder> felifetime_db;
  std::map<uint32_t, ELifetimeInfo> felifetime_cache;

  std::vector<float> fopdet_vuv_eff;
  std::vector<float> fopdet_vis_eff;
  std::vector<int>   fopdet_mask;

  std::unique_ptr<phot::SemiAnalyticalModel> fsemi_model;
  fhicl::ParameterSet fvuv_params;
  fhicl::ParameterSet fvis_params;
  std::shared_ptr<phot::OpticalPath> foptical_path_tool;

  sbndDB::PMTCalibrationDatabase const* fpmt_calib_db;

  opdet::sbndPDMapAlg opdetmap; //map for photon detector types
  unsigned int nchan = opdetmap.size();

  geo::GeometryCore const* geom;

  TTree* _tree;
  int _run, _subrun, _event;
  int _pfpid; // ID of the matched slice 
  double _opflash_time; // time of matched opflash 

  std::vector<double> _dep_pe; // vector of measured photo-electron (PE), one entry = one channel 
  std::vector<double> _rec_gamma; // vector of reconstructed photon count, one entry = one channel 
  std::vector<double> _visibility;

  double _slice_L; // reconstructed photon count 
  double _slice_Q; // reconstructed electron count 
  double _slice_E; // reconstructed deposited energy 
};


sbnd::LightCaloProducer::LightCaloProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fvuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  fvis_params = p.get<fhicl::ParameterSet>("VIVHits");
  foptical_path_tool = std::shared_ptr<phot::OpticalPath>(art::make_tool<phot::OpticalPath>(p.get<fhicl::ParameterSet>("OpticalPathTool")));
  fsemi_model = std::make_unique<phot::SemiAnalyticalModel>(fvuv_params, fvis_params, foptical_path_tool, true, false);
  fpmt_calib_db = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();

  fopflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  fopflash_ara_producer_v = p.get<std::vector<std::string>>("OpFlashAraProducers");
  fpd_types = p.get<std::vector<std::string>>("PDTypes"); 
  fslice_producer = p.get<std::string>("SliceProducer");
  fopt0_producer  = p.get<std::string>("OpT0FinderProducer");
  fbcfm_producer  = p.get<std::string>("BCFMProducer");
  fuse_opt0     = p.get<bool>("UseOpT0Finder");
  fuse_bcfm     = p.get<bool>("UseBCFM");

  fverbose = p.get<bool>("Verbose");
  ffill_tree = p.get<bool>("FillTree");

  fbcfmscore_cut = p.get<float>("bcfmScoreCut");
  fopt0score_cut = p.get<float>("opt0ScoreCut");
  fopt0_frac_diff_cut = p.get<float>("opt0FractionalCut");

  fopflash_min = p.get<float>("OpFlashMin");
  fopflash_max = p.get<float>("OpFlashMax");

  fflash_offset    = p.get<float>("FlashOffset");
  fnoise_thresh    = p.get<std::vector<float>>("FlashNoiseThreshold");
  fupper_thresh    = p.get<std::vector<float>>("OpDetMaxPEThreshold");

  fcal_area_const  = p.get<std::vector<float>>("CalAreaConstants");
  fopdet_mask      = p.get<std::vector<int>>("OpDetMask");

  fpmtcoated_viseff_tpc0  = p.get<double>("PMTCoatedVISEff_tpc0");
  fpmtcoated_vuveff_tpc0  = p.get<double>("PMTCoatedVUVEff_tpc0");
  fpmtuncoated_eff_tpc0   = p.get<double>("PMTUncoatedEff_tpc0");
  fpmtcoated_viseff_tpc1  = p.get<double>("PMTCoatedVISEff_tpc1");
  fpmtcoated_vuveff_tpc1  = p.get<double>("PMTCoatedVUVEff_tpc1");
  fpmtuncoated_eff_tpc1   = p.get<double>("PMTUncoatedEff_tpc1");

  fxarapucavuv_vuveff = p.get<double>("XArapucaVUVVUVEff");
  fxarapucavuv_viseff = p.get<double>("XArapucaVUVVISEff");
  fxarapucavis_eff    = p.get<double>("XArapucaVISEff");

  // Electron lifetime database configuration
  fuse_elifetime_db = p.get<bool>("UseELifetimeDB", false);
  if (fuse_elifetime_db) {
    felifetime_db_file = p.get<std::string>("ELifetimeDBFile");
    felifetime_db_tag  = p.get<std::string>("ELifetimeDBTag");
    felifetime_db = std::make_unique<lariov::DBFolder>(felifetime_db_file, "", "", felifetime_db_tag, true, false);
  }

  // fill efficiency vectors 
  fopdet_vuv_eff.resize(nchan,0.);
  fopdet_vis_eff.resize(nchan,0.);
  for (unsigned int ch = 0; ch < nchan; ++ch) {
    std::string pd_type = opdetmap.pdType(ch);
    if (pd_type == "pmt_coated") {
      if (ch%2==0){
        fopdet_vuv_eff[ch] = fpmtcoated_vuveff_tpc0;
        fopdet_vis_eff[ch] = fpmtcoated_viseff_tpc0;
      }
      else{
        fopdet_vuv_eff[ch] = fpmtcoated_vuveff_tpc1;
        fopdet_vis_eff[ch] = fpmtcoated_viseff_tpc1;
      }
    } else if (pd_type == "pmt_uncoated") {
      if (ch%2==0)
        fopdet_vis_eff[ch] = fpmtuncoated_eff_tpc0;
      else
        fopdet_vis_eff[ch] = fpmtuncoated_eff_tpc1;
    }
    else if (pd_type == "xarapuca_vuv") {
      fopdet_vuv_eff[ch] = fxarapucavuv_vuveff;
      fopdet_vis_eff[ch] = fxarapucavuv_viseff;
    }
    else if (pd_type == "xarapuca_vis") {
      fopdet_vis_eff[ch] = fxarapucavis_eff;
    }
  }

  geom = lar::providerFrom<geo::Geometry>();

  art::ServiceHandle<art::TFileService> fs;
  if (ffill_tree){
    _tree = fs->make<TTree>("lightcalo","");
    _tree->Branch("run",           &_run,          "run/I");
    _tree->Branch("subrun",        &_subrun,       "subrun/I");
    _tree->Branch("event",         &_event,        "event/I");
    _tree->Branch("pfpid",         &_pfpid,        "pfpid/I");
    _tree->Branch("opflash_time",  &_opflash_time, "opflash_time/D");

    _tree->Branch("rec_gamma",    "std::vector<double>", &_rec_gamma);
    _tree->Branch("dep_pe",       "std::vector<double>", &_dep_pe);
    _tree->Branch("visibility",   "std::vector<double>", &_visibility);

    _tree->Branch("slice_Q",       &_slice_Q,      "slice_Q/D");
    _tree->Branch("slice_L",       &_slice_L,      "slice_L/D");
    _tree->Branch("slice_E",       &_slice_E,      "slice_E/D");
  }
  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::LightCalo>>();
  produces<art::Assns<recob::Slice,   sbn::LightCalo>>();
  produces<art::Assns<recob::OpFlash, sbn::LightCalo>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::LightCaloProducer::produce(art::Event& e)
{

  std::unique_ptr<std::vector<sbn::LightCalo>> lightcalo_v (new std::vector<sbn::LightCalo>);
  std::unique_ptr< art::Assns<recob::Slice, sbn::LightCalo>> slice_assn_v (new art::Assns<recob::Slice, sbn::LightCalo>);
  std::unique_ptr< art::Assns<recob::OpFlash, sbn::LightCalo>> flash_assn_v (new art::Assns<recob::OpFlash, sbn::LightCalo>);

  CalculateCalorimetry(e, lightcalo_v, slice_assn_v, flash_assn_v);

  e.put(std::move(lightcalo_v));
  e.put(std::move(slice_assn_v));
  e.put(std::move(flash_assn_v));
}

void sbnd::LightCaloProducer::CalculateCalorimetry(art::Event& e, 
                                                   std::unique_ptr<std::vector<sbn::LightCalo>>& lightcalo_v,
                                                   std::unique_ptr< art::Assns<recob::Slice, sbn::LightCalo>>& slice_assn_v,
                                                   std::unique_ptr< art::Assns<recob::OpFlash, sbn::LightCalo>>& flash_assn_v)
{
  // services 
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
  
  art::ServiceHandle<sim::LArG4Parameters const> g4param;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // get slices 
  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(fslice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()){
    std::cout << "[LightCaloProducer] : " << fslice_producer << " doesn't have good slices!" << std::endl;
    return;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(fslice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    std::cout << "[LightCaloProducer] : " << fslice_producer << " doesn't have good PFParticles!" << std::endl;
    return;
  }

  ::art::Handle<std::vector<recob::SpacePoint>> spacepoint_h;
  e.getByLabel(fslice_producer, spacepoint_h);
  if(!spacepoint_h.isValid() || spacepoint_h->empty()) {
    std::cout << "[LightCaloProducer] : " << fslice_producer << " don't have good SpacePoints!" << std::endl;
    return;
  }
  
  std::vector<art::Ptr<recob::OpFlash>> flash0_v;
  std::vector<art::Ptr<recob::OpFlash>> flash1_v;

  for (size_t i=0; i<fopflash_producer_v.size(); i++){
    ::art::Handle<std::vector<recob::OpFlash>> flash_h;
    e.getByLabel(fopflash_producer_v[i], flash_h);
    if (!flash_h.isValid() || flash_h->empty()) {
      std::cout << "[LightCaloProducer] : " << "don't have good PMT flashes from producer " << fopflash_producer_v[i] << std::endl;
    }
    else{
      if (fopflash_producer_v[i].find("tpc0") != std::string::npos)
        art::fill_ptr_vector(flash0_v, flash_h);
      else if (fopflash_producer_v[i].find("tpc1") != std::string::npos)
        art::fill_ptr_vector(flash1_v, flash_h);
    }
  }

  art::FindManyP<recob::PFParticle> slice_to_pfp (slice_h, e, fslice_producer);
  art::FindManyP<recob::Hit>        slice_to_hit (slice_h, e, fslice_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoint(pfp_h, e, fslice_producer);
  art::FindManyP<recob::Hit>        spacepoint_to_hit(spacepoint_h, e, fslice_producer);

  std::vector<art::Ptr<recob::Slice>> match_slices_v; 
  std::vector<art::Ptr<recob::OpFlash>> match_op0; 
  std::vector<art::Ptr<recob::OpFlash>> match_op1; 

  std::map<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>> match_slice_opflash_map;

  // use templated member helper to fill match vectors
  if (fuse_bcfm) {
    ::art::Handle<std::vector<sbn::TPCPMTBarycenterMatch>> bcfm_h;
    e.getByLabel(fbcfm_producer, bcfm_h);
    if(!bcfm_h.isValid() || bcfm_h->empty()) {
      std::cout << "[LightCaloProducer] : " << "don't have good barycenter matches!" << std::endl;
      return;
    }
    std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> bcfm_v;
    art::fill_ptr_vector(bcfm_v, bcfm_h);

    CollectMatches(bcfm_h, bcfm_v, fbcfm_producer, e, match_slices_v, match_op0, match_op1,
                    [this](art::Ptr<sbn::TPCPMTBarycenterMatch> bcfm) {
                      if (bcfm->flashTime > fopflash_max || bcfm->flashTime < fopflash_min) return false;
                      if (bcfm->score < fbcfmscore_cut) return false;
                      return true;
                    });
  }
  else if (fuse_opt0){
    ::art::Handle<std::vector<sbn::OpT0Finder>> opt0_h;
    e.getByLabel(fopt0_producer, opt0_h);
    if(!opt0_h.isValid() || opt0_h->empty()) {
      std::cout << "[LightCaloProducer] : " << "don't have good OpT0Finder matches!" << std::endl;
      return;
    }
    std::vector<art::Ptr<sbn::OpT0Finder>> opt0_v;
    art::fill_ptr_vector(opt0_v, opt0_h);

    CollectMatches(opt0_h, opt0_v, fopt0_producer, e, match_slices_v, match_op0, match_op1,
                    [this](art::Ptr<sbn::OpT0Finder> opt0){
                      auto opt0_measPE = opt0->measPE;
                      auto opt0_hypoPE = opt0->hypoPE;
                      auto opt0_frac_diff = std::abs((opt0_hypoPE - opt0_measPE)/opt0_measPE);
                      if (opt0->time < fopflash_min ||  opt0->time > fopflash_max) return false;
                      if (opt0->score < fopt0score_cut) return false;
                      if (opt0_frac_diff > fopt0_frac_diff_cut) return false;
                      return true;
                    });
  }

  if (match_slices_v.empty()) return; 

  bool use_arapucas;
  if (std::find(fpd_types.begin(), fpd_types.end(), "xarapuca_vuv") != fpd_types.end() || 
      std::find(fpd_types.begin(), fpd_types.end(), "xarapuca_vis") != fpd_types.end()){
        use_arapucas = true; 
      }
      else use_arapucas = false; 

  std::vector<art::Ptr<recob::OpFlash>> flash0_ara_v;
  std::vector<art::Ptr<recob::OpFlash>> flash1_ara_v;
  if (use_arapucas){
    ::art::Handle<std::vector<recob::OpFlash>> flash0_ara_h;
    ::art::Handle<std::vector<recob::OpFlash>> flash1_ara_h;

    for (size_t i=0; i<fopflash_ara_producer_v.size(); i++){
      ::art::Handle<std::vector<recob::OpFlash>> flash_ara_h;
      e.getByLabel(fopflash_ara_producer_v[i], flash_ara_h);
      if (!flash_ara_h.isValid() || flash_ara_h->empty()) {
        std::cout << "[LightCaloProducer] : " << "don't have good X-ARAPUCA flashes from producer " << fopflash_ara_producer_v[i] << std::endl;
      }
      else{
        if (fopflash_ara_producer_v[i].find("tpc0") != std::string::npos)
          art::fill_ptr_vector(flash0_ara_v, flash_ara_h);
        else if (fopflash_ara_producer_v[i].find("tpc1") != std::string::npos)
          art::fill_ptr_vector(flash1_ara_v, flash_ara_h);
      }
    }
  }

  // Get electron lifetime (from database if enabled, otherwise from detector properties)
  double elifetime_tpc0, elifetime_tpc1;
  if (fuse_elifetime_db) {
    ELifetimeInfo elife_info = GetELifetimeFromDB(e.id().run());
    elifetime_tpc0 = elife_info.tau_tpc0;
    elifetime_tpc1 = elife_info.tau_tpc1;
  } else {
    elifetime_tpc0 = det_prop.ElectronLifetime();
    elifetime_tpc1 = det_prop.ElectronLifetime();
  }

  for (size_t n_slice=0; n_slice < match_slices_v.size(); n_slice++){
    // initialize tree variables 
    _pfpid = -1;
    _opflash_time=-1e9;
    _rec_gamma.clear(); 
    _dep_pe.clear();
    _visibility.clear();

    _slice_Q = 0; 
    _slice_L = 0; 
    _slice_E = 0;

    std::vector<geo::Point_t> sp_xyz;
    std::vector<double> sp_charge; // vector of charge info for charge-weighting

    auto slice = match_slices_v[n_slice];
    // sum charge information (without position info) for Q 
    // find which plane has the most integrated charge for this slice
    std::vector<art::Ptr<recob::Hit>> slice_hits_v = slice_to_hit.at(slice.key());
    std::vector<double> plane_charge{0.,0.,0.};
    std::vector<int>    plane_hits{0,0,0};
    bool has_sp0 = false;
    bool has_sp1 = false;

    for (size_t i=0; i < slice_hits_v.size(); i++){
      auto hit = slice_hits_v[i];
      auto drift_time = clock_data.TPCTick2TrigTime(hit->PeakTime());
      double elifetime = (hit->WireID().TPC == 0) ? elifetime_tpc0 : elifetime_tpc1;
      double atten_correction = std::exp(drift_time/elifetime); // exp(us/us)
      auto hit_plane = hit->View();
      plane_charge.at(hit_plane) += hit->Integral()*atten_correction*(1/fcal_area_const.at(hit_plane));
      plane_hits.at(hit_plane)++; 
    }

    int bestHits =  std::max_element(plane_hits.begin(), plane_hits.end()) - plane_hits.begin();

    _slice_Q = plane_charge.at(bestHits);

    // get charge information to create the weighted map 
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfp.at(slice.key());
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];
      if (pfp->IsPrimary()) _pfpid = pfp->Self();

      std::vector<art::Ptr<recob::SpacePoint>> sp_v = pfp_to_spacepoint.at(pfp.key());
      for (size_t n_sp=0; n_sp < sp_v.size(); n_sp++){
        auto sp = sp_v[n_sp];
        std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hit.at(sp.key());
        for (size_t n_hit=0; n_hit < hit_v.size(); n_hit++){
          auto hit = hit_v[n_hit];
          //
          if (hit->View() !=bestHits) continue;
          const auto &position(sp->XYZ());
          geo::Point_t xyz(position[0],position[1],position[2]);
          auto drift_time = clock_data.TPCTick2TrigTime(hit->PeakTime());
          double elifetime = (hit->WireID().TPC == 0) ? elifetime_tpc0 : elifetime_tpc1;
          double atten_correction = std::exp(drift_time/elifetime); // exp(us/us)
          double charge = (1/fcal_area_const.at(hit->View()))*atten_correction*hit->Integral();
          if (xyz.X() < 0) has_sp0 = true;
          else             has_sp1 = true;

          sp_xyz.push_back(xyz);
          sp_charge.push_back(charge);
        }
      } // end spacepoint loop 
    } // end pfp loop

    double flash_time = -999;
    auto opflash0 = (match_op0.at(n_slice));
    auto opflash1 = (match_op1.at(n_slice));
    bool flash_in_0 = false;
    bool flash_in_1 = false;
    // find matching flashes if slice has reco spacepoints in both TPCs
    if ((has_sp0 && has_sp1) && (opflash0.isNull() || opflash1.isNull())){
      if (opflash0.isNull() && opflash1.isNull()){
        if (fverbose) std::cout << "[LightCaloProducer] : " << "No opflashes found for slice with spacepoints in both TPCs." << std::endl;
        return;
      }
      else if (opflash0.isNull()) opflash0 = FindMatchingFlash(flash0_v, opflash1->Time());
      else if (opflash1.isNull()) opflash1 = FindMatchingFlash(flash1_v, opflash0->Time());
    }

    if (!opflash0.isNull() && opflash0->TotalPE() > fnoise_thresh[0]){
      flash_in_0 = true;
      flash_time = opflash0->Time();
    }
    if (!opflash1.isNull() && opflash1->TotalPE() > fnoise_thresh[0]){
      flash_in_1 = true; 
      flash_time = opflash1->Time(); 
    }

    if (flash_in_0==false && flash_in_1==false)
      return;

    // get total L count
    std::vector<std::vector<double>> visibility_maps = CalcVisibility(sp_xyz,sp_charge);
    auto dir_visibility_map = visibility_maps[0];
    auto ref_visibility_map = visibility_maps[1];

    std::vector<double> total_pe(nchan,0.); 
    std::vector<double> total_err(nchan,0.);
    std::vector<double> total_gamma(nchan, 0.);

    // combining flash PE information from separate TPCs into a single vector 
    for (int tpc=0; tpc<2; tpc++){
      bool found_flash = (tpc==0)? flash_in_0 : flash_in_1; 
      if (found_flash){
        auto flash_pe_v = (tpc==0)? opflash0->PEs() : opflash1->PEs(); 
        // if using arapucas, need to combine PMT and arapuca PE information into a single vector 
        if (use_arapucas){
          auto flash_ara_v = (tpc==0)? flash0_ara_v : flash1_ara_v;
          // for PMT flashes, the PE vector is shortened and don't include the last 6 entries for ARAPUCAs 
          if (flash_pe_v.size()!= nchan) flash_pe_v.resize(nchan,0);
          // find matching xarapuca flash
          auto flash_ara = FindMatchingFlash(flash_ara_v, flash_time);
          if (flash_ara.isNull()) continue;
          // add arapuca PEs to the main vector
          for (size_t ich=0; ich < (flash_ara->PEs()).size(); ich++){
            if ( (flash_ara->PEs()).at(ich) > fupper_thresh[1]) continue; // skip if above max PE threshold
            flash_pe_v.at(ich) += (flash_ara->PEs()).at(ich);
          }
        } // end of arapuca if 
        for (size_t ich=0; ich<flash_pe_v.size(); ich++){
          if (std::find(fpd_types.begin(), fpd_types.end(), opdetmap.pdType(ich) ) == fpd_types.end() ) continue;
          if (flash_pe_v[ich] > fupper_thresh[0]) continue; // skip if above max PE threshold
          total_pe[ich] += flash_pe_v[ich];
        }
      }
    } // end of TPC loop 

    // error is proportional to the amount of light
    // that actually reached the optical detector 
    for (size_t ich=0; ich<total_pe.size(); ich++){
      // if pmt channel is not on or reconstructable, set to zero
      if (opdetmap.pdType(ich) == "pmt_coated" || opdetmap.pdType(ich) == "pmt_uncoated"){
        if(!fpmt_calib_db->getReconstructChannel(ich) || !fpmt_calib_db->getOnPMT(ich))
          total_pe.at(ich) = 0;
      }
      // # if ich is in fopdet_mask, set to zero 
      if (fopdet_mask.end() != std::find(fopdet_mask.begin(), fopdet_mask.end(), ich))
        total_pe.at(ich) = 0;
      total_err.at(ich) = std::sqrt(total_pe.at(ich));
    }
    _visibility.resize(nchan,0);
    // calculate the photon estimates for every entry in total_pe
    CalcLight(total_pe, dir_visibility_map, ref_visibility_map, total_gamma);
    
    // fill tree variables 
    _opflash_time = flash_time;
    _dep_pe = total_pe;
    _rec_gamma = total_gamma;

    // calculate final light estimate     
    _slice_L  = CalcMean(total_gamma,total_err);
    _slice_E = (_slice_L + _slice_Q)*1e-9*g4param->Wph(); // GeV, Wph = 19.5 eV  

    if (fverbose){
      std::cout << "[LightCaloProducer] : " << "charge: " << _slice_Q << std::endl;
      std::cout << "[LightCaloProducer] : " << "light:  " << _slice_L << std::endl;
      std::cout << "[LightCaloProducer] : " << "energy: " << _slice_E << std::endl;
    }

    sbn::LightCalo lightcalo{_slice_Q,_slice_L,_slice_E,bestHits};
    lightcalo_v->push_back(lightcalo);
    util::CreateAssn(*this, e, *lightcalo_v, slice,    *slice_assn_v);
    util::CreateAssn(*this, e, *lightcalo_v, opflash0, *flash_assn_v);
    util::CreateAssn(*this, e, *lightcalo_v, opflash1, *flash_assn_v);

    if (ffill_tree)
      _tree->Fill();
  } // end slice loop
} // end produce 


// define functions 
template <typename MatchT, typename CheckFunc>
void sbnd::LightCaloProducer::CollectMatches(const art::Handle<std::vector<MatchT>> &handle,
                                             const std::vector<art::Ptr<MatchT>> &fm_v,
                                             const std::string &label,
                                             art::Event &e,
                                             std::vector<art::Ptr<recob::Slice>> &match_slices_v,
                                             std::vector<art::Ptr<recob::OpFlash>> &match_op0,
                                             std::vector<art::Ptr<recob::OpFlash>> &match_op1,
                                             CheckFunc check)
{
  std::map<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>> match_slice_opflash_map;
  art::FindManyP<recob::Slice> fm_to_slice(handle, e, label);
  art::FindManyP<recob::OpFlash> fm_to_flash(handle, e, label);

  for (size_t n_fm = 0; n_fm < fm_v.size(); ++n_fm) {
    auto fm = fm_v[n_fm];
    if (!check(fm)) continue;

    std::vector<art::Ptr<recob::Slice>> slice_v = fm_to_slice.at(fm.key());
    std::vector<art::Ptr<recob::OpFlash>> flash_v = fm_to_flash.at(fm.key());

    assert(slice_v.size() == 1);
    assert(flash_v.size() == 1);

    auto slice = slice_v.front();
    auto flash = flash_v.front();

    auto it = match_slice_opflash_map.find(slice);
    if (it == match_slice_opflash_map.end()){
      std::vector<art::Ptr<recob::OpFlash>> fv;
      fv.push_back(flash);
      match_slice_opflash_map.insert(std::pair<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>>(slice, fv));
    }
    else{
      it->second.push_back(flash);
    }
  }

  for (auto it = match_slice_opflash_map.begin(); it != match_slice_opflash_map.end(); ++it){
    auto slice = it->first;
    auto flash_v = it->second;
    if (flash_v.size() > 2){
      std::cout << "[LightCaloProducer] : " << "more than two opflash matched to this slice!" << std::endl;
      continue;
    }
    bool found_opflash0 = false;
    bool found_opflash1 = false;

    for (size_t n_flash=0; n_flash < flash_v.size(); n_flash++){
      auto flash = flash_v[n_flash];
      if (flash->XCenter() > 0){
        found_opflash1 = true;
        match_op1.push_back(flash);
      }
      else if (flash->XCenter() < 0){
        found_opflash0 = true;
        match_op0.push_back(flash);
      }
    } // end opflash loop
    if (found_opflash0 == false && found_opflash1 == false)
      continue;
    else if (found_opflash0 || found_opflash1){
      match_slices_v.push_back(slice);
      art::Ptr<recob::OpFlash> nullOpFlash;
      if (found_opflash0==false) {
        match_op0.push_back(nullOpFlash);
      }
      else if (found_opflash1==false){
        match_op1.push_back(nullOpFlash);
      }
    }
  }
}

art::Ptr<recob::OpFlash> sbnd::LightCaloProducer::FindMatchingFlash(const std::vector<art::Ptr<recob::OpFlash>> &flash_v,
                                                                    double ref_time)
{
  double best_dt = std::numeric_limits<double>::max();  
  art::Ptr<recob::OpFlash> best;
  for (size_t i = 0; i < flash_v.size(); ++i) {
    art::Ptr<recob::OpFlash> cand = flash_v[i];
    double dt = std::abs(cand->Time() - ref_time);
    if (dt < best_dt && dt < fflash_offset) { best_dt = dt; best = cand; }
  }
  return best;
}


std::vector<std::vector<double>> sbnd::LightCaloProducer::CalcVisibility(std::vector<geo::Point_t> xyz_v,
                                                                         std::vector<double> charge_v){
  // returns of two vectors (len is # of opdet) for the visibility for every opdet                                     
  std::vector<double> dir_visibility_map(nchan, 0);
  std::vector<double> ref_visibility_map(nchan, 0);
  double sum_charge0 = 0;
  double sum_charge1 = 0;

  for (size_t i=0; i<xyz_v.size(); i++){
    geo::Point_t const xyz = xyz_v[i];
    auto charge = charge_v[i];
    if (charge <=0) continue;

    if (xyz.X() < 0) sum_charge0+= charge;
    else sum_charge1 += charge; 

    std::vector<double> direct_visibility;
    std::vector<double> reflect_visibility;
    fsemi_model->detectedDirectVisibilities(direct_visibility, xyz);
    fsemi_model->detectedReflectedVisibilities(reflect_visibility, xyz);

    // weight by charge
    for (size_t ch=0; ch<direct_visibility.size(); ch++){
      dir_visibility_map[ch] += charge*direct_visibility[ch];
      ref_visibility_map[ch] += charge*reflect_visibility[ch];
    }    
  } // end spacepoint loop
  // normalize by the total charge in each TPC
  for (size_t ch=0; ch < dir_visibility_map.size(); ch++){
    if (ch%2 == 0){
      dir_visibility_map[ch] /= sum_charge0;
      ref_visibility_map[ch] /= sum_charge0;
    }
    if (ch%2 == 1){
      dir_visibility_map[ch] /= sum_charge1;
      ref_visibility_map[ch] /= sum_charge1;
    }
  }
  return {dir_visibility_map, ref_visibility_map};
}
void sbnd::LightCaloProducer::CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> dir_visibility,
                                   std::vector<double> ref_visibility,
                                   std::vector<double> &total_gamma_v){
  for (size_t ch = 0; ch < flash_pe_v.size(); ch++){
    auto pe = flash_pe_v[ch];
    auto vuv_eff = fopdet_vuv_eff.at(ch);
    auto vis_eff = fopdet_vis_eff.at(ch);
    auto tot_visibility = vuv_eff*dir_visibility[ch] + vis_eff*ref_visibility[ch];
    _visibility.at(ch) = tot_visibility;
    if((pe == 0) || std::isinf(1/tot_visibility))
      continue;
    // deposited light is inverse of visibility * PE count 
    total_gamma_v[ch] += (1/tot_visibility)*pe;
  }
}

double sbnd::LightCaloProducer::CalcMean(std::vector<double> total_gamma, std::vector<double> total_err){
  // calculates a weighted average, loops over per tpc and then per channel
  double total_mean=0;
  for (int tpc=0; tpc<2; tpc++){
    double wgt_num = 0;
    double wgt_denom = 0;

    for (size_t ich=0; ich<total_gamma.size(); ich++){
      if (int(ich%2) != tpc) continue;
      if (total_gamma[ich]<=0) continue;
      wgt_num   += total_gamma[ich]*(1./total_err[ich]);
      wgt_denom += (1./total_err[ich]);
    }
    if (wgt_denom!=0) total_mean += wgt_num/wgt_denom;
  }
  return total_mean; 
}

sbnd::LightCaloProducer::ELifetimeInfo sbnd::LightCaloProducer::GetELifetimeFromDB(uint64_t run) {
  // Check cache first
  if (felifetime_cache.count(run)) {
    return felifetime_cache.at(run);
  }

  // Query database - translate run into fake "timestamp" 
  // (same convention as NormalizeDriftSQLite)
  felifetime_db->UpdateData((run + 1000000000) * 1000000000);

  ELifetimeInfo info;
  double tau_E, tau_W;
  felifetime_db->GetNamedChannelData(0, "etau_sce_spatial_east", tau_E);
  felifetime_db->GetNamedChannelData(0, "etau_sce_spatial_west", tau_W);
  
  info.tau_tpc0 = tau_E*1e3; // the db value is in ms, convert to us
  info.tau_tpc1 = tau_W*1e3; // the db value is in ms, convert to us

  if (fverbose) {
    std::cout << "[LightCaloProducer] : Electron lifetime from DB for run " << run << std::endl;
    std::cout << "[LightCaloProducer] : TPC0 (East): " << info.tau_tpc0 << " us" << std::endl;
    std::cout << "[LightCaloProducer] : TPC1 (West): " << info.tau_tpc1 << " us" << std::endl;
  }

  // Cache the result
  felifetime_cache[run] = info;
  return info;
}

DEFINE_ART_MODULE(sbnd::LightCaloProducer)
