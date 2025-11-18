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

// LArSoft MC includes 
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/sim.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// SBND includes
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/Common/Reco/LightCaloInfo.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <numeric>
#include <memory>
#include <algorithm> // sort 
#include <cmath>
#include <functional>
#include <map>

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

  template <typename MatchT, typename CheckFunc>
  void CollectMatches(const art::Handle<std::vector<MatchT>> &handle,
                      const std::vector<art::Ptr<MatchT>> &fm_v,
                      const std::string &label,
                      art::Event &e,
                      std::vector<art::Ptr<recob::Slice>> &match_slices_v,
                      std::vector<art::Ptr<recob::OpFlash>> &match_op0,
                      std::vector<art::Ptr<recob::OpFlash>> &match_op1,
                      CheckFunc check);

  // Returns visibility vector for all opdets given charge/position information 
  std::vector<std::vector<double>> CalcVisibility(std::vector<geo::Point_t> xyz_v, 
                                                  std::vector<double> charge_v);

  // Fills reconstructed photon count vector (total_gamma_v) for all opdets given charge/position information 
  void CalcLight(std::vector<double>  flash_pe_v, 
                 std::vector<double>  dir_visibility,
                 std::vector<double>  ref_visibility, 
                 std::vector<double> &total_gamma_v);

  // Returns the median of the light vector 
  double CalcMedian(std::vector<double> total_light);

  // Returns the mean of the light vector 
  double CalcMean(std::vector<double> total_light);
  double CalcMean(std::vector<double> total_light, std::vector<double> total_err);

  // Performs truth validation, saves info to TTree 
  void TruthValidation(art::Event& e, art::ServiceHandle<cheat::ParticleInventoryService> piserv, double flash_time);

  // fcl parameters 
  std::vector<std::string> _opflash_producer_v;
  std::vector<std::string> _opflash_ara_producer_v;
  std::vector<std::string> _pd_types; 
  std::string _slice_producer;
  std::string _opt0_producer;
  std::string _bcfm_producer;

  bool _use_arapucas;
  bool _use_opt0;
  bool _use_bcfm;
  float _nuscore_cut; 
  float _fopt0score_cut;
  bool _verbose;

  float _fopflash_min;
  float _fopflash_max;
  float _fopt0_frac_diff_cut;

  float _pmt_ara_offset; 
  std::vector<float> _noise_thresh;

  std::vector<float> _cal_area_const; 
  std::vector<float> _opdet_vuv_eff;
  std::vector<float> _opdet_vis_eff;
  std::vector<int> _opdet_mask;
  float _scint_prescale;

  bool _truth_validation; 
  std::string _simenergy_producer;
  bool _truth_neutrino;

  std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;
  std::shared_ptr<phot::OpticalPath> _optical_path_tool;

  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();

  geo::GeometryCore const* geom;

  TTree* _tree;
  int _run, _subrun, _event;
  int _pfpid; // ID of the matched slice 
  double _opflash_time; // time of matched opflash 

  std::vector<double> _dep_pe; // vector of measured photo-electron (PE), one entry = one channel 
  std::vector<double> _rec_gamma; // vector of reconstructed photon count, one entry = one channel 

  double _true_gamma;  // true photon count from all energy depositions 
  double _true_charge; // true electron count from all energy depositions 
  double _true_energy; // true deposited energy 

  std::vector<double> _visibility;

  double _slice_L; // reconstructed photon count 
  double _slice_Q; // reconstructed electron count 
  double _slice_E; // reconstructed deposited energy 

};


sbnd::LightCaloProducer::LightCaloProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _optical_path_tool = std::shared_ptr<phot::OpticalPath>(art::make_tool<phot::OpticalPath>(p.get<fhicl::ParameterSet>("OpticalPathTool")));
  _semi_model = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, _optical_path_tool, true, false);

  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _opflash_ara_producer_v = p.get<std::vector<std::string>>("OpFlashAraProducers");
  _pd_types = p.get<std::vector<std::string>>("PDTypes"); 
  _slice_producer = p.get<std::string>("SliceProducer");
  _opt0_producer  = p.get<std::string>("OpT0FinderProducer");
  _bcfm_producer  = p.get<std::string>("BCFMProducer");
  _use_opt0     = p.get<bool>("UseOpT0Finder");
  _use_bcfm     = p.get<bool>("UseBCFM");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fopt0score_cut = p.get<float>("opt0ScoreCut");
  _verbose = p.get<bool>("Verbose");

  _fopflash_min = p.get<float>("OpFlashMin");
  _fopflash_max = p.get<float>("OpFlashMax");
  _fopt0_frac_diff_cut = p.get<float>("OpT0FractionalCut");

  _pmt_ara_offset  = p.get<float>("PMTARAFlashOffset");
  _noise_thresh    = p.get<std::vector<float>>("FlashNoiseThreshold");

  _cal_area_const  = p.get<std::vector<float>>("CalAreaConstants");
  _opdet_vuv_eff   = p.get<std::vector<float>>("OpDetVUVEfficiencies");
  _opdet_vis_eff   = p.get<std::vector<float>>("OpDetVISEfficiencies");
  _opdet_mask      = p.get<std::vector<int>>("OpDetMask");
  _scint_prescale  = p.get<float>("ScintPreScale");

  _truth_validation   = p.get<bool>("TruthValidation");
  _simenergy_producer = p.get<std::string>("SimEnergyProducer"); 
  _truth_neutrino     = p.get<bool>("TruthNeutrino");

  geom = lar::providerFrom<geo::Geometry>();

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("lightcalo","");
  _tree->Branch("run",           &_run,          "run/I");
  _tree->Branch("subrun",        &_subrun,       "subrun/I");
  _tree->Branch("event",         &_event,        "event/I");
  _tree->Branch("pfpid",         &_pfpid,        "pfpid/I");
  _tree->Branch("opflash_time",  &_opflash_time, "opflash_time/D");

  _tree->Branch("rec_gamma",    "std::vector<double>", &_rec_gamma);
  _tree->Branch("dep_pe",       "std::vector<double>", &_dep_pe);
  _tree->Branch("visibility",   "std::vector<double>", &_visibility);

  _tree->Branch("slice_L",       &_slice_L,      "slice_L/D");
  _tree->Branch("slice_Q",       &_slice_Q,      "slice_Q/D");
  _tree->Branch("slice_E",       &_slice_E,      "slice_E/D");

  _tree->Branch("true_gamma",    &_true_gamma,   "true_gamma/D");
  _tree->Branch("true_charge",   &_true_charge,  "true_charge/D");
  _tree->Branch("true_energy",   &_true_energy,  "true_energy/D");

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
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // get slices 
  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()){
    std::cout << "don't have good slices!" << std::endl;
    return;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_slice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    std::cout << "don't have good PFParticle!" << std::endl;
    return;
  }

  ::art::Handle<std::vector<recob::SpacePoint>> spacepoint_h;
  e.getByLabel(_slice_producer, spacepoint_h);
  if(!spacepoint_h.isValid() || spacepoint_h->empty()) {
    std::cout << "don't have good SpacePoints!" << std::endl;
    return;
  }

  auto const & flash0_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[0]);
  auto const & flash1_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[1]);
  if( (!flash0_h.isValid() || flash0_h->empty()) && (!flash1_h.isValid() || flash1_h->empty())) {
    std::cout << "don't have good PMT flashes from producer " << _opflash_producer_v[0] << " or "  << _opflash_producer_v[1] << std::endl;
    return;
  }

  art::FindManyP<recob::PFParticle> slice_to_pfp (slice_h, e, _slice_producer);
  art::FindManyP<recob::Hit>        slice_to_hit (slice_h, e, _slice_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoint(pfp_h, e, _slice_producer);
  art::FindManyP<recob::Hit>        spacepoint_to_hit(spacepoint_h, e, _slice_producer);

  std::vector<art::Ptr<recob::Slice>> match_slices_v; 
  std::vector<art::Ptr<recob::OpFlash>> match_op0; 
  std::vector<art::Ptr<recob::OpFlash>> match_op1; 

  std::map<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>> match_slice_opflash_map;

  // use templated member helper to fill match vectors
  if (_use_bcfm) {
    ::art::Handle<std::vector<sbn::TPCPMTBarycenterMatch>> bcfm_h;
    e.getByLabel(_bcfm_producer, bcfm_h);
    if(!bcfm_h.isValid() || bcfm_h->empty()) {
      std::cout << "don't have good barycenter matches!" << std::endl;
      return;
    }
    std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> bcfm_v;
    art::fill_ptr_vector(bcfm_v, bcfm_h);

    CollectMatches(bcfm_h, bcfm_v, _bcfm_producer, e, match_slices_v, match_op0, match_op1,
                    [this](art::Ptr<sbn::TPCPMTBarycenterMatch> bcfm) {
                      if (bcfm->flashTime > _fopflash_max || bcfm->flashTime < _fopflash_min) return false;
                      if (bcfm->score < 0.02) return false;
                      return true;
                    });
  }
  else {
    ::art::Handle<std::vector<sbn::OpT0Finder>> opt0_h;
    e.getByLabel(_opt0_producer, opt0_h);
    if(!opt0_h.isValid() || opt0_h->empty()) {
      std::cout << "don't have good OpT0Finder matches!" << std::endl;
      return;
    }
    std::vector<art::Ptr<sbn::OpT0Finder>> opt0_v;
    art::fill_ptr_vector(opt0_v, opt0_h);

    CollectMatches(opt0_h, opt0_v, _opt0_producer, e, match_slices_v, match_op0, match_op1,
                    [this](art::Ptr<sbn::OpT0Finder> opt0){
                      auto opt0_measPE = opt0->measPE;
                      auto opt0_hypoPE = opt0->hypoPE;
                      auto opt0_frac_diff = std::abs((opt0_hypoPE - opt0_measPE)/opt0_measPE);
                      if (opt0->time < _fopflash_min ||  opt0->time > _fopflash_max) return false;
                      if (opt0->score < _fopt0score_cut) return false;
                      if (opt0_frac_diff > _fopt0_frac_diff_cut) return false;
                      return true;
                    });
  }

  if (std::find(_pd_types.begin(), _pd_types.end(), "xarapuca_vuv") != _pd_types.end() || 
      std::find(_pd_types.begin(), _pd_types.end(), "xarapuca_vis") != _pd_types.end()){
        _use_arapucas = true; 
      }
      else
      _use_arapucas = false; 

  std::vector<art::Ptr<recob::OpFlash>> flash0_ara_v;
  std::vector<art::Ptr<recob::OpFlash>> flash1_ara_v;
  if (_use_arapucas){
    ::art::Handle<std::vector<recob::OpFlash>> flash0_ara_h;
    ::art::Handle<std::vector<recob::OpFlash>> flash1_ara_h;

    for (size_t i=0; i<2; i++){
      ::art::Handle<std::vector<recob::OpFlash>> flash_ara_h;
      e.getByLabel(_opflash_ara_producer_v[i], flash_ara_h);
      if (!flash_ara_h.isValid() || flash_ara_h->empty()) {
        std::cout << "don't have good X-ARAPUCA flashes from producer " << _opflash_ara_producer_v[i] << std::endl;
      }
      else art::fill_ptr_vector((i==0)? flash0_ara_v : flash1_ara_v, flash_ara_h);
    }
  }

  int nsuccessful_matches=0;
  for (size_t n_slice=0; n_slice < match_slices_v.size(); n_slice++){
    // initialize tree2 variables 
    _pfpid = -1; 

    _slice_Q = 0; // total amount of charge
    _slice_L = 0; // total amount of light
    _slice_E = 0;

    std::vector<geo::Point_t> sp_xyz;
    std::vector<double> sp_charge; // vector of charge info for charge-weighting
 
    double flash_time = -999;
    auto opflash0 = (match_op0.at(n_slice));
    auto opflash1 = (match_op1.at(n_slice));
    bool flash_in_0 = false;
    bool flash_in_1 = false;
    // set threshold above noise PE levels for the flash 
    float noise_thresh = (!_use_arapucas)? _noise_thresh[0] : _noise_thresh[1]; 

    if (!opflash0.isNull() &&  opflash1.isNull()){
      flash_time = opflash0->Time();
      if (opflash0->TotalPE() > noise_thresh)
        flash_in_0 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 0 (" << opflash0->TotalPE() << ") below noise threshold ... Skipping" << std::endl;
    } 
    else if ( opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash1->Time(); 
      if (opflash1->TotalPE() > noise_thresh)
        flash_in_1 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 1 (" << opflash1->TotalPE() << ") below noise threshold ... Skipping" << std::endl;
    }
    else if (!opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash0->Time();
      if (opflash0->TotalPE() > noise_thresh)
        flash_in_0 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 0 (" << opflash0->TotalPE() << ") below noise threshold ... Skipping" << std::endl;
      if (opflash1->TotalPE() > noise_thresh)
        flash_in_1 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 1 (" << opflash1->TotalPE() << ") below noise threshold ... Skipping" << std::endl;
    }
    else if  (opflash0.isNull() &&  opflash1.isNull()){
      std::cout << "No usable opflashes (none above threshold)." << std::endl;
      return;
    }
    auto slice = match_slices_v[n_slice];
    // sum charge information (without position info) for Q 
    // find which plane has the most integrated charge for this slice
    std::vector<art::Ptr<recob::Hit>> slice_hits_v = slice_to_hit.at(slice.key());
    std::vector<double> plane_charge{0.,0.,0.};
    std::vector<int>    plane_hits{0,0,0};
    for (size_t i=0; i < slice_hits_v.size(); i++){
      auto hit = slice_hits_v[i];
      auto drift_time = hit->PeakTime()*0.5 - clock_data.TriggerOffsetTPC(); 
      double atten_correction = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)
      auto hit_plane = hit->View();
      plane_charge.at(hit_plane) += hit->Integral()*atten_correction*(1/_cal_area_const.at(hit_plane));
      plane_hits.at(hit_plane)++; 
    }

    uint bestPlane = std::max_element(plane_charge.begin(), plane_charge.end()) - plane_charge.begin(); 
    uint bestHits =  std::max_element(plane_hits.begin(), plane_hits.end()) - plane_hits.begin();

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
          // correct for e- attenuation 
          auto drift_time = hit->PeakTime()*0.5 - clock_data.TriggerOffsetTPC(); 
          double atten_correction = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)
          double charge = (1/_cal_area_const.at(bestPlane))*atten_correction*hit->Integral();
          sp_xyz.push_back(xyz);
          sp_charge.push_back(charge);
        }
      } // end spacepoint loop 
    } // end pfp loop

    // get total L count
    std::vector<std::vector<double>> visibility_maps = CalcVisibility(sp_xyz,sp_charge);
    auto dir_visibility_map = visibility_maps[0];
    auto ref_visibility_map = visibility_maps[1];

    std::vector<double> total_pe(_nchan,0.); 
    std::vector<double> total_err(_nchan,0.);
    std::vector<double> total_gamma(_nchan, 0.);

    // combining flash PE information from separate TPCs into a single vector 
    for (int tpc=0; tpc<2; tpc++){
      bool found_flash = (tpc==0)? flash_in_0 : flash_in_1; 
      if (found_flash){
        auto flash_pe_v = (tpc==0)? opflash0->PEs() : opflash1->PEs(); 
        // if using arapucas, need to combine PMT and arapuca PE information into a single vector 
        if (_use_arapucas){
          auto flash_ara_v = (tpc==0)? flash0_ara_v : flash1_ara_v;
          // for PMT flashes, the PE vector is shortened and don't include the last 6 entries for ARAPUCAs 
          if (flash_pe_v.size()!= _nchan) flash_pe_v.insert(flash_pe_v.end(), {0,0,0,0,0,0}); 
          for (size_t nara=0; nara < flash_ara_v.size(); nara++){
            auto const &flash_ara = *flash_ara_v[nara];
            if (abs(flash_time-flash_ara.Time()) < _pmt_ara_offset){
              if (_verbose)
                std::cout << "Combining PMT+XARA Flashes with PMT time: " << flash_time << ", ARA time: " << flash_ara.Time() << std::endl;            
              for (size_t ich=0; ich < (flash_ara.PEs()).size(); ich++)
                flash_pe_v.at(ich) += (flash_ara.PEs()).at(ich);
              break;
            }     
          }
        } // end of arapuca if 
        for (size_t ich=0; ich<flash_pe_v.size(); ich++){
          if (std::find(_pd_types.begin(), _pd_types.end(), _opdetmap.pdType(ich) ) == _pd_types.end() ) continue;
          total_pe[ich] += flash_pe_v[ich];
        }
      }
    } // end of TPC loop 

    // mask out specific channels in opdet mask
    for (size_t imask=0; imask<_opdet_mask.size(); imask++){
      total_pe.at(_opdet_mask.at(imask)) = 0;
    }

    // error is proportional to the amount of light
    // that actually reached the optical detector 
    for (size_t ich=0; ich<total_pe.size(); ich++)
      total_err.at(ich) = std::sqrt(total_pe.at(ich));

    _visibility.resize(_nchan,0);
    // calculate the photon estimates for every entry in total_pe
    CalcLight(total_pe, dir_visibility_map, ref_visibility_map, total_gamma);
    
    // fill tree variables 
    _opflash_time = flash_time;
    _dep_pe = total_pe;
    _rec_gamma = total_gamma;

    // calculate final light estimate     
    _slice_L  = CalcMean(total_gamma,total_err);
    _slice_E = (_slice_L + _slice_Q)*1e-6*g4param->Wph(); // MeV, Wph = 19.5 eV  

    if (_verbose){
      std::cout << "charge: " << _slice_Q << std::endl;
      std::cout << "light:  " << _slice_L << std::endl;
      std::cout << "energy: " << _slice_E << std::endl;
    }

    sbn::LightCalo lightcalo(_slice_Q,_slice_L,_slice_E,bestHits,plane_charge);
    lightcalo_v->push_back(lightcalo);
    util::CreateAssn(*this, e, *lightcalo_v, slice,    *slice_assn_v);
    util::CreateAssn(*this, e, *lightcalo_v, opflash0, *flash_assn_v);
    util::CreateAssn(*this, e, *lightcalo_v, opflash1, *flash_assn_v);

    _true_gamma = 0; 
    _true_charge = 0;
    _true_energy = 0;

    if (_truth_validation){
      ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
      e.getByLabel(_simenergy_producer, energyDeps_h);
      std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
      
      if (!energyDeps_h.isValid() || energyDeps_h->empty())
        std::cout << "Don't have good SimEnergyDeposits!" << std::endl;
      else{ 
        art::fill_ptr_vector(energyDeps, energyDeps_h);
        for (size_t n_dep=0; n_dep < energyDeps.size(); n_dep++){
          auto energyDep = energyDeps[n_dep];
          const auto trackID = energyDep->TrackID();
          const double time = energyDep->Time() * 1e-3; // us 

          art::Ptr<simb::MCTruth> mctruth = piserv->TrackIdToMCTruth_P(trackID);

          if ((_truth_neutrino && mctruth->Origin()==simb::kBeamNeutrino ) || 
              (!_truth_neutrino && abs(time-flash_time) < 1)){
            // note: we divide by the prescale because NumPhotons() stored in simulation has the scint prescale applied 
            _true_gamma  += energyDep->NumPhotons()/_scint_prescale;
            _true_charge += energyDep->NumElectrons(); 
            _true_energy += energyDep->Energy(); 
          }       
        }
      }
    }
    nsuccessful_matches++;
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
      std::cout << "more than two opflash matched to this slice!" << std::endl;
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
    if (found_opflash0 == false && found_opflash1 == false){
      std::cout << "no opflashes matched to this slice" << std::endl;
      continue;
    }
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


std::vector<std::vector<double>> sbnd::LightCaloProducer::CalcVisibility(std::vector<geo::Point_t> xyz_v,
                                                                    std::vector<double> charge_v){
  // returns of two vectors (len is # of opdet) for the visibility for every opdet                                     
  if (xyz_v.size() != charge_v.size()) std::cout << "spacepoint coord and charge vector size mismatch" << std::endl;

  std::vector<double> dir_visibility_map(_nchan, 0);
  std::vector<double> ref_visibility_map(_nchan, 0);
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
    _semi_model->detectedDirectVisibilities(direct_visibility, xyz);
    _semi_model->detectedReflectedVisibilities(reflect_visibility, xyz);
    // if (dir_visibility_map.size() != direct_visibility.size()) std::cout << "mismatch of visibility vector size" << std::endl;

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
    auto vuv_eff = _opdet_vuv_eff.at(ch);
    auto vis_eff = _opdet_vis_eff.at(ch);
    auto tot_visibility = vuv_eff*dir_visibility[ch] + vis_eff*ref_visibility[ch];
    _visibility.at(ch) = tot_visibility;
    if((pe == 0) || std::isinf(1/tot_visibility))
      continue;
    // deposited light is inverse of visibility * PE count 
    total_gamma_v[ch] += (1/tot_visibility)*pe;
  }
}

double sbnd::LightCaloProducer::CalcMedian(std::vector<double> total_gamma){
  std::vector<double> tpc0_gamma; 
  std::vector<double> tpc1_gamma;
  // split into two TPCs 
  for (size_t i=0; i<total_gamma.size(); i++){
    if (total_gamma[i] <=0) continue;
    if (i%2==0) tpc0_gamma.push_back(total_gamma[i]);
    if (i%2==1) tpc1_gamma.push_back(total_gamma[i]);
  }
  double median_gamma=0; 

  for (int tpc=0; tpc < 2; tpc++){
    // make a copy to avoid changing in-place 
    std::vector<double> gamma_v = ((tpc==0)? std::vector<double>(tpc0_gamma) : std::vector<double>(tpc1_gamma));
    const auto median_it = gamma_v.begin() + gamma_v.size() / 2;
    std::nth_element(gamma_v.begin(), median_it , gamma_v.end());
    auto median = *median_it;
    median_gamma+=median;
  }
  return median_gamma;
}

double sbnd::LightCaloProducer::CalcMean(std::vector<double> total_gamma, std::vector<double> total_err){
  // calculates a weighted average, loops over per tpc and then per channel
  double total_mean=0;
  double wgt_num = 0;
  double wgt_denom = 0;

  for (size_t ich=0; ich<total_gamma.size(); ich++){
    if (total_gamma[ich]<=0) continue;
    wgt_num   += total_gamma[ich]*(1./total_err[ich]);
    wgt_denom += (1./total_err[ich]);
  }
  if (wgt_denom!=0) total_mean += wgt_num/wgt_denom;
  return total_mean; 
}

double sbnd::LightCaloProducer::CalcMean(std::vector<double> total_gamma){
  double total_mean=0;
  double wgt_num = 0;
  double wgt_denom = 0;

  for (size_t ich=0; ich<total_gamma.size(); ich++){
    if (total_gamma[ich]<=0) continue;
    wgt_num   += total_gamma[ich];
    wgt_denom += 1.;
  }
  if (wgt_denom!=0) total_mean += wgt_num/wgt_denom;
  return total_mean; 
}

DEFINE_ART_MODULE(sbnd::LightCaloProducer)
