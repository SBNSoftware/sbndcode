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
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/sim.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// SBND includes
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/Common/Reco/LightCaloInfo.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <numeric>
#include <memory>
#include <algorithm> // sort 

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

  // Performs truth validation, saves info to TTree 
  void TruthValidation(art::Event& e, art::ServiceHandle<cheat::ParticleInventoryService> piserv, double flash_time);

  // fcl parameters 
  std::vector<std::string> _opflash_producer_v;
  std::vector<std::string> _opflash_ara_producer_v;
  std::string _slice_producer;
  std::string _opt0_producer;
  bool _use_arapucas;
  bool _use_opt0;
  float _nuscore_cut; 
  float _fopt0score_cut;
  bool _verbose;

  float _fopt0_flash_min;
  float _fopt0_flash_max;
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

  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();

  geo::GeometryCore const* geom;

  TTree* _tree;
  int _run, _subrun, _event;
  int _match_type;
  // match_type key: 
  /// -1: no slices passed both the nuscore and flash match score cut (in the entire event)
  /// -2: no opflashes times found in coincidence with simpleflash time (in the entire event)
  /// -3: no opflashes in coincidence with simpleflash (for this slice)
  /// -4: opflashes are below the noise threshold (for this slice)
  /// 1: successful match 

  TTree* _tree2;
  int _nmatch=0; // number of matches in an event 
  int _pfpid; // ID of the matched slice 
  double _opflash_time; // time of matched opflash 

  std::vector<double> _dep_pe; // vector of measured photo-electron (PE), one entry = one channel 
  std::vector<double> _rec_gamma; // vector of reconstructed photon count, one entry = one channel 

  double _true_gamma;  // true photon count from all energy depositions 
  double _true_charge; // true electron count from all energy depositions 
  double _true_energy; // true deposited energy 

  std::vector<double> _visibility;

  double _median_gamma; // median of all reconstructed light estimates 
  double _mean_gamma;   // mean of all reconstructed light estimates

  double _mean_charge; // avg charge from all three planes 
  double _max_charge;  // charge from the plane with the highest amount of charge
  double _comp_charge; // charge from the plane with the highest number of hits, "highest completeness"

  double _slice_L; // reconstructed photon count 
  double _slice_Q; // reconstructed electron count 
  double _slice_E; // reconstructed deposited energy 

  std::vector<double> _charge = std::vector<double>(3); // reconstructed electron count per plane 
  std::vector<double> _light_med = std::vector<double>(3);  // median reconstructed photon count per plane
  std::vector<double> _light_avg = std::vector<double>(3);  // average reconstructed photon count per plane
  std::vector<double> _energy = std::vector<double>(3);
};


sbnd::LightCaloProducer::LightCaloProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);

  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _opflash_ara_producer_v = p.get<std::vector<std::string>>("OpFlashAraProducers");
  _slice_producer = p.get<std::string>("SliceProducer");
  _opt0_producer  = p.get<std::string>("OpT0FinderProducer");
  _use_arapucas = p.get<bool>("UseArapucas");
  _use_opt0     = p.get<bool>("UseOpT0Finder");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fopt0score_cut = p.get<float>("opt0ScoreCut");
  _verbose = p.get<bool>("Verbose");

  _fopt0_flash_min = p.get<float>("OpT0FlashMin");
  _fopt0_flash_max = p.get<float>("OpT0FlashMax");
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

  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::LightCalo>>();
  produces<art::Assns<recob::Slice, sbn::LightCalo>>();
  // Call appropriate consumes<>() for any products to be retrieved by this module.

}

void sbnd::LightCaloProducer::produce(art::Event& e)
{
  // services 
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
  
  art::ServiceHandle<sim::LArG4Parameters const> g4param;
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  std::cout << "run: " << _run <<  ", subrun: " << _subrun  << ", event: " <<  _event << std::endl;

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
  if (!_use_arapucas && _verbose)
    std::cout << "Using PMT OpFlash only..." << std::endl;  
  if( (!flash0_h.isValid() || flash0_h->empty()) && (!flash1_h.isValid() || flash1_h->empty())) {
    std::cout << "don't have good PMT flashes from producer " << _opflash_producer_v[0] << " or "  << _opflash_producer_v[1] << std::endl;
    return;
  }

  art::FindManyP<recob::PFParticle> slice_to_pfp (slice_h, e, _slice_producer);
  art::FindManyP<recob::Hit>        slice_to_hit (slice_h, e, _slice_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoint(pfp_h, e, _slice_producer);
  art::FindManyP<recob::Hit> spacepoint_to_hit(spacepoint_h, e, _slice_producer);

  std::vector<art::Ptr<recob::Slice>> match_slices_v; 
  std::vector<art::Ptr<recob::OpFlash>> match_op0; 
  std::vector<art::Ptr<recob::OpFlash>> match_op1; 

  if (_use_opt0){
    ::art::Handle<std::vector<sbn::OpT0Finder>> opt0_h;
    e.getByLabel(_opt0_producer, opt0_h);
    if(!opt0_h.isValid() || opt0_h->empty()) {
      std::cout << "don't have good OpT0Finder matches!" << std::endl;
      return;
    }
    std::vector<art::Ptr<sbn::OpT0Finder>> opt0_v;
    art::fill_ptr_vector(opt0_v, opt0_h);

    std::map<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>> match_slice_opflash_map;
    art::FindManyP<recob::Slice> opt0_to_slice(opt0_h, e, _opt0_producer);
    art::FindManyP<recob::OpFlash> opt0_to_flash(opt0_h, e, _opt0_producer);

    for (size_t n_opt0=0; n_opt0 < opt0_v.size(); n_opt0++){
      auto opt0 = opt0_v[n_opt0];
      std::vector<art::Ptr<recob::Slice>> slice_v = opt0_to_slice.at(opt0.key());
      std::vector<art::Ptr<recob::OpFlash>> flash_v = opt0_to_flash.at(opt0.key());
      
      assert(slice_v.size() == 1);
      assert(flash_v.size() == 1);

      auto slice = slice_v.front();
      auto flash = flash_v.front();

      auto opt0_score = opt0->score;
      auto opt0_time = opt0->time; 
      auto opt0_measPE = opt0->measPE;
      auto opt0_hypoPE = opt0->hypoPE;
      auto opt0_frac_diff = std::abs((opt0_hypoPE - opt0_measPE)/opt0_measPE);

      if (opt0_time < _fopt0_flash_min || opt0_time > _fopt0_flash_max) continue;
      if (opt0_score < _fopt0score_cut) continue;
      if (opt0_frac_diff > _fopt0_frac_diff_cut) continue;

      // check that slice is not already in the map 
      auto it = match_slice_opflash_map.find(slice);
      if (it == match_slice_opflash_map.end()){
        std::vector<art::Ptr<recob::OpFlash>> flash_v;
        flash_v.push_back(flash);
        match_slice_opflash_map.insert(std::pair<art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::OpFlash>>>(slice, flash_v));
      }
      else{
        it->second.push_back(flash);
      }
    }

    for (auto it = match_slice_opflash_map.begin(); it != match_slice_opflash_map.end(); ++it){
      auto slice = it->first;
      auto flash_v = it->second;
      if (flash_v.size() > 2){
        std::cout << "more than one opflash matched to this slice!" << std::endl;
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
    } // end opt0finder loop
  }
  //  else {}
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
    _nmatch++;
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
      _match_type = -3;
      _tree->Fill();
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

    _mean_charge = (plane_charge[0] + plane_charge[1] + plane_charge[2])/3; 
    _max_charge  = plane_charge.at(bestPlane);
    _comp_charge  = plane_charge.at(bestHits);

    _slice_Q = _comp_charge;

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
        for (size_t ich=0; ich<flash_pe_v.size(); ich++) total_pe[ich] += flash_pe_v[ich];
      }
    } // end of TPC loop 

    // mask out specific channels in opdet mask
    for (size_t imask=0; imask<_opdet_mask.size(); imask++){
      total_pe.at(_opdet_mask.at(imask)) = 0;
    }

    _visibility.resize(_nchan,0);
    // calculate the photon estimates for every entry in total_pe
    CalcLight(total_pe, dir_visibility_map, ref_visibility_map, total_gamma);
    
    // fill tree variables 
    _opflash_time = flash_time;
    _dep_pe = total_pe;
    _rec_gamma = total_gamma;

    // calculate final light estimate   
    // ! TODO: final light estimate should be weighted average 
    // ! where the weights are 1/poisson_err 
    _median_gamma = CalcMedian(total_gamma);
    _mean_gamma   = CalcMean(total_gamma);
    
    if (_median_gamma!=0 && !std::isnan(_median_gamma))
      _slice_L = _median_gamma;
    else
      _slice_L = _mean_gamma;
    
    _slice_E = (_slice_L + _slice_Q)*1e-6*g4param->Wph(); // MeV, Wph = 19.5 eV   

    _true_gamma = 0; 
    _true_charge = 0;
    _true_energy = 0;

    if (_truth_validation){
      ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
      e.getByLabel(_simenergy_producer, energyDeps_h);
      std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
      
      if (!energyDeps_h.isValid() || energyDeps_h->empty())
        std::cout << "Don't have good SimEnergyDeposits!" << std::endl;
      else 
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
    else{
      _true_gamma = -9999; 
      _true_charge = -9999;
      _true_energy = -9999;
    } 
    nsuccessful_matches++;
  } // end slice loop
  _match_type=nsuccessful_matches;
  _tree->Fill();
} // end analyze


// define functions 

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
    if (i%2==0) tpc0_gamma.push_back(total_gamma[i]);
    if (i%2==1) tpc1_gamma.push_back(total_gamma[i]);
  }
  double median_gamma=0; 
  for (int tpc=0; tpc < 2; tpc++){
    double median = 0;
    auto gamma_v = (tpc==0)? tpc0_gamma:tpc1_gamma;
    std::vector<int> idx(gamma_v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
            [&](int A, int B) -> bool {
                  return gamma_v[A] < gamma_v[B];
              });
    // count number of zero entries
    int zero_counter = 0;
    for (size_t i=0; i < gamma_v.size(); i++){
      if (gamma_v.at(i) <= 0) zero_counter++;
    }
    int med_idx=0;
    if (zero_counter != int(gamma_v.size())){
      med_idx = idx.at(int((gamma_v.size()-zero_counter))/2 + zero_counter); 
      median = gamma_v.at(med_idx);
    }
    median_gamma+=median;
  }
  return median_gamma;
}

double sbnd::LightCaloProducer::CalcMean(std::vector<double> total_gamma){
  std::vector<double> tpc0_gamma; 
  std::vector<double> tpc1_gamma;
  for (size_t i=0; i<total_gamma.size(); i++){
    if (i%2==0) tpc0_gamma.push_back(total_gamma[i]);
    if (i%2==1) tpc1_gamma.push_back(total_gamma[i]);
  }
  double mean_gamma=0;
  for (int tpc=0; tpc < 2; tpc++){
    auto gamma_v = (tpc==0)? tpc0_gamma:tpc1_gamma;

    double counter=0;
    double sum=0;
    for (size_t i=0; i< gamma_v.size(); i++){
      auto gamma = gamma_v[i];
      if (gamma>0){
        counter+=1.0;
        sum+=gamma;
      }
    }
    if (sum!=0) mean_gamma+= sum/counter;
  }
  return mean_gamma;
}

DEFINE_ART_MODULE(sbnd::LightCaloProducer)
