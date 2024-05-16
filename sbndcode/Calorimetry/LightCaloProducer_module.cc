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
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// SBND includes
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/Common/Reco/LightCaloInfo.h"

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

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Matches SimpleFlash time to OpFlashes, fills match_v with succesfully matched OpFlashes
  // returns true if match was found  
  bool MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                    std::vector<art::Ptr<recob::OpFlash>> flash_v,
                    std::vector<art::Ptr<recob::OpFlash>> &match_v);

  // Returns visibility vector for all opdets given charge/position information 
  void CalcVisibility(std::vector<geo::Point_t> xyz_v, 
                      std::vector<double> charge_v,
                      std::vector<double> &dir_vis_v,
                      std::vector<double> &ref_vis_v);

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
  std::string _flashmatch_producer;
  bool _use_arapucas;
  float _nuscore_cut; 
  float _fmscore_cut;
  bool _use_all_planes;
  bool _verbose;

  float _simple_op_offset;
  float _pmt_ara_offset; 
  float _noise_thresh;
  std::vector<float> _pmt_pe_range; 

  std::vector<float> _cal_area_const; 
  std::vector<float> _opdet_dir_eff;
  std::vector<float> _opdet_ref_eff;
  float _scint_prescale;

  bool _truth_validation; 
  std::string _simenergy_producer;
  bool _truth_neutrino;

  std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;

  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();

  TTree* _tree;
  int _run, _subrun, _event;
  int _match_type;
  // match_type key: 
  /// -1: no slices passed both the nuscore and flash match score cut (in the entire event)
  /// -2: no opflashes times found in coincidence with simpleflash time (in the entire event)
  /// -3: no opflashes in coincidence with simpleflash (for this slice)
  /// -4: opflashes are below the noise threshold (for this slice)
  /// 1: successful match 

  int _pfpid; // ID of the matched slice 
  double _opflash_time; // time of matched opflash 

  std::vector<double> _dep_pe; // vector of measured photo-electron (PE), one entry = one channel 
  std::vector<double> _rec_gamma; // vector of reconstructed photon counts, one entry = one channel (plane with most hits)

  double _true_gamma;  // true photon count from all energy depositions 
  double _true_charge; // true electron count from all energy depositions 
  double _true_energy; // true deposited energy 

  std::vector<double> _charge = std::vector<double>(3); // reconstructed electron count per plane 
  std::vector<double> _light_med = std::vector<double>(3);  // median reconstructed photon count per plane
  std::vector<double> _light_avg = std::vector<double>(3);  // average reconstructed photon count per plane
  std::vector<double> _energy = std::vector<double>(3);

  double _slice_L; // reconstructed photon count (plane with most hits)
  double _slice_Q; // reconstructed electron count for plane with the most hits 
  double _slice_E; // reconstructed deposited energy, sum of slice_L and slice_Q

  double _frac_L;  // light fractional difference:  (L_{true} - L_{reco})/(L_{true})
  double _frac_Q;  // charge fractional difference: (Q_{true} - Q_{reco})/(Q_{true})
  double _frac_E;  // energy fractional difference: (E_{true} - E_{reco})/(E_{true})
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
  _flashmatch_producer = p.get<std::string>("FlashMatchProducer");
  _use_arapucas = p.get<bool>("UseArapucas");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fmscore_cut = p.get<float>("fmScoreCut");
  _use_all_planes = p.get<bool>("UseAllPlanes");
  _verbose = p.get<bool>("Verbose");

  _simple_op_offset= p.get<float>("SimpleOpFlashOffset");
  _pmt_ara_offset  = p.get<float>("PMTARAFlashOffset");
  _noise_thresh    = p.get<float>("FlashNoiseThreshold");
  _pmt_pe_range    = p.get<std::vector<float>>("PMTPERange");

  _cal_area_const  = p.get<std::vector<float>>("CalAreaConstants");
  _opdet_dir_eff   = p.get<std::vector<float>>("OpDetDirectEff");
  _opdet_ref_eff   = p.get<std::vector<float>>("OpDetReflectEff");
  _scint_prescale  = p.get<float>("ScintPreScale");

  _truth_validation   = p.get<bool>("TruthValidation");
  _simenergy_producer = p.get<std::string>("SimEnergyProducer"); 
  _truth_neutrino     = p.get<bool>("TruthNeutrino");

  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::LightCalo>>();
  produces<art::Assns<recob::Slice, sbn::LightCalo>>();
  // Call appropriate consumes<>() for any products to be retrieved by this module.

}

void sbnd::LightCaloProducer::produce(art::Event& e)
{

  std::unique_ptr<std::vector<sbn::LightCalo>> lightcalo_v (new std::vector<sbn::LightCalo>);
  std::unique_ptr<art::Assns<recob::Slice, sbn::LightCalo>> slice_lightcalo_assn_v (new art::Assns<recob::Slice, sbn::LightCalo>);

  // services 
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
  
  // std::cout<<  "trigger offset TPC: " << clock_data.TriggerOffsetTPC() * 1.6 / 10 << std::endl;
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<sim::LArG4Parameters const> g4param;
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  std::cout << "run: " << _run <<  ", subrun: " << _subrun  << ", event: " <<  _event << std::endl;
  
  // initialize tree variables 
  _match_type=0; _pfpid = -1; _opflash_time=-9999;
  _dep_pe.clear(); _rec_gamma.clear();
  _true_gamma = -9999; _true_charge = -9999; _true_energy = -9999;
  _charge.assign(3,0); _light_med.assign(3,0); _light_avg.assign(3,0); _energy.assign(3,0);
  _slice_L = -1; _slice_Q = -1; _slice_E = -1; 
  _frac_L = -9999; _frac_Q = -9999; _frac_E = -9999;

  // get slices 
  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()){
    std::cout << "don't have good slices!" << std::endl;
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_slice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    std::cout << "don't have good PFParticle!" << std::endl;
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  ::art::Handle<std::vector<recob::SpacePoint>> spacepoint_h;
  e.getByLabel(_slice_producer, spacepoint_h);
  if(!spacepoint_h.isValid() || spacepoint_h->empty()) {
    std::cout << "don't have good SpacePoints!" << std::endl;
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  auto const & flash0_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[0]);
  auto const & flash1_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[1]);
  if (!_use_arapucas && _verbose)
    std::cout << "Using PMT OpFlash only..." << std::endl;  
  if( (!flash0_h.isValid() || flash0_h->empty()) && (!flash1_h.isValid() || flash1_h->empty())) {
    std::cout << "don't have good PMT flashes from producer " << _opflash_producer_v[0] << " or "  << _opflash_producer_v[1] << std::endl;
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  // Construct the vector of Slices
  std::vector<art::Ptr<recob::Slice>> slice_v;
  art::fill_ptr_vector(slice_v, slice_h);

  // Get associations 
  art::FindManyP<recob::PFParticle> slice_to_pfp (slice_h, e, _slice_producer);
  art::FindManyP<recob::Hit>        slice_to_hit (slice_h, e, _slice_producer);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_to_meta(pfp_h, e, _slice_producer);
  art::FindManyP<sbn::SimpleFlashMatch> pfp_to_sfm (pfp_h, e, _flashmatch_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoint(pfp_h, e, _slice_producer);
  art::FindManyP<recob::Hit> spacepoint_to_hit(spacepoint_h, e, _slice_producer);

  std::vector<art::Ptr<recob::Slice>> match_slices_v; 
  std::vector<int> match_pfpid_v;
  std::vector<art::Ptr<sbn::SimpleFlashMatch>> match_fm_v;

  //------------------------------//

  //* OBTAINING VALID SLICES BEGIN *// 

  for (size_t n_slice=0; n_slice < slice_v.size(); n_slice++){
    float nu_score = -9999;
    float fm_score = -9999;
    auto slice = slice_v[n_slice];
    bool found_fm = false;
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfp.at(n_slice);
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];

      // only select the PRIMARY pfp 
      if(!pfp->IsPrimary())
        continue;
      if(_truth_neutrino && !(abs(pfp->PdgCode()) == 12 || abs(pfp->PdgCode()) == 14|| abs(pfp->PdgCode()) == 16))
        continue;
      // if primary, get nu-score 
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpmeta_v = pfp_to_meta.at(pfp->Self());
      const art::Ptr<larpandoraobj::PFParticleMetadata> pfpmeta = pfpmeta_v.front();
      larpandoraobj::PFParticleMetadata::PropertiesMap propmap = pfpmeta->GetPropertiesMap();
      if (propmap.count("NuScore")) nu_score = propmap.at("NuScore");
      else nu_score = -1;
      
      // get fm-score 
      std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v = pfp_to_sfm.at(pfp.key());
      if (fm_v.empty()){
        std::cout << "No SimpleFlashMatch objects associated with this PFP!" << std::endl;
        continue;
      }
      if (fm_v.size() > 1) std::cout << "WARNING: more than one SimpleFlashMatch for one pfp?" << std::endl;
      for (size_t n_fm=0; n_fm < fm_v.size(); n_fm++){
        auto fm = fm_v.at(n_fm);
        fm_score = fm->score.total;
        if (nu_score > _nuscore_cut && fm_score < _fmscore_cut && fm_score > 0){
          found_fm = true;
          match_fm_v.push_back(fm);
        }
      } // end flashmatch loop
      if (found_fm ==true){
        match_slices_v.push_back(slice);
        match_pfpid_v.push_back(pfp->Self());
      }
    } // end pfp loop
  } // end slice loop
  if (match_slices_v.empty() && !slice_v.empty()){
    std::cout << "no slices passed the cuts" << std::endl;
    _match_type = -1; _pfpid = -1;
    _tree->Fill();
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  if (match_slices_v.size() != match_fm_v.size()){
    std::cout << "slice and flashmatch vector length mismatch!" << std::endl;
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }

  //* OBTAINING VALID SLICES END *// 

  //------------------------------//

  //* OBTAINING OPFLASH INFORMATION BEGIN*// 
  // tpc0
  std::vector<art::Ptr<recob::OpFlash>> flash0_v;
  art::fill_ptr_vector(flash0_v, flash0_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op0_v; 
  bool found_opflash0 = MatchOpFlash(match_fm_v,flash0_v, match_op0_v);

  std::vector<art::Ptr<recob::OpFlash>> flash0_ara_v;
  // if using arapucas
  if (_use_arapucas){
    ::art::Handle<std::vector<recob::OpFlash>> flash0_ara_h;
    e.getByLabel(_opflash_ara_producer_v[0],flash0_ara_h);
    if (_verbose) std::cout << "Using PMT OpFlash + X-ARAPUCA OpFlash..." << std::endl;
    if (!flash0_ara_h.isValid() || flash0_ara_h->empty()) {
      std::cout << "don't have good X-ARAPUCA flashes from producer " << _opflash_ara_producer_v[0] << std::endl;
    }
    else
      art::fill_ptr_vector(flash0_ara_v, flash0_ara_h);
  }
  // tpc1 
  std::vector<art::Ptr<recob::OpFlash>> flash1_v;
  art::fill_ptr_vector(flash1_v, flash1_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op1_v; 
  bool found_opflash1 = MatchOpFlash(match_fm_v,flash1_v, match_op1_v);

  std::vector<art::Ptr<recob::OpFlash>> flash1_ara_v;
  if (_use_arapucas){
    ::art::Handle<std::vector<recob::OpFlash>> flash1_ara_h;
    e.getByLabel(_opflash_ara_producer_v[1],flash1_ara_h);
    if (_verbose) std::cout << "Using PMT OpFlash + X-ARAPUCA OpFlash..." << std::endl;
    if (!flash1_ara_h.isValid() || flash1_ara_h->empty()) {
      std::cout << "don't have good X-ARAPUCA flashes from producer " << _opflash_ara_producer_v[1] << std::endl;
    }
    else 
      art::fill_ptr_vector(flash1_ara_v, flash1_ara_h);
  }
  // if no opflashes found in either TPC
  if (found_opflash0 == false && found_opflash1 == false){
    std::cout << "No OpFlashes matched to SimpleFlashes (for the entire event)" << std::endl;
    for (size_t n_slice=0; n_slice < match_pfpid_v.size(); n_slice++){
      _match_type = -2;
      _pfpid = match_pfpid_v.at(n_slice);
      _tree->Fill();
    }
    e.put(std::move(lightcalo_v));
    e.put(std::move(slice_lightcalo_assn_v));
    return;
  }
  //* OBTAINING OPFLASH INFORMATION END *// 

  for (size_t n_slice=0; n_slice < match_slices_v.size(); n_slice++){
    // initialize tree variables 
    _pfpid = match_pfpid_v.at(n_slice); 

    std::vector<std::vector<geo::Point_t>> sp_xyz(3); // position info of each spacepoint per plane 
    std::vector<std::vector<double>> sp_charge(3); // vector of charge info for charge-weighting for each plane

    if (_verbose)
      std::cout << "Reconstructing slice " << n_slice << std::endl;

    double flash_time = -999;
    auto opflash0 = (match_op0_v.at(n_slice));
    auto opflash1 = (match_op1_v.at(n_slice));
    bool flash_in_0 = false;
    bool flash_in_1 = false;
    if (!opflash0.isNull()){
      if (opflash0->TotalPE() > _noise_thresh)
        flash_in_0 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 0 (" << opflash0->TotalPE() << ") below noise threshold ..." << std::endl;
    } 
    if (!opflash1.isNull()){
      if (opflash1->TotalPE() > _noise_thresh)
        flash_in_1 = true;
      else if (_verbose)
        std::cout << "Flash Total PE in TPC 1 (" << opflash1->TotalPE() << ") below noise threshold ... Skipping" << std::endl;
    } 
    if  (!flash_in_0 && !flash_in_1){
      if (opflash0.isNull() || opflash1.isNull()){
        std::cout << "No OpFlashes matched with SimpleFlash objects (for this slice)" << std::endl;
        _match_type = -3;
      }
      else if (opflash0->TotalPE() < _noise_thresh|| opflash1->TotalPE() < _noise_thresh){
        std::cout << "All OpFlashes are below noise threshold..." << std::endl;
        _match_type = -4;
      }
       _opflash_time=-9999;
      _dep_pe.clear(); _rec_gamma.clear();
      _true_gamma = -9999; _true_charge = -9999; _true_energy = -9999;
      _charge.assign(3,0); _light_med.assign(3,0); _light_avg.assign(3,0); _energy.assign(3,0);
      _slice_L = -1; _slice_Q = -1; _slice_E = -1; 
      _frac_L = -9999; _frac_Q = -9999; _frac_E = -9999;
      _tree->Fill();
      continue;
    }

    if (flash_in_0) flash_time = opflash0->Time();
    else if (flash_in_1) flash_time = opflash1->Time();

    _slice_Q = 0; // total amount of charge
    _slice_L = 0; // total amount of light
    _slice_E = 0; 
    
    auto slice = match_slices_v[n_slice];
    // sum charge information (without position info) for Q, higher completeness 
    // find which plane has the most integrated charge for this slice
    std::vector<art::Ptr<recob::Hit>> slice_hits_v = slice_to_hit.at(slice.key());
    std::vector<double> plane_charge{0.,0.,0.};
    std::vector<int>    plane_hits{0,0,0};
    for (size_t i=0; i < slice_hits_v.size(); i++){
      auto hit = slice_hits_v[i];
      // TODO: fix hardcoded beam readout time
      auto drift_time = (hit->PeakTime() - 500)*0.5; // assuming TPC beam readout starts at 500 ticks, conversion = 0.5 us/tick  
      double atten_correction = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)
      auto hit_plane = hit->View();
      plane_charge.at(hit_plane) += hit->Integral()*atten_correction*(1/_cal_area_const.at(hit_plane));
      plane_hits.at(hit_plane)++; 
    }
    int bestPlane =  std::max_element(plane_hits.begin(), plane_hits.end()) - plane_hits.begin();
    _slice_Q  = plane_charge.at(bestPlane);
    _charge = plane_charge;
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
          auto hit_plane = hit->View();
          // if not the best plane and we're only using the best plane
          if ( hit_plane != bestPlane && !_use_all_planes)
              continue;
          else{
            // get position
            const auto &position(sp->XYZ());  
            geo::Point_t xyz(position[0],position[1],position[2]);
            // correct for e- attenuation 
            double drift_time = ((2.0*geom->DetHalfWidth()) - abs(position[0]))/(det_prop.DriftVelocity()); // cm / (cm/us) 
            double atten_correction = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)
            double charge = (1/_cal_area_const.at(hit_plane))*atten_correction*hit->Integral();

            sp_xyz.at(hit_plane).push_back(xyz);
            sp_charge.at(hit_plane).push_back(charge);
          }
        }
      } // end spacepoint loop 
    } // end pfp loop
    std::vector<double> total_pe(_nchan,0.);  // contains PE info from both TPCs: PMTs and ARA (if applicable)

    // combining flash PE information from separate TPCs into a single vector 
    for (int tpc=0; tpc<2; tpc++){
      bool found_flash = (tpc==0)? flash_in_0 : flash_in_1; 
      if (found_flash){
        auto opflash_pe = (tpc==0)? opflash0->PEs() : opflash1->PEs();
        for (size_t ich=0; ich < opflash_pe.size(); ich++) total_pe[ich] += opflash_pe[ich];
        // if using arapucas, need to combine PMT and arapuca PE information into a single vector 
        if (_use_arapucas){
          auto flash_ara_v = (tpc==0)? flash0_ara_v : flash1_ara_v;
          if (flash_ara_v.empty()) // if there are no valid arapuca flashes
            std::cout << "Using PMT information Only..." << std::endl;
          else{
            for (size_t nara=0; nara < flash_ara_v.size(); nara++){
              auto const &flash_ara = *flash_ara_v[nara];
              if (abs(flash_time-flash_ara.Time()) < _pmt_ara_offset){
                if (_verbose)
                  std::cout << "Combining PMT+XARA Flashes with PMT time: " << flash_time << ", ARA time: " << flash_ara.Time() << std::endl;            
                for (size_t ich=0; ich < (flash_ara.PEs()).size(); ich++)
                  total_pe.at(ich) += (flash_ara.PEs()).at(ich);
                break;
              } // end if arapuca and pmt flash time match
            } // end of arapuca flash loop
          }
        } // end of arapuca if 
      } // end of found flash if 
    } // end of TPC loop
    // fill tree variables 
    _opflash_time = flash_time;
    _dep_pe = total_pe;

    // create variables to push into the data product (not ordered by typical plane number)
    std::vector<double> lightcalo_charge;
    std::vector<double> lightcalo_light;
    std::vector<double> lightcalo_energy;
    std::vector<int> lightcalo_plane;

    std::vector<std::vector<double>> total_gamma(3, std::vector<double>(_nchan,0.));
    for (int nplane = 0; nplane<3; nplane++){
      // get total L count
      if (!sp_xyz.at(nplane).empty()){
        std::vector<double> dir_vis_v(_nchan, 0.); // direct visibility 
        std::vector<double> ref_vis_v(_nchan, 0.); // reflected visibility
        CalcVisibility(sp_xyz.at(nplane),sp_charge.at(nplane), dir_vis_v, ref_vis_v);
        // calculate the photon estimates for every entry in total_pe
        CalcLight(total_pe, dir_vis_v, ref_vis_v, total_gamma.at(nplane));
        _light_med.at(nplane) = CalcMedian(total_gamma.at(nplane));
        _light_avg.at(nplane) = CalcMedian(total_gamma.at(nplane));
        _energy.at(nplane) = (_charge.at(nplane) + _light_med.at(nplane))*1e-6*g4param->Wph();
        double gamma_avg = _light_avg.at(nplane);
        double gamma_med = _light_med.at(nplane);
        if (nplane==bestPlane){
          _rec_gamma = total_gamma.at(nplane);
          if (gamma_med!=0 && !std::isnan(gamma_med)) _slice_L = gamma_med;
          else _slice_L = gamma_avg;
          // if the plane is the bestPlane, insert it at the front
          lightcalo_charge.insert(lightcalo_charge.begin(),_charge.at(nplane));
          lightcalo_light.insert(lightcalo_light.begin(),gamma_med);
          lightcalo_energy.insert(lightcalo_energy.begin(),_energy.at(nplane));
          lightcalo_plane.insert(lightcalo_plane.begin(),nplane);
        }
        else{
          lightcalo_charge.push_back(_charge.at(nplane));
          lightcalo_light.push_back(gamma_med);
          lightcalo_energy.push_back(_energy.at(nplane));
          lightcalo_plane.push_back(nplane);
        }
      }
      else{
        _light_med.at(nplane) = -1.; _light_avg.at(nplane) = -1.; _energy.at(nplane) = -1.;
        lightcalo_charge.push_back(-1.); lightcalo_light.push_back(-1.); lightcalo_energy.push_back(-1.);
        lightcalo_plane.push_back(-1);  
      }
    }
    _slice_E = _energy.at(bestPlane);

    // sbn::LightCalo lightcalo(lightcalo_charge, lightcalo_light,lightcalo_energy,
    //                          lightcalo_plane,flash_time);

    // lightcalo_v->push_back(lightcalo);
    // util::CreateAssn(*this, e, *lightcalo_v, slice, *slice_lightcalo_assn_v);

    ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
    e.getByLabel(_simenergy_producer, energyDeps_h);
    std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
    
    if (_truth_validation) TruthValidation(e, piserv, flash_time); 
    else{
      _true_gamma = -9999; _true_charge = -9999; _true_energy = -9999;
      _frac_L = -9999; _frac_Q = -9999; _frac_E = -9999;
    }
    _match_type = 1;
    _tree->Fill();
  } // end slice loop
  e.put(std::move(lightcalo_v));
  e.put(std::move(slice_lightcalo_assn_v));
}

void sbnd::LightCaloProducer::beginJob()
{
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("slice_tree","");

  _tree->Branch("run",           &_run,          "run/I");
  _tree->Branch("subrun",        &_subrun,       "subrun/I");
  _tree->Branch("event",         &_event,        "event/I"); 
  _tree->Branch("match_type",    &_match_type,   "match_type/I");
  _tree->Branch("pfpid",         &_pfpid,        "pfpid/I");
  _tree->Branch("opflash_time",  &_opflash_time, "opflash_time/D");

  _tree->Branch("rec_gamma",    "std::vector<double>", &_rec_gamma);
  _tree->Branch("dep_pe",       "std::vector<double>", &_dep_pe);

  _tree->Branch("true_gamma",    &_true_gamma,   "true_gamma/D");
  _tree->Branch("true_charge",   &_true_charge,  "true_charge/D");
  _tree->Branch("true_energy",   &_true_energy,  "true_energy/D");

  _tree->Branch("charge_plane0", &_charge[0],    "charge_plane0/D");
  _tree->Branch("charge_plane1", &_charge[1],    "charge_plane1/D");
  _tree->Branch("charge_plane2", &_charge[2],    "charge_plane2/D");
  
  _tree->Branch("light_med_plane0", &_light_med[0] ,"light_med_plane0/D");
  _tree->Branch("light_med_plane1", &_light_med[1] ,"light_med_plane1/D");
  _tree->Branch("light_med_plane2", &_light_med[2] ,"light_med_plane2/D");

  _tree->Branch("light_avg_plane0", &_light_avg[0] ,"light_avg_plane0/D");
  _tree->Branch("light_avg_plane1", &_light_avg[1] ,"light_avg_plane1/D");
  _tree->Branch("light_avg_plane2", &_light_avg[2] ,"light_avg_plane2/D");

  _tree->Branch("energy_plane0", &_energy[0],    "energy_plane0/D");
  _tree->Branch("energy_plane1", &_energy[1],    "energy_plane1/D");
  _tree->Branch("energy_plane2", &_energy[2],    "energy_plane2/D");

  _tree->Branch("slice_L",       &_slice_L,      "slice_L/D");
  _tree->Branch("slice_Q",       &_slice_Q,      "slice_Q/D");
  _tree->Branch("slice_E",       &_slice_E,      "slice_E/D");
  _tree->Branch("frac_L",        &_frac_L,       "frac_L/D");
  _tree->Branch("frac_Q",        &_frac_Q,       "frac_Q/D");
  _tree->Branch("frac_E",        &_frac_E,       "frac_E/D");
}

void sbnd::LightCaloProducer::endJob()
{
  // Implementation of optional member function here.
}

// define functions 

bool sbnd::LightCaloProducer::MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                                      std::vector<art::Ptr<recob::OpFlash>> flash_v,
                                      std::vector<art::Ptr<recob::OpFlash>> &match_v){
  bool any_match = false; 
  for (size_t ifm=0; ifm<fm_v.size();ifm++){
    auto fm = fm_v[ifm];
    art::Ptr<recob::OpFlash> nullOpFlash;
    bool found_match = false;
    auto match_time = fm->time;
    for (size_t iop=0; iop<flash_v.size(); iop++){
      auto opflash = flash_v[iop];
      // only select opflashes that are coincident with the SimpleFlash t0
      if (abs(opflash->Time() - match_time) < _simple_op_offset){
        found_match = true;
        any_match = true;
        match_v.push_back(opflash);
        break;
      }
    }
    if (found_match == false) match_v.push_back(nullOpFlash);
  }
  if (match_v.size() != fm_v.size()) std::cout << "mismatched opflash and simpleflash vector sizes!" << std::endl;
  return any_match;
}


void sbnd::LightCaloProducer::CalcVisibility(std::vector<geo::Point_t> xyz_v,
                                             std::vector<double> charge_v,
                                             std::vector<double> &dir_vis_v,
                                             std::vector<double> &ref_vis_v){
  // returns of a vector (len is # of opdet) for the visibility for every opdet                                     
  if (xyz_v.size() != charge_v.size()) std::cout << "spacepoint coord and charge vector size mismatch" << std::endl;
  double sum_charge0 = 0;
  double sum_charge1 = 0;

  for (size_t i=0; i<xyz_v.size(); i++){
    geo::Point_t const xyz = xyz_v[i];
    auto charge = charge_v[i];

    if (xyz.X() < 0) sum_charge0+= charge;
    else sum_charge1 += charge; 

    std::vector<double> direct_visibility;
    std::vector<double> reflect_visibility;
    _semi_model->detectedDirectVisibilities(direct_visibility, xyz);
    _semi_model->detectedReflectedVisibilities(reflect_visibility, xyz);
    if (dir_vis_v.size() != direct_visibility.size() || ref_vis_v.size() != direct_visibility.size()  )
      std::cout << "mismatch of visibility vector size" << std::endl;

    // weight by charge
    for (size_t ch=0; ch<direct_visibility.size(); ch++){
      if (_opdetmap.isPDType(ch, "pmt_uncoated") || _opdetmap.isPDType(ch, "xarapuca_vis") || _opdetmap.isPDType(ch, "pmt_coated"))
        ref_vis_v.at(ch) += charge*reflect_visibility[ch];
      if (_opdetmap.isPDType(ch, "xarapuca_vuv") || _opdetmap.isPDType(ch, "pmt_coated"))
        dir_vis_v.at(ch) += charge*direct_visibility[ch];
    }    
  } // end spacepoint loop
  // normalize by the total charge in each TPC
  for (size_t ch=0; ch < dir_vis_v.size(); ch++){
    if (ch%2 == 0){ 
      if (sum_charge0 !=0) {
        dir_vis_v[ch] /= sum_charge0; 
        ref_vis_v[ch] /= sum_charge0;
      }
      else{ dir_vis_v[ch] = 0; ref_vis_v[ch] = 0;}
      }
    if (ch%2 == 1){ 
      if (sum_charge1 !=0){
        dir_vis_v[ch] /= sum_charge1; 
        ref_vis_v[ch] /= sum_charge1;
      }
      else { dir_vis_v[ch] = 0; ref_vis_v[ch] = 0;}
    }
  }
}
void sbnd::LightCaloProducer::CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> dir_visibility,
                                   std::vector<double> ref_visibility,
                                   std::vector<double> &total_gamma_v){
  for (size_t ch = 0; ch < flash_pe_v.size(); ch++){
    auto pe = flash_pe_v[ch];
    double efficiency = 0.03; 
    if((pe == 0) || std::isinf(abs(1/(dir_visibility[ch]+ref_visibility[ch]))))
      continue;
    if ( pe < _pmt_pe_range.at(0) || pe > _pmt_pe_range.at(1))
      continue;
    if (_opdetmap.isPDType(ch, "pmt_uncoated"))
      efficiency = _opdet_dir_eff[0]; 
    else if (_opdetmap.isPDType(ch, "pmt_coated"))
      efficiency = _opdet_dir_eff[1];
    else if (_opdetmap.isPDType(ch, "xarapuca_vis"))
      efficiency = _opdet_dir_eff[2];
    else if (_opdetmap.isPDType(ch, "xarapuca_vuv"))
      efficiency = _opdet_dir_eff[3];
    // deposited light is inverse of efficiency * inverse of visibility * PE count 
    total_gamma_v[ch] += (1/efficiency)*(1/(dir_visibility[ch]+ref_visibility[ch]))*pe; 
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
    median_gamma+=median; // sum the two tpc median light
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
    if (counter!=0) mean_gamma+= sum/counter;
  }
  return mean_gamma;
}

void sbnd::LightCaloProducer::TruthValidation(art::Event& e, art::ServiceHandle<cheat::ParticleInventoryService> piserv, double flash_time){
  _true_gamma = 0; 
  _true_charge = 0;
  _true_energy = 0;

  ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
  e.getByLabel(_simenergy_producer, energyDeps_h);
  std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
  
  if (!energyDeps_h.isValid() || energyDeps_h->empty()){
    std::cout << "Don't have good SimEnergyDeposits!" << std::endl;
    _true_gamma = -9999; _true_charge = -9999; _true_energy = -9999;
    _frac_L = -9999; _frac_Q = -9999; _frac_E = -9999;
    return;
  }
  else art::fill_ptr_vector(energyDeps, energyDeps_h);
    
  for (size_t n_dep=0; n_dep < energyDeps.size(); n_dep++){
    auto energyDep = energyDeps[n_dep];
    const auto trackID = energyDep->TrackID();
    const double time = energyDep->Time() * 1e-3; // us 

    art::Ptr<simb::MCTruth> mctruth = piserv->TrackIdToMCTruth_P(trackID);

    // if validating a neutrino event, check that the origin of the deposit is from a beam neutrino
    // if validating a non-neutrino event (e.g. cosmics), match the deposit time with the flash time (may not be complete)
    if ((_truth_neutrino && mctruth->Origin()==simb::kBeamNeutrino ) || 
        (!_truth_neutrino && abs(time-flash_time) < 1)){
      // note: we divide by the prescale because NumPhotons() stored in simulation has the scint prescale applied 
      _true_gamma  += energyDep->NumPhotons()/_scint_prescale;
      _true_charge += energyDep->NumElectrons(); 
      _true_energy += energyDep->Energy(); 
    }       
  }
  _frac_L = _true_gamma > 0? (_true_gamma - _slice_L)/_true_gamma : -9999; 
  _frac_Q = _true_charge > 0? (_true_charge - _slice_Q)/_true_charge : -9999;
  _frac_E = _true_energy > 0? (_true_energy - _slice_E)/_true_energy : -9999;

  if (_verbose){
    std::cout << "ratio of gamma (median/true):  " << _slice_L/_true_gamma << std::endl;
    std::cout << "ratio of electron (comp/true): " << _slice_Q/_true_charge << std::endl;
    std::cout << "ratio of energy (calc/true):   " << _slice_E/_true_energy << std::endl;
  }
}

DEFINE_ART_MODULE(sbnd::LightCaloProducer)
