////////////////////////////////////////////////////////////////////////
// Class:       LightCaloAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        LightCaloAna_module.cc
//
// Generated at Fri Oct 14 13:43:49 2022 by Lynn Tung using cetskelgen
// from  version .
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
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/sim.h"

// SBND includes
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <numeric>
#include <memory>
#include <algorithm> // sort 

namespace sbnd {
  class LightCaloAna;
}


class sbnd::LightCaloAna : public art::EDAnalyzer {
public:
  explicit LightCaloAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LightCaloAna(LightCaloAna const&) = delete;
  LightCaloAna(LightCaloAna&&) = delete;
  LightCaloAna& operator=(LightCaloAna const&) = delete;
  LightCaloAna& operator=(LightCaloAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // define functions 
  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();

  bool MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                    std::vector<art::Ptr<recob::OpFlash>> flash_v,
                    std::vector<art::Ptr<recob::OpFlash>> &match_v);
  std::vector<double> CalcVisibility(std::vector<geo::Point_t> xyz_v, std::vector<double> charge_v);
  void CalcLight(std::vector<double> flash_pe_v, std::vector<double> visibility, std::vector<double>&total_pe_v);
  double CalcMedian(std::vector<double> total_light);
  double CalcMean(std::vector<double> total_light);
  std::vector<std::string> _opflash_producer_v;
  std::vector<std::string> _opflash_ara_producer_v;
  std::string _slice_producer;
  std::string _flashmatch_producer;
  bool _use_arapucas;
  float _nuscore_cut; 
  float _fmscore_cut;
  std::vector<float> _cal_area_const; 

  std::unique_ptr<SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;


  int _run, _subrun, _event;

  TTree* _tree;
  int _match_type;

  TTree* _tree2;
  int _nmatch=0;
  int _pfpid; 
  double _opflash_time;

  std::vector<double> _dep_pe; 
  std::vector<double> _rec_gamma; 

  double _true_gamma;
  double _true_charge;
  double _true_energy;

  double _median_gamma;
  double _mean_gamma;

  double _mean_charge; 
  double _max_charge;
  double _comp_charge;

  double _slice_L; 
  double _slice_Q; 
  double _slice_E;

  double _frac_L; 
  double _frac_Q; 
  double _frac_E;

};


sbnd::LightCaloAna::LightCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);

  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _opflash_ara_producer_v = p.get<std::vector<std::string>>("OpFlashAraProducers");
  _slice_producer = p.get<std::string>("SliceProducer");
  _flashmatch_producer = p.get<std::string>("FlashMatchProducer");
  _use_arapucas = p.get<bool>("UseArapucas");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fmscore_cut = p.get<float>("fmScoreCut");
  _cal_area_const    = p.get<std::vector<float>>("CalAreaConstants");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("slice_tree","");
  _tree->Branch("run",             &_run,        "run/I");
  _tree->Branch("subrun",          &_subrun,     "subrun/I");
  _tree->Branch("event",           &_event,      "event/I");
  _tree->Branch("match_type",      &_match_type, "match_type/I");

  _tree2 = fs->make<TTree>("match_tree","");
  _tree2->Branch("run",           &_run,          "run/I");
  _tree2->Branch("subrun",        &_subrun,       "subrun/I");
  _tree2->Branch("event",         &_event,        "event/I");
  _tree2->Branch("nmatch",        &_nmatch,       "nmatch/I");
  _tree2->Branch("pfpid",         &_pfpid,        "pfpid/I");
  _tree2->Branch("opflash_time",  &_opflash_time, "opflash_time/D");

  _tree2->Branch("rec_gamma",    "std::vector<double>", &_rec_gamma);
  _tree2->Branch("dep_pe",       "std::vector<double>", &_dep_pe);

  _tree2->Branch("true_gamma",    &_true_gamma,   "true_gamma/D");
  _tree2->Branch("true_charge",   &_true_charge,  "true_charge/D");
  _tree2->Branch("true_energy",   &_true_energy,  "true_energy/D");

  _tree2->Branch("median_gamma",  &_median_gamma, "median_gamma/D");
  _tree2->Branch("mean_gamma",    &_mean_gamma,   "mean_gamma/D");
  _tree2->Branch("mean_charge",   &_mean_charge,  "mean_charge/D");
  _tree2->Branch("max_charge",    &_max_charge,   "max_charge/D");
  _tree2->Branch("comp_charge",   &_comp_charge,  "comp_charge/D");

  _tree2->Branch("slice_L",       &_slice_L,      "slice_L/D");
  _tree2->Branch("slice_Q",       &_slice_Q,      "slice_Q/D");
  _tree2->Branch("slice_E",       &_slice_E,      "slice_E/D");
  _tree2->Branch("frac_L",        &_frac_L,       "frac_L/D");
  _tree2->Branch("frac_Q",        &_frac_Q,       "frac_Q/D");
  _tree2->Branch("frac_E",        &_frac_E,       "frac_E/D");

}

void sbnd::LightCaloAna::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  std::cout << "run: " << _run <<  ", subrun: " << _subrun  << ", event: " <<  _event << std::endl;

  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));

  // get slices 
  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()){
    std::cout << "dont have good slices!" << std::endl;
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
    std::cout << "Don't have good SpacePoints!" << std::endl;
    return;
  }

  auto const & flash0_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[0]);
  auto const & flash1_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[1]);
  if( (!flash0_h.isValid() || flash0_h->empty()) && (!flash1_h.isValid() || flash1_h->empty())) {
    std::cout << "don't have good flashes from producer " << _opflash_producer_v[0] << " or "  << _opflash_producer_v[1] << std::endl;
    return;
  }
  auto const & flash0_ara_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_ara_producer_v[0]);
  auto const & flash1_ara_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_ara_producer_v[1]);
  if( _use_arapucas && (!flash0_ara_h.isValid() || flash0_ara_h->empty()) && (!flash1_ara_h.isValid() || flash1_ara_h->empty())) {
    std::cout << "don't have good flashes from producer " << _opflash_ara_producer_v[0] << " or "  << _opflash_ara_producer_v[1] << std::endl;
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
  std::vector<art::Ptr<sbn::SimpleFlashMatch>> match_fm_v;

  for (size_t n_slice=0; n_slice < slice_v.size(); n_slice++){
    float nu_score = -9999;
    float fm_score = -9999;
    // float fm_time  = -9999;
    auto slice = slice_v[n_slice];
    bool found_fm = false;
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfp.at(n_slice);
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];

      // only select the PRIMARY pfp 
      if(!pfp->IsPrimary() && !(abs(pfp->PdgCode()) == 12 || abs(pfp->PdgCode()) == 14|| abs(pfp->PdgCode()) == 16))
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
        continue;
      }
      if (fm_v.size() > 1)
        std::cout << "more than one match for one pfp?" << std::endl;
      for (size_t n_fm=0; n_fm < fm_v.size(); n_fm++){
        auto fm = fm_v.at(n_fm);
        fm_score = fm->score.total;
        // fm_time  = fm->time;
        if (nu_score > _nuscore_cut && fm_score < _fmscore_cut && fm_score > 0){
          found_fm = true;
          match_fm_v.push_back(fm);
        }
      } // end flashmatch loop
      if (found_fm ==true) match_slices_v.push_back(slice);
    } // end pfp loop
  } // end slice loop
  if (match_slices_v.empty() && !slice_v.empty()){
    std::cout << "no slices passed the cuts" << std::endl;
    _match_type = -1;
    _tree->Fill();
    return;
  }

  if (match_slices_v.size() != match_fm_v.size()){
    std::cout << "slice and flashmatch vector length mismatch!" << std::endl;
    return;
  }

  // get relevant opflashes
  // tpc0
  std::vector<art::Ptr<recob::OpFlash>> flash0_v;
  art::fill_ptr_vector(flash0_v, flash0_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op0; 
  bool found_opflash0 = MatchOpFlash(match_fm_v,flash0_v, match_op0);

  std::vector<art::Ptr<recob::OpFlash>> flash0_ara_v;
  art::fill_ptr_vector(flash0_ara_v, flash0_ara_h);

  // tpc1 
  std::vector<art::Ptr<recob::OpFlash>> flash1_v;
  art::fill_ptr_vector(flash1_v, flash1_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op1; 
  bool found_opflash1 = MatchOpFlash(match_fm_v,flash1_v, match_op1);

  std::vector<art::Ptr<recob::OpFlash>> flash1_ara_v;
  art::fill_ptr_vector(flash1_ara_v, flash1_ara_h);

  if (found_opflash0 == false && found_opflash1 == false){
    std::cout << "no opflashes matched to simpleflashes" << std::endl;
    _match_type = -2;
    _tree->Fill();
    return;
  }

  const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>&
    energyDeps(e.getValidHandle<std::vector<sim::SimEnergyDeposit>>("ionandscint"));
  
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

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
    double noise_thresh = (_use_arapucas)? 2000 : 1000; 

    if (!opflash0.isNull() &&  opflash1.isNull()){
      flash_time = opflash0->Time();
      if (opflash0->TotalPE() > noise_thresh)
        flash_in_0 = true;
      else
        std::cout << "Flash Total PE in TPC 0 too low ("<< opflash0->TotalPE() << ")... Skipping" << std::endl;
    } 
    else if ( opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash1->Time(); 
      if (opflash1->TotalPE() > noise_thresh)
        flash_in_1 = true;
      else
        std::cout << "Flash Total PE in TPC 1 too low(" << opflash1->TotalPE() << ")... Skipping" << std::endl;
    }
    else if (!opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash0->Time();
      if (opflash0->TotalPE() > noise_thresh)
        flash_in_0 = true;
      else 
        std::cout << "Flash Total PE in TPC 0 too low ("<< opflash0->TotalPE() << ")... Skipping" << std::endl;
      if (opflash1->TotalPE() > noise_thresh)
        flash_in_1 = true;
      else
        std::cout << "Flash Total PE in TPC 1 too low(" << opflash1->TotalPE() << ")... Skipping" << std::endl;
    }
    else if  (opflash0.isNull() &&  opflash1.isNull()){
      std::cout << "No opflashes matched with SimpleFlash objects" << std::endl;
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
      auto drift_time = (hit->PeakTime() - 500)*0.5; // us 
      double atten_correction = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms
      auto hit_plane = hit->View();
      plane_charge.at(hit_plane) += hit->Integral()*atten_correction*(1/_cal_area_const.at(hit_plane));
      plane_hits.at(hit_plane)++; 
    }
    std::cout << "charge: " << plane_charge[0] << ", " << plane_charge[1] << ", " << plane_charge[2] << std::endl;
    uint bestPlane = std::max_element(plane_charge.begin(), plane_charge.end()) - plane_charge.begin(); 
    uint bestHits =  std::max_element(plane_hits.begin(), plane_hits.end()) - plane_hits.begin();

    _mean_charge = (plane_charge[0] + plane_charge[1] + plane_charge[2])/3; 
    _max_charge  = plane_charge.at(bestPlane);
    _comp_charge  = plane_charge.at(bestHits);

    _slice_Q = _comp_charge;

    double sps_Q = 0;

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
          if (hit->View() !=bestPlane) continue;
          const auto &position(sp->XYZ());
          geo::Point_t xyz(position[0],position[1],position[2]);
          // correct for e- attenuation 
          double drift_time = abs(200.0 - abs(position[0]))/(0.16); // in us, drift velocity = 0.16 cm/us 
          double atten_correction = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms
          double charge = (1/.0201293)*atten_correction*hit->Integral();
          sp_xyz.push_back(xyz);
          sp_charge.push_back(charge);

          sps_Q += charge;
        }
      } // end spacepoint loop 
    } // end pfp loop

    // get total L count
    std::vector<double> visibility_map = CalcVisibility(sp_xyz,sp_charge);
    std::vector<double> total_pe(_nchan,0.); 
    std::vector<double> total_gamma(_nchan, 0.);

    if (flash_in_0){
      std::cout << "Using OpFlash in TPC 0..." << std::endl;
      auto flash_pe_v = opflash0->PEs();
      flash_pe_v.insert(flash_pe_v.end(), {0,0,0,0,0,0});
      if (_use_arapucas){
        for (size_t nara=0; nara < flash0_ara_v.size(); nara++){
          auto const &flash0_ara = *flash0_ara_v[nara];
          if (flash0_ara.Time() - flash_time < 0.05){
            for (size_t ich=0; ich < (flash0_ara.PEs()).size(); ich++)
              flash_pe_v.at(ich) += (flash0_ara.PEs()).at(ich);
            break;
          }  
        }
      }
      for (size_t ich=0; ich<flash_pe_v.size(); ich++) total_pe[ich] += flash_pe_v[ich];
      CalcLight(flash_pe_v, visibility_map, total_gamma);
    }
    if (flash_in_1){
      auto flash_pe_v = opflash1->PEs();
      std::cout << "Using OpFlash in TPC 1..." << std::endl;
      flash_pe_v.insert(flash_pe_v.end(), {0,0,0,0,0,0});
      if (_use_arapucas){
        for (size_t nara=0; nara < flash0_ara_v.size(); nara++){
          auto const &flash0_ara = *flash0_ara_v[nara];
          if (flash0_ara.Time() - flash_time < 0.05){
            for (size_t ich=0; ich < (flash0_ara.PEs()).size(); ich++)
              flash_pe_v.at(ich) += (flash0_ara.PEs()).at(ich);
            break;
          }  
        }
      }
      for (size_t ich=0; ich<flash_pe_v.size(); ich++) total_pe[ich] += flash_pe_v[ich];
      CalcLight(flash_pe_v, visibility_map, total_gamma);
    }

    _median_gamma = CalcMedian(total_gamma);
    _mean_gamma   = CalcMean(total_gamma);
    // protect against double flash outliers: 
  
      // fill tree variables 
    _dep_pe = total_pe;
    _rec_gamma = total_gamma;
    
    _slice_L = _median_gamma;
    _slice_E = (_slice_L + _slice_Q)*19.5*1e-6; 

    // get truth information 
    _true_gamma = 0; 
    _true_charge = 0;
    _true_energy = 0;

    double _true_gamma_0 = 0;
    double _true_gamma_1 = 0; 

    for (const sim::SimEnergyDeposit& energyDep:*energyDeps){
      auto trackID = energyDep.TrackID();
      auto start_x   = energyDep.StartX(); 

      art::Ptr<simb::MCTruth> mctruth = pi_serv->TrackIdToMCTruth_P(trackID);
      if (mctruth->Origin()==simb::kBeamNeutrino){
      
        _true_gamma  += energyDep.NumPhotons()/0.03; // scint-prescale = 0.03
        _true_charge += energyDep.NumElectrons(); 
        _true_energy += energyDep.Energy(); 
        if (start_x > 0) _true_gamma_1+= energyDep.NumPhotons()/0.03;
        if (start_x < 0) _true_gamma_0+= energyDep.NumPhotons()/0.03;
      }
    }

    // }
    std::cout << "ratio of gamma (median/true):  " << _median_gamma/_true_gamma << std::endl;
    std::cout << "ratio of gamma (mean/true):    " << _mean_gamma/_true_gamma << std::endl;
    if (flash_in_0 && flash_in_1){
      std::cout << "ratio of gamma (mean/true) TPC 0: " << mean_gamma0/_true_gamma_0 << std::endl;
      std::cout << "ratio of gamma (mean/true) TPC 1: " << mean_gamma1/_true_gamma_1 << std::endl;
    }

    std::cout << "ratio of electron (mean/true): " << _mean_charge/_true_charge << std::endl;
    std::cout << "ratio of electron (max/true):  " << _max_charge/_true_charge << std::endl;
    std::cout << "ratio of electron (comp/true): " << _comp_charge/_true_charge << std::endl;

    std::cout << "ratio of energy (calc/true):   " << _slice_E/_true_energy << std::endl;
    // std::cout << "fractional energy difference:  " << (_true_energy - _slice_E)/_true_energy << std::endl;
  
    _opflash_time = flash_time;

    _frac_L = (_true_gamma - _slice_L)/_true_gamma; 
    _frac_Q = (_true_charge - _slice_Q)/_true_charge;
    _frac_E = (_true_energy - _slice_E)/_true_energy;
    _tree2->Fill();
    nsuccessful_matches++;
  } // end slice loop
  _match_type=nsuccessful_matches;
  _tree->Fill();

} // end analyze


// define functions 

bool sbnd::LightCaloAna::MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
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
      // only select opflashes that match the time and have total PE above noise levels
      if (abs( opflash->Time() - match_time) < 0.17){
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


std::vector<double> sbnd::LightCaloAna::CalcVisibility(std::vector<geo::Point_t> xyz_v,
                                                       std::vector<double> charge_v){
  // returns of a vector (len is # of opdet) for the visibility for every opdet                                     
  if (xyz_v.size() != charge_v.size()) std::cout << "spacepoint coord and charge vector size mismatch" << std::endl;

  std::vector<double> visibility(_nchan, 0);
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
    if (visibility.size() != direct_visibility.size()) std::cout << "mismatch of visibility vector size" << std::endl;

    // weight by charge
    for (size_t ch=0; ch<direct_visibility.size(); ch++){
      if (_opdetmap.isPDType(ch, "pmt_uncoated") || _opdetmap.isPDType(ch, "xarapuca_vis"))
        visibility[ch] += charge*reflect_visibility[ch];
      else if (_opdetmap.isPDType(ch, "xarapuca_vuv"))
        visibility[ch] += charge*direct_visibility[ch];
      else if (_opdetmap.isPDType(ch, "pmt_coated"))
        visibility[ch] += charge*(direct_visibility[ch] + reflect_visibility[ch]);
      
      // std::cout << "x: " << xyz.X() << "vis: " << visibility[ch] << std::endl; 
    }    
  } // end spacepoint loop
  // normalize by the total charge in each TPC
  for (size_t ch=0; ch < visibility.size(); ch++){
    if (ch%2 == 0) visibility[ch] /= sum_charge0;
    if (ch%2 == 1) visibility[ch] /= sum_charge1;
  }
  // std::transform(visibility.begin(), visibility.end(), 
  //                visibility.begin(), [sum_charge](double &k){ return k/sum_charge; });
  return visibility;
}
void sbnd::LightCaloAna::CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> visibility,
                                   std::vector<double> &total_pe_v){
  for (size_t ch = 0; ch < flash_pe_v.size(); ch++){
    auto pe = flash_pe_v[ch];
    double efficiency = 0.03; 
    if((pe == 0) || std::isinf(1/visibility[ch]))
      continue;
    if (_opdetmap.isPDType(ch, "pmt_uncoated") || _opdetmap.isPDType(ch, "pmt_coated"))
      efficiency = 0.03;
    else if (_opdetmap.isPDType(ch, "xarapuca_vuv"))
      efficiency = 0.021;
    else if (_opdetmap.isPDType(ch, "xarapuca_vis"))
      efficiency = 0.014;
    total_pe_v[ch] += (1/efficiency)*pe*(1/visibility[ch]); 
  }
}

double sbnd::LightCaloAna::CalcMedian(std::vector<double> total_gamma){
  std::vector<double> tpc0_gamma; 
  std::vector<double> tpc1_gamma;
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
    // count number of zero entries: 
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

double sbnd::LightCaloAna::CalcMean(std::vector<double> total_gamma){
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


DEFINE_ART_MODULE(sbnd::LightCaloAna)
