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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// SBND includes
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <numeric>
#include <memory>

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
  // art::Ptr<recob::Slice> SelectSlice(std::vector<art::Ptr<recob::Slice>> slice_v, std::vector<float> nuscore_v, std::vector<float> fmscore_v
  //                                    std::vector<art::Ptr<recob::OpFlash>> opflash_v);
  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();

  bool MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                    std::vector<art::Ptr<recob::OpFlash>> flash_v,
                    std::vector<art::Ptr<recob::OpFlash>> &match_v);
  std::vector<double> CalcVisibility(std::vector<geo::Point_t> xyz_v, std::vector<double> charge_v);
  void CalcLight(std::vector<double> flash_pe_v, std::vector<double> visibility, std::vector<double>&total_pe_v);
  std::vector<std::string> _opflash_producer_v;
  std::string _slice_producer;
  std::string _flashmatch_producer;
  float _nuscore_cut; 
  float _fmscore_cut;

  std::unique_ptr<SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;

  TTree* _tree;
  // int _run, _subrun, _event;
  std::vector<double> _sp_x; 
  std::vector<double> _sp_x_true;
  std::vector<double> _sp_peakT;

  TTree* _tree2;
  int _nmatch=0;
  std::vector<double> _edep_time;
  double _opflash_time;

  double _true_gamma;
  double _true_charge;
  double _true_energy;
  double _slice_L; 
  double _slice_Q; 
  double _slice_E;
};


sbnd::LightCaloAna::LightCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);

  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _slice_producer = p.get<std::string>("SliceProducer");
  _flashmatch_producer = p.get<std::string>("FlashMatchProducer");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fmscore_cut = p.get<float>("fmScoreCut");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("slice_tree","");
  _tree->Branch("sp_x"    ,  "std::vector<double>", &_sp_x);
  _tree->Branch("sp_x_true", "std::vector<double>", &_sp_x_true);
  _tree->Branch("sp_peakT",  "std::vector<double>", &_sp_peakT);
  // _tree->Branch("run",             &_run,                             "run/I");
  // _tree->Branch("subrun",          &_subrun,                          "subrun/I");
  // _tree->Branch("event",           &_event,                           "event/I");


  _tree2 = fs->make<TTree>("match_tree","");
  _tree2->Branch("nmatch",        &_nmatch,       "nmatch/I");
  _tree2->Branch("edep_time",   "std::vector<double>", &_edep_time);
  _tree2->Branch("opflash_time",  &_opflash_time, "opflash_time/D");
  _tree2->Branch("true_gamma",    &_true_gamma,   "true_gamma/D");
  _tree2->Branch("true_charge",   &_true_charge,  "true_charge/D");
  _tree2->Branch("true_energy",   &_true_energy,  "true_energy/D");
  _tree2->Branch("slice_L",       &_slice_L,      "slice_L/D");
  _tree2->Branch("slice_Q",       &_slice_Q,      "slice_Q/D");
  _tree2->Branch("slice_E",       &_slice_E,      "slice_E/D");
}

void sbnd::LightCaloAna::analyze(art::Event const& e)
{
  int run    = e.id().run();
  int subrun = e.id().subRun();
  int event  = e.id().event();
  std::cout << "run: " << run <<  ", subrun: " << subrun  << ", event: " <<  event << std::endl;

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

  // tpc1 
  std::vector<art::Ptr<recob::OpFlash>> flash1_v;
  art::fill_ptr_vector(flash1_v, flash1_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op1; 
  bool found_opflash1 = MatchOpFlash(match_fm_v,flash1_v, match_op1);

  if (found_opflash0 == false && found_opflash1 == false){
    std::cout << "no opflashes matched to simpleflashes" << std::endl;
    return;
  }

  const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>&
    energyDeps(e.getValidHandle<std::vector<sim::SimEnergyDeposit>>("ionandscint"));
  
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  for (size_t n_slice=0; n_slice < match_slices_v.size(); n_slice++){
    _nmatch++;
    _sp_x.clear(); _sp_x_true.clear(); _sp_peakT.clear();
    _edep_time.clear();

    std::vector<geo::Point_t> sp_xyz;
    std::vector<double> sp_charge; // vector of charge info for charge-weighting
 
    double flash_time = -999;
    _slice_Q = 0; // total amount of charge
    _slice_L = 0; // total amount of light
    auto opflash0 = (match_op0.at(n_slice));
    auto opflash1 = (match_op1.at(n_slice));
    bool flash_in_0 = false;
    bool flash_in_1 = false;
    if (!opflash0.isNull() &&  opflash1.isNull()){
      flash_time = opflash0->Time();
      flash_in_0 = true;
    } 
    else if ( opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash1->Time(); 
      flash_in_1 = true;
    }
    else if (!opflash0.isNull() && !opflash1.isNull()){
      flash_time = opflash0->Time();
      flash_in_0 = true;
      flash_in_1 = true;
    }
    else if  (opflash0.isNull() &&  opflash1.isNull()){
      std::cout << "err! no opflashes :(" << std::endl;
      return;
    }
    
    auto slice = match_slices_v[n_slice];
    // find which plane has the most hits for this slice
    std::vector<art::Ptr<recob::Hit>> slice_hits_v = slice_to_hit.at(slice.key());
    int nhit0 = 0;
    int nhit1 = 0; 
    int nhit2 = 0;
    for (size_t i=0; i < slice_hits_v.size(); i++){
      auto hit = slice_hits_v[i];
      if (hit->View()==0) nhit0++;
      if (hit->View()==1) nhit1++;
      if (hit->View()==2) nhit2++;
    }
    uint plane = 2; 
    if ( nhit0 >= nhit1 && nhit0 >= nhit2) plane = 0;
    else if ( nhit1 >= nhit0 && nhit1 >= nhit2) plane = 1;
    else plane = 2;

    std::cout << "plane with most hits: " << plane << std::endl;

    // get charge information 
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfp.at(slice.key());
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];
      std::vector<art::Ptr<recob::SpacePoint>> sp_v = pfp_to_spacepoint.at(pfp.key());
      for (size_t n_sp=0; n_sp < sp_v.size(); n_sp++){
        auto sp = sp_v[n_sp];
        std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hit.at(sp.key());
        for (size_t n_hit=0; n_hit < hit_v.size(); n_hit++){
          auto hit = hit_v[n_hit];
          if (hit->View() !=plane) continue;
          const auto &position(sp->XYZ());
          // get true position information 
          const std::vector<double> xyz_true(bt_serv->HitToXYZ(clockData, hit));
          geo::Point_t xyz(position[0],position[1],position[2]);
          // std::cout << "spacepoint pos: (" << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
          sp_xyz.push_back(xyz);
          double drift_time = abs(200.0 - abs(position[0]))/(0.16); // in us, drift velocity = 0.16 cm/us 
          double atten_correction = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms
          double charge = (1/.0201293)*atten_correction*hit->Integral();
          sp_charge.push_back(charge);
          _slice_Q += charge;
          _sp_x.push_back(position[0]);
          _sp_x_true.push_back(xyz_true[0]);
          _sp_peakT.push_back(hit->PeakTime());
        }
      } // end spacepoint loop 
    } // end pfp loop
    // get total L count
    std::vector<double> visibility_map = CalcVisibility(sp_xyz,sp_charge);
    std::vector<double> total_pe(_nchan, 0);
    if (flash_in_0){
      auto flash_pe_v = opflash0->PEs();
      CalcLight(flash_pe_v, visibility_map, total_pe);
    }
    if (flash_in_1){
      auto flash_pe_v = opflash1->PEs();
      CalcLight(flash_pe_v, visibility_map, total_pe);
    }
    // calculate the average: 
    double sum_gamma = 0;
    double counter = 0.0;
    for (auto i : total_pe){
      if (i>0){
        counter+= 1.0;
        sum_gamma += i;
      }
    }
    // TO-DO: keep the amount of light as the average amount of light? 
    if (counter != 0)
      _slice_L = sum_gamma/counter;
    else
      _slice_L = 0;
    _true_gamma = 0; 
    _true_charge = 0;
    _true_energy = 0;
    std::cout << "flash time: " << flash_time << std::endl;
    for (const sim::SimEnergyDeposit& energyDep:*energyDeps){
      const double time = energyDep.Time() * 1e-3; 
      if (abs(time-flash_time) < 1){
        _true_gamma  += energyDep.NumPhotons()/0.03; 
        _true_charge += energyDep.NumElectrons(); 
        _true_energy += energyDep.Energy(); 
        _edep_time.push_back(time);
      }
    }
    _slice_E = (_slice_L + _slice_Q)*19.5*1e-6; 

    std::cout << "true gamma: " <<  _true_gamma << std::endl;
    std::cout << "calc gamma: " << _slice_L << std::endl;
    std::cout << "true electrons: " << _true_charge << std::endl;
    std::cout << "calc electrons: " << _slice_Q << std::endl;

    std::cout << "true deposited energy: " << _true_energy << std::endl;
    std::cout << "calc deposited energy: " << _slice_E << std::endl;

    std::cout << "ratio of gamma (calc/true):    " << _slice_L/_true_gamma << std::endl;
    std::cout << "ratio of electron (calc/true): " << _slice_Q/_true_charge << std::endl;
    std::cout << "ratio of energy (calc/true):   " << _slice_E/_true_energy << std::endl;
    std::cout << "fractional energy difference:  " << (_true_energy - _slice_E)/_true_energy << std::endl;
  
    _opflash_time = flash_time;
    
    _tree->Fill();
    _tree2->Fill();
  } // end slice loop
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
      if (abs( opflash->Time() - match_time) < 0.05){
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
  double sum_charge = std::accumulate(charge_v.begin(), charge_v.end(), 0.0);

  for (size_t i=0; i<xyz_v.size(); i++){
    geo::Point_t const xyz = xyz_v[i];
    auto charge = charge_v[i];
    std::vector<double> direct_visibility;
    std::vector<double> reflect_visibility;
    _semi_model->detectedDirectVisibilities(direct_visibility, xyz);
    _semi_model->detectedReflectedVisibilities(reflect_visibility, xyz);
    if (visibility.size() != direct_visibility.size()) std::cout << "mismatch of visibility vector size" << std::endl;

    // sum direct and reflected, weight by charge, and add to total visibility map
    for (size_t ivis=0; ivis<direct_visibility.size(); ivis++){
      visibility[ivis] += charge*(direct_visibility[ivis] + reflect_visibility[ivis]);
    }    
  } // end spacepoint loop
  // normalize by the total charge 
  std::transform(visibility.begin(), visibility.end(), 
                 visibility.begin(), [sum_charge](double &k){ return k/sum_charge; });
  return visibility;
}
void sbnd::LightCaloAna::CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> visibility,
                                   std::vector<double> &total_pe_v){
  for (size_t ichan = 0; ichan < flash_pe_v.size(); ichan++){
    auto pe = flash_pe_v[ichan];
    if((pe == 0) || std::isinf(1/visibility[ichan]))
      continue;
    total_pe_v[ichan] += (1/0.03)*pe*(1/visibility[ichan]); 
  }
}


DEFINE_ART_MODULE(sbnd::LightCaloAna)
