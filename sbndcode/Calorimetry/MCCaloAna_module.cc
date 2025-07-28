////////////////////////////////////////////////////////////////////////
// Class:       MCCaloAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        MCCaloAna_module.cc
//
// Generated at Tue Apr 29 12:33:40 2025 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
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

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "art_root_io/TFileService.h"
#include "TFile.h"
#include "TTree.h"

#include <algorithm>

namespace sbnd {
  class MCCaloAna;
}


class sbnd::MCCaloAna : public art::EDAnalyzer {
public:
  explicit MCCaloAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCCaloAna(MCCaloAna const&) = delete;
  MCCaloAna(MCCaloAna&&) = delete;
  MCCaloAna& operator=(MCCaloAna const&) = delete;
  MCCaloAna& operator=(MCCaloAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  TTree* _tree;

  int _run;
  int _subrun;
  int _event;

  std::vector<double> _nu_E;
  std::vector<int> _nu_CCNC;
  std::vector<double> _true_gamma;
  std::vector<double> _true_charge;
  std::vector<double> _true_energy;

  std::vector<double> _true_gamma_tpc;
  std::vector<double> _true_measphotons;
  std::vector<double> _true_resimphotons;

  std::vector<double> _chargecorr_gamma;
  std::vector<double> _chargecorr_visibility;
  std::vector<double> _lightcorr_gamma; 
  std::vector<double> _lightcorr_gamma_resim;
  std::vector<double> _lightcorr_visibility;

  std::vector<float> _opdet_vuv_eff;
  std::vector<float> _opdet_vis_eff;

  std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;
  opdet::sbndPDMapAlg _opdetmap; //map for photon detector types
  unsigned int _nchan = _opdetmap.size();


  std::vector<std::vector<double>> CalcVisibility(std::vector<geo::Point_t> xyz_v, 
                                                  std::vector<double> charge_v);
  std::vector<std::vector<double>> CalcPE(std::vector<geo::Point_t> xyz_v, 
                                                  std::vector<double> light_v);
  std::vector<std::vector<double>> FillSimPhotonsLite(std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> photonHandle_list);
  void CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> dir_visibility,
                                   std::vector<double> ref_visibility,
                                   std::vector<double> &total_gamma_v,
                                   std::vector<double> &total_visibility_v);
};

sbnd::MCCaloAna::MCCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);
  _opdet_vuv_eff   = p.get<std::vector<float>>("OpDetVUVEfficiencies");
  _opdet_vis_eff   = p.get<std::vector<float>>("OpDetVISEfficiencies");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("mccalo_tree","");
  _tree->Branch("run",             &_run,        "run/I");
  _tree->Branch("subrun",          &_subrun,     "subrun/I");
  _tree->Branch("event",           &_event,      "event/I");
  _tree->Branch("nu_E",           "std::vector<double>",  &_nu_E);
  _tree->Branch("nu_CCNC",        "std::vector<int>",     &_nu_CCNC);
  _tree->Branch("true_gamma",      "std::vector<double>",  &_true_gamma);
  _tree->Branch("true_charge",     "std::vector<double>",  &_true_charge);
  _tree->Branch("true_energy",     "std::vector<double>",  &_true_energy);
  _tree->Branch("true_gamma_tpc",  "std::vector<double>",  &_true_gamma_tpc);
  _tree->Branch("true_measphotons","std::vector<double>", &_true_measphotons);
  _tree->Branch("true_resimphotons","std::vector<double>", &_true_resimphotons);
  _tree->Branch("chargecorr_gamma", "std::vector<double>", &_chargecorr_gamma);
  _tree->Branch("chargecorr_visibility","std::vector<double>", &_chargecorr_visibility);
  _tree->Branch("lightcorr_gamma", "std::vector<double>", &_lightcorr_gamma);
  _tree->Branch("lightcorr_gamma_resim","std::vector<double>",&_lightcorr_gamma_resim);
  _tree->Branch("lightcorr_visibility","std::vector<double>", &_lightcorr_visibility);
}

void sbnd::MCCaloAna::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _nu_E.clear();
  _nu_CCNC.clear();
  _true_gamma.clear();
  _true_charge.clear();
  _true_energy.clear();

  _true_gamma_tpc.clear();
  _true_measphotons.clear();
  _true_resimphotons.clear();
  
  _chargecorr_gamma.clear();
  _chargecorr_visibility.clear();
  _lightcorr_gamma.clear();
  _lightcorr_gamma_resim.clear();
  _lightcorr_visibility.clear();

  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  art::Handle<std::vector<simb::MCTruth>> mctruth_handle;
  e.getByLabel("generator", mctruth_handle);

  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  if (mctruth_handle.isValid()) 
    art::fill_ptr_vector(mctruths, mctruth_handle);

  ::art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps_h;
  e.getByLabel("ionandscint", energyDeps_h);
  std::vector<art::Ptr<sim::SimEnergyDeposit>> energyDeps; 
  
  if (!energyDeps_h.isValid() || energyDeps_h->empty())
    std::cout << "Don't have good SimEnergyDeposits!" << std::endl;
  else 
    art::fill_ptr_vector(energyDeps, energyDeps_h);

  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> fLitePhotonHandle_list;
  fLitePhotonHandle_list = e.getMany<std::vector<sim::SimPhotonsLite>>();    
  auto opdet_simphotons = FillSimPhotonsLite(fLitePhotonHandle_list);

  std::vector<std::vector<double>> edep_photons;
  std::vector<std::vector<double>> edep_electrons;
  std::vector<std::vector<double>> edep_energy;
  std::vector<std::vector<geo::Point_t>> edep_points;
  // std::vector<std::vector<double>> edep_x, edep_y, edep_z; 

  _nu_E.resize(mctruths.size(), 0);
  _nu_CCNC.resize(mctruths.size(), 0);
  _true_gamma.resize(mctruths.size(), 0);
  _true_charge.resize(mctruths.size(), 0);
  _true_energy.resize(mctruths.size(), 0);
  _true_measphotons.resize(_nchan,0);
  _true_resimphotons.resize(_nchan,0);

  _true_gamma_tpc.resize(2,0);
  
  edep_photons.resize(mctruths.size());
  edep_electrons.resize(mctruths.size());
  edep_energy.resize(mctruths.size());
  edep_points.resize(mctruths.size());

  for (size_t n_dep=0; n_dep < energyDeps.size(); n_dep++){
    auto energyDep = energyDeps[n_dep];
    const auto trackID = energyDep->TrackID();

    art::Ptr<simb::MCTruth> mctruth = piserv->TrackIdToMCTruth_P(trackID);
    if (mctruth->Origin() != simb::kBeamNeutrino) continue;

    auto it = std::find(mctruths.begin(), mctruths.end(), mctruth);
    if (it == mctruths.end()) {
      std::cout << "No matching MCTruth found for trackID: " << trackID << std::endl;
      continue;
    }

    // get the index of the mctruth in the set 
    auto mctruth_index = std::distance(mctruths.begin(), it);

    _true_gamma[mctruth_index]  += energyDep->NumPhotons();
    _true_charge[mctruth_index] += energyDep->NumElectrons();
    _true_energy[mctruth_index] += energyDep->Energy();

    edep_photons.at(mctruth_index).push_back(energyDep->NumPhotons());
    edep_electrons.at(mctruth_index).push_back(energyDep->NumElectrons());
    edep_energy.at(mctruth_index).push_back(energyDep->Energy());
    edep_points.at(mctruth_index).push_back(energyDep->Start());

    if (energyDep->Start().X() < 0)
      _true_gamma_tpc.at(0) += energyDep->NumPhotons();
    else
      _true_gamma_tpc.at(1) += energyDep->NumPhotons();
  }


  auto resim_pe  = CalcPE(edep_points.at(0),edep_photons.at(0));

  for (size_t ich=0; ich < opdet_simphotons.at(0).size(); ich++){
    std::string pd_type = _opdetmap.pdType(ich);
    if(pd_type=="xarapuca_vis" || pd_type=="xarapuca_vuv") continue;
    _true_measphotons.at(ich) = opdet_simphotons.at(0).at(ich) + opdet_simphotons.at(1).at(ich);
    _true_resimphotons.at(ich) = resim_pe.at(0).at(ich) + resim_pe.at(1).at(ich);
    // std::cout << "ich " << ich <<  " (dir sim/resim): " << opdet_simphotons.at(0).at(ich) << ", " << resim_pe.at(0).at(ich) << std::endl;
    // std::cout << "ich " << ich <<  " (ref sim/resim): " << opdet_simphotons.at(1).at(ich) << ", " << resim_pe.at(1).at(ich) << std::endl;
  }

  for (size_t nnu; nnu < mctruths.size(); nnu++){
    auto true_vis_charge = CalcVisibility(edep_points.at(nnu),edep_electrons.at(nnu));
    auto true_vis_light  = CalcVisibility(edep_points.at(nnu),edep_photons.at(nnu));

    _chargecorr_gamma.resize(_nchan,0);
    _chargecorr_visibility.resize(_nchan,0);
    _lightcorr_gamma.resize(_nchan,0);
    _lightcorr_gamma_resim.resize(_nchan,0);
    _lightcorr_visibility.resize(_nchan,0);
    CalcLight(_true_measphotons,  true_vis_charge.at(0), true_vis_charge.at(1),_chargecorr_gamma,_chargecorr_visibility);
    CalcLight(_true_measphotons,  true_vis_light.at(0), true_vis_light.at(1),_lightcorr_gamma,_lightcorr_visibility);
    CalcLight(_true_resimphotons, true_vis_light.at(0), true_vis_light.at(1),_lightcorr_gamma_resim,_lightcorr_visibility);
  }

  _tree->Fill();
}

std::vector<std::vector<double>> sbnd::MCCaloAna::CalcVisibility(std::vector<geo::Point_t> xyz_v,
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

std::vector<std::vector<double>> sbnd::MCCaloAna::CalcPE(std::vector<geo::Point_t> xyz_v,
                                                                     std::vector<double> light_v){
  // returns of two vectors (len is # of opdet) for the visibility for every opdet                                     

  std::vector<std::vector<double>> pe_v(2);
  for (size_t i=0; i<2;i++){
    pe_v.at(i).resize(_nchan,0);
  }

  for (size_t i=0; i<xyz_v.size(); i++){
    geo::Point_t const xyz = xyz_v[i];
    auto light = light_v[i];

    std::vector<double> direct_visibility;
    std::vector<double> reflect_visibility;
    _semi_model->detectedDirectVisibilities(direct_visibility, xyz);
    _semi_model->detectedReflectedVisibilities(reflect_visibility, xyz);

    for (size_t ich=0;ich<_nchan;ich++){
      pe_v.at(0).at(ich) += light*direct_visibility.at(ich);
      pe_v.at(1).at(ich) += light*reflect_visibility.at(ich);
    }
  }
  return pe_v;
}

std::vector<std::vector<double>> sbnd::MCCaloAna::FillSimPhotonsLite(std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> photonHandle_list){
  std::vector<std::string> fPDTypes{"pmt_coated", "pmt_uncoated"};
  std::vector<std::vector<double>> opdet_simphotons;
  opdet_simphotons.resize(2);
  for (size_t i=0; i<2; i++){
    opdet_simphotons.at(i).resize(_nchan,0);
  }

  for ( const art::Handle<std::vector<sim::SimPhotonsLite>>& litePhotonHandle: (photonHandle_list) ){

    std::string spLabel = litePhotonHandle.provenance()->moduleLabel();
    std::vector<std::string> fSimPhotonsModuleLabel{"pdfastsim"};
    if(std::find(fSimPhotonsModuleLabel.begin(), fSimPhotonsModuleLabel.end(), spLabel)==fSimPhotonsModuleLabel.end()) continue;

    // Reflected light
    bool reflected = (litePhotonHandle.provenance()->productInstanceName() == "Reflected");

    // Loop over the SimPhotonsLite
    for ( auto const& fLitePhotons : (*litePhotonHandle) ){
      
      int opch=fLitePhotons.OpChannel;

      std::string pd_type = _opdetmap.pdType(opch);
      // Channels not sensitive to reflected light
      if(reflected && pd_type=="xarapuca_vuv") continue;
      // Channels not sensitive to direct light
      if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;
      
      // Only save the PD types specified in the fhicl list
      if(std::find(fPDTypes.begin(), fPDTypes.end(), pd_type ) != fPDTypes.end() ){
        
        std::map<int, int> fLitePhotons_map = fLitePhotons.DetectedPhotons;
        int nphotons=0;
        
        for(auto fphoton = fLitePhotons_map.begin(); fphoton!= fLitePhotons_map.end(); fphoton++){
          nphotons+=fphoton->second;
        }

        // Fill #photons per OpChannel
        if(reflected){
          opdet_simphotons.at(1).at(opch)+=nphotons;
        }
        else{
          opdet_simphotons.at(0).at(opch)+=nphotons;
        }
      }
    }
  }

  return opdet_simphotons;
}


void sbnd::MCCaloAna::CalcLight(std::vector<double> flash_pe_v,
                                   std::vector<double> dir_visibility,
                                   std::vector<double> ref_visibility,
                                   std::vector<double> &total_gamma_v,
                                   std::vector<double> &tot_visibility_v){
  for (size_t ch = 0; ch < flash_pe_v.size(); ch++){
    auto pe = flash_pe_v[ch];
    // auto vuv_eff = _opdet_vuv_eff.at(ch);
    // auto vis_eff = _opdet_vis_eff.at(ch);
    auto vuv_eff = 1.0;
    auto vis_eff = 1.0;
    auto tot_visibility = vuv_eff*dir_visibility[ch] + vis_eff*ref_visibility[ch];
    tot_visibility_v.at(ch) = tot_visibility;
    if((pe == 0) || std::isinf(1/tot_visibility))
      continue;
    // deposited light is inverse of visibility * PE count 
    total_gamma_v[ch] += (1/tot_visibility)*pe;
  }
}

DEFINE_ART_MODULE(sbnd::MCCaloAna)
