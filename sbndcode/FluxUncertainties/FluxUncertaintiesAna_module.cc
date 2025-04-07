////////////////////////////////////////////////////////////////////////
// Class:       FluxUncertaintiesAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        FluxUncertaintiesAna_module.cc
//
// Generated at Fri Mar  7 15:30:10 2025 by Francisco Nicolas-Arnaldos using cetskelgen
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
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

#include "TTree.h"
#include "TH1D.h"

namespace flux {
  class FluxUncertaintiesAna;
}


class flux::FluxUncertaintiesAna : public art::EDAnalyzer {
public:
  explicit FluxUncertaintiesAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FluxUncertaintiesAna(FluxUncertaintiesAna const&) = delete;
  FluxUncertaintiesAna(FluxUncertaintiesAna&&) = delete;
  FluxUncertaintiesAna& operator=(FluxUncertaintiesAna const&) = delete;
  FluxUncertaintiesAna& operator=(FluxUncertaintiesAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  std::string _mctruth_label;
  std::string _fluxeventweight_label;
  bool _save_tree;
  bool _save_histo;
  int _nUniverses;
  std::vector<std::string> _systematic_list;
  bool _save_total_sytematic;
  int _histoNBins;
  double _histoMin;
  double _histoMax;

  // Map of histograms for universes (all neutrino flavours)
  std::map<std::string, std::vector<TH1D*>> _histo_map;
  // Map of histograms for universes (by neutrino flavour)
  std::map<std::string, std::vector<TH1D*>> _histo_map_numu;
  std::map<std::string, std::vector<TH1D*>> _histo_map_numubar;
  std::map<std::string, std::vector<TH1D*>> _histo_map_nue;
  std::map<std::string, std::vector<TH1D*>> _histo_map_nuebar;

  // Central value histograms
  TH1D* _histo_central;
  TH1D* _histo_central_numu;
  TH1D* _histo_central_numubar;
  TH1D* _histo_central_nue;
  TH1D* _histo_central_nuebar;

  // TTree
  TTree* _tree;

  // Event info
  int _run, _subrun, _event;

  // Neutrino info
  double _nu_e, _nu_pdg;

  // Flux weights
  int _evtwgt_flux_nfunc; ///< Number of functions used for FLUX reweighting (multisim)
  std::vector<std::string> _evtwgt_flux_funcname; ///< Names of the functions used for FLUX reweighting (multisim)
  std::vector<std::vector<float>> _evtwgt_flux_weight; ///< Weights per function name used for FLUX reweighting (multisim)
  std::vector<float> _evtwgt_flux_oneweight; ///< Weights for FLUX reweighting (multisim) (combines all variations)
};



flux::FluxUncertaintiesAna::FluxUncertaintiesAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  _mctruth_label(p.get<std::string>("MCTruthLabel")),
  _fluxeventweight_label(p.get<std::string>("FluxEventWeightLabel")),
  _save_tree(p.get<bool>("SaveTree")),
  _save_histo(p.get<bool>("SaveHisto")),
  _nUniverses(p.get<int>("NUniverses")),
  _systematic_list(p.get<std::vector<std::string>>("SystematicList")),
  _save_total_sytematic(p.get<bool>("SaveTotalSystematic")),
  _histoNBins(p.get<int>("HistoNBins")),
  _histoMin(p.get<double>("HistoMin")),
  _histoMax(p.get<double>("HistoMax"))
  // ,
  // More initializers here.
{
  if(_save_total_sytematic){
    _systematic_list.push_back("total");
  }
}


void flux::FluxUncertaintiesAna::analyze(art::Event const& e)
{

  // Get the event info
  _run = e.run();
  _subrun = e.subRun();
  _event = e.event();

  // Get the MCTruth
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_label, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product for " << _mctruth_label << " not found." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

  // Loop over MCTruth
  for (size_t i = 0; i < mct_v.size(); i++) {

    // Get the FluxEventWeightMap
    art::FindManyP<sbn::evwgh::EventWeightMap> mct_to_fluxewm(mct_h, e, _fluxeventweight_label);
    
    std::cout<<"mct_to_fluxewm.size() " << mct_to_fluxewm.size() << std::endl;
    // Flux weights
    std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> flux_ewm_v = mct_to_fluxewm.at(i);
    
    if (flux_ewm_v.size() != 1) {
      std::cout << "EventWeightMap of " << _fluxeventweight_label << " bigger than 1?" << std::endl;
    }
    std::map<std::string, std::vector<float>> evtwgt_map = *(flux_ewm_v[0]);

    // Get the neutrino info
    _nu_e = mct_v[i]->GetNeutrino().Nu().E();
    _nu_pdg = mct_v[i]->GetNeutrino().Nu().PdgCode();

    // Fill the central value
    _histo_central->Fill(_nu_e);
    if(_nu_pdg == 14) _histo_central_numu->Fill(_nu_e);
    if(_nu_pdg == -14) _histo_central_numubar->Fill(_nu_e);
    if(_nu_pdg == 12) _histo_central_nue->Fill(_nu_e);
    if(_nu_pdg == -12) _histo_central_nuebar->Fill(_nu_e);

    std::cout << "Analyzing event " << _run << " " << _subrun << " " << _event << " with neutrino energy " << _nu_e << " and pdg " << _nu_pdg << std::endl;
    // Reset vars
    _evtwgt_flux_funcname.clear();
    _evtwgt_flux_weight.clear();
    _evtwgt_flux_oneweight.clear();

    std::vector<float> previous_weights;
    std::vector<float> final_weights;

    for(auto it : evtwgt_map) {

      std::string func_name = it.first;
      std::vector<float> weight_v = it.second;

      // Initialize the vectors to calculate the final weights
      if (previous_weights.size() == 0) {
        previous_weights.resize(weight_v.size(), 1.);
        final_weights.resize(weight_v.size(), 1.);
      }

      // Calculate the final weights
      std::transform(previous_weights.begin(), previous_weights.end(),
                       weight_v.begin(),
                       final_weights.begin(),
                       std::multiplies<float>());
      previous_weights = final_weights;

      // Fill the tree
      _evtwgt_flux_funcname.push_back(func_name);
      _evtwgt_flux_weight.push_back(weight_v);

      // Fill the histograms
      if (_save_histo) {
        for (size_t i = 0; i < weight_v.size(); i++) {
          _histo_map[func_name][i]->Fill(_nu_e, weight_v[i]);
          if(_nu_pdg == 14) _histo_map_numu[func_name][i]->Fill(_nu_e, weight_v[i]);
          if(_nu_pdg == -14) _histo_map_numubar[func_name][i]->Fill(_nu_e, weight_v[i]);
          if(_nu_pdg == 12) _histo_map_nue[func_name][i]->Fill(_nu_e, weight_v[i]);
          if(_nu_pdg == -12) _histo_map_nuebar[func_name][i]->Fill(_nu_e, weight_v[i]);
        }
      }

    }

    // Fill the final weights
    _evtwgt_flux_oneweight = final_weights;
    _evtwgt_flux_nfunc = _evtwgt_flux_funcname.size();

    // Fill the histograms for the total
    if (_save_histo && _save_total_sytematic) {
      for (size_t i = 0; i < final_weights.size(); i++) {
        _histo_map["total"][i]->Fill(_nu_e, final_weights[i]);
        if(_nu_pdg == 14) _histo_map_numu["total"][i]->Fill(_nu_e, final_weights[i]);
        if(_nu_pdg == -14) _histo_map_numubar["total"][i]->Fill(_nu_e, final_weights[i]);
        if(_nu_pdg == 12) _histo_map_nue["total"][i]->Fill(_nu_e, final_weights[i]);
        if(_nu_pdg == -12) _histo_map_nuebar["total"][i]->Fill(_nu_e, final_weights[i]);
      }
    }

    // Fill the tree
    if (_save_tree) _tree->Fill();
  
  }

}


void flux::FluxUncertaintiesAna::beginJob()
{

  art::ServiceHandle<art::TFileService> fs;

  if(_save_tree){
    _tree = fs->make<TTree>("FluxUncertTree","");

    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("subrun", &_subrun, "subrun/I");
    _tree->Branch("event", &_event, "event/I");

    _tree->Branch("nu_e", &_nu_e, "nu_e/D");
    _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/D");

    _tree->Branch("evtwgt_flux_nfunc", &_evtwgt_flux_nfunc, "evtwgt_flux_nfunc/I");
    _tree->Branch("evtwgt_flux_funcname", "std::vector<std::string>", &_evtwgt_flux_funcname);
    _tree->Branch("evtwgt_flux_weight", "std::vector<std::vector<float>>", &_evtwgt_flux_weight);
    _tree->Branch("evtwgt_flux_oneweight", "std::vector<float>", &_evtwgt_flux_oneweight);
  }

  // Initialize the histograms
  if (_save_histo) {

    // Central value
    _histo_central = fs->make<TH1D>("central_value", "central_value", _histoNBins, _histoMin, _histoMax);
    _histo_central_numu = fs->make<TH1D>("central_value_numu", "central_value_numu", _histoNBins, _histoMin, _histoMax);
    _histo_central_numubar = fs->make<TH1D>("central_value_numubar", "central_value_numubar", _histoNBins, _histoMin, _histoMax);
    _histo_central_nue = fs->make<TH1D>("central_value_nue", "central_value_nue", _histoNBins, _histoMin, _histoMax);
    _histo_central_nuebar = fs->make<TH1D>("central_value_nuebar", "central_value_nuebar", _histoNBins, _histoMin, _histoMax);

    // Initialize the histograms for the systematics
    for (auto syst : _systematic_list) {
      std::vector<TH1D*> hV;
      std::vector<TH1D*> hV_numu;
      std::vector<TH1D*> hV_numubar;
      std::vector<TH1D*> hV_nue;
      std::vector<TH1D*> hV_nuebar;

      art::TFileDirectory tfdir  = fs->mkdir( syst.c_str() );
      for (int i = 0; i < _nUniverses; i++) {
        std::string name;

        name = syst + "_universe_" + std::to_string(i);
        TH1D* h = tfdir.make<TH1D>(name.c_str(), name.c_str(), _histoNBins, _histoMin, _histoMax);
        hV.push_back(h);

        name = syst + "_numu_universe_" + std::to_string(i);
        TH1D* h_numu = tfdir.make<TH1D>(name.c_str(), name.c_str(), _histoNBins, _histoMin, _histoMax);
        hV_numu.push_back(h_numu);

        name = syst + "_numubar_universe_" + std::to_string(i);
        TH1D* h_numubar = tfdir.make<TH1D>(name.c_str(), name.c_str(), _histoNBins, _histoMin, _histoMax);
        hV_numubar.push_back(h_numubar);

        name = syst + "_nue_universe_" + std::to_string(i);
        TH1D* h_nue = tfdir.make<TH1D>(name.c_str(), name.c_str(), _histoNBins, _histoMin, _histoMax);
        hV_nue.push_back(h_nue);

        name = syst + "_nuebar_universe_" + std::to_string(i);
        TH1D* h_nuebar = tfdir.make<TH1D>(name.c_str(), name.c_str(), _histoNBins, _histoMin, _histoMax);
        hV_nuebar.push_back(h_nuebar);
  
      }
      _histo_map[syst] = hV;
      _histo_map_numu[syst] = hV_numu;
      _histo_map_numubar[syst] = hV_numubar;
      _histo_map_nue[syst] = hV_nue;
      _histo_map_nuebar[syst] = hV_nuebar;
    }

  }

  
}


DEFINE_ART_MODULE(flux::FluxUncertaintiesAna)