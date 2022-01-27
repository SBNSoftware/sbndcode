////////////////////////////////////////////////////////////////////////
// Class:       PrismAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        PrismAnalyzer_module.cc
//
// Generated at Tue Dec  1 11:40:17 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <functional>

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

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TString.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/CoreUtils/zip.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

class PrismAnalyzer;


class PrismAnalyzer : public art::EDAnalyzer {
public:
  explicit PrismAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PrismAnalyzer(PrismAnalyzer const&) = delete;
  PrismAnalyzer(PrismAnalyzer&&) = delete;
  PrismAnalyzer& operator=(PrismAnalyzer const&) = delete;
  PrismAnalyzer& operator=(PrismAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(art::SubRun const& sr) override;
  virtual void endSubRun(art::SubRun const& sr) override;
  void endJob();

private:
  /// Returns the off-axis angle in radians
  float GetOffAxisAngle(float x, float y, float z);

  /// Returns the reco energy using QE formula
  float GetEnergyQE(float muon_energy, float muon_px, float muon_py, float muon_pz);

  /// Returns the reco energy using calorometric approach
  float GetEnergyCalo(float & e_excit, float & k_recoil, float & k_neutron, bool argoneut = false);

  /// Returns the T recoil
  float GetRecoilEnergy(int n_protons, TVector3 muon_t_mom, std::vector<TVector3> proton_t_mom);

  /// Returns the separation + exictation energy
  float GetSeparationEnergy(int n_protons);

  std::string _mctruth_producer = "generator";
  std::string _flux_eventweight_multisim_producer;
  float _beam_origin_x = 73.78;
  float _beam_origin_y = 0.0;
  float _beam_origin_z = 11000.0;
  // std::vector<float> _beam_origin; // = {-0.457 * 100, 0 * 100, 110 * 100}; // cm

  float _p_th;
  float _pi_th;

  const TDatabasePDG *_pdg_db = TDatabasePDG::Instance();

  TTree* _tree;

  int _run, _subrun, _event;

  float _nu_e; ///< Neutrino energy
  int _nu_pdg; ///< Neutrino PDG code
  int _nu_ccnc; ///< 0: CC, 1: NC
  int _nu_mode; ///< Neutrino interaction mode
  int _nu_int_type; ///< Neutrino interaction type
  float _nu_vtx_x; ///< Neutrino vertex X
  float _nu_vtx_y; ///< Neutrino vertex Y
  float _nu_vtx_z; ///< Neutrino vertex Z
  float _nu_px; ///< Neutrino momentum along X
  float _nu_py; ///< Neutrino momentum along Y
  float _nu_pz; ///< Neutrino momentum along Z
  float _nu_oaa; ///< Off-Axis Angle (angle w.r.t. beam)
  float _nu_roaa; ///< Real Off-Axis Angle (angle between neutrino parent p and neutrino p)
  float _nu_l; ///< L [cm]: Distance the neutrino travelled from production to interaction point
  int _nu_decay; ///< Neutrino parent decay type (see dkproc_t in dk2nu)
  float _nu_lepton_e; ///< Final state lepton energy
  float _nu_lepton_px; ///< Final state lepton px
  float _nu_lepton_py; ///< Final state lepton py
  float _nu_lepton_pz; ///< Final state lepton pz
  float _nu_e_reco_kin; ///< Neutrino reconstructed energy using QE formula
  float _nu_e_reco_calo; ///< Neutrino reconstructed energy using calorimetry
  float _nu_e_reco_caloa; ///< Neutrino reconstructed energy using calorimetry (argoneut)
  float _nu_e_excit; ///< Excitation + separation energy
  float _nu_e_recoil; ///< Nucleus kinetic recoil energy
  float _nu_e_neutron; ///< Neutron kinetic energy
  std::vector<int> _pars_pdg; ///< All other particles produced - pdg code
  std::vector<float> _pars_e; ///< All other particles produced - energy
  std::vector<float> _pars_px; ///< All other particles produced - p_x
  std::vector<float> _pars_py; ///< All other particles produced - p_y
  std::vector<float> _pars_pz; ///< All other particles produced - p_z

  float _nu_gt_FShadSystP4;
  float _nu_gt_FSleptonP4;
  float _nu_gt_TgtP4;

  float _nu_prod_vtx_x; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_y; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_z; ///< Neutrino production vertex in detector coordinates
  float _nu_prod_vtx_x_beam; ///< Neutrino production vertex in beamline coordinates
  float _nu_prod_vtx_y_beam; ///< Neutrino production vertex in beamline coordinates
  float _nu_prod_vtx_z_beam; ///< Neutrino production vertex in beamline coordinates

  int _p_type; ///< Neutrino parent PDG code
  float _p_dpx; ///< Neutrino parent px at neutrino production vertex
  float _p_dpy; ///< Neutrino parent py at neutrino production vertex
  float _p_dpz; ///< Neutrino parent pz at neutrino production vertex

  int _nu_pip_mult; ///< Pi plus and minum multiplicity
  int _nu_pi0_mult; ///< Pi zero multiplicity
  int _nu_p_mult; ///< Proton multiplicity

  int _evtwgt_flux_nfunc; ///< Number of functions used for FLUX reweighting (multisim)
  std::vector<std::string> _evtwgt_flux_funcname; ///< Names of the functions used for FLUX reweighting (multisim)
  std::vector<int> _evtwgt_flux_nweight; ///< Number of weights per function name used for FLUX reweighting (multisim)
  std::vector<std::vector<float>> _evtwgt_flux_weight; ///< Weights per function name used for FLUX reweighting (multisim)
  std::vector<float> _evtwgt_flux_oneweight; ///< Weights for FLUX reweighting (multisim) (combines all variations)

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun

};


PrismAnalyzer::PrismAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _flux_eventweight_multisim_producer = p.get<std::string>("FluxEventWeightProducer", "fluxweight");

  _beam_origin_x = p.get<float>("BeamCenterX");
  _beam_origin_y = p.get<float>("BeamCenterY");
  _beam_origin_z = p.get<float>("BeamCenterZ");

  _p_th = p.get<float>("ProtonThreshold", 35e-3); // Gev
  _pi_th = p.get<float>("PionThreshold", 20e-3); // GeV

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree","");

  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");

  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
  _tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
  _tree->Branch("nu_int_type", &_nu_int_type, "nu_int_type/I");
  _tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
  _tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
  _tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
  _tree->Branch("nu_px", &_nu_px, "nu_px/F");
  _tree->Branch("nu_py", &_nu_py, "nu_px/F");
  _tree->Branch("nu_pz", &_nu_pz, "nu_px/F");
  _tree->Branch("nu_oaa", &_nu_oaa, "nu_oaa/F");
  _tree->Branch("nu_roaa", &_nu_roaa, "nu_roaa/F");
  _tree->Branch("nu_l", &_nu_l, "nu_l/F");
  _tree->Branch("nu_decay", &_nu_decay, "nu_decay/I");
  _tree->Branch("nu_lepton_e", &_nu_lepton_e, "nu_lepton_e/F");
  _tree->Branch("nu_lepton_px", &_nu_lepton_px, "nu_lepton_px/F");
  _tree->Branch("nu_lepton_py", &_nu_lepton_py, "nu_lepton_py/F");
  _tree->Branch("nu_lepton_pz", &_nu_lepton_pz, "nu_lepton_pz/F");
  _tree->Branch("nu_e_reco_kin", &_nu_e_reco_kin, "nu_e_reco_kin/F");
  _tree->Branch("nu_e_reco_calo", &_nu_e_reco_calo, "nu_e_reco_calo/F");
  _tree->Branch("nu_e_reco_caloa", &_nu_e_reco_caloa, "nu_e_reco_caloa/F");
  _tree->Branch("nu_e_excit", &_nu_e_excit, "nu_e_excit/F");
  _tree->Branch("nu_e_recoil", &_nu_e_recoil, "nu_e_recoil/F");
  _tree->Branch("nu_e_neutron", &_nu_e_neutron, "nu_e_neutron/F");

  _tree->Branch("nu_gt_FShadSystP4", &_nu_gt_FShadSystP4, "nu_gt_FShadSystP4/F");
  _tree->Branch("nu_gt_FSleptonP4", &_nu_gt_FSleptonP4, "nu_gt_FSleptonP4/F");
  _tree->Branch("nu_gt_TgtP4", &_nu_gt_TgtP4, "nu_gt_TgtP4/F");

  _tree->Branch("nu_prod_vtx_x", &_nu_prod_vtx_x, "nu_prod_vtx_x/F");
  _tree->Branch("nu_prod_vtx_y", &_nu_prod_vtx_y, "nu_prod_vtx_y/F");
  _tree->Branch("nu_prod_vtx_z", &_nu_prod_vtx_z, "nu_prod_vtx_z/F");
  _tree->Branch("nu_prod_vtx_x_beam", &_nu_prod_vtx_x_beam, "nu_prod_vtx_x_beam/F");
  _tree->Branch("nu_prod_vtx_y_beam", &_nu_prod_vtx_y_beam, "nu_prod_vtx_y_beam/F");
  _tree->Branch("nu_prod_vtx_z_beam", &_nu_prod_vtx_z_beam, "nu_prod_vtx_z_beam/F");

  _tree->Branch("pars_pdg", "std::vector<int>", &_pars_pdg);
  _tree->Branch("pars_e", "std::vector<float>", &_pars_e);

  _tree->Branch("p_type", &_p_type, "p_type/I");
  _tree->Branch("p_dpx", &_p_dpx, "p_dpx/F");
  _tree->Branch("p_dpy", &_p_dpy, "p_dpy/F");
  _tree->Branch("p_dpz", &_p_dpz, "p_dpz/F");

  _tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
  _tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
  _tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");

  // _tree->Branch("evtwgt_flux_nfunc", &_evtwgt_flux_nfunc, "evtwgt_flux_nfunc/I");
  // _tree->Branch("evtwgt_flux_funcname", "std::vector<std::string>", &_evtwgt_flux_funcname);
  // _tree->Branch("evtwgt_flux_nweight", "std::vector<int>", &_evtwgt_flux_nweight);
  // _tree->Branch("evtwgt_flux_weight", "std::vector<std::vector<float>>", &_evtwgt_flux_weight);
  _tree->Branch("evtwgt_flux_oneweight", "std::vector<float>", &_evtwgt_flux_oneweight);

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void PrismAnalyzer::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  //
  // Get the associated MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_producer, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product " << _mctruth_producer << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

  //
  // Get the associated MCFlux and Flux Event Weight Map (ewm)
  //
  art::FindManyP<simb::MCFlux> mct_to_mcf (mct_h, e, _mctruth_producer);
  art::FindManyP<simb::GTruth> mct_to_gt (mct_h, e, _mctruth_producer);
  art::FindManyP<sbn::evwgh::EventWeightMap> mct_to_fluxewm(mct_h, e, _flux_eventweight_multisim_producer);

  _pars_pdg.clear();
  _pars_e.clear();
  _nu_pi0_mult = _nu_pip_mult = _nu_p_mult = 0;

  //
  // Loop over the neutrino interactions in this event
  //
  for (size_t i = 0; i < mct_v.size(); i++) {
    if (mct_v.at(i)->Origin() != simb::kBeamNeutrino) {
      std::cout << "[PrismAnalyzer] MCTruth from generator does not have neutrino origin?!" << std::endl;
    }

    // Get the associated GTruth
    std::vector<art::Ptr<simb::GTruth>> gt_v = mct_to_gt.at(i);
    assert(gt_v.size() == 1);
    auto gt = gt_v[0];

    _nu_e = mct_v[i]->GetNeutrino().Nu().E();
    _nu_pdg = mct_v[i]->GetNeutrino().Nu().PdgCode();
    _nu_ccnc = mct_v[i]->GetNeutrino().CCNC();
    _nu_mode = mct_v[i]->GetNeutrino().Mode();
    _nu_int_type = mct_v[i]->GetNeutrino().InteractionType();
    _nu_vtx_x = mct_v[i]->GetNeutrino().Nu().Vx();
    _nu_vtx_y = mct_v[i]->GetNeutrino().Nu().Vy();
    _nu_vtx_z = mct_v[i]->GetNeutrino().Nu().Vz();
    _nu_px = mct_v[i]->GetNeutrino().Nu().Px();
    _nu_py = mct_v[i]->GetNeutrino().Nu().Py();
    _nu_pz = mct_v[i]->GetNeutrino().Nu().Pz();

    _nu_oaa = GetOffAxisAngle(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z);

    _nu_gt_FShadSystP4 = gt->fFShadSystP4.E();
    _nu_gt_FSleptonP4 = gt->fFSleptonP4.E();
    _nu_gt_TgtP4 = gt->fTgtP4.E();

    // Get the associated MCFlux
    std::vector<art::Ptr<simb::MCFlux>> mcf_v = mct_to_mcf.at(i);
    assert(mcf_v.size() == 1);
    auto mcf = mcf_v[0];
    assert(mcf->fntype == _nu_pdg);

    // Save the neutrino production vertex, this is in beam coordinates
    _nu_prod_vtx_x_beam = mcf->fvx;
    _nu_prod_vtx_y_beam = mcf->fvy;
    _nu_prod_vtx_z_beam = mcf->fvz;

    // Also convert it to detector coordinates
    _nu_prod_vtx_x = mcf->fvx - _beam_origin_x;
    _nu_prod_vtx_y = mcf->fvy - _beam_origin_y;
    _nu_prod_vtx_z = mcf->fvz - _beam_origin_z;

    auto diff_x = _nu_vtx_x - _nu_prod_vtx_x;
    auto diff_y = _nu_vtx_y - _nu_prod_vtx_y;
    auto diff_z = _nu_vtx_z - _nu_prod_vtx_z;
    _nu_l = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);

    _nu_decay = mcf->fndecay;
    _p_type = mcf->fptype;
    _p_dpx = mcf->fpdpx;
    _p_dpy = mcf->fpdpy;
    _p_dpz = mcf->fpdpz;
    _nu_roaa = TVector3(_nu_px, _nu_py, _nu_pz).Angle(TVector3(_p_dpx, _p_dpy, _p_dpz));

    _nu_lepton_e = -9999;
    _nu_lepton_px = -9999;
    _nu_lepton_py = -9999;
    _nu_lepton_pz = -9999;
    _nu_e_reco_kin = -9999;
    _nu_e_reco_calo = -9999;
    _nu_e_reco_caloa = -9999;

    for (int p = 0; p < mct_v[i]->NParticles(); p++) {
      auto const & mcp = mct_v[i]->GetParticle(p);

      if (mcp.StatusCode() != 1) continue;

      _pars_pdg.push_back(mcp.PdgCode());
      _pars_e.push_back(mcp.E());
      _pars_px.push_back(mcp.Px());
      _pars_py.push_back(mcp.Py());
      _pars_pz.push_back(mcp.Pz());

      if (mcp.PdgCode() == 111) {
        _nu_pi0_mult++;
      } else if (std::abs(mcp.PdgCode()) == 211) {
        _nu_pip_mult++;
      }
      else if (std::abs(mcp.PdgCode()) == 2112) {
        _nu_p_mult++;
      } else if (std::abs(mcp.PdgCode()) == 13 || std::abs(mcp.PdgCode()) == 11) {
        _nu_lepton_e = mcp.E();
        _nu_lepton_px = mcp.Px();
        _nu_lepton_py = mcp.Py();
        _nu_lepton_pz = mcp.Pz();
        _nu_e_reco_kin = GetEnergyQE(mcp.E(), mcp.Px(), mcp.Py(), mcp.Pz());
      }

    }

    std::cout << "******************************************************" << std::endl;
    _nu_e_reco_calo = GetEnergyCalo(_nu_e_excit, _nu_e_recoil, _nu_e_neutron, false);
    std::cout << "******************************************************" << std::endl;
    _nu_e_reco_caloa = GetEnergyCalo(_nu_e_excit, _nu_e_recoil, _nu_e_neutron, true);

    std::cout << "True energy: " << _nu_e << " Reco calo: " << _nu_e_reco_calo << std::endl;
    std::cout << "True energy: " << _nu_e << " Reco calo: " << _nu_e_reco_caloa << " (argoneut)" << std::endl;
    std::cout << "******************************************************" << std::endl;
    // _n_pi0 = 0;
    // for (int p = 0; p < mct_v[i]->NParticles(); p++) {
    //   auto const & mcp = mct_v[i]->GetParticle(p);
    //   if (mcp.StatusCode() != 1) continue;
    //   if (mcp.PdgCode() != 111) continue;
    //   _n_pi0++;
    // }
    // for (int p = 0; p < mct_v[i]->NParticles(); p++) {
    //   auto const & mcp = mct_v[i]->GetParticle(p);
    //   if (mcp.StatusCode() != 1) continue;
    //   if (mcp.PdgCode() != 2112) continue;
    //   std::cout << "Got neutron, mother is " << mcp.Mother() << ", energy is " << mcp.E() << std::endl;
    // }

    // Flux weights
    std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> flux_ewm_v = mct_to_fluxewm.at(i);
    if (flux_ewm_v.size() != 1) {
      std::cout << "[PrismAnalyzer] EventWeightMap of " << _flux_eventweight_multisim_producer << " bigger than 1?" << std::endl;
    }
    std::map<std::string, std::vector<float>> evtwgt_map = *(flux_ewm_v[0]);

    _evtwgt_flux_funcname.clear();
    _evtwgt_flux_weight.clear();
    _evtwgt_flux_nweight.clear();
    _evtwgt_flux_oneweight.clear();

    std::vector<float> previous_weights;
    std::vector<float> final_weights;

    int countFunc = 0;
    for(auto it : evtwgt_map) {
      std::string func_name = it.first;
      std::vector<float> weight_v = it.second;

      if (previous_weights.size() == 0) {
        previous_weights.resize(weight_v.size(), 1.);
        final_weights.resize(weight_v.size(), 1.);
      }

      _evtwgt_flux_funcname.push_back(func_name);
      _evtwgt_flux_weight.push_back(weight_v);
      _evtwgt_flux_nweight.push_back(weight_v.size());
      countFunc++;

      std::cout << "weight_v.at(0) " << weight_v.at(0) << ", final_weights.at(0) " << final_weights.at(0) << std::endl;

      // Construct a single weight
      std::transform(previous_weights.begin(), previous_weights.end(),
                     weight_v.begin(),
                     final_weights.begin(),
                     std::multiplies<float>());
      previous_weights = final_weights;
    }
    _evtwgt_flux_nfunc = countFunc;
    _evtwgt_flux_oneweight = final_weights;

    _tree->Fill();
  }

  // art::Handle<std::vector<simb::MCParticle>> mcp_h;
  // e.getByLabel("largeant", mcp_h);
  // if(!mcp_h.isValid()){
  //   std::cout << "MCParticle product not found..." << std::endl;
  //   throw std::exception();
  // }
  // std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  // art::fill_ptr_vector(mcp_v, mcp_h);
  // for (size_t p = 0; p < mcp_v.size(); p++) {
  //   auto mcp = mcp_v[p];
  //   if (mcp->StatusCode() != 1) continue;
  //   if (mcp->Mother() != 1) continue;
  //   std::cout << "MCParticle " << p
  //             << ", trackid " << mcp->TrackId()
  //             << ", pdg " << mcp->PdgCode()
  //             << ", E " << mcp->E()
  //             << ", mother " << mcp->Mother()
  //             << ", process " << mcp->Process() << std::endl;
  // }
}

float PrismAnalyzer::GetOffAxisAngle(float x, float y, float z) {
  TVector3 nu_vtx(x, y, z);
  TVector3 beam_center(_beam_origin_x, _beam_origin_y, _beam_origin_z);
  TVector3 nu_vtx_beam = nu_vtx + beam_center;
  TVector3 beam(0, 0, 1);

  return beam.Angle(nu_vtx_beam);
}

float PrismAnalyzer::GetEnergyQE(float muon_energy, float px, float py, float pz) {
  float p_mass = _pdg_db->GetParticle(2212)->Mass(); // proton mass
  float n_mass = _pdg_db->GetParticle(2112)->Mass(); // neutron mass
  float m_mass = _pdg_db->GetParticle(13)->Mass(); // muon mass
  float Eb = 0.030; // binding energy

  float p = std::sqrt(px * px + py * py + pz * pz);
  // float pxy = TMath::Sqrt( px*px + py*py )
  float theta = TMath::ACos(pz / p);

  float num = p_mass * p_mass;
  num -= m_mass * m_mass;
  num -= (n_mass - Eb) * (n_mass - Eb);
  num += 2 * (n_mass - Eb) * muon_energy;

  float den = 2. * (n_mass - Eb - muon_energy + p * std::cos(theta));

  return num / den;
}

float PrismAnalyzer::GetEnergyCalo(float & e_excit, float & k_recoil, float & k_neutron, bool argoneut) {

  float calo_e = 0;
  int n_protons = 0;
  std::vector<TVector3> proton_t_mom = std::vector<TVector3>(); // Transverse momenta of protons
  TVector3 muon_t_mom = TVector3(0, 0, 0); // Transverse momentum of muon
  e_excit = k_recoil = k_neutron = 0;

  // for (auto&& [pdg, e]: util::zip(pdgs, energies)) {
  for (size_t i = 0; i < _pars_pdg.size(); i++) {
    int pdg = _pars_pdg[i];
    float e = _pars_e[i];
    float px = _pars_px[i];
    float py = _pars_py[i];

    auto particle = _pdg_db->GetParticle(pdg);

    if (!particle) {
      continue;
    }

    float mass = particle->Mass(); // GeV
    float charge = particle->Charge() / 3;
    std::cout << std::endl;
    std::cout << "PDG: " << pdg << " - Mass: " << mass << " = Charge: " << charge << std::endl;

    if (charge == 0) {
      std::cout << "Charge is zero." << std::endl;
    }

    float k = e - mass;
    // float p_t = std::sqrt(px * px + py * py);
    TVector3 p_t = TVector3(px, py, 0.);

    TString p_class = TString(particle->ParticleClass());
    std::cout << "Class: " << p_class << std::endl;

    bool is_lepton = p_class.Contains("Lepton");
    bool is_meson = p_class.Contains("Meson");
    bool is_proton = p_class.Contains("Baryon") && (charge != 0);
    bool is_neutron = p_class.Contains("Baryon") && (charge == 0);

    if (argoneut)
    {
      if (is_meson and k > 8e-3) {
        calo_e += e;
        std::cout << "Meson, e = " << e << std::endl;
      }
      else if (is_proton and k > 21e-3) {
        calo_e += k;
        n_protons++;
        proton_t_mom.push_back(p_t);
        std::cout << "Proton, k = " << k << std::endl;
      }
      else if (is_neutron) {
        k_neutron += k;
      }
      else if (is_lepton) {
        calo_e += e;
        muon_t_mom = p_t;
        std::cout << "Lepton, e = " << e << std::endl;
      }
      else std::cout << "Unknown" << std::endl;
    } else
    {
      if (is_meson and k > _pi_th) calo_e += e;
      else if (is_proton and k > _p_th) calo_e += k;
      else if (is_lepton) calo_e += e;
      else std::cout << "Unknown" << std::endl;
    }
  }

  if (argoneut) {
    k_recoil = GetRecoilEnergy(n_protons, muon_t_mom, proton_t_mom);
    e_excit = GetSeparationEnergy(n_protons);

    calo_e += k_recoil;
    calo_e += e_excit;

    if (n_protons == 2) std::cout << "+++ 2 protons +++" << std::endl;
    else calo_e *= -1; // assign a negative value if not 2 protons, as the estimation is not correct
  }

  return calo_e;
}

float PrismAnalyzer::GetRecoilEnergy(int n_protons, TVector3 muon_t_mom, std::vector<TVector3> proton_t_mom) {

  // Estimate p_t miss
  TVector3 p_t_miss_v = muon_t_mom;
  for (auto t : proton_t_mom) { p_t_miss_v += t; };

  float p_t_miss = p_t_miss_v.Mag() * -1;

  std::cout << "Pt miss is " << p_t_miss << std::endl;


  // 40Ar mass: 39.9623831 amu * 931.5 MeV / amu = 37.22 GeV
  float Ar40_mass = 37.22; // GeV
  double Cl38_mass = 35.3669; // GeV
  float proton_mass = 0.938; // GeV
  float m_anp = Ar40_mass - n_protons * proton_mass;
  if (n_protons == 2) {
    m_anp = Cl38_mass;
  }

  std::cout << "M_{A-np} is " << m_anp << ", n protons: " << n_protons << std::endl;

  // Estimate recoil kinetic energy
  // float t_recoil = p_t_miss * p_t_miss / (2 * m_anp);
  float t_recoil = std::sqrt(p_t_miss * p_t_miss + m_anp * m_anp) - m_anp;

  std::cout << "T_{recoil} is " << t_recoil << std::endl;

  return t_recoil;

}

float PrismAnalyzer::GetSeparationEnergy(int n_protons) {

  // return 30.4*1e-3; // GeV

  if (n_protons == 2) {
    return 30.4*1e-3; // GeV
  }

  return 8.595e-3 * n_protons; // GeV
  // http://barwinski.net/isotopes/all_isotopes_query_results.php?element=18-Argon&order_isotopes_by_field=total_binding_energy/mass_number+DESC&all_radioactive_stable=ALL
}


void PrismAnalyzer::beginSubRun(art::SubRun const& sr) {

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(_mctruth_producer, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }
  std::cout << "POT for this subrun: " << _sr_pot << std::endl;

  _sr_tree->Fill();

}

void PrismAnalyzer::endJob() {

}

void PrismAnalyzer::endSubRun(art::SubRun const& sr) {

}

DEFINE_ART_MODULE(PrismAnalyzer)

















