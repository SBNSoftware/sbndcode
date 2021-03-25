////////////////////////////////////////////////////////////////////////
// Class:       PrismAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        PrismAnalyzer_module.cc
//
// Generated at Tue Dec  1 11:40:17 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
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

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"

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

  std::string _mctruth_producer = "generator";
  std::vector<float> _beam_origin = {-0.457 * 100, 0 * 100, 110 * 100}; // cm

  const TDatabasePDG *_pdg_db = TDatabasePDG::Instance();

  TTree* _tree;

  int _run, _subrun, _event;

  float _nu_e;
  int _nu_pdg;
  int _nu_ccnc;
  int _nu_mode;
  int _nu_int_type;
  float _nu_vtx_x;
  float _nu_vtx_y;
  float _nu_vtx_z;
  float _nu_px;
  float _nu_py;
  float _nu_pz;
  float _nu_oaa; ///< Off-Axis Angle (angle w.r.t. beam)
  float _nu_roaa; ///< Real Off-Axis Angle (angle between neutrino parent p and neutrino p)
  float _nu_l; ///< L [cm]: Distance the neutrino travelled from production to interaction point
  int _nu_decay; ///< Neutrino parent decay type (see dkproc_t in dk2nu)
  float _nu_mu_e; ///< Final state muon (if any) energy
  float _nu_mu_px; ///< Final state muon (if any) px
  float _nu_mu_py; ///< Final state muon (if any) py
  float _nu_mu_pz; ///< Final state muon (if any) pz
  float _nu_e_reco; ///< Neutrino reconstructe energy using QE formula
  std::vector<int> _pars_pdg; ///< All other particles produced - pdg code
  std::vector<float> _pars_e; ///< All other particles produced - energy

  float _nu_prod_vtx_x;
  float _nu_prod_vtx_y;
  float _nu_prod_vtx_z;
  float _nu_prod_vtx_x_beam;
  float _nu_prod_vtx_y_beam;
  float _nu_prod_vtx_z_beam;

  int _p_type; ///< Neutrino parent PDG code
  float _p_dpx; ///< Neutrino parent px at neutrino production vertex
  float _p_dpy; ///< Neutrino parent py at neutrino production vertex
  float _p_dpz; ///< Neutrino parent pz at neutrino production vertex

  int _nu_pip_mult; ///< Pi0 multiplicity
  int _nu_pi0_mult; ///< Pi plus multiplicity
  int _nu_p_mult; ///< Proton multiplicity

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

};


PrismAnalyzer::PrismAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
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
  _tree->Branch("nu_mu_e", &_nu_mu_e, "nu_mu_e/F");
  _tree->Branch("nu_mu_px", &_nu_mu_px, "nu_mu_px/F");
  _tree->Branch("nu_mu_py", &_nu_mu_py, "nu_mu_py/F");
  _tree->Branch("nu_mu_pz", &_nu_mu_pz, "nu_mu_pz/F");
  _tree->Branch("nu_e_reco", &_nu_e_reco, "nu_e_reco/F");

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
  // Get the associated MCFlux
  //
  art::FindManyP<simb::MCFlux> mct_to_mcf (mct_h, e, _mctruth_producer);

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
    _nu_prod_vtx_x = mcf->fvx - _beam_origin[0];
    _nu_prod_vtx_y = mcf->fvy - _beam_origin[1];
    _nu_prod_vtx_z = mcf->fvz - _beam_origin[2];

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

    _nu_mu_e = -9999;
    _nu_mu_px = -9999;
    _nu_mu_py = -9999;
    _nu_mu_pz = -9999;
    _nu_e_reco = -9999;

    for (int p = 0; p < mct_v[i]->NParticles(); p++) {
      auto const & mcp = mct_v[i]->GetParticle(p);

      if (mcp.StatusCode() != 1) continue;

      _pars_pdg.push_back(mcp.PdgCode());
      _pars_e.push_back(mcp.E());

      if (mcp.PdgCode() == 111) {
        _nu_pi0_mult++;
      } else if (std::abs(mcp.PdgCode()) == 211) {
        _nu_pip_mult++;
      }
      else if (std::abs(mcp.PdgCode()) == 2112) {
        _nu_p_mult++;
      } else if (std::abs(mcp.PdgCode()) == 13) {
        _nu_mu_e = mcp.E();
        _nu_mu_px = mcp.Px();
        _nu_mu_py = mcp.Py();
        _nu_mu_pz = mcp.Pz();
        _nu_e_reco = GetEnergyQE(mcp.E(), mcp.Px(), mcp.Py(), mcp.Pz());
      }
    }

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
  TVector3 beam_center(_beam_origin[0], _beam_origin[1], _beam_origin[2]);
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

















