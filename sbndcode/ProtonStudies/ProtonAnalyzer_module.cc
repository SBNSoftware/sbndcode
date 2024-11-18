////////////////////////////////////////////////////////////////////////
// Class:       ProtonAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        ProtonAnalyzer_module.cc
//
// Generated at Fri Nov  8 13:50:20 2024 by Marco Del Tutto using cetskelgen
// from cetlib version 3.18.02.
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
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/CoreUtils/zip.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "larcore/Geometry/Geometry.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

class ProtonAnalyzer;


class ProtonAnalyzer : public art::EDAnalyzer {
public:
  explicit ProtonAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtonAnalyzer(ProtonAnalyzer const&) = delete;
  ProtonAnalyzer(ProtonAnalyzer&&) = delete;
  ProtonAnalyzer& operator=(ProtonAnalyzer const&) = delete;
  ProtonAnalyzer& operator=(ProtonAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(art::SubRun const& sr) override;
  virtual void endSubRun(art::SubRun const& sr) override;
  void endJob();

private:

  /// Check if the point is in the detector
  bool InDetector(const double& x, const double& y, const double& z) const;

  /// Check if an MCP passes through the detector, sets step to the step index when particle is in detector
  bool InDetector(const art::Ptr<simb::MCParticle>, int & step);


  std::string _mctruth_producer = "generator";
  std::string _mcpar_producer = "largeant";
  std::string _sed_producer = "largeant:LArG4DetectorServicevolTPCActive";

  const TDatabasePDG *_pdg_db = TDatabasePDG::Instance();

  std::map<unsigned int, simb::MCParticle> _trackid_to_mcparticle;


  double _x_max; //!< x-max of volume box used to determine whether to save track information
  double _x_min; //!< x-min of volume box used to determine whether to save track information
  double _y_max; //!< y-max of volume box used to determine whether to save track information
  double _y_min; //!< y-min of volume box used to determine whether to save track information
  double _z_max; //!< z-max of volume box used to determine whether to save track information
  double _z_min; //!< z-min of volume box used to determine whether to save track information


  std::vector<std::string> _flash_label_v;
  std::vector<std::string> _ophit_label_v;
  // std::vector<std::string> _mcflash_label_v;
  std::vector<int> _tpc_v;
  std::vector<float> _flash_time_window;

  opdet::sbndPDMapAlg _pds_map;


  TTree* _protontree;
  float _nu_edep_from_all; ///< Energy deposit from all particles with neutrino origin
  float _nu_edep_from_p; ///< Energy deposit from protons with neutrino origin
  float _nu_edep_from_p_nin; ///< Energy deposit from protons coming from neutronInelastic with neutrino origin
  float _n_nu_np_av; ///< Number of protons in AV from neutrino origin
  float _n_nu_in_av; ///< Number of neutrino events in AV

  float _nu_in_av_e; ///< Neutrino energy
  int _nu_in_av_pdg; ///< Neutrino PDG code
  int _nu_in_av_ccnc; ///< 0: CC, 1: NC
  int _nu_in_av_mode; ///< Neutrino interaction mode
  int _nu_in_av_int_type; ///< Neutrino interaction type
  float _nu_in_av_vtx_x; ///< Neutrino vertex X
  float _nu_in_av_vtx_y; ///< Neutrino vertex Y
  float _nu_in_av_vtx_z; ///< Neutrino vertex Z

  float _proton_ke; ///< Highest-energy proton kinetic energy
  std::string _proton_process; ///< Highest-energy proton process
  float _proton_vtx_x; ///< Highest-energy proton vertex x
  float _proton_vtx_y; ///< Highest-energy proton vertex y
  float _proton_vtx_z; ///< Highest-energy proton vertex z
  float _proton_end_x; ///< Highest-energy proton end x
  float _proton_end_y; ///< Highest-energy proton end y
  float _proton_end_z; ///< Highest-energy proton end z

  double _flash_tpc0_time, _flash_tpc1_time;
  double _flash_tpc0_total_pe, _flash_tpc1_total_pe;
  std::vector<double> _flash_tpc0_pe_v, _flash_tpc1_pe_v;
  double _flash_tpc0_y, _flash_tpc1_y;
  double _flash_tpc0_yerr, _flash_tpc1_yerr;
  double _flash_tpc0_z, _flash_tpc1_z;
  double _flash_tpc0_zerr, _flash_tpc1_zerr;


  TTree* _tree;

  int _run, _subrun, _event;

  float _nu_e; ///< Neutrino energy
  int _nu_pdg; ///< Neutrino PDG code
  int _nu_ccnc; ///< 0: CC, 1: NC
  int _nu_mode; ///< Neutrino interaction mode
  int _nu_int_type; ///< Neutrino interaction type
  float _nu_x; ///< Neutrino Bjorken x
  float _nu_y; ///< Neutrino inelasticity y
  float _nu_w; ///< Neutrino hadronic invariant mass w
  float _nu_q2; ///< Neutrino q2
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

  int _nu_pip_mult; ///< Pi plus and minus multiplicity
  int _nu_pi0_mult; ///< Pi zero multiplicity
  int _nu_p_mult; ///< Proton multiplicity

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun

};


ProtonAnalyzer::ProtonAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  // _beam_origin_x = p.get<float>("BeamCenterX");

  _flash_label_v     = p.get<std::vector<std::string>>("FlashProducer",   _flash_label_v);
  // _mcflash_label_v = p.get<std::vector<std::string>>("MCFlashProducer", _mcflash_label_v);
  _ophit_label_v     = p.get<std::vector<std::string>>("OpHitProducer", _ophit_label_v);
  _tpc_v             = p.get<std::vector<int>>("TPC", _tpc_v);
  _flash_time_window = p.get<std::vector<float>>("FlashTimeWindow", _flash_time_window);


  art::ServiceHandle<art::TFileService> fs;

  //
  // Proton TTree
  //
  _protontree = fs->make<TTree>("protontree","");

  _protontree->Branch("run", &_run, "run/I");
  _protontree->Branch("subrun", &_subrun, "subrun/I");
  _protontree->Branch("event", &_event, "event/I");

  _protontree->Branch("nu_edep_from_all", &_nu_edep_from_all, "nu_edep_from_all/F");
  _protontree->Branch("nu_edep_from_p", &_nu_edep_from_p, "nu_edep_from_p/F");
  _protontree->Branch("nu_edep_from_p_nin", &_nu_edep_from_p_nin, "nu_edep_from_p_nin/F");
  _protontree->Branch("n_nu_np_av", &_n_nu_np_av, "n_nu_np_av/F");
  _protontree->Branch("n_nu_in_av", &_n_nu_in_av, "n_nu_in_av/F");

  _protontree->Branch("nu_in_av_e", &_nu_in_av_e, "nu_in_av_e/F");
  _protontree->Branch("nu_in_av_pdg", &_nu_in_av_pdg, "nu_in_av_pdg/I");
  _protontree->Branch("nu_in_av_ccnc", &_nu_in_av_ccnc, "nu_in_av_ccnc/I");
  _protontree->Branch("nu_in_av_mode", &_nu_in_av_mode, "nu_in_av_mode/I");
  _protontree->Branch("nu_in_av_int_type", &_nu_in_av_int_type, "nu_in_av_int_type/I");
  _protontree->Branch("nu_in_av_vtx_x", &_nu_in_av_vtx_x, "nu_in_av_vtx_x/F");
  _protontree->Branch("nu_in_av_vtx_y", &_nu_in_av_vtx_y, "nu_in_av_vtx_y/F");
  _protontree->Branch("nu_in_av_vtx_z", &_nu_in_av_vtx_z, "nu_in_av_vtx_z/F");

  _protontree->Branch("proton_ke", &_proton_ke, "proton_ke/F");
  _protontree->Branch("proton_process", "std::string", &_proton_process);
  _protontree->Branch("proton_vtx_x", &_proton_vtx_x, "proton_vtx_x/F");
  _protontree->Branch("proton_vtx_y", &_proton_vtx_y, "proton_vtx_y/F");
  _protontree->Branch("proton_vtx_z", &_proton_vtx_z, "proton_vtx_z/F");
  _protontree->Branch("proton_end_x", &_proton_end_x, "proton_end_x/F");
  _protontree->Branch("proton_end_y", &_proton_end_y, "proton_end_y/F");
  _protontree->Branch("proton_end_z", &_proton_end_z, "proton_end_z/F");

  _protontree->Branch("flash_tpc0_time", &_flash_tpc0_time, "flash_tpc0_time/D");
  _protontree->Branch("flash_tpc0_total_pe", &_flash_tpc0_total_pe, "flash_tpc0_total_pe/D");
  _protontree->Branch("flash_tpc0_y", &_flash_tpc0_y, "flash_tpc0_y/D");
  _protontree->Branch("flash_tpc0_yerr", &_flash_tpc0_yerr, "flash_tpc0_yerr/D");
  _protontree->Branch("flash_tpc0_z", &_flash_tpc0_z, "flash_tpc0_z/D");
  _protontree->Branch("flash_tpc0_zerr", &_flash_tpc0_zerr, "flash_tpc0_zerr/D");
  _protontree->Branch("flash_tpc0_pe_v", "std::vector<double>", &_flash_tpc0_pe_v);
  _protontree->Branch("flash_tpc1_time", &_flash_tpc1_time, "flash_tpc1_time/D");
  _protontree->Branch("flash_tpc1_total_pe", &_flash_tpc1_total_pe, "flash_tpc1_total_pe/D");
  _protontree->Branch("flash_tpc1_y", &_flash_tpc1_y, "flash_tpc1_y/D");
  _protontree->Branch("flash_tpc1_yerr", &_flash_tpc1_yerr, "flash_tpc1_yerr/D");
  _protontree->Branch("flash_tpc1_z", &_flash_tpc1_z, "flash_tpc1_z/D");
  _protontree->Branch("flash_tpc1_zerr", &_flash_tpc1_zerr, "flash_tpc1_zerr/D");
  _protontree->Branch("flash_tpc1_pe_v", "std::vector<double>", &_flash_tpc1_pe_v);


  //
  // Neutrino TTree
  //
  _tree = fs->make<TTree>("tree","");

  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");

  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
  _tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
  _tree->Branch("nu_int_type", &_nu_int_type, "nu_int_type/I");
  _tree->Branch("nu_x", &_nu_x, "nu_x/D");
  _tree->Branch("nu_y", &_nu_y, "nu_x/D");
  _tree->Branch("nu_w", &_nu_w, "nu_x/D");
  _tree->Branch("nu_q2", &_nu_q2, "nu_x/D");
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


  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");





  // Iterate over all TPC's to get bounding box that covers volumes of each individual TPC in the detector
  // art::ServiceHandle<geo::Geometry const> geo;
  _x_min = -200;
  _y_min = -200;
  _z_min = 0;
  _x_max = 200;
  _y_max = 200;
  _z_max = 500;


}

void ProtonAnalyzer::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();
  std::cout << "RUN: " << _run << ", SUBRUN: " << _subrun << ", EVENT:" << _event << std::endl;

  TParticlePDG * p_pdg_par = _pdg_db->GetParticle(2212);

  


  //
  // Get the MCTracks
  //
  art::Handle<std::vector<sim::MCTrack> > mc_track_h;
  e.getByLabel("mcreco", mc_track_h);
  if(!mc_track_h.isValid()){
    std::cout << "MCTrack product " << "mcreco" << " not found..." << std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<sim::MCTrack>> mc_track_v;
  art::fill_ptr_vector(mc_track_v, mc_track_h);

  for (size_t i = 0; i < mc_track_v.size(); i++) {
    auto mc_track = mc_track_v.at(i);

    std::cout << "MCTrack " << i << ": ancestor pdg " << mc_track->AncestorPdgCode()
                                  << ", trackid " << mc_track->TrackID()
                                  << ", ancestor process " << mc_track->AncestorProcess()
                                  << ", mother pdg " << mc_track->MotherPdgCode()
                                  << ", mother process " << mc_track->MotherProcess()
                                  << ", PDG " << mc_track->PdgCode()
                                  << ", process " << mc_track->Process()
                                  << ", origin " << mc_track->Origin()
                                  << ", size " << mc_track->size()
                                  << std::endl;
  }


  //
  // Get the MCParticles
  //
  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  e.getByLabel(_mcpar_producer, mcp_h);
  if(!mcp_h.isValid()){
    std::cout << "MCParticle product not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcp_h);
  std::cout << "There are N MCParticle = " << mcp_v.size() << std::endl;

  _trackid_to_mcparticle.clear();

  for (size_t p = 0; p < mcp_v.size(); p++) {
    auto mcp = mcp_v[p];

    _trackid_to_mcparticle[mcp->TrackId()] = *mcp;

    if (mcp->StatusCode() != 1) continue;
    // if (mcp->Mother() != 1) continue;


    TParticlePDG * pdg_par = _pdg_db->GetParticle(mcp->PdgCode());

    std::cout << "MCParticle " << p
              << ", trackid " << mcp->TrackId()
              << ", pdg " << mcp->PdgCode()
              << ", E " << mcp->E()
              << ", Vx " << mcp->Vx()
              << ", Vy " << mcp->Vy()
              << ", Vz " << mcp->Vz()
              << ", mother " << mcp->Mother()
              << ", process " << mcp->Process()
              << ", end process " << mcp->EndProcess()
              << ", n traj " << mcp->NumberTrajectoryPoints();

    if (pdg_par) {
      std::cout << ", KE " << mcp->E() - pdg_par->Mass()
                << ", charge " << pdg_par->Charge() / 3;
    }
    std::cout << std::endl;
  }


  _nu_edep_from_all = 0;
  _nu_edep_from_p = 0;
  _nu_edep_from_p_nin = 0;
  _n_nu_in_av = 0;
  _n_nu_np_av = 0;
  _proton_ke = 0;
  _proton_vtx_x = 0;
  _proton_vtx_y = 0;
  _proton_vtx_z = 0;


  //
  // Get the SimEnergyDeposit
  //
  const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>&
      seds(e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(_sed_producer));

  // Loop over all the SimEnergyDeposit in the event
  for (const sim::SimEnergyDeposit& sed : *seds)
  {

    // Tracks with 20000000 are cosmics, with 10000000 are from neutrino interactions
    if (sed.TrackID() <= -20000000 or sed.TrackID() >= 20000000) continue;

    // Negative track IDs belong to particles dropped in the simulation
    if (sed.TrackID() < 0) continue;

    // Find the MCParticle associated with this energy deposit
    auto iter = _trackid_to_mcparticle.find(sed.TrackID());
    if (iter == _trackid_to_mcparticle.end()) {
      continue;
    }
    auto mcp = iter->second;

    std::cout << "SED - track ID: " << sed.TrackID()
              << ", PDG " << sed.PdgCode()
              << ", Energy " << sed.Energy() 
              << ", Electrons " << sed.NumElectrons()
              << ", X " << sed.MidPointX()
              << ", Z " << sed.MidPointZ()
              << ", Process " << mcp.Process() << " " << mcp.EndProcess() << std::endl;
 

    // Accumulate the energy deposited
    _nu_edep_from_all += sed.Energy();

    // Select energy deposits from protons only
    if(std::abs(sed.PdgCode()) == 2212) {
      _nu_edep_from_p += sed.Energy();

      if (mcp.Process() == "neutronInelastic") {
        _nu_edep_from_p_nin += sed.Energy();
      }

      // Save the highest energy proton
      float ke = mcp.E() - p_pdg_par->Mass();
      if (ke > _proton_ke) {
        _proton_ke = ke;
        _proton_process = mcp.Process();
        _proton_vtx_x = mcp.Vx();
        _proton_vtx_y = mcp.Vy();
        _proton_vtx_z = mcp.Vz();
        _proton_end_x = mcp.EndX();
        _proton_end_y = mcp.EndY();
        _proton_end_z = mcp.EndZ();
      }
    }

    // Check particle time is within the beam time
    // const double time(sed.Time() * 1e-3); // [ns] -> [us]

    sed.Time();

  }
  std::cout << "SED from nu all [MeV]: " << _nu_edep_from_all << std::endl;
  std::cout << "SED from nu p [MeV]: " << _nu_edep_from_p << std::endl;
  std::cout << "SED from nu p from neutronInelastic [MeV]: " << _nu_edep_from_p_nin << std::endl;


  // _n_nu_np_av?
  // _n_nu_in_av?


  _nu_in_av_e = -9999;
  _nu_in_av_pdg = -9999;
  _nu_in_av_ccnc = -9999;
  _nu_in_av_mode = -9999;
  _nu_in_av_int_type = -9999;
  _nu_in_av_vtx_x = -9999;
  _nu_in_av_vtx_y = -9999;
  _nu_in_av_vtx_z = -9999;


  //
  // Get the MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_producer, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product " << _mctruth_producer << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

  // Get the associated MCFlux and GTruth
  art::FindManyP<simb::MCFlux> mct_to_mcf (mct_h, e, _mctruth_producer);
  art::FindManyP<simb::GTruth> mct_to_gt (mct_h, e, _mctruth_producer);

  _pars_pdg.clear();
  _pars_e.clear();
  _nu_pi0_mult = _nu_pip_mult = _nu_p_mult = 0;

  //
  // Loop over the neutrino interactions in this event
  //
  for (size_t i = 0; i < mct_v.size(); i++) {
    if (mct_v.at(i)->Origin() != simb::kBeamNeutrino) {
      std::cout << "[ProtonAnalyzer] MCTruth from generator does not have neutrino origin?!" << std::endl;
    }

    // Get the associated GTruth
    std::vector<art::Ptr<simb::GTruth>> gt_v = mct_to_gt.at(i);
    assert(gt_v.size() == 1);
    auto gt = gt_v[0];

    std::cout << "Neutrino info: "
              << "Pdg " << mct_v[i]->GetNeutrino().Nu().PdgCode()
              << ", CCNC " << mct_v[i]->GetNeutrino().CCNC() 
              << ", Mode " << mct_v[i]->GetNeutrino().Mode()
              << ", X " << mct_v[i]->GetNeutrino().Nu().Vx()
              << ", Y " << mct_v[i]->GetNeutrino().Nu().Vy()
              << ", Z " << mct_v[i]->GetNeutrino().Nu().Vz()
              << std::endl;

    if (InDetector(mct_v[i]->GetNeutrino().Nu().Vx(), mct_v[i]->GetNeutrino().Nu().Vy(), mct_v[i]->GetNeutrino().Nu().Vz())) {
      _n_nu_in_av++;

      _nu_in_av_e = mct_v[i]->GetNeutrino().Nu().E();
      _nu_in_av_pdg = mct_v[i]->GetNeutrino().Nu().PdgCode();
      _nu_in_av_ccnc = mct_v[i]->GetNeutrino().CCNC();
      _nu_in_av_mode = mct_v[i]->GetNeutrino().Mode();
      _nu_in_av_int_type = mct_v[i]->GetNeutrino().InteractionType();
      _nu_in_av_vtx_x = mct_v[i]->GetNeutrino().Nu().Vx();
      _nu_in_av_vtx_y = mct_v[i]->GetNeutrino().Nu().Vy();
      _nu_in_av_vtx_z = mct_v[i]->GetNeutrino().Nu().Vz();
    }


    _nu_e = mct_v[i]->GetNeutrino().Nu().E();
    _nu_pdg = mct_v[i]->GetNeutrino().Nu().PdgCode();
    _nu_ccnc = mct_v[i]->GetNeutrino().CCNC();
    _nu_mode = mct_v[i]->GetNeutrino().Mode();
    _nu_int_type = mct_v[i]->GetNeutrino().InteractionType();
    _nu_x = mct_v[i]->GetNeutrino().X();
    _nu_y = mct_v[i]->GetNeutrino().Y();
    _nu_w = mct_v[i]->GetNeutrino().W();
    _nu_q2 = mct_v[i]->GetNeutrino().QSqr();
    _nu_vtx_x = mct_v[i]->GetNeutrino().Nu().Vx();
    _nu_vtx_y = mct_v[i]->GetNeutrino().Nu().Vy();
    _nu_vtx_z = mct_v[i]->GetNeutrino().Nu().Vz();
    _nu_px = mct_v[i]->GetNeutrino().Nu().Px();
    _nu_py = mct_v[i]->GetNeutrino().Nu().Py();
    _nu_pz = mct_v[i]->GetNeutrino().Nu().Pz();

    // _nu_oaa = GetOffAxisAngle(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z);

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
    // _nu_prod_vtx_x = mcf->fvx - _beam_origin_x;
    // _nu_prod_vtx_y = mcf->fvy - _beam_origin_y;
    // _nu_prod_vtx_z = mcf->fvz - _beam_origin_z;

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

    std::cout << "Pars out: " << std::endl;

    for (int p = 0; p < mct_v[i]->NParticles(); p++) {
      auto const & mcp = mct_v[i]->GetParticle(p);

      if (mcp.StatusCode() != 1) continue;

      std::cout << "\t"
              << "Pdg " << mcp.PdgCode()
              << ", E " << mcp.E()
              << ", X " << mcp.Vx()
              << ", Y " << mcp.Vy()
              << ", Z " << mcp.Vz()
              << std::endl;

      if (InDetector(mcp.Vx(), mcp.Vy(), mcp.Vz()) and std::abs(mcp.PdgCode()) == 2212) {
        _n_nu_np_av++;
      }

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
        // _nu_e_reco_kin = GetEnergyQE(mcp.E(), mcp.Px(), mcp.Py(), mcp.Pz());
      }

    }


    _tree->Fill();
  }



  //
  // OpFlash
  //
  // art::Handle<std::vector<recob::OpHit>> ophit_h;
  // e.getByLabel("ophitpmt", ophit_h);
  // if(!ophit_h.isValid()){
  //   std::cout << "Invalid producer for recob::OpHit: " << "ophitpmt" << ". Ignoring." << std::endl;
  // }
  // std::vector<art::Ptr<recob::OpHit>> ophit_v;
  // art::fill_ptr_vector(ophit_v, ophit_h);


  _flash_tpc0_time = -1;
  _flash_tpc0_total_pe = -1;
  // _flash_tpc0_pe_v = -1;
  _flash_tpc0_y = -1;
  _flash_tpc0_yerr = -1;
  _flash_tpc0_z = -1;
  _flash_tpc0_zerr = -1;
  _flash_tpc1_time = -1;
  _flash_tpc1_total_pe = -1;
  // _flash_tpc1_pe_v = -1;
  _flash_tpc1_y = -1;
  _flash_tpc1_yerr = -1;
  _flash_tpc1_z = -1;
  _flash_tpc1_zerr = -1;

  for (size_t l = 0; l < _flash_label_v.size(); l++) {
    art::Handle<std::vector<recob::OpFlash>> flash_h;
    e.getByLabel(_flash_label_v[l], flash_h);
    if(!flash_h.isValid()){
      std::cout << "Invalid producer for recob::OpFlash: " << _flash_label_v[l] << ". Ignoring." << std::endl;
    }

    art::FindManyP<recob::OpHit> flashToOpHitAssns(flash_h, e, _ophit_label_v[l]);

    for (size_t i = 0; i < flash_h->size(); i++) {
      auto const& f = (*flash_h)[i];
      std::cout << "Flash " << i << ", time " << f.AbsTime() << std::endl;

      if (f.AbsTime() < _flash_time_window[0] or f.AbsTime() > _flash_time_window[1]) {
        continue;
      }

      if (_tpc_v[l] == 0) {
        std::cout << "\t Keeping this one as tpc " << _tpc_v[l] << std::endl;
        _flash_tpc0_time = f.AbsTime();
        _flash_tpc0_total_pe = f.TotalPE();
        _flash_tpc0_pe_v = f.PEs();
        _flash_tpc0_y = f.YCenter();
        _flash_tpc0_yerr = f.YWidth();
        _flash_tpc0_z = f.ZCenter();
        _flash_tpc0_zerr = f.ZWidth();
      } else {
        _flash_tpc1_time = f.AbsTime();
        _flash_tpc1_total_pe = f.TotalPE();
        _flash_tpc1_pe_v = f.PEs();
        _flash_tpc1_y = f.YCenter();
        _flash_tpc1_yerr = f.YWidth();
        _flash_tpc1_z = f.ZCenter();
        _flash_tpc1_zerr = f.ZWidth();
      }   
    }
  }

  //
  // Fill the protontree
  //
  _protontree->Fill();

}



void ProtonAnalyzer::beginSubRun(art::SubRun const& sr) {

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

void ProtonAnalyzer::endJob() {

}

void ProtonAnalyzer::endSubRun(art::SubRun const& sr) {

}


bool ProtonAnalyzer::InDetector(const double& x,
                               const double& y,
                               const double& z) const {
  return !( x > _x_max || x < _x_min ||
            z > _z_max || z < _z_min ||
            y > _y_max || y < _y_min );
}

bool ProtonAnalyzer::InDetector(art::Ptr<simb::MCParticle> mcp, int & step) {
  auto t = mcp->Trajectory();
  for (size_t i = 0; i < t.size(); i++) {
    if (InDetector(t.X(i), t.Y(i), t.Z(i))) {
      step = i;
      return true;
    }
  }
  return false;
}


DEFINE_ART_MODULE(ProtonAnalyzer)
















