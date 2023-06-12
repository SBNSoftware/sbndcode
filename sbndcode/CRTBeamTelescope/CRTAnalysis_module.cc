////////////////////////////////////////////////////////////////////////
// Class:       CRTAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTAnalysis_module.cc
//
// Generated at Tue Feb  1 17:05:06 2022 by Marco Del Tutto using cetskelgen
// from  version .
// Modified by Jiaoyan Li/李 娇瑒 (jiaoyang.li@ed.ac.uk)
////////////////////////////////////////////////////////////////////////

// Search for **To-do**
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
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"

#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

class CRTAnalysis;

class CRTAnalysis : public art::EDAnalyzer {
public:
  explicit CRTAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTAnalysis(CRTAnalysis const&) = delete;
  CRTAnalysis(CRTAnalysis&&) = delete;
  CRTAnalysis& operator=(CRTAnalysis const&) = delete;
  CRTAnalysis& operator=(CRTAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(art::SubRun const& sr) override;

private:

  art::ServiceHandle<cheat::ParticleInventoryService> _particleInventory;

  sbnd::CRTBackTracker _crt_back_tracker;

  std::string _feb_data_producer;

  std::string _mctruth_label;
  std::string _g4_label;
  std::string _auxdethit_label;
  // std::string _crtdata_label;
  std::string _crthit_label;
  std::string _crttrack_label;
  std::string _febdata_label;
  std::string _pot_label;

  bool _debug;
  bool _keep_mcp_from_adh; ///< Whether or not to keep only MCParticles that produce energy deposit on the CRT
  bool _data_mode;

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

  float _mct_sp_pdg; ///< Single particle PDG
  float _mct_sp_e;   ///< Single particle energy
  float _mct_sp_px;  ///< Single particle momentum X
  float _mct_sp_py;  ///< Single particle momentum Y
  float _mct_sp_pz;  ///< Single particle momentum Z
  float _mct_sp_vx;  ///< Single particle vertex X
  float _mct_sp_vy;  ///< Single particle vertex Y
  float _mct_sp_vz;  ///< Single particle vertex Z

  float _mct_sp_2_pdg; ///< For the double-muon studys, the second particle PDG
  float _mct_sp_2_e;   ///< For the double-muon studys, the second particle energy
  float _mct_sp_2_px;  ///< For the double-muon studys, the second particle momentum X
  float _mct_sp_2_py;  ///< For the double-muon studys, the second particle momentum Y
  float _mct_sp_2_pz;  ///< For the double-muon studys, the second particle momentum Z
  float _mct_sp_2_vx;  ///< For the double-muon studys, the second particle vertex X
  float _mct_sp_2_vy;  ///< For the double-muon studys, the second particle vertex Y
  float _mct_sp_2_vz;  ///< For the double-muon studys, the second particle vertex Z

  float _mct_darkNeutrino_e;
  float _weight; //the weight which store as the vertex of dark neutrino.

  std::vector<int> _mcp_pdg; ///< G4 MCParticle PDG
  std::vector<double> _mcp_e; ///< G4 MCParticle Energy
  //std::vector<std::vector<double>> _mcp_e_vec; ///< G4 MCParticle Energy per trajectory point.
  std::vector<double> _mcp_px; ///< G4 MCParticle Momentum along X
  std::vector<double> _mcp_py; ///< G4 MCParticle Momentum along Y
  std::vector<double> _mcp_pz; ///< G4 MCParticle Momentum along Z
  std::vector<double> _mcp_startx; ///< G4 MCParticle start point X
  std::vector<double> _mcp_starty; ///< G4 MCParticle start point Y
  std::vector<double> _mcp_startz; ///< G4 MCParticle start point Z
  std::vector<double> _mcp_endx; ///< G4 MCParticle end point X
  std::vector<double> _mcp_endy; ///< G4 MCParticle end point Y
  std::vector<double> _mcp_endz; ///< G4 MCParticle end point Z
  std::vector<int> _mcp_isprimary; ///< G4 MCParticle, true if primary
  std::vector<int> _mcp_trackid; ///< G4 MCParticle, track ID
  std::vector<int> _mcp_makes_adh; ///< G4 MCParticle, true if this MCP deposits energy in CRT

  std::vector<double> _adh_t; ///< AuxDetHit time
  std::vector<double> _adh_e; ///< AuxDetHit deposited energy
  std::vector<double> _adh_x; ///< AuxDetHit X
  std::vector<double> _adh_y; ///< AuxDetHit Y
  std::vector<double> _adh_z; ///< AuxDetHit Z
  std::vector<double> _adh_p_exit; ///< AuxDetHit exit momentum
  std::vector<int> _adh_trackid; ///< AuxDetHit trackID of MCParticle

  std::vector<double> _chit_x; ///< CRT hit x
  std::vector<double> _chit_y; ///< CRT hit y
  std::vector<double> _chit_z; ///< CRT hit z
  std::vector<double> _chit_ex; ///< CRT hit error on x
  std::vector<double> _chit_ey; ///< CRT hit error on y
  std::vector<double> _chit_ez; ///< CRT hit error on z
  std::vector<double> _chit_t0; ///< CRT hit t0
  std::vector<double> _chit_t1; ///< CRT hit t1
  std::vector<double> _chit_t1_diff; ///< CRT hit difference between the two 1D hit t1s
  std::vector<uint64_t> _chit_unix_s; ///< CRT hit unix timestamp
  std::vector<double> _chit_h1_t0; ///< CRT hit t0 (1DHit 1)
  std::vector<double> _chit_h2_t0; ///< CRT hit t0 (1DHit 2)
  std::vector<double> _chit_h1_t1; ///< CRT hit t1 (1DHit 1)
  std::vector<double> _chit_h2_t1; ///< CRT hit t1 (1DHit 2)
  std::vector<double> _chit_pes; ///< CRT hit PEs
  std::vector<int> _chit_plane; ///< CRT hit plane
  std::vector<float> _chit_true_t; ///< CRT hit true time (from sim energy dep)
  std::vector<float> _chit_true_e; ///< CRT hit true energy (from sim energy dep)
  std::vector<float> _chit_true_x; ///< CRT hit true x position (from sim energy dep)
  std::vector<float> _chit_true_y; ///< CRT hit true y position (from sim energy dep)
  std::vector<float> _chit_true_z; ///< CRT hit true z position (from sim energy dep)
  // std::vector<std::vector<int> > _chit_true_mcp_trackids; ///< CRT hit contributing true trackIDs
  // std::vector<std::vector<int> > _chit_true_mcp_pdg; ///< CRT hit contributing PDG codes
  // std::vector<std::vector<double> > _chit_true_mcp_e; ///< CRT hit contributing MCP's Energy
  // std::vector<std::vector<double> > _chit_true_mcp_px; ///< CRT hit contributing MCP's Momentum along X
  // std::vector<std::vector<double> > _chit_true_mcp_py; ///< CRT hit contributing MCP's Momentum along Y
  // std::vector<std::vector<double> > _chit_true_mcp_pz; ///< CRT hit contributing MCP's Momentum along Z
  // std::vector<std::vector<double> > _chit_true_mcp_startx; ///< CRT hit contributing MCP's start point X
  // std::vector<std::vector<double> > _chit_true_mcp_starty; ///< CRT hit contributing MCP's start point Y
  // std::vector<std::vector<double> > _chit_true_mcp_startz; ///< CRT hit contributing MCP's start point Z
  // std::vector<std::vector<double> > _chit_true_mcp_endx; ///< CRT hit contributing MCP's end point X
  // std::vector<std::vector<double> > _chit_true_mcp_endy; ///< CRT hit contributing MCP's end point Y
  // std::vector<std::vector<double> > _chit_true_mcp_endz; ///< CRT hit contributing MCP's end point Z
  // std::vector<std::vector<int> > _chit_true_mcp_isprimary; ///< CRT hit contributing MCP's, true if primary
  std::vector<std::vector<uint16_t> > _chit_sipm_raw_adc; ///< The 4 contributing ADC values to this CRTHit before corrections are made or pedestals subtracted
  std::vector<std::vector<uint16_t> > _chit_sipm_adc; ///< The 4 contributing ADC values to this CRTHit before corrections are made (but after pedestal subtraction)
  std::vector<std::vector<uint16_t> > _chit_sipm_corr_adc; ///< The 4 contributing ADC values to this CRTHit
  std::vector<std::vector<uint16_t> > _chit_sipm_channel_id; ///< The IDs of the four SiPMs that were used to make this CRTHit
  std::vector<std::vector<uint16_t> > _chit_sipm_feb_mac5; ///< The IDs of the two FEBs that were used in the making of this hit.

  std::vector<double> _ct_time; ///< CRT track time
  std::vector<double> _ct_pes; ///< CRT track PEs
  std::vector<double> _ct_length; ///< CRT track length
  std::vector<double> _ct_tof; ///< CRT track time of flight
  std::vector<double> _ct_true_tof;  ///< CRT track time of flight, true
  std::vector<double> _ct_hit1_x; ///< CRT track x position of hit 1
  std::vector<double> _ct_hit1_y; ///< CRT track y position of hit 1
  std::vector<double> _ct_hit1_z; ///< CRT track z position of hit 1
  std::vector<double> _ct_hit2_x; ///< CRT track x position of hit 2
  std::vector<double> _ct_hit2_y; ///< CRT track y position of hit 2
  std::vector<double> _ct_hit2_z; ///< CRT track z position of hit 2
  std::vector<double> _ct_hit1_ex; ///< CRT track x position error of hit 1
  std::vector<double> _ct_hit1_ey; ///< CRT track y position error of hit 1
  std::vector<double> _ct_hit1_ez; ///< CRT track z position error of hit 1
  std::vector<double> _ct_hit2_ex; ///< CRT track x position error of hit 2
  std::vector<double> _ct_hit2_ey; ///< CRT track y position error of hit 2
  std::vector<double> _ct_hit2_ez; ///< CRT track z position error of hit 2
  std::vector<double> _ct_hit1_t0; ///< CRT track t0 of hit 1
  std::vector<double> _ct_hit1_t1; ///< CRT track t1 of hit 1
  std::vector<double> _ct_hit2_t0; ///< CRT track t0 of hit 2
  std::vector<double> _ct_hit2_t1; ///< CRT track t1 of hit 2
  std::vector<uint16_t> _ct_hit1_nhits; ///< CRT track, number of hits merged for 'hit 1'
  std::vector<uint16_t> _ct_hit2_nhits; ///< CRT track, number of hits merged for 'hit 2'
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_raw_adc; ///< CRT track, 4 ADC values from hit 1 before corrections are made or pedestals subtracted
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_adc; ///< CRT track, 4 ADC values from hit 1 before corrections are made (but with pedestals subtracted)
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_corr_adc; ///< CRT track, 4 ADC values from hit 1
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_raw_adc; ///< CRT track, 4 ADC values from hit 2 before corrections are made or pedestals subtracted
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_adc; ///< CRT track, 4 ADC values from hit 2 before corrections are made (but with pedestals subtracted)
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_corr_adc; ///< CRT track, 4 ADC values from hit 2
  // truth info for crt_track
  std::vector<int> _ct_pdg; ///< CRT track, truth information of the pdg code 
  std::vector<double> _ct_energy; ///< CRT track, truth information of the particle energy
  std::vector<double> _ct_deposited_energy; ///< CRT track, truth information of the deposited energy for both upstream and downstream
  std::vector<double> _ct_purity; ///< CRT track, truth information of selection purity

  std::vector<uint16_t> _feb_mac5; ///< FEBData Mac5 ID
  std::vector<uint16_t> _feb_flags; ///< FEBData Flags
  std::vector<uint32_t> _feb_ts0; ///< FEBData Ts0
  std::vector<uint32_t> _feb_ts1; ///< FEBData Ts1
  std::vector<uint32_t> _feb_unixs; ///< FEBData UnixS
  std::vector<std::vector<uint16_t>> _feb_adc; ///< FEBData 32 ADC values
  std::vector<uint32_t> _feb_coinc; ///< FEBData Coinc

  // std::vector<uint16_t> _crt_channel; ///< crtData channel ID
  // std::vector<uint32_t> _crt_t0; ///< crtData T0
  // std::vector<uint32_t> _crt_t1; ///< crtData T1
  // std::vector<uint32_t> _crt_adc; ///< crtData 32 ADC values
  // int _crt_n_up; ///< crtData Number of SiPMs in the upstream tagger
  // int _crt_n_dw; ///< crtData Number of SiPMs in the downstream tagger
  // std::vector<int> _crt_true_nide; ///< crtData, assns truth, number of IDEs
  // std::vector<double> _crt_true_t; ///< crtData, assns truth, IDE time
  // std::vector<double> _crt_true_e; ///< crtData, assns truth, IDE deposited energy
  // std::vector<double> _crt_true_x; ///< crtData, assns truth, IDE x
  // std::vector<double>_crt_true_y; ///< crtData, assns truth, IDE y
  // std::vector<double>_crt_true_z; ///< crtData, assns truth, IDE z
  // std::vector<double>_crt_true_exitp; ///< crtData, assns truth, IDE exit momentum
  // std::vector<int> _crt_true_mcp_pdg; ///< crtData, assns truth, MCP PDG
  // std::vector<int> _crt_true_mcp_isprimary; ///< crtData, assns truth, MCP is primary?
  // std::vector<std::string> _crt_true_mcp_process; ///< crtData, assns truth, MCP process
  // std::vector<int> _crt_true_mcp_mother; ///< crtData, assns truth, MCP mother ID
  // std::vector<int> _crt_true_mcp_mother_pdg; ///< crtData, assns truth, MCP mother PDG
  // std::vector<double> _crt_true_mcp_time; ///< crtData, assns truth, MCP time
  // std::vector<double> _crt_true_mcp_k_start; ///< crtData, assns truth, MCP kinetic energy at start
  // std::vector<double> _crt_true_mcp_k_end; ///< crtData, assns truth, MCP kinetic energy at end
  // std::vector<double> _crt_true_mcp_vx; ///< crtData, assns truth, MCP start X
  // std::vector<double> _crt_true_mcp_vy; ///< crtData, assns truth, MCP start Y
  // std::vector<double> _crt_true_mcp_vz; ///< crtData, assns truth, MCP start Z

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun
  double _sr_spills; ///< Number of Spills per subrun
};


CRTAnalysis::CRTAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  //_feb_data_producer = p.get<std::string>("FEBDataProducer", fhicl::ParameterSet());
  _mctruth_label     = p.get<std::string>("MCTruthLabel", "generator");
  _g4_label          = p.get<std::string>("G4Label", "largeant");
  _auxdethit_label   = p.get<std::string>("AuxDetHitLabel", "largeant:LArG4DetectorServicevolAuxDetSensitiveCRTStripBERN");
  // _crtdata_label  = p.get<std::string>("CRTDataLabel", "crt");
  _febdata_label     = p.get<std::string>("FEBDataLabel", "crtsim");
  _crthit_label      = p.get<std::string>("CRTHitLabel", "crthit");
  _crttrack_label    = p.get<std::string>("CRTTrackLabel", "crttrack");
  _debug             = p.get<bool>("Debug", false);
  _keep_mcp_from_adh = p.get<bool>("KeepMCPFromADH", false);
  _data_mode         = p.get<bool>("DataMode", false);
  _pot_label         = p.get<std::string>("POTLabel", "generator");
  _crt_back_tracker  = p.get<fhicl::ParameterSet>("CRTBackTracker", fhicl::ParameterSet());

  art::ServiceHandle<art::TFileService> fs;

  _tree = fs->make<TTree>("tree","");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");

  if(!_data_mode)
  {
    _tree->Branch("nu_e", &_nu_e, "nu_e/F");
    _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
    _tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
    _tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
    _tree->Branch("nu_int_type", &_nu_int_type, "nu_int_type/I");
    _tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
    _tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
    _tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
    _tree->Branch("nu_px", &_nu_px, "nu_px/F");
    _tree->Branch("nu_py", &_nu_py, "nu_py/F");
    _tree->Branch("nu_pz", &_nu_pz, "nu_pz/F");

    // MCTruth info for the signal (dark neutrino)
    _tree->Branch("mct_sp_pdg", &_mct_sp_pdg, "mct_sp_pdg/F");
    _tree->Branch("mct_sp_e", &_mct_sp_e, "mct_sp_e/F");
    _tree->Branch("mct_sp_px", &_mct_sp_px, "mct_sp_px/F");
    _tree->Branch("mct_sp_py", &_mct_sp_py, "mct_sp_py/F");
    _tree->Branch("mct_sp_pz", &_mct_sp_pz, "mct_sp_pz/F");
    _tree->Branch("mct_sp_vx", &_mct_sp_vx, "mct_sp_vx/F");
    _tree->Branch("mct_sp_vy", &_mct_sp_vy, "mct_sp_vy/F");
    _tree->Branch("mct_sp_vz", &_mct_sp_vz, "mct_sp_vz/F");

    _tree->Branch("mct_sp_2_pdg", &_mct_sp_2_pdg, "mct_sp_2_pdg/F");
    _tree->Branch("mct_sp_2_e", &_mct_sp_2_e, "mct_sp_2_e/F");
    _tree->Branch("mct_sp_2_py", &_mct_sp_2_py, "mct_sp_2_py/F");
    _tree->Branch("mct_sp_2_pz", &_mct_sp_2_pz, "mct_sp_2_pz/F");
    _tree->Branch("mct_sp_2_vx", &_mct_sp_2_vx, "mct_sp_2_vx/F");
    _tree->Branch("mct_sp_2_vy", &_mct_sp_2_vy, "mct_sp_2_vy/F");
    _tree->Branch("mct_sp_2_vz", &_mct_sp_2_vz, "mct_sp_2_vz/F");

    _tree->Branch("mct_darkNeutrino_e", &_mct_darkNeutrino_e, "mct_darkNeutrino_e/F");
    _tree->Branch("weight", &_weight, "weight/F");

    // MCParticle truth info. 
    _tree->Branch("mcp_pdg", "std::vector<int>", &_mcp_pdg);
    _tree->Branch("mcp_e", "std::vector<double>", &_mcp_e);
    //_tree->Branch("mcp_e_vec", "std::vector<std::vector<double>>", &_mcp_e_vec);
    _tree->Branch("mcp_px", "std::vector<double>", &_mcp_px);
    _tree->Branch("mcp_py", "std::vector<double>", &_mcp_py);
    _tree->Branch("mcp_pz", "std::vector<double>", &_mcp_pz);
    _tree->Branch("mcp_startx", "std::vector<double>", &_mcp_startx);
    _tree->Branch("mcp_starty", "std::vector<double>", &_mcp_starty);
    _tree->Branch("mcp_startz", "std::vector<double>", &_mcp_startz);
    _tree->Branch("mcp_endx", "std::vector<double>", &_mcp_endx);
    _tree->Branch("mcp_endy", "std::vector<double>", &_mcp_endy);
    _tree->Branch("mcp_endz", "std::vector<double>", &_mcp_endz);
    _tree->Branch("mcp_isprimary", "std::vector<int>", &_mcp_isprimary);
    _tree->Branch("mcp_trackid", "std::vector<int>", &_mcp_trackid);
    _tree->Branch("mcp_makes_adh", "std::vector<int>", &_mcp_makes_adh);

    _tree->Branch("adh_t", "std::vector<double>", &_adh_t);
    _tree->Branch("adh_e", "std::vector<double>", &_adh_e);
    _tree->Branch("adh_x", "std::vector<double>", &_adh_x);
    _tree->Branch("adh_y", "std::vector<double>", &_adh_y);
    _tree->Branch("adh_z", "std::vector<double>", &_adh_z);
    _tree->Branch("adh_p_exit", "std::vector<double>", &_adh_p_exit);
    _tree->Branch("adh_trackid", "std::vector<int>", &_adh_trackid);

  }

  _tree->Branch("chit_x", "std::vector<double>", &_chit_x);
  _tree->Branch("chit_y", "std::vector<double>", &_chit_y);
  _tree->Branch("chit_z", "std::vector<double>", &_chit_z);
  _tree->Branch("chit_ex", "std::vector<double>", &_chit_ex);
  _tree->Branch("chit_ey", "std::vector<double>", &_chit_ey);
  _tree->Branch("chit_ez", "std::vector<double>", &_chit_ez);
  _tree->Branch("chit_t0", "std::vector<double>", &_chit_t0);
  _tree->Branch("chit_t1", "std::vector<double>", &_chit_t1);
  _tree->Branch("chit_t1_diff", "std::vector<double>", &_chit_t1_diff);
  _tree->Branch("chit_unix_s", "std::vector<uint64_t>", &_chit_unix_s);
  _tree->Branch("chit_h1_t0", "std::vector<double>", &_chit_h1_t0);
  _tree->Branch("chit_h2_t0", "std::vector<double>", &_chit_h2_t0);
  _tree->Branch("chit_h1_t1", "std::vector<double>", &_chit_h1_t1);
  _tree->Branch("chit_h2_t1", "std::vector<double>", &_chit_h2_t1);
  _tree->Branch("chit_pes", "std::vector<double>", &_chit_pes);
  _tree->Branch("chit_plane", "std::vector<int>", &_chit_plane);
  if(!_data_mode)
  {
    _tree->Branch("chit_true_t", "std::vector<float>", &_chit_true_t);
    _tree->Branch("chit_true_e", "std::vector<float>", &_chit_true_e);
    _tree->Branch("chit_true_x", "std::vector<float>", &_chit_true_x);
    _tree->Branch("chit_true_y", "std::vector<float>", &_chit_true_y);
    _tree->Branch("chit_true_z", "std::vector<float>", &_chit_true_z);
    // _tree->Branch("chit_true_mcp_trackids", "std::vector<std::vector<int> >", &_chit_true_mcp_trackids);
    // _tree->Branch("chit_true_mcp_pdg", "std::vector<std::vector<int> >", &_chit_true_mcp_pdg);
    // _tree->Branch("chit_true_mcp_e", "std::vector<std::vector<double> >", &_chit_true_mcp_e);
    // _tree->Branch("chit_true_mcp_px", "std::vector<std::vector<double> >", &_chit_true_mcp_px);
    // _tree->Branch("chit_true_mcp_py", "std::vector<std::vector<double> >", &_chit_true_mcp_py);
    // _tree->Branch("chit_true_mcp_pz", "std::vector<std::vector<double> >", &_chit_true_mcp_pz);
    // _tree->Branch("chit_true_mcp_startx", "std::vector<std::vector<double> >", &_chit_true_mcp_startx);
    // _tree->Branch("chit_true_mcp_starty", "std::vector<std::vector<double> >", &_chit_true_mcp_starty);
    // _tree->Branch("chit_true_mcp_startz", "std::vector<std::vector<double> >", &_chit_true_mcp_startz);
    // _tree->Branch("chit_true_mcp_endx", "std::vector<std::vector<double> >", &_chit_true_mcp_endx);
    // _tree->Branch("chit_true_mcp_endy", "std::vector<std::vector<double> >", &_chit_true_mcp_endy);
    // _tree->Branch("chit_true_mcp_endz", "std::vector<std::vector<double> >", &_chit_true_mcp_endz);
    // _tree->Branch("chit_true_mcp_isprimary", "std::vector<std::vector<int> >", &_chit_true_mcp_isprimary);
  }
  _tree->Branch("chit_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_raw_adc);
  _tree->Branch("chit_sipm_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_adc);
  _tree->Branch("chit_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_corr_adc);
  _tree->Branch("chit_sipm_channel_id", "std::vector<std::vector<uint16_t> >", &_chit_sipm_channel_id);
  _tree->Branch("chit_sipm_feb_mac5", "std::vector<std::vector<uint16_t> >", &_chit_sipm_feb_mac5);

  _tree->Branch("ct_time", "std::vector<double>", &_ct_time);
  _tree->Branch("ct_pes", "std::vector<double>", &_ct_pes);
  _tree->Branch("ct_length", "std::vector<double>", &_ct_length);
  _tree->Branch("ct_tof", "std::vector<double>", &_ct_tof);
  if(!_data_mode) _tree->Branch("ct_true_tof", "std::vector<double>", &_ct_true_tof);
  _tree->Branch("ct_hit1_x", "std::vector<double>", &_ct_hit1_x);
  _tree->Branch("ct_hit1_y", "std::vector<double>", &_ct_hit1_y);
  _tree->Branch("ct_hit1_z", "std::vector<double>", &_ct_hit1_z);
  _tree->Branch("ct_hit2_x", "std::vector<double>", &_ct_hit2_x);
  _tree->Branch("ct_hit2_y", "std::vector<double>", &_ct_hit2_y);
  _tree->Branch("ct_hit2_z", "std::vector<double>", &_ct_hit2_z);
  _tree->Branch("ct_hit1_ex", "std::vector<double>", &_ct_hit1_ex);
  _tree->Branch("ct_hit1_ey", "std::vector<double>", &_ct_hit1_ey);
  _tree->Branch("ct_hit1_ez", "std::vector<double>", &_ct_hit1_ez);
  _tree->Branch("ct_hit2_ex", "std::vector<double>", &_ct_hit2_ex);
  _tree->Branch("ct_hit2_ey", "std::vector<double>", &_ct_hit2_ey);
  _tree->Branch("ct_hit2_ez", "std::vector<double>", &_ct_hit2_ez);
  _tree->Branch("ct_hit1_t0", "std::vector<double>", &_ct_hit1_t0);
  _tree->Branch("ct_hit1_t1", "std::vector<double>", &_ct_hit1_t1);
  _tree->Branch("ct_hit2_t0", "std::vector<double>", &_ct_hit2_t0);
  _tree->Branch("ct_hit2_t1", "std::vector<double>", &_ct_hit2_t1);
  _tree->Branch("ct_hit1_nhits", "std::vector<uint16_t>", &_ct_hit1_nhits);
  _tree->Branch("ct_hit2_nhits", "std::vector<uint16_t>", &_ct_hit2_nhits);
  _tree->Branch("ct_hit1_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_raw_adc);
  _tree->Branch("ct_hit1_sipm_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_adc);
  _tree->Branch("ct_hit1_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_corr_adc);
  _tree->Branch("ct_hit2_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_raw_adc);
  _tree->Branch("ct_hit2_sipm_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_adc);
  _tree->Branch("ct_hit2_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_corr_adc);
  // truth information
  if(!_data_mode){
    _tree->Branch("ct_pdg", "std::vector<int>", &_ct_pdg);
    _tree->Branch("ct_energy", "std::vector<double>", &_ct_energy);
    _tree->Branch("ct_deposited_energy", "std::vector<double>", &_ct_deposited_energy);
    _tree->Branch("ct_purity", "std::vector<double>", &_ct_purity);
  }

  _tree->Branch("feb_mac5", "std::vector<uint16_t>", &_feb_mac5);
  _tree->Branch("feb_flags", "std::vector<uint16_t>", &_feb_flags);
  _tree->Branch("feb_ts0", "std::vector<uint32_t>", &_feb_ts0);
  _tree->Branch("feb_ts1", "std::vector<uint32_t>", &_feb_ts1);
  _tree->Branch("feb_unixs", "std::vector<uint32_t>", &_feb_unixs);
  _tree->Branch("feb_adc", "std::vector<std::vector<uint16_t>>", &_feb_adc);
  _tree->Branch("feb_coinc", "std::vector<uint32_t>", &_feb_coinc);

  // _tree->Branch("crt_channel", "std::vector<uint16_t>", &_crt_channel);
  // _tree->Branch("crt_t0", "std::vector<uint32_t>", &_crt_t0);
  // _tree->Branch("crt_t1", "std::vector<uint32_t>", &_crt_t1);
  // _tree->Branch("crt_adc", "std::vector<uint32_t>", &_crt_adc);
  // _tree->Branch("crt_n_up", &_crt_n_up, "crt_n_up/I");
  // _tree->Branch("crt_n_dw", &_crt_n_dw, "crt_n_dw/I");
  // _tree->Branch("crt_true_nide", "std::vector<int>", &_crt_true_nide);
  // _tree->Branch("crt_true_t", "std::vector<double>", &_crt_true_t);
  // _tree->Branch("crt_true_e", "std::vector<double>", &_crt_true_e);
  // _tree->Branch("crt_true_x", "std::vector<double>", &_crt_true_x);
  // _tree->Branch("crt_true_y", "std::vector<double>", &_crt_true_y);
  // _tree->Branch("crt_true_z", "std::vector<double>", &_crt_true_z);
  // _tree->Branch("crt_true_exitp", "std::vector<double>", &_crt_true_exitp);
  // _tree->Branch("crt_true_mcp_pdg", "std::vector<int>", &_crt_true_mcp_pdg);
  // _tree->Branch("crt_true_mcp_isprimary", "std::vector<int>", &_crt_true_mcp_isprimary);
  // _tree->Branch("crt_true_mcp_process", "std::vector<std::string>", &_crt_true_mcp_process);
  // _tree->Branch("crt_true_mcp_mother", "std::vector<int>", &_crt_true_mcp_mother);
  // _tree->Branch("crt_true_mcp_mother_pdg", "std::vector<int>", &_crt_true_mcp_mother_pdg);
  // _tree->Branch("crt_true_mcp_time", "std::vector<double>", &_crt_true_mcp_time);
  // _tree->Branch("crt_true_mcp_k_start", "std::vector<double>", &_crt_true_mcp_k_start);
  // _tree->Branch("crt_true_mcp_k_end", "std::vector<double>", &_crt_true_mcp_k_end);
  // _tree->Branch("crt_true_mcp_vx", "std::vector<double>", &_crt_true_mcp_vx);
  // _tree->Branch("crt_true_mcp_vy", "std::vector<double>", &_crt_true_mcp_vy);
  // _tree->Branch("crt_true_mcp_vz", "std::vector<double>", &_crt_true_mcp_vz);

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void CRTAnalysis::analyze(art::Event const& e)
{ 

  _crt_back_tracker.Initialize(e); // Initialise the backtrack alg. 

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if(_debug) std::cout << "This is event " << _event << std::endl;

  art::Handle<std::vector<simb::MCTruth>> mct_h;
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  std::map<int, simb::MCParticle> trackid_to_mcp;
  art::Handle<std::vector<sim::AuxDetHit>> adh_h;
  std::vector<art::Ptr<sim::AuxDetHit>> adh_v;

  //
  // Get the MCTruth
  //
  if(!_data_mode)
  {
    e.getByLabel(_mctruth_label, mct_h);
    if(!mct_h.isValid()){
      std::cout << "MCTruth product " << _mctruth_label << " not found..." << std::endl;
      throw std::exception();
    }
    art::fill_ptr_vector(mct_v, mct_h);

    //
    // Get the MCParticles from G4
    //
    e.getByLabel(_g4_label, mcp_h);
    if(!mcp_h.isValid()){
      std::cout << "MCTruth product " << _g4_label << " not found..." << std::endl;
      throw std::exception();
    }
    art::fill_ptr_vector(mcp_v, mcp_h);

    std::map<int, simb::MCParticle> trackid_to_mcp;
    for (auto mcp : mcp_v) { trackid_to_mcp[mcp->TrackId()] = *mcp; }

    //
    // Get the AuxDetHits from G4
    //
    e.getByLabel(_auxdethit_label, adh_h);
    if(!adh_h.isValid()){
      std::cout << "AuxDetHit product " << _auxdethit_label << " not found..." << std::endl;
      throw std::exception();
    }
    art::fill_ptr_vector(adh_v, adh_h);
  }

  //
  // Get the CRT Hits
  //
  art::Handle<std::vector<sbn::crt::CRTHit>> crt_hit_h;
  e.getByLabel(_crthit_label, crt_hit_h);
  if(!crt_hit_h.isValid()){
    std::cout << "CRTHit product " << _crthit_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTHit>> crt_hit_v;
  art::fill_ptr_vector(crt_hit_v, crt_hit_h);

  // Get the CRTHit to FEBData association
  art::FindManyP<sbnd::crt::FEBData> crt_hit_to_feb_data(crt_hit_h, e, _crthit_label);


  //
  // Get the FEB Data
  //
  art::Handle<std::vector<sbnd::crt::FEBData>> feb_data_h;
  e.getByLabel(_febdata_label, feb_data_h);
  if(!feb_data_h.isValid()){
    std::cout << "FEBData product " << _febdata_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbnd::crt::FEBData>> feb_data_v;
  art::fill_ptr_vector(feb_data_v, feb_data_h);

  // Get the FEBData to AuxDetIDE & FEBTruthInfo association
  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> *febdata_to_ides;
  if(!_data_mode) {
    febdata_to_ides = new art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo>(feb_data_h, e, _febdata_label);
  }
  

  // // Get the CRT Hits to CRTData association
  // art::FindManyP<sbnd::crt::CRTData> crt_hit_to_data (crt_hit_h, e, _crthit_label);

  // // Get the CRTData to AuxDetIDE association
  // art::Handle<std::vector<sbnd::crt::CRTData>> crt_data_h;
  // e.getByLabel(_crtdata_label, crt_data_h);
  // art::FindManyP<sim::AuxDetIDE> crt_data_to_ides (crt_data_h, e, _crtdata_label);


  // //
  // // Get the CRT Data
  // //
  // if(!crt_data_h.isValid()){
  //   std::cout << "CRTData product " << _crtdata_label << " not found..." << std::endl;
  //   throw std::exception();
  // }
  // std::vector<art::Ptr<sbnd::crt::CRTData>> crt_data_v;
  // art::fill_ptr_vector(crt_data_v, crt_data_h);


  //
  // Get the CRT Tracks
  //
  art::Handle<std::vector<sbn::crt::CRTTrack>> crt_track_h;
  e.getByLabel(_crttrack_label, crt_track_h);
  if(!crt_track_h.isValid()){
    std::cout << "CRTTrack product " << _crttrack_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTTrack>> crt_track_v;
  art::fill_ptr_vector(crt_track_v, crt_track_h);

  // Get the CRT Tracks to Hits association
  art::FindManyP<sbn::crt::CRTHit> crt_track_to_hit (crt_track_h, e, _crttrack_label);


  if(!_data_mode)
  {
    //
    // Fill the MCTruth in the tree
    //
    assert(mct_v.size() <= 3);
    auto mct = mct_v[0];
    if (mct->Origin() == simb::kBeamNeutrino)
    {
      _nu_e = mct->GetNeutrino().Nu().E();
      _nu_pdg = mct->GetNeutrino().Nu().PdgCode();
      _nu_ccnc = mct->GetNeutrino().CCNC();
      _nu_mode = mct->GetNeutrino().Mode();
      _nu_int_type = mct->GetNeutrino().InteractionType();
      _nu_vtx_x = mct->GetNeutrino().Nu().Vx();
      _nu_vtx_y = mct->GetNeutrino().Nu().Vy();
      _nu_vtx_z = mct->GetNeutrino().Nu().Vz();
      _nu_px = mct->GetNeutrino().Nu().Px();
      _nu_py = mct->GetNeutrino().Nu().Py();
      _nu_pz = mct->GetNeutrino().Nu().Pz();

      // In the neutrino case, save the outgoing lepton in the first mct_sp object
      for (int p = 0; p < mct->NParticles(); p++)
      {
        auto particle = mct->GetParticle(p);
        if (std::abs(particle.PdgCode()) == 13 || std::abs(particle.PdgCode()) == 11) {
          _mct_sp_pdg = particle.PdgCode();
          _mct_sp_e = particle.E();
          _mct_sp_px = particle.Px();
          _mct_sp_py = particle.Py();
          _mct_sp_pz = particle.Pz();
          _mct_sp_vx = particle.Vx();
          _mct_sp_vy = particle.Vy();
          _mct_sp_vz = particle.Vz();
        }
      }
      if (_debug) std::cout << "Neutrino E = " << _nu_e
                            << ", x = " << _nu_vtx_x
                            << ", y = " << _nu_vtx_y
                            << ", z = " << _nu_vtx_z << std::endl;

    }

    // To-do (?) comment by Jiaoyang:
    // Removed the if-condition for now as in we are using TextFileGen instead of SingleGen.
    // Thus the if-condiction is no longer hold.
    //else if (mct->Origin() == simb::kSingleParticle) {
    else {
      assert(mct->NParticles() <= 3); // die if mct->NPartickes() > 3
      auto particle = mct->GetParticle(0);
      _mct_sp_pdg = particle.PdgCode();
      _mct_sp_e = particle.E();
      _mct_sp_px = particle.Px();
      _mct_sp_py = particle.Py();
      _mct_sp_pz = particle.Pz();
      _mct_sp_vx = particle.Vx();
      _mct_sp_vy = particle.Vy();
      _mct_sp_vz = particle.Vz();

      if (mct->NParticles() > 1){ // modified to include double MIPs.
        auto particle_2 = mct->GetParticle(1);
        _mct_sp_2_pdg = particle_2.PdgCode();
        _mct_sp_2_e = particle_2.E();
        _mct_sp_2_px = particle_2.Px();
        _mct_sp_2_py = particle_2.Py();
        _mct_sp_2_pz = particle_2.Pz();
        _mct_sp_2_vx = particle_2.Vx();
        _mct_sp_2_vy = particle_2.Vy();
        _mct_sp_2_vz = particle_2.Vz();
      }

      if (mct->NParticles() > 2){
        auto particle_3 = mct->GetParticle(2);
        _mct_darkNeutrino_e = particle_3.E();
        _weight = particle_3.Vx();
      }

      if (_debug) std::cout<<mct->NParticles()<<" "<<_mct_sp_pdg<<" "<<_mct_sp_2_pdg<<std::endl;
    }


    //
    // Fill the AuxDetHits in the tree
    //
    size_t n_adh = adh_v.size();
    // std::cout << "n_adh " << n_adh << std::endl;
    _adh_t.resize(n_adh);
    _adh_e.resize(n_adh);
    _adh_x.resize(n_adh);
    _adh_y.resize(n_adh);
    _adh_z.resize(n_adh);
    _adh_p_exit.resize(n_adh);
    _adh_trackid.resize(n_adh);

    std::set<unsigned int> trackids_from_adh;

    for (size_t i = 0; i < n_adh; i++) {
      auto auxdethit = adh_v[i];
      // if (0.5 * (auxdethit->GetEntryZ() + auxdethit->GetExitZ()) > -164) continue;
      _adh_t[i] = 0.5 * (auxdethit->GetEntryT() + auxdethit->GetExitT());
      _adh_e[i] = auxdethit->GetEnergyDeposited();
      _adh_x[i] = 0.5 * (auxdethit->GetEntryX() + auxdethit->GetExitX());
      _adh_y[i] = 0.5 * (auxdethit->GetEntryY() + auxdethit->GetExitY());
      _adh_z[i] = 0.5 * (auxdethit->GetEntryZ() + auxdethit->GetExitZ());
      _adh_p_exit[i] = std::sqrt(auxdethit->GetExitMomentumX()*auxdethit->GetExitMomentumX() +
                                 auxdethit->GetExitMomentumY()*auxdethit->GetExitMomentumY() +
                                 auxdethit->GetExitMomentumZ()*auxdethit->GetExitMomentumZ());
      _adh_trackid[i] = auxdethit->GetTrackID();
      trackids_from_adh.insert(auxdethit->GetTrackID());
      if (_debug) std::cout << "Adding adh with track ID " << auxdethit->GetTrackID() << std::endl;
    }


    //
    // Fill the MCParticles in the tree
    //
    size_t n_mcp = 0; // mcp_v.size();
    // Count mcp that have status code different than zero
    for (auto mcp : mcp_v) { if (mcp->StatusCode() == 1) n_mcp++; }
    // If we only keep MCP that make energy deposits in the CRT, count only those
    if (_keep_mcp_from_adh) { n_mcp = trackids_from_adh.size(); }
    _mcp_pdg.resize(n_mcp);
    _mcp_e.resize(n_mcp);

    /*
    for (auto mcp : mcp_v){ 
      if (mcp->StatusCode() == 1) 
      {
        size_t n_trajectoryPoints = mcp->NumberTrajectoryPoints(); 
        _mcp_e_vec.resize(n_mcp, std::vector<double>(n_trajectoryPoints)); 
      }
    }*/
    _mcp_px.resize(n_mcp);
    _mcp_py.resize(n_mcp);
    _mcp_pz.resize(n_mcp);
    _mcp_startx.resize(n_mcp);
    _mcp_starty.resize(n_mcp);
    _mcp_startz.resize(n_mcp);
    _mcp_endx.resize(n_mcp);
    _mcp_endy.resize(n_mcp);
    _mcp_endz.resize(n_mcp);
    _mcp_isprimary.resize(n_mcp);
    _mcp_trackid.resize(n_mcp);
    _mcp_makes_adh.resize(n_mcp);

    size_t counter = 0;
    for (size_t i = 0; i < mcp_v.size(); i++) {
      auto particle = mcp_v[i];
      // Exclude particles that are not propagated
      if (particle->StatusCode() != 1) continue;
      // Exclude particles that don't make energy deposit in the CRT, if we only keep those
      if (_keep_mcp_from_adh and
           trackids_from_adh.find(particle->TrackId()) == trackids_from_adh.end()) {
        continue;
      }
      // if (particle->Process() != "primary" || particle->StatusCode() != 1) continue;
      _mcp_pdg[counter] = particle->PdgCode();
      _mcp_e[counter] = particle->E();

      //for (size_t j = 0; j < _mcp_e_vec[counter].size(); j++)
        //_mcp_e_vec[counter].at(j) = particle->E(j);

      _mcp_px[counter] = particle->Px();
      _mcp_py[counter] = particle->Py();
      _mcp_pz[counter] = particle->Pz();
      _mcp_startx[counter] = particle->Vx();
      _mcp_starty[counter] = particle->Vy();
      _mcp_startz[counter] = particle->Vz();
      _mcp_endx[counter] = particle->EndX();
      _mcp_endy[counter] = particle->EndY();
      _mcp_endz[counter] = particle->EndZ();
      _mcp_isprimary[counter] = particle->Process() == "primary";
      _mcp_trackid[counter] = particle->TrackId();
      _mcp_makes_adh[counter] = trackids_from_adh.find(particle->TrackId()) != trackids_from_adh.end();


      if (_debug) {
        std::cout << counter << ", " << i << "MCP " << _mcp_pdg[counter] << ", trackID = " << _mcp_trackid[counter]
                  << ", p = " << particle->P()
                  << ", (" << particle->Px() << "," << particle->Py() << "," << particle->Pz() << ")"
                  << ", process = "
                  << particle->Process() << ": start point = ("
                  << _mcp_startx[counter] << ", " << _mcp_starty[counter] << ", " << _mcp_startz[counter]
                  << ") - start process: " << particle->Process()
                  << " --- end point = (" 
                  << _mcp_endx[counter] << ", " << _mcp_endy[counter] << ", " << _mcp_endz[counter] << "); "<< std::endl;
      }
      counter++;
    }
  }



  //
  // Fill the CRT Hits in the tree
  //
  size_t n_hits = crt_hit_v.size();
  _chit_x.resize(n_hits);
  _chit_y.resize(n_hits);
  _chit_z.resize(n_hits);
  _chit_ex.resize(n_hits);
  _chit_ey.resize(n_hits);
  _chit_ez.resize(n_hits);
  _chit_t0.resize(n_hits);
  _chit_t1.resize(n_hits);
  _chit_t1_diff.resize(n_hits);
  _chit_unix_s.resize(n_hits);
  _chit_h1_t0.resize(n_hits);
  _chit_h2_t0.resize(n_hits);
  _chit_h1_t1.resize(n_hits);
  _chit_h2_t1.resize(n_hits);
  _chit_pes.resize(n_hits);
  _chit_plane.resize(n_hits);
  _chit_true_t.resize(n_hits);
  _chit_true_e.resize(n_hits);
  _chit_true_x.resize(n_hits);
  _chit_true_y.resize(n_hits);
  _chit_true_z.resize(n_hits);
  // _chit_true_mcp_trackids.resize(n_hits);
  // _chit_true_mcp_pdg.resize(n_hits);
  // _chit_true_mcp_e.resize(n_hits);
  // _chit_true_mcp_px.resize(n_hits);
  // _chit_true_mcp_py.resize(n_hits);
  // _chit_true_mcp_pz.resize(n_hits);
  // _chit_true_mcp_startx.resize(n_hits);
  // _chit_true_mcp_starty.resize(n_hits);
  // _chit_true_mcp_startz.resize(n_hits);
  // _chit_true_mcp_endx.resize(n_hits);
  // _chit_true_mcp_endy.resize(n_hits);
  // _chit_true_mcp_endz.resize(n_hits);
  // _chit_true_mcp_px.resize(n_hits);
  // _chit_true_mcp_py.resize(n_hits);
  // _chit_true_mcp_isprimary.resize(n_hits);
  _chit_sipm_raw_adc.resize(n_hits);
  _chit_sipm_adc.resize(n_hits);
  _chit_sipm_corr_adc.resize(n_hits);
  _chit_sipm_channel_id.resize(n_hits);
  _chit_sipm_feb_mac5.resize(n_hits);

  for (size_t i = 0; i < n_hits; i++) {

    auto hit = crt_hit_v[i];

    _chit_x[i] = hit->x_pos;
    _chit_y[i] = hit->y_pos;
    _chit_z[i] = hit->z_pos;
    _chit_ex[i] = hit->x_err;
    _chit_ey[i] = hit->y_err;
    _chit_ez[i] = hit->z_err;
    _chit_t0[i] = hit->ts0_ns;
    _chit_t1[i] = hit->ts1_ns;
    _chit_t1_diff[i] = hit->ts0_ns_corr;   // the variable name in the object is old and is just a placeholder for diff, don't worry!
    _chit_unix_s[i] = hit->ts0_s;
    _chit_pes[i] = hit->peshit;

    if (hit->tagger == "volTaggerNorth_0") {
      _chit_plane[i] = 0; // upstream
    } else {
      _chit_plane[i] = 1; // downstream
    }

    _chit_sipm_raw_adc[i].resize(4);
    _chit_sipm_adc[i].resize(4);
    _chit_sipm_corr_adc[i].resize(4);
    const std::array<uint16_t,4> raw_adcs  = hit->raw_adcs;
    const std::array<uint16_t,4> adcs      = hit->adcs;
    const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;

    for(uint adc_i = 0; adc_i < 4; ++adc_i)
    {
      _chit_sipm_raw_adc[i][adc_i]  = raw_adcs[adc_i];
      _chit_sipm_adc[i][adc_i]      = adcs[adc_i];
      _chit_sipm_corr_adc[i][adc_i] = corr_adcs[adc_i];
    }

    // CRTHit matches to AuxDetIDE.  **To-do** need to find the correct channel for CRT Hits. 
    auto feb_datas = crt_hit_to_feb_data.at(hit.key()); // Use the Assn to find the AuxDetIDE from FEBData. 
    if(feb_datas.size() != 2) std::cout << "ERROR: CRTHit associated to " << feb_datas.size() << " FEBDatas" << std::endl;

    _chit_sipm_feb_mac5[i].resize(feb_datas.size());

    _chit_sipm_channel_id[i].resize(2);
    _chit_sipm_channel_id[i][0] = hit->channel0;
    _chit_sipm_channel_id[i][1] = hit->channel1;

    for(unsigned ii = 0; ii < feb_datas.size(); ++ii)
    {
      const auto& feb_data = feb_datas[ii];
      _chit_sipm_feb_mac5[i][ii] = feb_data->Mac5();
    }

    // From the hit, get the associated FEBData,
    // then the associated AuxDetIDE, so we can
    // retrieve the truth info
    if(!_data_mode) {
      _chit_true_t[i] = 0;
      _chit_true_e[i] = 0;
      size_t n_ides = 0;
      for (auto feb_data : feb_datas) {
        auto ide_v = febdata_to_ides->at(feb_data.key());
        for (auto ide : ide_v) {
          _chit_true_t[i] += 0.5 * (ide->entryT + ide->exitT);
          _chit_true_e[i] += ide->energyDeposited;
          _chit_true_x[i] += 0.5 * (ide->entryX + ide->exitX);
          _chit_true_y[i] += 0.5 * (ide->entryY + ide->exitY);
          _chit_true_z[i] += 0.5 * (ide->entryZ + ide->exitZ);
          n_ides++;
        }
      }
      _chit_true_t[i] /= n_ides;
      _chit_true_e[i] /= n_ides;
      _chit_true_x[i] /= n_ides;
      _chit_true_y[i] /= n_ides;
      _chit_true_z[i] /= n_ides;
    }

    if (_debug) std::cout << "CRT hit, z = " << _chit_z[i] << ", h1 time " << _chit_h1_t1[i] << ", h2 time " << _chit_h2_t1[i] << ", hit time " << _chit_t1[i] << std::endl;
  }


  //
  // Fill the CRT Tracks in the tree
  //
  size_t n_tracks = crt_track_v.size();
  _ct_pes.resize(n_tracks);
  _ct_time.resize(n_tracks);
  _ct_length.resize(n_tracks);
  _ct_tof.resize(n_tracks);
  _ct_true_tof.resize(n_tracks);
  _ct_hit1_x.resize(n_tracks);
  _ct_hit1_y.resize(n_tracks);
  _ct_hit1_z.resize(n_tracks);
  _ct_hit2_x.resize(n_tracks);
  _ct_hit2_y.resize(n_tracks);
  _ct_hit2_z.resize(n_tracks);
  _ct_hit1_ex.resize(n_tracks);
  _ct_hit1_ey.resize(n_tracks);
  _ct_hit1_ez.resize(n_tracks);
  _ct_hit2_ex.resize(n_tracks);
  _ct_hit2_ey.resize(n_tracks);
  _ct_hit2_ez.resize(n_tracks);
  _ct_hit1_t0.resize(n_tracks);
  _ct_hit1_t1.resize(n_tracks);
  _ct_hit2_t0.resize(n_tracks);
  _ct_hit2_t1.resize(n_tracks);
  _ct_hit1_nhits.resize(n_tracks);
  _ct_hit2_nhits.resize(n_tracks);
  _ct_hit1_sipm_raw_adc.resize(n_tracks);
  _ct_hit1_sipm_adc.resize(n_tracks);
  _ct_hit1_sipm_corr_adc.resize(n_tracks);
  _ct_hit2_sipm_raw_adc.resize(n_tracks);
  _ct_hit2_sipm_adc.resize(n_tracks);
  _ct_hit2_sipm_corr_adc.resize(n_tracks);

  if(!_data_mode){
    _ct_pdg.resize(n_tracks);
    _ct_energy.resize(n_tracks);
    _ct_deposited_energy.resize(n_tracks);
    _ct_purity.resize(n_tracks);
  }



  for (size_t i = 0; i < n_tracks; ++i){

    auto track = crt_track_v[i];

    _ct_pes[i] = track->peshit;
    _ct_time[i] = track->ts1_ns;
    _ct_length[i] = track->length;
    _ct_tof[i] = track->ts0_ns_h2 - track->ts0_ns_h1;
    _ct_hit1_x[i] = track->x1_pos;
    _ct_hit1_y[i] = track->y1_pos;
    _ct_hit1_z[i] = track->z1_pos;
    _ct_hit2_x[i] = track->x2_pos;
    _ct_hit2_y[i] = track->y2_pos;
    _ct_hit2_z[i] = track->z2_pos;
    _ct_hit1_ex[i] = track->x1_err;
    _ct_hit1_ey[i] = track->y1_err;
    _ct_hit1_ez[i] = track->z1_err;
    _ct_hit2_ex[i] = track->x2_err;
    _ct_hit2_ey[i] = track->y2_err;
    _ct_hit2_ez[i] = track->z2_err;

    if (_debug) std::cout << "CRT track, z1 = " << _ct_hit1_z[i] << ", z2 = " << _ct_hit2_z[i] << ", ts0_ns_h1 = " << track->ts0_ns_h1 << ", ts0_ns_h2 = " << track->ts0_ns_h2 << ", tof = " << _ct_tof[i] << std::endl;

    auto hit_v = crt_track_to_hit.at(track.key());
    _ct_hit1_nhits[i] = 0;
    _ct_hit2_nhits[i] = 0;

    for(auto hit : hit_v)
    {
    	if(std::signbit(hit->z_pos) == std::signbit(track->z1_pos))
    	  ++_ct_hit1_nhits[i];
    	else if(std::signbit(hit->z_pos) == std::signbit(track->z2_pos))
    	  ++_ct_hit2_nhits[i];
    }

    _ct_hit1_sipm_raw_adc[i].resize(4 * _ct_hit1_nhits[i]);
    _ct_hit1_sipm_adc[i].resize(4 * _ct_hit1_nhits[i]);
    _ct_hit1_sipm_corr_adc[i].resize(4 * _ct_hit1_nhits[i]);
    _ct_hit2_sipm_raw_adc[i].resize(4 * _ct_hit2_nhits[i]);
    _ct_hit2_sipm_adc[i].resize(4 * _ct_hit2_nhits[i]);
    _ct_hit2_sipm_corr_adc[i].resize(4 * _ct_hit2_nhits[i]);

    if(!_data_mode) _ct_true_tof[i] = 0;

    unsigned used_hits_1 = 0, used_hits_2 = 0;

    for (size_t i_hit = 0; i_hit < hit_v.size(); i_hit++)
    {
    	auto hit = hit_v[i_hit];
    	  
    	if(std::signbit(hit->z_pos) == std::signbit(track->z1_pos))
  	  {
  	    _ct_hit1_t0[i] = hit->ts0_ns;
  	    _ct_hit1_t1[i] = hit->ts1_ns;

  	    const std::array<uint16_t,4> raw_adcs  = hit->raw_adcs;
  	    const std::array<uint16_t,4> adcs      = hit->adcs;
  	    const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;

  	    for(uint adc_i = 0; adc_i < 4; ++adc_i)
	      {
      		_ct_hit1_sipm_raw_adc[i][4 * used_hits_1 + adc_i]  = raw_adcs[adc_i];
      		_ct_hit1_sipm_adc[i][4 * used_hits_1 + adc_i]      = adcs[adc_i];
      		_ct_hit1_sipm_corr_adc[i][4 * used_hits_1 + adc_i] = corr_adcs[adc_i];
  	     }
  	    ++used_hits_1;
  	  }
    	else if(std::signbit(hit->z_pos) == std::signbit(track->z2_pos))
  	  {
  	    _ct_hit2_t0[i] = hit->ts0_ns;
  	    _ct_hit2_t1[i] = hit->ts1_ns;

  	    const std::array<uint16_t,4> raw_adcs  = hit->raw_adcs;
  	    const std::array<uint16_t,4> adcs      = hit->adcs;
  	    const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;

  	    for(uint adc_i = 0; adc_i < 4; ++adc_i)
	      {
      		_ct_hit2_sipm_raw_adc[i][4 * used_hits_2 + adc_i]  = raw_adcs[adc_i];
      		_ct_hit2_sipm_adc[i][4 * used_hits_2 + adc_i]      = adcs[adc_i];
      		_ct_hit2_sipm_corr_adc[i][4 * used_hits_2 + adc_i] = corr_adcs[adc_i];
	      }
  	    ++used_hits_2;
  	  }
    }
    // truth information
    if (!_data_mode){
      const sbnd::CRTBackTracker::TruthMatchMetrics truthMatch = _crt_back_tracker.TruthMatrixFromTotalEnergy(e, track);
      
      _ct_pdg[i]              = truthMatch.pdg;
      _ct_energy[i]           = truthMatch.particle_energy;
      _ct_deposited_energy[i] = truthMatch.depEnergy_total;
      _ct_purity[i]           = truthMatch.purity;
      if (_debug) std::cout<<"Number of tracks: "<<n_tracks<<" This track "<< i <<", id: "<<truthMatch.trackid<<", PDG code: "<<truthMatch.pdg<<", energy: "<<truthMatch.particle_energy
                           <<", total deposited energy: "<<truthMatch.depEnergy_total<<", selection purity: "<<truthMatch.purity<<std::endl;
    }
  }





  //
  // Fill the FEBData objects in the tree
  //
  size_t n_febdata = feb_data_v.size();

  _feb_mac5.resize(n_febdata);
  _feb_flags.resize(n_febdata);
  _feb_ts0.resize(n_febdata);
  _feb_ts1.resize(n_febdata);
  _feb_unixs.resize(n_febdata);
  _feb_adc.resize(n_febdata, std::vector<uint16_t>(32));
  _feb_coinc.resize(n_febdata);

  for (size_t i = 0; i < n_febdata; ++i){

    auto feb_data = feb_data_v[i];

    _feb_mac5[i] = feb_data->Mac5();
    _feb_flags[i] = feb_data->Flags();
    _feb_ts0[i] = feb_data->Ts0();
    _feb_ts1[i] = feb_data->Ts1();
    _feb_unixs[i] = feb_data->UnixS();
    for (size_t j = 0; j < 32; j++) {
      _feb_adc[i][j] = feb_data->ADC(j);
    }
    _feb_coinc[i] = feb_data->Coinc();
  }


  // //
  // // Fill the CRTData objects in the tree
  // //
  // size_t n_crtdata = crt_data_v.size();

  // _crt_channel.resize(n_crtdata);
  // _crt_t0.resize(n_crtdata);
  // _crt_t1.resize(n_crtdata);
  // _crt_adc.resize(n_crtdata);
  // _crt_true_nide.resize(n_crtdata);
  // _crt_true_t.resize(n_crtdata);
  // _crt_true_e.resize(n_crtdata);
  // _crt_true_x.resize(n_crtdata);
  // _crt_true_y.resize(n_crtdata);
  // _crt_true_z.resize(n_crtdata);
  // _crt_true_exitp.resize(n_crtdata);
  // _crt_true_mcp_pdg.resize(n_crtdata);
  // _crt_true_mcp_isprimary.resize(n_crtdata);
  // _crt_true_mcp_process.resize(n_crtdata);
  // _crt_true_mcp_mother.resize(n_crtdata);
  // _crt_true_mcp_mother_pdg.resize(n_crtdata);
  // _crt_true_mcp_time.resize(n_crtdata);
  // _crt_true_mcp_k_start.resize(n_crtdata);
  // _crt_true_mcp_k_end.resize(n_crtdata);
  // _crt_true_mcp_vx.resize(n_crtdata);
  // _crt_true_mcp_vy.resize(n_crtdata);
  // _crt_true_mcp_vz.resize(n_crtdata);

  // int n_up = 0;
  // int n_dw = 0;

  // for (size_t i = 0; i < n_crtdata; ++i){

  //   auto crt_data = crt_data_v[i];

  //   _crt_channel[i] = crt_data->Channel();
  //   _crt_t0[i] = crt_data->T0();
  //   _crt_t1[i] = crt_data->T1();
  //   _crt_adc[i]= crt_data->ADC();

  //   if (_crt_channel[i] <= 191) {
  //     n_up++;
  //   } else {
  //     n_dw++;
  //   }

  //   auto ide_v = crt_data_to_ides.at(crt_data.key());
  //   auto ide = ide_v[0]; // TODO: saving only the first IDE, may be enough?
  //   _crt_true_nide[i] = ide_v.size();
  //   _crt_true_t[i] = 0.5 * (ide->entryT + ide->exitT);
  //   _crt_true_e[i] = ide->energyDeposited;
  //   _crt_true_x[i] = 0.5 * (ide->entryX + ide->exitX);
  //   _crt_true_y[i] = 0.5 * (ide->entryY + ide->exitY);
  //   _crt_true_z[i] = 0.5 * (ide->entryZ + ide->exitZ);
  //   _crt_true_exitp[i] = std::sqrt(ide->exitMomentumX*ide->exitMomentumX +
  //                                   ide->exitMomentumY*ide->exitMomentumY +
  //                                   ide->exitMomentumZ*ide->exitMomentumZ);
  //   // Get the MCParticle that caused this CRTData
  //   auto mcp = trackid_to_mcp[ide->trackID];
  //   _crt_true_mcp_pdg[i] = mcp.PdgCode();
  //   _crt_true_mcp_isprimary[i] = mcp.Process() == "primary";
  //   _crt_true_mcp_process[i] = mcp.Process();
  //   _crt_true_mcp_mother[i] = mcp.Mother();
  //   _crt_true_mcp_mother_pdg[i] = trackid_to_mcp[mcp.Mother()].PdgCode();;
  //   if (mcp.NumberTrajectoryPoints()) {
  //     _crt_true_mcp_time[i] = mcp.T();
  //     _crt_true_mcp_k_start[i] = mcp.E() - mcp.Mass();
  //     _crt_true_mcp_k_end[i] = mcp.E(mcp.Trajectory().size() - 1) - mcp.Mass();
  //     _crt_true_mcp_vx[i] = mcp.Vx();
  //     _crt_true_mcp_vy[i] = mcp.Vy();
  //     _crt_true_mcp_vz[i] = mcp.Vz();
  //   } else {
  //     std::cout << "Don't have trajectory points for track ID " << ide->trackID << std::endl;
  //     _crt_true_mcp_time[i] = -9999;
  //     _crt_true_mcp_k_start[i] = -9999;
  //     _crt_true_mcp_k_end[i] = -9999;
  //     _crt_true_mcp_vx[i] = -9999;
  //     _crt_true_mcp_vy[i] = -9999;
  //     _crt_true_mcp_vz[i] = -9999;
  //   }
  // }

  // _crt_n_up = n_up;
  // _crt_n_dw = n_dw;


  //
  // Fill the Tree
  //
  _tree->Fill();
}


void CRTAnalysis::beginSubRun(art::SubRun const& sr) {

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;

  sr.getByLabel(_pot_label, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
    _sr_spills = pot_handle->totspills;
  } else {
    _sr_pot = 0.;
    _sr_spills = 0.;
  }

  std::cout << "POT for this subrun: " << _sr_pot << " (" << _sr_spills << " spills)" << std::endl;

  _sr_tree->Fill();

}

DEFINE_ART_MODULE(CRTAnalysis)
