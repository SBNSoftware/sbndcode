////////////////////////////////////////////////////////////////////////
// Class:       CRTAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTAnalysis_module.cc
//
// Generated at Tue Feb  1 17:05:06 2022 by Marco Del Tutto using cetskelgen
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
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

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

  std::map<int, int> TrackIDToPrimaryIDMapMaker(art::Event const& e, std::vector<art::Ptr<simb::MCParticle> > &mcp_v);

private:

  sbnd::CRTGeoAlg _crt_geo;

  std::string _mctruth_label;
  std::string _g4_label;
  std::string _auxdethit_label;
  std::string _crtdata_label;
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
  float _mct_sp_e; ///< Single particle energy
  float _mct_sp_px; ///< Single particle momentum X
  float _mct_sp_py; ///< Single particle momentum Y
  float _mct_sp_pz; ///< Single particle momentum Z
  float _mct_sp_vx; ///< Single particle vertex X
  float _mct_sp_vy; ///< Single particle vertex Y
  float _mct_sp_vz; ///< Single particle vertex Z

  float _mct_sp_2_pdg; ///< For the double-muon studys, the second particle PDG
  float _mct_sp_2_e; ///< For the double-muon studys, the second particle energy
  float _mct_sp_2_px; ///< For the double-muon studys, the second particle momentum X
  float _mct_sp_2_py; ///< For the double-muon studys, the second particle momentum Y
  float _mct_sp_2_pz; ///< For the double-muon studys, the second particle momentum Z
  float _mct_sp_2_vx; ///< For the double-muon studys, the second particle vertex X
  float _mct_sp_2_vy; ///< For the double-muon studys, the second particle vertex Y
  float _mct_sp_2_vz; ///< For the double-muon studys, the second particle vertex Z

  float _mct_darkNeutrino_e;
  float _weight; //the weight which store as the vertex of dark neutrino.

  std::vector<int> _mcp_pdg; ///< G4 MCParticle PDG
  std::vector<double> _mcp_e; ///< G4 MCParticle Energy
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
  std::vector<int> _mcp_primaryid; ///< G4 MCParticle, primary ID
  std::vector<int> _mcp_makes_adh; ///< G4 MCParticle, true if this MCP deposits energy in CRT
  std::vector<int> _mcp_primary_makes_adh; ///< G4 MCParticle, true if this MCP's primary deposits energy in CRT
  std::vector<double> _mcp_upstream_dist; ///< Distance from the centre of the upstream tagger at the point the particle would intersect that plane
  std::vector<double> _mcp_downstream_dist; ///< Distance from the centre of the downstream tagger at the point the particle would intersect that plane
  std::vector<int> _mcp_intersects_upstream; ///< Does the extrapolation of the particle's trajectory intercept the upstream tagger (with a buffer zone)
  std::vector<int> _mcp_intersects_downstream; ///< Does the extrapolation of the particle's trajectory intercept the downstream tagger (with a buffer zone)
  std::vector<double> _mcp_cl_ap_upstream; ///< How close does the extrapolation of the particle's trajectory come to the upstream tagger?
  std::vector<double> _mcp_cl_ap_downstream; ///< How close does the extrapolation of the particle's trajectory come to the downstream tagger?
  std::vector<int> _mcp_intersects_upstream1000; ///< Does the extrapolation of the particle's trajectory intercept the upstream tagger (with a buffer zone)
  std::vector<int> _mcp_intersects_downstream1000; ///< Does the extrapolation of the particle's trajectory intercept the downstream tagger (with a buffer zone)

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
  std::vector<std::vector<int> > _chit_true_mcp_trackids; ///< CRT hit contributing true trackIDs
  std::vector<std::vector<int> > _chit_true_mcp_pdg; ///< CRT hit contributing PDG codes
  std::vector<std::vector<double> > _chit_true_mcp_e; ///< CRT hit contributing MCP's Energy
  std::vector<std::vector<double> > _chit_true_mcp_px; ///< CRT hit contributing MCP's Momentum along X
  std::vector<std::vector<double> > _chit_true_mcp_py; ///< CRT hit contributing MCP's Momentum along Y
  std::vector<std::vector<double> > _chit_true_mcp_pz; ///< CRT hit contributing MCP's Momentum along Z
  std::vector<std::vector<double> > _chit_true_mcp_startx; ///< CRT hit contributing MCP's start point X
  std::vector<std::vector<double> > _chit_true_mcp_starty; ///< CRT hit contributing MCP's start point Y
  std::vector<std::vector<double> > _chit_true_mcp_startz; ///< CRT hit contributing MCP's start point Z
  std::vector<std::vector<double> > _chit_true_mcp_endx; ///< CRT hit contributing MCP's end point X
  std::vector<std::vector<double> > _chit_true_mcp_endy; ///< CRT hit contributing MCP's end point Y
  std::vector<std::vector<double> > _chit_true_mcp_endz; ///< CRT hit contributing MCP's end point Z
  std::vector<std::vector<int> > _chit_true_mcp_isprimary; ///< CRT hit contributing MCP's, true if primary
  std::vector<std::vector<uint16_t> > _chit_sipm_adc; ///< The 4 contributing ADC values to this CRTHit
  std::vector<std::vector<uint16_t> > _chit_sipm_channel_id; ///< The IDs of the four SiPMs that were used to make this CRTHit
  std::vector<std::vector<uint16_t> > _chit_sipm_feb_mac5; ///< The IDs of the two FEBs that were used in the making of this hit.

  std::vector<double> _ct_time; ///< CRT track time
  std::vector<double> _ct_pes; ///< CRT track PEs
  std::vector<double> _ct_length; ///< CRT track length
  std::vector<double> _ct_tof; ///< CRT track time of flight
  std::vector<double> _ct_true_tof;  ///< CRT track time of flight, true
  std::vector<double> _ct_x1; ///< CRT track x1
  std::vector<double> _ct_y1; ///< CRT track y1
  std::vector<double> _ct_z1; ///< CRT track z1
  std::vector<double> _ct_x2; ///< CRT track x2
  std::vector<double> _ct_y2; ///< CRT track y2
  std::vector<double> _ct_z2; ///< CRT track z2

  std::vector<uint16_t> _feb_mac5; ///< FEBData Mac5 ID
  std::vector<uint16_t> _feb_flags; ///< FEBData Flags
  std::vector<uint32_t> _feb_ts0; ///< FEBData Ts0
  std::vector<uint32_t> _feb_ts1; ///< FEBData Ts1
  std::vector<uint32_t> _feb_unixs; ///< FEBData UnixS
  std::vector<std::vector<uint16_t>> _feb_adc; ///< FEBData 32 ADC values
  std::vector<uint32_t> _feb_coinc; ///< FEBData Coinc

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
  _mctruth_label = p.get<std::string>("MCTruthLabel", "generator");
  _g4_label = p.get<std::string>("G4Label", "largeant");
  _auxdethit_label = p.get<std::string>("AuxDetHitLabel", "largeant:LArG4DetectorServicevolAuxDetSensitiveCRTStripBERN");
  _crthit_label = p.get<std::string>("CRTHitLabel", "crthit");
  _crttrack_label = p.get<std::string>("CRTTrackLabel", "crttrack");
  _febdata_label = p.get<std::string>("FEBDataLabel", "crtsim");
  _debug = p.get<bool>("Debug", false);
  _keep_mcp_from_adh = p.get<bool>("KeepMCPFromADH", false);
  _data_mode = p.get<bool>("DataMode", false);
  _pot_label = p.get<std::string>("POTLabel", "generator");


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
      _tree->Branch("nu_py", &_nu_py, "nu_px/F");
      _tree->Branch("nu_pz", &_nu_pz, "nu_px/F");

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

      _tree->Branch("mcp_pdg", "std::vector<int>", &_mcp_pdg);
      _tree->Branch("mcp_e", "std::vector<double>", &_mcp_e);
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
      _tree->Branch("mcp_primaryid", "std::vector<int>", &_mcp_primaryid);
      _tree->Branch("mcp_makes_adh", "std::vector<int>", &_mcp_makes_adh);
      _tree->Branch("mcp_primary_makes_adh", "std::vector<int>", &_mcp_primary_makes_adh);
      _tree->Branch("mcp_upstream_dist", "std::vector<double>", &_mcp_upstream_dist);
      _tree->Branch("mcp_downstream_dist", "std::vector<double>", &_mcp_downstream_dist);
      _tree->Branch("mcp_intercepts_upstream", "std::vector<int>", &_mcp_intersects_upstream);
      _tree->Branch("mcp_intercepts_downstream", "std::vector<int>", &_mcp_intersects_downstream);
      _tree->Branch("mcp_cl_ap_upstream", "std::vector<double>" , &_mcp_cl_ap_upstream);
      _tree->Branch("mcp_cl_ap_downstream", "std::vector<double>" , &_mcp_cl_ap_downstream);
      _tree->Branch("mcp_intercepts_upstream1000", "std::vector<int>", &_mcp_intersects_upstream1000);
      _tree->Branch("mcp_intercepts_downstream1000", "std::vector<int>", &_mcp_intersects_downstream1000);

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
      _tree->Branch("chit_true_mcp_trackids", "std::vector<std::vector<int> >", &_chit_true_mcp_trackids);
      _tree->Branch("chit_true_mcp_pdg", "std::vector<std::vector<int> >", &_chit_true_mcp_pdg);
      _tree->Branch("chit_true_mcp_e", "std::vector<std::vector<double> >", &_chit_true_mcp_e);
      _tree->Branch("chit_true_mcp_px", "std::vector<std::vector<double> >", &_chit_true_mcp_px);
      _tree->Branch("chit_true_mcp_py", "std::vector<std::vector<double> >", &_chit_true_mcp_py);
      _tree->Branch("chit_true_mcp_pz", "std::vector<std::vector<double> >", &_chit_true_mcp_pz);
      _tree->Branch("chit_true_mcp_startx", "std::vector<std::vector<double> >", &_chit_true_mcp_startx);
      _tree->Branch("chit_true_mcp_starty", "std::vector<std::vector<double> >", &_chit_true_mcp_starty);
      _tree->Branch("chit_true_mcp_startz", "std::vector<std::vector<double> >", &_chit_true_mcp_startz);
      _tree->Branch("chit_true_mcp_endx", "std::vector<std::vector<double> >", &_chit_true_mcp_endx);
      _tree->Branch("chit_true_mcp_endy", "std::vector<std::vector<double> >", &_chit_true_mcp_endy);
      _tree->Branch("chit_true_mcp_endz", "std::vector<std::vector<double> >", &_chit_true_mcp_endz);
      _tree->Branch("chit_true_mcp_isprimary", "std::vector<std::vector<int> >", &_chit_true_mcp_isprimary);
    }
  _tree->Branch("chit_sipm_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_adc);
  _tree->Branch("chit_sipm_channel_id", "std::vector<std::vector<uint16_t> >", &_chit_sipm_channel_id);
  _tree->Branch("chit_sipm_feb_mac5", "std::vector<std::vector<uint16_t> >", &_chit_sipm_feb_mac5);

  _tree->Branch("ct_time", "std::vector<double>", &_ct_time);
  _tree->Branch("ct_pes", "std::vector<double>", &_ct_pes);
  _tree->Branch("ct_length", "std::vector<double>", &_ct_length);
  _tree->Branch("ct_tof", "std::vector<double>", &_ct_tof);
  if(!_data_mode) _tree->Branch("ct_true_tof", "std::vector<double>", &_ct_true_tof);
  _tree->Branch("ct_x1", "std::vector<double>", &_ct_x1);
  _tree->Branch("ct_y1", "std::vector<double>", &_ct_y1);
  _tree->Branch("ct_z1", "std::vector<double>", &_ct_z1);
  _tree->Branch("ct_x2", "std::vector<double>", &_ct_x2);
  _tree->Branch("ct_y2", "std::vector<double>", &_ct_y2);
  _tree->Branch("ct_z2", "std::vector<double>", &_ct_z2);

  _tree->Branch("feb_mac5", "std::vector<uint16_t>", &_feb_mac5);
  _tree->Branch("feb_flags", "std::vector<uint16_t>", &_feb_flags);
  _tree->Branch("feb_ts0", "std::vector<uint32_t>", &_feb_ts0);
  _tree->Branch("feb_ts1", "std::vector<uint32_t>", &_feb_ts1);
  _tree->Branch("feb_unixs", "std::vector<uint32_t>", &_feb_unixs);
  _tree->Branch("feb_adc", "std::vector<std::vector<uint16_t>>", &_feb_adc);
  _tree->Branch("feb_coinc", "std::vector<uint32_t>", &_feb_coinc);

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
  _sr_tree->Branch("spills", &_sr_spills, "spills/D");
}

void CRTAnalysis::analyze(art::Event const& e)
{
  if(_debug)
    {
      for(auto const &[name, tagger] : _crt_geo.GetTaggers())
	{
	  std::cout << "Tagger:  " << tagger.name << '\n'
		    << "X - Min: " << tagger.minX << " Max: " << tagger.maxX << '\n'
		    << "Y - Min: " << tagger.minY << " Max: " << tagger.maxY << '\n'
		    << "Z - Min: " << tagger.minZ << " Max: " << tagger.maxZ << '\n' << std::endl;
	}
    }
  
  
  std::vector<double> upstream_limits = {999999., -999999., 999999., -999999., 999999., -999999.};
  auto north_tagger = _crt_geo.GetTagger("volTaggerNorth_0");
  std::vector<double> downstream_limits = {north_tagger.minX, north_tagger.maxX, north_tagger.minY, north_tagger.maxY, north_tagger.minZ, north_tagger.maxZ};

  for(std::string name : {"volTaggerSouthOne_0", "volTaggerSouthTwo_0", "volTaggerSouthThree_0"})
    {
      auto tagger = _crt_geo.GetTagger(name);
      if(tagger.minX < upstream_limits[0]) upstream_limits[0] = tagger.minX;
      if(tagger.maxX > upstream_limits[1]) upstream_limits[1] = tagger.maxX;
      if(tagger.minY < upstream_limits[2]) upstream_limits[2] = tagger.minY;
      if(tagger.maxY > upstream_limits[3]) upstream_limits[3] = tagger.maxY;
      if(tagger.minZ < upstream_limits[4]) upstream_limits[4] = tagger.minZ;
      if(tagger.maxZ > upstream_limits[5]) upstream_limits[5] = tagger.maxZ;
    }

  if(_debug)
    {
      std::cout << "- Upstream Tagger -\n" 
		<< "X - Min: " << upstream_limits[0] << " Max: " << upstream_limits[1] << '\n'
		<< "Y - Min: " << upstream_limits[2] << " Max: " << upstream_limits[3] << '\n'
		<< "Z - Min: " << upstream_limits[4] << " Max: " << upstream_limits[5] << '\n' << std::endl;
      
      std::cout << "- Downstream Tagger -\n" 
		<< "X - Min: " << downstream_limits[0] << " Max: " << downstream_limits[1] << '\n'
		<< "Y - Min: " << downstream_limits[2] << " Max: " << downstream_limits[3] << '\n'
		<< "Z - Min: " << downstream_limits[4] << " Max: " << downstream_limits[5] << '\n' << std::endl;
    }

  for(unsigned i = 0; i < 6; ++i)
    {
      if(i%2 == 0) {
	upstream_limits[i] -= 200.; downstream_limits[i] -= 200.;
      }
      if(i%2 == 1) {
	upstream_limits[i] += 200.; downstream_limits[i] += 200.;
      }
    }
  
  if(_debug)
    {
      std::cout << "- Upstream Tagger -\n" 
		<< "X - Min: " << upstream_limits[0] << " Max: " << upstream_limits[1] << '\n'
		<< "Y - Min: " << upstream_limits[2] << " Max: " << upstream_limits[3] << '\n'
		<< "Z - Min: " << upstream_limits[4] << " Max: " << upstream_limits[5] << '\n' << std::endl;
      
      std::cout << "- Downstream Tagger -\n" 
		<< "X - Min: " << downstream_limits[0] << " Max: " << downstream_limits[1] << '\n'
		<< "Y - Min: " << downstream_limits[2] << " Max: " << downstream_limits[3] << '\n'
		<< "Z - Min: " << downstream_limits[4] << " Max: " << downstream_limits[5] << '\n' << std::endl;
    }

  const double upstream_x   = 0.5 * (upstream_limits[0] + upstream_limits[1]);
  const double upstream_y   = 0.5 * (upstream_limits[2] + upstream_limits[3]);
  const double upstream_z   = 0.5 * (upstream_limits[4] + upstream_limits[5]);
  const double downstream_x = 0.5 * (downstream_limits[0] + downstream_limits[1]);
  const double downstream_y = 0.5 * (downstream_limits[2] + downstream_limits[3]);
  const double downstream_z = 0.5 * (downstream_limits[4] + downstream_limits[5]);

  const TVector3 upstream_centre   = TVector3(upstream_x, upstream_y, upstream_z);
  const TVector3 downstream_centre = TVector3(downstream_x, downstream_y, downstream_z);

  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(_debug) std::cout << "This is event " << _event << std::endl;

  art::Handle<std::vector<simb::MCTruth>> mct_h;
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  std::map<int, simb::MCParticle> trackid_to_mcp;
  art::Handle<std::vector<sim::AuxDetHit>> adh_h;
  std::vector<art::Ptr<sim::AuxDetHit>> adh_v;
  
  if(!_data_mode)
    {
      //
      // Get the MCTruth
      //
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
  /*  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> *feb_data_to_ides;
  if(!_data_mode) {
    feb_data_to_ides = new art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo>(feb_data_h, e, _crtdata_label);
  }
  */
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

  if(!_data_mode) {
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
      assert(mct->NParticles() <= 3); // die if mct->NPartickes() != 1
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

    std::map<int, int> trackid_to_primaryid_map = TrackIDToPrimaryIDMapMaker(e, mcp_v);
    if(trackid_to_primaryid_map.size() != mcp_v.size())
      {
	std::cout << "TrackIDToPrimaryIDMap is wrong size (" << trackid_to_primaryid_map.size() << " vs. " << mcp_v.size() << ")" << std::endl;
	//      throw std::exception();
      }


    //
    // Fill the AuxDetHits in the tree
    //
    size_t n_adh = adh_v.size();
    if (_debug) std::cout << "n_adh " << n_adh << std::endl;
    _adh_t.resize(n_adh);
    _adh_e.resize(n_adh);
    _adh_x.resize(n_adh);
    _adh_y.resize(n_adh);
    _adh_z.resize(n_adh);
    _adh_p_exit.resize(n_adh);
    _adh_trackid.resize(n_adh);

    std::set<unsigned int> trackids_from_adh;
    std::set<unsigned int> primaryids_from_adh;

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
      primaryids_from_adh.insert(trackid_to_primaryid_map[auxdethit->GetTrackID()]);
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
    _mcp_primaryid.resize(n_mcp);
    _mcp_makes_adh.resize(n_mcp);
    _mcp_primary_makes_adh.resize(n_mcp);
    _mcp_upstream_dist.resize(n_mcp);
    _mcp_downstream_dist.resize(n_mcp);
    _mcp_intersects_upstream.resize(n_mcp);
    _mcp_intersects_downstream.resize(n_mcp);
    _mcp_cl_ap_upstream.resize(n_mcp);
    _mcp_cl_ap_downstream.resize(n_mcp);
    _mcp_intersects_upstream1000.resize(n_mcp);
    _mcp_intersects_downstream1000.resize(n_mcp);

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
      _mcp_primaryid[counter] = trackid_to_primaryid_map[particle->TrackId()];
      _mcp_makes_adh[counter] = trackids_from_adh.find(particle->TrackId()) != trackids_from_adh.end();
      if(_mcp_isprimary[counter])
	_mcp_primary_makes_adh[counter] = primaryids_from_adh.find(particle->TrackId()) != primaryids_from_adh.end();
      else
	_mcp_primary_makes_adh[counter] = -1;

      const double t_upstream   = (upstream_z -_mcp_startz[counter]) / _mcp_pz[counter];
      const double t_downstream = (downstream_z -_mcp_startz[counter]) / _mcp_pz[counter];

      const double mcp_upstream_x = _mcp_startx[counter] + _mcp_px[counter] * t_upstream;
      const double mcp_upstream_y = _mcp_starty[counter] + _mcp_py[counter] * t_upstream;

      const double mcp_downstream_x = _mcp_startx[counter] + _mcp_px[counter] * t_downstream;
      const double mcp_downstream_y = _mcp_starty[counter] + _mcp_py[counter] * t_downstream;

      const TVector3 mcp_upstream   = TVector3(mcp_upstream_x, mcp_upstream_y, upstream_z);
      const TVector3 mcp_downstream = TVector3(mcp_downstream_x, mcp_downstream_y, downstream_z);

      _mcp_upstream_dist[counter]   = (mcp_upstream - upstream_centre).Mag();
      _mcp_downstream_dist[counter] = (mcp_downstream - downstream_centre).Mag();

      std::vector<double> upstream_wide_limits = upstream_limits;
      std::vector<double> downstream_wide_limits = downstream_limits;
      std::vector<double> buffer_box = {-0., 0., -0., 0., -0., 0};

      for (unsigned j = 0; j < 6; ++j) {
	upstream_wide_limits[j] += buffer_box[j];
	downstream_wide_limits[j] += buffer_box[j];
      }

      const std::vector<double> x0 = {_mcp_startx[counter], _mcp_starty[counter], _mcp_startz[counter]};
      const std::vector<double> dx = {_mcp_px[counter], _mcp_py[counter], _mcp_pz[counter]};

      bool intersects_upstream = false;
      double closest_approach_upstream = std::numeric_limits<double>::max();

      for (int bnd = 0; bnd != 6; ++bnd) {
	if (bnd<2) {
	  double p2[3] = {upstream_wide_limits[bnd],  x0[1] + (dx[1]/dx[0])*(upstream_wide_limits[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(upstream_wide_limits[bnd] - x0[0])};
	  if ( p2[1] >= upstream_wide_limits[2] && p2[1] <= upstream_wide_limits[3] &&
	       p2[2] >= upstream_wide_limits[4] && p2[2] <= upstream_wide_limits[5] ) {
	    intersects_upstream = true;
	    closest_approach_upstream = 0.;
	    break;
	  }
	  else if (p2[1] >= upstream_wide_limits[2] && p2[1] <= upstream_wide_limits[3])
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[4] - p2[2]), std::abs(upstream_wide_limits[5] - p2[2])});
	  else
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[2] - p2[1]), std::abs(upstream_wide_limits[3] - p2[1])});
	}
	else if (bnd>=2 && bnd<4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[1])*(upstream_wide_limits[bnd] - x0[1]), upstream_wide_limits[bnd], x0[2] + (dx[2]/dx[1])*(upstream_wide_limits[bnd] - x0[1])};
	  if ( p2[0] >= upstream_wide_limits[0] && p2[0] <= upstream_wide_limits[1] &&
	       p2[2] >= upstream_wide_limits[4] && p2[2] <= upstream_wide_limits[5] ) {
	    intersects_upstream = true;
	    closest_approach_upstream = 0.;
	    break;
	  }
	  else if (p2[0] >= upstream_wide_limits[0] && p2[0] <= upstream_wide_limits[1])
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[4] - p2[2]), std::abs(upstream_wide_limits[5] - p2[2])});
	  else
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[0] - p2[0]), std::abs(upstream_wide_limits[1] - p2[0])});
	}
	else if (bnd>=4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[2])*(upstream_wide_limits[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(upstream_wide_limits[bnd] - x0[2]), upstream_wide_limits[bnd]};
	  if ( p2[0] >= upstream_wide_limits[0] && p2[0] <= upstream_wide_limits[1] &&
	       p2[1] >= upstream_wide_limits[2] && p2[1] <= upstream_wide_limits[3] ) {
	    intersects_upstream = true;
	    closest_approach_upstream = 0.;
	    break;
	  }
	  else if (p2[0] >= upstream_wide_limits[0] && p2[0] <= upstream_wide_limits[1])
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[2] - p2[1]), std::abs(upstream_wide_limits[3] - p2[1])});
	  else
	    closest_approach_upstream = std::min({closest_approach_upstream, std::abs(upstream_wide_limits[0] - p2[0]), std::abs(upstream_wide_limits[1] - p2[0])});
	}
      }

      bool intersects_downstream = false;
      double closest_approach_downstream = std::numeric_limits<double>::max();

      for (int bnd = 0; bnd != 6; ++bnd) {
	if (bnd<2) {
	  double p2[3] = {downstream_wide_limits[bnd],  x0[1] + (dx[1]/dx[0])*(downstream_wide_limits[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(downstream_wide_limits[bnd] - x0[0])};
	  if ( p2[1] >= downstream_wide_limits[2] && p2[1] <= downstream_wide_limits[3] &&
	       p2[2] >= downstream_wide_limits[4] && p2[2] <= downstream_wide_limits[5] ) {
	    intersects_downstream = true;
	    closest_approach_downstream = 0.;
	    break;
	  }
	  else if (p2[1] >= downstream_wide_limits[2] && p2[1] <= downstream_wide_limits[3])
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[4] - p2[2]), std::abs(downstream_wide_limits[5] - p2[2])});
	  else
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[2] - p2[1]), std::abs(downstream_wide_limits[3] - p2[1])});
	}
	else if (bnd>=2 && bnd<4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[1])*(downstream_wide_limits[bnd] - x0[1]), downstream_wide_limits[bnd], x0[2] + (dx[2]/dx[1])*(downstream_wide_limits[bnd] - x0[1])};
	  if ( p2[0] >= downstream_wide_limits[0] && p2[0] <= downstream_wide_limits[1] &&
	       p2[2] >= downstream_wide_limits[4] && p2[2] <= downstream_wide_limits[5] ) {
	    intersects_downstream = true;
	    closest_approach_downstream = 0.;
	    break;
	  }
	  else if (p2[0] >= downstream_wide_limits[0] && p2[0] <= downstream_wide_limits[1])
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[4] - p2[2]), std::abs(downstream_wide_limits[5] - p2[2])});
	  else
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[0] - p2[0]), std::abs(downstream_wide_limits[1] - p2[0])});
	}
	else if (bnd>=4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[2])*(downstream_wide_limits[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(downstream_wide_limits[bnd] - x0[2]), downstream_wide_limits[bnd]};
	  if ( p2[0] >= downstream_wide_limits[0] && p2[0] <= downstream_wide_limits[1] &&
	       p2[1] >= downstream_wide_limits[2] && p2[1] <= downstream_wide_limits[3] ) {
	    intersects_downstream = true;
	    closest_approach_downstream = 0.;
	    break;
	  }
	  else if (p2[0] >= downstream_wide_limits[0] && p2[0] <= downstream_wide_limits[1])
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[2] - p2[1]), std::abs(downstream_wide_limits[3] - p2[1])});
	  else
	    closest_approach_downstream = std::min({closest_approach_downstream, std::abs(downstream_wide_limits[0] - p2[0]), std::abs(downstream_wide_limits[1] - p2[0])});
	}
      }

      _mcp_intersects_upstream[counter] = intersects_upstream;
      _mcp_intersects_downstream[counter] = intersects_downstream;
      _mcp_cl_ap_upstream[counter] = closest_approach_upstream;
      _mcp_cl_ap_downstream[counter] = closest_approach_downstream;

      std::vector<double> upstream_wide_limits1000 = upstream_limits;
      std::vector<double> downstream_wide_limits1000 = downstream_limits;
      std::vector<double> buffer_box1000 = {-1000., 1000., -1000., 1000., -1000., 1000.};

      for (unsigned j = 0; j < 6; ++j) {
	upstream_wide_limits1000[j] += buffer_box1000[j];
	downstream_wide_limits1000[j] += buffer_box1000[j];
      }

      bool intersects_upstream1000 = false;

      for (int bnd = 0; bnd != 6; ++bnd) {
	if (bnd<2) {
	  double p2[3] = {upstream_wide_limits1000[bnd],  x0[1] + (dx[1]/dx[0])*(upstream_wide_limits1000[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(upstream_wide_limits1000[bnd] - x0[0])};
	  if ( p2[1] >= upstream_wide_limits1000[2] && p2[1] <= upstream_wide_limits1000[3] &&
	       p2[2] >= upstream_wide_limits1000[4] && p2[2] <= upstream_wide_limits1000[5] ) {
	    intersects_upstream1000 = true;
	    break;
	  }
	}
	else if (bnd>=2 && bnd<4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[1])*(upstream_wide_limits1000[bnd] - x0[1]), upstream_wide_limits1000[bnd], x0[2] + (dx[2]/dx[1])*(upstream_wide_limits1000[bnd] - x0[1])};
	  if ( p2[0] >= upstream_wide_limits1000[0] && p2[0] <= upstream_wide_limits1000[1] &&
	       p2[2] >= upstream_wide_limits1000[4] && p2[2] <= upstream_wide_limits1000[5] ) {
	    intersects_upstream1000 = true;
	    break;
	  }
	}
	else if (bnd>=4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[2])*(upstream_wide_limits1000[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(upstream_wide_limits1000[bnd] - x0[2]), upstream_wide_limits1000[bnd]};
	  if ( p2[0] >= upstream_wide_limits1000[0] && p2[0] <= upstream_wide_limits1000[1] &&
	       p2[1] >= upstream_wide_limits1000[2] && p2[1] <= upstream_wide_limits1000[3] ) {
	    intersects_upstream1000 = true;
	    break;
	  }
	}
      }

      bool intersects_downstream1000 = false;

      for (int bnd = 0; bnd != 6; ++bnd) {
	if (bnd<2) {
	  double p2[3] = {downstream_wide_limits1000[bnd],  x0[1] + (dx[1]/dx[0])*(downstream_wide_limits1000[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(downstream_wide_limits1000[bnd] - x0[0])};
	  if ( p2[1] >= downstream_wide_limits1000[2] && p2[1] <= downstream_wide_limits1000[3] &&
	       p2[2] >= downstream_wide_limits1000[4] && p2[2] <= downstream_wide_limits1000[5] ) {
	    intersects_downstream1000 = true;
	    break;
	  }
	}
	else if (bnd>=2 && bnd<4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[1])*(downstream_wide_limits1000[bnd] - x0[1]), downstream_wide_limits1000[bnd], x0[2] + (dx[2]/dx[1])*(downstream_wide_limits1000[bnd] - x0[1])};
	  if ( p2[0] >= downstream_wide_limits1000[0] && p2[0] <= downstream_wide_limits1000[1] &&
	       p2[2] >= downstream_wide_limits1000[4] && p2[2] <= downstream_wide_limits1000[5] ) {
	    intersects_downstream1000 = true;
	    break;
	  }
	}
	else if (bnd>=4) {
	  double p2[3] = {x0[0] + (dx[0]/dx[2])*(downstream_wide_limits1000[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(downstream_wide_limits1000[bnd] - x0[2]), downstream_wide_limits1000[bnd]};
	  if ( p2[0] >= downstream_wide_limits1000[0] && p2[0] <= downstream_wide_limits1000[1] &&
	       p2[1] >= downstream_wide_limits1000[2] && p2[1] <= downstream_wide_limits1000[3] ) {
	    intersects_downstream1000 = true;
	    break;
	  }
	}
      }

      _mcp_intersects_upstream1000[counter] = intersects_upstream1000;
      _mcp_intersects_downstream1000[counter] = intersects_downstream1000;

      if(_debug)
	{
	  std::cout << "MCP - " << counter << '\n'
		    << "Start:  " << _mcp_startx[counter] << ", " << _mcp_starty[counter] << ", " << _mcp_startz[counter] << '\n'
		    << "Mom:    " << _mcp_px[counter] << ", " << _mcp_py[counter] << ", " << _mcp_pz[counter] << '\n'
		    << "t_up:   " << t_upstream << '\n'
		    << "t_down: " << t_downstream << '\n' 
		    << "UP:     " << mcp_upstream.x() << ", " << mcp_upstream.y() << ", " << mcp_upstream.z() << '\n'
		    << "DOWN:   " << mcp_downstream.x() << ", " << mcp_downstream.y() << ", " << mcp_downstream.z() << '\n'
		    << "DIST-U: " << _mcp_upstream_dist[counter] << '\n'
		    << "DIST-D: " << _mcp_downstream_dist[counter] << '\n'
		    << std::endl;
	}


      if (_debug) {
	std::cout << "MCP " << _mcp_pdg[counter] << ", trackID = " << _mcp_trackid[counter]
		  << ", p = " << particle->P()
		  << "(" << particle->Px() << "," << particle->Py() << "," << particle->Pz() << ")"
		  << ", process = "
		  << particle->Process() << ": start point = ("
		  << _mcp_startx[counter] << ", " << _mcp_starty[counter] << ", " << _mcp_startz[counter]
		  << ") - start process: " << particle->Process()
		  << " --- end point = (" << std::endl;
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
  _chit_true_mcp_trackids.resize(n_hits);
  _chit_true_mcp_pdg.resize(n_hits);
  _chit_true_mcp_e.resize(n_hits);
  _chit_true_mcp_px.resize(n_hits);
  _chit_true_mcp_py.resize(n_hits);
  _chit_true_mcp_pz.resize(n_hits);
  _chit_true_mcp_startx.resize(n_hits);
  _chit_true_mcp_starty.resize(n_hits);
  _chit_true_mcp_startz.resize(n_hits);
  _chit_true_mcp_endx.resize(n_hits);
  _chit_true_mcp_endy.resize(n_hits);
  _chit_true_mcp_endz.resize(n_hits);
  _chit_true_mcp_px.resize(n_hits);
  _chit_true_mcp_py.resize(n_hits);
  _chit_true_mcp_isprimary.resize(n_hits);
  _chit_sipm_adc.resize(n_hits);
  _chit_sipm_channel_id.resize(n_hits);
  _chit_sipm_feb_mac5.resize(n_hits);

  for (size_t i = 0; i < n_hits; i++) {

    auto hit = crt_hit_v[i];

    _chit_x[i] = hit->x_pos;
    _chit_y[i] = hit->y_pos;
    _chit_z[i] = hit->z_pos;
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

    auto feb_datas = crt_hit_to_feb_data.at(hit.key());
    if(feb_datas.size() != 2) std::cout << "ERROR: CRTHit associated to " << feb_datas.size() << " FEBDatas" << std::endl;

    _chit_sipm_adc[i].resize(2 * feb_datas.size());
    _chit_sipm_feb_mac5[i].resize(feb_datas.size());

    _chit_sipm_channel_id[i].resize(2);
    _chit_sipm_channel_id[i][0] = hit->channel0;
    _chit_sipm_channel_id[i][1] = hit->channel1;

    for(unsigned ii = 0; ii < feb_datas.size(); ++ii)
      {
	const auto& feb_data = feb_datas[ii];
	_chit_sipm_feb_mac5[i][ii] = feb_data->Mac5();
	_chit_sipm_adc[i][2*ii] = 0;
	_chit_sipm_adc[i][2*ii+1] = 0;

	if(_chit_sipm_channel_id[i][0] / 32 == _chit_sipm_feb_mac5[i][ii]) 
	  {
	    const uint otherchannel = (_chit_sipm_channel_id[i][0] % 2 == 0) ? _chit_sipm_channel_id[i][0] + 1 : _chit_sipm_channel_id[i][0] - 1;
	    _chit_sipm_adc[i][2*ii]   = feb_data->ADC()[_chit_sipm_channel_id[i][0] % 32];
	    _chit_sipm_adc[i][2*ii+1] = feb_data->ADC()[otherchannel % 32];
	  }
	if(_chit_sipm_channel_id[i][1] / 32 == _chit_sipm_feb_mac5[i][ii])
	  {
	    const uint otherchannel = (_chit_sipm_channel_id[i][1] % 2 == 0) ? _chit_sipm_channel_id[i][1] + 1 : _chit_sipm_channel_id[i][1] - 1;
	    _chit_sipm_adc[i][2*ii]   = feb_data->ADC()[_chit_sipm_channel_id[i][1] % 32];
	    _chit_sipm_adc[i][2*ii+1] = feb_data->ADC()[otherchannel % 32];
	  }
      }

    /*
    size_t n_ides = 0;
    if(!_data_mode) {
      
      // From the hit, get the associated CRTData,
      // then the associated AuxDetIDE, so we can
      // retrieve the truth info
      _chit_true_t[i] = 0;
      _chit_true_e[i] = 0;
      _chit_true_x[i] = 0;
      _chit_true_y[i] = 0;
      _chit_true_z[i] = 0;

      auto feb_datas = crt_hit_to_feb_data.at(hit.key());
      _chit_h1_t0[i] = crt_data_v[0]->T0();
      _chit_h2_t0[i] = crt_data_v[2]->T0();
      _chit_h1_t1[i] = crt_data_v[0]->T1();
      _chit_h2_t1[i] = crt_data_v[2]->T1();
      
      for(auto const &feb_data : feb_datas)
	{
	  auto ides       = feb_data_to_ides->at(feb_data.key());
	  auto feb_truths = feb_data_to_ides->at(feb_data.key());
	  
	  for(unsigned ii = 0; ii < ides.size(); ++ii)
	    {
	    if(feb_truths[ii]->GetChannel());
	    
	    n_ides += ;
	    }
	}
	
      _chit_true_mcp_trackids[i].resize(n_ides);
      _chit_true_mcp_pdg[i].resize(n_ides);
      _chit_true_mcp_e[i].resize(n_ides);
      _chit_true_mcp_px[i].resize(n_ides);
      _chit_true_mcp_py[i].resize(n_ides);
      _chit_true_mcp_pz[i].resize(n_ides);
      _chit_true_mcp_startx[i].resize(n_ides);
      _chit_true_mcp_starty[i].resize(n_ides);
      _chit_true_mcp_startz[i].resize(n_ides);
      _chit_true_mcp_endx[i].resize(n_ides);
      _chit_true_mcp_endy[i].resize(n_ides);
      _chit_true_mcp_endz[i].resize(n_ides);
      _chit_true_mcp_px[i].resize(n_ides);
      _chit_true_mcp_py[i].resize(n_ides);
      _chit_true_mcp_isprimary[i].resize(n_ides);

      _chit_sipm_adc[i].resize(crt_data_v.size());
      _chit_sipm_channel_id[i].resize(crt_data_v.size());
      _chit_sipm_feb_mac5[i].resize(crt_data_v.size());

      size_t data_counter = 0, ide_counter = 0;
      for (auto crt_data : crt_data_v) {
	auto feb_v = crt_data_to_feb_data->at(crt_data.key());

	if(feb_v.size() != 1) std::cout << "========== ERROR: Found " << feb_v.size() << " FEBDatas for one CRTData?" << std::endl;

	_chit_sipm_adc[i][data_counter] = crt_data->ADC();
	_chit_sipm_channel_id[i][data_counter] = crt_data->Channel();
	_chit_sipm_feb_mac5[i][data_counter] = feb_v[0]->Mac5();

	auto ide_v = crt_data_to_ides->at(crt_data.key());
	for (auto ide : ide_v) {
	  _chit_true_t[i] += 0.5 * (ide->entryT + ide->exitT);
	  _chit_true_e[i] += ide->energyDeposited;
	  _chit_true_x[i] += 0.5 * (ide->entryX + ide->exitX);
	  _chit_true_y[i] += 0.5 * (ide->entryY + ide->exitY);
	  _chit_true_z[i] += 0.5 * (ide->entryZ + ide->exitZ);
	
	  if(_debug)
	    {
	      std::cout << "IDE Track ID: " << ide->trackID << std::endl;
	      if(trackid_to_mcp.find(ide->trackID) != trackid_to_mcp.end())
		std::cout << "Opened map: Success! " << ide->energyDeposited << std::endl;
	      else
		std::cout << "Opened map: Failure! " << ide->energyDeposited << std::endl;
	    }

	  auto mcp = trackid_to_mcp[ide->trackID];
	  _chit_true_mcp_trackids[i][ide_counter] = mcp.TrackId();
	  _chit_true_mcp_pdg[i][ide_counter] = mcp.PdgCode();

	  if(mcp.NumberTrajectoryPoints())
	    {
	      _chit_true_mcp_e[i][ide_counter] = mcp.E();
	      _chit_true_mcp_px[i][ide_counter] = mcp.Px();
	      _chit_true_mcp_py[i][ide_counter] = mcp.Py();
	      _chit_true_mcp_pz[i][ide_counter] = mcp.Pz();
	      _chit_true_mcp_startx[i][ide_counter] = mcp.Vx();
	      _chit_true_mcp_starty[i][ide_counter] = mcp.Vy();
	      _chit_true_mcp_startz[i][ide_counter] = mcp.Vz();
	      _chit_true_mcp_endx[i][ide_counter] = mcp.EndX();
	      _chit_true_mcp_endy[i][ide_counter] = mcp.EndY();
	      _chit_true_mcp_endz[i][ide_counter] = mcp.EndZ();
	    }
	  else
	    {
	      _chit_true_mcp_e[i][ide_counter] = -9999.;
	      _chit_true_mcp_px[i][ide_counter] = -9999.;
	      _chit_true_mcp_py[i][ide_counter] = -9999.;
	      _chit_true_mcp_pz[i][ide_counter] = -9999.;
	      _chit_true_mcp_startx[i][ide_counter] = -9999.;
	      _chit_true_mcp_starty[i][ide_counter] = -9999.;
	      _chit_true_mcp_startz[i][ide_counter] = -9999.;
	      _chit_true_mcp_endx[i][ide_counter] = -9999.;
	      _chit_true_mcp_endy[i][ide_counter] = -9999.;
	      _chit_true_mcp_endz[i][ide_counter] = -9999.;
	    }

	  _chit_true_mcp_isprimary[i][ide_counter] = mcp.Process() == "primary";

	  ++ide_counter;
	}
	++data_counter;
      }
      _chit_true_t[i] /= n_ides;
      _chit_true_e[i] /= n_ides;
      _chit_true_x[i] /= n_ides;
      _chit_true_y[i] /= n_ides;
      _chit_true_z[i] /= n_ides;
    }
    */    
    if (_debug) std::cout << "CRT hit, z = " << _chit_z[i] << ", h1 time " << _chit_h1_t1[i] << ", h2 time " << _chit_h2_t1[i] << ", hit time " << _chit_t1[i] << std::endl;
  }


  //
  // Fill the CRT Tracks in the tree
  //
  size_t n_tracks = crt_track_v.size();

  //std::cout<<crt_track_v.size()<<std::endl;

  _ct_pes.resize(n_tracks);
  _ct_time.resize(n_tracks);
  _ct_length.resize(n_tracks);
  _ct_tof.resize(n_tracks);
  _ct_true_tof.resize(n_tracks);
  _ct_x1.resize(n_tracks);
  _ct_y1.resize(n_tracks);
  _ct_z1.resize(n_tracks);
  _ct_x2.resize(n_tracks);
  _ct_y2.resize(n_tracks);
  _ct_z2.resize(n_tracks);

  for (size_t i = 0; i < n_tracks; ++i){

    auto track = crt_track_v[i];

    _ct_pes[i] = track->peshit;
    _ct_time[i] = track->ts1_ns;
    _ct_length[i] = track->length;
    _ct_tof[i] = track->ts0_ns_h2 - track->ts0_ns_h1;
    _ct_x1[i] = track->x1_pos;
    _ct_y1[i] = track->y1_pos;
    _ct_z1[i] = track->z1_pos;
    _ct_x2[i] = track->x2_pos;
    _ct_y2[i] = track->y2_pos;
    _ct_z2[i] = track->z2_pos;
    if (_debug) std::cout << "CRT track, z1 = " << _ct_z1[i] << ", z2 = " << _ct_z2[i] << ", ts0_ns_h1 = " << track->ts0_ns_h1 << ", ts0_ns_h2 = " << track->ts0_ns_h2 << ", tof = " << _ct_tof[i] << std::endl;

    if(!_data_mode) {
      // From the hit, get the associated CRTData,
      // then the associated AuxDetIDE, so we can
      // retrieve the truth info
      _ct_true_tof[i] = 0;
      // 1. Get the hits
      auto hit_v = crt_track_to_hit.at(track.key());
      assert(hit_v.size == 2); // 2 hits per track
      for (size_t i_hit = 0; i_hit < hit_v.size(); i_hit++) {
	/*	auto hit = hit_v[i_hit];
	float hit_time = 0;
	size_t n_ides = 0;
	  auto ide_v = feb_data_to_ides->at(crt_.key());
	  for (auto ide : ide_v) {
	    hit_time += 0.5 * (ide->entryT + ide->exitT);
	    n_ides++;
	  }
	}
	hit_time /= n_ides;
	if (i_hit == 0) _ct_true_tof[i] -= hit_time;
	if (i_hit == 1) _ct_true_tof[i] += hit_time;
	*/
      }
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
    _sr_pot    = pot_handle->totpot;
    _sr_spills = pot_handle->totspills;
  } else {
    _sr_pot    = 0.;
    _sr_spills = 0.;
  }
  if (_debug) std::cout << "POT for this subrun: " << _sr_pot << " (" << _sr_spills << " spills)" << std::endl;

  _sr_tree->Fill();

}

std::map<int, int> CRTAnalysis::TrackIDToPrimaryIDMapMaker(art::Event const& e, std::vector<art::Ptr<simb::MCParticle> > &mcp_v) 
{
  std::map<int, art::Ptr<simb::MCParticle> > particle_map;

  for(auto const & mcp : mcp_v)
    {
      particle_map[mcp->TrackId()] = mcp;
    }

  std::map<int, int> map;

  for(auto const &mcp : mcp_v)
    {
      if(mcp->StatusCode() != 1) continue;

      if(mcp->Process() == "primary")
	map[mcp->TrackId()] = mcp->TrackId();
      else 
	{
	  bool filled = false;
	  art::Ptr<simb::MCParticle> mother = mcp;
	  
	  while(!filled)
	    {
	      if(map.find(mother->Mother()) != map.end())
		{
		  map[mcp->TrackId()] = map[mother->Mother()];
		  filled = true;
		}
	      mother = particle_map[mother->Mother()];
	    }
	}
    }

  return map;

}
      
DEFINE_ART_MODULE(CRTAnalysis)
