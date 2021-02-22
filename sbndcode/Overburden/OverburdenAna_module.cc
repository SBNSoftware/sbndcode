////////////////////////////////////////////////////////////////////////
// Class:       OverburdenAna
// Plugin Type: producer (art v3_04_00)
// File:        OverburdenAna_module.cc
//
// Generated at Tue Jan 28 20:46:15 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include <memory>
#include <map>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcorealg/CoreUtils/ParticleFilters.h" // util::PositionInVolumeFilter
#include "larcore/Geometry/Geometry.h"

#include "TTree.h"

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

class OverburdenAna;



class OverburdenAna : public art::EDProducer {
public:
  explicit OverburdenAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OverburdenAna(OverburdenAna const&) = delete;
  OverburdenAna(OverburdenAna&&) = delete;
  OverburdenAna& operator=(OverburdenAna const&) = delete;
  OverburdenAna& operator=(OverburdenAna&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  virtual void beginSubRun(art::SubRun & sr) override;

  using Point_t = std::array<double, 3>;

private:

  /// Clear vectors
  void clear_vectors();

  /// Finds the ancestors that was created in the Ovrburden, if any
  int FindMotherInOverburden(simb::MCParticle);

  /// Check if the point is in the detector
  bool InDetector(const double& x, const double& y, const double& z) const;

  /// Check if an MCP passes through the detector
  bool InDetector(const art::Ptr<simb::MCParticle>);

  /// Saves the pi0 info to a separate tree, give the pi0 track id
  void SavePi0ShowerInfo(int pi0_track_id);

  /// Configures and returns a particle filter
  std::unique_ptr<util::PositionInVolumeFilter> CreateParticleVolumeFilter
      (std::set<std::string> const& vol_names) const;

  std::map<unsigned int, simb::MCParticle> _trackid_to_mcparticle;

  std::unique_ptr<util::PositionInVolumeFilter> _part_filter;

  std::string _mctruth_producer = "generator"; // For storing POT information
  std::string _mcparticle_producer = "largeant";
  std::string _mctrack_producer = "mcreco";
  std::string _mcshower_producer = "mcreco";
  bool _save_pi0_tree = true;

  std::vector<std::string> _overburden_volumes = {"volShieldingLid", "volShieldingTop", "volMezzanineLid"};
  std::vector<unsigned int> _pi0_ids;

  double _x_max; //!< x-max of volume box used to determine whether to save track information
  double _x_min; //!< x-min of volume box used to determine whether to save track information
  double _y_max; //!< y-max of volume box used to determine whether to save track information
  double _y_min; //!< y-min of volume box used to determine whether to save track information
  double _z_max; //!< z-max of volume box used to determine whether to save track information
  double _z_min; //!< z-min of volume box used to determine whether to save track information

  boost::uuids::uuid _uuid; ///< A unique ID to identify different events in files with same event number
  std::string _uuid_str; ///< Same as uuid, but converted to string

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

  std::vector<float> _mcp_px;
  std::vector<float> _mcp_py;
  std::vector<float> _mcp_pz;
  std::vector<float> _mcp_e;
  std::vector<float> _mcp_vx;
  std::vector<float> _mcp_vy;
  std::vector<float> _mcp_vz;
  std::vector<float> _mcp_endx;
  std::vector<float> _mcp_endy;
  std::vector<float> _mcp_endz;
  std::vector<float> _mcp_pdg;
  std::vector<float> _mcp_mother;
  std::vector<float> _mcp_status_code;
  std::vector<std::string> _mcp_process;
  std::vector<std::string> _mcp_end_process;
  std::vector<bool> _mcp_intpc;

  std::vector<float> _mcs_pdg;
  std::vector<std::string> _mcs_process;
  std::vector<float> _mcs_start_x;
  std::vector<float> _mcs_start_y;
  std::vector<float> _mcs_start_z;
  std::vector<float> _mcs_end_x;
  std::vector<float> _mcs_end_y;
  std::vector<float> _mcs_end_z;
  std::vector<float> _mcs_start_px;
  std::vector<float> _mcs_start_py;
  std::vector<float> _mcs_start_pz;
  std::vector<float> _mcs_start_e;
  std::vector<float> _mcs_mother_pdg;
  std::vector<std::string> _mcs_mother_process;
  std::vector<float> _mcs_mother_start_x;
  std::vector<float> _mcs_mother_start_y;
  std::vector<float> _mcs_mother_start_z;
  std::vector<float> _mcs_mother_start_e;
  std::vector<float> _mcs_mother_end_x;
  std::vector<float> _mcs_mother_end_y;
  std::vector<float> _mcs_mother_end_z;
  std::vector<float> _mcs_mother_end_e;
  std::vector<float> _mcs_ancestor_pdg;
  std::vector<std::string> _mcs_ancestor_process;
  std::vector<float> _mcs_ancestor_start_e;
  std::vector<float> _mcs_ancestor_end_e;
  std::vector<float> _mcs_mother_in_ob_trackid;
  int _n_mcs_lt1; ///< Number of MC showers with energy less than 10 MeV
  int _n_mcs_lt1_from_ob; ///< Number of MC showers with energy less than 10 MeV and coming from OB

  std::vector<float> _mct_pdg;
  std::vector<std::string> _mct_process;
  std::vector<float> _mct_start_x;
  std::vector<float> _mct_start_y;
  std::vector<float> _mct_start_z;
  std::vector<float> _mct_end_x;
  std::vector<float> _mct_end_y;
  std::vector<float> _mct_end_z;
  std::vector<float> _mct_start_px;
  std::vector<float> _mct_start_py;
  std::vector<float> _mct_start_pz;
  std::vector<float> _mct_start_e;
  std::vector<float> _mct_mother_pdg;
  std::vector<std::string> _mct_mother_process;
  std::vector<float> _mct_ancestor_pdg;
  std::vector<std::string> _mct_ancestor_process;
  std::vector<float> _mct_mother_in_ob_trackid;
  int _n_mct_lt1; ///< Number of MC tracks with energy less than 10 MeV
  int _n_mct_lt1_from_ob; ///< Number of MC tracks with energy less than 10 MeV and coming from OB


  TTree* _pi0_tree; ///< A pi0 TTree, one entry per pi0

  float _pi0_par_e; ///< The pi0 energy
  float _pi0_par_start_x; ///< The pi0 start x
  float _pi0_par_start_y; ///< The pi0 start y
  float _pi0_par_start_z; ///< The pi0 start z
  float _pi0_par_end_x; ///< The pi0 end x
  float _pi0_par_end_y; ///< The pi0 end y
  float _pi0_par_end_z; ///< The pi0 end z
  int _pi0_par_mother_pdg; ///< The pi0 mother pdg
  float _pi0_par_mother_e; ///< The pi0 mother energy
  int _pi0_par_ancestor_trackid; ///< The pi0 primary particle ancestor track id (allows easy pi0 clustering by cosmic interaction)
  std::string _pi0_par_ancestor_uuid; ///< The pi0 unique ID for the event

  std::vector<int> _pi0_daughters_pdg; ///< All the pi0 daughters (usually two photons) (pdg)
  std::vector<float> _pi0_daughters_e; ///< All the pi0 daughters (usually two photons) (energy)
  // std::vector<std::string> _pi0_daughters_startprocess; ///< All the pi0 daughters (usually two photons) (energy) (start process)
  // std::vector<std::string> _pi0_daughters_endprocess; ////< All the pi0 daughters (usually two photons) (energy) (end process)
  std::vector<float> _pi0_daughters_start_x; ///< All the pi0 daughters (usually two photons) (start x)
  std::vector<float> _pi0_daughters_start_y; ///< All the pi0 daughters (usually two photons) (start y)
  std::vector<float> _pi0_daughters_start_z; ///< All the pi0 daughters (usually two photons) (start z)
  std::vector<float> _pi0_daughters_end_x; ///< All the pi0 daughters (usually two photons) (end x)
  std::vector<float> _pi0_daughters_end_y; ///< All the pi0 daughters (usually two photons) (end y)
  std::vector<float> _pi0_daughters_end_z; ///< All the pi0 daughters (usually two photons) (end z)

  std::vector<int> _pi0_event_particles_pdg; ///< All the particles produced together with the pi0
  std::vector<float> _pi0_event_particles_e; ///< All the particles produced together with the pi0

  std::vector<int> _pi0_genealogy_pdg; ///< The full pi0 genealogy, from its mother all the way to the ancestor (pdg)
  std::vector<std::string> _pi0_genealogy_startprocess; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start process)
  std::vector<std::string> _pi0_genealogy_endprocess; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end process)
  std::vector<int> _pi0_genealogy_mother; ///< The full pi0 genealogy, from its mother all the way to the ancestor (mother track id)
  std::vector<float> _pi0_genealogy_e; ///< The full pi0 genealogy, from its mother all the way to the ancestor (energy)
  std::vector<float> _pi0_genealogy_start_x; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start x)
  std::vector<float> _pi0_genealogy_start_y; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start y)
  std::vector<float> _pi0_genealogy_start_z; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start z)
  std::vector<float> _pi0_genealogy_end_x; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end x)
  std::vector<float> _pi0_genealogy_end_y; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end y)
  std::vector<float> _pi0_genealogy_end_z; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end z)
  std::vector<float> _pi0_genealogy_trackid; ///< The full pi0 genealogy, from its mother all the way to the ancestor (trackid)



  TTree* _sr_tree; ///< TTree filled per subrun
  int _sr_run; ///< Subrun Run number
  int _sr_subrun; ///< Subrun Subrun number
  double _sr_begintime; ///< Subrun start time
  double _sr_endtime; ///< Subrun end time
  double _sr_pot; ///< Subrun POT
};


OverburdenAna::OverburdenAna(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("OBTree","");

  _tree->Branch("run",     &_run,     "run/I");
  _tree->Branch("subrun",     &_subrun,     "subrun/I");
  _tree->Branch("event",     &_event,     "event/I");

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

  _tree->Branch("mcp_px", "std::vector<float>", &_mcp_px);
  _tree->Branch("mcp_py", "std::vector<float>", &_mcp_py);
  _tree->Branch("mcp_pz", "std::vector<float>", &_mcp_pz);
  _tree->Branch("mcp_e", "std::vector<float>", &_mcp_e);
  _tree->Branch("mcp_vx", "std::vector<float>", &_mcp_vx);
  _tree->Branch("mcp_vy", "std::vector<float>", &_mcp_vy);
  _tree->Branch("mcp_vz", "std::vector<float>", &_mcp_vz);
  _tree->Branch("mcp_endx", "std::vector<float>", &_mcp_endx);
  _tree->Branch("mcp_endy", "std::vector<float>", &_mcp_endy);
  _tree->Branch("mcp_endz", "std::vector<float>", &_mcp_endz);
  _tree->Branch("mcp_pdg", "std::vector<float>", &_mcp_pdg);
  _tree->Branch("mcp_mother", "std::vector<float>", &_mcp_mother);
  _tree->Branch("mcp_status_code", "std::vector<float>", &_mcp_status_code);
  _tree->Branch("mcp_process", "std::vector<std::string>", &_mcp_process);
  _tree->Branch("mcp_end_process", "std::vector<std::string>", &_mcp_end_process);
  _tree->Branch("mcp_intpc", "std::vector<bool>", &_mcp_intpc);

  _tree->Branch("mct_pdg", "std::vector<float>", &_mct_pdg);
  _tree->Branch("mct_process", "std::vector<std::string>", &_mct_process);
  _tree->Branch("mct_start_x", "std::vector<float>", &_mct_start_x);
  _tree->Branch("mct_start_y", "std::vector<float>", &_mct_start_y);
  _tree->Branch("mct_start_z", "std::vector<float>", &_mct_start_z);
  _tree->Branch("mct_end_x", "std::vector<float>", &_mct_end_x);
  _tree->Branch("mct_end_y", "std::vector<float>", &_mct_end_y);
  _tree->Branch("mct_end_z", "std::vector<float>", &_mct_end_z);
  _tree->Branch("mct_start_px", "std::vector<float>", &_mct_start_px);
  _tree->Branch("mct_start_py", "std::vector<float>", &_mct_start_py);
  _tree->Branch("mct_start_pz", "std::vector<float>", &_mct_start_pz);
  _tree->Branch("mct_start_e", "std::vector<float>", &_mct_start_e);
  _tree->Branch("mct_mother_pdg", "std::vector<float>", &_mct_mother_pdg);
  _tree->Branch("mct_ancestor_pdg", "std::vector<float>", &_mct_ancestor_pdg);
  _tree->Branch("mct_mother_process", "std::vector<std::string>>", &_mct_mother_process);
  _tree->Branch("mct_ancestor_process", "std::vector<std::string>>", &_mct_ancestor_process);
  _tree->Branch("mct_mother_in_ob_trackid", "std::vector<float>", &_mct_mother_in_ob_trackid);
  _tree->Branch("mct_n_lt1",     &_n_mct_lt1,     "mct_n_lt1/I");
  _tree->Branch("mct_n_lt1_from_ob",     &_n_mct_lt1_from_ob,     "mct_n_lt1_from_ob/I");

  _tree->Branch("mcs_pdg", "std::vector<float>", &_mcs_pdg);
  _tree->Branch("mcs_process", "std::vector<std::string>", &_mcs_process);
  _tree->Branch("mcs_start_x", "std::vector<float>", &_mcs_start_x);
  _tree->Branch("mcs_start_y", "std::vector<float>", &_mcs_start_y);
  _tree->Branch("mcs_start_z", "std::vector<float>", &_mcs_start_z);
  _tree->Branch("mcs_end_x", "std::vector<float>", &_mcs_end_x);
  _tree->Branch("mcs_end_y", "std::vector<float>", &_mcs_end_y);
  _tree->Branch("mcs_end_z", "std::vector<float>", &_mcs_end_z);
  _tree->Branch("mcs_start_px", "std::vector<float>", &_mcs_start_px);
  _tree->Branch("mcs_start_py", "std::vector<float>", &_mcs_start_py);
  _tree->Branch("mcs_start_pz", "std::vector<float>", &_mcs_start_pz);
  _tree->Branch("mcs_start_e", "std::vector<float>", &_mcs_start_e);
  _tree->Branch("mcs_mother_pdg", "std::vector<float>", &_mcs_mother_pdg);
  _tree->Branch("mcs_ancestor_pdg", "std::vector<float>", &_mcs_ancestor_pdg);
  _tree->Branch("mcs_mother_start_x", "std::vector<float>", &_mcs_mother_start_x);
  _tree->Branch("mcs_mother_start_y", "std::vector<float>", &_mcs_mother_start_y);
  _tree->Branch("mcs_mother_start_z", "std::vector<float>", &_mcs_mother_start_z);
  _tree->Branch("mcs_mother_start_e", "std::vector<float>", &_mcs_mother_start_e);
  _tree->Branch("mcs_mother_end_x", "std::vector<float>", &_mcs_mother_end_x);
  _tree->Branch("mcs_mother_end_y", "std::vector<float>", &_mcs_mother_end_y);
  _tree->Branch("mcs_mother_end_z", "std::vector<float>", &_mcs_mother_end_z);
  _tree->Branch("mcs_mother_end_e", "std::vector<float>", &_mcs_mother_end_e);
  _tree->Branch("mcs_mother_process", "std::vector<std::string>>", &_mcs_mother_process);
  _tree->Branch("mcs_ancestor_process", "std::vector<std::string>>", &_mcs_ancestor_process);
  _tree->Branch("mcs_ancestor_start_e", "std::vector<float>", &_mcs_ancestor_start_e);
  _tree->Branch("mcs_ancestor_end_e", "std::vector<float>", &_mcs_ancestor_end_e);
  _tree->Branch("mcs_mother_in_ob_trackid", "std::vector<float>", &_mcs_mother_in_ob_trackid);
  _tree->Branch("mcs_n_lt1",     &_n_mcs_lt1,     "mcs_n_lt1/I");
  _tree->Branch("mcs_n_lt1_from_ob",     &_n_mcs_lt1_from_ob,     "mcs_n_lt1_from_ob/I");

  if (_save_pi0_tree) {
    _pi0_tree = fs->make<TTree>("Pi0Tree","");

    _pi0_tree->Branch("run", &_run, "run/I");
    _pi0_tree->Branch("subrun", &_subrun, "subrun/I");
    _pi0_tree->Branch("event", &_event, "event/I");

    _pi0_tree->Branch("pi0_par_e", &_pi0_par_e, "pi0_par_e/F");
    _pi0_tree->Branch("pi0_par_start_x", &_pi0_par_start_x, "pi0_par_start_x/F");
    _pi0_tree->Branch("pi0_par_start_y", &_pi0_par_start_y, "pi0_par_start_y/F");
    _pi0_tree->Branch("pi0_par_start_z", &_pi0_par_start_z, "pi0_par_start_z/F");
    _pi0_tree->Branch("pi0_par_end_x", &_pi0_par_end_x, "pi0_par_end_x/F");
    _pi0_tree->Branch("pi0_par_end_y", &_pi0_par_end_y, "pi0_par_end_y/F");
    _pi0_tree->Branch("pi0_par_end_z", &_pi0_par_end_z, "pi0_par_end_z/F");
    _pi0_tree->Branch("pi0_par_mother_pdg", &_pi0_par_mother_pdg, "pi0_par_mother_pdg/I");
    _pi0_tree->Branch("pi0_par_mother_e", &_pi0_par_mother_e, "pi0_par_mother_e/F");
    _pi0_tree->Branch("pi0_par_ancestor_trackid", &_pi0_par_ancestor_trackid, "pi0_par_ancestor_trackid/I");
    _pi0_tree->Branch("pi0_par_ancestor_uuid", "std::string", &_pi0_par_ancestor_uuid);

    _pi0_tree->Branch("pi0_event_particles_pdg", "std::vector<int>", &_pi0_event_particles_pdg);
    _pi0_tree->Branch("pi0_event_particles_e", "std::vector<float>", &_pi0_event_particles_e);

    _pi0_tree->Branch("pi0_daughters_pdg", "std::vector<int>", &_pi0_daughters_pdg);
    _pi0_tree->Branch("pi0_daughters_e", "std::vector<float>", &_pi0_daughters_e);
    // _pi0_tree->Branch("pi0_daughters_startprocess", "std::vector<std::string>", &_pi0_daughters_startprocess);
    // _pi0_tree->Branch("pi0_daughters_endprocess", "std::vector<std::string>", &_pi0_daughters_endprocess);
    _pi0_tree->Branch("pi0_daughters_start_x", "std::vector<float>", &_pi0_daughters_start_x);
    _pi0_tree->Branch("pi0_daughters_start_y", "std::vector<float>", &_pi0_daughters_start_y);
    _pi0_tree->Branch("pi0_daughters_start_z", "std::vector<float>", &_pi0_daughters_start_z);
    _pi0_tree->Branch("pi0_daughters_end_x", "std::vector<float>", &_pi0_daughters_end_x);
    _pi0_tree->Branch("pi0_daughters_end_y", "std::vector<float>", &_pi0_daughters_end_y);
    _pi0_tree->Branch("pi0_daughters_end_z", "std::vector<float>", &_pi0_daughters_end_z);

    _pi0_tree->Branch("pi0_genealogy_pdg", "std::vector<int>", &_pi0_genealogy_pdg);
    _pi0_tree->Branch("pi0_genealogy_startprocess", "std::vector<std::string>", &_pi0_genealogy_startprocess);
    _pi0_tree->Branch("pi0_genealogy_endprocess", "std::vector<std::string>", &_pi0_genealogy_endprocess);
    _pi0_tree->Branch("pi0_genealogy_mother", "std::vector<int>", &_pi0_genealogy_mother);
    _pi0_tree->Branch("pi0_genealogy_e", "std::vector<float>", &_pi0_genealogy_e);
    _pi0_tree->Branch("pi0_genealogy_start_x", "std::vector<float>", &_pi0_genealogy_start_x);
    _pi0_tree->Branch("pi0_genealogy_start_y", "std::vector<float>", &_pi0_genealogy_start_y);
    _pi0_tree->Branch("pi0_genealogy_start_z", "std::vector<float>", &_pi0_genealogy_start_z);
    _pi0_tree->Branch("pi0_genealogy_end_x", "std::vector<float>", &_pi0_genealogy_end_x);
    _pi0_tree->Branch("pi0_genealogy_end_y", "std::vector<float>", &_pi0_genealogy_end_y);
    _pi0_tree->Branch("pi0_genealogy_end_z", "std::vector<float>", &_pi0_genealogy_end_z);
    _pi0_tree->Branch("pi0_genealogy_trackid", "std::vector<float>", &_pi0_genealogy_trackid);

  }

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");

  // Iterate over all TPC's to get bounding box that covers volumes of each individual TPC in the detector
  art::ServiceHandle<geo::Geometry const> geo;
  _x_min = std::min_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MinX() < rhs.BoundingBox().MinX();})->MinX();
  _y_min = std::min_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MinY() < rhs.BoundingBox().MinY();})->MinY();
  _z_min = std::min_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MinZ() < rhs.BoundingBox().MinZ();})->MinZ();
  _x_max = std::max_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MaxX() < rhs.BoundingBox().MaxX();})->MaxX();
  _y_max = std::max_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MaxY() < rhs.BoundingBox().MaxY();})->MaxY();
  _z_max = std::max_element(geo->begin_TPC(), geo->end_TPC(), [](auto const &lhs, auto const &rhs){ return lhs.BoundingBox().MaxZ() < rhs.BoundingBox().MaxZ();})->MaxZ();

  std::cout << "TPC limits: " << std::endl;
  std::cout << "\tx_max" << _x_max << std::endl;
  std::cout << "\tx_min" << _x_min << std::endl;
  std::cout << "\ty_max" << _y_max << std::endl;
  std::cout << "\ty_min" << _y_min << std::endl;
  std::cout << "\tz_max" << _z_max << std::endl;
  std::cout << "\tz_min" << _z_min << std::endl;
}

void OverburdenAna::produce(art::Event& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  _uuid = boost::uuids::random_generator()();
  _uuid_str = boost::lexical_cast<std::string>(_uuid) + "-" + _event;
  std::cout << "Event uuid " << _uuid_str << std::endl;

  std::set<std::string> volnameset(_overburden_volumes.begin(), _overburden_volumes.end());
  _part_filter = CreateParticleVolumeFilter(volnameset);

  bool in_ob = _part_filter->mustKeep(Point_t{{0, 800, 0}});
  std::cout << "Is {0, 800, 0} in the OB (should be yes)? " << (in_ob ? "YES" : "NO") << std::endl;
  in_ob = _part_filter->mustKeep(Point_t{{0, 300, 0}});
  std::cout << "Is {0, 300, 0} in the OB (should be no)? " << (in_ob ? "YES" : "NO") << std::endl;

  // std::vector<std::string> names = {"volCryostat"};
  // std::set<std::string> set(names.begin(), names.end());
  // auto filter = CreateParticleVolumeFilter(set);
  // std::cout << "NEW Is {0, 0, 100} in the volCryostat? " << (filter->mustKeep(Point_t{{0, 0, 100}}) ? "YES" : "NO") << std::endl;
  // std::cout << "NEW Is {0, 800, 100} in the volCryostat? " << (filter->mustKeep(Point_t{{0, 800, 100}}) ? "YES" : "NO") << std::endl;


  //
  // MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_producer, mct_h);
  if(mct_h.isValid()){

    std::vector<art::Ptr<simb::MCTruth>> mct_v;
    art::fill_ptr_vector(mct_v, mct_h);

    //
    // Loop over the neutrino interactions in this event
    //
    for (size_t i = 0; i < mct_v.size(); i++) {
      // if (mct_v.at(i)->Origin() != simb::kBeamNeutrino) {
      //   std::cout << "[OverburdenAna] MCTruth from generator does not have neutrino origin?!" << std::endl;
      // }

      if(!mct_v[i]->NeutrinoSet()) {
        break;
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
    }
  } else {
    std::cout << "MCTruth product " << _mctruth_producer << " not found..." << std::endl;
  }

  //
  // MCParticle
  //
  art::Handle<std::vector<simb::MCParticle> > mcp_h;
  e.getByLabel(_mcparticle_producer, mcp_h);
  if(!mcp_h.isValid()){
    std::cout << "MCParticle product " << _mcparticle_producer << " not found..." << std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcp_h);

  clear_vectors();

  for (size_t i = 0; i < mcp_v.size(); i++) {
    auto mcp = mcp_v.at(i);

    _trackid_to_mcparticle[mcp->TrackId()] = *mcp;

    // bool in_det = InDetector(mcp);

    // Only save the MCP if it's a primary
    // if (mcp->Process() == "primary" || in_det) {
    if (mcp->Process() == "primary") {

      _mcp_px.push_back(mcp->Px());
      _mcp_py.push_back(mcp->Py());
      _mcp_pz.push_back(mcp->Pz());
      _mcp_e.push_back(mcp->E());

      _mcp_vx.push_back(mcp->Vx());
      _mcp_vy.push_back(mcp->Vy());
      _mcp_vz.push_back(mcp->Vz());
      _mcp_endx.push_back(mcp->EndX());
      _mcp_endy.push_back(mcp->EndY());
      _mcp_endz.push_back(mcp->EndZ());

      _mcp_pdg.push_back(mcp->PdgCode());
      _mcp_mother.push_back(mcp->Mother());
      _mcp_status_code.push_back(mcp->StatusCode());
      _mcp_process.push_back(mcp->Process());
      _mcp_end_process.push_back(mcp->EndProcess());

      _mcp_intpc.push_back(InDetector(mcp));
    }
  }



  //
  // MCTrack
  //
  art::Handle<std::vector<sim::MCTrack> > mc_track_h;
  e.getByLabel(_mctrack_producer, mc_track_h);
  if(!mc_track_h.isValid()){
    std::cout << "MCTrack product " << _mctrack_producer << " not found..." << std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<sim::MCTrack>> mc_track_v;
  art::fill_ptr_vector(mc_track_v, mc_track_h);

  for (size_t i = 0; i < mc_track_v.size(); i++) {
    auto mc_track = mc_track_v.at(i);

    // std::cout << "MCTrack " << i << ": ancestor pdg " << mc_track->AncestorPdgCode()
    //                               << ", ancestor process " << mc_track->AncestorProcess()
    //                               << ", mother pdg " << mc_track->MotherPdgCode()
    //                               << ", mother process " << mc_track->MotherProcess()
    //                              << "| PDG " << mc_track->PdgCode()
    //                               << ", process " << mc_track->Process()
    //                               << std::endl;

    auto iter = _trackid_to_mcparticle.find(mc_track->TrackID());
    int mother_in_ob = -1;
    if (iter != _trackid_to_mcparticle.end()) {
      mother_in_ob = FindMotherInOverburden(iter->second);
    }

    // Don't save MCT with energy less than 1 MeV
    if (mc_track->Start().E() < 1) { // MeV
      _n_mct_lt1 ++;
      if (mother_in_ob != -1){
        _n_mct_lt1_from_ob ++;
      }
      continue;
    }

    // Don't save MCS that are not in the TPCs
    if (mc_track->size() == 0) {
      continue;
    }


    _mct_pdg.push_back(mc_track->PdgCode());
    _mct_process.push_back(mc_track->Process());

    _mct_start_x.push_back(mc_track->Start().X());
    _mct_start_y.push_back(mc_track->Start().Y());
    _mct_start_z.push_back(mc_track->Start().Z());

    _mct_end_x.push_back(mc_track->End().X());
    _mct_end_y.push_back(mc_track->End().Y());
    _mct_end_z.push_back(mc_track->End().Z());

    _mct_start_px.push_back(mc_track->Start().Px());
    _mct_start_py.push_back(mc_track->Start().Py());
    _mct_start_pz.push_back(mc_track->Start().Pz());
    _mct_start_e.push_back(mc_track->Start().E());

    _mct_mother_pdg.push_back(mc_track->MotherPdgCode());
    _mct_mother_process.push_back(mc_track->MotherProcess());
    _mct_ancestor_pdg.push_back(mc_track->AncestorPdgCode());
    _mct_ancestor_process.push_back(mc_track->AncestorProcess());

    _mct_mother_in_ob_trackid.push_back(mother_in_ob);

  }


  //
  // MCShower
  //
  art::Handle<std::vector<sim::MCShower> > mc_shower_h;
  e.getByLabel(_mcshower_producer, mc_shower_h);
  if(!mc_shower_h.isValid()){
    std::cout << "MCShower product " << _mcshower_producer << " not found..." << std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<sim::MCShower>> mc_shower_v;
  art::fill_ptr_vector(mc_shower_v, mc_shower_h);

  for (size_t i = 0; i < mc_shower_v.size(); i++) {
    auto mc_shower = mc_shower_v.at(i);

    // std::cout << "MCShower " << i << ": ancestor pdg " << mc_shower->AncestorPdgCode()
    //                               << ", ancestor process " << mc_shower->AncestorProcess()
    //                               << ", mother pdg " << mc_shower->MotherPdgCode()
    //                               << ", mother process " << mc_shower->MotherProcess()
    //                               << "| PDG " << mc_shower->PdgCode()
    //                               << ", process " << mc_shower->Process()
    //                               << std::endl;

    auto iter = _trackid_to_mcparticle.find(mc_shower->TrackID());
    int mother_in_ob = -1;
    if (iter != _trackid_to_mcparticle.end()) {
      mother_in_ob = FindMotherInOverburden(iter->second);
    }

    // Don't save MCS with energy less than 1 MeV
    if (mc_shower->Start().E() < 1) { // MeV
      _n_mcs_lt1 ++;
      if (mother_in_ob != -1){
        _n_mcs_lt1_from_ob ++;
      }
      continue;
    }

    // Don't save MCS that are not in the TPCs
    // Special case for photon showers, which can start outside
    // bool start_in_det = InDetector(mc_shower->Start().X(), mc_shower->Start().Y(), mc_shower->Start().Z());
    bool end_in_det = InDetector(mc_shower->End().X(), mc_shower->End().Y(), mc_shower->End().Z());
    // if (!(start_in_det || end_in_det)) {
    //   continue;
    // }
    // if (mc_shower->PdgCode() == 22 && !end_in_det) {
    //   continue;
    // } else if (mc_shower->PdgCode() != 22 && !(start_in_det || end_in_det)) {
    //   continue;
    // }
    if (!end_in_det) {
      continue;
    }

    _mcs_pdg.push_back(mc_shower->PdgCode());
    _mcs_process.push_back(mc_shower->Process());

    _mcs_start_x.push_back(mc_shower->Start().X());
    _mcs_start_y.push_back(mc_shower->Start().Y());
    _mcs_start_z.push_back(mc_shower->Start().Z());

    _mcs_end_x.push_back(mc_shower->End().X());
    _mcs_end_y.push_back(mc_shower->End().Y());
    _mcs_end_z.push_back(mc_shower->End().Z());

    _mcs_start_px.push_back(mc_shower->Start().Px());
    _mcs_start_py.push_back(mc_shower->Start().Py());
    _mcs_start_pz.push_back(mc_shower->Start().Pz());
    _mcs_start_e.push_back(mc_shower->Start().E());

    _mcs_mother_pdg.push_back(mc_shower->MotherPdgCode());
    _mcs_mother_process.push_back(mc_shower->MotherProcess());
    _mcs_mother_start_x.push_back(mc_shower->MotherStart().X());
    _mcs_mother_start_y.push_back(mc_shower->MotherStart().Y());
    _mcs_mother_start_z.push_back(mc_shower->MotherStart().Z());
    _mcs_mother_start_e.push_back(mc_shower->MotherStart().E());
    _mcs_mother_end_x.push_back(mc_shower->MotherEnd().X());
    _mcs_mother_end_y.push_back(mc_shower->MotherEnd().Y());
    _mcs_mother_end_z.push_back(mc_shower->MotherEnd().Z());
    _mcs_mother_end_e.push_back(mc_shower->MotherEnd().E());

    _mcs_ancestor_pdg.push_back(mc_shower->AncestorPdgCode());
    _mcs_ancestor_process.push_back(mc_shower->AncestorProcess());
    _mcs_ancestor_start_e.push_back(mc_shower->AncestorStart().E());
    _mcs_ancestor_end_e.push_back(mc_shower->AncestorEnd().E());

    _mcs_mother_in_ob_trackid.push_back(mother_in_ob);

    if (mc_shower->MotherPdgCode() == 111 && _save_pi0_tree) {
      SavePi0ShowerInfo(mc_shower->MotherTrackID());
    }


  }


   _tree->Fill();
}

void OverburdenAna::SavePi0ShowerInfo(int pi0_track_id) {
  std::cout << "SavePi0ShowerInfo***, pi0_track_id: " << pi0_track_id << std::endl;
  auto it = std::find(_pi0_ids.begin(), _pi0_ids.end(), pi0_track_id);
  if (it != _pi0_ids.end()) {
    return;
  }
  std::cout << "SavePi0ShowerInfo***, got it" << std::endl;
  _pi0_ids.push_back(pi0_track_id);

  // Get the pi0 MCParticle
  auto iter = _trackid_to_mcparticle.find(pi0_track_id);
  if (iter == _trackid_to_mcparticle.end()) {
    return;
  }
  simb::MCParticle pi0_mcp = iter->second;
  std::cout << "Pi0 MCP, PDG = " << pi0_mcp.PdgCode() << ", E = " << pi0_mcp.E() << ", n daughters " << pi0_mcp.NumberDaughters() << std::endl;

  // Save the information on the pi0 itself
  _pi0_par_e = pi0_mcp.E();
  _pi0_par_start_x = pi0_mcp.Vx();
  _pi0_par_start_y = pi0_mcp.Vy();
  _pi0_par_start_z = pi0_mcp.Vz();
  _pi0_par_end_x = pi0_mcp.EndX();
  _pi0_par_end_y = pi0_mcp.EndY();
  _pi0_par_end_z = pi0_mcp.EndZ();

  // Get the pi0 mother MCParticle
  iter = _trackid_to_mcparticle.find(pi0_mcp.Mother());
  if (iter == _trackid_to_mcparticle.end()) {
    return;
  }
  simb::MCParticle pi0_mother_mcp = iter->second;
  _pi0_par_mother_pdg = pi0_mother_mcp.PdgCode();
  _pi0_par_mother_e = pi0_mother_mcp.E();
  std::cout << "Pi0 MCP Mother, PDG = " << pi0_mother_mcp.PdgCode() << ", E = " << pi0_mother_mcp.E() << std::endl;

  // Get the daughters of the pi0 mother
  for (int d = 0; d < pi0_mother_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(pi0_mother_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _pi0_event_particles_pdg.push_back(daughter.PdgCode());
      _pi0_event_particles_e.push_back(daughter.E());
    }
  }

  // Save the mother first...
  _pi0_genealogy_pdg.push_back(pi0_mother_mcp.PdgCode());
  _pi0_genealogy_startprocess.push_back(pi0_mother_mcp.Process());
  _pi0_genealogy_endprocess.push_back(pi0_mother_mcp.EndProcess());
  _pi0_genealogy_mother.push_back(pi0_mother_mcp.Mother());
  _pi0_genealogy_e.push_back(pi0_mother_mcp.E());
  _pi0_genealogy_start_x.push_back(pi0_mother_mcp.Vx());
  _pi0_genealogy_start_y.push_back(pi0_mother_mcp.Vy());
  _pi0_genealogy_start_z.push_back(pi0_mother_mcp.Vy());
  _pi0_genealogy_end_x.push_back(pi0_mother_mcp.EndX());
  _pi0_genealogy_end_y.push_back(pi0_mother_mcp.EndY());
  _pi0_genealogy_end_z.push_back(pi0_mother_mcp.EndZ());
  _pi0_genealogy_trackid.push_back(pi0_mother_mcp.TrackId());


  // ... then save all the other ancestors
  simb::MCParticle mcp = pi0_mother_mcp;
  while(true) {
    iter = _trackid_to_mcparticle.find(mcp.Mother());
    if (iter == _trackid_to_mcparticle.end()) {
      break;
    }
    mcp = iter->second;
    _pi0_genealogy_pdg.push_back(mcp.PdgCode());
    _pi0_genealogy_startprocess.push_back(mcp.Process());
    _pi0_genealogy_endprocess.push_back(mcp.EndProcess());
    _pi0_genealogy_mother.push_back(mcp.Mother());
    _pi0_genealogy_e.push_back(mcp.E());
    _pi0_genealogy_start_x.push_back(mcp.Vx());
    _pi0_genealogy_start_y.push_back(mcp.Vy());
    _pi0_genealogy_start_z.push_back(mcp.Vy());
    _pi0_genealogy_end_x.push_back(mcp.EndX());
    _pi0_genealogy_end_y.push_back(mcp.EndY());
    _pi0_genealogy_end_z.push_back(mcp.EndZ());
    _pi0_genealogy_trackid.push_back(mcp.TrackId());
  }

  _pi0_par_ancestor_trackid = _pi0_genealogy_trackid.back();
  _pi0_par_ancestor_uuid = _uuid_str + "-" + std::to_string(_pi0_par_ancestor_trackid);
  std::cout << "Pi0 uuid " << _pi0_par_ancestor_uuid << std::endl;

  // Get the dautghers of this pi0
  for (int d = 0; d < pi0_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(pi0_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _pi0_daughters_pdg.push_back(daughter.PdgCode());
      _pi0_daughters_e.push_back(daughter.E());
      // _pi0_daughters_startprocess.push_back(daughter.Process());
      // _pi0_daughters_endprocess.push_back(daughter.EndProcess());
      _pi0_daughters_start_x.push_back(daughter.Vx());
      _pi0_daughters_start_y.push_back(daughter.Vy());
      _pi0_daughters_start_z.push_back(daughter.Vz());
      _pi0_daughters_end_x.push_back(daughter.EndX());
      _pi0_daughters_end_y.push_back(daughter.EndY());
      _pi0_daughters_end_z.push_back(daughter.EndZ());
    }
  }

  // Fill the tree and reset the variables
  _pi0_tree->Fill();

  _pi0_event_particles_pdg.clear();
  _pi0_event_particles_e.clear();

  _pi0_daughters_pdg.clear();
  _pi0_daughters_e.clear();
  // _pi0_daughters_startprocess.clear();
  // _pi0_daughters_endprocess.clear();
  _pi0_daughters_start_x.clear();
  _pi0_daughters_start_y.clear();
  _pi0_daughters_start_z.clear();
  _pi0_daughters_end_x.clear();
  _pi0_daughters_end_y.clear();
  _pi0_daughters_end_z.clear();

  _pi0_genealogy_pdg.clear();
  _pi0_genealogy_startprocess.clear();
  _pi0_genealogy_endprocess.clear();
  _pi0_genealogy_mother.clear();
  _pi0_genealogy_e.clear();
  _pi0_genealogy_start_x.clear();
  _pi0_genealogy_start_y.clear();
  _pi0_genealogy_start_z.clear();
  _pi0_genealogy_end_x.clear();
  _pi0_genealogy_end_y.clear();
  _pi0_genealogy_end_z.clear();
  _pi0_genealogy_trackid.clear();

}

int OverburdenAna::FindMotherInOverburden(simb::MCParticle mcp) {

  if (mcp.Process() == "primary") {
    return -1;
  }

  // We use this filter not to actualy filter, but to check if
  // the vertex of this particle is in the OB
  bool vtx_in_ob = _part_filter->mustKeep(Point_t{{ mcp.Vx(),
                                                    mcp.Vy(),
                                                    mcp.Vz() }});

  auto iter = _trackid_to_mcparticle.find(mcp.Mother());
  if (iter == _trackid_to_mcparticle.end()) {
    return -1;
  }
  auto mother = iter->second;


  if (vtx_in_ob                          // If this particle has a vertex in the OB
    && mcp.Process() != "primary"        // and this particle is not a primary one
    // && mother.Process() == "primary"     // and the mother of it is a primary
    ) {
    return mcp.TrackId();                // Then return it, as is something created by a primary in the OB
  }

  // Otherwise, keep looking
  return FindMotherInOverburden(mother);

}

std::unique_ptr<util::PositionInVolumeFilter> OverburdenAna::CreateParticleVolumeFilter
    (std::set<std::string> const& vol_names) const
  {

    // if we don't have favourite volumes, don't even bother creating a filter
    if (vol_names.empty()) return {};

    auto const& geom = *art::ServiceHandle<geo::Geometry const>();

    std::vector<std::vector<TGeoNode const*>> node_paths
      = geom.FindAllVolumePaths(vol_names);
    std::cout << "Found " << node_paths.size() << " node paths." << std::endl;

    // collection of interesting volumes
    util::PositionInVolumeFilter::AllVolumeInfo_t GeoVolumePairs;
    GeoVolumePairs.reserve(node_paths.size()); // because we are obsessed

    //for each interesting volume, follow the node path and collect
    //total rotations and translations
    for (size_t iVolume = 0; iVolume < node_paths.size(); ++iVolume) {
      std::vector<TGeoNode const*> path = node_paths[iVolume];

      TGeoTranslation* pTransl = new TGeoTranslation(0.,0.,0.);
      TGeoRotation* pRot = new TGeoRotation();
      for (TGeoNode const* node: path) {
        TGeoTranslation thistranslate(*node->GetMatrix());
        TGeoRotation thisrotate(*node->GetMatrix());
        pTransl->Add(&thistranslate);
        *pRot=*pRot * thisrotate;
      }

      //for some reason, pRot and pTransl don't have tr and rot bits set correctly
      //make new translations and rotations so bits are set correctly
      TGeoTranslation* pTransl2 = new TGeoTranslation(pTransl->GetTranslation()[0],
                                                        pTransl->GetTranslation()[1],
                                                      pTransl->GetTranslation()[2]);
      double phi=0.,theta=0.,psi=0.;
      pRot->GetAngles(phi,theta,psi);
      TGeoRotation* pRot2 = new TGeoRotation();
      pRot2->SetAngles(phi,theta,psi);

      TGeoCombiTrans* pTransf = new TGeoCombiTrans(*pTransl2,*pRot2);

      GeoVolumePairs.emplace_back(node_paths[iVolume].back()->GetVolume(), pTransf);

    }

    return std::make_unique<util::PositionInVolumeFilter>(std::move(GeoVolumePairs));

  } // CreateParticleVolumeFilter()

void OverburdenAna::clear_vectors() {

  _n_mct_lt1 = 0;
  _n_mct_lt1_from_ob = 0;
  _n_mcs_lt1 = 0;
  _n_mcs_lt1_from_ob = 0;

  _mcp_px.clear();
  _mcp_py.clear();
  _mcp_pz.clear();
  _mcp_e.clear();
  _mcp_vx.clear();
  _mcp_vy.clear();
  _mcp_vz.clear();
  _mcp_endx.clear();
  _mcp_endy.clear();
  _mcp_endz.clear();
  _mcp_pdg.clear();
  _mcp_mother.clear();
  _mcp_status_code.clear();
  _mcp_process.clear();
  _mcp_end_process.clear();
  _mcp_intpc.clear();

  _mcs_pdg.clear();
  _mcs_process.clear();
  _mcs_start_x.clear();
  _mcs_start_y.clear();
  _mcs_start_z.clear();
  _mcs_end_x.clear();
  _mcs_end_y.clear();
  _mcs_end_z.clear();
  _mcs_start_px.clear();
  _mcs_start_py.clear();
  _mcs_start_pz.clear();
  _mcs_start_e.clear();
  _mcs_mother_pdg.clear();
  _mcs_mother_process.clear();
  _mcs_mother_start_x.clear();
  _mcs_mother_start_y.clear();
  _mcs_mother_start_z.clear();
  _mcs_mother_start_e.clear();
  _mcs_mother_end_x.clear();
  _mcs_mother_end_y.clear();
  _mcs_mother_end_z.clear();
  _mcs_mother_end_e.clear();
  _mcs_ancestor_pdg.clear();
  _mcs_ancestor_process.clear();
  _mcs_ancestor_start_e.clear();
  _mcs_ancestor_end_e.clear();
  _mcs_mother_in_ob_trackid.clear();

  _mct_pdg.clear();
  _mct_process.clear();
  _mct_start_x.clear();
  _mct_start_y.clear();
  _mct_start_z.clear();
  _mct_end_x.clear();
  _mct_end_y.clear();
  _mct_end_z.clear();
  _mct_start_px.clear();
  _mct_start_py.clear();
  _mct_start_pz.clear();
  _mct_start_e.clear();
  _mct_mother_pdg.clear();
  _mct_mother_process.clear();
  _mct_ancestor_pdg.clear();
  _mct_ancestor_process.clear();
  _mct_mother_in_ob_trackid.clear();

  _pi0_ids.clear();
}

bool OverburdenAna::InDetector(const double& x,
                               const double& y,
                               const double& z) const {
  return !( x > _x_max || x < _x_min ||
            z > _z_max || z < _z_min ||
            y > _y_max || y < _y_min );
}

bool OverburdenAna::InDetector(art::Ptr<simb::MCParticle> mcp) {
  auto t = mcp->Trajectory();
  for (size_t i = 0; i < t.size(); i++) {
    if (InDetector(t.X(i), t.Y(i), t.Z(i))) return true;
  }
  return false;
}


void OverburdenAna::beginSubRun(art::SubRun & sr) {

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



DEFINE_ART_MODULE(OverburdenAna)
