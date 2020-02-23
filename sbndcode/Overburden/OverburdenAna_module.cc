////////////////////////////////////////////////////////////////////////
// Class:       OverburdenAna
// Plugin Type: producer (art v3_04_00)
// File:        OverburdenAna_module.cc
//
// Generated at Tue Jan 28 20:46:15 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

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

#include <memory>
#include <map>

#include "TTree.h"
// #include "TGeoVolume.h"
// #include "TGeoMatrix.h" // TGeoCombiTrans
// #include "TLorentzVector.h"
// #include "TVector3.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

#include "larcorealg/CoreUtils/ParticleFilters.h" // util::PositionInVolumeFilter
#include "larcore/Geometry/Geometry.h"

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

  /// Configures and returns a particle filter
  std::unique_ptr<util::PositionInVolumeFilter> CreateParticleVolumeFilter
      (std::set<std::string> const& vol_names) const;

  std::map<unsigned int, simb::MCParticle> _trackid_to_mcparticle;

  std::unique_ptr<util::PositionInVolumeFilter> _part_filter;

  std::string _mcparticle_producer = "largeant";
  std::string _mctrack_producer = "mcreco";
  std::string _mcshower_producer = "mcreco";

  std::vector<std::string> _overburden_volumes = {"volShieldingLid", "volShieldingTop", "volMezzanineLid"};

  double _x_max; //!< x-max of volume box used to determine whether to save track information
  double _x_min; //!< x-min of volume box used to determine whether to save track information
  double _y_max; //!< y-max of volume box used to determine whether to save track information
  double _y_min; //!< y-min of volume box used to determine whether to save track information
  double _z_max; //!< z-max of volume box used to determine whether to save track information
  double _z_min; //!< z-min of volume box used to determine whether to save track information

  TTree* _tree;

  int _run, _subrun, _event;

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
  std::vector<float> _mcs_ancestor_pdg;
  std::vector<std::string> _mcs_ancestor_process;
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
};


OverburdenAna::OverburdenAna(fhicl::ParameterSet const& p)
  : EDProducer{p} 
{
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("OBTree","");

  _tree->Branch("run",     &_run,     "run/I");
  _tree->Branch("subrun",     &_subrun,     "subrun/I");
  _tree->Branch("event",     &_event,     "event/I");

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
  _tree->Branch("mcs_mother_process", "std::vector<std::string>>", &_mcs_mother_process);
  _tree->Branch("mcs_ancestor_process", "std::vector<std::string>>", &_mcs_ancestor_process);
  _tree->Branch("mcs_mother_in_ob_trackid", "std::vector<float>", &_mcs_mother_in_ob_trackid);
  _tree->Branch("mcs_n_lt1",     &_n_mcs_lt1,     "mcs_n_lt1/I");
  _tree->Branch("mcs_n_lt1_from_ob",     &_n_mcs_lt1_from_ob,     "mcs_n_lt1_from_ob/I");

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

    bool in_det = InDetector(mcp);

    // Only save the MCP if it's a primary, or if it crosses the det
    if (mcp->Process() == "primary" || in_det) { 
    
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
    bool start_in_det = InDetector(mc_shower->Start().X(), mc_shower->Start().Y(), mc_shower->Start().Z());
    bool end_in_det = InDetector(mc_shower->End().X(), mc_shower->End().Y(), mc_shower->End().Z());
    if (!(start_in_det || end_in_det)) {
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
    _mcs_ancestor_pdg.push_back(mc_shower->AncestorPdgCode());
    _mcs_ancestor_process.push_back(mc_shower->AncestorProcess());
    
    _mcs_mother_in_ob_trackid.push_back(mother_in_ob);
   
    
  }
   

   _tree->Fill();
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
  _mcs_ancestor_pdg.clear();
  _mcs_ancestor_process.clear();
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



DEFINE_ART_MODULE(OverburdenAna)
