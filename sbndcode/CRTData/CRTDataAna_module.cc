////////////////////////////////////////////////////////////////////////
// Class:       CRTDataAna
// Plugin Type: producer (Unknown Unknown)
// File:        CRTDataAna_module.cc
//
// Generated at Wed Oct 20 11:57:57 2021 by Marco Del Tutto using cetskelgen
// from  version .
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
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

class CRTDataAna;


class CRTDataAna : public art::EDProducer {
public:
  explicit CRTDataAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTDataAna(CRTDataAna const&) = delete;
  CRTDataAna(CRTDataAna&&) = delete;
  CRTDataAna& operator=(CRTDataAna const&) = delete;
  CRTDataAna& operator=(CRTDataAna&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  virtual void beginSubRun(art::SubRun & sr) override;

  using Point_t = std::array<double, 3>;

private:

  /// Clear vectors
  void clear_vectors();

  /// Check if the point is in the detector
  bool InDetector(const double& x, const double& y, const double& z) const;

  /// Check if an MCP passes through the detector, sets step to the step index when particle is in detector
  bool InDetector(const art::Ptr<simb::MCParticle>, int & step);

  std::string _mctruth_producer = "generator"; // For storing POT information
  std::string _mcparticle_producer = "largeant";

  std::map<unsigned int, simb::MCParticle> _trackid_to_mcparticle;


  double _x_max; //!< x-max of volume box used to determine whether to save track information
  double _x_min; //!< x-min of volume box used to determine whether to save track information
  double _y_max; //!< y-max of volume box used to determine whether to save track information
  double _y_min; //!< y-min of volume box used to determine whether to save track information
  double _z_max; //!< z-max of volume box used to determine whether to save track information
  double _z_min; //!< z-min of volume box used to determine whether to save track information


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
  int _nu_pip_mult; ///< Pi0 multiplicity
  int _nu_pi0_mult; ///< Pi plus multiplicity
  int _nu_p_mult; ///< Proton multiplicity
  std::vector<int> _pars_pdg; ///< All other particles produced - pdg code
  std::vector<float> _pars_e; ///< All other particles produced - energy

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
  std::vector<float> _mcp_intpc_e;

  TTree* _sr_tree; ///< TTree filled per subrun
  int _sr_run; ///< Subrun Run number
  int _sr_subrun; ///< Subrun Subrun number
  double _sr_begintime; ///< Subrun start time
  double _sr_endtime; ///< Subrun end time
  double _sr_pot; ///< Subrun POT

};


CRTDataAna::CRTDataAna(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
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
  _tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
  _tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
  _tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");
  _tree->Branch("pars_pdg", "std::vector<int>", &_pars_pdg);
  _tree->Branch("pars_e", "std::vector<float>", &_pars_e);

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
  _tree->Branch("mcp_intpc_e", "std::vector<float>", &_mcp_intpc_e);

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
  std::cout << "\tx_max: " << _x_max << std::endl;
  std::cout << "\tx_min: " << _x_min << std::endl;
  std::cout << "\ty_max: " << _y_max << std::endl;
  std::cout << "\ty_min: " << _y_min << std::endl;
  std::cout << "\tz_max: " << _z_max << std::endl;
  std::cout << "\tz_min: " << _z_min << std::endl;
}

void CRTDataAna::produce(art::Event& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

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

      // // Do not save neutrinos interacting in the detector if we are using a dirt sample
      // if (_simulating_dirt && InDetector(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z)) {
      //   return;
      // }

      _pars_pdg.clear();
      _pars_e.clear();
      _nu_pip_mult = 0;
      _nu_pi0_mult = 0;
      _nu_p_mult = 0; 

      for (int p = 0; p < mct_v[i]->NParticles(); p++) {
        auto const & mcp = mct_v[i]->GetParticle(p);

        if (mcp.StatusCode() != 1) continue;

        _pars_pdg.push_back(mcp.PdgCode());
        _pars_e.push_back(mcp.E());

        if (mcp.PdgCode() == 111) {
          _nu_pi0_mult++;
          std::cout << "There is a GENIE pi0 with energy " << mcp.E() << std::endl;
        } else if (std::abs(mcp.PdgCode()) == 211) {
          _nu_pip_mult++;
        }
        else if (std::abs(mcp.PdgCode()) == 2112) {
          _nu_p_mult++;
        }
      }
      std::cout << "_nu_e " << _nu_e << std::endl;
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
    _tree->Fill();
    return;
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

      _mcp_cross_crt = IntersectsCRTs(mcp);

      int step = 0;
      bool in_det = InDetector(mcp, step);
      _mcp_intpc.push_back(in_det);
      if (in_det) {
        _mcp_intpc_e.push_back(mcp->E(step));
      } else {
        _mcp_intpc_e.push_back(-9999.);
      }
    }
  }

}

bool CRTDataAna::IntersectsCRTs(art::Ptr<simb::MCParticle> mcp) {

  auto int_up = GetIntersection(TVector3(mcp->Vx(),mcp->Vy(),mcp->Vz()),
                                TVector3(mcp->Px(),mcp->Py(),mcp->Pz()),
                                -200);

  auto int_down = GetIntersection(TVector3(mcp->Vx(),mcp->Vy(),mcp->Vz()),
                                  TVector3(mcp->Px(),mcp->Py(),mcp->Pz()),
                                  -200);

}


TVector3 CRTDataAna::GetIntersection(TVector3 pos, TVector3 dir, float z_location) {

  TVector3 plane_point(0, 0, z_location);
  TVector3 plane_normal(0, 0, 1);

  TVector3 diff = pos - plane_point;
  double prod1 = diff.Dot(plane_normal);
  double prod2 = dir.Dot(plane_normal);
  double prod3 = prod1 / prod2;
  return pos - dir * prod3;

}


bool CRTDataAna::InDetector(const double& x,
                            const double& y,
                            const double& z) const {
  return !( x > _x_max || x < _x_min ||
            z > _z_max || z < _z_min ||
            y > _y_max || y < _y_min );
}

bool CRTDataAna::InDetector(art::Ptr<simb::MCParticle> mcp, int & step) {
  auto t = mcp->Trajectory();
  for (size_t i = 0; i < t.size(); i++) {
    if (InDetector(t.X(i), t.Y(i), t.Z(i))) {
      step = i - 1;
      if (step < 0) {
        step = 0;
      }
      return true;
    }
  }
  return false;
}

void CRTDataAna::clear_vectors() {


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
  _mcp_intpc_e.clear();

}


void CRTDataAna::beginSubRun(art::SubRun & sr) {

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


DEFINE_ART_MODULE(CRTDataAna)
