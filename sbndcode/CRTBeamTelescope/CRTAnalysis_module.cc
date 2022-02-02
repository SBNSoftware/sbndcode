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
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "larcoreobj/SummaryData/POTSummary.h"

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

  std::string _mctruth_label;
  std::string _crthit_label;
  std::string _crttrack_label;

  TTree* _tree;

  int _run, _subrun, _event;

  float _mct_sp_pdg; ///< Single particle PDG
  float _mct_sp_e; ///< Single particle energy
  float _mct_sp_px; ///< Single particle momentum X
  float _mct_sp_py; ///< Single particle momentum Y
  float _mct_sp_pz; ///< Single particle momentum Z
  float _mct_sp_vx; ///< Single particle vertex X
  float _mct_sp_vy; ///< Single particle vertex Y
  float _mct_sp_vz; ///< Single particle vertex Z

  std::vector<double> _chit_x; ///< CRT hit x
  std::vector<double> _chit_y; ///< CRT hit y
  std::vector<double> _chit_z; ///< CRT hit z
  std::vector<double> _chit_time; ///< CRT hit time
  std::vector<double> _chit_pes; ///< CRT hit PEs
  std::vector<int> _chit_plane; ///< CRT hit plane

  std::vector<double> _ct_time; ///< CRT track time
  std::vector<double> _ct_pes; ///< CRT track PEs
  std::vector<double> _ct_x1; ///< CRT track x1
  std::vector<double> _ct_y1; ///< CRT track y1
  std::vector<double> _ct_z1; ///< CRT track z1
  std::vector<double> _ct_x2; ///< CRT track x2
  std::vector<double> _ct_y2; ///< CRT track y2
  std::vector<double> _ct_z2; ///< CRT track z2


  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun

};


CRTAnalysis::CRTAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _mctruth_label = p.get<std::string>("MCTruthLabel", "generator");
  _crthit_label = p.get<std::string>("CRTHitLabel", "crthit");
  _crttrack_label = p.get<std::string>("CRTTrackLabel", "crttrack");

  art::ServiceHandle<art::TFileService> fs;

  _tree = fs->make<TTree>("tree","");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");

  _tree->Branch("mct_sp_pdg", &_mct_sp_pdg, "mct_sp_pdg/F");
  _tree->Branch("mct_sp_e", &_mct_sp_e, "mct_sp_e/F");
  _tree->Branch("mct_sp_px", &_mct_sp_px, "mct_sp_px/F");
  _tree->Branch("mct_sp_py", &_mct_sp_py, "mct_sp_py/F");
  _tree->Branch("mct_sp_pz", &_mct_sp_pz, "mct_sp_pz/F");
  _tree->Branch("mct_sp_vx", &_mct_sp_vx, "mct_sp_vx/F");
  _tree->Branch("mct_sp_vy", &_mct_sp_vy, "mct_sp_vy/F");
  _tree->Branch("mct_sp_vz", &_mct_sp_vz, "mct_sp_vz/F");

  _tree->Branch("chit_x", "std::vector<double>", &_chit_x);
  _tree->Branch("chit_y", "std::vector<double>", &_chit_y);
  _tree->Branch("chit_z", "std::vector<double>", &_chit_z);
  _tree->Branch("chit_time", "std::vector<double>", &_chit_time);
  _tree->Branch("chit_pes", "std::vector<double>", &_chit_pes);
  _tree->Branch("chit_plane", "std::vector<int>", &_chit_plane);

  _tree->Branch("ct_time", "std::vector<double>", &_ct_time);
  _tree->Branch("ct_pes", "std::vector<double>", &_ct_pes);
  _tree->Branch("ct_x1", "std::vector<double>", &_ct_x1);
  _tree->Branch("ct_y1", "std::vector<double>", &_ct_y1);
  _tree->Branch("ct_z1", "std::vector<double>", &_ct_z1);
  _tree->Branch("ct_x2", "std::vector<double>", &_ct_x2);
  _tree->Branch("ct_y2", "std::vector<double>", &_ct_y2);
  _tree->Branch("ct_z2", "std::vector<double>", &_ct_z2);


  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void CRTAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  //
  // Get the MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_label, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product " << _mctruth_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

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


  //
  // Fill the MCTruth in the tree
  //
  assert(mct_v.size() == 1);
  auto mct = mct_v[0];
  if (mct->Origin() == simb::kBeamNeutrino) {

  }
  else if (mct->Origin() == simb::kSingleParticle) {
    assert(mct->NParticles() == 1);
    auto particle = mct->GetParticle(0);
    _mct_sp_pdg = particle.PdgCode();
    _mct_sp_e = particle.E();
    _mct_sp_px = particle.Px();
    _mct_sp_py = particle.Py();
    _mct_sp_pz = particle.Pz();
    _mct_sp_vx = particle.Vx();
    _mct_sp_vy = particle.Vy();
    _mct_sp_vz = particle.Vz();

  }

  //
  // Fill the CRT Hits in the tree
  //
  size_t n_hits = crt_hit_v.size();
  _chit_x.resize(n_hits);
  _chit_y.resize(n_hits);
  _chit_z.resize(n_hits);
  _chit_time.resize(n_hits);
  _chit_pes.resize(n_hits);
  _chit_plane.resize(n_hits);

  for (size_t i = 0; i < n_hits; i++) {

    auto hit = crt_hit_v[i];

    _chit_x[i] = hit->x_pos;
    _chit_y[i] = hit->y_pos;
    _chit_z[i] = hit->z_pos;
    _chit_time[i] = hit->ts1_ns * 0.001;
    _chit_pes[i] = hit->peshit;

    if (hit->tagger == "volTaggerNorth_0") {
      _chit_plane[i] = 0;
    } else {
      _chit_plane[i] = 1;
    }

  }

  //
  // Fill the CRT Tracks in the tree
  //
  size_t n_tracks = crt_track_v.size();
  _ct_pes.resize(n_tracks);
  _ct_time.resize(n_tracks);
  _ct_x1.resize(n_tracks);
  _ct_y1.resize(n_tracks);
  _ct_z1.resize(n_tracks);
  _ct_x2.resize(n_tracks);
  _ct_y2.resize(n_tracks);
  _ct_z2.resize(n_tracks);

  for (size_t i = 0; i < n_tracks; ++i){

    auto track = crt_track_v[i];

    _ct_pes[i] = track->peshit;
    _ct_time[i] = track->ts1_ns * 0.001;
    _ct_x1[i] = track->x1_pos;
    _ct_y1[i] = track->y1_pos;
    _ct_z1[i] = track->z1_pos;
    _ct_x2[i] = track->x2_pos;
    _ct_y2[i] = track->y2_pos;
    _ct_z2[i] = track->z2_pos;
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
  sr.getByLabel(_mctruth_label, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }
  std::cout << "POT for this subrun: " << _sr_pot << std::endl;

  _sr_tree->Fill();

}

DEFINE_ART_MODULE(CRTAnalysis)
