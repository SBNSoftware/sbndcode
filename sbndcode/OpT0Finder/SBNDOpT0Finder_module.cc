////////////////////////////////////////////////////////////////////////
// Class:       SBNDOpT0Finder
// Plugin Type: producer (art v3_05_01)
// File:        SBNDOpT0Finder_module.cc
//
// Generated at Sun Oct  4 17:27:36 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderTypes.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashMatchManager.h"
#include "sbncode/OpT0Finder/flashmatch/Algorithms/PhotonLibHypothesis.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/sim.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>



class SBNDOpT0Finder;


class SBNDOpT0Finder : public art::EDProducer {
public:
  explicit SBNDOpT0Finder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDOpT0Finder(SBNDOpT0Finder const&) = delete;
  SBNDOpT0Finder(SBNDOpT0Finder&&) = delete;
  SBNDOpT0Finder& operator=(SBNDOpT0Finder const&) = delete;
  SBNDOpT0Finder& operator=(SBNDOpT0Finder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  /// Performs the matching in a specified tpc
  void DoMatch(art::Event& e,
               int tpc,
               std::unique_ptr<std::vector<anab::T0>> & t0_v,
               std::unique_ptr< art::Assns<recob::Slice, anab::T0>> & slice_t0_assn_v,
               std::unique_ptr< art::Assns<recob::OpFlash, anab::T0>> & flash_t0_assn_v);

  /// Constructs all the LightClusters (TPC Objects) in a specified TPC
  bool ConstructLightClusters(art::Event& e, unsigned int tpc);

  /// Returns the number of photons given charge and PFParticle
  // float GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &pfp, art::ServiceHandle<sim::LArG4Parameters const> g4param);

  /// Convert from a list of PDS names to a list of op channels
  std::vector<int> PDNamesToList(std::vector<std::string>);

  /// Returns a list of uncoated PMTs that are a subset of those in ch_to_use
  std::vector<int> GetUncoatedPTMList(std::vector<int> ch_to_use);

  std::unique_ptr<SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;

  ::flashmatch::FlashMatchManager _mgr; ///< The flash matching manager
  std::vector<flashmatch::FlashMatch_t> _result_v; ///< Matching result will be stored here

  std::vector<std::string> _opflash_producer_v; ///< The OpFlash producers (to be set)
  std::vector<unsigned int> _tpc_v; ///< TPC number per OpFlash producer (to be set)
  std::string _slice_producer; ///< The Slice producer (to be set)
  std::string _trk_producer;
  std::string _shw_producer;
  std::string _calo_producer;

  double _flash_trange_start; ///< The time start from where to include flashes (to be set)
  double _flash_trange_end; ///< The time stop from where to stop including flashes (to be set)

  float _calibration_const;  /// conversion from (ADC*time ticks) to e- (to be set), given in units of (ADC*time ticks)/e-
  float _dQdx_limit;
  float _pitch_limit;
  bool  _exclude_outliers;

  float _track_to_photons; ///< The conversion factor betweeen hit integral and photons (to be set)
  float _shower_to_photons; ///< The conversion factor betweeen hit integral and photons (to be set)

  // bool _calc_correlated;
  // art::ServiceHandle<sim::LArG4Parameters const> g4param;

  std::vector<std::string> _photo_detectors; ///< The photodetector to use (to be set)
  std::vector<int> _opch_to_use; ///< List of opch to use (will be infered from _photo_detectors)
  std::vector<int> _uncoated_pmts; ///< List of uncoated opch to use (will be infered from _opch_to_use)

  opdet::sbndPDMapAlg _pds_map; ///< map for photon detector types
  // std::unique_ptr<opdet::sbndPDMapAlg> _pds_map;

  bool _select_nus;

  std::vector<flashmatch::QCluster_t> _light_cluster_v; ///< Vector that contains all the TPC objects

  std::map<int, art::Ptr<recob::Slice>> _clusterid_to_slice; /// Will contain map tpc object id -> Slice
  std::map<int, art::Ptr<recob::OpFlash>> _flashid_to_opflash; /// Will contain map flash id -> OpFlash

  TTree* _tree1;
  int _run, _subrun, _event;
  int _tpc;
  int _matchid, _flashid, _tpcid, _sliceid, _pfpid; 
  double _t0, _score;
  double _tpc_xmin, _qll_xmin;
  double _hypo_pe, _flash_pe;
  std::vector<double> _flash_spec;
  std::vector<double> _hypo_spec;

  TTree* _tree2;
  std::vector<float> _dep_x, _dep_y, _dep_z, _dep_E, _dep_charge, _dep_photons, _dep_pitch;
  std::vector<int> _dep_slice;
  std::vector<int> _dep_pfpid;
  std::vector<int> _dep_trk;
};


SBNDOpT0Finder::SBNDOpT0Finder(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  produces<std::vector<anab::T0>>();
  produces<art::Assns<recob::Slice, anab::T0>>();
  produces<art::Assns<recob::OpFlash, anab::T0>>();

  ::art::ServiceHandle<geo::Geometry> geo;

  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);

  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _tpc_v = p.get<std::vector<unsigned int>>("TPCs");
  _slice_producer = p.get<std::string>("SliceProducer");
  _trk_producer   = p.get<std::string>("TrackProducer");
  _shw_producer   = p.get<std::string>("ShowerProducer");
  _calo_producer  = p.get<std::string>("CaloProducer");

  _flash_trange_start = p.get<double>("FlashVetoTimeStart", 0);
  _flash_trange_end = p.get<double>("FlashVetoTimeEnd", 2);

  _photo_detectors = p.get<std::vector<std::string>>("PhotoDetectors");
  _opch_to_use = this->PDNamesToList(_photo_detectors);
  _uncoated_pmts = this->GetUncoatedPTMList(_opch_to_use);

  _calibration_const = p.get<float>("CalibrationConst");
  _dQdx_limit        = p.get<float>("dQdxLimit");
  _pitch_limit       = p.get<float>("PitchLimit");
  _exclude_outliers  = p.get<bool>("ExcludeOutliers");

  _track_to_photons = p.get<float>("ChargeToNPhotonsTrack");
  _shower_to_photons = p.get<float>("ChargeToNPhotonsShower");

  _select_nus = p.get<bool>("SelectNeutrino", true);

  if (_tpc_v.size() != _opflash_producer_v.size()) {
    throw cet::exception("SBNDOpT0Finder")
      << "TPC vector and OpFlash producer vector don't have the same size, check your fcl params.";
  }

  _mgr.Configure(p.get<flashmatch::Config_t>("FlashMatchConfig"));

  _mgr.SetChannelMask(_opch_to_use);

  _mgr.SetUncoatedPMTs(_uncoated_pmts);

  _mgr.SetSemiAnalyticalModel(std::move(_semi_model));

  _flash_spec.resize(geo->NOpDets(), 0.);
  _hypo_spec.resize(geo->NOpDets(), 0.);

  art::ServiceHandle<art::TFileService> fs;
  // art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  _tree1 = fs->make<TTree>("slice_deposition_tree","");
  _tree1->Branch("run",             &_run,                             "run/I");
  _tree1->Branch("subrun",          &_subrun,                          "subrun/I");
  _tree1->Branch("event",           &_event,                           "event/I");
  _tree1->Branch("dep_slice", "std::vector<int>", &_dep_slice);
  _tree1->Branch("dep_pfpid", "std::vector<int>", &_dep_pfpid);
  _tree1->Branch("dep_x", "std::vector<float>", &_dep_x);
  _tree1->Branch("dep_y", "std::vector<float>", &_dep_y);
  _tree1->Branch("dep_z", "std::vector<float>", &_dep_z);
  _tree1->Branch("dep_E", "std::vector<float>", &_dep_E);
  _tree1->Branch("dep_charge", "std::vector<float>", &_dep_charge);
  _tree1->Branch("dep_photons", "std::vector<float>", &_dep_photons);
  _tree1->Branch("dep_pitch"  , "std::vector<float>", &_dep_pitch);
  _tree1->Branch("dep_trk", "std::vector<int>", &_dep_trk);

  _tree2 = fs->make<TTree>("flash_match_tree","");
  _tree2->Branch("run",             &_run,                             "run/I");
  _tree2->Branch("subrun",          &_subrun,                          "subrun/I");
  _tree2->Branch("event",           &_event,                           "event/I");
  _tree2->Branch("tpc",             &_tpc,                             "tpc/I");
  _tree2->Branch("matchid",         &_matchid,                         "matchid/I");
  _tree2->Branch("tpcid",           &_tpcid,                           "tpcid/I");
  _tree2->Branch("sliceid",         &_sliceid,                         "sliceid/I");
  _tree2->Branch("pfpid",           &_pfpid,                           "pfpid/I");
  _tree2->Branch("flashid",         &_flashid,                         "flashid/I");
  _tree2->Branch("tpc_xmin",        &_tpc_xmin,                        "tpc_xmin/D");
  _tree2->Branch("qll_xmin",        &_qll_xmin,                        "qll_xmin/D");
  _tree2->Branch("t0",              &_t0,                              "t0/D");
  _tree2->Branch("score",           &_score,                           "score/D");
  _tree2->Branch("hypo_pe",         &_hypo_pe,                         "hypo_pe/D");
  _tree2->Branch("flash_pe",        &_flash_pe,                        "flash_pe/D");
  _tree2->Branch("hypo_spec",       "std::vector<double>",             &_hypo_spec);
  _tree2->Branch("flash_spec",      "std::vector<double>",             &_flash_spec);
}

void SBNDOpT0Finder::produce(art::Event& e)
{
  std::unique_ptr<std::vector<anab::T0>> t0_v (new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Slice, anab::T0>> slice_t0_assn_v (new art::Assns<recob::Slice, anab::T0>);
  std::unique_ptr< art::Assns<recob::OpFlash, anab::T0>> flash_t0_assn_v (new art::Assns<recob::OpFlash, anab::T0>);

  // _mgr.PrintConfig();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  std::cout << "run: " << _run <<  std::endl;
  std::cout << "subrun: " << _subrun  << std::endl;
  std::cout << "event: " <<  _event << std::endl;

  // Loop over the specified TPCs
  for (auto tpc : _tpc_v) {

    mf::LogInfo("SBNDOpT0Finder") << "Performing matching in TPC " << tpc << std::endl;

    // Reset the manager and the result vector
    _mgr.Reset();
    _result_v.clear();
    _tpc = tpc;

    // Tell the manager what TPC and cryostat we are going to be doing
    // the matching in. For SBND, the cryostat is always zero.
    _mgr.SetTPCCryo(tpc, 0);

    // Perform the matching in the specified TPC
    DoMatch(e, tpc, t0_v, slice_t0_assn_v, flash_t0_assn_v);

  }

  // Finally, place the anab::T0 vector and the associations in the Event
  e.put(std::move(t0_v));
  e.put(std::move(slice_t0_assn_v));
  e.put(std::move(flash_t0_assn_v));

  return;
}

void SBNDOpT0Finder::DoMatch(art::Event& e,
                             int tpc,
                             std::unique_ptr<std::vector<anab::T0>> & t0_v,
                             std::unique_ptr< art::Assns<recob::Slice, anab::T0>> & slice_t0_assn_v,
                             std::unique_ptr< art::Assns<recob::OpFlash, anab::T0>> & flash_t0_assn_v) {

  _flashid_to_opflash.clear();
  _clusterid_to_slice.clear();

  auto const & flash_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[tpc]);
  if(!flash_h.isValid() || flash_h->empty()) {
    mf::LogInfo("SBNDOpT0Finder") << "Don't have good flashes from producer "
                                  << _opflash_producer_v[tpc] << std::endl;
    return;
  }

  // Construct the vector of OpFlashes
  std::vector<art::Ptr<recob::OpFlash>> flash_v;
  art::fill_ptr_vector(flash_v, flash_h);

  ::art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<sim::LArG4Parameters const> g4param;

  int n_flashes = 0;
  std::vector<::flashmatch::Flash_t> all_flashes;

  for (size_t n = 0; n < flash_v.size(); n++) {

    auto const& flash = *flash_v[n];

    mf::LogDebug("SBNDOpT0Finder") << "Flash time from " << _opflash_producer_v[tpc] << ": " << flash.Time() << std::endl;

    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }

    _flashid_to_opflash[n_flashes] = flash_v[n];

    n_flashes++;

    // Construct a Flash_t
    ::flashmatch::Flash_t f;
    f.x = f.x_err = 0;
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int op_ch = 0; op_ch < f.pe_v.size(); op_ch++) {
      unsigned int opdet = geo->OpDetFromOpChannel(op_ch);
      if (std::find(_opch_to_use.begin(), _opch_to_use.end(), op_ch) == _opch_to_use.end()) {
        f.pe_v[opdet] = 0;
        f.pe_err_v[opdet] = 0;
      } else {
        f.pe_v[opdet] = flash.PE(op_ch);
        f.pe_err_v[opdet] = sqrt(flash.PE(op_ch));
      }
    }
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.time = flash.Time();
    f.idx = n_flashes-1;
    all_flashes.resize(n_flashes);
    all_flashes[n_flashes-1] = f;

  } // flash loop

  // Don't waste time if there are no flashes
  if (n_flashes == 0) {
    mf::LogInfo("SBNDOpT0Finder") << "Zero good flashes in this event." << std::endl;
    return;
  }

  // Get all the ligh clusters
  // auto light_cluster_v = GetLighClusters(e);
  if (!ConstructLightClusters(e, tpc)) {
    mf::LogInfo("SBNDOpT0Finder") << "Cannot construct Light Clusters." << std::endl;
    return;
  }

  // Don't waste time if there are no clusters
  if (!_light_cluster_v.size()) {
    mf::LogInfo("SBNDOpT0Finder") << "No slices to work with in TPC " << tpc << "." << std::endl;
    return;
  }

  // Emplace flashes to Flash Matching Manager
  for (auto f : all_flashes) {
    _mgr.Emplace(std::move(f));
  }

  // Emplace clusters to Flash Matching Manager
  for (auto lc : _light_cluster_v) {
    _mgr.Emplace(std::move(lc));
  }

  // Run the matching
  _result_v = _mgr.Match();

  // Loop over the matching results
  for(_matchid = 0; _matchid < (int)(_result_v.size()); ++_matchid) {

    auto const& match = _result_v[_matchid];

    _tpcid    = match.tpc_id;
    _flashid  = match.flash_id;
    _score    = match.score;
    _qll_xmin = match.tpc_point.x;

    mf::LogInfo("SBNDOpT0Finder") << "Matched TPC object " << _tpcid
                                  << " with flash number " << _flashid
                                  << " in TPC " << tpc
                                  << " -> score: " << _score
                                  << ", qll xmin: " << _qll_xmin << std::endl;

    // Get the minimum x position of the TPC Object
    _tpc_xmin = 1.e4;
    for(auto const& pt : _mgr.QClusterArray()[_tpcid]) {
      if(pt.x < _tpc_xmin) _tpc_xmin = pt.x;
    }

    // Get the matched flash time, the t0
    auto const& flash = _mgr.FlashArray()[_flashid];
    _t0 = flash.time;

    // Save the reconstructed flash and hypothesis flash PE spectrum
    if(_hypo_spec.size() != match.hypothesis.size()) {
      throw cet::exception("SBNDOpT0Finder") << "Hypothesis size mismatch!";
    }
    for(size_t pmt=0; pmt<_hypo_spec.size(); ++pmt) _hypo_spec[pmt]  = match.hypothesis[pmt];
    for(size_t pmt=0; pmt<_hypo_spec.size(); ++pmt) _flash_spec[pmt] = flash.pe_v[pmt];

    // Also save the total number of photoelectrons
    _flash_pe = 0.;
    _hypo_pe  = 0.;
    for(auto const& v : _hypo_spec) _hypo_pe += v;
    for(auto const& v : _flash_spec) _flash_pe += v;


    // Construct the anab::T0 dataproduc to put in the Event
    auto t0 = anab::T0(_t0,        // "Time": The recontructed flash time, or t0
                       _flash_pe,  // "TriggerType": placing the reconstructed total PE instead
                       _tpcid,     // "TriggerBits": placing the tpc id instead
                       _flashid,   // "ID": placing the flash id instead
                       _score);    // "TriggerConfidence": Matching score

    t0_v->push_back(t0);
    util::CreateAssn(*this, e, *t0_v, _clusterid_to_slice[_tpcid], *slice_t0_assn_v);
    util::CreateAssn(*this, e, *t0_v, _flashid_to_opflash[_flashid], *flash_t0_assn_v);

    art::Ptr<recob::Slice> ptr_slice = _clusterid_to_slice[_tpcid];
    int slice_id = ptr_slice->ID();
    _sliceid = slice_id;
    // std::cout << "slice_id: " << slice_id << std::endl;
    // std::cout << "_sliceid: " << _sliceid << std::endl;

    ::art::Handle<std::vector<recob::Slice>> slice_h;
    e.getByLabel(_slice_producer, slice_h);
    if(!slice_h.isValid() || slice_h->empty()) {
      mf::LogWarning("SBNDOpT0Finder") << "Don't have good Slices." << std::endl;
          }
    // Construct the vector of Slices
    std::vector<art::Ptr<recob::Slice>> slice_v;
    art::fill_ptr_vector(slice_v, slice_h);
    art::FindManyP<recob::PFParticle> slice_to_pfps (slice_h, e, _slice_producer);
    for (size_t n_slice = 0; n_slice < slice_h->size(); n_slice++) {
      auto slice = slice_v[n_slice];
      if (slice->ID() != _sliceid) continue;

      std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfps.at(n_slice);
      for (size_t n_pfp = 0; n_pfp < pfp_v.size(); n_pfp++) {
        auto pfp = pfp_v[n_pfp];
        if (!pfp->IsPrimary()) continue;
        // std::cout << "pfp id: " << pfp->Self() << std::endl;
        _pfpid = pfp->Self();
      }
    }
    double new_score = t0.TriggerConfidence(); 
    // std::cout << "new_score: " << new_score << std::endl;
    _score = new_score; 
    // std::cout << "score: " << _score << std::endl;
    // std::cout << "alt score:" << t0.TriggerConfidence() << std::endl;
    // _score    = t0.TriggerConfidence();

    _tree2->Fill();
  }

}

bool SBNDOpT0Finder::ConstructLightClusters(art::Event& e, unsigned int tpc) {
  // One slice is one QCluster_t.
  // Start from a slice, get all the PFParticles, from there get all the spacepoints, from
  // there get all the hits on the collection plane.
  // Use the charge on the collection plane to estimate the light, and the 3D spacepoint
  // position for the 3D location.

  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  _light_cluster_v.clear();

  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()) {
    mf::LogWarning("SBNDOpT0Finder") << "Don't have good Slices." << std::endl;
    return false;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_slice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    mf::LogWarning("SBNDOpT0Finder") << "Don't have good PFParticle." << std::endl;
    return false;
  }

  ::art::Handle<std::vector<recob::SpacePoint>> spacepoint_h;
  e.getByLabel(_slice_producer, spacepoint_h);
  if(!spacepoint_h.isValid() || spacepoint_h->empty()) {
    mf::LogWarning("SBNDOpT0Finder") << "Don't have good SpacePoint." << std::endl;
    return false;
  }

  ::art::Handle<std::vector<recob::Track>> trk_h; 
  e.getByLabel(_trk_producer, trk_h);

  ::art::Handle<std::vector<recob::Shower>> shw_h;
  e.getByLabel(_shw_producer, shw_h);

  // Construct the vector of Slices
  std::vector<art::Ptr<recob::Slice>> slice_v;
  art::fill_ptr_vector(slice_v, slice_h);

  // Get the associations between slice->pfp->spacepoint->hit
  art::FindManyP<recob::PFParticle> slice_to_pfps (slice_h, e, _slice_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoints (pfp_h, e, _slice_producer);
  art::FindManyP<recob::Hit> spacepoint_to_hits (spacepoint_h, e, _slice_producer);
  // For using track calorimetry objects, get slice->pfp->track->calo 
  art::FindManyP<recob::Track>  pfp_to_trks (pfp_h, e, _trk_producer);
  art::FindManyP<recob::Shower> pfp_to_shws (pfp_h, e, _shw_producer);
  art::FindManyP<anab::Calorimetry> trk_to_calo (trk_h, e, _calo_producer);
  // for truth information 
  art::FindManyP<recob::Hit> trk_to_hits (trk_h, e, _trk_producer);
  art::FindManyP<recob::Hit> shw_to_hits (shw_h, e, _shw_producer);

  // Loop over the Slices
  for (size_t n_slice = 0; n_slice < slice_h->size(); n_slice++) {
    // std::cout <<"n_slice: " << n_slice << std::endl;

    flashmatch::QCluster_t light_cluster;

    _dep_slice.clear();
    _dep_pfpid.clear();
    _dep_x.clear();
    _dep_y.clear();
    _dep_z.clear();
    _dep_E.clear();
    _dep_charge.clear();
    _dep_photons.clear();
    _dep_pitch.clear();
    _dep_trk.clear();

    // Get the associated PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfps.at(n_slice);

    if (_select_nus){
      bool nu_pfp = false;
      for (size_t n_pfp = 0; n_pfp < pfp_v.size(); n_pfp++) {
        auto pfp = pfp_v[n_pfp];
        unsigned pfpPDGC = std::abs(pfp->PdgCode());
        if ((pfpPDGC == 12) || (pfpPDGC == 14) || (pfpPDGC == 16))
          nu_pfp = true;
      }
      if (nu_pfp == false)
        break;
    }

    for (size_t n_pfp = 0; n_pfp < pfp_v.size(); n_pfp++) {

      auto pfp = pfp_v[n_pfp];      
      auto pfpistrack = ::lar_pandora::LArPandoraHelper::IsTrack(pfp);
      auto pfpisshower = ::lar_pandora::LArPandoraHelper::IsShower(pfp);

      if (pfpistrack){
        std::vector<art::Ptr<recob::Track>> track_v = pfp_to_trks.at(pfp.key());
        for (size_t n_trk = 0; n_trk < track_v.size(); n_trk++){
          auto track = track_v[n_trk];
          // section to analyze calo hits vs. track hits 
          std::vector<art::Ptr<anab::Calorimetry>> calo_v = trk_to_calo.at(track.key());
          std::cout << "nhits in calo plane 0: " << (calo_v[0])->dQdx().size() 
                    << ", plane 1: " << (calo_v[1])->dQdx().size()
                    << ", plane 2: " << (calo_v[2])->dQdx().size() 
                    << std::endl;

          float calo_integral = 0.0;
          // end section

          // Get the associate calo objects, we expect one for each plane... 

          // auto calo = calo_v[0]; // choose collection plane (calo object only stores pitch for collection plane)
          // choose plane with the most hits: 
          const unsigned int maxHits(std::max({ calo_v[0]->dEdx().size(), calo_v[1]->dEdx().size(), calo_v[2]->dEdx().size() }));
          const int bestPlane((calo_v[0]->dEdx().size() == maxHits) ? 0 : (calo_v[2]->dEdx().size() == maxHits) ? 2 : (calo_v[1]->dEdx().size() == maxHits) ? 1 : -1);
          if (bestPlane == -1)
            continue;
          auto calo   = calo_v[bestPlane];
          auto dEdx_v = calo->dEdx(); // assuming units in MeV/cm
          auto dADCdx_v = calo->dQdx(); // this is in ADC/cm!!!!!!
          auto pitch_v = calo->TrkPitchVec(); // assuming units in cm 
          auto pos_v   = calo->XYZ();

          // create vector of e- instead of ADC units 
          std::vector<float> dQdx_v(dADCdx_v.size(),0);
          for (size_t n_calo = 0; n_calo < dADCdx_v.size(); n_calo++){
            dQdx_v[n_calo] = dADCdx_v[n_calo]*(1/_calibration_const);
            calo_integral += dADCdx_v[n_calo]*pitch_v[n_calo]; 
          }
          // std::cout << "number of hits in calo in track: " << dADCdx_v.size() << std::endl;
          for (size_t n_calo = 0; n_calo < calo->dEdx().size(); n_calo++){
            // only select steps that are in the right TPC
            auto &position = pos_v[n_calo];
            auto x_calo = position.X();
            if ((x_calo < 0 && tpc==1 ) || (x_calo > 0 && tpc==0)) continue; // skip if not in the correct TPC 
            // variables to fill the tree
            float dE; 
            float dQ;
            float pitch; 
            float nphotons; 
            int   trk_val; 

            // lifetime correction: TODO: remove hardcoded detector lim, drift velocity, and electron lifetime 
            double drift_time = abs(200.0 - abs(x_calo))/(0.16); // in us, drift velocity = 0.16 cm/us 
            double atten_corr = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms

            // steps that do not contain an outlier: 
            if (pitch_v[n_calo] < _pitch_limit && dQdx_v[n_calo] < _dQdx_limit){
              pitch = pitch_v[n_calo];
              dQ = dQdx_v[n_calo] * pitch * atten_corr;
              dE = dEdx_v[n_calo] * pitch; // this value *is* lifetime corrected
              nphotons = dE/(19.5*1e-6) - dQ;
              trk_val = 1;
            } 
            // steps that do contain an outlier: 
            else if(_exclude_outliers == true) // skip outlier values 
              continue;
            else if(_exclude_outliers == false){ // use outlier values 
              pitch = -1.;
              dQ = dQdx_v[n_calo] * pitch_v[n_calo] * atten_corr;  
              dE = -1.;
              nphotons = dQ*_track_to_photons; 
              trk_val = 0;
            }
            // Fill tree variables 
            _dep_slice.push_back(n_slice);
            _dep_pfpid.push_back(pfp->Self());
            _dep_x.push_back(position.X());
            _dep_y.push_back(position.Y());
            _dep_z.push_back(position.Z());
            _dep_E.push_back(dE);
            _dep_charge.push_back(dQ);
            _dep_photons.push_back(nphotons);
            _dep_pitch.push_back(pitch);
            _dep_trk.push_back(trk_val);

            // emplace this point into the light cluster 
            light_cluster.emplace_back(position.X(),
                                      position.Y(),
                                      position.Z(),
                                      nphotons,
                                      trk_val);
          } // end loop over calo steps 
          // get truth information 
          std::vector<art::Ptr<recob::Hit>> hit_v = trk_to_hits.at(track.key());
          int n_hit_plane0 =0;
          int n_hit_plane1 =0;
          int n_hit_plane2 =0;
          // int hit_best_plane;

          // double hit_integral = 0.;
          for (size_t n_hit=0; n_hit<hit_v.size(); n_hit++){
            auto hit = hit_v[n_hit];
            if (hit->View() == 0)  n_hit_plane0++;
            if (hit->View() == 1)  n_hit_plane1++;
            if (hit->View() == 2)  n_hit_plane2++;
            // if (hit->View() == bestPlane){
            //   hit_integral += hit->Integral();
            // }
          }
          std::cout << "nhits in trk plane 0: " << n_hit_plane0
                    << ", plane 1: " << n_hit_plane1
                    << ", plane 2: " << n_hit_plane2
                    << std::endl;    

          // std::cout << "calo integral tot: " << calo_integral << std::endl;    
          // std::cout << "hits integral tot: " << hit_integral << std::endl;
          // std::cout << "ratio of calo/hit: " << calo_integral/hit_integral << std::endl;
          // std::cout << "ratio of hit/calo: " << hit_integral/calo_integral << std::endl;
        } // end loop over tracks 
      } // end calo track section 
      else if (pfpisshower){
        // look inside shower
        std::vector<art::Ptr<recob::Shower>> shower_v = pfp_to_shws.at(pfp.key());
        for (size_t n_shw = 0; n_shw < shower_v.size(); n_shw++){
          std::cout << "shower section start" << std::endl;
          auto shower = shower_v[n_shw];
          int bestPlane = shower->best_plane();
          const std::vector<double> showerEnergy = shower->Energy();
          std::cout << "bestPlane: " << bestPlane << std::endl;
          std::cout << "energy: (" << showerEnergy[0] << ", " << showerEnergy[1] << ", " << showerEnergy[2] << ")" << std::endl;
          std::vector<art::Ptr<recob::Hit>> hit_v = shw_to_hits.at(shower.key());
          int nhit0=0; int nhit1=0; int nhit2=0;
          for (size_t n_hit=0; n_hit<hit_v.size(); n_hit++){
            auto hit = hit_v[n_hit];
            if (hit->View() == 0)  nhit0++;
            if (hit->View() == 1)  nhit1++;
            if (hit->View() == 2)  nhit2++;
          }
          const std::vector<int> shw_nhit{nhit0,nhit1,nhit2};
          std::cout << "nhits in shw: (" << shw_nhit[0] << ", " << shw_nhit[1] << ", " << shw_nhit[2] << ")" << std::endl;
          std::cout << "avg energy per hit (bestplane): " << showerEnergy[bestPlane] / shw_nhit[bestPlane] << std::endl;
          std::cout << "shower section end" << std::endl;
        }
      }
      else{ // if pfp is not a track (belongs to shower or neither track nor shower)
        // Get the associated SpacePoints
        std::vector<art::Ptr<recob::SpacePoint>> spacepoint_v = pfp_to_spacepoints.at(pfp.key());

        for (size_t n_spacepoint = 0; n_spacepoint < spacepoint_v.size(); n_spacepoint++) {
          auto spacepoint = spacepoint_v[n_spacepoint];

          // Get the associated hits
          std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hits.at(spacepoint.key());

          for (size_t n_hit = 0; n_hit < hit_v.size(); n_hit++) {

            auto hit = hit_v[n_hit];

            // Only select hits from the collection plane
            if (hit->View() != geo::kZ) {
              continue;
            }

            // Only use hits (and so spacepoints) that are in the specified TPC
            if (hit->WireID().TPC != tpc) {
              continue;
            }

            const auto &position(spacepoint->XYZ());
            double drift_time = abs(200.0 - abs(position[0]))/(0.16); // in us, drift velocity = 0.16 cm/us 
            double atten_corr = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms

            const auto charge((1/_calibration_const)*hit->Integral()*atten_corr);
            auto nphotons = charge*_shower_to_photons;
            // obtain truth info
            float ide_energy = 0.0;
            float ide_charge = 0.0;
            if ((charge > 0) & (hit->Channel() != 8294) & (hit->Channel() != 8327)){
              std::vector<const sim::IDE*> ide_v = bt_serv->HitToSimIDEs_Ps(clockData, hit);
              for (auto const ide: ide_v){
                ide_energy += ide->energy;
                ide_charge += ide->numElectrons; // electrons at the wires 
              }
            }

            // std::cout << "ide energy: " << ide_energy << std::endl;
            // std::cout << "ide charge: " << ide_charge*atten_corr << std::endl;
            // std::cout << "ide photon: " << ide_energy/(19.5*1e-6) - ide_charge*atten_corr << std::endl;
            
            // if (nphotons<0) std::cout << "nphotons val: " << nphotons << std::endl;
            // Emplace this point with charge to the light cluster
            light_cluster.emplace_back(position[0],
                                       position[1],
                                       position[2],
                                       nphotons,
                                       0);

            // Also save the quantites for the output tree
            _dep_slice.push_back(n_slice);
            _dep_pfpid.push_back(pfp->Self());
            _dep_x.push_back(position[0]);
            _dep_y.push_back(position[1]);
            _dep_z.push_back(position[2]);
            _dep_E.push_back(-1.);
            _dep_charge.push_back(charge);
            _dep_photons.push_back(nphotons);
            _dep_pitch.push_back(-1.);
            _dep_trk.push_back(0);
            
            // _tru_charge.push_back(ide_charge*atten_corr);
            // _tru_photon.push_back(ide_energy/(19.5*1e-6) - ide_charge*atten_corr);
            // std::cout << "true charge/dep charge: " << (ide_charge*atten_corr)/charge << std::endl;
            // std::cout << "true photons/dep photons: " << (ide_energy/(19.5*1e-6) - ide_charge*atten_corr)/nphotons << std::endl;
          }
        } // End loop over Spacepoints
      } // end non-track loop
    } // End loop over PFParticle

    _tree1->Fill();

    // Don't include clusters with zero points
    if (!light_cluster.size()) {
      continue;
    }

    // Save the light cluster, and remember the correspondance from index to slice
    _clusterid_to_slice[_light_cluster_v.size()] = slice_v.at(n_slice);
    _light_cluster_v.emplace_back(light_cluster);

  } // End loop over Slices

  return true;
}

std::vector<int> SBNDOpT0Finder::PDNamesToList(std::vector<std::string> pd_names) {

  std::vector<int> out_ch_v;

  for (auto name : pd_names) {
    auto ch_v = _pds_map.getChannelsOfType(name);
    out_ch_v.insert(out_ch_v.end(), ch_v.begin(), ch_v.end());
  }

  return out_ch_v;

}

std::vector<int> SBNDOpT0Finder::GetUncoatedPTMList(std::vector<int> ch_to_use) {
  std::vector<int> out_v;

  for (auto ch : ch_to_use) {
    if (_pds_map.isPDType(ch, "pmt_uncoated")) {
      out_v.push_back(ch);
    }
  }

  return out_v;
}



DEFINE_ART_MODULE(SBNDOpT0Finder)


