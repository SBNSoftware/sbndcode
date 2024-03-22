////////////////////////////////////////////////////////////////////////
// Class:       SBNDOpT0Finder
// Plugin Type: producer (art v3_05_01)
// File:        SBNDOpT0Finder_module.cc
//
// Authors: Marco del Tutto and Lynn Tung 

// Flash Matching using charge-to-light likelihood scoring 

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
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <algorithm>
#include <cmath>



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

  void DoMatch(art::Event& e,
               int tpc,
               std::unique_ptr<std::vector<sbn::OpT0Finder>> & t0_v,
               std::unique_ptr< art::Assns<recob::Slice, sbn::OpT0Finder>> & slice_t0_assn_v,
               std::unique_ptr< art::Assns<recob::OpFlash, sbn::OpT0Finder>> & flash_t0_assn_v);

  /// Constructs all the LightClusters (TPC Objects) in a specified TPC
  bool ConstructLightClusters(art::Event& e, unsigned int tpc);

  /// Convert from a list of PDS names to a list of op channels
  std::vector<int> PDNamesToList(std::vector<std::string> pd_names);

  /// return if "xarapucas" are included in the list of photodetectors
  bool UseArapucas(std::vector<std::string> pd_names);

  /// returns if a channel corresponds to PMT (0) or X-Arapuca (1)
  std::vector<int> GetChannelTypes(int nopdets);

  /// Returns a list of uncoated PMTs that are a subset of those in ch_to_use
  std::vector<int> GetUncoatedPMTList(std::vector<int> ch_to_use);

  std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;

  ::flashmatch::FlashMatchManager _mgr; ///< The flash matching manager
  std::vector<flashmatch::FlashMatch_t> _result_v; ///< Matching result will be stored here

  std::vector<std::string> _opflash_producer_v; ///< The OpFlash producers (to be set)
  std::vector<std::string> _opflash_ara_producer_v;
  bool _use_arapucas;
  std::vector<unsigned int> _tpc_v; ///< TPC number per OpFlash producer (to be set)
  std::string _slice_producer; ///< The Slice producer (to be set)
  std::string _trk_producer;
  std::string _shw_producer;
  std::string _calo_producer;

  double _flash_trange_start; ///< The time start from where to include flashes (to be set)
  double _flash_trange_end; ///< The time stop from where to stop including flashes (to be set)

  bool _select_nus;
  bool _collection_only;
  std::vector<float> _cal_area_const; 
  std::vector<int>   _opch_to_skip;
  float _dQdx_limit;
  float _pitch_limit;

  bool  _exclude_exiting;
  bool  _track_const_conv;
  bool  _shower_const_conv; 

  float _track_to_photons; ///< The conversion factor betweeen hit integral and photons (to be set)
  float _shower_to_photons; ///< The conversion factor betweeen hit integral and photons (to be set)

  std::vector<std::string> _photo_detectors; ///< The photodetector to use (to be set)
  std::vector<int> _opch_to_use; ///< List of opch to use (will be infered from _photo_detectors)
  std::vector<int> _opch_types; ///< List of opch types, 0 for PMT and 1 for xARAPUCAs, -1 for other 
  std::vector<int> _uncoated_pmts; ///< List of uncoated opch to use (will be infered from _opch_to_use)
  std::vector<geo::Point_t> _opch_centers; ///< List of opch cneter coordinates 

  opdet::sbndPDMapAlg _pds_map; ///< map for photon detector types

  std::vector<flashmatch::QCluster_t> _light_cluster_v; ///< Vector that contains all the TPC objects

  std::map<int, art::Ptr<recob::Slice>> _clusterid_to_slice; /// Will contain map tpc object id -> Slice
  std::map<int, art::Ptr<recob::OpFlash>> _flashid_to_opflash; /// Will contain map flash id -> OpFlash

  TTree* _tree1;
  int _run, _subrun, _event;
  int _tpc;
  int _matchid, _flashid, _tpcid, _sliceid, _pfpid; 
  double _t0, _score;
  double _hypo_pe, _flash_pe;
  std::vector<double> _flash_spec;
  std::vector<double> _hypo_spec;
  int _nopdets_masked;

  TTree* _tree2;
  std::vector<float> _dep_x, _dep_y, _dep_z, _dep_E, _dep_charge, _dep_photons, _dep_pitch;
  std::vector<int> _dep_slice;
  std::vector<int> _dep_pfpid;
  std::vector<int> _dep_trk;
};


SBNDOpT0Finder::SBNDOpT0Finder(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  produces<std::vector<sbn::OpT0Finder>>();
  produces<art::Assns<recob::Slice, sbn::OpT0Finder>>();
  produces<art::Assns<recob::OpFlash, sbn::OpT0Finder>>();

  ::art::ServiceHandle<geo::Geometry> geo;

  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _semi_model = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);

  _opflash_producer_v =     p.get<std::vector<std::string>>("OpFlashProducers");
  _opflash_ara_producer_v = p.get<std::vector<std::string>>("OpFlashAraProducers");
  _tpc_v = p.get<std::vector<unsigned int>>("TPCs");
  _slice_producer = p.get<std::string>("SliceProducer");
  _trk_producer   = p.get<std::string>("TrackProducer");
  _shw_producer   = p.get<std::string>("ShowerProducer");
  _calo_producer  = p.get<std::string>("CaloProducer");

  _flash_trange_start = p.get<double>("FlashVetoTimeStart", 0);
  _flash_trange_end = p.get<double>("FlashVetoTimeEnd", 2);

  _photo_detectors = p.get<std::vector<std::string>>("PhotoDetectors");
  _use_arapucas = this->UseArapucas(_photo_detectors);
  _opch_to_use = this->PDNamesToList(_photo_detectors);
  _opch_types  = this->GetChannelTypes(geo->NOpDets());
  _uncoated_pmts = this->GetUncoatedPMTList(_opch_to_use);

  _select_nus        = p.get<bool>("SelectNeutrino");
  _collection_only   = p.get<bool>("CollectionPlaneOnly");
  _cal_area_const    = p.get<std::vector<float>>("CalAreaConstants");
  _opch_to_skip      = p.get<std::vector<int>>("OpChannelsToSkip");
  _dQdx_limit        = p.get<float>("dQdxLimit");
  _pitch_limit       = p.get<float>("PitchLimit");

  _exclude_exiting   = p.get<bool>("ExcludeExitingOpDets");
  _track_const_conv  = p.get<bool>("TrackConstantConversion");
  _shower_const_conv = p.get<bool>("ShowerConstantConversion");

  _track_to_photons  = p.get<float>("ChargeToNPhotonsTrack");
  _shower_to_photons = p.get<float>("ChargeToNPhotonsShower");


  if (_tpc_v.size() != _opflash_producer_v.size()) {
    throw cet::exception("SBNDOpT0Finder")
      << "TPC vector and OpFlash producer vector don't have the same size, check your fcl params.";
  }

  _mgr.Configure(p.get<flashmatch::Config_t>("FlashMatchConfig"));

  _mgr.SetSemiAnalyticalModel(std::move(_semi_model));

  _flash_spec.resize(geo->NOpDets(), 0.);
  _hypo_spec.resize(geo->NOpDets(), 0.);
  _opch_centers.resize(geo->NOpDets());

  art::ServiceHandle<art::TFileService> fs;

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
  _tree2->Branch("t0",              &_t0,                              "t0/D");
  _tree2->Branch("score",           &_score,                           "score/D");
  _tree2->Branch("hypo_pe",         &_hypo_pe,                         "hypo_pe/D");
  _tree2->Branch("flash_pe",        &_flash_pe,                        "flash_pe/D");
  _tree2->Branch("hypo_spec",       "std::vector<double>",             &_hypo_spec);
  _tree2->Branch("flash_spec",      "std::vector<double>",             &_flash_spec);
  _tree2->Branch("nopdets_masked",  &_nopdets_masked,                  "nmaked_opdets/I");
}

void SBNDOpT0Finder::produce(art::Event& e)
{
  std::unique_ptr<std::vector<sbn::OpT0Finder>> opt0_result_v (new std::vector<sbn::OpT0Finder>);
  std::unique_ptr< art::Assns<recob::Slice, sbn::OpT0Finder>> slice_opt0_assn_v (new art::Assns<recob::Slice, sbn::OpT0Finder>);
  std::unique_ptr< art::Assns<recob::OpFlash, sbn::OpT0Finder>> flash_opt0_assn_v (new art::Assns<recob::OpFlash, sbn::OpT0Finder>);

  // set default masks at the beginning of every event 
  _mgr.SetChannelMask(_opch_to_use);
  _mgr.SetChannelType(_opch_types);
  _mgr.SetUncoatedPMTs(_uncoated_pmts);
  // _mgr.PrintConfig();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

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
    DoMatch(e, tpc, opt0_result_v, slice_opt0_assn_v, flash_opt0_assn_v);

  }

  // Finally, place the anab::T0 vector and the associations in the Event
  e.put(std::move(opt0_result_v));
  e.put(std::move(slice_opt0_assn_v));
  e.put(std::move(flash_opt0_assn_v));

  return;
}

void SBNDOpT0Finder::DoMatch(art::Event& e,
                             int tpc,
                             std::unique_ptr<std::vector<sbn::OpT0Finder>> & opt0_result_v,
                             std::unique_ptr< art::Assns<recob::Slice, sbn::OpT0Finder>> & slice_opt0_assn_v,
                             std::unique_ptr< art::Assns<recob::OpFlash, sbn::OpT0Finder>> & flash_opt0_assn_v) {

  _flashid_to_opflash.clear();
  _clusterid_to_slice.clear();
  _tpc = tpc;
  _sliceid = -1;
  _pfpid = -1; 
  _flashid = -1;
  _t0 = -9999;
  _score = -1; 
  _hypo_pe = -1;
  _flash_pe = -1;
  std::fill(_hypo_spec.begin(), _hypo_spec.end(), 0);
  std::fill(_flash_spec.begin(), _flash_spec.end(), 0);
  _nopdets_masked=0;

  auto const & flash_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[tpc]);
  if(!flash_h.isValid() || flash_h->empty()) {
    mf::LogInfo("SBNDOpT0Finder") << "Don't have good flashes from producer "
                                  << _opflash_producer_v[tpc] << std::endl;
    _matchid = -1; 
    _tree2->Fill();
    return;
  }

  // Construct the vector of OpFlashes
  std::vector<art::Ptr<recob::OpFlash>> flash_pmt_v;
  art::fill_ptr_vector(flash_pmt_v, flash_h);
  
  // if using arapucas: 
  std::vector<art::Ptr<recob::OpFlash>> flash_ara_v;
  if (_use_arapucas){
    auto const & flash_ara_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_ara_producer_v[tpc]);
    if(!flash_ara_h.isValid() || flash_ara_h->empty()) {
      mf::LogInfo("SBNDOpT0Finder") << "Don't have good flashes from producer "
                                    << _opflash_ara_producer_v[tpc] << std::endl;
      _matchid = -1; 
      _tree2->Fill();
      return;
    }
    art::fill_ptr_vector(flash_ara_v, flash_ara_h);
  }

  ::art::ServiceHandle<geo::Geometry> geo;

  int n_flashes = 0;
  std::vector<::flashmatch::Flash_t> all_flashes;

  std::vector<recob::OpFlash> flash_comb_v;
  std::vector<bool> combine_v;

  for (size_t i=0; i < flash_pmt_v.size(); i++){
    auto const& flash_pmt = *flash_pmt_v[i];
    // combine xara + pmt PE information
    if (_use_arapucas){
      std::vector<double> combined_pe(geo->NOpDets(), 0.0); 
      bool combine = false;
      for (size_t j=0; j < flash_ara_v.size(); j++){
        auto const& flash_ara = *flash_ara_v[j];
        // if the ara and pmt flashes match: 
        if (abs(flash_pmt.Time() - flash_ara.Time()) < 0.05){
          combine = true;
          if (flash_pmt.Time() > _flash_trange_start  && flash_pmt.Time() < _flash_trange_end)
            mf::LogInfo("SBNDOpT0Finder") << "Combining PMT OpFlash (time: " 
                                          << flash_pmt.Time() 
                                          << "), with ARA OpFlash (time: " 
                                          << flash_ara.Time() << ")" << std::endl;
          // add the arapuca flash PE to the pmt flash PE (the PEs vector are different sizes for PMT and xARAPUCAs)
          for(unsigned int pmt_ch = 0; pmt_ch < flash_pmt.PEs().size(); pmt_ch++)
            combined_pe.at(pmt_ch) += flash_pmt.PEs().at(pmt_ch);
          for(unsigned int ara_ch = 0; ara_ch < flash_ara.PEs().size(); ara_ch++)
            combined_pe.at(ara_ch) += flash_ara.PEs().at(ara_ch);
  
          // create new flash with combined PE information and pmt flash information
          recob::OpFlash new_flash(flash_pmt.Time(), flash_pmt.TimeWidth(), flash_pmt.AbsTime(),
            flash_pmt.Frame(), combined_pe, flash_pmt.InBeamFrame(), flash_pmt.OnBeamTime(), 
            flash_pmt.FastToTotal(), flash_pmt.XCenter(), flash_pmt.XWidth(), 
            flash_pmt.YCenter(), flash_pmt.YWidth(), flash_pmt.ZCenter(), flash_pmt.ZWidth());
        
          flash_comb_v.push_back(new_flash);
          break;
        }
      }
      // if no arapuca flashes are found
      if (combine == false){
        flash_comb_v.push_back(flash_pmt);
        if (flash_pmt.Time() > _flash_trange_start  && flash_pmt.Time() < _flash_trange_end)
          mf::LogInfo("SBNDOpT0Finder") << "Unable to combine XARAPUCA OpFlash with PMT OpFlash... Using PMT Only" << std::endl;
      }
      combine_v.push_back(combine);
    }
  }
  int nflashes_tot = (_use_arapucas)? flash_comb_v.size():flash_pmt_v.size(); 
  auto xara_opch = PDNamesToList({"xarapuca_vis","xarapuca_vuv"}); 

  for (int n = 0; n < nflashes_tot; n++) {

    auto const& flash = (_use_arapucas)? flash_comb_v.at(n) : *flash_pmt_v[n];

    mf::LogDebug("SBNDOpT0Finder") << "Flash time from " << _opflash_producer_v[tpc] << ": " << flash.Time() << std::endl;

    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }

    _flashid_to_opflash[n_flashes] = flash_pmt_v[n];

    n_flashes++;

    // Construct a Flash_t
    ::flashmatch::Flash_t f;
    f.x = f.x_err = 0;
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    f.pds_mask_v.resize(geo->NOpDets(), 0);

    for (unsigned int op_ch = 0; op_ch < f.pe_v.size(); op_ch++) {
      bool skip = std::find(_opch_to_skip.begin(), _opch_to_skip.end(), op_ch) != _opch_to_skip.end(); 
      bool skip_ara = ((_use_arapucas == false) && (std::find(xara_opch.begin(), xara_opch.end(), op_ch) != xara_opch.end()));
      bool skip_combine = ((_use_arapucas) && (combine_v.at(n) == false) && (std::find(xara_opch.begin(), xara_opch.end(), op_ch) != xara_opch.end())); 
      if (skip || skip_ara || skip_combine ){
        f.pds_mask_v.at(op_ch) = 1;
        f.pe_v[op_ch] = 0.;
        f.pe_err_v[op_ch] = 0.;
      }
      else {
        f.pe_v[op_ch] = flash.PE(op_ch);
        f.pe_err_v[op_ch] = sqrt(flash.PE(op_ch));
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
    _matchid = -2; 
    _tree2->Fill();
    return;
  }

  // Fill vector of opch centers 
  for (size_t opch=0; opch < geo->NOpDets(); opch++){
    _opch_centers[opch] = geo->OpDetGeoFromOpChannel(opch).GetCenter();
  }

  // Get all the light clusters
  // auto light_cluster_v = GetLighClusters(e);
  if (!ConstructLightClusters(e, tpc)) {
    mf::LogInfo("SBNDOpT0Finder") << "Cannot construct Light Clusters." << std::endl;
    _matchid = -3; 
    _tree2->Fill();
    return;
  }

  // Don't waste time if there are no clusters
  if (!_light_cluster_v.size()) {
    mf::LogInfo("SBNDOpT0Finder") << "No slices to work with in TPC " << tpc << "." << std::endl;
    _matchid = -3; 
    _tree2->Fill();
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

    // Get the matched flash time, the t0
    auto const& flash = _mgr.FlashArray()[_flashid];
    auto const& cluster = _mgr.QClusterArray()[_tpcid];
    _t0 = flash.time;

    // Save the reconstructed flash and hypothesis flash PE spectrum
    if(_hypo_spec.size() != match.hypothesis.size()) {
      throw cet::exception("SBNDOpT0Finder") << "Hypothesis size mismatch!";
    }

    _nopdets_masked = 0;

    std::vector<int> result_opch_v(geo->NOpDets(), 0);
    for(size_t opch=0; opch<_hypo_spec.size(); ++opch){
      if (flash.pds_mask_v.at(opch)!=0 || cluster.tpc_mask_v.at(opch)!=0){
        _hypo_spec[opch] = 0;
        _flash_spec[opch] = 0;
        _nopdets_masked++;
      }
      else{
        _hypo_spec[opch]  = match.hypothesis[opch];
        _flash_spec[opch] = flash.pe_v[opch];
        result_opch_v[opch] = 1;
      }
    }
    // Also save the total number of photoelectrons
    _flash_pe = 0.;
    _hypo_pe  = 0.;
    for(auto const& v : _hypo_spec) _hypo_pe += v;
    for(auto const& v : _flash_spec) _flash_pe += v;

    mf::LogInfo("SBNDOpT0Finder") << "Matched TPC object " << _tpcid
                                  << " with flash number " << _flashid
                                  << " in TPC " << tpc
                                  << " at time " << _t0
                                  << " -> score: " << _score << std::endl;

    sbn::OpT0Finder opt0_result(tpc, _t0, _score, _flash_pe, _hypo_pe,
                                _flash_spec, _hypo_spec, result_opch_v);

    opt0_result_v->push_back(opt0_result);
    util::CreateAssn(*this, e, *opt0_result_v, _clusterid_to_slice[_tpcid], *slice_opt0_assn_v);
    util::CreateAssn(*this, e, *opt0_result_v, _flashid_to_opflash[_flashid], *flash_opt0_assn_v);

    art::Ptr<recob::Slice> ptr_slice = _clusterid_to_slice[_tpcid];
    int slice_id = ptr_slice->ID();
    _sliceid = slice_id;

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
        _pfpid = pfp->Self();
      }
    }
    _tree2->Fill();
  }

}

bool SBNDOpT0Finder::ConstructLightClusters(art::Event& e, unsigned int tpc) {
  // One slice is one QCluster_t.
  // Start from a slice, get all the PFParticles, from there get all the spacepoints, from
  // there get all the hits on the collection plane.
  // Use the charge on the collection plane to estimate the light, and the 3D spacepoint
  // position for the 3D location.

  _light_cluster_v.clear();

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
  art::ServiceHandle<sim::LArG4Parameters const> g4param;
  ::art::ServiceHandle<geo::Geometry> geo;

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
  // For using track calorimetry objects, get slice->pfp->track->calo 
  art::FindManyP<recob::Track>  pfp_to_trks (pfp_h, e, _trk_producer);
  art::FindManyP<anab::Calorimetry> trk_to_calo (trk_h, e, _calo_producer);
  // For using track constant objects 
  art::FindManyP<recob::SpacePoint> trk_to_spacepoints(trk_h, e, _trk_producer);
  // Get spacepoint to Shower
  art::FindManyP<recob::Shower> pfp_to_shws (pfp_h, e, _shw_producer);
  art::FindManyP<recob::SpacePoint> shw_to_spacepoints(shw_h, e, _shw_producer);
  art::FindManyP<recob::Hit> spacepoint_to_hits (spacepoint_h, e, _slice_producer);

  // Loop over the Slices
  for (size_t n_slice = 0; n_slice < slice_h->size(); n_slice++) {
    flashmatch::QCluster_t light_cluster;
    light_cluster.tpc_mask_v.resize(geo->NOpDets(), 0);

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

    std::vector<int> exit_opch; // mask of opch near the exit point for uncontained tracks

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

          // ** exiting track section ** 
          // find if the track is uncontained and intersects the wire planes
          bool uncontained = false;
          auto const trk_start = track->Start();
          auto const trk_end   = track->End();
          geo::Point_t exit_pt; 

          if (abs(trk_start.X()) >= 2.0*geo->DetHalfWidth()-3.0) {exit_pt = trk_start; uncontained = true;}
          else if (abs(trk_end.X()) >= 2.0*geo->DetHalfWidth()-3.0) {exit_pt = trk_end; uncontained = true;}
          if (uncontained && _exclude_exiting){
            mf::LogInfo("SBNDOpT0Finder") << "Found particle with exit point: " 
                                          << exit_pt.X() << ", " 
                                          << exit_pt.Y() << ", " 
                                          << exit_pt.Z() << std::endl;

            int tpc = (exit_pt.X() > 0)? 1 : 0; 
            for (size_t opch=0; opch < geo->NOpDets(); opch++){
              if (int(opch)%2 != tpc) continue;
              // only coated PMTs and vuv arapucas will be affected by direct light
              if (_pds_map.isPDType(opch, "pmt_uncoated") || _pds_map.isPDType(opch, "xarapuca_vis")) continue;
              if (_use_arapucas && _pds_map.isPDType(opch, "xarapuca_vuv")) continue;
              auto center = _opch_centers.at(opch);
              // find which optical detectors are within range of an exiting particle
              // ** uses the projection in the beam direction ** 
              if ((abs(center.Z() - (exit_pt.Z() + 50*std::cos(track->Theta()))) <= 75) && 
                  (abs(center.Y() - (exit_pt.Y() + 50*std::cos(track->ZenithAngle()))) <= 75)){
                exit_opch.push_back(opch);
                light_cluster.tpc_mask_v.at(opch) = 1;
              }
            }
          }
          // ** end exiting section ** 

          if(!_track_const_conv){       
            // ** calo section ** 
            // access the vector from the association, **not necessarily ordered by planes** 
            std::vector<art::Ptr<anab::Calorimetry>> calo_assn_v = trk_to_calo.at(track.key());
            // fill a different vector that is **correctly ordered** by plane 
            std::vector<art::Ptr<anab::Calorimetry>> calo_v(3);
            for (auto calo : calo_assn_v){
              if (calo->PlaneID().Plane == 0) calo_v.at(0) = calo;
              if (calo->PlaneID().Plane == 1) calo_v.at(1) = calo;
              if (calo->PlaneID().Plane == 2) calo_v.at(2) = calo;
            }
            // choose the plane that we want 
            int bestPlane_trk = 2;
            if (!_collection_only){
              const unsigned int maxHits(std::max({ calo_v[0]->dEdx().size(), calo_v[1]->dEdx().size(), calo_v[2]->dEdx().size() }));
              bestPlane_trk = ((calo_v[2]->dEdx().size() == maxHits) ? 2 : (calo_v[1]->dEdx().size() == maxHits) ? 1 : (calo_v[0]->dEdx().size() == maxHits) ? 0 : -1);
              if (bestPlane_trk == -1)
                continue;
            }

            auto calo = calo_v[bestPlane_trk];
            auto dEdx_v = calo->dEdx(); // assuming units in MeV/cm
            auto dADCdx_v = calo->dQdx(); // this is in ADC/cm!!!!!!
            auto pitch_v = calo->TrkPitchVec(); // assuming units in cm 
            auto pos_v   = calo->XYZ();

            // create vector of e- instead of ADC units 
            std::vector<float> dQdx_v(dADCdx_v.size(),0);
            for (size_t n_calo = 0; n_calo < dADCdx_v.size(); n_calo++){
              dQdx_v[n_calo] = dADCdx_v[n_calo]*(1/_cal_area_const.at(bestPlane_trk));
            }

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

              double drift_time = ((2.0*geo->DetHalfWidth()) - abs(x_calo))/(det_prop.DriftVelocity()); // cm / (cm/us) 
              double atten_corr = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)

              // steps that do not contain an outlier: 
              if (pitch_v[n_calo] < _pitch_limit && dQdx_v[n_calo] < _dQdx_limit){
                pitch = pitch_v[n_calo];
                dQ = dQdx_v[n_calo] * pitch * atten_corr;
                dE = dEdx_v[n_calo] * pitch; // this value *is already* lifetime corrected
                nphotons = dE/(g4param->Wph()*1e-6) - dQ;
                trk_val = 1;
              } 
              // steps that do contain an outlier: 
              else{
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
                                        nphotons);
            } // end loop over calo steps }
          } // end calo (not using constant conversion)
          if (_track_const_conv){
            int bestPlane_trk=2; 
            std::vector<art::Ptr<recob::SpacePoint>> spacepoint_v = trk_to_spacepoints.at(track.key());
            for (size_t n_spacepoint = 0; n_spacepoint < spacepoint_v.size(); n_spacepoint++) {
              auto spacepoint = spacepoint_v[n_spacepoint];
              std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hits.at(spacepoint.key());
              if (!_collection_only){ // find best track planes if other planes are allowed 
                int nhit0_trk=0; int nhit1_trk=0; int nhit2_trk=0;
                for (size_t n_hit = 0; n_hit < hit_v.size(); n_hit++) {
                  auto hit = hit_v[n_hit];
                  if (hit->View() == 0 && hit->WireID().TPC == tpc) nhit0_trk++; 
                  if (hit->View() == 1 && hit->WireID().TPC == tpc) nhit1_trk++; 
                  if (hit->View() == 2 && hit->WireID().TPC == tpc) nhit2_trk++;
                 } 
                const int maxHits = std::max({ nhit0_trk, nhit1_trk, nhit2_trk});
                bestPlane_trk = ((nhit2_trk == maxHits) ? 2 : (nhit1_trk == maxHits) ? 1 : (nhit0_trk == maxHits) ? 0 : -1);
              } 
              for (size_t n_hit = 0; n_hit < hit_v.size(); n_hit++) {
                auto hit = hit_v[n_hit];
                // Only select hits from the collection plane/best plane and in the specified TPC
                if (hit->View() != bestPlane_trk) continue;
                if (hit->WireID().TPC != tpc) continue; 

                const auto &position(spacepoint->XYZ());
                double drift_time = ((2.0*geo->DetHalfWidth()) - abs(position[0]))/(det_prop.DriftVelocity()); // cm / (cm/us) 
                double atten_corr = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)

                const auto charge((1/_cal_area_const.at(bestPlane_trk))*hit->Integral()*atten_corr);
                double nphotons = charge*_track_to_photons;

                // Emplace this point with charge to the light cluster
                light_cluster.emplace_back(position[0],
                                          position[1],
                                          position[2],
                                          nphotons);

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
              }
            }  // End loop over Spacepoints
          } // end trk const conversion 
        } // end loop over tracks 
      } // end if pfpistrack

      else if (pfpisshower){
        std::vector<art::Ptr<recob::Shower>> shower_v = pfp_to_shws.at(pfp.key());
        for (size_t n_shw = 0; n_shw < shower_v.size(); n_shw++){
          auto shower = shower_v[n_shw];
          int bestPlane_shw=2; 
          std::vector<art::Ptr<recob::SpacePoint>> spacepoint_v = shw_to_spacepoints.at(shower.key());
          for (size_t n_spacepoint = 0; n_spacepoint < spacepoint_v.size(); n_spacepoint++) {
            auto spacepoint = spacepoint_v[n_spacepoint];
            std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hits.at(spacepoint.key());
            
            if (!_collection_only){ // find best shower plane if other planes are allowed
              bestPlane_shw = shower->best_plane();
              if ( (shower->Energy()).at(bestPlane_shw) == -999){
                int nhit0_shw=0; int nhit1_shw=0; int nhit2_shw=0;
                for (size_t n_hit=0; n_hit<hit_v.size(); n_hit++){
                  auto hit = hit_v[n_hit];
                  if (hit->View() == 0 && hit->WireID().TPC == tpc) nhit0_shw++; 
                  if (hit->View() == 1 && hit->WireID().TPC == tpc) nhit1_shw++; 
                  if (hit->View() == 2 && hit->WireID().TPC == tpc) nhit2_shw++;
                }
                const int maxHits = std::max({ nhit0_shw, nhit1_shw, nhit2_shw});
                bestPlane_shw = ((nhit2_shw == maxHits) ? 2 : (nhit1_shw == maxHits) ? 1 : (nhit0_shw == maxHits) ? 0 : -1);
              }
            }

            for (size_t n_hit = 0; n_hit < hit_v.size(); n_hit++) {
              auto hit = hit_v[n_hit];
              // Only select hits from the collection plane/best plane and in the specified TPC
              if (hit->View() != bestPlane_shw) continue;
              if (hit->WireID().TPC != tpc) continue; 

              const auto &position(spacepoint->XYZ());
              double drift_time = ((2.0*geo->DetHalfWidth()) - abs(position[0]))/(det_prop.DriftVelocity()); // cm / (cm/us) 
              double atten_corr = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)

              const auto charge((1/_cal_area_const.at(bestPlane_shw))*hit->Integral()*atten_corr);
              double nphotons = 0;
 
              if ( _shower_const_conv) nphotons = charge*_shower_to_photons;
              else{
                std::cout << "Only have shower constant conversion calculation... using constant conversion" << std::endl;
                nphotons = charge*_shower_to_photons;
              }
              // Emplace this point with charge to the light cluster
              light_cluster.emplace_back(position[0],
                                        position[1],
                                        position[2],
                                        nphotons);

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
              _dep_trk.push_back(2);
            }
          } // End loop over Spacepoints
        } // end shower loop
      } // end if pfpisshower
    } // End loop over PFParticle

    _tree1->Fill();

    // Don't include clusters with zero points
    if (!light_cluster.size()) {
      continue;
    }

    // Save the light cluster, and remember the correspondance from index to slice
    _clusterid_to_slice[_light_cluster_v.size()] = slice_v.at(n_slice);

    _light_cluster_v.emplace_back(light_cluster);

    if (!exit_opch.empty()){
      std::cout << "Not evaluating the following OpDets due to exiting particle: { ";
      for (auto opch : exit_opch)
          std::cout << opch << ' ';
      std::cout << "}\n";
    }
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
 
bool SBNDOpT0Finder::UseArapucas(std::vector<std::string> pd_names){
  bool found_ara = false;
  for (auto name : pd_names){
    if ((name == "xarapuca_vis") || ( name =="xarapuca_vis")){
      found_ara = true;
      break;
    }
  }
  return found_ara;
}

std::vector<int> SBNDOpT0Finder::GetChannelTypes(int nopdets){
  std::vector<int> out_v(nopdets,-1); 
  for (size_t ch = 0; ch < out_v.size(); ch ++){
    if (_pds_map.isPDType(ch, "pmt_uncoated") || _pds_map.isPDType(ch, "pmt_coated"))
      out_v.at(ch) = 0;
    else if (_pds_map.isPDType(ch, "xarapuca_vis") || _pds_map.isPDType(ch, "xarapuca_vuv"))
      out_v.at(ch) = 1;
  }
  return out_v;
}

std::vector<int> SBNDOpT0Finder::GetUncoatedPMTList(std::vector<int> ch_to_use) {
  std::vector<int> out_v;

  for (auto ch : ch_to_use) {
    if (_pds_map.isPDType(ch, "pmt_uncoated")) {
      out_v.push_back(ch);
    }
    else if (_pds_map.isPDType(ch, "xarapuca_vis")){
      out_v.push_back(ch);
    }
  }

  return out_v;
}



DEFINE_ART_MODULE(SBNDOpT0Finder)
