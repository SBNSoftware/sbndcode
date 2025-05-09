////////////////////////////////////////////////////////////////////////
// Class:       CRTAnalysis
// Plugin Type: analyzer
// File:        CRTAnalysis_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
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

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/ChannelMaps/CRT/CRTChannelMapService.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

namespace sbnd::crt {
  class CRTAnalysis;
}

class sbnd::crt::CRTAnalysis : public art::EDAnalyzer {
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

  void AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec);

  void AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec);

  void AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec);

  void AnalyseSimDeposits(std::vector<art::Ptr<sim::AuxDetSimChannel>> &SimDepositVec);

  void AnalyseFEBDatas(std::vector<art::Ptr<FEBData>> &FEBDataVec);

  void AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec);

  void AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec, const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints,
                          const art::FindManyP<CRTStripHit> &clustersToStripHits);

  void AnalyseTrueDepositsPerTagger(const std::map<CRTBackTrackerAlg::Category, bool> &recoStatusMap);

  void AnalyseTrueDeposits(const std::map<int, std::pair<bool, bool>> &recoStatusMap);

  void AnalyseCRTTracks(const art::Event &e, const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec, const art::FindManyP<CRTSpacePoint> &tracksToSpacePoints,
                        const art::FindOneP<CRTCluster> &spacePointsToClusters, const art::FindManyP<CRTStripHit> &clustersToStripHits);

  void AnalyseTPCMatching(const art::Event &e, const art::Handle<std::vector<recob::Track>> &TPCTrackHandle,
                          const art::Handle<std::vector<CRTSpacePoint>> &CRTSpacePointHandle, const art::Handle<std::vector<CRTCluster>> &CRTClusterHandle,
                          const art::Handle<std::vector<recob::PFParticle>> &PFPHandle,
                          const std::map<CRTBackTrackerAlg::Category, bool> &spacePointRecoStatusMap = {},
                          const std::map<int, std::pair<bool, bool>> &trackRecoStatusMap = {});

  double SPTimeDelta(std::vector<std::pair<int, double>> &time_set);

private:

  CRTGeoAlg fCRTGeoAlg;
  TPCGeoAlg fTPCGeoAlg;
  CRTBackTrackerAlg fCRTBackTrackerAlg;

  std::string fMCParticleModuleLabel, fSimDepositModuleLabel, fFEBDataModuleLabel, fCRTStripHitModuleLabel,
    fCRTClusterModuleLabel, fCRTSpacePointModuleLabel, fCRTTrackModuleLabel, fTPCTrackModuleLabel,
    fCRTSpacePointMatchingModuleLabel, fCRTTrackMatchingModuleLabel, fPFPModuleLabel, fPTBModuleLabel,
    fTDCModuleLabel, fTimingReferenceModuleLabel;
  bool fDebug, fDataMode, fNoTPC, fHasPTB, fHasTDC, fTruthMatch;
  //! Adding some of the reco parameters to save corrections
  double fPEAttenuation, fTimeWalkNorm, fTimeWalkScale, fPropDelay;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;
  int _crt_timing_reference_type;
  int _crt_timing_reference_channel;

  //mc truth
  std::vector<int16_t>              _mc_trackid;
  std::vector<int16_t>              _mc_pdg;
  std::vector<int16_t>              _mc_status;
  std::vector<uint16_t>             _mc_ndaughters;
  std::vector<std::vector<int16_t>> _mc_daughters;
  std::vector<double>               _mc_vx;
  std::vector<double>               _mc_vy;
  std::vector<double>               _mc_vz;
  std::vector<double>               _mc_vt;
  std::vector<double>               _mc_endx;
  std::vector<double>               _mc_endy;
  std::vector<double>               _mc_endz;
  std::vector<double>               _mc_endt;
  std::vector<double>               _mc_vpx;
  std::vector<double>               _mc_vpy;
  std::vector<double>               _mc_vpz;
  std::vector<double>               _mc_ve;
  std::vector<double>               _mc_endpx;
  std::vector<double>               _mc_endpy;
  std::vector<double>               _mc_endpz;
  std::vector<double>               _mc_ende;

  //G4 detector id
  std::vector<int16_t> _ide_trackid;
  std::vector<float>   _ide_e;
  std::vector<float>   _ide_entryx;
  std::vector<float>   _ide_entryy;
  std::vector<float>   _ide_entryz;
  std::vector<float>   _ide_entryt;
  std::vector<float>   _ide_exitx;
  std::vector<float>   _ide_exity;
  std::vector<float>   _ide_exitz;
  std::vector<float>   _ide_exitt;

  //front end mother board
  std::vector<uint16_t>              _feb_mac5;
  std::vector<int16_t>               _feb_tagger;
  std::vector<uint16_t>              _feb_flags;
  std::vector<uint32_t>              _feb_ts0;
  std::vector<uint32_t>              _feb_ts1;
  std::vector<uint32_t>              _feb_unixs;
  std::vector<std::vector<uint16_t>> _feb_adc;
  std::vector<uint32_t>              _feb_coinc;

  //strip hit to select the strip which has ADC above threshold
  std::vector<uint32_t> _sh_channel;
  std::vector<int16_t>  _sh_tagger;
  std::vector<double>   _sh_ts0;
  std::vector<double>   _sh_ts1;
  std::vector<uint32_t> _sh_unixs;
  std::vector<double>   _sh_pos;
  std::vector<double>   _sh_err;
  std::vector<uint16_t> _sh_adc1;
  std::vector<uint16_t> _sh_adc2;
  std::vector<bool>     _sh_saturated1;
  std::vector<bool>     _sh_saturated2;
  std::vector<int>      _sh_truth_trackid;
  std::vector<double>   _sh_truth_completeness;
  std::vector<double>   _sh_truth_purity;
  std::vector<double>   _sh_truth_pos;
  std::vector<double>   _sh_truth_energy;
  std::vector<double>   _sh_truth_time;

  //cluster from x-y coincidence for CRTSpacePoint , this is what we normally call a CRT hit
  std::vector<double>                _cl_ts0;
  std::vector<double>                _cl_ts1;
  std::vector<uint32_t>              _cl_unixs;
  std::vector<uint16_t>              _cl_nhits;
  std::vector<int16_t>               _cl_tagger;
  std::vector<uint8_t>               _cl_composition;
  std::vector<std::vector<uint32_t>> _cl_channel_set;
  std::vector<std::vector<uint16_t>> _cl_adc_set;
  std::vector<std::vector<double>>   _cl_sh_ts0_set; //! To store t0 from x-y coincidences
  std::vector<std::vector<double>>   _cl_sh_ts1_set; //! To store t1 from x-y coincidences
  std::vector<std::vector<uint16_t>> _cl_sh_feb_mac5_set; //! MAC5 addresses of StripHit FEBs
  std::vector<std::vector<double>>   _cl_sh_time_walk_set; //! Time walk correction
  std::vector<std::vector<double>>   _cl_sh_prop_delay_set; //! Light propagation correction
  std::vector<int>                   _cl_truth_trackid;
  std::vector<double>                _cl_truth_completeness;
  std::vector<double>                _cl_truth_purity;
  std::vector<double>                _cl_truth_hit_completeness;
  std::vector<double>                _cl_truth_hit_purity;
  std::vector<int>                   _cl_truth_pdg;
  std::vector<double>                _cl_truth_energy;
  std::vector<double>                _cl_truth_time;
  std::vector<double>                _cl_truth_x;
  std::vector<double>                _cl_truth_y;
  std::vector<double>                _cl_truth_z;
  std::vector<double>                _cl_truth_core_energy;
  std::vector<double>                _cl_truth_core_time;
  std::vector<double>                _cl_truth_core_x;
  std::vector<double>                _cl_truth_core_y;
  std::vector<double>                _cl_truth_core_z;
  std::vector<bool>                  _cl_has_sp;
  std::vector<double>                _cl_sp_x;
  std::vector<double>                _cl_sp_ex;
  std::vector<double>                _cl_sp_y;
  std::vector<double>                _cl_sp_ey;
  std::vector<double>                _cl_sp_z;
  std::vector<double>                _cl_sp_ez;
  std::vector<double>                _cl_sp_pe;
  std::vector<double>                _cl_sp_ts0;
  std::vector<double>                _cl_sp_ets0;
  std::vector<double>                _cl_sp_ts1;
  std::vector<double>                _cl_sp_ets1;
  std::vector<double>                _cl_sp_dts0;
  std::vector<double>                _cl_sp_dts1;
  std::vector<bool>                  _cl_sp_complete;

  //backtrack truth information from reco level
  std::vector<int>     _td_tag_trackid;
  std::vector<int>     _td_tag_pdg;
  std::vector<int16_t> _td_tag_tagger;
  std::vector<double>  _td_tag_energy;
  std::vector<double>  _td_tag_time;
  std::vector<double>  _td_tag_x;
  std::vector<double>  _td_tag_y;
  std::vector<double>  _td_tag_z;
  std::vector<bool>    _td_tag_reco_status;

  std::vector<int>     _td_trackid;
  std::vector<int>     _td_pdg;
  std::vector<double>  _td_energy;
  std::vector<double>  _td_time;
  std::vector<bool>    _td_reconstructable;
  std::vector<bool>    _td_reco_status;
  std::vector<bool>    _td_reco_triple;

  //track level information
  std::vector<double>                _tr_start_x;
  std::vector<double>                _tr_start_y;
  std::vector<double>                _tr_start_z;
  std::vector<double>                _tr_end_x;
  std::vector<double>                _tr_end_y;
  std::vector<double>                _tr_end_z;
  std::vector<double>                _tr_dir_x;
  std::vector<double>                _tr_dir_y;
  std::vector<double>                _tr_dir_z;
  std::vector<double>                _tr_ts0;
  std::vector<double>                _tr_ets0;
  std::vector<double>                _tr_ts1;
  std::vector<double>                _tr_ets1;
  std::vector<double>                _tr_pe;
  std::vector<double>                _tr_length;
  std::vector<double>                _tr_tof;
  std::vector<double>                _tr_theta;
  std::vector<double>                _tr_phi;
  std::vector<bool>                  _tr_triple;
  std::vector<int16_t>               _tr_tagger1;
  std::vector<int16_t>               _tr_tagger2;
  std::vector<int16_t>               _tr_tagger3;
  std::vector<std::vector<uint32_t>> _tr_channel_set;
  std::vector<std::vector<uint16_t>> _tr_adc_set;
  std::vector<int>                   _tr_truth_trackid;
  std::vector<double>                _tr_truth_completeness;
  std::vector<double>                _tr_truth_purity;
  std::vector<int>                   _tr_truth_pdg;
  std::vector<double>                _tr_truth_energy;
  std::vector<double>                _tr_truth_time;
  std::vector<double>                _tr_truth_start_x;
  std::vector<double>                _tr_truth_start_y;
  std::vector<double>                _tr_truth_start_z;
  std::vector<double>                _tr_truth_end_x;
  std::vector<double>                _tr_truth_end_y;
  std::vector<double>                _tr_truth_end_z;
  std::vector<double>                _tr_truth_dir_x;
  std::vector<double>                _tr_truth_dir_y;
  std::vector<double>                _tr_truth_dir_z;
  std::vector<double>                _tr_truth_particle_energy;
  std::vector<double>                _tr_truth_length;
  std::vector<double>                _tr_truth_tof;
  std::vector<double>                _tr_truth_theta;
  std::vector<double>                _tr_truth_phi;

  std::vector<double>                _tpc_start_x;
  std::vector<double>                _tpc_start_y;
  std::vector<double>                _tpc_start_z;
  std::vector<double>                _tpc_end_x;
  std::vector<double>                _tpc_end_y;
  std::vector<double>                _tpc_end_z;
  std::vector<double>                _tpc_start_dir_x;
  std::vector<double>                _tpc_start_dir_y;
  std::vector<double>                _tpc_start_dir_z;
  std::vector<double>                _tpc_end_dir_x;
  std::vector<double>                _tpc_end_dir_y;
  std::vector<double>                _tpc_end_dir_z;
  std::vector<double>                _tpc_length;
  std::vector<double>                _tpc_track_score;
  std::vector<int>                   _tpc_truth_trackid;
  std::vector<int>                   _tpc_truth_pdg;
  std::vector<double>                _tpc_truth_energy;
  std::vector<double>                _tpc_truth_time;
  std::vector<bool>                  _tpc_sp_matchable;
  std::vector<bool>                  _tpc_sp_matched;
  std::vector<std::vector<uint32_t>> _tpc_sp_channel_set;
  std::vector<bool>                  _tpc_sp_good_match;
  std::vector<double>                _tpc_sp_xshift;
  std::vector<double>                _tpc_sp_ts0;
  std::vector<double>                _tpc_sp_ts1;
  std::vector<int>                   _tpc_sp_tagger;
  std::vector<int>                   _tpc_sp_nhits;
  std::vector<double>                _tpc_sp_x;
  std::vector<double>                _tpc_sp_y;
  std::vector<double>                _tpc_sp_z;
  std::vector<double>                _tpc_sp_score;
  std::vector<bool>                  _tpc_tr_matchable;
  std::vector<bool>                  _tpc_tr_matched;
  std::vector<bool>                  _tpc_tr_good_match;
  std::vector<double>                _tpc_tr_ts0;
  std::vector<double>                _tpc_tr_ts1;
  std::vector<std::vector<int>>      _tpc_tr_taggers;
  std::vector<double>                _tpc_tr_start_x;
  std::vector<double>                _tpc_tr_start_y;
  std::vector<double>                _tpc_tr_start_z;
  std::vector<double>                _tpc_tr_end_x;
  std::vector<double>                _tpc_tr_end_y;
  std::vector<double>                _tpc_tr_end_z;
  std::vector<double>                _tpc_tr_score;

  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;

  std::vector<uint64_t> _ptb_llt_trigger;
  std::vector<uint64_t> _ptb_llt_timestamp;

  std::vector<uint32_t>    _tdc_channel;
  std::vector<uint64_t>    _tdc_timestamp;
  std::vector<uint64_t>    _tdc_offset;
  std::vector<std::string> _tdc_name;
};

sbnd::crt::CRTAnalysis::CRTAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  , fCRTBackTrackerAlg(p.get<fhicl::ParameterSet>("CRTBackTrackerAlg", fhicl::ParameterSet()))
{
  fMCParticleModuleLabel            = p.get<std::string>("MCParticleModuleLabel", "largeant");
  fSimDepositModuleLabel            = p.get<std::string>("SimDepositModuleLabel", "genericcrt");
  fFEBDataModuleLabel               = p.get<std::string>("FEBDataModuleLabel", "crtsim");
  fCRTStripHitModuleLabel           = p.get<std::string>("CRTStripHitModuleLabel", "crtstrips");
  fCRTClusterModuleLabel            = p.get<std::string>("CRTClusterModuleLabel", "crtclustering");
  fCRTSpacePointModuleLabel         = p.get<std::string>("CRTSpacePointModuleLabel", "crtspacepoints");
  fCRTTrackModuleLabel              = p.get<std::string>("CRTTrackModuleLabel", "crttracks");
  fTPCTrackModuleLabel              = p.get<std::string>("TPCTrackModuleLabel", "pandoraSCETrack");
  fCRTSpacePointMatchingModuleLabel = p.get<std::string>("CRTSpacePointMatchingModuleLabel", "crtspacepointmatchingSCE");
  fCRTTrackMatchingModuleLabel      = p.get<std::string>("CRTTrackMatchingModuleLabel", "crttrackmatchingSCE");
  fPFPModuleLabel                   = p.get<std::string>("PFPModuleLabel", "pandora");
  fPTBModuleLabel                   = p.get<std::string>("PTBModuleLabel", "ptbdecoder");
  fTDCModuleLabel                   = p.get<std::string>("TDCModuleLabel", "tdcdecoder");
  fTimingReferenceModuleLabel       = p.get<std::string>("TimingReferenceModuleLabel", "crtstrips");
  fDebug                            = p.get<bool>("Debug", false);
  fDataMode                         = p.get<bool>("DataMode", false);
  fNoTPC                            = p.get<bool>("NoTPC", false);
  fHasPTB                           = p.get<bool>("HasPTB", false);
  fHasTDC                           = p.get<bool>("HasTDC", false);
  fTruthMatch                       = p.get<bool>("TruthMatch", true);
  //! Adding some of the reco parameters to save corrections
  fPEAttenuation                    = p.get<double>("PEAttenuation", 1.0);
  fTimeWalkNorm                     = p.get<double>("TimeWalkNorm",  0.0);
  fTimeWalkScale                    = p.get<double>("TimeWalkScale", 0.0);
  fPropDelay                        = p.get<double>("PropDelay",     0.0);

  if(!fDataMode && fTruthMatch)
    fCRTBackTrackerAlg = CRTBackTrackerAlg(p.get<fhicl::ParameterSet>("CRTBackTrackerAlg", fhicl::ParameterSet()));

  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);
  fTree->Branch("event", &_event);
  fTree->Branch("crt_timing_reference_type", &_crt_timing_reference_type);
  fTree->Branch("crt_timing_reference_channel", &_crt_timing_reference_channel);

  if(!fDataMode)
    {
      fTree->Branch("mc_trackid", "std::vector<int16_t>", &_mc_trackid);
      fTree->Branch("mc_pdg", "std::vector<int16_t>", &_mc_pdg);
      fTree->Branch("mc_status", "std::vector<int16_t>", &_mc_status);
      fTree->Branch("mc_ndaughters", "std::vector<uint16_t>", &_mc_ndaughters);
      fTree->Branch("mc_daughters", "std::vector<std::vector<int16_t>>", &_mc_daughters);
      fTree->Branch("mc_vx", "std::vector<double>", &_mc_vx);
      fTree->Branch("mc_vy", "std::vector<double>", &_mc_vy);
      fTree->Branch("mc_vz", "std::vector<double>", &_mc_vz);
      fTree->Branch("mc_vt", "std::vector<double>", &_mc_vt);
      fTree->Branch("mc_endx", "std::vector<double>", &_mc_endx);
      fTree->Branch("mc_endy", "std::vector<double>", &_mc_endy);
      fTree->Branch("mc_endz", "std::vector<double>", &_mc_endz);
      fTree->Branch("mc_endt", "std::vector<double>", &_mc_endt);
      fTree->Branch("mc_vpx", "std::vector<double>", &_mc_vpx);
      fTree->Branch("mc_vpy", "std::vector<double>", &_mc_vpy);
      fTree->Branch("mc_vpz", "std::vector<double>", &_mc_vpz);
      fTree->Branch("mc_ve", "std::vector<double>", &_mc_ve);
      fTree->Branch("mc_endpx", "std::vector<double>", &_mc_endpx);
      fTree->Branch("mc_endpy", "std::vector<double>", &_mc_endpy);
      fTree->Branch("mc_endpz", "std::vector<double>", &_mc_endpz);
      fTree->Branch("mc_ende", "std::vector<double>", &_mc_ende);

      fTree->Branch("ide_trackid", "std::vector<int16_t>", &_ide_trackid);
      fTree->Branch("ide_e", "std::vector<float>", &_ide_e);
      fTree->Branch("ide_entryx", "std::vector<float>", &_ide_entryx);
      fTree->Branch("ide_entryy", "std::vector<float>", &_ide_entryy);
      fTree->Branch("ide_entryz", "std::vector<float>", &_ide_entryz);
      fTree->Branch("ide_entryt", "std::vector<float>", &_ide_entryt);
      fTree->Branch("ide_exitx", "std::vector<float>", &_ide_exitx);
      fTree->Branch("ide_exity", "std::vector<float>", &_ide_exity);
      fTree->Branch("ide_exitz", "std::vector<float>", &_ide_exitz);
      fTree->Branch("ide_exitt", "std::vector<float>", &_ide_exitt);
    }

  fTree->Branch("feb_mac5", "std::vector<uint16_t>", &_feb_mac5);
  fTree->Branch("feb_tagger", "std::vector<int16_t>", &_feb_tagger);
  fTree->Branch("feb_flags", "std::vector<uint16_t>", &_feb_flags);
  fTree->Branch("feb_ts0", "std::vector<uint32_t>", &_feb_ts0);
  fTree->Branch("feb_ts1", "std::vector<uint32_t>", &_feb_ts1);
  fTree->Branch("feb_unixs", "std::vector<uint32_t>", &_feb_unixs);
  fTree->Branch("feb_adc", "std::vector<std::vector<uint16_t>>", &_feb_adc);
  fTree->Branch("feb_coinc", "std::vector<uint32_t>", &_feb_coinc);

  fTree->Branch("sh_channel", "std::vector<uint32_t>", &_sh_channel);
  fTree->Branch("sh_tagger", "std::vector<int16_t>", &_sh_tagger);
  fTree->Branch("sh_ts0", "std::vector<double>", &_sh_ts0);
  fTree->Branch("sh_ts1", "std::vector<double>", &_sh_ts1);
  fTree->Branch("sh_unixs", "std::vector<uint32_t>", &_sh_unixs);
  fTree->Branch("sh_pos", "std::vector<double>", &_sh_pos);
  fTree->Branch("sh_err", "std::vector<double>", &_sh_err);
  fTree->Branch("sh_adc1", "std::vector<uint16_t>", &_sh_adc1);
  fTree->Branch("sh_adc2", "std::vector<uint16_t>", &_sh_adc2);
  fTree->Branch("sh_saturated1", "std::vector<bool>", &_sh_saturated1);
  fTree->Branch("sh_saturated2", "std::vector<bool>", &_sh_saturated2);
  if(!fDataMode && fTruthMatch)
    {
      fTree->Branch("sh_truth_trackid", "std::vector<int>", &_sh_truth_trackid);
      fTree->Branch("sh_truth_completeness", "std::vector<double>", &_sh_truth_completeness);
      fTree->Branch("sh_truth_purity", "std::vector<double>", &_sh_truth_purity);
      fTree->Branch("sh_truth_pos", "std::vector<double>", &_sh_truth_pos);
      fTree->Branch("sh_truth_energy", "std::vector<double>", &_sh_truth_energy);
      fTree->Branch("sh_truth_time", "std::vector<double>", &_sh_truth_time);
    }

  fTree->Branch("cl_ts0", "std::vector<double>", &_cl_ts0);
  fTree->Branch("cl_ts1", "std::vector<double>", &_cl_ts1);
  fTree->Branch("cl_unixs", "std::vector<uint32_t>", &_cl_unixs);
  fTree->Branch("cl_nhits", "std::vector<uint16_t>", &_cl_nhits);
  fTree->Branch("cl_tagger", "std::vector<int16_t>", &_cl_tagger);
  fTree->Branch("cl_composition", "std::vector<uint8_t>", &_cl_composition);
  fTree->Branch("cl_channel_set", "std::vector<std::vector<uint32_t>>", &_cl_channel_set);
  fTree->Branch("cl_adc_set", "std::vector<std::vector<uint16_t>>", &_cl_adc_set);
  fTree->Branch("cl_sh_ts0_set", "std::vector<std::vector<double>>", &_cl_sh_ts0_set);
  fTree->Branch("cl_sh_ts1_set", "std::vector<std::vector<double>>", &_cl_sh_ts1_set);
  fTree->Branch("cl_sh_feb_mac5_set", "std::vector<std::vector<uint16_t>>", &_cl_sh_feb_mac5_set);
  fTree->Branch("cl_sh_time_walk_set", "std::vector<std::vector<double>>", &_cl_sh_time_walk_set);
  fTree->Branch("cl_sh_prop_delay_set", "std::vector<std::vector<double>>", &_cl_sh_prop_delay_set);
  if(!fDataMode && fTruthMatch)
    {
      fTree->Branch("cl_truth_trackid", "std::vector<int>", &_cl_truth_trackid);
      fTree->Branch("cl_truth_completeness", "std::vector<double>", &_cl_truth_completeness);
      fTree->Branch("cl_truth_purity", "std::vector<double>", &_cl_truth_purity);
      fTree->Branch("cl_truth_hit_completeness", "std::vector<double>", &_cl_truth_hit_completeness);
      fTree->Branch("cl_truth_hit_purity", "std::vector<double>", &_cl_truth_hit_purity);
      fTree->Branch("cl_truth_pdg", "std::vector<int>", &_cl_truth_pdg);
      fTree->Branch("cl_truth_energy", "std::vector<double>", &_cl_truth_energy);
      fTree->Branch("cl_truth_time", "std::vector<double>", &_cl_truth_time);
      fTree->Branch("cl_truth_x", "std::vector<double>", &_cl_truth_x);
      fTree->Branch("cl_truth_y", "std::vector<double>", &_cl_truth_y);
      fTree->Branch("cl_truth_z", "std::vector<double>", &_cl_truth_z);
      fTree->Branch("cl_truth_core_energy", "std::vector<double>", &_cl_truth_core_energy);
      fTree->Branch("cl_truth_core_time", "std::vector<double>", &_cl_truth_core_time);
      fTree->Branch("cl_truth_core_x", "std::vector<double>", &_cl_truth_core_x);
      fTree->Branch("cl_truth_core_y", "std::vector<double>", &_cl_truth_core_y);
      fTree->Branch("cl_truth_core_z", "std::vector<double>", &_cl_truth_core_z);
    }
  fTree->Branch("cl_has_sp", "std::vector<bool>", &_cl_has_sp);
  fTree->Branch("cl_sp_x", "std::vector<double>", &_cl_sp_x);
  fTree->Branch("cl_sp_ex", "std::vector<double>", &_cl_sp_ex);
  fTree->Branch("cl_sp_y", "std::vector<double>", &_cl_sp_y);
  fTree->Branch("cl_sp_ey", "std::vector<double>", &_cl_sp_ey);
  fTree->Branch("cl_sp_z", "std::vector<double>", &_cl_sp_z);
  fTree->Branch("cl_sp_ez", "std::vector<double>", &_cl_sp_ez);
  fTree->Branch("cl_sp_pe", "std::vector<double>", &_cl_sp_pe);
  fTree->Branch("cl_sp_ts0", "std::vector<double>", &_cl_sp_ts0);
  fTree->Branch("cl_sp_ets0", "std::vector<double>", &_cl_sp_ets0);
  fTree->Branch("cl_sp_ts1", "std::vector<double>", &_cl_sp_ts1);
  fTree->Branch("cl_sp_ets1", "std::vector<double>", &_cl_sp_ets1);
  fTree->Branch("cl_sp_dts0", "std::vector<double>", &_cl_sp_dts0);
  fTree->Branch("cl_sp_dts1", "std::vector<double>", &_cl_sp_dts1);
  fTree->Branch("cl_sp_complete", "std::vector<bool>", &_cl_sp_complete);

  if(!fDataMode && fTruthMatch)
    {
      fTree->Branch("td_tag_trackid", "std::vector<int>", &_td_tag_trackid);
      fTree->Branch("td_tag_pdg", "std::vector<int>", &_td_tag_pdg);
      fTree->Branch("td_tag_tagger", "std::vector<int16_t>", &_td_tag_tagger);
      fTree->Branch("td_tag_energy", "std::vector<double>", &_td_tag_energy);
      fTree->Branch("td_tag_time", "std::vector<double>", &_td_tag_time);
      fTree->Branch("td_tag_x", "std::vector<double>", &_td_tag_x);
      fTree->Branch("td_tag_y", "std::vector<double>", &_td_tag_y);
      fTree->Branch("td_tag_z", "std::vector<double>", &_td_tag_z);
      fTree->Branch("td_tag_reco_status", "std::vector<bool>", &_td_tag_reco_status);

      fTree->Branch("td_trackid", "std::vector<int>", &_td_trackid);
      fTree->Branch("td_pdg", "std::vector<int>", &_td_pdg);
      fTree->Branch("td_energy", "std::vector<double>", &_td_energy);
      fTree->Branch("td_time", "std::vector<double>", &_td_time);
      fTree->Branch("td_reconstructable", "std::vector<bool>", &_td_reconstructable);
      fTree->Branch("td_reco_status", "std::vector<bool>", &_td_reco_status);
      fTree->Branch("td_reco_triple", "std::vector<bool>", &_td_reco_triple);
    }

  fTree->Branch("tr_start_x", "std::vector<double>", &_tr_start_x);
  fTree->Branch("tr_start_y", "std::vector<double>", &_tr_start_y);
  fTree->Branch("tr_start_z", "std::vector<double>", &_tr_start_z);
  fTree->Branch("tr_end_x", "std::vector<double>", &_tr_end_x);
  fTree->Branch("tr_end_y", "std::vector<double>", &_tr_end_y);
  fTree->Branch("tr_end_z", "std::vector<double>", &_tr_end_z);
  fTree->Branch("tr_dir_x", "std::vector<double>", &_tr_dir_x);
  fTree->Branch("tr_dir_y", "std::vector<double>", &_tr_dir_y);
  fTree->Branch("tr_dir_z", "std::vector<double>", &_tr_dir_z);
  fTree->Branch("tr_ts0", "std::vector<double>", &_tr_ts0);
  fTree->Branch("tr_ets0", "std::vector<double>", &_tr_ets0);
  fTree->Branch("tr_ts1", "std::vector<double>", &_tr_ts1);
  fTree->Branch("tr_ets1", "std::vector<double>", &_tr_ets1);
  fTree->Branch("tr_pe", "std::vector<double>", &_tr_pe);
  fTree->Branch("tr_length", "std::vector<double>", &_tr_length);
  fTree->Branch("tr_tof", "std::vector<double>", &_tr_tof);
  fTree->Branch("tr_theta", "std::vector<double>", &_tr_theta);
  fTree->Branch("tr_phi", "std::vector<double>", &_tr_phi);
  fTree->Branch("tr_triple", "std::vector<bool>", &_tr_triple);
  fTree->Branch("tr_tagger1", "std::vector<int16_t>", &_tr_tagger1);
  fTree->Branch("tr_tagger2", "std::vector<int16_t>", &_tr_tagger2);
  fTree->Branch("tr_tagger3", "std::vector<int16_t>", &_tr_tagger3);
  fTree->Branch("tr_channel_set", "std::vector<std::vector<uint32_t>>", &_tr_channel_set);
  fTree->Branch("tr_adc_set", "std::vector<std::vector<uint16_t>>", &_tr_adc_set);
  if(!fDataMode && fTruthMatch)
    {
      fTree->Branch("tr_truth_trackid", "std::vector<int>", &_tr_truth_trackid);
      fTree->Branch("tr_truth_completeness", "std::vector<double>", &_tr_truth_completeness);
      fTree->Branch("tr_truth_purity", "std::vector<double>", &_tr_truth_purity);
      fTree->Branch("tr_truth_pdg", "std::vector<int>", &_tr_truth_pdg);
      fTree->Branch("tr_truth_energy", "std::vector<double>", &_tr_truth_energy);
      fTree->Branch("tr_truth_time", "std::vector<double>", &_tr_truth_time);
      fTree->Branch("tr_truth_start_x", "std::vector<double>", &_tr_truth_start_x);
      fTree->Branch("tr_truth_start_y", "std::vector<double>", &_tr_truth_start_y);
      fTree->Branch("tr_truth_start_z", "std::vector<double>", &_tr_truth_start_z);
      fTree->Branch("tr_truth_end_x", "std::vector<double>", &_tr_truth_end_x);
      fTree->Branch("tr_truth_end_y", "std::vector<double>", &_tr_truth_end_y);
      fTree->Branch("tr_truth_end_z", "std::vector<double>", &_tr_truth_end_z);
      fTree->Branch("tr_truth_dir_x", "std::vector<double>", &_tr_truth_dir_x);
      fTree->Branch("tr_truth_dir_y", "std::vector<double>", &_tr_truth_dir_y);
      fTree->Branch("tr_truth_dir_z", "std::vector<double>", &_tr_truth_dir_z);
      fTree->Branch("tr_truth_particle_energy", "std::vector<double>", &_tr_truth_particle_energy);
      fTree->Branch("tr_truth_length", "std::vector<double>", &_tr_truth_length);
      fTree->Branch("tr_truth_tof", "std::vector<double>", &_tr_truth_tof);
      fTree->Branch("tr_truth_theta", "std::vector<double>", &_tr_truth_theta);
      fTree->Branch("tr_truth_phi", "std::vector<double>", &_tr_truth_phi);
    }

  if(!fNoTPC)
    {
      fTree->Branch("tpc_start_x", "std::vector<double>", &_tpc_start_x);
      fTree->Branch("tpc_start_y", "std::vector<double>", &_tpc_start_y);
      fTree->Branch("tpc_start_z", "std::vector<double>", &_tpc_start_z);
      fTree->Branch("tpc_end_x", "std::vector<double>", &_tpc_end_x);
      fTree->Branch("tpc_end_y", "std::vector<double>", &_tpc_end_y);
      fTree->Branch("tpc_end_z", "std::vector<double>", &_tpc_end_z);
      fTree->Branch("tpc_start_dir_x", "std::vector<double>", &_tpc_start_dir_x);
      fTree->Branch("tpc_start_dir_y", "std::vector<double>", &_tpc_start_dir_y);
      fTree->Branch("tpc_start_dir_z", "std::vector<double>", &_tpc_start_dir_z);
      fTree->Branch("tpc_end_dir_x", "std::vector<double>", &_tpc_end_dir_x);
      fTree->Branch("tpc_end_dir_y", "std::vector<double>", &_tpc_end_dir_y);
      fTree->Branch("tpc_end_dir_z", "std::vector<double>", &_tpc_end_dir_z);
      fTree->Branch("tpc_length", "std::vector<double>", &_tpc_length);
      fTree->Branch("tpc_track_score", "std::vector<double>", &_tpc_track_score);
      fTree->Branch("tpc_sp_matched", "std::vector<bool>", &_tpc_sp_matched);
      fTree->Branch("tpc_sp_channel_set", "std::vector<std::vector<uint32_t>>", &_tpc_sp_channel_set);
      fTree->Branch("tpc_sp_xshift", "std::vector<double>", &_tpc_sp_xshift);
      fTree->Branch("tpc_sp_ts0", "std::vector<double>", &_tpc_sp_ts0);
      fTree->Branch("tpc_sp_ts1", "std::vector<double>", &_tpc_sp_ts1);
      fTree->Branch("tpc_sp_tagger", "std::vector<int>", &_tpc_sp_tagger);
      fTree->Branch("tpc_sp_nhits", "std::vector<int>", &_tpc_sp_nhits);
      fTree->Branch("tpc_sp_x", "std::vector<double>", &_tpc_sp_x);
      fTree->Branch("tpc_sp_y", "std::vector<double>", &_tpc_sp_y);
      fTree->Branch("tpc_sp_z", "std::vector<double>", &_tpc_sp_z);
      fTree->Branch("tpc_sp_score", "std::vector<double>", &_tpc_sp_score);
      fTree->Branch("tpc_tr_matched", "std::vector<bool>", &_tpc_tr_matched);
      fTree->Branch("tpc_tr_ts0", "std::vector<double>", &_tpc_tr_ts0);
      fTree->Branch("tpc_tr_ts1", "std::vector<double>", &_tpc_tr_ts1);
      fTree->Branch("tpc_tr_taggers", "std::vector<std::vector<int>>", &_tpc_tr_taggers);
      fTree->Branch("tpc_tr_start_x", "std::vector<double>", &_tpc_tr_start_x);
      fTree->Branch("tpc_tr_start_y", "std::vector<double>", &_tpc_tr_start_y);
      fTree->Branch("tpc_tr_start_z", "std::vector<double>", &_tpc_tr_start_z);
      fTree->Branch("tpc_tr_end_x", "std::vector<double>", &_tpc_tr_end_x);
      fTree->Branch("tpc_tr_end_y", "std::vector<double>", &_tpc_tr_end_y);
      fTree->Branch("tpc_tr_end_z", "std::vector<double>", &_tpc_tr_end_z);
      fTree->Branch("tpc_tr_score", "std::vector<double>", &_tpc_tr_score);
      if(!fDataMode && fTruthMatch)
        {
          fTree->Branch("tpc_truth_trackid", "std::vector<int>", &_tpc_truth_trackid);
          fTree->Branch("tpc_truth_pdg", "std::vector<int>", &_tpc_truth_pdg);
          fTree->Branch("tpc_truth_energy", "std::vector<double>", &_tpc_truth_energy);
          fTree->Branch("tpc_truth_time", "std::vector<double>", &_tpc_truth_time);
          fTree->Branch("tpc_sp_matchable", "std::vector<bool>", &_tpc_sp_matchable);
          fTree->Branch("tpc_sp_good_match", "std::vector<bool>", &_tpc_sp_good_match);
          fTree->Branch("tpc_tr_matchable", "std::vector<bool>", &_tpc_tr_matchable);
          fTree->Branch("tpc_tr_good_match", "std::vector<bool>", &_tpc_tr_good_match);
        }
    }

  if(fHasPTB)
    {
      fTree->Branch("ptb_hlt_trigger", "std::vector<uint64_t>", &_ptb_hlt_trigger);
      fTree->Branch("ptb_hlt_timestamp", "std::vector<uint64_t>", &_ptb_hlt_timestamp);
      fTree->Branch("ptb_llt_trigger", "std::vector<uint64_t>", &_ptb_llt_trigger);
      fTree->Branch("ptb_llt_timestamp", "std::vector<uint64_t>", &_ptb_llt_timestamp);
    }

  if(fHasTDC)
    {
      fTree->Branch("tdc_channel", "std::vector<uint32_t>", &_tdc_channel);
      fTree->Branch("tdc_timestamp", "std::vector<uint64_t>", &_tdc_timestamp);
      fTree->Branch("tdc_offset", "std::vector<uint64_t>", &_tdc_offset);
      fTree->Branch("tdc_name", "std::vector<std::string>", &_tdc_name);
    }

  if(fDebug)
    {
      for(auto const &[name, tagger] : fCRTGeoAlg.GetTaggers())
        {
          std::cout << "Tagger:  " << tagger.name << '\n'
                    << "X - Min: " << tagger.minX << " Max: " << tagger.maxX << '\n'
                    << "Y - Min: " << tagger.minY << " Max: " << tagger.maxY << '\n'
                    << "Z - Min: " << tagger.minZ << " Max: " << tagger.maxZ << '\n' << std::endl;
        }

      std::cout << std::endl;

      for(auto const &[name, module] : fCRTGeoAlg.GetModules())
        {
          std::cout << "Module:  " << module.name << " (" << module.taggerName << ")" << '\n';
          if(module.minos)
            std::cout << "MINOS module" << std::endl;
          std::cout << "X - Min: " << module.minX << " Max: " << module.maxX << " Diff: " << module.maxX - module.minX << '\n' 
                    << "Y - Min: " << module.minY << " Max: " << module.maxY << " Diff: " << module.maxY - module.minY << '\n' 
                    << "Z - Min: " << module.minZ << " Max: " << module.maxZ << " Diff: " << module.maxZ - module.minZ << '\n' 
                    << "Orientation: " << module.orientation << '\n' << std::endl;
        }

      std::cout << std::endl;

      for(auto const &[name, sipm] : fCRTGeoAlg.GetSiPMs())
        {
          std::cout << "SiPM:  " << sipm.channel << " (" << sipm.channel/32 << " - " << sipm.channel%32 << ")" << '\n'
                    << "x: " << sipm.x << " y: " << sipm.y << " z: " << sipm.z << std::endl;
        }
    }
}

void sbnd::crt::CRTAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(fDebug)
    std::cout << "CRTAnalysis -- This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  if(!fDataMode && fTruthMatch)
    {
      if(fDebug)
        std::cout << "CRTAnalysis -- Setting up back track maps" << std::endl;

      fCRTBackTrackerAlg.SetupMaps(e);
      fCRTBackTrackerAlg.RunSpacePointRecoStatusChecks(e);
      fCRTBackTrackerAlg.RunTrackRecoStatusChecks(e);
    }

  _crt_timing_reference_type    = -1;
  _crt_timing_reference_channel = -1;

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing Timing Reference Info" << std::endl;

  art::Handle<raw::TimingReferenceInfo> TimingReferenceHandle;
  e.getByLabel(fTimingReferenceModuleLabel, TimingReferenceHandle);
  if(TimingReferenceHandle.isValid())
    {
      _crt_timing_reference_type    = TimingReferenceHandle->timingType;
      _crt_timing_reference_channel = TimingReferenceHandle->timingChannel;
    }

  if(fHasPTB)
    {
      if(fDebug)
        std::cout << "CRTAnalysis -- Processing PTB" << std::endl;

      // Get PTBs
      art::Handle<std::vector<raw::ptb::sbndptb>> PTBHandle;
      e.getByLabel(fPTBModuleLabel, PTBHandle);
      if(!PTBHandle.isValid()){
        std::cout << "PTB product " << fPTBModuleLabel << " not found..." << std::endl;
        throw std::exception();
      }
      std::vector<art::Ptr<raw::ptb::sbndptb>> PTBVec;
      art::fill_ptr_vector(PTBVec, PTBHandle);

      // Fill PTB variables
      AnalysePTBs(PTBVec);
    }

  if(fHasTDC)
    {
      if(fDebug)
        std::cout << "CRTAnalysis -- Processing TDC" << std::endl;

      // Get TDCs
      art::Handle<std::vector<sbnd::timing::DAQTimestamp>> TDCHandle;
      e.getByLabel(fTDCModuleLabel, TDCHandle);
      if(!TDCHandle.isValid()){
        std::cout << "TDC product " << fTDCModuleLabel << " not found..." << std::endl;
        throw std::exception();
      }
      std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> TDCVec;
      art::fill_ptr_vector(TDCVec, TDCHandle);

      // Fill TDC variables
      AnalyseTDCs(TDCVec);
    }

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing FEBDatas" << std::endl;

  // Get FEBDatas
  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  if(!FEBDataHandle.isValid()){
    std::cout << "FEBData product " << fFEBDataModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Fill FEBData variables
  AnalyseFEBDatas(FEBDataVec);

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing StripHits" << std::endl;

  // Get CRTStripHits
  art::Handle<std::vector<CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  if(!CRTStripHitHandle.isValid()){
    std::cout << "CRTStripHit product " << fCRTStripHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  // Fill CRTStripHit variables
  AnalyseCRTStripHits(e, CRTStripHitVec);

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing Clusters" << std::endl;

  // Get CRTClusters
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTCluster>> CRTClusterVec;
  art::fill_ptr_vector(CRTClusterVec, CRTClusterHandle);

  // Get CRTCluster to CRTSpacePoint Assns
  art::FindManyP<CRTSpacePoint> clustersToSpacePoints(CRTClusterHandle, e, fCRTSpacePointModuleLabel);

  // Get CRTCluster to CRTStripHit Assns
  art::FindManyP<CRTStripHit> clustersToStripHits(CRTClusterHandle, e, fCRTClusterModuleLabel);

  // Fill CRTCluster variables
  AnalyseCRTClusters(e, CRTClusterVec, clustersToSpacePoints, clustersToStripHits);

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing Tracks" << std::endl;

  // Get CRTTracks
  art::Handle<std::vector<CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  if(!CRTTrackHandle.isValid()){
    std::cout << "CRTTrack product " << fCRTTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTTrack>> CRTTrackVec;
  art::fill_ptr_vector(CRTTrackVec, CRTTrackHandle);

  // Get CRTTrack to CRTSpacePoint Assns
  art::FindManyP<CRTSpacePoint> tracksToSpacePoints(CRTTrackHandle, e, fCRTTrackModuleLabel);

  // Get CRTSpacePoints
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get CRTSpacePoint to CRTCluster Assns
  art::FindOneP<CRTCluster> spacepointsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  // Fill CRTTrack variables
  AnalyseCRTTracks(e, CRTTrackVec, tracksToSpacePoints, spacepointsToClusters, clustersToStripHits);

  if(fNoTPC)
    {
      // Fill the Tree
      fTree->Fill();

      return;
    }

  // Get TPC tracks
  art::Handle<std::vector<recob::Track>> TPCTrackHandle;
  e.getByLabel(fTPCTrackModuleLabel, TPCTrackHandle);
  if(!TPCTrackHandle.isValid()){
    std::cout << "TPCTrack product " << fTPCTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get PFPs
  art::Handle<std::vector<recob::PFParticle>> PFPHandle;
  e.getByLabel(fPFPModuleLabel, PFPHandle);
  if(!PFPHandle.isValid()){
    std::cout << "PFP product " << fPFPModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  if(fDataMode)
    {
      if(fDebug)
        std::cout << "CRTAnalysis -- Processing TPC" << std::endl;

      // Fill TPC matching variables
      AnalyseTPCMatching(e, TPCTrackHandle, CRTSpacePointHandle, CRTClusterHandle, PFPHandle);

      // Fill the Tree
      fTree->Fill();

      return;
    }

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing MCParticles" << std::endl;

  // Get MCParticles
  art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
  e.getByLabel(fMCParticleModuleLabel, MCParticleHandle);
  if(!MCParticleHandle.isValid()){
    std::cout << "MCParticle product " << fMCParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCParticle>> MCParticleVec;
  art::fill_ptr_vector(MCParticleVec, MCParticleHandle);

  // Fill MCParticle variables
  AnalyseMCParticles(MCParticleVec);

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing Sim Deposits" << std::endl;

  // Get SimDeposits
  art::Handle<std::vector<sim::AuxDetSimChannel>> SimDepositHandle;
  e.getByLabel(fSimDepositModuleLabel, SimDepositHandle);
  if(!SimDepositHandle.isValid()){
    std::cout << "SimDeposit product " << fSimDepositModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sim::AuxDetSimChannel>> SimDepositVec;
  art::fill_ptr_vector(SimDepositVec, SimDepositHandle);

  // Fill SimDeposit variables
  AnalyseSimDeposits(SimDepositVec);

  if(fDebug)
    std::cout << "CRTAnalysis -- Processing CRT Back Tracking" << std::endl;

  if(fTruthMatch)
    {
      // Get Map of TrueDeposits per tagger from BackTracker
      std::map<CRTBackTrackerAlg::Category, bool> spacePointRecoStatusMap = fCRTBackTrackerAlg.GetSpacePointRecoStatusMap();

      // Fill TrueDeposit variables
      AnalyseTrueDepositsPerTagger(spacePointRecoStatusMap);

      // Get Map of TrueDeposits from BackTracker
      std::map<int, std::pair<bool, bool>> trackRecoStatusMap = fCRTBackTrackerAlg.GetTrackRecoStatusMap();

      // Fill TrueDeposit variables
      AnalyseTrueDeposits(trackRecoStatusMap);

      if(fDebug)
        std::cout << "CRTAnalysis -- Processing TPC" << std::endl;

      // Fill TPC matching variables
      AnalyseTPCMatching(e, TPCTrackHandle, CRTSpacePointHandle, CRTClusterHandle, PFPHandle, spacePointRecoStatusMap, trackRecoStatusMap);
    }
  else
    {
      if(fDebug)
        std::cout << "CRTAnalysis -- Processing TPC" << std::endl;

      // Fill TPC matching variables
      AnalyseTPCMatching(e, TPCTrackHandle, CRTSpacePointHandle, CRTClusterHandle, PFPHandle);
    }
  // Fill the Tree
  fTree->Fill();
}

void sbnd::crt::CRTAnalysis::AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec)
{
  unsigned nHLTs = 0;

  for(auto const& ptb : PTBVec)
    nHLTs += ptb->GetNHLTriggers();

  _ptb_hlt_trigger.resize(nHLTs);
  _ptb_hlt_timestamp.resize(nHLTs);

  unsigned hlt_i = 0;

  for(auto const& ptb : PTBVec)
    {
      for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i)
        {
          _ptb_hlt_trigger[hlt_i]   = ptb->GetHLTrigger(i).trigger_word;
          _ptb_hlt_timestamp[hlt_i] = ptb->GetHLTrigger(i).timestamp;

          ++hlt_i;
        }
    }

  unsigned nLLTs = 0;

  for(auto const& ptb : PTBVec)
    nLLTs += ptb->GetNLLTriggers();

  _ptb_llt_trigger.resize(nLLTs);
  _ptb_llt_timestamp.resize(nLLTs);

  unsigned llt_i = 0;

  for(auto const& ptb : PTBVec)
    {
      for(unsigned i = 0; i < ptb->GetNLLTriggers(); ++i)
        {
          _ptb_llt_trigger[llt_i]   = ptb->GetLLTrigger(i).trigger_word;
          _ptb_llt_timestamp[llt_i] = ptb->GetLLTrigger(i).timestamp;

          ++llt_i;
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec)
{
  const unsigned nTDCs = TDCVec.size();

  _tdc_channel.resize(nTDCs);
  _tdc_timestamp.resize(nTDCs);
  _tdc_offset.resize(nTDCs);
  _tdc_name.resize(nTDCs);

  unsigned tdc_i = 0;

  for(auto const& tdc : TDCVec)
    {
      _tdc_channel[tdc_i]   = tdc->Channel();
      _tdc_timestamp[tdc_i] = tdc->Timestamp();
      _tdc_offset[tdc_i]    = tdc->Offset();
      _tdc_name[tdc_i]      = tdc->Name();

      ++tdc_i;
    }
}

void sbnd::crt::CRTAnalysis::AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec)
{
  const unsigned nMCParticles = MCParticleVec.size();

  _mc_trackid.resize(nMCParticles);
  _mc_pdg.resize(nMCParticles);
  _mc_status.resize(nMCParticles);
  _mc_ndaughters.resize(nMCParticles);
  _mc_daughters.resize(nMCParticles);
  _mc_vx.resize(nMCParticles);
  _mc_vy.resize(nMCParticles);
  _mc_vz.resize(nMCParticles);
  _mc_vt.resize(nMCParticles);
  _mc_endx.resize(nMCParticles);
  _mc_endy.resize(nMCParticles);
  _mc_endz.resize(nMCParticles);
  _mc_endt.resize(nMCParticles);
  _mc_vpx.resize(nMCParticles);
  _mc_vpy.resize(nMCParticles);
  _mc_vpz.resize(nMCParticles);
  _mc_ve.resize(nMCParticles);
  _mc_endpx.resize(nMCParticles);
  _mc_endpy.resize(nMCParticles);
  _mc_endpz.resize(nMCParticles);
  _mc_ende.resize(nMCParticles);
  
  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      _mc_trackid[i]    = mcp->TrackId();
      _mc_pdg[i]        = mcp->PdgCode();
      _mc_status[i]     = mcp->StatusCode();
      _mc_ndaughters[i] = mcp->NumberDaughters();
      _mc_vx[i]         = mcp->Vx();
      _mc_vy[i]         = mcp->Vy();
      _mc_vz[i]         = mcp->Vz();
      _mc_vt[i]         = mcp->T();
      _mc_endx[i]       = mcp->EndX();
      _mc_endy[i]       = mcp->EndY();
      _mc_endz[i]       = mcp->EndZ();
      _mc_endt[i]       = mcp->EndT();
      _mc_vpx[i]        = mcp->Px();
      _mc_vpy[i]        = mcp->Py();
      _mc_vpz[i]        = mcp->Pz();
      _mc_ve[i]         = mcp->E();
      _mc_endpx[i]      = mcp->EndPx();
      _mc_endpy[i]      = mcp->EndPy();
      _mc_endpz[i]      = mcp->EndPz();
      _mc_ende[i]       = mcp->EndE();

      _mc_daughters[i].resize(_mc_ndaughters[i]);

      for(unsigned ii = 0; ii < _mc_ndaughters[i]; ++ii)
        _mc_daughters[i][ii] = mcp->Daughter(ii);
    }
}

void sbnd::crt::CRTAnalysis::AnalyseSimDeposits(std::vector<art::Ptr<sim::AuxDetSimChannel>> &SimDepositVec)
{
  const unsigned nAuxDetSimChannels = SimDepositVec.size();
  unsigned nIDEs = 0;

  for(unsigned i = 0; i < nAuxDetSimChannels; ++i)
    {
      const auto auxDetSimChannel = SimDepositVec[i];
      nIDEs += auxDetSimChannel->AuxDetIDEs().size();
    }

  _ide_trackid.resize(nIDEs);
  _ide_e.resize(nIDEs);
  _ide_entryx.resize(nIDEs);
  _ide_entryy.resize(nIDEs);
  _ide_entryz.resize(nIDEs);
  _ide_entryt.resize(nIDEs);
  _ide_exitx.resize(nIDEs);
  _ide_exity.resize(nIDEs);
  _ide_exitz.resize(nIDEs);
  _ide_exitt.resize(nIDEs);

  unsigned ide_counter = 0;

  for(unsigned i = 0; i < nAuxDetSimChannels; ++i)
    {
      const auto auxDetSimChannel = SimDepositVec[i];
      const auto ideVec           = auxDetSimChannel->AuxDetIDEs();

      for(unsigned ii = 0; ii < ideVec.size(); ++ii)
        {
          const auto ide = ideVec[ii];

          _ide_trackid[ide_counter] = ide.trackID;
          _ide_e[ide_counter]       = ide.energyDeposited;
          _ide_entryx[ide_counter]  = ide.entryX;
          _ide_entryy[ide_counter]  = ide.entryY;
          _ide_entryz[ide_counter]  = ide.entryZ;
          _ide_entryt[ide_counter]  = ide.entryT;
          _ide_exitx[ide_counter]   = ide.exitX;
          _ide_exity[ide_counter]   = ide.exitY;
          _ide_exitz[ide_counter]   = ide.exitZ;
          _ide_exitt[ide_counter]   = ide.exitT;
          
          ++ide_counter;
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseFEBDatas(std::vector<art::Ptr<FEBData>> &FEBDataVec)
{
  const unsigned nFEBData = FEBDataVec.size();

  _feb_mac5.resize(nFEBData);
  _feb_tagger.resize(nFEBData);
  _feb_flags.resize(nFEBData);
  _feb_ts0.resize(nFEBData);
  _feb_ts1.resize(nFEBData);
  _feb_unixs.resize(nFEBData);
  _feb_adc.resize(nFEBData, std::vector<uint16_t>(32));
  _feb_coinc.resize(nFEBData);

  for(unsigned i = 0; i < nFEBData; ++i)
    {
      const auto data = FEBDataVec[i];
      
      _feb_mac5[i]   = data->Mac5();
      _feb_tagger[i] = fCRTGeoAlg.AuxDetIndexToTaggerEnum(data->Mac5());
      _feb_flags[i]  = data->Flags();
      _feb_ts0[i]    = data->Ts0();
      _feb_ts1[i]    = data->Ts1();
      _feb_unixs[i]  = data->UnixS();
      _feb_coinc[i]  = data->Coinc();

      for(unsigned j = 0; j < 32; ++j)
        _feb_adc[i][j] = data->ADC(j);
    }
}

void sbnd::crt::CRTAnalysis::AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec)
{
  const unsigned nStripHits = CRTStripHitVec.size();
  
  _sh_channel.resize(nStripHits);
  _sh_tagger.resize(nStripHits);
  _sh_ts0.resize(nStripHits);
  _sh_ts1.resize(nStripHits);
  _sh_unixs.resize(nStripHits);
  _sh_pos.resize(nStripHits);
  _sh_err.resize(nStripHits);
  _sh_adc1.resize(nStripHits);
  _sh_adc2.resize(nStripHits);
  _sh_saturated1.resize(nStripHits);
  _sh_saturated2.resize(nStripHits);
  _sh_truth_trackid.resize(nStripHits);
  _sh_truth_completeness.resize(nStripHits);
  _sh_truth_purity.resize(nStripHits);
  _sh_truth_pos.resize(nStripHits);
  _sh_truth_energy.resize(nStripHits);
  _sh_truth_time.resize(nStripHits);

  for(unsigned i = 0; i < nStripHits; ++i)
    {
      const auto hit = CRTStripHitVec[i];

      _sh_channel[i]    = hit->Channel();
      _sh_tagger[i]     = fCRTGeoAlg.ChannelToTaggerEnum(hit->Channel());
      _sh_ts0[i]        = hit->Ts0();
      _sh_ts1[i]        = hit->Ts1();
      _sh_unixs[i]      = hit->UnixS();
      _sh_pos[i]        = hit->Pos();
      _sh_err[i]        = hit->Error();
      _sh_adc1[i]       = hit->ADC1();
      _sh_adc2[i]       = hit->ADC2();
      _sh_saturated1[i] = hit->Saturated1();
      _sh_saturated2[i] = hit->Saturated2();

      if(!fDataMode && fTruthMatch)
        {
          const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, hit);
          const std::vector<double> localpos = fCRTGeoAlg.StripWorldToLocalPos(hit->Channel(), truthMatch.deposit.x, truthMatch.deposit.y, truthMatch.deposit.z);
          const double width = fCRTGeoAlg.GetStrip(hit->Channel()).width;

          _sh_truth_trackid[i]      = truthMatch.trackid;
          _sh_truth_completeness[i] = truthMatch.completeness;
          _sh_truth_purity[i]       = truthMatch.purity;
          _sh_truth_pos[i]          = localpos[1] + width / 2.;
          _sh_truth_energy[i]       = truthMatch.deposit.energy;
          _sh_truth_time[i]         = truthMatch.deposit.time;
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec, const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints,
                                                const art::FindManyP<CRTStripHit> &clustersToStripHits)
{
  const unsigned nClusters = CRTClusterVec.size();

  _cl_ts0.resize(nClusters);
  _cl_ts1.resize(nClusters);
  _cl_unixs.resize(nClusters);
  _cl_nhits.resize(nClusters);
  _cl_tagger.resize(nClusters);
  _cl_composition.resize(nClusters);
  _cl_channel_set.resize(nClusters);
  _cl_adc_set.resize(nClusters);
  _cl_sh_ts0_set.resize(nClusters);
  _cl_sh_ts1_set.resize(nClusters);
  _cl_sh_feb_mac5_set.resize(nClusters);
  _cl_sh_time_walk_set.resize(nClusters);
  _cl_sh_prop_delay_set.resize(nClusters);
  _cl_truth_trackid.resize(nClusters);
  _cl_truth_completeness.resize(nClusters);
  _cl_truth_purity.resize(nClusters);
  _cl_truth_hit_completeness.resize(nClusters);
  _cl_truth_hit_purity.resize(nClusters);
  _cl_truth_pdg.resize(nClusters);
  _cl_truth_energy.resize(nClusters);
  _cl_truth_time.resize(nClusters);
  _cl_truth_x.resize(nClusters);
  _cl_truth_y.resize(nClusters);
  _cl_truth_z.resize(nClusters);
  _cl_truth_core_energy.resize(nClusters);
  _cl_truth_core_time.resize(nClusters);
  _cl_truth_core_x.resize(nClusters);
  _cl_truth_core_y.resize(nClusters);
  _cl_truth_core_z.resize(nClusters);
  _cl_has_sp.resize(nClusters);
  _cl_sp_x.resize(nClusters);
  _cl_sp_ex.resize(nClusters);
  _cl_sp_y.resize(nClusters);
  _cl_sp_ey.resize(nClusters);
  _cl_sp_z.resize(nClusters);
  _cl_sp_ez.resize(nClusters);
  _cl_sp_pe.resize(nClusters);
  _cl_sp_ts0.resize(nClusters);
  _cl_sp_ets0.resize(nClusters);
  _cl_sp_ts1.resize(nClusters);
  _cl_sp_ets1.resize(nClusters);
  _cl_sp_dts0.resize(nClusters);
  _cl_sp_dts1.resize(nClusters);
  _cl_sp_complete.resize(nClusters);

  for(unsigned i = 0; i < nClusters; ++i)
    {
      const auto cluster = CRTClusterVec[i];

      _cl_ts0[i]         = cluster->Ts0();
      _cl_ts1[i]         = cluster->Ts1();
      _cl_unixs[i]       = cluster->UnixS();
      _cl_nhits[i]       = cluster->NHits();
      _cl_tagger[i]      = cluster->Tagger();
      _cl_composition[i] = cluster->Composition();

      const auto spacepoints = clustersToSpacePoints.at(cluster.key());

      const auto striphits = clustersToStripHits.at(cluster.key());
      _cl_channel_set[i].resize(_cl_nhits[i]);
      _cl_adc_set[i].resize(2 * _cl_nhits[i]);
      _cl_sh_ts0_set[i].resize(_cl_nhits[i]);
      _cl_sh_ts1_set[i].resize(_cl_nhits[i]);
      _cl_sh_feb_mac5_set[i].resize(_cl_nhits[i]);
      _cl_sh_time_walk_set[i].resize(_cl_nhits[i]);
      _cl_sh_prop_delay_set[i].resize(_cl_nhits[i]);

      std::vector<std::pair<int, double>> ts0_set, ts1_set;

      for(unsigned ii = 0; ii < _cl_nhits[i]; ++ii)
        {
          const auto striphit = striphits[ii];
          _cl_channel_set[i][ii] = striphit->Channel();
          _cl_adc_set[i][2*ii]   = striphit->ADC1();
          _cl_adc_set[i][2*ii+1] = striphit->ADC2();
          _cl_sh_ts0_set[i][ii]  = striphit->Ts0();
          _cl_sh_ts1_set[i][ii]  = striphit->Ts1();

          if(fDataMode)
            {
              art::ServiceHandle<SBND::CRTChannelMapService> ChannelMapService;
              SBND::CRTChannelMapService::ModuleInfo_t module_info = ChannelMapService->GetModuleInfoFromOfflineID( striphit->Channel() / 32 );
              _cl_sh_feb_mac5_set[i][ii] = ( module_info.valid ) ? module_info.feb_mac5 : 0;
            }
          else
            _cl_sh_feb_mac5_set[i][ii] = striphit->Channel() / 32;

          /*
           * The below segment reimplements the CorrectTime() method
           * from CRTReco/CRTClusterCharacterisationAlg.cc .
           * Because the Ts0(), Ts1() getters invoked in _cl_sp_ts*, _cl_sh_ts*_set are raw T0/1
           * counters, the time walk and propagation delay are saved as explicit branches here.
           */
          if(spacepoints.size() == 1)
            {
              double pe0 = fCRTGeoAlg.GetSiPM( striphit->Channel() ).gain * striphit->ADC1();
              double pe1 = fCRTGeoAlg.GetSiPM( striphit->Channel() + 1 ).gain * striphit->ADC2();
              double pe  = pe0 + pe1;

              double dist = fCRTGeoAlg.DistanceDownStrip( spacepoints[0]->Pos(), striphit->Channel() );

              double corr = std::pow( dist - fPEAttenuation, 2.0 ) / std::pow( fPEAttenuation, 2.0 );
              double tw_pe = pe * corr;

              _cl_sh_time_walk_set[i][ii]  = fTimeWalkNorm * std::exp( -fTimeWalkScale * tw_pe );
              _cl_sh_prop_delay_set[i][ii] = fPropDelay * dist;

              ts0_set.push_back({_cl_sh_feb_mac5_set[i][ii], _cl_sh_ts0_set[i][ii] - _cl_sh_time_walk_set[i][ii] - _cl_sh_prop_delay_set[i][ii]});
              ts1_set.push_back({_cl_sh_feb_mac5_set[i][ii], _cl_sh_ts1_set[i][ii] - _cl_sh_time_walk_set[i][ii] - _cl_sh_prop_delay_set[i][ii]});
            }
          else
            {
              _cl_sh_time_walk_set[i][ii]  = -999999.;
              _cl_sh_prop_delay_set[i][ii] = -999999.;
            }
        }

      if(!fDataMode && fTruthMatch)
        {
          const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, cluster);
          _cl_truth_trackid[i]          = truthMatch.trackid;
          _cl_truth_completeness[i]     = truthMatch.completeness;
          _cl_truth_purity[i]           = truthMatch.purity;
          _cl_truth_hit_completeness[i] = truthMatch.hitcompleteness;
          _cl_truth_hit_purity[i]       = truthMatch.hitpurity;
          _cl_truth_pdg[i]              = truthMatch.deposit.pdg;
          _cl_truth_energy[i]           = truthMatch.deposit.energy;
          _cl_truth_time[i]             = truthMatch.deposit.time;
          _cl_truth_x[i]                = truthMatch.deposit.x;
          _cl_truth_y[i]                = truthMatch.deposit.y;
          _cl_truth_z[i]                = truthMatch.deposit.z;
          _cl_truth_core_energy[i]      = truthMatch.deposit.coreEnergy;
          _cl_truth_core_time[i]        = truthMatch.deposit.coreTime;
          _cl_truth_core_x[i]           = truthMatch.deposit.coreX;
          _cl_truth_core_y[i]           = truthMatch.deposit.coreY;
          _cl_truth_core_z[i]           = truthMatch.deposit.coreZ;
        }

      if(spacepoints.size() == 1)
        { 
          const auto spacepoint = spacepoints[0];
          
          _cl_has_sp[i]      = true;
          _cl_sp_x[i]        = spacepoint->X();
          _cl_sp_ex[i]       = spacepoint->XErr();
          _cl_sp_y[i]        = spacepoint->Y();
          _cl_sp_ey[i]       = spacepoint->YErr();
          _cl_sp_z[i]        = spacepoint->Z();
          _cl_sp_ez[i]       = spacepoint->ZErr();
          _cl_sp_pe[i]       = spacepoint->PE();
          _cl_sp_ts0[i]      = spacepoint->Ts0();
          _cl_sp_ets0[i]     = spacepoint->Ts0Err();
          _cl_sp_ts1[i]      = spacepoint->Ts1();
          _cl_sp_ets1[i]     = spacepoint->Ts1Err();
          _cl_sp_dts0[i]     = SPTimeDelta(ts0_set);
          _cl_sp_dts1[i]     = SPTimeDelta(ts1_set);
          _cl_sp_complete[i] = spacepoint->Complete();
        }
      else
        {
          _cl_has_sp[i]      = false;
          _cl_sp_x[i]        = -999999.;
          _cl_sp_ex[i]       = -999999.;
          _cl_sp_y[i]        = -999999.;
          _cl_sp_ey[i]       = -999999.;
          _cl_sp_z[i]        = -999999.;
          _cl_sp_ez[i]       = -999999.;
          _cl_sp_pe[i]       = -999999.;
          _cl_sp_ts0[i]      = -999999.;
          _cl_sp_ets0[i]     = -999999.;
          _cl_sp_ts1[i]      = -999999.;
          _cl_sp_ets1[i]     = -999999.;
          _cl_sp_dts0[i]     = -999999.;
          _cl_sp_dts1[i]     = -999999.;
          _cl_sp_complete[i] = false;
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseTrueDepositsPerTagger(const std::map<CRTBackTrackerAlg::Category, bool> &recoStatusMap)
{
  const unsigned nTrueDeposits = recoStatusMap.size();

  _td_tag_trackid.resize(nTrueDeposits);
  _td_tag_pdg.resize(nTrueDeposits);
  _td_tag_tagger.resize(nTrueDeposits);
  _td_tag_energy.resize(nTrueDeposits);
  _td_tag_time.resize(nTrueDeposits);
  _td_tag_x.resize(nTrueDeposits);
  _td_tag_y.resize(nTrueDeposits);
  _td_tag_z.resize(nTrueDeposits);
  _td_tag_reco_status.resize(nTrueDeposits);

  unsigned entry = 0;
  for(auto const& [category, status] : recoStatusMap)
    {
      const CRTBackTrackerAlg::TrueDeposit deposit = fCRTBackTrackerAlg.GetTrueDeposit(category);

      _td_tag_trackid[entry]     = deposit.trackid;
      _td_tag_pdg[entry]         = deposit.pdg;
      _td_tag_tagger[entry]      = deposit.tagger;
      _td_tag_energy[entry]      = deposit.energy;
      _td_tag_time[entry]        = deposit.time;
      _td_tag_x[entry]           = deposit.x;
      _td_tag_y[entry]           = deposit.y;
      _td_tag_z[entry]           = deposit.z;
      _td_tag_reco_status[entry] = status;

      ++entry;
    }
}

void sbnd::crt::CRTAnalysis::AnalyseTrueDeposits(const std::map<int, std::pair<bool, bool>> &recoStatusMap)
{
  const unsigned nTrueDeposits = recoStatusMap.size();

  _td_trackid.resize(nTrueDeposits);
  _td_pdg.resize(nTrueDeposits);
  _td_energy.resize(nTrueDeposits);
  _td_time.resize(nTrueDeposits);
  _td_reconstructable.resize(nTrueDeposits);
  _td_reco_status.resize(nTrueDeposits);
  _td_reco_triple.resize(nTrueDeposits);

  unsigned entry = 0;
  for(auto const& [trackid, status] : recoStatusMap)
    {
      const CRTBackTrackerAlg::TrueDeposit deposit = fCRTBackTrackerAlg.GetTrueDeposit(trackid);

      _td_trackid[entry]         = deposit.trackid;
      _td_pdg[entry]             = deposit.pdg;
      _td_energy[entry]          = deposit.energy;
      _td_time[entry]            = deposit.time;
      _td_reconstructable[entry] = deposit.reconstructable;
      _td_reco_status[entry]     = status.first;
      _td_reco_triple[entry]     = status.second;

      ++entry;
    }
}

void sbnd::crt::CRTAnalysis::AnalyseCRTTracks(const art::Event &e, const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec, const art::FindManyP<CRTSpacePoint> &tracksToSpacePoints,
                                              const art::FindOneP<CRTCluster> &spacePointsToClusters, const art::FindManyP<CRTStripHit> &clustersToStripHits)
{
  const unsigned nTracks = CRTTrackVec.size();

  _tr_start_x.resize(nTracks);
  _tr_start_y.resize(nTracks);
  _tr_start_z.resize(nTracks);
  _tr_end_x.resize(nTracks);
  _tr_end_y.resize(nTracks);
  _tr_end_z.resize(nTracks);
  _tr_dir_x.resize(nTracks);
  _tr_dir_y.resize(nTracks);
  _tr_dir_z.resize(nTracks);
  _tr_ts0.resize(nTracks);
  _tr_ets0.resize(nTracks);
  _tr_ts1.resize(nTracks);
  _tr_ets1.resize(nTracks);
  _tr_pe.resize(nTracks);
  _tr_length.resize(nTracks);
  _tr_tof.resize(nTracks);
  _tr_theta.resize(nTracks);
  _tr_phi.resize(nTracks);
  _tr_triple.resize(nTracks);
  _tr_tagger1.resize(nTracks);
  _tr_tagger2.resize(nTracks);
  _tr_tagger3.resize(nTracks);
  _tr_channel_set.resize(nTracks);
  _tr_adc_set.resize(nTracks);
  _tr_truth_trackid.resize(nTracks);
  _tr_truth_completeness.resize(nTracks);
  _tr_truth_purity.resize(nTracks);
  _tr_truth_pdg.resize(nTracks);
  _tr_truth_energy.resize(nTracks);
  _tr_truth_time.resize(nTracks);
  _tr_truth_start_x.resize(nTracks);
  _tr_truth_start_y.resize(nTracks);
  _tr_truth_start_z.resize(nTracks);
  _tr_truth_end_x.resize(nTracks);
  _tr_truth_end_y.resize(nTracks);
  _tr_truth_end_z.resize(nTracks);
  _tr_truth_dir_x.resize(nTracks);
  _tr_truth_dir_y.resize(nTracks);
  _tr_truth_dir_z.resize(nTracks);
  _tr_truth_particle_energy.resize(nTracks);
  _tr_truth_length.resize(nTracks);
  _tr_truth_tof.resize(nTracks);
  _tr_truth_theta.resize(nTracks);
  _tr_truth_phi.resize(nTracks);

  for(unsigned i = 0; i < nTracks; ++i)
    {
      const auto track = CRTTrackVec[i];

      const geo::Point_t start = track->Start();
      _tr_start_x[i] = start.X();
      _tr_start_y[i] = start.Y();
      _tr_start_z[i] = start.Z();

      const geo::Point_t end = track->End();
      _tr_end_x[i] = end.X();
      _tr_end_y[i] = end.Y();
      _tr_end_z[i] = end.Z();

      const geo::Vector_t dir = track->Direction();
      _tr_dir_x[i] = dir.X();
      _tr_dir_y[i] = dir.Y();
      _tr_dir_z[i] = dir.Z();

      _tr_ts0[i]     = track->Ts0();
      _tr_ets0[i]    = track->Ts0Err();
      _tr_ts1[i]     = track->Ts1();
      _tr_ets1[i]    = track->Ts1Err();
      _tr_pe[i]      = track->PE();
      _tr_length[i]  = track->Length();
      _tr_tof[i]     = track->ToF();
      _tr_theta[i]   = TMath::RadToDeg() * track->Theta();
      _tr_phi[i]     = TMath::RadToDeg() * track->Phi();
      _tr_triple[i]  = track->Triple();

      unsigned tag_i = 0;

      for(auto const &tagger : track->Taggers())
        {
          if(tag_i == 0)
            _tr_tagger1[i] = tagger;
          else if(tag_i == 1)
            _tr_tagger2[i] = tagger;
          else if(tag_i == 2)
            _tr_tagger3[i] = tagger;

          ++tag_i;
        }

      _tr_channel_set[i].clear();
      _tr_adc_set[i].clear();

      const auto spacepoints = tracksToSpacePoints.at(track.key());

      for(auto const& spacepoint : spacepoints)
        {
          const auto cluster = spacePointsToClusters.at(spacepoint.key());

          const auto striphits = clustersToStripHits.at(cluster.key());

          for(auto const& striphit : striphits)
            {
              _tr_channel_set[i].push_back(striphit->Channel());
              _tr_adc_set[i].push_back(striphit->ADC1());
              _tr_adc_set[i].push_back(striphit->ADC2());
            }
        }

      if(!fDataMode && fTruthMatch)
        {
          const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, track);
          _tr_truth_trackid[i]          = truthMatch.trackid;
          _tr_truth_completeness[i]     = truthMatch.completeness;
          _tr_truth_purity[i]           = truthMatch.purity;
          _tr_truth_pdg[i]              = truthMatch.deposit.pdg;
          _tr_truth_energy[i]           = truthMatch.deposit.energy;
          _tr_truth_time[i]             = truthMatch.deposit.time;

          _tr_truth_particle_energy[i] = truthMatch.trackinfo.energy;

          const geo::Point_t true_start(truthMatch.trackinfo.deposit1.x, truthMatch.trackinfo.deposit1.y, truthMatch.trackinfo.deposit1.z);
          const geo::Point_t true_end(truthMatch.trackinfo.deposit2.x, truthMatch.trackinfo.deposit2.y, truthMatch.trackinfo.deposit2.z);
          const geo::Vector_t true_dir = (true_end - true_start).Unit();

          _tr_truth_start_x[i] = true_start.X();
          _tr_truth_start_y[i] = true_start.Y();
          _tr_truth_start_z[i] = true_start.Z();

          _tr_truth_end_x[i] = true_end.X();
          _tr_truth_end_y[i] = true_end.Y();
          _tr_truth_end_z[i] = true_end.Z();

          _tr_truth_dir_x[i] = true_dir.X();
          _tr_truth_dir_y[i] = true_dir.Y();
          _tr_truth_dir_z[i] = true_dir.Z();

          _tr_truth_length[i] = (true_end - true_start).R();
          _tr_truth_tof[i]    = truthMatch.trackinfo.deposit2.time - truthMatch.trackinfo.deposit1.time;
          _tr_truth_theta[i]  = TMath::RadToDeg() * true_dir.Theta();
          _tr_truth_phi[i]    = TMath::RadToDeg() * true_dir.Phi();
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseTPCMatching(const art::Event &e, const art::Handle<std::vector<recob::Track>> &TPCTrackHandle,
                                                const art::Handle<std::vector<CRTSpacePoint>> &CRTSpacePointHandle, const art::Handle<std::vector<CRTCluster>> &CRTClusterHandle,
                                                const art::Handle<std::vector<recob::PFParticle>> &PFPHandle,
                                                const std::map<CRTBackTrackerAlg::Category, bool> &spacePointRecoStatusMap,
                                                const std::map<int, std::pair<bool, bool>> &trackRecoStatusMap)
{
  std::vector<art::Ptr<recob::Track>> TPCTrackVec;
  art::fill_ptr_vector(TPCTrackVec, TPCTrackHandle);

  const unsigned nTracks = TPCTrackVec.size();

  _tpc_start_x.resize(nTracks);
  _tpc_start_y.resize(nTracks);
  _tpc_start_z.resize(nTracks);
  _tpc_end_x.resize(nTracks);
  _tpc_end_y.resize(nTracks);
  _tpc_end_z.resize(nTracks);
  _tpc_start_dir_x.resize(nTracks);
  _tpc_start_dir_y.resize(nTracks);
  _tpc_start_dir_z.resize(nTracks);
  _tpc_end_dir_x.resize(nTracks);
  _tpc_end_dir_y.resize(nTracks);
  _tpc_end_dir_z.resize(nTracks);
  _tpc_length.resize(nTracks);
  _tpc_track_score.resize(nTracks);
  _tpc_truth_trackid.resize(nTracks);
  _tpc_truth_pdg.resize(nTracks);
  _tpc_truth_energy.resize(nTracks);
  _tpc_truth_time.resize(nTracks);
  _tpc_sp_matchable.resize(nTracks);
  _tpc_sp_matched.resize(nTracks);
  _tpc_sp_channel_set.resize(nTracks);
  _tpc_sp_good_match.resize(nTracks);
  _tpc_sp_xshift.resize(nTracks);
  _tpc_sp_ts0.resize(nTracks);
  _tpc_sp_ts1.resize(nTracks);
  _tpc_sp_tagger.resize(nTracks);
  _tpc_sp_nhits.resize(nTracks);
  _tpc_sp_x.resize(nTracks);
  _tpc_sp_y.resize(nTracks);
  _tpc_sp_z.resize(nTracks);
  _tpc_sp_score.resize(nTracks);
  _tpc_tr_matchable.resize(nTracks);
  _tpc_tr_matched.resize(nTracks);
  _tpc_tr_good_match.resize(nTracks);
  _tpc_tr_ts0.resize(nTracks);
  _tpc_tr_ts1.resize(nTracks);
  _tpc_tr_taggers.resize(nTracks);
  _tpc_tr_start_x.resize(nTracks);
  _tpc_tr_start_y.resize(nTracks);
  _tpc_tr_start_z.resize(nTracks);
  _tpc_tr_end_x.resize(nTracks);
  _tpc_tr_end_y.resize(nTracks);
  _tpc_tr_end_z.resize(nTracks);
  _tpc_tr_score.resize(nTracks);

  art::FindOneP<recob::PFParticle>                 tracksToPFPs(TPCTrackHandle, e, fTPCTrackModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpsToMetadata(PFPHandle, e, fPFPModuleLabel);
  art::FindManyP<recob::Hit>                       tracksToHits(TPCTrackHandle, e, fTPCTrackModuleLabel);
  art::FindOneP<CRTSpacePoint, anab::T0>           tracksToSPMatches(TPCTrackHandle, e, fCRTSpacePointMatchingModuleLabel);
  art::FindOneP<CRTCluster>                        spsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);
  art::FindManyP<CRTStripHit>                      clustersToStripHits(CRTClusterHandle, e, fCRTClusterModuleLabel);
  art::FindOneP<CRTTrack, anab::T0>                tracksToTrackMatches(TPCTrackHandle, e, fCRTTrackMatchingModuleLabel);

  const detinfo::DetectorClocksData clockData   = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const detinfo::DetectorPropertiesData detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  geo::GeometryCore const* geometryService      = lar::providerFrom<geo::Geometry>();

  int nActualTracks = 0;

  for(unsigned i = 0; i < nTracks; ++i)
    {
      const auto track = TPCTrackVec[i];
      const auto pfp   = tracksToPFPs.at(track.key());
      const auto meta  = pfpsToMetadata.at(pfp.key());
      const auto props = meta->GetPropertiesMap();

      if(pfp->PdgCode() != 13)
        continue;

      const geo::Point_t start = track->Start();
      _tpc_start_x[nActualTracks] = start.X();
      _tpc_start_y[nActualTracks] = start.Y();
      _tpc_start_z[nActualTracks] = start.Z();

      const geo::Point_t end = track->End();
      _tpc_end_x[nActualTracks] = end.X();
      _tpc_end_y[nActualTracks] = end.Y();
      _tpc_end_z[nActualTracks] = end.Z();

      const std::pair<geo::Vector_t, geo::Vector_t> dirs = CRTCommonUtils::AverageTrackDirections(track, 0.2);

      const geo::Vector_t start_dir = dirs.first;
      _tpc_start_dir_x[nActualTracks] = start_dir.X();
      _tpc_start_dir_y[nActualTracks] = start_dir.Y();
      _tpc_start_dir_z[nActualTracks] = start_dir.Z();

      const geo::Vector_t end_dir = dirs.first;
      _tpc_end_dir_x[nActualTracks] = end_dir.X();
      _tpc_end_dir_y[nActualTracks] = end_dir.Y();
      _tpc_end_dir_z[nActualTracks] = end_dir.Z();

      _tpc_length[nActualTracks] = track->Length();

      const auto trackscore = props.find("TrackScore");
      if(trackscore != props.end())
        _tpc_track_score[nActualTracks] = trackscore->second;
      else
        _tpc_track_score[nActualTracks] = -std::numeric_limits<double>::max();

      const art::Ptr<CRTSpacePoint> spacepoint = tracksToSPMatches.at(track.key());
      const art::Ptr<CRTTrack> crttrack = tracksToTrackMatches.at(track.key());

      const std::vector<art::Ptr<recob::Hit>> trackHits = tracksToHits.at(track.key());

      if(spacepoint.isNonnull())
        {
          const anab::T0 spMatch                             = tracksToSPMatches.data(track.key()).ref();
          const art::Ptr<CRTCluster> cluster                 = spsToClusters.at(spacepoint.key());
          const std::vector<art::Ptr<CRTStripHit>> striphits = clustersToStripHits.at(cluster.key());

          const int driftDirection     = TPCGeoUtil::DriftDirectionFromHits(geometryService, trackHits);
          const double crtShiftingTime = fDataMode ? spacepoint->Ts0() * 1e-3 : spacepoint->Ts1() * 1e-3;

          _tpc_sp_matched[nActualTracks] = true;
          _tpc_sp_xshift[nActualTracks]  = driftDirection * crtShiftingTime * detProp.DriftVelocity();
          _tpc_sp_ts0[nActualTracks]     = spacepoint->Ts0();
          _tpc_sp_ts1[nActualTracks]     = spacepoint->Ts1();
          _tpc_sp_tagger[nActualTracks]  = cluster->Tagger();
          _tpc_sp_nhits[nActualTracks]   = cluster->NHits();
          _tpc_sp_x[nActualTracks]       = spacepoint->X();
          _tpc_sp_y[nActualTracks]       = spacepoint->Y();
          _tpc_sp_z[nActualTracks]       = spacepoint->Z();
          _tpc_sp_score[nActualTracks]   = spMatch.TriggerConfidence();

          _tpc_sp_channel_set[nActualTracks] = std::vector<uint32_t>();

          for(auto const &striphit : striphits)
            _tpc_sp_channel_set[nActualTracks].push_back(striphit->Channel());
        }
      else
        {
          _tpc_sp_matched[nActualTracks] = false;
          _tpc_sp_ts0[nActualTracks]     = -std::numeric_limits<double>::max();
          _tpc_sp_ts1[nActualTracks]     = -std::numeric_limits<double>::max();
          _tpc_sp_tagger[nActualTracks]  = -std::numeric_limits<int>::max();
          _tpc_sp_x[nActualTracks]       = -std::numeric_limits<double>::max();
          _tpc_sp_y[nActualTracks]       = -std::numeric_limits<double>::max();
          _tpc_sp_z[nActualTracks]       = -std::numeric_limits<double>::max();
          _tpc_sp_score[nActualTracks]   = -std::numeric_limits<double>::max();
        }

      if(crttrack.isNonnull())
        {
          const anab::T0 trackMatch = tracksToTrackMatches.data(track.key()).ref();

          _tpc_tr_matched[nActualTracks] = true;
          _tpc_tr_ts0[nActualTracks]     = crttrack->Ts0();
          _tpc_tr_ts1[nActualTracks]     = crttrack->Ts1();
          _tpc_tr_score[nActualTracks]   = trackMatch.TriggerConfidence();

          const std::set<CRTTagger> taggers = crttrack->Taggers();
          const geo::Point_t start          = crttrack->Start();
          const geo::Point_t end            = crttrack->End();

          _tpc_tr_taggers[nActualTracks] = std::vector<int>();

          for(auto const& tagger : taggers)
            _tpc_tr_taggers[nActualTracks].push_back(tagger);

          _tpc_tr_start_x[nActualTracks] = start.X();
          _tpc_tr_start_y[nActualTracks] = start.Y();
          _tpc_tr_start_z[nActualTracks] = start.Z();
          _tpc_tr_end_x[nActualTracks]   = end.X();
          _tpc_tr_end_y[nActualTracks]   = end.Y();
          _tpc_tr_end_z[nActualTracks]   = end.Z();
        }
      else
        {
          _tpc_tr_matched[nActualTracks] = false;
          _tpc_tr_ts0[nActualTracks]     = -std::numeric_limits<double>::max();
          _tpc_tr_ts1[nActualTracks]     = -std::numeric_limits<double>::max();
          _tpc_tr_score[nActualTracks]   = -std::numeric_limits<double>::max();
          _tpc_tr_taggers[nActualTracks] = std::vector<int>();
          _tpc_tr_start_x[nActualTracks] = -std::numeric_limits<double>::max();
          _tpc_tr_start_y[nActualTracks] = -std::numeric_limits<double>::max();
          _tpc_tr_start_z[nActualTracks] = -std::numeric_limits<double>::max();
          _tpc_tr_end_x[nActualTracks]   = -std::numeric_limits<double>::max();
          _tpc_tr_end_y[nActualTracks]   = -std::numeric_limits<double>::max();
          _tpc_tr_end_z[nActualTracks]   = -std::numeric_limits<double>::max();
        }

      if(!fDataMode && fTruthMatch)
        {
          const int trackid = fCRTBackTrackerAlg.RollUpID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true));

          _tpc_truth_trackid[nActualTracks] = trackid;

          int pdg;
          double energy, time;
          fCRTBackTrackerAlg.TrueParticlePDGEnergyTime(trackid, pdg, energy, time);

          _tpc_truth_pdg[nActualTracks]    = pdg;
          _tpc_truth_energy[nActualTracks] = energy;
          _tpc_truth_time[nActualTracks]   = time;

          bool sp_matchable = false, tr_matchable = false;

          for(auto const& [ category, status ] : spacePointRecoStatusMap)
            sp_matchable |= (category.trackid == trackid && status);

          tr_matchable |= trackRecoStatusMap.count(trackid) != 0;
          if(tr_matchable) tr_matchable = trackRecoStatusMap.at(trackid).first;

          _tpc_sp_matchable[nActualTracks] = sp_matchable;
          _tpc_tr_matchable[nActualTracks] = tr_matchable;

          if(spacepoint.isNonnull())
            {
              const art::Ptr<CRTCluster> cluster                    = spsToClusters.at(spacepoint.key());
              const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, cluster);

              _tpc_sp_good_match[nActualTracks] = truthMatch.trackid == trackid;
            }
          else
            {
              _tpc_sp_good_match[nActualTracks] = false;
            }

          if(crttrack.isNonnull())
            {
              const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, crttrack);

              _tpc_tr_good_match[nActualTracks] = truthMatch.trackid == trackid;
            }
          else
            {
              _tpc_tr_good_match[nActualTracks] = false;
            }
        }

      ++nActualTracks;
    }

  _tpc_start_x.resize(nActualTracks);
  _tpc_start_y.resize(nActualTracks);
  _tpc_start_z.resize(nActualTracks);
  _tpc_end_x.resize(nActualTracks);
  _tpc_end_y.resize(nActualTracks);
  _tpc_end_z.resize(nActualTracks);
  _tpc_start_dir_x.resize(nActualTracks);
  _tpc_start_dir_y.resize(nActualTracks);
  _tpc_start_dir_z.resize(nActualTracks);
  _tpc_end_dir_x.resize(nActualTracks);
  _tpc_end_dir_y.resize(nActualTracks);
  _tpc_end_dir_z.resize(nActualTracks);
  _tpc_length.resize(nActualTracks);
  _tpc_track_score.resize(nActualTracks);
  _tpc_truth_trackid.resize(nActualTracks);
  _tpc_truth_pdg.resize(nActualTracks);
  _tpc_truth_energy.resize(nActualTracks);
  _tpc_truth_time.resize(nActualTracks);
  _tpc_sp_matchable.resize(nActualTracks);
  _tpc_sp_matched.resize(nActualTracks);
  _tpc_sp_tagger.resize(nActualTracks);
  _tpc_sp_nhits.resize(nActualTracks);
  _tpc_sp_x.resize(nActualTracks);
  _tpc_sp_y.resize(nActualTracks);
  _tpc_sp_z.resize(nActualTracks);
  _tpc_sp_good_match.resize(nActualTracks);
  _tpc_sp_ts0.resize(nActualTracks);
  _tpc_sp_ts1.resize(nActualTracks);
  _tpc_sp_score.resize(nActualTracks);
  _tpc_tr_matchable.resize(nActualTracks);
  _tpc_tr_matched.resize(nActualTracks);
  _tpc_tr_taggers.resize(nActualTracks);
  _tpc_tr_start_x.resize(nActualTracks);
  _tpc_tr_start_y.resize(nActualTracks);
  _tpc_tr_start_z.resize(nActualTracks);
  _tpc_tr_end_x.resize(nActualTracks);
  _tpc_tr_end_y.resize(nActualTracks);
  _tpc_tr_end_z.resize(nActualTracks);
  _tpc_tr_good_match.resize(nActualTracks);
  _tpc_tr_ts0.resize(nActualTracks);
  _tpc_tr_ts1.resize(nActualTracks);
  _tpc_tr_score.resize(nActualTracks);
}

double sbnd::crt::CRTAnalysis::SPTimeDelta(std::vector<std::pair<int, double>> &time_set)
{
  std::sort(time_set.begin(), time_set.end(), [](auto &a, auto &b)
  { return a.first < b.first; });

  double sum = 0;

  for(unsigned int i = 0; i < time_set.size(); ++i)
    {
      for(unsigned int ii = i + 1; ii < time_set.size(); ++ii)
        sum += (time_set[i].second - time_set[ii].second);
    }

  sum /= (time_set.size() - 1);

  return sum;
}

DEFINE_ART_MODULE(sbnd::crt::CRTAnalysis)
