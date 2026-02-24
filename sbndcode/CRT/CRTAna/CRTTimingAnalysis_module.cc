////////////////////////////////////////////////////////////////////////
// Class:       CRTTimingAnalysis
// Plugin Type: analyzer
// File:        CRTTimingAnalysis_module.cc
// Author:      Henry Lay (h.lay@sheffield.ac.uk)
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
#include "canvas/Persistency/Common/FindOneP.h"

#include "TTree.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/Common/Reco/CorrectedOpFlashTiming.h"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoService.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

namespace sbnd::crt {
  class CRTTimingAnalysis;
}

class sbnd::crt::CRTTimingAnalysis : public art::EDAnalyzer {
public:
  explicit CRTTimingAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTimingAnalysis(CRTTimingAnalysis const&) = delete;
  CRTTimingAnalysis(CRTTimingAnalysis&&) = delete;
  CRTTimingAnalysis& operator=(CRTTimingAnalysis const&) = delete;
  CRTTimingAnalysis& operator=(CRTTimingAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec);

  void AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec);

  void SortReferencing();

  void AnalyseCRTSpacePoints(const std::vector<art::Ptr<CRTSpacePoint>> &CRTSpacePointVec,
                             const art::FindOneP<CRTCluster> &spacepointsToClusters,
                             const art::FindManyP<CRTStripHit> &clustersToStripHits);

  void AnalyseCRTTracks(const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec,
                        const art::FindManyP<CRTSpacePoint> &tracksToSpacePoints,
                        const art::FindOneP<CRTCluster> &spacepointsToClusters,
                        const art::FindManyP<CRTStripHit> &clustersToStripHits);
    
  void AnalyseTPCSlices(const std::vector<art::Ptr<recob::Slice>> &TPCSliceVec,
                        const art::FindManyP<sbn::CorrectedOpFlashTiming> &sliceToCorrectedOpFlashes,
                        const art::FindManyP<recob::PFParticle> &sliceToPFPs,
                        const art::FindOneP<recob::Track> &pfpToTrack,
                        const art::FindOneP<CRTSpacePoint, anab::T0> &trackToCRTSpacePoint,
                        const art::FindOneP<CRTCluster> &spacepointsToClusters,
                        const art::FindManyP<CRTStripHit> &clustersToStripHits);

private:

  art::ServiceHandle<CRTGeoService> fCRTGeoService;

  // fcl Controlled Variables
  std::string fCRTClusterModuleLabel, fCRTSpacePointModuleLabel, fCRTTrackModuleLabel,
    fPTBModuleLabel, fTDCModuleLabel, fTimingReferenceModuleLabel, fTPCSliceModuleLabel,
    fCorrectedOpFlashModuleLabel, fTPCTrackModuleLabel, fCRTSpacePointMatchingModuleLabel;
  std::vector<uint32_t> fAllowedPTBHLTs;

  // Global Storage
  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;

  std::vector<uint64_t> _ptb_llt_trigger;
  std::vector<uint64_t> _ptb_llt_timestamp;

  std::vector<uint32_t>    _tdc_channel;
  std::vector<uint64_t>    _tdc_timestamp;
  std::vector<uint64_t>    _tdc_offset;
  std::vector<std::string> _tdc_name;

  // Trees
  TTree *fSPTree, *fTrTree, *fTPCTree;

  // Tree Variables
  int _run;
  int _subrun;
  int _event;
  int _crt_timing_reference_type;
  int _crt_timing_reference_channel;

  bool   _etrig_good;
  bool   _rwm_good;
  bool   _ptb_hlt_beam_gate_good;
  bool   _crt_t1_reset_good;
  double _rwm_etrig_diff;
  double _ptb_hlt_beam_gate_etrig_diff;
  double _rwm_crt_t1_reset_diff;
  double _ptb_hlt_beam_gate_crt_t1_reset_diff;
  double _rwm_ptb_hlt_beam_gate_diff;

  uint16_t             _sp_nhits;
  int16_t              _sp_tagger;
  double               _sp_x;
  double               _sp_y;
  double               _sp_z;
  double               _sp_pe;
  double               _sp_ts0;
  double               _sp_ts0_rwm_ref;
  double               _sp_ts0_ptb_hlt_beam_gate_ref;
  double               _sp_dts0;
  double               _sp_ts1;
  double               _sp_ts1_rwm_ref;
  double               _sp_ts1_ptb_hlt_beam_gate_ref;
  double               _sp_dts1;
  bool                 _sp_single_timing_chain;
  int16_t              _sp_timing_chain;
  std::vector<int32_t> _sp_sh_channel_set;
  std::vector<int16_t> _sp_sh_mac5_set;
  std::vector<int16_t> _sp_sh_timing_chain_set;
  std::vector<double>  _sp_sh_ts0_set;
  std::vector<double>  _sp_sh_ts1_set;
  std::vector<double>  _sp_sh_time_walk_set;
  std::vector<double>  _sp_sh_prop_delay_set;
  std::vector<double>  _sp_sh_cable_length_set;
  std::vector<double>  _sp_sh_calib_offset_ts0_set;
  std::vector<double>  _sp_sh_calib_offset_ts1_set;

  double  _tr_start_x;
  double  _tr_start_y;
  double  _tr_start_z;
  double  _tr_end_x;
  double  _tr_end_y;
  double  _tr_end_z;
  double  _tr_dir_x;
  double  _tr_dir_y;
  double  _tr_dir_z;
  double  _tr_ts0;
  double  _tr_ts0_rwm_ref;
  double  _tr_ts0_ptb_hlt_beam_gate_ref;
  double  _tr_ts1;
  double  _tr_ts1_rwm_ref;
  double  _tr_ts1_ptb_hlt_beam_gate_ref;
  double  _tr_pe;
  double  _tr_length;
  double  _tr_length_tof;
  double  _tr_tof_ts0;
  double  _tr_tof_diff_ts0;
  double  _tr_tof_ts1;
  double  _tr_tof_diff_ts1;
  double  _tr_theta;
  double  _tr_phi;
  bool    _tr_triple;
  int16_t _tr_tagger1;
  int16_t _tr_tagger2;
  int16_t _tr_tagger3;
  int16_t _tr_start_tagger;
  double  _tr_start_dts0;
  bool    _tr_start_single_timing_chain;
  int16_t _tr_start_timing_chain;
  int16_t _tr_end_tagger;
  double  _tr_end_dts0;
  bool    _tr_end_single_timing_chain;
  int16_t _tr_end_timing_chain;

  bool   _tpc_has_corrected_opflash;
  bool   _tpc_has_crt_sp_match;
  double _tpc_opflash_t0;
  double _tpc_opflash_nutof_light;
  double _tpc_opflash_nutof_charge;
  double _tpc_opflash_t0_corrected;
  double _tpc_opflash_t0_corrected_rwm;
  double _tpc_crt_sp_score;
  double _tpc_crt_sp_ts0;
  double _tpc_crt_sp_ts0_rwm_ref;
  double _tpc_crt_sp_ts0_ptb_hlt_beam_gate_ref;
  double _tpc_crt_sp_dts0;
  double _tpc_crt_sp_ts1;
  double _tpc_crt_sp_ts1_rwm_ref;
  double _tpc_crt_sp_ts1_ptb_hlt_beam_gate_ref;
  double _tpc_crt_sp_dts1;
};

sbnd::crt::CRTTimingAnalysis::CRTTimingAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fCRTClusterModuleLabel            = p.get<std::string>("CRTClusterModuleLabel");
  fCRTSpacePointModuleLabel         = p.get<std::string>("CRTSpacePointModuleLabel");
  fCRTTrackModuleLabel              = p.get<std::string>("CRTTrackModuleLabel");
  fPTBModuleLabel                   = p.get<std::string>("PTBModuleLabel");
  fTDCModuleLabel                   = p.get<std::string>("TDCModuleLabel");
  fTimingReferenceModuleLabel       = p.get<std::string>("TimingReferenceModuleLabel");
  fTPCSliceModuleLabel              = p.get<std::string>("TPCSliceModuleLabel");
  fCorrectedOpFlashModuleLabel      = p.get<std::string>("CorrectedOpFlashModuleLabel");
  fTPCTrackModuleLabel              = p.get<std::string>("TPCTrackModuleLabel");
  fCRTSpacePointMatchingModuleLabel = p.get<std::string>("CRTSpacePointMatchingModuleLabel");
  fAllowedPTBHLTs                   = p.get<std::vector<uint32_t>>("AllowedPTBHLTs");

  art::ServiceHandle<art::TFileService> fs;

  fSPTree = fs->make<TTree>("spacepoints","");
  fSPTree->Branch("run", "int", &_run);
  fSPTree->Branch("subrun", "int", &_subrun);
  fSPTree->Branch("event", "int", &_event);
  fSPTree->Branch("crt_timing_reference_type", "int", &_crt_timing_reference_type);
  fSPTree->Branch("crt_timing_reference_channel", "int", &_crt_timing_reference_channel);

  fSPTree->Branch("etrig_good", "bool", &_etrig_good);
  fSPTree->Branch("rwm_good", "bool", &_rwm_good);
  fSPTree->Branch("ptb_hlt_beam_gate_good", "bool", &_ptb_hlt_beam_gate_good);
  fSPTree->Branch("crt_t1_reset_good", "bool", &_crt_t1_reset_good);
  fSPTree->Branch("rwm_etrig_diff", "double", &_rwm_etrig_diff);
  fSPTree->Branch("ptb_hlt_beam_gate_etrig_diff", "double", &_ptb_hlt_beam_gate_etrig_diff);
  fSPTree->Branch("rwm_crt_t1_reset_diff", "double", &_rwm_crt_t1_reset_diff);
  fSPTree->Branch("ptb_hlt_beam_gate_crt_t1_reset_diff", "double", &_ptb_hlt_beam_gate_crt_t1_reset_diff);
  fSPTree->Branch("rwm_ptb_hlt_beam_gate_diff", "double", &_rwm_ptb_hlt_beam_gate_diff);

  fSPTree->Branch("sp_nhits", "uint16_t", &_sp_nhits);
  fSPTree->Branch("sp_tagger", "int16_t", &_sp_tagger);
  fSPTree->Branch("sp_x", "double", &_sp_x);
  fSPTree->Branch("sp_y", "double", &_sp_y);
  fSPTree->Branch("sp_z", "double", &_sp_z);
  fSPTree->Branch("sp_pe", "double", &_sp_pe);
  fSPTree->Branch("sp_ts0", "double", &_sp_ts0);
  fSPTree->Branch("sp_ts0_rwm_ref", "double", &_sp_ts0_rwm_ref);
  fSPTree->Branch("sp_ts0_ptb_hlt_beam_gate_ref", "double", &_sp_ts0_ptb_hlt_beam_gate_ref);
  fSPTree->Branch("sp_dts0", "double", &_sp_dts0);
  fSPTree->Branch("sp_ts1", "double", &_sp_ts1);
  fSPTree->Branch("sp_ts1_rwm_ref", "double", &_sp_ts1_rwm_ref);
  fSPTree->Branch("sp_ts1_ptb_hlt_beam_gate_ref", "double", &_sp_ts1_ptb_hlt_beam_gate_ref);
  fSPTree->Branch("sp_dts1", "double", &_sp_dts1);
  fSPTree->Branch("sp_single_timing_chain", "bool", &_sp_single_timing_chain);
  fSPTree->Branch("sp_timing_chain", "int16_t", &_sp_timing_chain);
  fSPTree->Branch("sp_sh_channel_set", "std::vector<int32_t>", &_sp_sh_channel_set);
  fSPTree->Branch("sp_sh_mac5_set", "std::vector<int16_t>", &_sp_sh_mac5_set);
  fSPTree->Branch("sp_sh_timing_chain_set", "std::vector<int16_t>", &_sp_sh_timing_chain_set);
  fSPTree->Branch("sp_sh_ts0_set", "std::vector<double>", &_sp_sh_ts0_set);
  fSPTree->Branch("sp_sh_ts1_set", "std::vector<double>", &_sp_sh_ts1_set);
  fSPTree->Branch("sp_sh_time_walk_set", "std::vector<double>", &_sp_sh_time_walk_set);
  fSPTree->Branch("sp_sh_prop_delay_set", "std::vector<double>", &_sp_sh_prop_delay_set);
  fSPTree->Branch("sp_sh_cable_length_set", "std::vector<double>", &_sp_sh_cable_length_set);
  fSPTree->Branch("sp_sh_calib_offset_ts0_set", "std::vector<double>", &_sp_sh_calib_offset_ts0_set);
  fSPTree->Branch("sp_sh_calib_offset_ts1_set", "std::vector<double>", &_sp_sh_calib_offset_ts1_set);

  fTrTree = fs->make<TTree>("tracks","");
  fTrTree->Branch("run", "int", &_run);
  fTrTree->Branch("subrun", "int", &_subrun);
  fTrTree->Branch("event", "int", &_event);
  fTrTree->Branch("crt_timing_reference_type", "int", &_crt_timing_reference_type);
  fTrTree->Branch("crt_timing_reference_channel", "int", &_crt_timing_reference_channel);

  fTrTree->Branch("etrig_good", "bool", &_etrig_good);
  fTrTree->Branch("rwm_good", "bool", &_rwm_good);
  fTrTree->Branch("ptb_hlt_beam_gate_good", "bool", &_ptb_hlt_beam_gate_good);
  fTrTree->Branch("crt_t1_reset_good", "bool", &_crt_t1_reset_good);
  fTrTree->Branch("rwm_etrig_diff", "double", &_rwm_etrig_diff);
  fTrTree->Branch("ptb_hlt_beam_gate_etrig_diff", "double", &_ptb_hlt_beam_gate_etrig_diff);
  fTrTree->Branch("rwm_crt_t1_reset_diff", "double", &_rwm_crt_t1_reset_diff);
  fTrTree->Branch("ptb_hlt_beam_gate_crt_t1_reset_diff", "double", &_ptb_hlt_beam_gate_crt_t1_reset_diff);
  fTrTree->Branch("rwm_ptb_hlt_beam_gate_diff", "double", &_rwm_ptb_hlt_beam_gate_diff);

  fTrTree->Branch("tr_start_x", "double", &_tr_start_x);
  fTrTree->Branch("tr_start_y", "double", &_tr_start_y);
  fTrTree->Branch("tr_start_z", "double", &_tr_start_z);
  fTrTree->Branch("tr_end_x", "double", &_tr_end_x);
  fTrTree->Branch("tr_end_y", "double", &_tr_end_y);
  fTrTree->Branch("tr_end_z", "double", &_tr_end_z);
  fTrTree->Branch("tr_dir_x", "double", &_tr_dir_x);
  fTrTree->Branch("tr_dir_y", "double", &_tr_dir_y);
  fTrTree->Branch("tr_dir_z", "double", &_tr_dir_z);
  fTrTree->Branch("tr_ts0", "double", &_tr_ts0);
  fTrTree->Branch("tr_ts0_rwm_ref", "double", &_tr_ts0_rwm_ref);
  fTrTree->Branch("tr_ts0_ptb_hlt_beam_gate_ref", "double", &_tr_ts0_ptb_hlt_beam_gate_ref);
  fTrTree->Branch("tr_ts1", "double", &_tr_ts1);
  fTrTree->Branch("tr_ts1_rwm_ref", "double", &_tr_ts1_rwm_ref);
  fTrTree->Branch("tr_ts1_ptb_hlt_beam_gate_ref", "double", &_tr_ts1_ptb_hlt_beam_gate_ref);
  fTrTree->Branch("tr_pe", "double", &_tr_pe);
  fTrTree->Branch("tr_length", "double", &_tr_length);
  fTrTree->Branch("tr_length_tof", "double", &_tr_length_tof);
  fTrTree->Branch("tr_tof_ts0", "double", &_tr_tof_ts0);
  fTrTree->Branch("tr_tof_diff_ts0", "double", &_tr_tof_diff_ts0);
  fTrTree->Branch("tr_tof_ts1", "double", &_tr_tof_ts1);
  fTrTree->Branch("tr_tof_diff_ts1", "double", &_tr_tof_diff_ts1);
  fTrTree->Branch("tr_theta", "double", &_tr_theta);
  fTrTree->Branch("tr_phi", "double", &_tr_phi);
  fTrTree->Branch("tr_triple", "bool", &_tr_triple);
  fTrTree->Branch("tr_tagger1", "int16_t", &_tr_tagger1);
  fTrTree->Branch("tr_tagger2", "int16_t", &_tr_tagger2);
  fTrTree->Branch("tr_tagger3", "int16_t", &_tr_tagger3);
  fTrTree->Branch("tr_start_tagger", "int16_t", &_tr_start_tagger);
  fTrTree->Branch("tr_start_dts0", "double", &_tr_start_dts0);
  fTrTree->Branch("tr_start_single_timing_chain", "bool", &_tr_start_single_timing_chain);
  fTrTree->Branch("tr_start_timing_chain", "int16_t", &_tr_start_timing_chain);
  fTrTree->Branch("tr_end_tagger", "int16_t", &_tr_end_tagger);
  fTrTree->Branch("tr_end_dts0", "double", &_tr_end_dts0);
  fTrTree->Branch("tr_end_single_timing_chain", "bool", &_tr_end_single_timing_chain);
  fTrTree->Branch("tr_end_timing_chain", "int16_t", &_tr_end_timing_chain);

  fTPCTree = fs->make<TTree>("slices","");
  fTPCTree->Branch("run", "int", &_run);
  fTPCTree->Branch("subrun", "int", &_subrun);
  fTPCTree->Branch("event", "int", &_event);
  fTPCTree->Branch("crt_timing_reference_type", "int", &_crt_timing_reference_type);
  fTPCTree->Branch("crt_timing_reference_channel", "int", &_crt_timing_reference_channel);

  fTPCTree->Branch("etrig_good", "bool", &_etrig_good);
  fTPCTree->Branch("rwm_good", "bool", &_rwm_good);
  fTPCTree->Branch("ptb_hlt_beam_gate_good", "bool", &_ptb_hlt_beam_gate_good);
  fTPCTree->Branch("crt_t1_reset_good", "bool", &_crt_t1_reset_good);
  fTPCTree->Branch("rwm_etrig_diff", "double", &_rwm_etrig_diff);
  fTPCTree->Branch("ptb_hlt_beam_gate_etrig_diff", "double", &_ptb_hlt_beam_gate_etrig_diff);
  fTPCTree->Branch("rwm_crt_t1_reset_diff", "double", &_rwm_crt_t1_reset_diff);
  fTPCTree->Branch("ptb_hlt_beam_gate_crt_t1_reset_diff", "double", &_ptb_hlt_beam_gate_crt_t1_reset_diff);
  fTPCTree->Branch("rwm_ptb_hlt_beam_gate_diff", "double", &_rwm_ptb_hlt_beam_gate_diff);

  fTPCTree->Branch("tpc_has_corrected_opflash", "bool", &_tpc_has_corrected_opflash);
  fTPCTree->Branch("tpc_has_crt_sp_match", "bool", &_tpc_has_crt_sp_match);
  fTPCTree->Branch("tpc_opflash_t0", "double", &_tpc_opflash_t0);
  fTPCTree->Branch("tpc_opflash_nutof_light", "double", &_tpc_opflash_nutof_light);
  fTPCTree->Branch("tpc_opflash_nutof_charge", "double", &_tpc_opflash_nutof_charge);
  fTPCTree->Branch("tpc_opflash_t0_corrected", "double", &_tpc_opflash_t0_corrected);
  fTPCTree->Branch("tpc_opflash_t0_corrected_rwm", "double", &_tpc_opflash_t0_corrected_rwm);
  fTPCTree->Branch("tpc_crt_sp_score", "double", &_tpc_crt_sp_score);
  fTPCTree->Branch("tpc_crt_sp_ts0", "double", &_tpc_crt_sp_ts0);
  fTPCTree->Branch("tpc_crt_sp_ts0_rwm_ref", "double", &_tpc_crt_sp_ts0_rwm_ref);
  fTPCTree->Branch("tpc_crt_sp_ts0_ptb_hlt_beam_gate_ref", "double", &_tpc_crt_sp_ts0_ptb_hlt_beam_gate_ref);
  fTPCTree->Branch("tpc_crt_sp_dts0", "double", &_tpc_crt_sp_dts0);
  fTPCTree->Branch("tpc_crt_sp_ts1", "double", &_tpc_crt_sp_ts1);
  fTPCTree->Branch("tpc_crt_sp_ts1_rwm_ref", "double", &_tpc_crt_sp_ts1_rwm_ref);
  fTPCTree->Branch("tpc_crt_sp_ts1_ptb_hlt_beam_gate_ref", "double", &_tpc_crt_sp_ts1_ptb_hlt_beam_gate_ref);
  fTPCTree->Branch("tpc_crt_sp_dts1", "double", &_tpc_crt_sp_dts1);
}

void sbnd::crt::CRTTimingAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  _crt_timing_reference_type    = -1;
  _crt_timing_reference_channel = -1;

  art::Handle<raw::TimingReferenceInfo> TimingReferenceHandle;
  e.getByLabel(fTimingReferenceModuleLabel, TimingReferenceHandle);
  if(TimingReferenceHandle.isValid())
    {
      _crt_timing_reference_type    = TimingReferenceHandle->timingType;
      _crt_timing_reference_channel = TimingReferenceHandle->timingChannel;
    }

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

  SortReferencing();
  
  // Get CRTSpacePoints
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  // Get CRTSpacePoint to CRTCluster Assns
  art::FindOneP<CRTCluster> spacepointsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  // Get CRTClusters
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get CRTCluster to CRTStripHit Assns
  art::FindManyP<CRTStripHit> clustersToStripHits(CRTClusterHandle, e, fCRTClusterModuleLabel);

  // Fill CRTSpacePoint variables
  AnalyseCRTSpacePoints(CRTSpacePointVec, spacepointsToClusters, clustersToStripHits);

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

  // Fill CRTTrack variables
  AnalyseCRTTracks(CRTTrackVec, tracksToSpacePoints, spacepointsToClusters, clustersToStripHits);

  // Get TPCSlices
  art::Handle<std::vector<recob::Slice>> TPCSliceHandle;
  e.getByLabel(fTPCSliceModuleLabel, TPCSliceHandle);
  if(!TPCSliceHandle.isValid()){
    std::cout << "TPCSlice product " << fTPCSliceModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<recob::Slice>> TPCSliceVec;
  art::fill_ptr_vector(TPCSliceVec, TPCSliceHandle);

  // Get TPCSlice to CorrectedOpFlash Assns
  art::FindManyP<sbn::CorrectedOpFlashTiming> sliceToCorrectedOpFlashes(TPCSliceHandle, e, fCorrectedOpFlashModuleLabel);

  // Get TPCSlice to PFP Assns
  art::FindManyP<recob::PFParticle> sliceToPFPs(TPCSliceHandle, e, fTPCSliceModuleLabel);
  
  // Get TPCPFPs
  art::Handle<std::vector<recob::PFParticle>> TPCPFPHandle;
  e.getByLabel(fTPCSliceModuleLabel, TPCPFPHandle);
  if(!TPCPFPHandle.isValid()){
    std::cout << "TPCPFP product " << fTPCSliceModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get PFP to Track Assns
  art::FindOneP<recob::Track> pfpToTrack(TPCPFPHandle, e, fTPCTrackModuleLabel);

  // Get TPCTracks
  art::Handle<std::vector<recob::Track>> TPCTrackHandle;
  e.getByLabel(fTPCTrackModuleLabel, TPCTrackHandle);
  if(!TPCTrackHandle.isValid()){
    std::cout << "TPCTrack product " << fTPCTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Track to CRTSpacePoint Assns
  art::FindOneP<CRTSpacePoint, anab::T0> trackToCRTSpacePoint(TPCTrackHandle, e, fCRTSpacePointMatchingModuleLabel);

  AnalyseTPCSlices(TPCSliceVec, sliceToCorrectedOpFlashes, sliceToPFPs, pfpToTrack, trackToCRTSpacePoint, spacepointsToClusters, clustersToStripHits);
}

void sbnd::crt::CRTTimingAnalysis::AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec)
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
          _ptb_hlt_timestamp[hlt_i] = ptb->GetHLTrigger(i).timestamp * 20;

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
          _ptb_llt_timestamp[llt_i] = ptb->GetLLTrigger(i).timestamp * 20;

          ++llt_i;
        }
    }
}

void sbnd::crt::CRTTimingAnalysis::AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec)
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

void sbnd::crt::CRTTimingAnalysis::SortReferencing()
{
  _etrig_good = false; _rwm_good = false; _ptb_hlt_beam_gate_good = false; _crt_t1_reset_good = false;
  _rwm_etrig_diff = std::numeric_limits<double>::max(); _ptb_hlt_beam_gate_etrig_diff = std::numeric_limits<double>::max();
  _rwm_crt_t1_reset_diff = std::numeric_limits<double>::max(); _ptb_hlt_beam_gate_crt_t1_reset_diff = std::numeric_limits<double>::max();
  _rwm_ptb_hlt_beam_gate_diff = std::numeric_limits<double>::max();

  int etrig_count = 0, etrig_id = -1, rwm_count = 0, rwm_id = -1, crt_t1_reset_count = 0, crt_t1_reset_id = -1;

  for(unsigned int tdc_i = 0; tdc_i < _tdc_channel.size(); ++tdc_i)
    {
      if(_tdc_channel[tdc_i] == 4)
        {
          ++etrig_count;
          etrig_id = tdc_i;
        }
      else if(_tdc_channel[tdc_i] == 2)
        {
          ++rwm_count;
          rwm_id = tdc_i;
        }
      else if(_tdc_channel[tdc_i] == 0)
        {
          ++crt_t1_reset_count;
          crt_t1_reset_id = tdc_i;
        }
    }

  uint64_t etrig = std::numeric_limits<uint64_t>::max(), rwm = std::numeric_limits<uint64_t>::max(),
    hlt = std::numeric_limits<uint64_t>::max(), crt_t1_reset = std::numeric_limits<uint64_t>::max();

  if(etrig_count == 1)
    {
      _etrig_good  = true;
      etrig = _tdc_timestamp[etrig_id];
    }

  if(rwm_count == 1)
    {
      _rwm_good  = true;
      rwm = _tdc_timestamp[rwm_id];
    }

  if(etrig_count == 1)
    {
      double closest_diff = std::numeric_limits<double>::max();

      for(unsigned int ptb_i = 0; ptb_i < _ptb_hlt_trigger.size(); ++ptb_i)
        {
          std::bitset<32> hlt_bitmask = std::bitset<32>(_ptb_hlt_trigger[ptb_i]);

          for(uint32_t allowed_hlt : fAllowedPTBHLTs)
            {
              if(hlt_bitmask[allowed_hlt])
                {
                  _ptb_hlt_beam_gate_good = true;

                  uint64_t temp_hlt = _ptb_hlt_timestamp[ptb_i];
                  double diff  = etrig > temp_hlt ? etrig - temp_hlt : -1. * (temp_hlt - etrig);

                  if(std::abs(diff) < closest_diff)
                    {
                      closest_diff = diff;
                      hlt          = temp_hlt;
                    }
                }
            }
        }
    }

  if(crt_t1_reset_count == 1)
    {
      _crt_t1_reset_good  = true;
      crt_t1_reset = _tdc_timestamp[crt_t1_reset_id];
    }

  if(_etrig_good && _rwm_good)
    _rwm_etrig_diff = etrig > rwm ? etrig - rwm : -1. * (rwm - etrig);

  if(_etrig_good && _ptb_hlt_beam_gate_good)
    _ptb_hlt_beam_gate_etrig_diff = etrig > hlt ? etrig - hlt : -1. * (hlt - etrig);

  if(_crt_t1_reset_good && _rwm_good)
    _rwm_crt_t1_reset_diff = crt_t1_reset > rwm ? crt_t1_reset - rwm : -1. * (rwm - crt_t1_reset);

  if(_etrig_good && _crt_t1_reset_good && _ptb_hlt_beam_gate_good)
    _ptb_hlt_beam_gate_crt_t1_reset_diff = crt_t1_reset > hlt ? crt_t1_reset - hlt : -1. * (hlt - crt_t1_reset);

  if(_etrig_good && _rwm_good && _ptb_hlt_beam_gate_good)
    _rwm_ptb_hlt_beam_gate_diff = hlt > rwm ? hlt - rwm : -1. * (rwm - hlt);
}

void sbnd::crt::CRTTimingAnalysis::AnalyseCRTSpacePoints(const std::vector<art::Ptr<CRTSpacePoint>> &CRTSpacePointVec,
                                                         const art::FindOneP<CRTCluster> &spacepointsToClusters,
                                                         const art::FindManyP<CRTStripHit> &clustersToStripHits)
{
}

void sbnd::crt::CRTTimingAnalysis::AnalyseCRTTracks(const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec,
                                                    const art::FindManyP<CRTSpacePoint> &tracksToSpacePoints,
                                                    const art::FindOneP<CRTCluster> &spacepointsToClusters,
                                                    const art::FindManyP<CRTStripHit> &clustersToStripHits)
{
}

void sbnd::crt::CRTTimingAnalysis::AnalyseTPCSlices(const std::vector<art::Ptr<recob::Slice>> &TPCSliceVec,
                                                    const art::FindManyP<sbn::CorrectedOpFlashTiming> &sliceToCorrectedOpFlashes,
                                                    const art::FindManyP<recob::PFParticle> &sliceToPFPs,
                                                    const art::FindOneP<recob::Track> &pfpToTrack,
                                                    const art::FindOneP<CRTSpacePoint, anab::T0> &trackToCRTSpacePoint,
                                                    const art::FindOneP<CRTCluster> &spacepointsToClusters,
                                                    const art::FindManyP<CRTStripHit> &clustersToStripHits)
{
}

DEFINE_ART_MODULE(sbnd::crt::CRTTimingAnalysis)
