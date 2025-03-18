////////////////////////////////////////////////////////////////////////
// Class:       CRTTopHatAnalysis
// Plugin Type: analyzer
// File:        CRTTopHatAnalysis_module.cc
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

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

namespace sbnd::crt {
  class CRTTopHatAnalysis;
}

class sbnd::crt::CRTTopHatAnalysis : public art::EDAnalyzer {
public:
  explicit CRTTopHatAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTopHatAnalysis(CRTTopHatAnalysis const&) = delete;
  CRTTopHatAnalysis(CRTTopHatAnalysis&&) = delete;
  CRTTopHatAnalysis& operator=(CRTTopHatAnalysis const&) = delete;
  CRTTopHatAnalysis& operator=(CRTTopHatAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec);

  void AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec);

  void AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec);

  void AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec,
                          const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints);
  void AnalyseCRTTracks(const art::Event &e, const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec);

private:

  CRTGeoAlg fCRTGeoAlg;
  TPCGeoAlg fTPCGeoAlg;

  std::string fCRTStripHitModuleLabel, fCRTClusterModuleLabel, fCRTSpacePointModuleLabel,
    fCRTTrackModuleLabel, fPTBModuleLabel, fTDCModuleLabel, fTimingReferenceModuleLabel;
  bool fDebug, fCutT0;
  double fMinT0, fMaxT0;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;
  int _crt_timing_reference_type;
  int _crt_timing_reference_channel;

  //strip hit to select the strip which has ADC above threshold
  std::vector<uint32_t> _sh_channel;
  std::vector<int16_t>  _sh_tagger;
  std::vector<int64_t>  _sh_ts0;
  std::vector<int64_t>  _sh_ts1;
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

  //cluster from x-y coincidence for CRTSpacePoint
  std::vector<int64_t>  _cl_ts0;
  std::vector<int64_t>  _cl_ts1;
  std::vector<uint32_t> _cl_unixs;
  std::vector<uint16_t> _cl_nhits;
  std::vector<int16_t>  _cl_tagger;
  std::vector<uint8_t>  _cl_composition;
  std::vector<bool>     _cl_has_sp;
  std::vector<double>   _cl_sp_x;
  std::vector<double>   _cl_sp_ex;
  std::vector<double>   _cl_sp_y;
  std::vector<double>   _cl_sp_ey;
  std::vector<double>   _cl_sp_z;
  std::vector<double>   _cl_sp_ez;
  std::vector<double>   _cl_sp_pe;
  std::vector<double>   _cl_sp_ts0;
  std::vector<double>   _cl_sp_ets0;
  std::vector<double>   _cl_sp_ts1;
  std::vector<double>   _cl_sp_ets1;
  std::vector<bool>     _cl_sp_complete;

  //track level information
  std::vector<double>  _tr_start_x;
  std::vector<double>  _tr_start_y;
  std::vector<double>  _tr_start_z;
  std::vector<double>  _tr_end_x;
  std::vector<double>  _tr_end_y;
  std::vector<double>  _tr_end_z;
  std::vector<double>  _tr_dir_x;
  std::vector<double>  _tr_dir_y;
  std::vector<double>  _tr_dir_z;
  std::vector<double>  _tr_ts0;
  std::vector<double>  _tr_ets0;
  std::vector<double>  _tr_ts1;
  std::vector<double>  _tr_ets1;
  std::vector<double>  _tr_pe;
  std::vector<double>  _tr_length;
  std::vector<double>  _tr_tof;
  std::vector<double>  _tr_theta;
  std::vector<double>  _tr_phi;
  std::vector<bool>    _tr_triple;
  std::vector<int16_t> _tr_tagger1;
  std::vector<int16_t> _tr_tagger2;
  std::vector<int16_t> _tr_tagger3;

  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;

  std::vector<uint64_t> _ptb_llt_trigger;
  std::vector<uint64_t> _ptb_llt_timestamp;

  std::vector<uint32_t>    _tdc_channel;
  std::vector<uint64_t>    _tdc_timestamp;
  std::vector<uint64_t>    _tdc_offset;
  std::vector<std::string> _tdc_name;
};

sbnd::crt::CRTTopHatAnalysis::CRTTopHatAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
{
  fCRTStripHitModuleLabel     = p.get<std::string>("CRTStripHitModuleLabel", "crtstrips");
  fCRTClusterModuleLabel      = p.get<std::string>("CRTClusterModuleLabel", "crtclustering");
  fCRTSpacePointModuleLabel   = p.get<std::string>("CRTSpacePointModuleLabel", "crtspacepoints");
  fCRTTrackModuleLabel        = p.get<std::string>("CRTTrackModuleLabel", "crttracks");
  fPTBModuleLabel             = p.get<std::string>("PTBModuleLabel", "ptbdecoder");
  fTDCModuleLabel             = p.get<std::string>("TDCModuleLabel", "tdcdecoder");
  fTimingReferenceModuleLabel = p.get<std::string>("TimingReferenceModuleLabel", "crtstrips");
  fDebug                      = p.get<bool>("Debug", false);
  fCutT0                      = p.get<bool>("CutT0", false);
  fMinT0                      = p.get<double>("MinT0", std::numeric_limits<double>::min());
  fMaxT0                      = p.get<double>("MaxT0", std::numeric_limits<double>::max());

  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);
  fTree->Branch("event", &_event);
  fTree->Branch("crt_timing_reference_type", &_crt_timing_reference_type);
  fTree->Branch("crt_timing_reference_channel", &_crt_timing_reference_channel);

  fTree->Branch("sh_channel", "std::vector<uint32_t>", &_sh_channel);
  fTree->Branch("sh_tagger", "std::vector<int16_t>", &_sh_tagger);
  fTree->Branch("sh_ts0", "std::vector<int64_t>", &_sh_ts0);
  fTree->Branch("sh_ts1", "std::vector<int64_t>", &_sh_ts1);
  fTree->Branch("sh_unixs", "std::vector<uint32_t>", &_sh_unixs);
  fTree->Branch("sh_pos", "std::vector<double>", &_sh_pos);
  fTree->Branch("sh_err", "std::vector<double>", &_sh_err);
  fTree->Branch("sh_adc1", "std::vector<uint16_t>", &_sh_adc1);
  fTree->Branch("sh_adc2", "std::vector<uint16_t>", &_sh_adc2);
  fTree->Branch("sh_saturated1", "std::vector<bool>", &_sh_saturated1);
  fTree->Branch("sh_saturated2", "std::vector<bool>", &_sh_saturated2);

  fTree->Branch("cl_ts0", "std::vector<int64_t>", &_cl_ts0);
  fTree->Branch("cl_ts1", "std::vector<int64_t>", &_cl_ts1);
  fTree->Branch("cl_unixs", "std::vector<uint32_t>", &_cl_unixs);
  fTree->Branch("cl_nhits", "std::vector<uint16_t>", &_cl_nhits);
  fTree->Branch("cl_tagger", "std::vector<int16_t>", &_cl_tagger);
  fTree->Branch("cl_composition", "std::vector<uint8_t>", &_cl_composition);
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
  fTree->Branch("cl_sp_complete", "std::vector<bool>", &_cl_sp_complete);

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

  fTree->Branch("ptb_hlt_trigger", "std::vector<uint64_t>", &_ptb_hlt_trigger);
  fTree->Branch("ptb_hlt_timestamp", "std::vector<uint64_t>", &_ptb_hlt_timestamp);
  fTree->Branch("ptb_llt_trigger", "std::vector<uint64_t>", &_ptb_llt_trigger);
  fTree->Branch("ptb_llt_timestamp", "std::vector<uint64_t>", &_ptb_llt_timestamp);

  fTree->Branch("tdc_channel", "std::vector<uint32_t>", &_tdc_channel);
  fTree->Branch("tdc_timestamp", "std::vector<uint64_t>", &_tdc_timestamp);
  fTree->Branch("tdc_offset", "std::vector<uint64_t>", &_tdc_offset);
  fTree->Branch("tdc_name", "std::vector<std::string>", &_tdc_name);  
}

void sbnd::crt::CRTTopHatAnalysis::analyze(art::Event const& e)
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

  // Fill CRTCluster variables
  AnalyseCRTClusters(e, CRTClusterVec, clustersToSpacePoints);

  // Get CRTTracks
  art::Handle<std::vector<CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  if(!CRTTrackHandle.isValid()){
    std::cout << "CRTTrack product " << fCRTTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTTrack>> CRTTrackVec;
  art::fill_ptr_vector(CRTTrackVec, CRTTrackHandle);

  // Fill CRTTrack variables
  AnalyseCRTTracks(e, CRTTrackVec);

  fTree->Fill();
}

void sbnd::crt::CRTTopHatAnalysis::AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec)
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

void sbnd::crt::CRTTopHatAnalysis::AnalyseTDCs(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &TDCVec)
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

void sbnd::crt::CRTTopHatAnalysis::AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec)
{
  _sh_channel.clear();
  _sh_tagger.clear();
  _sh_ts0.clear();
  _sh_ts1.clear();
  _sh_unixs.clear();
  _sh_pos.clear();
  _sh_err.clear();
  _sh_adc1.clear();
  _sh_adc2.clear();
  _sh_saturated1.clear();
  _sh_saturated2.clear();

  for(auto const &hit : CRTStripHitVec)
    {
      if(fCutT0 && (hit->Ts0() < fMinT0 || hit->Ts0() > fMaxT0))
        continue;

      _sh_channel.push_back(hit->Channel());
      _sh_tagger.push_back(fCRTGeoAlg.ChannelToTaggerEnum(hit->Channel()));
      _sh_ts0.push_back(hit->Ts0());
      _sh_ts1.push_back(hit->Ts1());
      _sh_unixs.push_back(hit->UnixS());
      _sh_pos.push_back(hit->Pos());
      _sh_err.push_back(hit->Error());
      _sh_adc1.push_back(hit->ADC1());
      _sh_adc2.push_back(hit->ADC2());
      _sh_saturated1.push_back(hit->Saturated1());
      _sh_saturated2.push_back(hit->Saturated2());
    }
}

void sbnd::crt::CRTTopHatAnalysis::AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec,
                                                      const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints)
{
  _cl_ts0.clear();
  _cl_ts1.clear();
  _cl_unixs.clear();
  _cl_nhits.clear();
  _cl_tagger.clear();
  _cl_composition.clear();
  _cl_has_sp.clear();
  _cl_sp_x.clear();
  _cl_sp_ex.clear();
  _cl_sp_y.clear();
  _cl_sp_ey.clear();
  _cl_sp_z.clear();
  _cl_sp_ez.clear();
  _cl_sp_pe.clear();
  _cl_sp_ts0.clear();
  _cl_sp_ets0.clear();
  _cl_sp_ts1.clear();
  _cl_sp_ets1.clear();
  _cl_sp_complete.clear();

  for(auto const &cluster : CRTClusterVec)
    {
      if(fCutT0 && (cluster->Ts0() < fMinT0 || cluster->Ts0() > fMaxT0))
        continue;

      _cl_ts0.push_back(cluster->Ts0());
      _cl_ts1.push_back(cluster->Ts1());
      _cl_unixs.push_back(cluster->UnixS());
      _cl_nhits.push_back(cluster->NHits());
      _cl_tagger.push_back(cluster->Tagger());
      _cl_composition.push_back(cluster->Composition());

      const auto spacepoints = clustersToSpacePoints.at(cluster.key());
      if(spacepoints.size() == 1)
        { 
          const auto spacepoint = spacepoints[0];
          
          _cl_has_sp.push_back(true);
          _cl_sp_x.push_back(spacepoint->X());
          _cl_sp_ex.push_back(spacepoint->XErr());
          _cl_sp_y.push_back(spacepoint->Y());
          _cl_sp_ey.push_back(spacepoint->YErr());
          _cl_sp_z.push_back(spacepoint->Z());
          _cl_sp_ez.push_back(spacepoint->ZErr());
          _cl_sp_pe.push_back(spacepoint->PE());
          _cl_sp_ts0.push_back(spacepoint->Ts0());
          _cl_sp_ets0.push_back(spacepoint->Ts0Err());
          _cl_sp_ts1.push_back(spacepoint->Ts1());
          _cl_sp_ets1.push_back(spacepoint->Ts1Err());
          _cl_sp_complete.push_back(spacepoint->Complete());
        }
      else
        {
          _cl_has_sp.push_back(false);
          _cl_sp_x.push_back(-999999.);
          _cl_sp_ex.push_back(-999999.);
          _cl_sp_y.push_back(-999999.);
          _cl_sp_ey.push_back(-999999.);
          _cl_sp_z.push_back(-999999.);
          _cl_sp_ez.push_back(-999999.);
          _cl_sp_pe.push_back(-999999.);
          _cl_sp_ts0.push_back(-999999.);
          _cl_sp_ets0.push_back(-999999.);
          _cl_sp_ts1.push_back(-999999.);
          _cl_sp_ets1.push_back(-999999.);
          _cl_sp_complete.push_back(false);
        }
    }
}
void sbnd::crt::CRTTopHatAnalysis::AnalyseCRTTracks(const art::Event &e, const std::vector<art::Ptr<CRTTrack>> &CRTTrackVec)
{
  _tr_start_x.clear();
  _tr_start_y.clear();
  _tr_start_z.clear();
  _tr_end_x.clear();
  _tr_end_y.clear();
  _tr_end_z.clear();
  _tr_dir_x.clear();
  _tr_dir_y.clear();
  _tr_dir_z.clear();
  _tr_ts0.clear();
  _tr_ets0.clear();
  _tr_ts1.clear();
  _tr_ets1.clear();
  _tr_pe.clear();
  _tr_length.clear();
  _tr_tof.clear();
  _tr_theta.clear();
  _tr_phi.clear();
  _tr_triple.clear();
  _tr_tagger1.clear();
  _tr_tagger2.clear();
  _tr_tagger3.clear();

  for(auto const& track : CRTTrackVec)
    {
      if(fCutT0 && (track->Ts0() < fMinT0 || track->Ts0() > fMaxT0))
        continue;

      const geo::Point_t start = track->Start();
      _tr_start_x.push_back(start.X());
      _tr_start_y.push_back(start.Y());
      _tr_start_z.push_back(start.Z());

      const geo::Point_t end = track->End();
      _tr_end_x.push_back(end.X());
      _tr_end_y.push_back(end.Y());
      _tr_end_z.push_back(end.Z());

      const geo::Vector_t dir = track->Direction();
      _tr_dir_x.push_back(dir.X());
      _tr_dir_y.push_back(dir.Y());
      _tr_dir_z.push_back(dir.Z());

      _tr_ts0.push_back(track->Ts0());
      _tr_ets0.push_back(track->Ts0Err());
      _tr_ts1.push_back(track->Ts1());
      _tr_ets1.push_back(track->Ts1Err());
      _tr_pe.push_back(track->PE());
      _tr_length.push_back(track->Length());
      _tr_tof.push_back(track->ToF());
      _tr_theta.push_back(TMath::RadToDeg() * track->Theta());
      _tr_phi.push_back(TMath::RadToDeg() * track->Phi());
      _tr_triple.push_back(track->Triple());

      unsigned tag_i = 0;

      for(auto const &tagger : track->Taggers())
        {
          if(tag_i == 0)
            _tr_tagger1.push_back(tagger);
          else if(tag_i == 1)
            _tr_tagger2.push_back(tagger);
          else if(tag_i == 2)
            _tr_tagger3.push_back(tagger);

          ++tag_i;
        }
    }
}

DEFINE_ART_MODULE(sbnd::crt::CRTTopHatAnalysis)
