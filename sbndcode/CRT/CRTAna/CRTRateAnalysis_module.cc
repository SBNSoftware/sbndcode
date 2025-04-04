////////////////////////////////////////////////////////////////////////
// Class:       CRTRateAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTRateAnalysis_module.cc
//
// Generated at Tue Jan  7 05:28:06 2025 by Henry Lay using cetskelgen
// from cetlib version 3.18.02.
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

#include "TTree.h"

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "artdaq-core/Data/RawEvent.hh"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTBlob.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

namespace sbnd {
  namespace crt {
    class CRTRateAnalysis;
  }
}


class sbnd::crt::CRTRateAnalysis : public art::EDAnalyzer {
public:
  explicit CRTRateAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTRateAnalysis(CRTRateAnalysis const&) = delete;
  CRTRateAnalysis(CRTRateAnalysis&&) = delete;
  CRTRateAnalysis& operator=(CRTRateAnalysis const&) = delete;
  CRTRateAnalysis& operator=(CRTRateAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void ResetEventVars();
  void ResetRawVars();
  void ResetSpacePointVars();
  void ResetBlobVars();

private:

  CRTGeoAlg fCRTGeoAlg;

  std::string fFEBDataModuleLabel, fCRTSpacePointModuleLabel,
    fCRTBlobModuleLabel, fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel;

  TTree *fEventTree, *fRawTree, *fSpacePointTree, *fBlobTree;

  int32_t  _run, _subrun, _event;
  uint32_t _event_header_ts;

  int16_t  _tagger, _module;
  uint16_t _max_channel;

  int16_t _n_hits;
  double _x, _y, _z, _pe;

  int16_t _n_total_sps, _n_bottom_sps, _n_bottom_sps_twohits, _n_south_sps,
    _n_north_sps, _n_west_sps, _n_east_sps, _n_toplow_sps, _n_tophigh_sps;
};


sbnd::crt::CRTRateAnalysis::CRTRateAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fCRTBlobModuleLabel(p.get<std::string>("CRTBlobModuleLabel"))
  , fDAQHeaderModuleLabel(p.get<std::string>("DAQHeaderModuleLabel"))
  , fDAQHeaderInstanceLabel(p.get<std::string>("DAQHeaderInstanceLabel"))
{
  art::ServiceHandle<art::TFileService> fs;

  fEventTree = fs->make<TTree>("events", "");
  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);
  fEventTree->Branch("event_header_ts", &_event_header_ts);

  fRawTree = fs->make<TTree>("readouts", "");
  fRawTree->Branch("run", &_run);
  fRawTree->Branch("subrun", &_subrun);
  fRawTree->Branch("event", &_event);
  fRawTree->Branch("event_header_ts", &_event_header_ts);
  fRawTree->Branch("tagger", &_tagger);
  fRawTree->Branch("module", &_module);
  fRawTree->Branch("max_channel", &_max_channel);

  fSpacePointTree = fs->make<TTree>("spacepoints", "");
  fSpacePointTree->Branch("run", &_run);
  fSpacePointTree->Branch("subrun", &_subrun);
  fSpacePointTree->Branch("event", &_event);
  fSpacePointTree->Branch("event_header_ts", &_event_header_ts);
  fSpacePointTree->Branch("tagger", &_tagger);
  fSpacePointTree->Branch("n_hits", &_n_hits);
  fSpacePointTree->Branch("x", &_x);
  fSpacePointTree->Branch("y", &_y);
  fSpacePointTree->Branch("z", &_z);
  fSpacePointTree->Branch("pe", &_pe);

  fBlobTree = fs->make<TTree>("blobs", "");
  fBlobTree->Branch("run", &_run);
  fBlobTree->Branch("subrun", &_subrun);
  fBlobTree->Branch("event", &_event);
  fBlobTree->Branch("event_header_ts", &_event_header_ts);
  fBlobTree->Branch("n_total_sps", &_n_total_sps);
  fBlobTree->Branch("n_bottom_sps", &_n_bottom_sps);
  fBlobTree->Branch("n_bottom_sps_twohits", &_n_bottom_sps_twohits);
  fBlobTree->Branch("n_south_sps", &_n_south_sps);
  fBlobTree->Branch("n_north_sps", &_n_north_sps);
  fBlobTree->Branch("n_west_sps", &_n_west_sps);
  fBlobTree->Branch("n_east_sps", &_n_east_sps);
  fBlobTree->Branch("n_toplow_sps", &_n_toplow_sps);
  fBlobTree->Branch("n_tophigh_sps", &_n_tophigh_sps);
  fBlobTree->Branch("pe", &_pe);
}

void sbnd::crt::CRTRateAnalysis::analyze(art::Event const& e)
{
  ResetEventVars();

  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  art::Handle<artdaq::detail::RawEventHeader> DAQHeaderHandle;
  e.getByLabel(fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel, DAQHeaderHandle);

  if(DAQHeaderHandle.isValid())
    {
      artdaq::RawEvent rawHeaderEvent = artdaq::RawEvent(*DAQHeaderHandle);
      uint64_t raw_ts = rawHeaderEvent.timestamp();
      _event_header_ts = raw_ts / static_cast<uint32_t>(1e9);
    }

  fEventTree->Fill();

  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);

  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  for(auto const& data : FEBDataVec)
    {
      ResetRawVars();
      _tagger = fCRTGeoAlg.AuxDetIndexToTaggerEnum(data->Mac5());
      _module = data->Mac5();

      int max_adc = -1, max_ch = -1;

      for(int ch = 0; ch < 32; ++ch)
        {
          int adc = data->ADC(ch);

          if(adc > max_adc)
            {
              max_adc = adc;
              max_ch  = ch;
            }
        }

      _max_channel = _module * 32 + max_ch;

      fRawTree->Fill();
    }

  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);

  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::FindOneP<CRTCluster> spacePointsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  for(auto const& spacePoint : CRTSpacePointVec)
    {
      ResetSpacePointVars();

      const art::Ptr<CRTCluster> cluster = spacePointsToClusters.at(spacePoint.key());

      _tagger = cluster->Tagger();
      _n_hits = cluster->NHits();
      _x      = spacePoint->X();
      _y      = spacePoint->Y();
      _z      = spacePoint->Z();
      _pe     = spacePoint->PE();

      fSpacePointTree->Fill();
    }

  art::Handle<std::vector<CRTBlob>> CRTBlobHandle;
  e.getByLabel(fCRTBlobModuleLabel, CRTBlobHandle);

  std::vector<art::Ptr<CRTBlob>> CRTBlobVec;
  art::fill_ptr_vector(CRTBlobVec, CRTBlobHandle);

  art::FindManyP<CRTSpacePoint> blobsToSpacePoints(CRTBlobHandle, e, fCRTBlobModuleLabel);

  for(auto const& blob : CRTBlobVec)
    {
      ResetBlobVars();

      _n_total_sps   = blob->TotalSpacePoints();
      _n_bottom_sps  = blob->SpacePointsInTagger(kBottomTagger);
      _n_south_sps   = blob->SpacePointsInTagger(kSouthTagger);
      _n_north_sps   = blob->SpacePointsInTagger(kNorthTagger);
      _n_west_sps    = blob->SpacePointsInTagger(kWestTagger);
      _n_east_sps    = blob->SpacePointsInTagger(kEastTagger);
      _n_toplow_sps  = blob->SpacePointsInTagger(kTopLowTagger);
      _n_tophigh_sps = blob->SpacePointsInTagger(kTopHighTagger);
      _pe            = blob->PE();

      _n_bottom_sps_twohits = 0;

      const std::vector<art::Ptr<CRTSpacePoint>> spacePoints = blobsToSpacePoints.at(blob.key());

      for(auto const& spacePoint : spacePoints)
        {
          const art::Ptr<CRTCluster> cluster = spacePointsToClusters.at(spacePoint.key());

          if(cluster->Tagger() == kBottomTagger && cluster->NHits() > 1)
            ++_n_bottom_sps_twohits;
        }

      fBlobTree->Fill();
    }
}

void sbnd::crt::CRTRateAnalysis::ResetEventVars()
{
  _run = -1; _subrun = -1; _event = -1;
  _event_header_ts = std::numeric_limits<uint32_t>::max();
}

void sbnd::crt::CRTRateAnalysis::ResetRawVars()
{
  _tagger = -1; _module = -1;
  _max_channel = std::numeric_limits<uint16_t>::max();
}

void sbnd::crt::CRTRateAnalysis::ResetSpacePointVars()
{
  _tagger = -1;

  _n_hits = -1;

  _x  = std::numeric_limits<double>::lowest();
  _y  = std::numeric_limits<double>::lowest();
  _z  = std::numeric_limits<double>::lowest();
  _pe = std::numeric_limits<double>::lowest();
}

void sbnd::crt::CRTRateAnalysis::ResetBlobVars()
{
  _n_total_sps          = -1;
  _n_bottom_sps         = -1;
  _n_bottom_sps_twohits = -1;
  _n_south_sps          = -1;
  _n_north_sps          = -1;
  _n_west_sps           = -1;
  _n_east_sps           = -1;
  _n_toplow_sps         = -1;
  _n_tophigh_sps        = -1;

  _pe = std::numeric_limits<double>::lowest();
}

DEFINE_ART_MODULE(sbnd::crt::CRTRateAnalysis)
