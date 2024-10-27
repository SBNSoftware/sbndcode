////////////////////////////////////////////////////////////////////////
// Class:       ADRIFT
// Plugin Type: analyzer
// File:        ADRIFT_module.cc
//
// Author:      Henry Lay (h.lay@sheffield.ac.uk)
// 
// ADC
// Dynamic
// Range
// Improvement
// Facilitation
// Trees
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
#include "TH1D.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

namespace sbnd {
  namespace crt {
    class ADRIFT;
  }
}


class sbnd::crt::ADRIFT : public art::EDAnalyzer {
public:
  explicit ADRIFT(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ADRIFT(ADRIFT const&) = delete;
  ADRIFT(ADRIFT&&) = delete;
  ADRIFT& operator=(ADRIFT const&) = delete;
  ADRIFT& operator=(ADRIFT&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  CRTGeoAlg fCRTGeoAlg;

  std::string fFEBDataModuleLabel, fCRTStripHitModuleLabel, fCRTClusterModuleLabel,
    fCRTSpacePointModuleLabel, fCRTTrackModuleLabel;

  bool fOnly2HitSpacePoints;

  double fTrackAngleLimit;

  TTree* fChannelTree;

  int _channel, _gdml_id, _mac5, _raw_channel, _tagger;
  double _area, _y_average, _ped_calib, _ped_fit, _ped_peak,
    _sh_rate, _sp_rate, _tr_rate,
    _sh_peak_fit, _sh_peak_peak, _sh_sat_ratio_total, _sh_sat_ratio_peak,
    _sp_peak_fit, _sp_peak_peak, _sp_sat_ratio_total, _sp_sat_ratio_peak,
    _tr_peak_fit, _tr_peak_peak, _tr_sat_ratio_total, _tr_sat_ratio_peak,
    _tr_lim_angle_peak_fit, _tr_lim_angle_peak_peak, _tr_lim_angle_sat_ratio_total,
    _tr_lim_angle_sat_ratio_peak;
  bool _horizontal;

  std::map<uint, TH1D*> hADCPed, hADCSH, hADCSP, hADCTr, hADCTrLA;
};


sbnd::crt::ADRIFT::ADRIFT(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  {
    fFEBDataModuleLabel       = p.get<std::string>("FEBDataModuleLabel");
    fCRTStripHitModuleLabel   = p.get<std::string>("CRTStripHitModuleLabel");
    fCRTClusterModuleLabel    = p.get<std::string>("CRTClusterModuleLabel");
    fCRTSpacePointModuleLabel = p.get<std::string>("CRTSpacePointModuleLabel");
    fCRTTrackModuleLabel      = p.get<std::string>("CRTTrackModuleLabel");
    fOnly2HitSpacePoints      = p.get<bool>("Only2HitSpacePoints");
    fTrackAngleLimit          = p.get<double>("TrackAngleLimit");

    art::ServiceHandle<art::TFileService> fs;

    fChannelTree = fs->make<TTree>("channel_tree", "");
    fChannelTree->Branch("channel", &_channel);
    fChannelTree->Branch("gdml_id", &_gdml_id);
    fChannelTree->Branch("mac5", &_mac5);
    fChannelTree->Branch("raw_channel", &_raw_channel);
    fChannelTree->Branch("tagger", &_tagger);
    fChannelTree->Branch("area", &_area);
    fChannelTree->Branch("y_average", &_y_average);
    fChannelTree->Branch("horizontal", &_horizontal);
    fChannelTree->Branch("ped_calib", &_ped_calib);
    fChannelTree->Branch("ped_fit", &_ped_fit);
    fChannelTree->Branch("ped_peak", &_ped_peak);
    fChannelTree->Branch("sh_rate", &_sh_rate);
    fChannelTree->Branch("sp_rate", &_sp_rate);
    fChannelTree->Branch("tr_rate", &_tr_rate);
    fChannelTree->Branch("sh_peak_fit", &_sh_peak_fit);
    fChannelTree->Branch("sh_peak_peak", &_sh_peak_peak);
    fChannelTree->Branch("sh_sat_ratio_total", &_sh_sat_ratio_total);
    fChannelTree->Branch("sh_sat_ratio_peak", &_sh_sat_ratio_peak);
    fChannelTree->Branch("sp_peak_fit", &_sp_peak_fit);
    fChannelTree->Branch("sp_peak_peak", &_sp_peak_peak);
    fChannelTree->Branch("sp_sat_ratio_total", &_sp_sat_ratio_total);
    fChannelTree->Branch("sp_sat_ratio_peak", &_sp_sat_ratio_peak);
    fChannelTree->Branch("tr_peak_fit", &_tr_peak_fit);
    fChannelTree->Branch("tr_peak_peak", &_tr_peak_peak);
    fChannelTree->Branch("tr_sat_ratio_total", &_tr_sat_ratio_total);
    fChannelTree->Branch("tr_sat_ratio_peak", &_tr_sat_ratio_peak);
    fChannelTree->Branch("tr_lim_angle_peak_fit", &_tr_lim_angle_peak_fit);
    fChannelTree->Branch("tr_lim_angle_peak_peak", &_tr_lim_angle_peak_peak);
    fChannelTree->Branch("tr_lim_angle_sat_ratio_total", &_tr_lim_angle_sat_ratio_total);
    fChannelTree->Branch("tr_lim_angle_sat_ratio_peak", &_tr_lim_angle_sat_ratio_peak);

    for(int ch = 0; ch < 4480; ++ch)
      {
        hADCPed[ch]  = new TH1D(Form("hADCPed_Channel%i", ch), ";ADC;Readouts", 4300, -.5, 4299.5);
        hADCSH[ch]   = new TH1D(Form("hADCSH_Channel%i", ch), ";ADC;Strip Hits", 4300, -.5, 4299.5);
        hADCSP[ch]   = new TH1D(Form("hADCSP_Channel%i", ch), ";ADC;Space Points", 4300, -.5, 4299.5);
        hADCTr[ch]   = new TH1D(Form("hADCTr_Channel%i", ch), ";ADC;Tracks", 4300, -.5, 4299.5);
        hADCTrLA[ch] = new TH1D(Form("hADCTrLA_Channel%i", ch), ";ADC;Tracks", 4300, -.5, 4299.5);
      }
  }

void sbnd::crt::ADRIFT::analyze(art::Event const& e)
{
  // Get FEBDatas
  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  if(!FEBDataHandle.isValid()){
    std::cout << "FEBData product " << fFEBDataModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Get FEBData to CRTStripHit Assns
  art::FindManyP<CRTStripHit> febDatasToStripHits(FEBDataHandle, e, fCRTStripHitModuleLabel);

  for(const art::Ptr<FEBData> &feb_data : FEBDataVec)
    {
      const std::vector<art::Ptr<CRTStripHit>> strip_hits = febDatasToStripHits.at(feb_data.key());
      const uint m5 = feb_data->Mac5();

      std::set<int> mask_channels;

      for(const art::Ptr<CRTStripHit> &strip_hit : strip_hits)
        {
          mask_channels.insert(strip_hit->Channel());
          mask_channels.insert(strip_hit->Channel() + 1);

          if(strip_hit->Channel() % 2)
            std::cout << "ODD Strip Hit Channel Number" << std::endl;
        }

      for(int i = 0; i < 32; ++i)
        {
          const int ch = m5 * 32 + i;

          if(mask_channels.count(ch))
            continue;

          hADCPed[ch]->Fill(feb_data->ADC(i));
        }
    }

  // Get CRTStripHits
  art::Handle<std::vector<CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  if(!CRTStripHitHandle.isValid()){
    std::cout << "CRTStripHit product " << fCRTStripHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  for(const art::Ptr<CRTStripHit> &strip_hit : CRTStripHitVec)
    {
      if(strip_hit->Channel() % 2)
        std::cout << "ODD Strip Hit Channel Number" << std::endl;
      
      CRTSiPMGeo sipm1 = fCRTGeoAlg.GetSiPM(strip_hit->Channel());
      CRTSiPMGeo sipm2 = fCRTGeoAlg.GetSiPM(strip_hit->Channel() + 1);

      hADCSH[strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
      hADCSH[strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);
    }

  // Get CRTSpacePoints
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  // Get CRTClusters
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get CRTSpacePoint to CRTCluster Assns
  art::FindOneP<CRTCluster> spacepointsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  // Get CRTCluster to CRTStripHit Assns
  art::FindManyP<CRTStripHit> clustersToStripHits(CRTClusterHandle, e, fCRTClusterModuleLabel);

  for(const art::Ptr<CRTSpacePoint> &space_point : CRTSpacePointVec)
    {
      const art::Ptr<CRTCluster> cluster = spacepointsToClusters.at(space_point.key());

      if(fOnly2HitSpacePoints && cluster->NHits() > 2)
        continue;

      const std::vector<art::Ptr<CRTStripHit>> strip_hits = clustersToStripHits.at(cluster.key());

      for(const art::Ptr<CRTStripHit> &strip_hit : strip_hits)
        {
          if(strip_hit->Channel() % 2)
            std::cout << "ODD Strip Hit Channel Number" << std::endl;
      
          CRTSiPMGeo sipm1 = fCRTGeoAlg.GetSiPM(strip_hit->Channel());
          CRTSiPMGeo sipm2 = fCRTGeoAlg.GetSiPM(strip_hit->Channel() + 1);
          
          hADCSP[strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
          hADCSP[strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);
        }
    }

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

  for(const art::Ptr<CRTTrack> &track : CRTTrackVec)
    {
      const std::vector<art::Ptr<CRTSpacePoint>> space_points = tracksToSpacePoints.at(track.key());

      const geo::Vector_t tr_dir_vec = track->Direction();
      const TVector3 tr_dir = TVector3(tr_dir_vec.X(), tr_dir_vec.Y(), tr_dir_vec.Z());

      for(const art::Ptr<CRTSpacePoint> &space_point : space_points)
        {
          const art::Ptr<CRTCluster> cluster = spacepointsToClusters.at(space_point.key());

          const std::vector<art::Ptr<CRTStripHit>> strip_hits = clustersToStripHits.at(cluster.key());
          const CRTTagger tagger = cluster->Tagger();
          TVector3 normal;
          if(tagger == kBottomTagger || tagger == kTopLowTagger || tagger == kTopHighTagger)
            normal = TVector3(0, 1, 0);
          else if(tagger == kWestTagger || tagger == kEastTagger)
            normal = TVector3(1, 0, 0);
          else if(tagger == kSouthTagger || tagger == kNorthTagger)
            normal = TVector3(0, 0, 1);

          const double angle = TMath::RadToDeg() * normal.Angle(tr_dir);

          for(const art::Ptr<CRTStripHit> &strip_hit : strip_hits)
            {
              if(strip_hit->Channel() % 2)
                std::cout << "ODD Strip Hit Channel Number" << std::endl;
      
              CRTSiPMGeo sipm1 = fCRTGeoAlg.GetSiPM(strip_hit->Channel());
              CRTSiPMGeo sipm2 = fCRTGeoAlg.GetSiPM(strip_hit->Channel() + 1);
          
              hADCTr[strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
              hADCTr[strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);

              if(angle < fTrackAngleLimit)
                {
                  hADCTrLA[strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
                  hADCTrLA[strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);
                }
            }
        }
    }
}

void sbnd::crt::ADRIFT::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::crt::ADRIFT::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::crt::ADRIFT)
