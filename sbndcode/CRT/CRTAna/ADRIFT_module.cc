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
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/ChannelMaps/CRT/CRTChannelMapService.h"

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

  void PedestalFit(TH1D* hADC, double &fit, double &chi2, bool &converged);
  void PedestalPeak(TH1D* hADC, double &peak);
  void Rate(TH1D* hADC, double &rate);
  void PeakPeak(TH1D* hADC, const double &ped, double &peak);
  void PeakFit(TH1D* hADC, const double &peak, const double &ped, double &fit, double &chi2, bool &converged);
  double Saturation(TH1D* hADC);
  void ResetVars();
  void SaveHist(TH1D *hADC, std::string &saveDir, std::string saveName, int rebin = 1);
  static double LanGau(double *x, double *par);

  CRTGeoAlg fCRTGeoAlg;

  // Inputs
  std::string fFEBDataModuleLabel, fCRTStripHitModuleLabel, fCRTClusterModuleLabel,
    fCRTSpacePointModuleLabel, fCRTTrackModuleLabel, fTopSaveDirectory;

  bool fOnly2HitSpacePoints, fSaveAllFits, fSaveBadFits, fSaveSubset, fTrackLA;

  double fTrackAngleLimit, fPullWindow;

  // Other Global Parameters
  int fNEvents;

  TTree* fChannelTree;

  int _channel, _gdml_id, _mac5, _raw_channel, _tagger;
  double _area, _y_average, _ped_calib, _ped_fit, _ped_fit_chi2, _ped_peak,
    _sh_rate, _sp_rate, _tr_rate,
    _sh_peak_fit, _sh_peak_fit_chi2, _sh_peak_peak, _sh_sat_ratio_total, _sh_sat_ratio_peak,
    _sp_peak_fit, _sp_peak_fit_chi2, _sp_peak_peak, _sp_sat_ratio_total, _sp_sat_ratio_peak,
    _tr_peak_fit, _tr_peak_fit_chi2, _tr_peak_peak, _tr_sat_ratio_total, _tr_sat_ratio_peak,
    _tr_lim_angle_peak_fit, _tr_lim_angle_peak_fit_chi2, _tr_lim_angle_peak_peak, _tr_lim_angle_sat_ratio_total,
    _tr_lim_angle_sat_ratio_peak;
  bool _horizontal, _ped_fit_converged, _sh_peak_fit_converged, _sp_peak_fit_converged,
    _tr_peak_fit_converged, _tr_lim_angle_peak_fit_converged;

  std::map<uint, TH1D*> hADCPed, hADCSH, hADCSP, hADCTr, hADCTrLA;

  std::string fPedestalSaveDirectory, fPeakSaveDirectory, fBadPedestalSaveDirectory, fBadPeakSaveDirectory,
    fPedestalSubsetSaveDirectory, fStripHitSubsetSaveDirectory, fSpacePointSubsetSaveDirectory, fTrackSubsetSaveDirectory,
    fTrackLASubsetSaveDirectory;
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
  fTopSaveDirectory         = p.get<std::string>("TopSaveDirectory");
  fOnly2HitSpacePoints      = p.get<bool>("Only2HitSpacePoints");
  fSaveAllFits              = p.get<bool>("SaveAllFits");
  fSaveBadFits              = p.get<bool>("SaveBadFits");
  fSaveSubset               = p.get<bool>("SaveSubset");
  fTrackLA                  = p.get<bool>("TrackLA");
  fTrackAngleLimit          = p.get<double>("TrackAngleLimit");
  fPullWindow               = p.get<double>("PullWindow");

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
  fChannelTree->Branch("ped_fit_chi2", &_ped_fit_chi2);
  fChannelTree->Branch("ped_fit_converged", &_ped_fit_converged);
  fChannelTree->Branch("ped_peak", &_ped_peak);
  fChannelTree->Branch("sh_rate", &_sh_rate);
  fChannelTree->Branch("sp_rate", &_sp_rate);
  fChannelTree->Branch("tr_rate", &_tr_rate);
  fChannelTree->Branch("sh_peak_fit", &_sh_peak_fit);
  fChannelTree->Branch("sh_peak_fit_chi2", &_sh_peak_fit_chi2);
  fChannelTree->Branch("sh_peak_fit_converged", &_sh_peak_fit_converged);
  fChannelTree->Branch("sh_peak_peak", &_sh_peak_peak);
  fChannelTree->Branch("sh_sat_ratio_total", &_sh_sat_ratio_total);
  fChannelTree->Branch("sh_sat_ratio_peak", &_sh_sat_ratio_peak);
  fChannelTree->Branch("sp_peak_fit", &_sp_peak_fit);
  fChannelTree->Branch("sp_peak_fit_chi2", &_sp_peak_fit_chi2);
  fChannelTree->Branch("sp_peak_fit_converged", &_sp_peak_fit_converged);
  fChannelTree->Branch("sp_peak_peak", &_sp_peak_peak);
  fChannelTree->Branch("sp_sat_ratio_total", &_sp_sat_ratio_total);
  fChannelTree->Branch("sp_sat_ratio_peak", &_sp_sat_ratio_peak);
  fChannelTree->Branch("tr_peak_fit", &_tr_peak_fit);
  fChannelTree->Branch("tr_peak_fit_chi2", &_tr_peak_fit_chi2);
  fChannelTree->Branch("tr_peak_fit_converged", &_tr_peak_fit_converged);
  fChannelTree->Branch("tr_peak_peak", &_tr_peak_peak);
  fChannelTree->Branch("tr_sat_ratio_total", &_tr_sat_ratio_total);
  fChannelTree->Branch("tr_sat_ratio_peak", &_tr_sat_ratio_peak);
  if(fTrackLA)
    {
      fChannelTree->Branch("tr_lim_angle_peak_fit", &_tr_lim_angle_peak_fit);
      fChannelTree->Branch("tr_lim_angle_peak_fit_chi2", &_tr_lim_angle_peak_fit_chi2);
      fChannelTree->Branch("tr_lim_angle_peak_fit_converged", &_tr_lim_angle_peak_fit_converged);
      fChannelTree->Branch("tr_lim_angle_peak_peak", &_tr_lim_angle_peak_peak);
      fChannelTree->Branch("tr_lim_angle_sat_ratio_total", &_tr_lim_angle_sat_ratio_total);
      fChannelTree->Branch("tr_lim_angle_sat_ratio_peak", &_tr_lim_angle_sat_ratio_peak);
    }

  for(int ch = 0; ch < 4480; ++ch)
    {
      hADCPed[ch]  = new TH1D(Form("hADCPed_Channel%i", ch), ";ADC;Readouts", 1000, -.5, 999.5);
      hADCSH[ch]   = new TH1D(Form("hADCSH_Channel%i", ch), ";ADC;Strip Hits", 4300, -.5, 4299.5);
      hADCSP[ch]   = new TH1D(Form("hADCSP_Channel%i", ch), ";ADC;Space Points", 4300, -.5, 4299.5);
      hADCTr[ch]   = new TH1D(Form("hADCTr_Channel%i", ch), ";ADC;Tracks", 4300, -.5, 4299.5);

      if(fTrackLA)
        hADCTrLA[ch] = new TH1D(Form("hADCTrLA_Channel%i", ch), ";ADC;Tracks", 4300, -.5, 4299.5);
    }
}

void sbnd::crt::ADRIFT::analyze(art::Event const& e)
{
  if((fSaveAllFits || fSaveBadFits) && fNEvents == 0)
    {
      gSystem->Exec(Form("mkdir -p %s", fTopSaveDirectory.c_str()));

      fPedestalSaveDirectory = Form("%s/run%i/pedestals", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fPedestalSaveDirectory.c_str()));

      fPeakSaveDirectory = Form("%s/run%i/peaks", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fPeakSaveDirectory.c_str()));

      fBadPedestalSaveDirectory = Form("%s/run%i/pedestals/bad", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fBadPedestalSaveDirectory.c_str()));

      fBadPeakSaveDirectory = Form("%s/run%i/peaks/bad", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fBadPeakSaveDirectory.c_str()));
    }

  if(fSaveSubset && fNEvents == 0)
    {
      gSystem->Exec(Form("mkdir -p %s", fTopSaveDirectory.c_str()));

      fPedestalSubsetSaveDirectory = Form("%s/run%i/subset/pedestals", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fPedestalSubsetSaveDirectory.c_str()));

      fStripHitSubsetSaveDirectory = Form("%s/run%i/subset/striphits", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fStripHitSubsetSaveDirectory.c_str()));

      fSpacePointSubsetSaveDirectory = Form("%s/run%i/subset/spacepoints", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fSpacePointSubsetSaveDirectory.c_str()));

      fTrackSubsetSaveDirectory = Form("%s/run%i/subset/tracks", fTopSaveDirectory.c_str(), e.run());
      gSystem->Exec(Form("mkdir -p %s", fTrackSubsetSaveDirectory.c_str()));

      if(fTrackLA)
        {
          fTrackLASubsetSaveDirectory = Form("%s/run%i/subset/tracks_limited_angle", fTopSaveDirectory.c_str(), e.run());
          gSystem->Exec(Form("mkdir -p %s", fTrackLASubsetSaveDirectory.c_str()));
        }
    }

  if((fSaveAllFits || fSaveBadFits || fSaveSubset) && fNEvents == 0)
    {
      gStyle->SetOptStat(0);
      gStyle->SetFrameLineWidth(2);
      gStyle->SetTextFont(62);
      gStyle->SetTextSize(0.07);
      gStyle->SetLabelFont(62, "xyz");
      gStyle->SetLabelSize(0.05, "xyz");
      gStyle->SetTitleSize(0.05, "xyz");
      gStyle->SetTitleFont(62, "xyz");
    }

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

              if(fTrackLA && angle < fTrackAngleLimit)
                {
                  hADCTrLA[strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
                  hADCTrLA[strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);
                }
            }
        }
    }

  ++fNEvents;
}

void sbnd::crt::ADRIFT::beginJob()
{
  fNEvents = 0;
}

void sbnd::crt::ADRIFT::endJob()
{
  art::ServiceHandle<SBND::CRTChannelMapService> ChannelMapService;

  for(int gdml_i = 0; gdml_i < 140; ++gdml_i)
    {
      ResetVars();

      std::cout << "=== Processing module " << gdml_i << std::endl;

      SBND::CRTChannelMapService::ModuleInfo_t moduleInfo = ChannelMapService->GetModuleInfoFromOfflineID(gdml_i);
      uint mac5 = moduleInfo.valid ? moduleInfo.feb_mac5 : 0;
      bool invert = moduleInfo.valid ? moduleInfo.channel_order_swapped : false;

      for(int ch_i = 0; ch_i < 32; ++ch_i)
        {
          const int ch = gdml_i * 32 + ch_i;

          if(fSaveSubset && (ch > 1471 && ch < 1728))
            {
              SaveHist(hADCPed[ch], fPedestalSubsetSaveDirectory, Form("pedestal_channel_%i", ch), 2);
              SaveHist(hADCSH[ch], fStripHitSubsetSaveDirectory, Form("strip_hit_channel_%i", ch), 20);
              SaveHist(hADCSP[ch], fSpacePointSubsetSaveDirectory, Form("space_point_channel_%i", ch), 20);
              SaveHist(hADCTr[ch], fTrackSubsetSaveDirectory, Form("track_channel_%i", ch), 20);

              if(fTrackLA)
                SaveHist(hADCTrLA[ch], fTrackLASubsetSaveDirectory, Form("track_limited_angle_channel_%i", ch), 20);
            }

          _channel     = ch;
          _gdml_id     = gdml_i;
          _mac5        = mac5;
          _raw_channel = invert ? 31 - ch_i : ch_i;
          _tagger      = fCRTGeoAlg.ChannelToTaggerEnum(ch);
          _area        = fCRTGeoAlg.StripArea(ch);
          _y_average   = fCRTGeoAlg.StripAverageY(ch);
          _ped_calib   = fCRTGeoAlg.GetSiPM(ch).pedestal;
          _horizontal  = _tagger == kBottomTagger || _tagger == kTopLowTagger || _tagger == kTopHighTagger;

          PedestalFit(hADCPed[ch], _ped_fit, _ped_fit_chi2, _ped_fit_converged);
          PedestalPeak(hADCPed[ch], _ped_peak);

          Rate(hADCSH[ch], _sh_rate);
          Rate(hADCSP[ch], _sp_rate);

          if(_tagger != kBottomTagger)
            Rate(hADCTr[ch], _tr_rate);

          PeakPeak(hADCSH[ch], _ped_calib, _sh_peak_peak);
          PeakFit(hADCSH[ch], _sh_peak_peak, _ped_calib, _sh_peak_fit, _sh_peak_fit_chi2, _sh_peak_fit_converged);
          PeakPeak(hADCSP[ch], _ped_calib, _sp_peak_peak);
          PeakFit(hADCSP[ch], _sp_peak_peak, _ped_calib, _sp_peak_fit, _sp_peak_fit_chi2, _sp_peak_fit_converged);

          if(_tagger != kBottomTagger)
            {
              PeakPeak(hADCTr[ch], _ped_calib, _tr_peak_peak);
              PeakFit(hADCTr[ch], _tr_peak_peak, _ped_calib, _tr_peak_fit, _tr_peak_fit_chi2, _tr_peak_fit_converged);

              if(fTrackLA)
                {
                  PeakPeak(hADCTrLA[ch], _ped_calib, _tr_lim_angle_peak_peak);
                  PeakFit(hADCTrLA[ch], _tr_lim_angle_peak_peak, _ped_calib, _tr_lim_angle_peak_fit, _tr_lim_angle_peak_fit_chi2, _tr_lim_angle_peak_fit_converged);
                }
            }

          const double sh_sat = Saturation(hADCSH[ch]);
          _sh_sat_ratio_total = sh_sat / hADCSH[ch]->GetEntries();
          _sh_sat_ratio_peak  = sh_sat / hADCSH[ch]->GetBinContent(hADCSH[ch]->FindBin(_sh_peak_peak));

          const double sp_sat = Saturation(hADCSP[ch]);
          _sp_sat_ratio_total = sp_sat / hADCSP[ch]->GetEntries();
          _sp_sat_ratio_peak  = sp_sat / hADCSP[ch]->GetBinContent(hADCSP[ch]->FindBin(_sp_peak_peak));

          if(_tagger != kBottomTagger)
            {
              const double tr_sat = Saturation(hADCTr[ch]);
              _tr_sat_ratio_total = tr_sat / hADCTr[ch]->GetEntries();
              _tr_sat_ratio_peak  = tr_sat / hADCTr[ch]->GetBinContent(hADCTr[ch]->FindBin(_tr_peak_peak));

              if(fTrackLA)
                {
                  const double tr_lim_angle_sat = Saturation(hADCTrLA[ch]);
                  _tr_lim_angle_sat_ratio_total = tr_lim_angle_sat / hADCTrLA[ch]->GetEntries();
                  _tr_lim_angle_sat_ratio_peak  = tr_lim_angle_sat / hADCTrLA[ch]->GetBinContent(hADCTrLA[ch]->FindBin(_tr_lim_angle_peak_peak));
                }
            }

          fChannelTree->Fill();
        }
    }
}

void sbnd::crt::ADRIFT::PedestalFit(TH1D* hADC, double &fit, double &chi2, bool &converged)
{
  TF1 *gaus = new TF1("gaus", "gaus", 0, 500);
  const TString name = hADC->GetName();
  TString ch_name = name;
  ch_name.Remove(0,15);

  TH1D* hADC2 = (TH1D*) hADC->Clone(name + "_for_fit");
  hADC2->Rebin(5);

  TFitResultPtr fit_result = hADC2->Fit(gaus, "QRS");
  converged = !(bool)(int(fit_result));

  if(!converged)
    std::cout << "Pedestal fit has not converged - " << hADC->GetName() << std::endl;

  fit  = gaus->GetParameter("Mean");
  chi2 = gaus->GetChisquare() / gaus->GetNDF();

  if(fSaveAllFits || (fSaveBadFits && !converged))
    {
      TCanvas *c = new TCanvas(Form("c%s", name.Data()), Form("c%s", name.Data()));
      c->cd();

      hADC2->SetLineColor(kOrange+2);
      hADC2->SetLineWidth(2);
      hADC2->Draw("histe");
      gaus->SetLineColor(kSpring-6);
      gaus->Draw("same");

      if(converged)
        {
          c->SaveAs(Form("%s/pedestal_fit_channel_%s.png", fPedestalSaveDirectory.c_str(), ch_name.Data()));
          c->SaveAs(Form("%s/pedestal_fit_channel_%s.pdf", fPedestalSaveDirectory.c_str(), ch_name.Data()));
        }
      else
        {
          c->SaveAs(Form("%s/pedestal_fit_channel_%s.png", fBadPedestalSaveDirectory.c_str(), ch_name.Data()));
          c->SaveAs(Form("%s/pedestal_fit_channel_%s.pdf", fBadPedestalSaveDirectory.c_str(), ch_name.Data()));
        }
    }
}

void sbnd::crt::ADRIFT::PedestalPeak(TH1D* hADC, double &peak)
{
  const int bin = hADC->GetMaximumBin();

  peak = hADC->GetBinCenter(bin);
}

void sbnd::crt::ADRIFT::Rate(TH1D* hADC, double &rate)
{
  rate = hADC->GetEntries() / (fNEvents * fPullWindow);
}

void sbnd::crt::ADRIFT::PeakPeak(TH1D* hADC, const double &ped, double &peak)
{
  int bin = 0;
  double max = std::numeric_limits<double>::lowest();

  for(int i = (int)ped + 20; i < 4000; ++i)
    {
      if(hADC->GetBinContent(i) > max)
        {
          max = hADC->GetBinContent(i);
          bin = i;
        }
    }

  peak = hADC->GetBinCenter(bin);
}

void sbnd::crt::ADRIFT::PeakFit(TH1D* hADC, const double &peak, const double &ped, double &fit, double &chi2, bool &converged)
{
  TF1 *langau = new TF1("langau", LanGau, ped + 20, 4000, 4);
  double params[4] = { 10, peak, hADC->GetEntries(), 50 };
  langau->SetParameters(params);
  const TString name = hADC->GetName();
  TString ch_name    = name;
  TString type       = "";

  if(name.Contains("SH"))
    {
      ch_name.Remove(0,14);
      type = "strip_hit";
    }
  else if(name.Contains("SP"))
    {
      ch_name.Remove(0,14);
      type = "space_point";
    }
  else if(name.Contains("TrLA"))
    {
      ch_name.Remove(0,16);
      type = "track_limited_angle";
    }
  else if(name.Contains("Tr"))
    {
      ch_name.Remove(0,14);
      type = "track";
    }
  else
    std::cout << "Cannot identify histogram type (" << name << ")" << std::endl;

  TH1D* hADC2 = (TH1D*) hADC->Clone(name + "_for_fit");
  hADC2->Rebin(20);

  TFitResultPtr fit_result = hADC2->Fit(langau, "QRS");
  converged = !(bool)(int(fit_result));

  if(!converged)
    std::cout << "Peak fit has not converged - " << hADC->GetName() << std::endl;

  fit  = langau->GetParameter(1);
  chi2 = langau->GetChisquare() / langau->GetNDF();

  if(fSaveAllFits || (fSaveBadFits && !converged))
    {
      TCanvas *c = new TCanvas(Form("c%s", name.Data()), Form("c%s", name.Data()));
      c->cd();

      hADC2->SetLineColor(kOrange+2);
      hADC2->SetLineWidth(2);
      hADC2->Draw("histe");
      langau->SetLineColor(kSpring-6);
      langau->Draw("same");

      if(converged)
        {
          c->SaveAs(Form("%s/peak_fit_%s_channel_%s.png", fPeakSaveDirectory.c_str(), type.Data(), ch_name.Data()));
          c->SaveAs(Form("%s/peak_fit_%s_channel_%s.pdf", fPeakSaveDirectory.c_str(), type.Data(), ch_name.Data()));
        }
      else
        {
          c->SaveAs(Form("%s/peak_fit_%s_channel_%s.png", fBadPeakSaveDirectory.c_str(), type.Data(), ch_name.Data()));
          c->SaveAs(Form("%s/peak_fit_%s_channel_%s.pdf", fBadPeakSaveDirectory.c_str(), type.Data(), ch_name.Data()));
        }
    }
}

double sbnd::crt::ADRIFT::Saturation(TH1D* hADC)
{
  const double sat = hADC->GetBinContent(4090);

  for(int i = 4080; i < 4101; ++i)
    {
      if(i == 4090)
        continue;

      if(hADC->GetBinContent(i) > sat)
        std::cout << "Saturation appears to be in bin " << i << " (" << hADC->GetBinContent(i)
                  << ") not bin 4090 (" << sat << ")" << std::endl;
    }

  return sat;
}

void sbnd::crt::ADRIFT::ResetVars()
{
  _channel     = std::numeric_limits<int>::lowest();
  _gdml_id     = std::numeric_limits<int>::lowest();
  _mac5        = std::numeric_limits<int>::lowest();
  _raw_channel = std::numeric_limits<int>::lowest();
  _tagger      = std::numeric_limits<int>::lowest();

  _area                         = std::numeric_limits<double>::lowest();
  _y_average                    = std::numeric_limits<double>::lowest();
  _ped_calib                    = std::numeric_limits<double>::lowest();
  _ped_fit                      = std::numeric_limits<double>::lowest();
  _ped_fit_chi2                 = std::numeric_limits<double>::lowest();
  _ped_peak                     = std::numeric_limits<double>::lowest();
  _sh_rate                      = std::numeric_limits<double>::lowest();
  _sp_rate                      = std::numeric_limits<double>::lowest();
  _tr_rate                      = std::numeric_limits<double>::lowest();
  _sh_peak_fit                  = std::numeric_limits<double>::lowest();
  _sh_peak_fit_chi2             = std::numeric_limits<double>::lowest();
  _sh_peak_peak                 = std::numeric_limits<double>::lowest();
  _sh_sat_ratio_total           = std::numeric_limits<double>::lowest();
  _sh_sat_ratio_peak            = std::numeric_limits<double>::lowest();
  _sp_peak_fit                  = std::numeric_limits<double>::lowest();
  _sp_peak_fit_chi2             = std::numeric_limits<double>::lowest();
  _sp_peak_peak                 = std::numeric_limits<double>::lowest();
  _sp_sat_ratio_total           = std::numeric_limits<double>::lowest();
  _sp_sat_ratio_peak            = std::numeric_limits<double>::lowest();
  _tr_peak_fit                  = std::numeric_limits<double>::lowest();
  _tr_peak_fit_chi2             = std::numeric_limits<double>::lowest();
  _tr_peak_peak                 = std::numeric_limits<double>::lowest();
  _tr_sat_ratio_total           = std::numeric_limits<double>::lowest();
  _tr_sat_ratio_peak            = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_fit        = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_fit_chi2   = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_peak       = std::numeric_limits<double>::lowest();
  _tr_lim_angle_sat_ratio_total = std::numeric_limits<double>::lowest();
  _tr_lim_angle_sat_ratio_peak  = std::numeric_limits<double>::lowest();

  _horizontal                      = false;
  _ped_fit_converged               = false;
  _sh_peak_fit_converged           = false;
  _sp_peak_fit_converged           = false;
  _tr_peak_fit_converged           = false;
  _tr_lim_angle_peak_fit_converged = false;
}

void sbnd::crt::ADRIFT::SaveHist(TH1D *hADC, std::string &saveDir, std::string saveName, int rebin)
{
  TCanvas *c = new TCanvas(Form("c%s", saveName.c_str()), Form("c%s", saveName.c_str()));
  c->cd();

  const TString name = hADC->GetName();
  TH1D* hADC2 = (TH1D*) hADC->Clone(name + "_for_save");

  hADC2->SetLineColor(kOrange+2);
  hADC2->SetLineWidth(2);
  hADC2->Rebin(rebin);
  hADC2->Draw("histe");

  c->SaveAs(Form("%s/%s.png", saveDir.c_str(), saveName.c_str()));
  c->SaveAs(Form("%s/%s.pdf", saveDir.c_str(), saveName.c_str()));
}

double sbnd::crt::ADRIFT::LanGau(double *x, double *par)
{
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  double np = 100.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  double xx;
  double mpc;
  double fland;
  double sum = 0.0;
  double xlow,xupp;
  double step;
  double i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

DEFINE_ART_MODULE(sbnd::crt::ADRIFT)
