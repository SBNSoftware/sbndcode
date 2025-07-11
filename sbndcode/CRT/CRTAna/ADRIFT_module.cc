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
#include "TText.h"

#include "artdaq-core/Data/RawEvent.hh"

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

  void MakeSaveDirectories(art::Event const &e);
  void AnalyseFEBDatas(art::Event const &e, const int window);
  void AnalyseStripHits(art::Event const &e, const int window);
  void AnalyseSpacePoints(art::Event const &e, const int window);
  void AnalyseTracks(art::Event const &e, const int window);
  void ProcessEntry(const int ch, const int window);
  void PedestalFit(TH1D* hADC, double &fit, double &std, double &chi2, bool &converged,
                   bool badChannel, const int window);
  void PedestalPeak(TH1D* hADC, double &peak);
  void Rate(TH1D* hADC, double &rate);
  void PeakPeak(TH1D* hADC, const double &ped, double &peak);
  void PeakFit(TH1D* hADC, const double &peak, const double &ped, double &fit,
               double &chi2, bool &converged, bool badChannel, const int window);
  double Saturation(TH1D* hADC);
  void ResetVars();
  void SaveHist(TH1D *hADC, std::string &saveDir, std::string saveName, int rebin, bool badChannel);
  static double LanGau(double *x, double *par);

  CRTGeoAlg fCRTGeoAlg;

  // Inputs
  std::string fFEBDataModuleLabel, fCRTStripHitModuleLabel, fCRTClusterModuleLabel,
    fCRTSpacePointModuleLabel, fCRTTrackModuleLabel, fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel,
    fTopSaveDirectory;

  bool fOnly2HitSpacePoints, fSaveAllFits, fSaveBadFits, fSaveSubset, fFEBs, fStripHits, fSpacePoints,
    fTracks, fTrackLA, fSaveROOTHists, fAnalysePE;

  double fTrackAngleLimit, fPullWindow;

  uint32_t fRawTSCorrection;

  std::vector<std::pair<uint64_t, uint64_t>> fUnixWindows;

  // Other Global Parameters
  int fNEvents;

  TTree* fChannelTree;

  uint64_t _unix_start, _unix_end;
  int _channel, _gdml_id, _mac5, _raw_channel, _tagger, _channel_status;
  double _area, _y_average, _ped_calib, _gain_calib, _ped_fit, _ped_fit_std, _ped_fit_chi2, _ped_peak,
    _ped_reset_fit, _ped_reset_fit_std, _ped_reset_fit_chi2, _ped_reset_peak, _raw_max_chan_rate, _sh_rate, _sp_rate, _tr_rate,
    _sh_peak_fit, _sh_peak_fit_chi2, _sh_peak_peak, _sh_pe_peak_fit, _sh_pe_peak_fit_chi2, _sh_pe_peak_peak,
    _sh_sat_rate, _sh_sat_ratio_total, _sh_sat_ratio_peak, _sp_peak_fit, _sp_peak_fit_chi2, _sp_peak_peak,
    _sp_pe_peak_fit, _sp_pe_peak_fit_chi2, _sp_pe_peak_peak, _sp_sat_rate, _sp_sat_ratio_total, _sp_sat_ratio_peak,
    _tr_peak_fit, _tr_peak_fit_chi2, _tr_peak_peak, _tr_pe_peak_fit, _tr_pe_peak_fit_chi2, _tr_pe_peak_peak,
    _tr_sat_rate, _tr_sat_ratio_total, _tr_sat_ratio_peak, _tr_lim_angle_peak_fit, _tr_lim_angle_peak_fit_chi2, _tr_lim_angle_peak_peak,
    _tr_lim_angle_pe_peak_fit, _tr_lim_angle_pe_peak_fit_chi2, _tr_lim_angle_pe_peak_peak, _tr_lim_angle_sat_rate,
    _tr_lim_angle_sat_ratio_total, _tr_lim_angle_sat_ratio_peak, _sh_ped_to_peak_fit, _sh_ped_reset_to_peak_fit,
    _tr_by_length_peak_fit, _tr_by_length_peak_fit_chi2, _tr_by_length_peak_peak,
    _tr_by_length_pe_peak_fit, _tr_by_length_pe_peak_fit_chi2, _tr_by_length_pe_peak_peak;
  bool _horizontal, _ped_fit_converged, _ped_reset_fit_converged, _sh_peak_fit_converged, _sh_pe_peak_fit_converged,
    _sp_peak_fit_converged, _sp_pe_peak_fit_converged, _tr_peak_fit_converged, _tr_pe_peak_fit_converged, _tr_lim_angle_peak_fit_converged,
    _tr_lim_angle_pe_peak_fit_converged, _tr_by_length_peak_fit_converged, _tr_by_length_pe_peak_fit_converged;

  std::map<uint, std::map<uint, TH1D*>> hADCPed, hADCPedReset, hADCMaxChan, hADCSH, hPESH, hADCSP, hPESP, hADCTr, hPETr,
    hADCTrLA, hPETrLA, hADCTrByLength, hPETrByLength;

  std::map<uint, std::string> fPedestalSaveDirectory, fPedestalResetSaveDirectory, fPeakSaveDirectory, fBadPedestalSaveDirectory,
    fBadPedestalResetSaveDirectory, fBadPeakSaveDirectory, fPedestalSubsetSaveDirectory, fPedestalResetSubsetSaveDirectory,
    fStripHitSubsetSaveDirectory, fStripHitPESubsetSaveDirectory, fSpacePointSubsetSaveDirectory, fSpacePointPESubsetSaveDirectory,
    fTrackSubsetSaveDirectory, fTrackPESubsetSaveDirectory, fTrackLASubsetSaveDirectory, fTrackLAPESubsetSaveDirectory,
    fTrackByLengthSubsetSaveDirectory, fTrackByLengthPESubsetSaveDirectory;
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
  fDAQHeaderModuleLabel     = p.get<std::string>("DAQHeaderModuleLabel");
  fDAQHeaderInstanceLabel   = p.get<std::string>("DAQHeaderInstanceLabel");
  fTopSaveDirectory         = p.get<std::string>("TopSaveDirectory");
  fOnly2HitSpacePoints      = p.get<bool>("Only2HitSpacePoints");
  fSaveAllFits              = p.get<bool>("SaveAllFits");
  fSaveBadFits              = p.get<bool>("SaveBadFits");
  fSaveSubset               = p.get<bool>("SaveSubset");
  fFEBs                     = p.get<bool>("FEBs");
  fStripHits                = p.get<bool>("StripHits");
  fSpacePoints              = p.get<bool>("SpacePoints");
  fTracks                   = p.get<bool>("Tracks");
  fTrackLA                  = p.get<bool>("TrackLA");
  fSaveROOTHists            = p.get<bool>("SaveROOTHists");
  fAnalysePE                = p.get<bool>("AnalysePE");
  fTrackAngleLimit          = p.get<double>("TrackAngleLimit");
  fPullWindow               = p.get<double>("PullWindow");
  fRawTSCorrection          = p.get<uint32_t>("RawTSCorrection");
  fUnixWindows              = p.get<std::vector<std::pair<uint64_t, uint64_t>>>("UnixWindows", { {0, std::numeric_limits<uint64_t>::max()} });

  art::ServiceHandle<art::TFileService> fs;

  fChannelTree = fs->make<TTree>("channel_tree", "");
  fChannelTree->Branch("unix_start", &_unix_start);
  fChannelTree->Branch("unix_end", &_unix_end);
  fChannelTree->Branch("nevents", &fNEvents);
  fChannelTree->Branch("channel", &_channel);
  fChannelTree->Branch("gdml_id", &_gdml_id);
  fChannelTree->Branch("mac5", &_mac5);
  fChannelTree->Branch("raw_channel", &_raw_channel);
  fChannelTree->Branch("tagger", &_tagger);
  fChannelTree->Branch("channel_status", &_channel_status);
  fChannelTree->Branch("area", &_area);
  fChannelTree->Branch("y_average", &_y_average);
  fChannelTree->Branch("horizontal", &_horizontal);
  fChannelTree->Branch("ped_calib", &_ped_calib);
  fChannelTree->Branch("gain_calib", &_gain_calib);
  if(fFEBs)
    {
      fChannelTree->Branch("ped_fit", &_ped_fit);
      fChannelTree->Branch("ped_fit_std", &_ped_fit_std);
      fChannelTree->Branch("ped_fit_chi2", &_ped_fit_chi2);
      fChannelTree->Branch("ped_fit_converged", &_ped_fit_converged);
      fChannelTree->Branch("ped_peak", &_ped_peak);
      fChannelTree->Branch("ped_reset_fit", &_ped_reset_fit);
      fChannelTree->Branch("ped_reset_fit_std", &_ped_reset_fit_std);
      fChannelTree->Branch("ped_reset_fit_chi2", &_ped_reset_fit_chi2);
      fChannelTree->Branch("ped_reset_fit_converged", &_ped_reset_fit_converged);
      fChannelTree->Branch("ped_reset_peak", &_ped_reset_peak);
      fChannelTree->Branch("raw_max_chan_rate", &_raw_max_chan_rate);
    }
  if(fStripHits)
    {
      fChannelTree->Branch("sh_rate", &_sh_rate);
      fChannelTree->Branch("sh_peak_fit", &_sh_peak_fit);
      fChannelTree->Branch("sh_peak_fit_chi2", &_sh_peak_fit_chi2);
      fChannelTree->Branch("sh_peak_fit_converged", &_sh_peak_fit_converged);
      fChannelTree->Branch("sh_peak_peak", &_sh_peak_peak);
      fChannelTree->Branch("sh_sat_rate", &_sh_sat_rate);
      fChannelTree->Branch("sh_sat_ratio_total", &_sh_sat_ratio_total);
      fChannelTree->Branch("sh_sat_ratio_peak", &_sh_sat_ratio_peak);
      fChannelTree->Branch("sh_ped_to_peak_fit", &_sh_ped_to_peak_fit);
      fChannelTree->Branch("sh_ped_reset_to_peak_fit", &_sh_ped_reset_to_peak_fit);

      if(fAnalysePE)
        {
          fChannelTree->Branch("sh_pe_peak_fit", &_sh_pe_peak_fit);
          fChannelTree->Branch("sh_pe_peak_fit_chi2", &_sh_pe_peak_fit_chi2);
          fChannelTree->Branch("sh_pe_peak_fit_converged", &_sh_pe_peak_fit_converged);
          fChannelTree->Branch("sh_pe_peak_peak", &_sh_pe_peak_peak);
        }
    }
  if(fSpacePoints)
    {
      fChannelTree->Branch("sp_rate", &_sp_rate);
      fChannelTree->Branch("sp_peak_fit", &_sp_peak_fit);
      fChannelTree->Branch("sp_peak_fit_chi2", &_sp_peak_fit_chi2);
      fChannelTree->Branch("sp_peak_fit_converged", &_sp_peak_fit_converged);
      fChannelTree->Branch("sp_peak_peak", &_sp_peak_peak);
      fChannelTree->Branch("sp_sat_rate", &_sp_sat_rate);
      fChannelTree->Branch("sp_sat_ratio_total", &_sp_sat_ratio_total);
      fChannelTree->Branch("sp_sat_ratio_peak", &_sp_sat_ratio_peak);

      if(fAnalysePE)
        {
          fChannelTree->Branch("sp_pe_peak_fit", &_sp_pe_peak_fit);
          fChannelTree->Branch("sp_pe_peak_fit_chi2", &_sp_pe_peak_fit_chi2);
          fChannelTree->Branch("sp_pe_peak_fit_converged", &_sp_pe_peak_fit_converged);
          fChannelTree->Branch("sp_pe_peak_peak", &_sp_pe_peak_peak);
        }
    }
  if(fTracks)
    {
      fChannelTree->Branch("tr_rate", &_tr_rate);
      fChannelTree->Branch("tr_peak_fit", &_tr_peak_fit);
      fChannelTree->Branch("tr_peak_fit_chi2", &_tr_peak_fit_chi2);
      fChannelTree->Branch("tr_peak_fit_converged", &_tr_peak_fit_converged);
      fChannelTree->Branch("tr_peak_peak", &_tr_peak_peak);
      fChannelTree->Branch("tr_sat_rate", &_tr_sat_rate);
      fChannelTree->Branch("tr_sat_ratio_total", &_tr_sat_ratio_total);
      fChannelTree->Branch("tr_sat_ratio_peak", &_tr_sat_ratio_peak);
      fChannelTree->Branch("tr_by_length_peak_fit", &_tr_by_length_peak_fit);
      fChannelTree->Branch("tr_by_length_peak_fit_chi2", &_tr_by_length_peak_fit_chi2);
      fChannelTree->Branch("tr_by_length_peak_fit_converged", &_tr_by_length_peak_fit_converged);
      fChannelTree->Branch("tr_by_length_peak_peak", &_tr_by_length_peak_peak);

      if(fAnalysePE)
        {
          fChannelTree->Branch("tr_pe_peak_fit", &_tr_pe_peak_fit);
          fChannelTree->Branch("tr_pe_peak_fit_chi2", &_tr_pe_peak_fit_chi2);
          fChannelTree->Branch("tr_pe_peak_fit_converged", &_tr_pe_peak_fit_converged);
          fChannelTree->Branch("tr_pe_peak_peak", &_tr_pe_peak_peak);
          fChannelTree->Branch("tr_by_length_pe_peak_fit", &_tr_by_length_pe_peak_fit);
          fChannelTree->Branch("tr_by_length_pe_peak_fit_chi2", &_tr_by_length_pe_peak_fit_chi2);
          fChannelTree->Branch("tr_by_length_pe_peak_fit_converged", &_tr_by_length_pe_peak_fit_converged);
          fChannelTree->Branch("tr_by_length_pe_peak_peak", &_tr_by_length_pe_peak_peak);
        }
    }
  if(fTrackLA)
    {
      fChannelTree->Branch("tr_lim_angle_peak_fit", &_tr_lim_angle_peak_fit);
      fChannelTree->Branch("tr_lim_angle_peak_fit_chi2", &_tr_lim_angle_peak_fit_chi2);
      fChannelTree->Branch("tr_lim_angle_peak_fit_converged", &_tr_lim_angle_peak_fit_converged);
      fChannelTree->Branch("tr_lim_angle_peak_peak", &_tr_lim_angle_peak_peak);
      fChannelTree->Branch("tr_lim_angle_sat_rate", &_tr_lim_angle_sat_rate);
      fChannelTree->Branch("tr_lim_angle_sat_ratio_total", &_tr_lim_angle_sat_ratio_total);
      fChannelTree->Branch("tr_lim_angle_sat_ratio_peak", &_tr_lim_angle_sat_ratio_peak);

      if(fAnalysePE)
        {
          fChannelTree->Branch("tr_lim_angle_pe_peak_fit", &_tr_lim_angle_pe_peak_fit);
          fChannelTree->Branch("tr_lim_angle_pe_peak_fit_chi2", &_tr_lim_angle_pe_peak_fit_chi2);
          fChannelTree->Branch("tr_lim_angle_pe_peak_fit_converged", &_tr_lim_angle_pe_peak_fit_converged);
          fChannelTree->Branch("tr_lim_angle_pe_peak_peak", &_tr_lim_angle_pe_peak_peak);
        }
    }

  for(uint window = 0; window < fUnixWindows.size(); ++window)
    {
      if(fSaveROOTHists)
        {
          for(int ch = 0; ch < 4480; ++ch)
            {
              if(fFEBs)
                {
                  hADCPed[window][ch]      = fs->make<TH1D>(Form("hADCPed_Window%i_Channel%i", window, ch), ";ADC;Readouts", 1000, -.5, 999.5);
                  hADCPedReset[window][ch] = fs->make<TH1D>(Form("hADCPedReset_Window%i_Channel%i", window, ch), ";ADC;Readouts", 1000, -.5, 999.5);
                  hADCMaxChan[window][ch]  = fs->make<TH1D>(Form("hADCMaxChan_Window%i_Channel%i", window, ch), ";ADC;Readouts - Max Channel", 4300, -.5, 4299.5);
                }

              if(fStripHits)
                {
                  hADCSH[window][ch] = fs->make<TH1D>(Form("hADCSH_Window%i_Channel%i", window, ch), ";ADC;Strip Hits", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPESH[window][ch] = fs->make<TH1D>(Form("hPESH_Window%i_Channel%i", window, ch), ";PE;Strip Hits", 4000, 0, 200);
                }

              if(fSpacePoints)
                {
                  hADCSP[window][ch] = fs->make<TH1D>(Form("hADCSP_Window%i_Channel%i", window, ch), ";ADC;Space Points", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPESP[window][ch] = fs->make<TH1D>(Form("hPESP_Window%i_Channel%i", window, ch), ";PE;Space Points", 4000, 0, 200);
                }

              if(fTracks)
                {
                  hADCTr[window][ch]         = fs->make<TH1D>(Form("hADCTr_Window%i_Channel%i", window, ch), ";ADC;Tracks", 4300, -.5, 4299.5);
                  hADCTrByLength[window][ch] = fs->make<TH1D>(Form("hADCTrByLength_Window%i_Channel%i", window, ch), ";ADC/Path Length;Tracks", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    {
                      hPETr[window][ch]         = fs->make<TH1D>(Form("hPETr_Window%i_Channel%i", window, ch), ";PE;Tracks", 4000, 0, 200);
                      hPETrByLength[window][ch] = fs->make<TH1D>(Form("hPETrByLength_Window%i_Channel%i", window, ch), ";PE/Path Length;Tracks", 4000, 0, 200);
                    }
                }

              if(fTrackLA)
                {
                  hADCTrLA[window][ch] = fs->make<TH1D>(Form("hADCTrLA_Window%i_Channel%i", window, ch), ";ADC;Tracks", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPETrLA[window][ch] = fs->make<TH1D>(Form("hPETrLA_Window%i_Channel%i", window, ch), ";PE;Tracks", 4000, 0, 200);
                }
            }
        }
      else
        {
          for(int ch = 0; ch < 4480; ++ch)
            {
              if(fFEBs)
                {
                  hADCPed[window][ch]      = new TH1D(Form("hADCPed_Window%i_Channel%i", window, ch), ";ADC;Readouts", 1000, -.5, 999.5);
                  hADCPedReset[window][ch] = new TH1D(Form("hADCPedReset_Window%i_Channel%i", window, ch), ";ADC;Readouts", 1000, -.5, 999.5);
                  hADCMaxChan[window][ch]  = new TH1D(Form("hADCMaxChan_Window%i_Channel%i", window, ch), ";ADC;Readouts - Max Channel", 4300, -.5, 4299.5);
                }

              if(fStripHits)
                {
                  hADCSH[window][ch] = new TH1D(Form("hADCSH_Window%i_Channel%i", window, ch), ";ADC;Strip Hits", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPESH[window][ch] = new TH1D(Form("hPESH_Window%i_Channel%i", window, ch), ";PE;Strip Hits", 4000, 0, 200);
                }

              if(fSpacePoints)
                {
                  hADCSP[window][ch] = new TH1D(Form("hADCSP_Window%i_Channel%i", window, ch), ";ADC;Space Points", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPESP[window][ch] = new TH1D(Form("hPESP_Window%i_Channel%i", window, ch), ";PE;Space Points", 4000, 0, 200);
                }

              if(fTracks)
                {
                  hADCTr[window][ch]         = new TH1D(Form("hADCTr_Window%i_Channel%i", window, ch), ";ADC;Tracks", 4300, -.5, 4299.5);
                  hADCTrByLength[window][ch] = new TH1D(Form("hADCTrByLength_Window%i_Channel%i", window, ch), ";ADC/Path Length;Tracks", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    {
                      hPETr[window][ch]         = new TH1D(Form("hPETr_Window%i_Channel%i", window, ch), ";PE;Tracks", 4000, 0, 200);
                      hPETrByLength[window][ch] = new TH1D(Form("hPETrByLength_Window%i_Channel%i", window, ch), ";PE/Path Length;Tracks", 4000, 0, 200);
                    }
                }

              if(fTrackLA)
                {
                  hADCTrLA[window][ch] = new TH1D(Form("hADCTrLA_Window%i_Channel%i", window, ch), ";ADC;Tracks", 4300, -.5, 4299.5);

                  if(fAnalysePE)
                    hPETrLA[window][ch] = new TH1D(Form("hPETrLA_Window%i_Channel%i", window, ch), ";PE;Tracks", 4000, 0, 200);
                }
            }
        }
    }

  for(uint i = 0; i < fUnixWindows.size(); ++i)
    {
      for(uint ii = i + 1; ii < fUnixWindows.size(); ++ii)
        {
          if(fUnixWindows[i].second > fUnixWindows[ii].first)
            {
              std::cout << "Unix windows overlap!" << std::endl;
              throw std::exception();
            }
        }
    }
}

void sbnd::crt::ADRIFT::analyze(art::Event const& e)
{
  // Get Raw TS
  art::Handle<artdaq::detail::RawEventHeader> DAQHeaderHandle;
  e.getByLabel(fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel, DAQHeaderHandle);
  if(!DAQHeaderHandle.isValid()){
    std::cout << "RawEventHeader product " << fDAQHeaderModuleLabel << " - " << fDAQHeaderInstanceLabel << " not found..." << std::endl;
    throw std::exception();
  }
  artdaq::RawEvent rawHeaderEvent = artdaq::RawEvent(*DAQHeaderHandle);
  uint64_t raw_ts = rawHeaderEvent.timestamp() - fRawTSCorrection;

  uint window = 0;
  bool found  = false;

  for(uint i = 0; i < fUnixWindows.size(); ++i)
    {
      if(raw_ts > fUnixWindows[i].first && raw_ts < fUnixWindows[i].second)
        {
          window = i;
          found  = true;
          break;
        }
    }

  if(!found)
    return;

  if(fNEvents == 0)
    MakeSaveDirectories(e);

  if(fFEBs)
    AnalyseFEBDatas(e, window);

  if(fStripHits)
    AnalyseStripHits(e, window);

  if(fSpacePoints)
    AnalyseSpacePoints(e, window);

  if(fTracks || fTrackLA)
    AnalyseTracks(e, window);

  ++fNEvents;
}

void sbnd::crt::ADRIFT::MakeSaveDirectories(art::Event const &e)
{
  for(uint window = 0; window < fUnixWindows.size(); ++window)
    {
      std::string directory_structure = fUnixWindows.size() == 1 ? Form("run%i", e.run()) : Form("run%i/window%i", e.run(), window);

      if((fSaveAllFits || fSaveBadFits))
        {
          gSystem->Exec(Form("mkdir -p %s", fTopSaveDirectory.c_str()));

          if(fFEBs)
            {
              fPedestalSaveDirectory[window] = Form("%s/%s/pedestals", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fPedestalSaveDirectory[window].c_str()));

              fPedestalResetSaveDirectory[window] = Form("%s/%s/pedestals_reset", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fPedestalResetSaveDirectory[window].c_str()));

              fBadPedestalSaveDirectory[window] = Form("%s/%s/pedestals/bad", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fBadPedestalSaveDirectory[window].c_str()));

              fBadPedestalResetSaveDirectory[window] = Form("%s/%s/pedestals_reset/bad", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fBadPedestalResetSaveDirectory[window].c_str()));
            }

          if(fStripHits || fSpacePoints || fTracks || fTrackLA)
            {
              fPeakSaveDirectory[window] = Form("%s/%s/peaks", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fPeakSaveDirectory[window].c_str()));

              fBadPeakSaveDirectory[window] = Form("%s/%s/peaks/bad", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fBadPeakSaveDirectory[window].c_str()));
            }
        }

      if(fSaveSubset)
        {
          gSystem->Exec(Form("mkdir -p %s", fTopSaveDirectory.c_str()));

          if(fFEBs)
            {
              fPedestalSubsetSaveDirectory[window] = Form("%s/%s/subset/pedestals", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fPedestalSubsetSaveDirectory[window].c_str()));

              fPedestalResetSubsetSaveDirectory[window] = Form("%s/%s/subset/pedestals_reset", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fPedestalResetSubsetSaveDirectory[window].c_str()));
            }

          if(fStripHits)
            {
              fStripHitSubsetSaveDirectory[window] = Form("%s/%s/subset/striphits", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fStripHitSubsetSaveDirectory[window].c_str()));

              if(fAnalysePE)
                {
                  fStripHitPESubsetSaveDirectory[window] = Form("%s/%s/subset/striphits_pe", fTopSaveDirectory.c_str(), directory_structure.c_str());
                  gSystem->Exec(Form("mkdir -p %s", fStripHitPESubsetSaveDirectory[window].c_str()));
                }
            }

          if(fSpacePoints)
            {
              fSpacePointSubsetSaveDirectory[window] = Form("%s/%s/subset/spacepoints", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fSpacePointSubsetSaveDirectory[window].c_str()));

              if(fAnalysePE)
                {
                  fSpacePointPESubsetSaveDirectory[window] = Form("%s/%s/subset/spacepoints_pe", fTopSaveDirectory.c_str(), directory_structure.c_str());
                  gSystem->Exec(Form("mkdir -p %s", fSpacePointPESubsetSaveDirectory[window].c_str()));
                }
            }

          if(fTracks)
            {
              fTrackSubsetSaveDirectory[window] = Form("%s/%s/subset/tracks", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fTrackSubsetSaveDirectory[window].c_str()));

              fTrackByLengthSubsetSaveDirectory[window] = Form("%s/%s/subset/tracks_by_length", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fTrackByLengthSubsetSaveDirectory[window].c_str()));

              if(fAnalysePE)
                {
                  fTrackPESubsetSaveDirectory[window] = Form("%s/%s/subset/tracks_pe", fTopSaveDirectory.c_str(), directory_structure.c_str());
                  gSystem->Exec(Form("mkdir -p %s", fTrackPESubsetSaveDirectory[window].c_str()));

                  fTrackByLengthPESubsetSaveDirectory[window] = Form("%s/%s/subset/tracks_by_length_pe", fTopSaveDirectory.c_str(), directory_structure.c_str());
                  gSystem->Exec(Form("mkdir -p %s", fTrackByLengthPESubsetSaveDirectory[window].c_str()));
                }
            }

          if(fTrackLA)
            {
              fTrackLASubsetSaveDirectory[window] = Form("%s/%s/subset/tracks_limited_angle", fTopSaveDirectory.c_str(), directory_structure.c_str());
              gSystem->Exec(Form("mkdir -p %s", fTrackLASubsetSaveDirectory[window].c_str()));

              if(fAnalysePE)
                {
                  fTrackLAPESubsetSaveDirectory[window] = Form("%s/%s/subset/tracks_limited_angle_pe", fTopSaveDirectory.c_str(), directory_structure.c_str());
                  gSystem->Exec(Form("mkdir -p %s", fTrackLAPESubsetSaveDirectory[window].c_str()));
                }
            }
        }
    }

  if((fSaveAllFits || fSaveBadFits || fSaveSubset))
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
}

void sbnd::crt::ADRIFT::AnalyseFEBDatas(art::Event const &e, const int window)
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
      const uint m5    = feb_data->Mac5();
      const uint flags = feb_data->Flags();

      std::set<int> mask_channels;

      for(const art::Ptr<CRTStripHit> &strip_hit : strip_hits)
        {
          mask_channels.insert(strip_hit->Channel());
          mask_channels.insert(strip_hit->Channel() + 1);

          if(strip_hit->Channel() % 2)
            std::cout << "ODD Strip Hit Channel Number" << std::endl;
        }

      int max_adc = -1, max_chan = -1;

      for(int i = 0; i < 32; ++i)
        {
          if(feb_data->ADC(i) > max_adc)
            {
              max_adc  = feb_data->ADC(i);
              max_chan = i;
            }

          const int ch = m5 * 32 + i;

          if(!mask_channels.count(ch))
            hADCPed[window][ch]->Fill(feb_data->ADC(i));

          if(flags != 1 && flags != 3)
            hADCPedReset[window][ch]->Fill(feb_data->ADC(i));
        }

      if((flags == 1 || flags == 3) && max_chan != -1)
        hADCMaxChan[window][m5 * 32 + max_chan]->Fill(feb_data->ADC(max_chan));
    }
}

void sbnd::crt::ADRIFT::AnalyseStripHits(art::Event const &e, const int window)
{
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

      hADCSH[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
      hADCSH[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);

      if(fAnalysePE)
        {
          if((strip_hit->ADC1() + sipm1.pedestal) < 4089)
            hPESH[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() * sipm1.gain);
          if((strip_hit->ADC2() + sipm2.pedestal) < 4089)
            hPESH[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() * sipm2.gain);
        }
    }
}

void sbnd::crt::ADRIFT::AnalyseSpacePoints(art::Event const &e, const int window)
{
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
          
          hADCSP[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
          hADCSP[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);

          if(fAnalysePE)
            {
              if((strip_hit->ADC1() + sipm1.pedestal) < 4089)
                hPESP[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() * sipm1.gain);
              if((strip_hit->ADC2() + sipm2.pedestal) < 4089)
                hPESP[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() * sipm2.gain);
            }
        }
    }
}

void sbnd::crt::ADRIFT::AnalyseTracks(art::Event const &e, const int window)
{
  // Get CRTTracks
  art::Handle<std::vector<CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  if(!CRTTrackHandle.isValid()){
    std::cout << "CRTTrack product " << fCRTTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTTrack>> CRTTrackVec;
  art::fill_ptr_vector(CRTTrackVec, CRTTrackHandle);

  // Get CRTSpacePoints
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get CRTClusters
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get CRTTrack to CRTSpacePoint Assns
  art::FindManyP<CRTSpacePoint> tracksToSpacePoints(CRTTrackHandle, e, fCRTTrackModuleLabel);

  // Get CRTSpacePoint to CRTCluster Assns
  art::FindOneP<CRTCluster> spacepointsToClusters(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  // Get CRTCluster to CRTStripHit Assns
  art::FindManyP<CRTStripHit> clustersToStripHits(CRTClusterHandle, e, fCRTClusterModuleLabel);

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
          const double path_length = 1 / TMath::Cos(normal.Angle(tr_dir));

          for(const art::Ptr<CRTStripHit> &strip_hit : strip_hits)
            {
              if(strip_hit->Channel() % 2)
                std::cout << "ODD Strip Hit Channel Number" << std::endl;
      
              CRTSiPMGeo sipm1 = fCRTGeoAlg.GetSiPM(strip_hit->Channel());
              CRTSiPMGeo sipm2 = fCRTGeoAlg.GetSiPM(strip_hit->Channel() + 1);
          
              if(fTracks)
                {
                  hADCTr[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
                  hADCTr[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);

                  if(fAnalysePE)
                    {
                      if((strip_hit->ADC1() + sipm1.pedestal) < 4089)
                        hPETr[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() * sipm1.gain);
                      if((strip_hit->ADC2() + sipm2.pedestal) < 4089)
                        hPETr[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() * sipm2.gain);
                    }

                  if((strip_hit->ADC1() + sipm1.pedestal) < 4089)
                    {
                      hADCTrByLength[window][strip_hit->Channel()]->Fill((strip_hit->ADC1() + sipm1.pedestal) / path_length);

                      if(fAnalysePE && (strip_hit->ADC1() + sipm1.pedestal) < 4089)
                        hPETrByLength[window][strip_hit->Channel()]->Fill((strip_hit->ADC1() * sipm1.gain) / path_length);
                    }
                  if((strip_hit->ADC2() + sipm2.pedestal) < 4089)
                    {
                      hADCTrByLength[window][strip_hit->Channel() + 1]->Fill((strip_hit->ADC2() + sipm2.pedestal) / path_length);

                      if(fAnalysePE && (strip_hit->ADC2() + sipm2.pedestal) < 4089)
                        hPETrByLength[window][strip_hit->Channel() + 1]->Fill((strip_hit->ADC2() * sipm2.gain) / path_length);
                    }
                }

              if(fTrackLA && angle < fTrackAngleLimit)
                {
                  hADCTrLA[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() + sipm1.pedestal);
                  hADCTrLA[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() + sipm2.pedestal);

                  if(fAnalysePE)
                    {
                      if((strip_hit->ADC1() + sipm1.pedestal) < 4089)
                        hPETrLA[window][strip_hit->Channel()]->Fill(strip_hit->ADC1() * sipm1.gain);
                      if((strip_hit->ADC2() + sipm2.pedestal) < 4089)
                        hPETrLA[window][strip_hit->Channel() + 1]->Fill(strip_hit->ADC2() * sipm2.gain);
                    }
                }
            }
        }
    }
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

          _channel        = ch;
          _gdml_id        = gdml_i;
          _mac5           = mac5;
          _raw_channel    = invert ? 31 - ch_i : ch_i;
          _tagger         = fCRTGeoAlg.ChannelToTaggerEnum(ch);
          _channel_status = fCRTGeoAlg.GetSiPM(ch).status;
          _area           = fCRTGeoAlg.StripArea(ch);
          _y_average      = fCRTGeoAlg.StripAverageY(ch);
          _ped_calib      = fCRTGeoAlg.GetSiPM(ch).pedestal;
          _gain_calib     = fCRTGeoAlg.GetSiPM(ch).gain;
          _horizontal     = _tagger == kBottomTagger || _tagger == kTopLowTagger || _tagger == kTopHighTagger;

          for(uint window = 0; window < fUnixWindows.size(); ++window)
            ProcessEntry(ch, window);
        }
    }
}

void sbnd::crt::ADRIFT::ProcessEntry(const int ch, const int window)
{
  if(fSaveSubset && (ch > 1471 && ch < 1728))
    {
      if(fFEBs)
        {
          SaveHist(hADCPed[window][ch], fPedestalSubsetSaveDirectory[window], Form("pedestal_channel_%i", ch), 2, _channel_status);
          SaveHist(hADCPedReset[window][ch], fPedestalResetSubsetSaveDirectory[window], Form("pedestal_reset_channel_%i", ch), 2, _channel_status);
        }

      if(fStripHits)
        {
          SaveHist(hADCSH[window][ch], fStripHitSubsetSaveDirectory[window], Form("strip_hit_channel_%i", ch), 20, _channel_status);

          if(fAnalysePE)
            SaveHist(hPESH[window][ch], fStripHitPESubsetSaveDirectory[window], Form("strip_hit_pe_channel_%i", ch), 20, _channel_status);
        }

      if(fSpacePoints)
        {
          SaveHist(hADCSP[window][ch], fSpacePointSubsetSaveDirectory[window], Form("space_point_channel_%i", ch), 20, _channel_status);
          if(fAnalysePE)
            SaveHist(hPESP[window][ch], fSpacePointPESubsetSaveDirectory[window], Form("space_point_pe_channel_%i", ch), 20, _channel_status);
        }

      if(fTracks)
        {
          SaveHist(hADCTr[window][ch], fTrackSubsetSaveDirectory[window], Form("track_channel_%i", ch), 20, _channel_status);
          SaveHist(hADCTrByLength[window][ch], fTrackByLengthSubsetSaveDirectory[window], Form("track_by_length_channel_%i", ch), 20, _channel_status);
          if(fAnalysePE)
            {
              SaveHist(hPETr[window][ch], fTrackPESubsetSaveDirectory[window], Form("track_pe_channel_%i", ch), 20, _channel_status);
              SaveHist(hPETrByLength[window][ch], fTrackByLengthPESubsetSaveDirectory[window], Form("track_by_length_pe_channel_%i", ch), 20, _channel_status);
            }
        }

      if(fTrackLA)
        {
          SaveHist(hADCTrLA[window][ch], fTrackLASubsetSaveDirectory[window], Form("track_limited_angle_channel_%i", ch), 20, _channel_status);
          if(fAnalysePE)
            SaveHist(hPETrLA[window][ch], fTrackLAPESubsetSaveDirectory[window], Form("track_limited_angle_pe_channel_%i", ch), 20, _channel_status);
        }

    }

  if(fFEBs)
    {
      PedestalFit(hADCPed[window][ch], _ped_fit, _ped_fit_std, _ped_fit_chi2, _ped_fit_converged, _channel_status, window);
      PedestalPeak(hADCPed[window][ch], _ped_peak);

      PedestalFit(hADCPedReset[window][ch], _ped_reset_fit, _ped_reset_fit_std, _ped_reset_fit_chi2, _ped_reset_fit_converged, _channel_status, window);
      PedestalPeak(hADCPedReset[window][ch], _ped_reset_peak);

      Rate(hADCMaxChan[window][ch], _raw_max_chan_rate);
    }

  if(fStripHits)
    {
      Rate(hADCSH[window][ch], _sh_rate);

      PeakPeak(hADCSH[window][ch], _ped_calib, _sh_peak_peak);
      PeakFit(hADCSH[window][ch], _sh_peak_peak, _ped_calib, _sh_peak_fit, _sh_peak_fit_chi2, _sh_peak_fit_converged, _channel_status, window);
      _sh_ped_to_peak_fit       = _sh_peak_fit - _ped_fit;
      _sh_ped_reset_to_peak_fit = _sh_peak_fit - _ped_reset_fit;

      if(fAnalysePE)
        {
          PeakPeak(hPESH[window][ch], _ped_calib, _sh_pe_peak_peak);
          PeakFit(hPESH[window][ch], _sh_pe_peak_peak, _ped_calib, _sh_pe_peak_fit, _sh_pe_peak_fit_chi2, _sh_pe_peak_fit_converged, _channel_status, window);
        }

      const double sh_sat = Saturation(hADCSH[window][ch]);
      _sh_sat_rate        = sh_sat / (fNEvents * fPullWindow);
      _sh_sat_ratio_total = sh_sat / hADCSH[window][ch]->GetEntries();
      _sh_sat_ratio_peak  = sh_sat / hADCSH[window][ch]->GetBinContent(hADCSH[window][ch]->FindBin(_sh_peak_peak));
    }

  if(fSpacePoints)
    {
      Rate(hADCSP[window][ch], _sp_rate);

      PeakPeak(hADCSP[window][ch], _ped_calib, _sp_peak_peak);
      PeakFit(hADCSP[window][ch], _sp_peak_peak, _ped_calib, _sp_peak_fit, _sp_peak_fit_chi2, _sp_peak_fit_converged, _channel_status, window);

      if(fAnalysePE)
        {
          PeakPeak(hPESP[window][ch], _ped_calib, _sp_pe_peak_peak);
          PeakFit(hPESP[window][ch], _sp_pe_peak_peak, _ped_calib, _sp_pe_peak_fit, _sp_pe_peak_fit_chi2, _sp_pe_peak_fit_converged, _channel_status, window);
        }

      const double sp_sat = Saturation(hADCSP[window][ch]);
      _sp_sat_rate        = sp_sat / (fNEvents * fPullWindow);
      _sp_sat_ratio_total = sp_sat / hADCSP[window][ch]->GetEntries();
      _sp_sat_ratio_peak  = sp_sat / hADCSP[window][ch]->GetBinContent(hADCSP[window][ch]->FindBin(_sp_peak_peak));
    }

  if(fTracks && _tagger != kBottomTagger)
    {
      Rate(hADCTr[window][ch], _tr_rate);

      PeakPeak(hADCTr[window][ch], _ped_calib, _tr_peak_peak);
      PeakFit(hADCTr[window][ch], _tr_peak_peak, _ped_calib, _tr_peak_fit, _tr_peak_fit_chi2, _tr_peak_fit_converged, _channel_status, window);

      PeakPeak(hADCTrByLength[window][ch], _ped_calib, _tr_by_length_peak_peak);
      PeakFit(hADCTrByLength[window][ch], _tr_by_length_peak_peak, _ped_calib, _tr_by_length_peak_fit, _tr_by_length_peak_fit_chi2, _tr_by_length_peak_fit_converged, _channel_status, window);

      if(fAnalysePE)
        {
          PeakPeak(hPETr[window][ch], _ped_calib, _tr_pe_peak_peak);
          PeakFit(hPETr[window][ch], _tr_pe_peak_peak, _ped_calib, _tr_pe_peak_fit, _tr_pe_peak_fit_chi2, _tr_pe_peak_fit_converged, _channel_status, window);

          PeakPeak(hPETrByLength[window][ch], _ped_calib, _tr_by_length_pe_peak_peak);
          PeakFit(hPETrByLength[window][ch], _tr_by_length_pe_peak_peak, _ped_calib, _tr_by_length_pe_peak_fit, _tr_by_length_pe_peak_fit_chi2, _tr_by_length_pe_peak_fit_converged, _channel_status, window);
        }

      const double tr_sat = Saturation(hADCTr[window][ch]);
      _tr_sat_rate        = tr_sat / (fNEvents * fPullWindow);
      _tr_sat_ratio_total = tr_sat / hADCTr[window][ch]->GetEntries();
      _tr_sat_ratio_peak  = tr_sat / hADCTr[window][ch]->GetBinContent(hADCTr[window][ch]->FindBin(_tr_peak_peak));
    }

  if(fTrackLA && _tagger != kBottomTagger)
    {
      PeakPeak(hADCTrLA[window][ch], _ped_calib, _tr_lim_angle_peak_peak);
      PeakFit(hADCTrLA[window][ch], _tr_lim_angle_peak_peak, _ped_calib, _tr_lim_angle_peak_fit, _tr_lim_angle_peak_fit_chi2, _tr_lim_angle_peak_fit_converged, _channel_status, window);

      if(fAnalysePE)
        {
          PeakPeak(hPETrLA[window][ch], _ped_calib, _tr_lim_angle_pe_peak_peak);
          PeakFit(hPETrLA[window][ch], _tr_lim_angle_pe_peak_peak, _ped_calib, _tr_lim_angle_pe_peak_fit, _tr_lim_angle_pe_peak_fit_chi2, _tr_lim_angle_pe_peak_fit_converged, _channel_status, window);
        }

      const double tr_lim_angle_sat = Saturation(hADCTrLA[window][ch]);
      _tr_lim_angle_sat_rate        = tr_lim_angle_sat / (fNEvents * fPullWindow);
      _tr_lim_angle_sat_ratio_total = tr_lim_angle_sat / hADCTrLA[window][ch]->GetEntries();
      _tr_lim_angle_sat_ratio_peak  = tr_lim_angle_sat / hADCTrLA[window][ch]->GetBinContent(hADCTrLA[window][ch]->FindBin(_tr_lim_angle_peak_peak));
    }

  fChannelTree->Fill();
}

void sbnd::crt::ADRIFT::PedestalFit(TH1D* hADC, double &fit, double &std, double &chi2, bool &converged,
                                    bool badChannel, const int window)
{
  TF1 *gaus = new TF1("gaus", "gaus", 0, 500);
  const TString name = hADC->GetName();
  TString ch_name = name;
  TString type    = "";

  const int window_log = window == 0 ? 0 : std::floor(std::log10(window));

  if(name.Contains("Reset"))
    {
      ch_name.Remove(0, 27 + window_log + 1);
      type = "reset";
    }
  else
    ch_name.Remove(0,22 + window_log + 1);

  int ch = std::stoi(ch_name.Data());

  TH1D* hADC2 = (TH1D*) hADC->Clone(name + "_for_fit");
  hADC2->Rebin(5);

  TFitResultPtr fit_result = hADC2->Fit(gaus, "QRS");
  converged = !(bool)(int(fit_result));

  if(!converged && !badChannel)
    std::cout << "Pedestal fit has not converged - " << hADC->GetName() << std::endl;

  fit  = gaus->GetParameter("Mean");
  std  = gaus->GetParameter("Sigma");
  chi2 = gaus->GetChisquare() / gaus->GetNDF();

  if(fSaveAllFits || (fSaveBadFits && !converged) || (fSaveSubset && ch > 1471 && ch < 1728))
    {
      TCanvas *c = new TCanvas(Form("c%s", name.Data()), Form("c%s", name.Data()));
      c->cd();

      hADC2->SetLineColor(kOrange+2);
      hADC2->SetLineWidth(2);
      hADC2->Draw("histe");
      gaus->SetLineColor(kSpring-6);
      gaus->Draw("same");

      if(badChannel)
        {
          TText *t = new TText(.5, .75, "Bad Channel");
          t->SetTextSize(0.1);
          t->SetTextColor(kRed);
          t->Draw();
        }

      if(converged)
        {
          if(type.Contains("reset"))
            {
              c->SaveAs(Form("%s/pedestal_reset_fit_channel_%s.png", fPedestalResetSaveDirectory[window].c_str(), ch_name.Data()));
              c->SaveAs(Form("%s/pedestal_reset_fit_channel_%s.pdf", fPedestalResetSaveDirectory[window].c_str(), ch_name.Data()));
            }
          else
            {
              c->SaveAs(Form("%s/pedestal_fit_channel_%s.png", fPedestalSaveDirectory[window].c_str(), ch_name.Data()));
              c->SaveAs(Form("%s/pedestal_fit_channel_%s.pdf", fPedestalSaveDirectory[window].c_str(), ch_name.Data()));
            }
        }
      else
        {
          if(type.Contains("reset"))
            {
              c->SaveAs(Form("%s/pedestal_reset_fit_channel_%s.png", fBadPedestalResetSaveDirectory[window].c_str(), ch_name.Data()));
              c->SaveAs(Form("%s/pedestal_reset_fit_channel_%s.pdf", fBadPedestalResetSaveDirectory[window].c_str(), ch_name.Data()));
            }
          else
            {
              c->SaveAs(Form("%s/pedestal_fit_channel_%s.png", fBadPedestalSaveDirectory[window].c_str(), ch_name.Data()));
              c->SaveAs(Form("%s/pedestal_fit_channel_%s.pdf", fBadPedestalSaveDirectory[window].c_str(), ch_name.Data()));
            }
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

void sbnd::crt::ADRIFT::PeakPeak(TH1D* hist, const double &ped, double &peak)
{
  int bin = 0;
  double max = std::numeric_limits<double>::lowest();

  const TString name = hist->GetName();

  const int start = name.Contains("PE") ? 0 : (int)ped + 20;

  for(int i = start; i < 4000; ++i)
    {
      if(hist->GetBinContent(i) > max)
        {
          max = hist->GetBinContent(i);
          bin = i;
        }
    }

  peak = hist->GetBinCenter(bin);
}

void sbnd::crt::ADRIFT::PeakFit(TH1D* hist, const double &peak, const double &ped, double &fit,
                                double &chi2, bool &converged, bool badChannel, const int window)
{
  const TString name = hist->GetName();
  const bool pe      = name.Contains("PE");
  TString ch_name    = name;
  TString type       = "";

  const double start = pe ? 0 : ped + 20;
  const double end   = pe ? 200 : 4000;
  const double width = pe ? 0.25 : 10;
  const double sigma = pe ? 1 : 50;

  TF1 *langau = new TF1("langau", LanGau, start, end, 4);
  double params[4] = { width, peak, hist->GetEntries(), sigma };
  langau->SetParameters(params);

  const int window_log = window == 0 ? 0 : std::floor(std::log10(window));

  if(name.Contains("SH"))
    {
      const int end = pe ? 20 + window_log + 1 : 21 + window_log + 1;
      ch_name.Remove(0, end);
      type = "strip_hit";
    }
  else if(name.Contains("SP"))
    {
      const int end = pe ? 20 + window_log + 1 : 21 + window_log + 1;
      ch_name.Remove(0, end);
      type = "space_point";
    }
  else if(name.Contains("TrLA"))
    {
      const int end = pe ? 22 + window_log + 1 : 23 + window_log + 1;
      ch_name.Remove(0, end);
      type = "track_limited_angle";
    }
  else if(name.Contains("TrByLength"))
    {
      const int end = pe ? 28 + window_log + 1 : 29 + window_log + 1;
      ch_name.Remove(0, end);
      type = "track_by_length";
    }
  else if(name.Contains("Tr"))
    {
      const int end = pe ? 20 + window_log + 1 : 21  + window_log + 1;
      ch_name.Remove(0, end);
      type = "track";
    }
  else
    std::cout << "Cannot identify histogram type (" << name << ")" << std::endl;

  int ch = std::stoi(ch_name.Data());

  TH1D* hist2 = (TH1D*) hist->Clone(name + "_for_fit");
  hist2->Rebin(20);

  TFitResultPtr fit_result = hist2->Fit(langau, "QRS");
  converged = !(bool)(int(fit_result));

  if(!converged && !badChannel)
    std::cout << "Peak fit has not converged - " << hist->GetName() << std::endl;

  fit  = langau->GetParameter(1);
  chi2 = langau->GetChisquare() / langau->GetNDF();

  if(fSaveAllFits || (fSaveBadFits && !converged) || (fSaveSubset && ch > 1471 && ch < 1728))
    {
      TCanvas *c = new TCanvas(Form("c%s", name.Data()), Form("c%s", name.Data()));
      c->cd();

      hist2->SetLineColor(kOrange+2);
      hist2->SetLineWidth(2);
      hist2->Draw("histe");
      langau->SetLineColor(kSpring-6);
      langau->Draw("same");

      if(badChannel)
        {
          TText *t = new TText(.5, .75, "Bad Channel");
          t->SetTextSize(0.1);
          t->SetTextColor(kRed);
          t->Draw();
        }

      if(pe)
        {
          if(converged)
            {
              c->SaveAs(Form("%s/peak_fit_%s_pe_channel_%s.png", fPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
              c->SaveAs(Form("%s/peak_fit_%s_pe_channel_%s.pdf", fPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
            }
          else
            {
              c->SaveAs(Form("%s/peak_fit_%s_pe_channel_%s.png", fBadPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
              c->SaveAs(Form("%s/peak_fit_%s_pe_channel_%s.pdf", fBadPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
            }
        }
      else
        {
          if(converged)
            {
              c->SaveAs(Form("%s/peak_fit_%s_channel_%s.png", fPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
              c->SaveAs(Form("%s/peak_fit_%s_channel_%s.pdf", fPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
            }
          else
            {
              c->SaveAs(Form("%s/peak_fit_%s_channel_%s.png", fBadPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
              c->SaveAs(Form("%s/peak_fit_%s_channel_%s.pdf", fBadPeakSaveDirectory[window].c_str(), type.Data(), ch_name.Data()));
            }
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
  _channel        = std::numeric_limits<int>::lowest();
  _gdml_id        = std::numeric_limits<int>::lowest();
  _mac5           = std::numeric_limits<int>::lowest();
  _raw_channel    = std::numeric_limits<int>::lowest();
  _tagger         = std::numeric_limits<int>::lowest();
  _channel_status = std::numeric_limits<int>::lowest();

  _area                          = std::numeric_limits<double>::lowest();
  _y_average                     = std::numeric_limits<double>::lowest();
  _ped_calib                     = std::numeric_limits<double>::lowest();
  _gain_calib                    = std::numeric_limits<double>::lowest();
  _ped_fit                       = std::numeric_limits<double>::lowest();
  _ped_fit_std                   = std::numeric_limits<double>::lowest();
  _ped_fit_chi2                  = std::numeric_limits<double>::lowest();
  _ped_peak                      = std::numeric_limits<double>::lowest();
  _ped_reset_fit                 = std::numeric_limits<double>::lowest();
  _ped_reset_fit_std             = std::numeric_limits<double>::lowest();
  _ped_reset_fit_chi2            = std::numeric_limits<double>::lowest();
  _ped_reset_peak                = std::numeric_limits<double>::lowest();
  _sh_rate                       = std::numeric_limits<double>::lowest();
  _sp_rate                       = std::numeric_limits<double>::lowest();
  _tr_rate                       = std::numeric_limits<double>::lowest();
  _sh_peak_fit                   = std::numeric_limits<double>::lowest();
  _sh_peak_fit_chi2              = std::numeric_limits<double>::lowest();
  _sh_peak_peak                  = std::numeric_limits<double>::lowest();
  _sh_sat_rate                   = std::numeric_limits<double>::lowest();
  _sh_sat_ratio_total            = std::numeric_limits<double>::lowest();
  _sh_sat_ratio_peak             = std::numeric_limits<double>::lowest();
  _sh_ped_to_peak_fit            = std::numeric_limits<double>::lowest();
  _sh_ped_reset_to_peak_fit      = std::numeric_limits<double>::lowest();
  _sp_peak_fit                   = std::numeric_limits<double>::lowest();
  _sp_peak_fit_chi2              = std::numeric_limits<double>::lowest();
  _sp_peak_peak                  = std::numeric_limits<double>::lowest();
  _sp_sat_rate                   = std::numeric_limits<double>::lowest();
  _sp_sat_ratio_total            = std::numeric_limits<double>::lowest();
  _sp_sat_ratio_peak             = std::numeric_limits<double>::lowest();
  _tr_peak_fit                   = std::numeric_limits<double>::lowest();
  _tr_peak_fit_chi2              = std::numeric_limits<double>::lowest();
  _tr_peak_peak                  = std::numeric_limits<double>::lowest();
  _tr_sat_rate                   = std::numeric_limits<double>::lowest();
  _tr_sat_ratio_total            = std::numeric_limits<double>::lowest();
  _tr_sat_ratio_peak             = std::numeric_limits<double>::lowest();
  _tr_by_length_peak_fit         = std::numeric_limits<double>::lowest();
  _tr_by_length_peak_fit_chi2    = std::numeric_limits<double>::lowest();
  _tr_by_length_peak_peak        = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_fit         = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_fit_chi2    = std::numeric_limits<double>::lowest();
  _tr_lim_angle_peak_peak        = std::numeric_limits<double>::lowest();
  _tr_lim_angle_sat_rate         = std::numeric_limits<double>::lowest();
  _tr_lim_angle_sat_ratio_total  = std::numeric_limits<double>::lowest();
  _tr_lim_angle_sat_ratio_peak   = std::numeric_limits<double>::lowest();
  _sh_pe_peak_fit                = std::numeric_limits<double>::lowest();
  _sh_pe_peak_fit_chi2           = std::numeric_limits<double>::lowest();
  _sh_pe_peak_peak               = std::numeric_limits<double>::lowest();
  _sp_pe_peak_fit                = std::numeric_limits<double>::lowest();
  _sp_pe_peak_fit_chi2           = std::numeric_limits<double>::lowest();
  _sp_pe_peak_peak               = std::numeric_limits<double>::lowest();
  _tr_pe_peak_fit                = std::numeric_limits<double>::lowest();
  _tr_pe_peak_fit_chi2           = std::numeric_limits<double>::lowest();
  _tr_pe_peak_peak               = std::numeric_limits<double>::lowest();
  _tr_by_length_pe_peak_fit      = std::numeric_limits<double>::lowest();
  _tr_by_length_pe_peak_fit_chi2 = std::numeric_limits<double>::lowest();
  _tr_by_length_pe_peak_peak     = std::numeric_limits<double>::lowest();
  _tr_lim_angle_pe_peak_fit      = std::numeric_limits<double>::lowest();
  _tr_lim_angle_pe_peak_fit_chi2 = std::numeric_limits<double>::lowest();
  _tr_lim_angle_pe_peak_peak     = std::numeric_limits<double>::lowest();

  _horizontal                         = false;
  _ped_fit_converged                  = false;
  _ped_reset_fit_converged            = false;
  _sh_peak_fit_converged              = false;
  _sp_peak_fit_converged              = false;
  _tr_peak_fit_converged              = false;
  _tr_by_length_peak_fit_converged    = false;
  _tr_lim_angle_peak_fit_converged    = false;
  _sh_pe_peak_fit_converged           = false;
  _sp_pe_peak_fit_converged           = false;
  _tr_pe_peak_fit_converged           = false;
  _tr_by_length_pe_peak_fit_converged = false;
  _tr_lim_angle_pe_peak_fit_converged = false;
}

void sbnd::crt::ADRIFT::SaveHist(TH1D *hist, std::string &saveDir, std::string saveName, int rebin, bool badChannel)
{
  TCanvas *c = new TCanvas(Form("c%s", saveName.c_str()), Form("c%s", saveName.c_str()));
  c->cd();

  const TString name = hist->GetName();
  TH1D* hist2 = (TH1D*) hist->Clone(name + "_for_save");

  hist2->SetLineColor(kOrange+2);
  hist2->SetLineWidth(2);
  hist2->Rebin(rebin);
  hist2->Draw("histe");

  if(badChannel)
    {
      TText *t = new TText(.5, .75, "Bad Channel");
      t->SetTextSize(0.1);
      t->SetTextColor(kRed);
      t->Draw();
    }

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
