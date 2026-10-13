////////////////////////////////////////////////////////////////////////
// Class:       XARAPUCATopHatPlot
// Plugin Type: analyzer (Unknown Unknown)
// File:        XARAPUCATopHatPlot_module.cc
// Author:      Alicia Vázquez-Ramos (aliciavr@ugr.es)
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

#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "art/Framework/Core/FileBlock.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TLinearFitter.h"
#include "TParameter.h"
#include "TFile.h"

#include <limits>
#include <algorithm>

namespace thp {
  class XARAPUCATopHatPlot;
  struct RisingEdgeFit;
  struct MaxADCInfo;
};

struct thp::RisingEdgeFit {
  float slope = 0.0f;                  // Slope (m)
  float y_intercept = 0.0f;            // Y-intercept (b) [ADCs]
  float x_intercept = 0.0f;            // X-intercept (-b/m)
  float flash_time = 0.0f;             // Flash time [us]
  float chi2 = 0.0f;                   // Chi-square value
  float reduced_chi2 = 0.0f;           // Reduced chi-square value
  int ndf = -1;                        // Number of degrees of freedom
  std::vector<float> fit_points_x;     // X values of the fit points
  std::vector<float> fit_points_y;     // Y values of the fit points
};

struct thp::MaxADCInfo {
  int max_ADC = 0;                     // Maximum ADC value
  size_t max_ADC_pos = 0;              // Position of the maximum ADC value
  float max_ADC_time = -1.0;           // Time of the maximum ADC value [us]
};

class thp::XARAPUCATopHatPlot : public art::EDAnalyzer {
public:
  explicit XARAPUCATopHatPlot(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  XARAPUCATopHatPlot(XARAPUCATopHatPlot const&) = delete;
  XARAPUCATopHatPlot(XARAPUCATopHatPlot&&) = delete;
  XARAPUCATopHatPlot& operator=(XARAPUCATopHatPlot const&) = delete;
  XARAPUCATopHatPlot& operator=(XARAPUCATopHatPlot&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void reset_event_variables();
  void get_rising_edge(std::vector<int> const& summed_wvfm, float time_stamp, MaxADCInfo& max_ADC_info, RisingEdgeFit& rising_edge_fit);

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  virtual void respondToOpenInputFile(art::FileBlock const &fb) override;


private:
  // TTrees
  TTree* fTree_events; // Event information tree.
  TTree* fTree_thp;    // Top-hat plot tree.
  TTree* fTree_config; // Configuration tree.

  // Tags
  art::InputTag fTDCDecodeLabel;
  art::InputTag fWaveformsModuleLabel;
  art::InputTag fTimingReferenceModuleLabel;

  // Verbose and debug options
  bool fVerbose;
  bool fSaveJittering;
  bool fDebugWaveforms;
  bool fDebugTDC;
  bool fDebugRisingEdgeFit;
  bool fDebugJittering;

  // Useful constants
  constexpr static float NS_TO_US = 1E-3;             // Conversion factor.
  constexpr static float US_TO_NS = 1E+3;             // Conversion factor.
  constexpr static int DEFAULT_INT = -999999;         // Default integer value.
  // Top-hat plot configuration constants
  constexpr static int NUM_FIT_POINTS = 7;           // Number of points used for the linear fit.
  constexpr static int BOARD_SLOT = 7;                // Board slot (BOARD_SLOT == 7 for APSAIA, BOARD_SLOT != 8 for ARARA. If all XA do not use it in: if (time_stamp > MIN_WINDOW && time_stamp < MAX_WINDOW) { // us.).
  constexpr static float MIN_WINDOW = -10.0;           //-2.5; // us.
  constexpr static float MAX_WINDOW = 10.0;           //-2.4; // us.
  constexpr static float NS_PER_SAMPLE = 16.0;        // ns.
  constexpr static float BASELINE_PERCENTAGE = 0.1;   // (10 %) 

  // Top-hat plot configuration variables
  int fNumFitPoints;
  int fBoardSlot;
  float fMinWindow;
  float fMaxWindow;
  float fNsPerSample;
  float fNsToUs;
  float fBaselinePercentage;

  // Filepath
  std::string filepath;

  // Event information variables
  int fRun;
  int fSubrun;
  int fEvent;

  // SPEC-TDC SPEC variables
  std::vector<uint64_t> fTDCCh0;
  std::vector<uint64_t> fTDCCh1;
  std::vector<uint64_t> fTDCCh2;
  std::vector<uint64_t> fTDCCh3;
  std::vector<uint64_t> fTDCCh4;
  std::vector<uint64_t> fTDCCh0UTC;
  std::vector<uint64_t> fTDCCh1UTC;
  std::vector<uint64_t> fTDCCh2UTC;
  std::vector<uint64_t> fTDCCh3UTC;
  std::vector<uint64_t> fTDCCh4UTC;

  // Top-hat plot variables
  std::vector<int> fRuns;                               // Run number for each summed waveform.
  std::vector<int> fSubruns;                            // Subrun number for each summed waveform.
  std::vector<int> fEvents;                             // Event number for each summed waveform.
  std::vector<uint16_t> fTimingRefTypes;                // Timing reference type for each summed waveform.
  std::vector<uint16_t> fTimingRefChannels;             // Timing reference channel for each summed waveform.
  std::vector<float> fFlashTimes;                       // Flash time for each summed waveform.
  std::vector<float> fFitSlopes;                        // Slope of the linear fit for each summed waveform.
  std::vector<float> fFitYIntercepts;                   // Y-intercept of the linear fit for each summed waveform.
  std::vector<float> fFitXIntercepts;                   // X-intercept of the linear fit for each summed waveform.
  std::vector<float> fFitChi2s;                         // Chi-square value for each summed waveform.
  std::vector<float> fFitNdfs;                          // Number of degrees of freedom for each summed waveform.
  std::vector<float> fFitReducedChi2s;                  // Reduced chi-square value for each summed waveform.
  std::vector <std::vector <float> > fFitPointsX;       // X values of the fit points for each summed waveform.
  std::vector <std::vector <float> > fFitPointsY;       // Y values of the fit points for each summed waveform.
  std::vector<int> fMaxADCs;                            // ADC max for each summed waveform.
  std::vector<float> fMaxADCsTimes;                     // ADC max for each summed waveform.
  std::vector<size_t> fMaxADCsPositions;                // ADC max position for each summed waveform.
  std::vector<float> fInitTimeStamps;                   // Init time stamp for each summed waveform.
  std::vector<int> fTimeStampNumbers;                   // Time stamp ID number (for this event).
  std::vector<int> fWvfmSizes;                          // Size of the summed waveform.
  std::vector<int> fSummedWfmsBaselines;                // Baseline of each summed waveform.
  std::vector <std::vector <int> > fSummedWvfms;        // Summed waveform constructed (N samples).
  std::vector <std::vector <int> > fRmBaselines;        // Vector with the baselines removed before constructing the summed waveform (64 baselines).
  // Jittering information
  std::vector<float> fJittering;                        // Trigger jittering information.
  std::vector<int> fTDiff;                              // Time difference between a waveforms and the previous one (ns).
  std::vector<int> fWvfmLengths;                        // Length of each waveform (ns).
  std::vector<int> fPrevWvfmLengths;                    // Length of the previous waveform (ns).
};


thp::XARAPUCATopHatPlot::XARAPUCATopHatPlot(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fTDCDecodeLabel = p.get<art::InputTag>("spectdc_product_name", "tdcdecoder");
  fWaveformsModuleLabel = p.get<art::InputTag>("waveforms_module_label", "xarapucadecoder:XARAPUCAChannels");
  fTimingReferenceModuleLabel = p.get<art::InputTag>("timing_reference_module_label", "xarapucadecoder");
  fVerbose = p.get<bool> ("verbose", false);
  fSaveJittering = p.get<bool> ("save_jittering", false);
  fDebugWaveforms = p.get<bool> ("debug_waveforms", false);
  fDebugTDC = p.get<bool> ("debug_TDC", false);
  fDebugRisingEdgeFit = p.get<bool> ("debug_rising_edge_fit", false);
  fDebugJittering = p.get<bool> ("debug_jittering", false);
}

void thp::XARAPUCATopHatPlot::analyze(art::Event const& e)
{
  // Reset event variables
  reset_event_variables();

  // Get event information
  fRun    =  e.id().run();
  fSubrun =  e.id().subRun();
  fEvent  =  e.id().event();
  fTree_events->Fill();

  // Start analyzing the event
  if (fVerbose | fDebugJittering) std::cout << "\nANALYZING " << e.id() << " (" << fRun << "," << fSubrun << ", " << fEvent << ")" << std::endl;

   // ========== SPEC-TDC INFORMATION ==========
   art::Handle<std::vector<sbnd::timing::DAQTimestamp>> TDC_handle;
   e.getByLabel(fTDCDecodeLabel, TDC_handle);

   if (!TDC_handle.isValid() || TDC_handle->size() == 0){
       if (fDebugTDC) std::cout << "No SPECTDC products found. Skip this event." << std::endl;
       return;
   } else{
    if (fDebugTDC) std::cout << "Decoded TDC-SPEC products found: " << TDC_handle->size() << " products." << std::endl;
    const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*TDC_handle);
    
    for (size_t i=0; i<tdc_v.size(); i++){
      auto tdc = tdc_v[i];
      const uint32_t  ch = tdc.Channel();
      const uint64_t  ts = tdc.Timestamp();
      //const uint64_t  offset = tdc.Offset();
      const std::string name  = tdc.Name();

      if (fDebugTDC) std::cout << "Ch " << ch << " " << name
                          << ", ts (ns) = " << ts%uint64_t(1e9)
                          << ", sec (s) = " << ts/uint64_t(1e9)
                          << std::endl;

      
      if(ch == 0){
          fTDCCh0UTC.push_back(ts);
          fTDCCh0.push_back(ts%uint64_t(1e9));
      }
      if(ch == 1){
          fTDCCh1UTC.push_back(ts);
          fTDCCh1.push_back(ts%uint64_t(1e9));
      }
      if(ch == 2){
          fTDCCh2UTC.push_back(ts);
          fTDCCh2.push_back(ts%uint64_t(1e9));
      }
      if(ch == 3){
          fTDCCh3UTC.push_back(ts);
          fTDCCh3.push_back(ts%uint64_t(1e9));
      }
      if(ch == 4){
          fTDCCh4UTC.push_back(ts);
          fTDCCh4.push_back(ts%uint64_t(1e9));
      }
    } 
  }

  // ========== TIMING REFERENCE INFO ========= 
  if (fVerbose) std::cout << "\nAccessing to " << fTimingReferenceModuleLabel << " label." << std::endl;

  uint16_t timing_ref_type = 0;
  uint16_t timing_ref_channel = 0;

  art::Handle<raw::TimingReferenceInfo> timing_ref_handle;
  e.getByLabel(fTimingReferenceModuleLabel, timing_ref_handle);
  
  if (timing_ref_handle.isValid()) {
    timing_ref_type = timing_ref_handle->timingType;
    timing_ref_channel = timing_ref_handle->timingChannel;
    if (fVerbose) std::cout << "Timing reference info: type = " << timing_ref_type << ", channel = " << timing_ref_channel << std::endl;
  }

  // ========== GETS THE WAVEFORMS FROM THE EVENT ==========
  if (fVerbose) std::cout << "\nAccessing to " << fWaveformsModuleLabel << " label." << std::endl;
  
  art::Handle <std::vector <raw::OpDetWaveform> > waveforms_handle;
  e.getByLabel(fWaveformsModuleLabel, waveforms_handle);

  // The art::Handle object is not valid.
  if (!waveforms_handle.isValid()) {
    if (fVerbose) std::cout << "\nWaveforms handle not valid for " << fWaveformsModuleLabel << "." << std::endl;
    if (fVerbose) std::cout << "\n" << waveforms_handle.whyFailed() << std::endl;
  
  // The art::Handle object is empty.
  } else if (waveforms_handle->empty()) {
    if (fVerbose) std::cout << "\nWaveforms handle is empty." << std::endl;
  
  // The art::Handle object is valid and not empty
  } else {
    size_t num_opdetwvfm = waveforms_handle->size();
    if (fVerbose) std::cout << "Valid ( " << num_opdetwvfm << " objects)." << std::endl;

    std::map <float, std::vector<int> > summed_wvfms_map;
    std::map <float, std::vector<int> > summed_wvfm_rm_baselines_map;
  
    // ====== Jittering variables ======
    float t_1 = 0.0;                // us.
    int32_t wvfm_length_1 = 0;      // ns.
    float t_2 = 0.0;                // us.
    int32_t wvfm_length_2 = 0;      // ns.
    int32_t jitt = 0;               // ns.
    int32_t t_diff = 0;             // ns.
    // ===== End jittering variables ======

    // Navigates over all the waveforms of this event.
    for (size_t w = 0; w < num_opdetwvfm; w++) {
      raw::OpDetWaveform wvfm = waveforms_handle->at(w);
      unsigned int channel_number = wvfm.ChannelNumber();
      unsigned int board_slot = channel_number / 100;
      float time_stamp = wvfm.TimeStamp();

      // ARARA reversing waveforms
      if (board_slot != 7) {
        for (size_t i = 0; i < wvfm.size(); i++) {
          wvfm[i] = -wvfm[i];
        }
      }

      // Filters waveforms near to the zero.
      if (time_stamp > MIN_WINDOW && time_stamp < MAX_WINDOW) { // us.
        
        if (fDebugWaveforms | (fDebugJittering && channel_number == 700)) std::cout << "[WVFM_ID: " << w << "]\tChannel number: " << channel_number << "\tTimeStamp: " << time_stamp << " us\tBoard slot: " << board_slot << std::endl;
        
        // Jittering computation (for APSAIA only).
        if (fSaveJittering) {
          if (channel_number == 700) {
            t_2 = time_stamp;
            wvfm_length_2 = static_cast<int32_t>(wvfm.size() * NS_PER_SAMPLE);   // ns.

            if (wvfm_length_1 == 0) {
              t_diff = DEFAULT_INT;
              jitt = DEFAULT_INT;
            } else {
              t_diff = static_cast<int32_t>(std::lround((t_2 - t_1) * US_TO_NS));
              jitt = t_diff - wvfm_length_1; // ns.
            }

            fJittering.push_back(jitt);
            fTDiff.push_back(t_diff);
            fWvfmLengths.push_back(wvfm_length_2);
            fPrevWvfmLengths.push_back(wvfm_length_1);

            if (fDebugJittering) {
              std::cout << "\t\t JITTERING COMPUTATION: t_1 = " << t_1 << " us, t_2 = " << t_2 << " us --> t_diff = " << t_diff << " ns, wvfm_length_1 = " << wvfm_length_1 << " ns, wvfm_length_2 = " << wvfm_length_2 << " ns --> jitt = " << jitt << " ns." << std::endl;
            }
            t_1 = t_2;
            wvfm_length_1 = wvfm_length_2;
          }
        }

        // End jittering computation.

        if (summed_wvfms_map.find(time_stamp) == summed_wvfms_map.end()) {
          summed_wvfms_map[time_stamp] = std::vector<int>(wvfm.size(), 0);
          summed_wvfm_rm_baselines_map[time_stamp] = std::vector<int>();
        }

        // Compute baseline.
        float rm_baseline = 0;
        size_t num_baseline_samples = wvfm.size() * BASELINE_PERCENTAGE;
        for (size_t s = 0; s < num_baseline_samples; s++) {
          rm_baseline += wvfm[s];
        }
        rm_baseline /= num_baseline_samples;
        if (fDebugWaveforms) std::cout << "REMOVING BASELINE [" << num_baseline_samples << " samples (" << BASELINE_PERCENTAGE << "%)]: " << rm_baseline << " ADCs." << std::endl; 
        summed_wvfm_rm_baselines_map[time_stamp].push_back(rm_baseline);

        // Sum waveforms element by element removing the baseline.
        for (size_t s = 0; s < wvfm.size(); s++) {
          summed_wvfms_map[time_stamp][s] += std::round(wvfm[s] - rm_baseline);
          //if (fdebug) std::cout << "\t\t\t" << summed_wvfms_map[time_stamp][s] << ":::\t\t" << std::round(wvfm[s] - rm_baseline) <<  "\t\t = " << wvfm[s] << "\t - \t" << rm_baseline << std::endl;
        }
      }
    } // End sum waveforms loop. 
    
    // ========== COMPUTES STATISTICS FOR EACH SUMMED WAVEFORM ==========
    int num_time_stamp = 0;
    for (const auto& [time_stamp, summed_wvfm] : summed_wvfms_map) {

      // Computes the baseline for the summed waveform.
      float summed_wvfm_baseline = 0;
      size_t num_baseline_samples = summed_wvfm.size() * BASELINE_PERCENTAGE;
      for (size_t s = 0; s < num_baseline_samples; s++) {
          summed_wvfm_baseline += summed_wvfm[s];
      }
      summed_wvfm_baseline /= num_baseline_samples;
      if (fVerbose) std::cout << "SUMMED BASELINE [" << num_baseline_samples << " samples (" << BASELINE_PERCENTAGE << "%)]: " << std::round(summed_wvfm_baseline) << " ADCs." << std::endl; 
      
      // Find the rising edge using a linear fit and the maximum ADC value information.
      thp::MaxADCInfo max_ADC_info;
      thp::RisingEdgeFit rising_edge_fit;
      get_rising_edge(summed_wvfm, time_stamp, max_ADC_info, rising_edge_fit);

      // Store the values.
      // - Event information.
      fRuns.push_back(fRun);
      fSubruns.push_back(fSubrun);
      fEvents.push_back(fEvent);

      // - Timing reference information.
      fTimingRefTypes.push_back(timing_ref_type);
      fTimingRefChannels.push_back(timing_ref_channel);

      // - Maximum ADC information.
      fMaxADCs.push_back(max_ADC_info.max_ADC);
      fMaxADCsTimes.push_back(max_ADC_info.max_ADC_time);
      fMaxADCsPositions.push_back(max_ADC_info.max_ADC_pos);

      // - Rising edge fit information.
      fFlashTimes.push_back(rising_edge_fit.flash_time);
      fFitSlopes.push_back(rising_edge_fit.slope);
      fFitYIntercepts.push_back(rising_edge_fit.y_intercept);
      fFitXIntercepts.push_back(rising_edge_fit.x_intercept);
      fFitChi2s.push_back(rising_edge_fit.chi2);
      fFitNdfs.push_back(rising_edge_fit.ndf);
      fFitReducedChi2s.push_back(rising_edge_fit.reduced_chi2);
      fFitPointsX.push_back(rising_edge_fit.fit_points_x);
      fFitPointsY.push_back(rising_edge_fit.fit_points_y);
      
      // - Waveform information.
      fInitTimeStamps.push_back(time_stamp);
      fTimeStampNumbers.push_back(num_time_stamp);
      fWvfmSizes.push_back(summed_wvfm.size());
      fSummedWfmsBaselines.push_back(summed_wvfm_baseline);
      fSummedWvfms.push_back(summed_wvfm);
      fRmBaselines.push_back(summed_wvfm_rm_baselines_map[time_stamp]);

      if (fVerbose)  {
        std::cout << "DATA: ["<< num_time_stamp << "]" << std::endl;
        std::cout << "\tMaximum ADC time: " << max_ADC_info.max_ADC_time << " us\tMaximum ADC pos: " << max_ADC_info.max_ADC_pos << ".\tMaximum ADC: " << max_ADC_info.max_ADC << std::endl;
        std::cout << "\tFlash time: " << rising_edge_fit.flash_time << " us\tX-Intercept: " << rising_edge_fit.x_intercept << "\tSlope: " << rising_edge_fit.slope << "\tY-Intercept: " << rising_edge_fit.y_intercept << std::endl;
        if (fDebugRisingEdgeFit) {
          std::cout << "\t\tFitting points: " << std::endl;
          for (size_t j = 0; j < NUM_FIT_POINTS; j++) {
            std::cout << "\t\t\tX: " << fFitPointsX.at(num_time_stamp)[j] << "\tY: " << fFitPointsY.at(num_time_stamp)[j] << std::endl;
          }
        }
      }
      num_time_stamp++;
      }
  }
  
  if (fDebugTDC) {
    std::cout << "Number of ch2 found: " << fTDCCh2UTC.size() << std::endl;
    std::cout << "Number of ch4 found: " << fTDCCh4UTC.size() << std::endl;
  }
  
  fTree_thp->Fill();
}

void thp::XARAPUCATopHatPlot::get_rising_edge(std::vector<int> const& summed_wvfm, float time_stamp, MaxADCInfo& max_ADC_info, RisingEdgeFit& rising_edge_fit)
{
  // Find the maximum ADC value and its position.
  max_ADC_info.max_ADC_pos = std::distance(summed_wvfm.begin(), std::max_element(summed_wvfm.begin(), summed_wvfm.end()));
  max_ADC_info.max_ADC = summed_wvfm[max_ADC_info.max_ADC_pos];
  max_ADC_info.max_ADC_time = max_ADC_info.max_ADC_pos * NS_PER_SAMPLE * NS_TO_US + time_stamp; // us.

  if (fDebugRisingEdgeFit) std::cout << "Maximum ADC: " << max_ADC_info.max_ADC << "\tPosition: " << max_ADC_info.max_ADC_pos << "\tTime: " << max_ADC_info.max_ADC_time << " us." << std::endl;

  // Find the flash time using a linear fit.
  rising_edge_fit.slope = -std::numeric_limits<float>::infinity();
  rising_edge_fit.fit_points_x.resize(NUM_FIT_POINTS);
  rising_edge_fit.fit_points_y.resize(NUM_FIT_POINTS);
  std::vector<float> fit_points_x(NUM_FIT_POINTS, 0.0);
  std::vector<float> fit_points_y(NUM_FIT_POINTS, 0.0);
  for (int s = max_ADC_info.max_ADC_pos; s >= NUM_FIT_POINTS; s--) {
    TLinearFitter fitter(1, "pol1");
    for (int i = s-NUM_FIT_POINTS, j = 0; i < s; i++, j++) {
      Double_t x = i;
      Double_t y = summed_wvfm[i];
      fitter.AddPoint(&x, y);
      fit_points_x[j] = x;
      fit_points_y[j] = y;
    }

    fitter.Eval();
    TVectorD params(2);
    fitter.GetParameters(params);

    float m = params[1];                                                // Slope
    float b = params[0];                                                // Y-intercept
    int ndf = fitter.GetNpoints() - fitter.GetNumberFreeParameters();   // Number of degrees of freedom (ndf)
    float chi2 = fitter.GetChisquare();                                 // Chi-square value
    float reduced_chi2 = (ndf > 0) ? chi2 / ndf : 0.0;                  // Reduced chi-square

    if (fDebugRisingEdgeFit) {
      std::cout << "Slope: " << m << "\tY-intercept: " << b << "\tChi2: " << chi2 << "\tNDF: " << ndf << "\tReduced Chi2: " << reduced_chi2 << std::endl;
    }

    if (m == 0) {
      if (fDebugRisingEdgeFit) std::cout << "Zero slope, cannot calculate the X-axis intercept." << std::endl;
    } else if (m > rising_edge_fit.slope) {
      rising_edge_fit.slope = m;
      rising_edge_fit.y_intercept = b;
      rising_edge_fit.x_intercept = -b / m;               // X-intercept (-b/m) is the actual flash time [us] 
      rising_edge_fit.chi2 = chi2;
      rising_edge_fit.ndf = ndf;
      rising_edge_fit.reduced_chi2 = reduced_chi2;
      for (int j = 0; j < NUM_FIT_POINTS; j++) {
        rising_edge_fit.fit_points_x[j] = fit_points_x[j];
        rising_edge_fit.fit_points_y[j] = fit_points_y[j];
      }
      if (fDebugRisingEdgeFit) {
        std::cout << "New X-intercept: " << rising_edge_fit.x_intercept << std::endl;
        std::cout << "Fitting points: " << std::endl;
        for (size_t j = 0; j < NUM_FIT_POINTS; j++) {
          std::cout << "\tX: " << rising_edge_fit.fit_points_x[j] << "\tY: " << rising_edge_fit.fit_points_y[j] << std::endl;
        }
      }
    }
  }
  rising_edge_fit.flash_time = rising_edge_fit.x_intercept * NS_PER_SAMPLE * NS_TO_US + time_stamp; // us.
}


void thp::XARAPUCATopHatPlot::beginJob()
{
  // Implementation of optional member function here.
  if (fVerbose) std::cout << "Time window of " << MAX_WINDOW - MIN_WINDOW << " us: [" << MIN_WINDOW << ", " << MAX_WINDOW << "] us." << std::endl;
  if (fVerbose) std::cout << "Nanoseconds per sample: " << NS_PER_SAMPLE << " ns." << std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  fTree_config = tfs->make<TTree>("thp_config", "Top-hat Plot configuration");
  fTree_events = tfs->make<TTree>("thp_events", "Top-hat Plot events");
  fTree_thp = tfs->make<TTree>("thp", "Top-hat Plot data");

  /********** CONFIGURATION TREE **********/
  fTree_config->Branch("NUM_FIT_POINTS", &fNumFitPoints, "NUM_FIT_POINTS/I");
  fTree_config->Branch("MIN_WINDOW", &fMinWindow, "MIN_WINDOW/F");
  fTree_config->Branch("MAX_WINDOW", &fMaxWindow, "MAX_WINDOW/F");
  fTree_config->Branch("NS_PER_SAMPLE", &fNsPerSample, "NS_PER_SAMPLE/F");
  fTree_config->Branch("NS_TO_US", &fNsToUs, "NS_TO_US/F");
  fTree_config->Branch("BASELINE_PERCENTAGE", &fBaselinePercentage, "BASELINE_PERCENTAGE/F");
  fNumFitPoints = NUM_FIT_POINTS;
  fBoardSlot = BOARD_SLOT;
  fMinWindow = MIN_WINDOW;
  fMaxWindow = MAX_WINDOW;
  fNsPerSample = NS_PER_SAMPLE;
  fNsToUs = NS_TO_US;
  fBaselinePercentage = BASELINE_PERCENTAGE;
  fTree_config->Fill();

  /********** EVENTS TREE **********/
  fTree_events->Branch("run", &fRun);
  fTree_events->Branch("subrun", &fSubrun);
  fTree_events->Branch("event", &fEvent);
  fTree_events->Branch("filepath", &filepath);

  /********** TOP-HAT PLOT TREE **********/
  // Event information.
  fTree_thp->Branch("run", &fRuns);
  fTree_thp->Branch("subrun", &fSubrun);
  fTree_thp->Branch("event", &fEvents);
  // TDC SPEC information.
  fTree_thp->Branch("tdc_ch0", &fTDCCh0);
  fTree_thp->Branch("tdc_ch1", &fTDCCh1);
  fTree_thp->Branch("tdc_ch2", &fTDCCh2);
  fTree_thp->Branch("tdc_ch3", &fTDCCh3);
  fTree_thp->Branch("tdc_ch4", &fTDCCh4);
  fTree_thp->Branch("tdc_ch0_UTC", &fTDCCh0UTC);
  fTree_thp->Branch("tdc_ch1_UTC", &fTDCCh1UTC);
  fTree_thp->Branch("tdc_ch2_UTC", &fTDCCh2UTC);
  fTree_thp->Branch("tdc_ch3_UTC", &fTDCCh3UTC);
  fTree_thp->Branch("tdc_ch4_UTC", &fTDCCh4UTC);
  // Timing reference information.
  fTree_thp->Branch("timing_ref_type", &fTimingRefTypes);
  fTree_thp->Branch("timing_ref_channel", &fTimingRefChannels);
  // Maximum ADC information.
  fTree_thp->Branch("max_ADC", &fMaxADCs);
  fTree_thp->Branch("max_ADC_time", &fMaxADCsTimes);
  fTree_thp->Branch("max_ADC_pos", &fMaxADCsPositions);
  // Rising edge fit information.
  fTree_thp->Branch("flash_time", &fFlashTimes);
  fTree_thp->Branch("fit_slope", &fFitSlopes);
  fTree_thp->Branch("fit_y_intercept", &fFitYIntercepts);
  fTree_thp->Branch("fit_x_intercept", &fFitXIntercepts);
  fTree_thp->Branch("fit_chi2", &fFitChi2s);
  fTree_thp->Branch("fit_ndf", &fFitNdfs);
  fTree_thp->Branch("fit_reduced_chi2", &fFitReducedChi2s);
  fTree_thp->Branch("fit_points_x", &fFitPointsX);
  fTree_thp->Branch("fit_points_y", &fFitPointsY);
  // Waveform information.
  fTree_thp->Branch("init_time_stamp", &fInitTimeStamps);
  fTree_thp->Branch("wvfm_size", &fWvfmSizes);
  fTree_thp->Branch("num_time_stamp", &fTimeStampNumbers);
  fTree_thp->Branch("summed_wvfm_baseline", &fSummedWfmsBaselines);
  fTree_thp->Branch("summed_wvfm", &fSummedWvfms);
  fTree_thp->Branch("rm_baselines", &fRmBaselines);
  // Jittering information.
  fTree_thp->Branch("jittering", &fJittering);
  fTree_thp->Branch("t_diff", &fTDiff);
  fTree_thp->Branch("wvfm_length", &fWvfmLengths);
  fTree_thp->Branch("prev_wvfm_length", &fPrevWvfmLengths);
  // Filepath
  fTree_thp->Branch("filepath", &filepath);
}

void thp::XARAPUCATopHatPlot::endJob() {

}

void thp::XARAPUCATopHatPlot::respondToOpenInputFile(art::FileBlock const &fb) {
  filepath = fb.fileName();
}

void thp::XARAPUCATopHatPlot::reset_event_variables() {
  /********** EVENTS TREE **********/
  fRun = -1;
  fSubrun = -1;
  fEvent = -1;

  /********** TOP-HAT PLOT TREE **********/
  // Event information.
  fRuns.clear();
  fSubruns.clear();
  fEvents.clear();
  // TDC SPEC information.
  fTDCCh0.clear();
  fTDCCh1.clear();
  fTDCCh2.clear();
  fTDCCh3.clear();
  fTDCCh4.clear();
  fTDCCh0UTC.clear();
  fTDCCh1UTC.clear();
  fTDCCh2UTC.clear();
  fTDCCh3UTC.clear();
  fTDCCh4UTC.clear();
  // Timing reference information.
  fTimingRefTypes.clear();
  fTimingRefChannels.clear();
  // Maximum ADC information.
  fMaxADCs.clear();
  fMaxADCsPositions.clear();
  fMaxADCsTimes.clear();
  // Rising edge fit information.
  fFlashTimes.clear();
  fFitSlopes.clear();
  fFitYIntercepts.clear();
  fFitXIntercepts.clear();
  fFitChi2s.clear();
  fFitNdfs.clear();
  fFitReducedChi2s.clear();
  fFitPointsX.clear();
  fFitPointsY.clear();
  // Waveform information.
  fInitTimeStamps.clear();
  fWvfmSizes.clear();
  fTimeStampNumbers.clear();
  fSummedWfmsBaselines.clear();
  fSummedWvfms.clear();
  fRmBaselines.clear();
  // Jittering information.
  fJittering.clear();
  fTDiff.clear();
  fWvfmLengths.clear();
  fPrevWvfmLengths.clear();
}

DEFINE_ART_MODULE(thp::XARAPUCATopHatPlot)
