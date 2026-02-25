////////////////////////////////////////////////////////////////////////
// Class:       CALLOS
// Plugin Type: analyzer (Unknown Unknown)
// File:        CALLOS_module.cc
// 
// Generated at Fri Feb  9 19:47:11 2024 by Rodrigo Alvarez Garrote using cetskelgen
// from  version .
// 
// This module is a Calibration Analyzer for Low Light Optical Signals
// CALLOS is made to read raw waveforms from raw and deconvolved stages and 
// compute the charges and average waveform for each optical channel selected.
// The resulting histograms and NTuples are saved as output.
// A SPE and a vector of charges is produced for each channel and for ALL events processed.
// ROI tools and PEAK finders are defined as Tools with templates to allow for
// different implementations to be used. 
// 
// The aim of this module is to be as fast as possible, hence the conscious decision of not using 
// if checks when possible. Specially for analyze() function(run on each event).
// 
// Note: the module assumes raw::waveform channels go from 0 to NChannels-1. 
// If needed, change from std::vector to std::map.
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
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Utilities/make_tool.h"
#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <TROOT.h>
// #include <TVirtualFFT.h>
#include <TMath.h>
#include <algorithm> 

#include "TTree.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "Tools/ROIFinderAlg.hh"
#include "Evd.hh"
#include "AverageWaveform.h"

namespace callos {
  class CALLOS;
}

class callos::CALLOS : public art::EDAnalyzer {
public:
  explicit CALLOS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CALLOS(CALLOS const&) = delete;
  CALLOS(CALLOS&&) = delete;
  CALLOS& operator=(CALLOS const&) = delete;
  CALLOS& operator=(CALLOS&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // placeholders
  // void beginJob(art::Run& run) override ;
  // void beginRun(art::Run& run) override ;
  // void beginSubRun(art::Run& run) override ;
  
  // void endSubRun(art::SubRun& sr) override ;
  // void endRun(art::Run& run) override ;
  void endJob() override ;

  // Tools shall be initialized here 
private:

  // Declare member data here.

  // fhicl parameters
  std::string fInputLabel;
  //std::vector<std::string> fPDTypes;
  bool fDEBUG;
  bool fFull_AverageWaveform;
  bool fPEMode;
  // Channel number
  std::vector<int> fChNumbers;
  // Selected channels
  std::vector<int> fSelectedChannels;
  int number_of_events_analyzed = 0; // Variable to store the number of events analyzed
  int number_of_events_analyzed_total;
  int number_waveforms = 0;
  int number_waveforms_total;
  int number_of_signals = 0;
  int number_of_signals_total;
  std::vector<int> waveform_wrong_baseline_channel;
  std::vector<int> waveform_no_signals_channel;
  // Number of selected channels
  int fNSelectedChannels;
  int fMaxADC;
  std::vector<int> fDictionary; //Need to assign a number to each channel
  std::vector<int> fChannels; // Max charge expected for each channel
  // Size of the Region of Interest(ROI) of each signal, must be smaller than RawWaveforms size. rodrigoa: maybe should be a tool parameter
  int fROISamples;

  // Number of samples from start of the ROI to peak. rodrigoa: maybe should be a tool parameter
  int fStartToPeak;

  // Specific channels to be analyzed (only used if not empty)
  //std::vector<int> fSpecificChannels;

  std::vector<TH1F*> hAverageWaveformHistograms_fewPEs;


  TH1F* hEventNumber;
  TH1F* hWaveformsNumber;
  TH1F* hSignalsNumber;

  std::vector<TH1F*> hFull_AverageWaveformHistograms;
  std::vector<TH1F*> hAverageWaveformHistograms_highPEs;

  std::vector<TH1F*> hAverageWaveformHistograms_1PE;

  //Charge container
  std::vector<std::vector<float>> Charge_SelectedChannels;
  std::vector<std::vector<float>> Charge_highPEs_SelectedChannels;
  std::vector<std::vector<float>> Amplitude_highPEs_SelectedChannels;
  std::vector<std::vector<float>> TimeStamp_SelectedChannels;
   std::vector<std::vector<float>> Event_number_SelectedChannels;
  //Pedestal charge container
  std::vector<std::vector<float>> PedestalCharge_SelectedChannels;
  //Peak container
  std::vector<std::vector<float>> Peak_SelectedChannels;
  //PeakTick container
  std::vector<std::vector<int>> PeakTick_SelectedChannels;
  //Baseline container
  std::vector<std::vector<float>> Baseline_SelectedChannels;
  //BaselineSTD container
  std::vector<std::vector<float>> BaselineSTD_SelectedChannels;
  // Average Waveforms container.
  std::vector<std::vector<std::vector<float>>> AverageWaveform_fewPEs_SelectedChannels;
  std::vector<std::vector<std::vector<float>>> Full_AverageWaveform_SelectedChannels;
  std::vector<std::vector<std::vector<float>>> AverageWaveform_highPEs_SelectedChannels;

  //PE Mode Avg Waveform 
  std::vector<std::vector<std::vector<float>>> PE_1_AverageWaveform_SelectedChannels;
  std::vector<std::vector<float>> PE_1_Peak_SelectedChannels;

  float EstimateBaseline(std::vector<float>& wf);
  float BaselineSTD(std::vector<float>& wf, float baseline);
  bool  SingleinPretrigger(std::vector<float>& wf);

  // Tool pointers
  std::unique_ptr<callos::ROIFINDERALG> fROIFinderAlgPtr;

  // le  tree
  TTree *fTree;
  TTree* fSummaryTree;
};


callos::CALLOS::CALLOS(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fInputLabel = p.get<std::string>("InputLabel","opdaq");
  //fPDTypes = p.get<std::vector<std::string>>("PDTypes",{"xarapuca_vuv"});
  fROISamples = p.get<int>("ROI_samples", 1000);
  fStartToPeak = p.get<int>("StartToPeak",200);
  //fSpecificChannels = p.get<std::vector<int>>("SpecificChannels", {});
  fDEBUG = p.get<bool>("DEBUG", true);
  fFull_AverageWaveform = p.get<bool>("Full_AverageWaveform", false);
  // get map info
  // for each pdtype in fPDTypes, get the channels and add them to the list of selected channels
  
  // fSelectedChannels.reserve(64); //hardcoded for now

  //   for (int i = 700; i <= 763; i++) {
  //   fSelectedChannels.push_back(i);
  // }
  fSelectedChannels = p.get<std::vector<int>>("Channels", {});
  // fSelectedChannels.reserve(32); //hardcoded for now

  //   for (int i = 700; i <= 731; i++) {
  //   fSelectedChannels.push_back(i);
  // }

  fNSelectedChannels = fSelectedChannels.size();

  // Initialize the average waveforms container
  if (fFull_AverageWaveform){
    Full_AverageWaveform_SelectedChannels.resize(fNSelectedChannels);
    Baseline_SelectedChannels.resize(fNSelectedChannels);
    Peak_SelectedChannels.resize(fNSelectedChannels);
  }
  else{
    PedestalCharge_SelectedChannels.resize(fNSelectedChannels);
    Charge_SelectedChannels.resize(fNSelectedChannels);
    Baseline_SelectedChannels.resize(fNSelectedChannels);
    BaselineSTD_SelectedChannels.resize(fNSelectedChannels);
    AverageWaveform_fewPEs_SelectedChannels.resize(fNSelectedChannels);
    Peak_SelectedChannels.resize(fNSelectedChannels);
    PeakTick_SelectedChannels.resize(fNSelectedChannels);
    AverageWaveform_highPEs_SelectedChannels.resize(fNSelectedChannels);
    Charge_highPEs_SelectedChannels.resize(fNSelectedChannels);
    Event_number_SelectedChannels.resize(fNSelectedChannels);
    Amplitude_highPEs_SelectedChannels.resize(fNSelectedChannels);
    TimeStamp_SelectedChannels.resize(fNSelectedChannels);
  }
  if (fPEMode){
    PE_1_AverageWaveform_SelectedChannels.resize(fNSelectedChannels);
    PE_1_Peak_SelectedChannels.resize(fNSelectedChannels);
  }
  // Initialize the ROIFinderAlg tool
  fROIFinderAlgPtr = art::make_tool<callos::ROIFINDERALG>(p.get<fhicl::ParameterSet>("ROIFinderAlg"));

  // Now that we know the number of selected channels, we can create the tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("callos","callos");

  fSummaryTree = tfs->make<TTree>("summary", "summary");
  fSummaryTree->Branch("EventNumber", &number_of_events_analyzed_total, "EventNumber/I");
  fSummaryTree->Branch("SignalsNumber", &number_of_signals_total, "SignalsNumber/I");
  fSummaryTree->Branch("WaveformsNumber", &number_waveforms_total, "WaveformsNumber/I");

  if (fFull_AverageWaveform)
  {
    art::TFileDirectory fullAverageDir = tfs->mkdir("Full_AverageWaveform");
    for (int i=0; i<fNSelectedChannels; i++)
    {
      hFull_AverageWaveformHistograms.push_back(fullAverageDir.make<TH1F>(
          ("hFull_AverageWaveform_" + std::to_string(fSelectedChannels[i])).c_str(),
          ("Full Average Waveform Channel " + std::to_string(fSelectedChannels[i])).c_str(),
          fROISamples, 0, fROISamples));
    }
    for (int i=0; i<fNSelectedChannels; i++)
    {
      fTree->Branch(("Baseline_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Baseline_SelectedChannels[i]);
      fTree->Branch(("Peak_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Peak_SelectedChannels[i]);
    }
  }
  else{
    for (int i = 0; i < fNSelectedChannels; ++i) 
    {
      hAverageWaveformHistograms_fewPEs.push_back(tfs->make<TH1F>(
          ("hAverageWaveform_" + std::to_string(i)).c_str(),
          ("Average Waveform Channel " + std::to_string(i)).c_str(),
          fROISamples, 0, fROISamples));
      hAverageWaveformHistograms_highPEs.push_back(tfs->make<TH1F>(
          ("hAverageWaveform_highPEs_" + std::to_string(i)).c_str(),
          ("Average Waveform High PEs Channel " + std::to_string(i)).c_str(),
          fROISamples, 0, fROISamples));
      if (fPEMode)
      {
        hAverageWaveformHistograms_1PE.push_back(tfs->make<TH1F>(
            ("hAverageWaveform_1PE_" + std::to_string(i)).c_str(),
            ("Average Waveform Channel " + std::to_string(i)).c_str(),
            fROISamples, 0, fROISamples));
      }
      
    }
    fTree->Branch("Waveform_wrong_baseline_channel", &waveform_wrong_baseline_channel);
    fTree->Branch("Waveform_no_signals_channel", &waveform_no_signals_channel);
    for (int i=0; i<fNSelectedChannels; i++)
    {

      fTree->Branch(("Charge_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Charge_SelectedChannels[i]);
      fTree->Branch(("Peak_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Peak_SelectedChannels[i]);
      fTree->Branch(("PedestalCharge_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PedestalCharge_SelectedChannels[i]);
      fTree->Branch(("Signal_time_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PeakTick_SelectedChannels[i]);
      fTree->Branch(("Charge_highPEs_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Charge_highPEs_SelectedChannels[i]);
      fTree->Branch(("EventNumber_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Event_number_SelectedChannels[i]);
      fTree->Branch(("Amplitude_highPEs_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Amplitude_highPEs_SelectedChannels[i]);
      fTree->Branch(("TimeStamp_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &TimeStamp_SelectedChannels[i]);
      fTree->Branch(("Baseline_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Baseline_SelectedChannels[i]);
      fTree->Branch(("BaselineSTD_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &BaselineSTD_SelectedChannels[i]);
      if (fPEMode)
      {
        fTree->Branch(("Peak_1PE_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PE_1_Peak_SelectedChannels[i]);
      }

    }

  }
}

// compute the baseline of the waveform
float callos::CALLOS::EstimateBaseline(std::vector<float> &wf) {
  return std::accumulate(wf.begin(), wf.end(), 0.0f) / wf.size();
}

// compute the standard deviation of the baseline
float callos::CALLOS::BaselineSTD(std::vector<float> &wf, float baseline) {
  double mean = baseline;
  double variance = 0.0;
  for (const auto& value : wf) {
    variance += (value - mean) * (value - mean);
  }

  if (wf.size() > 0)
  { 
    variance /= wf.size();
  }
  else
  {
    // std::cout << "SimpleROIAlg::BaselineSTD: wf.size() is 0"<<std::endl;
    return 0.0;
  }
  return std::sqrt(variance);
}
// bool callos::CALLOS::SingleinPretrigger(std::vector<float> &wf) {
//   std::vector<float> wf0;
//   wf0.assign(wf.begin(), wf.begin()+0.16 * wf.size());
//   float baselina0 = EstimateBaseline(wf0);
//   float baselinestd0 = BaselineSTD(wf0, baselina0);
//   std::vector<float> wf1;
//   wf1.assign(wf.begin(), wf.begin()+0.04 * wf.size());
//   float baseline1 = EstimateBaseline(wf1);
//   float baselinestd1 = BaselineSTD(wf1, baseline1);
//   std::vector<float> wf2;
//   float baseline2 = EstimateBaseline(wf2);
//   float baselinestd2 = BaselineSTD(wf2, baseline2);
//   wf2.assign(wf.begin()+0.04 * wf.size(), wf.begin()+0.08 * wf.size());
//   std::vector<float> wf3;
//   float baseline3 = EstimateBaseline(wf3);
//   float baselinestd3 = BaselineSTD(wf3, baseline3);
//   wf3.assign(wf.begin()+0.08 * wf.size(), wf.begin()+0.12 * wf.size());
//   std::vector<float> wf4;
//   float baseline4 = EstimateBaseline(wf4);
//   float baselinestd4 = BaselineSTD(wf4 baseline4);
//   wf4.assign(wf.begin()+0.12 * wf.size(), wf.begin()+0.16 * wf.size());

//   float min_baseline = std::min({baselina0, baseline1, baseline2, baseline3, baseline4});

//   double max_position = std::distance(wf0.begin(), std::max_element(wf0.begin(), wf0.end()));
//   std::vector<float> wf_max;
//   wf_max.assign(wf0.begin()+max_position-3, wf0.begin()+max_position+8);
//   float baseline_max = EstimateBaseline(wf_max);
//   float baselinestd_max = BaselineSTD(wf_max, baseline_max);
//   float diff_baseline_max = std::abs(baseline_max - min_baseline);
//   if (diff_baseline_max > 1){
//     return false; // Skip noisy waveforms
//   }
//   return true;
// }

bool callos::CALLOS::SingleinPretrigger(std::vector<float> &wf)
{
    if (wf.size() < 20) return false;   // protección básica

    const size_t N = wf.size();

    // Fracciones de ventana
    constexpr float PRE_FRAC   = 0.16f;
    constexpr float SLICE_FRAC = 0.04f;

    const size_t preN   = static_cast<size_t>(PRE_FRAC   * N);
    const size_t sliceN = static_cast<size_t>(SLICE_FRAC * N);

    // --- Helper lambda para extraer subventanas ---
    auto make_slice = [&](size_t i0, size_t i1)
    {
        i0 = std::min(i0, N);
        i1 = std::min(i1, N);
        return std::vector<float>(wf.begin() + i0, wf.begin() + i1);
    };

    // --- Ventana completa del pretrigger ---
    auto wf0 = make_slice(0, preN);
    float baseline0 = EstimateBaseline(wf0);
    //const float std0      = BaselineSTD(wf0, baseline0);

    // --- Subventanas ---
    std::vector<float> baselines;

    for (size_t i = 0; i < 4; ++i)
    {
        size_t i0 = i * sliceN;
        size_t i1 = (i + 1) * sliceN;

        auto slice = make_slice(i0, i1);
        float b    = EstimateBaseline(slice);
        baselines.push_back(b);
    }

    // Baseline mínimo
    float min_baseline =
        std::min({ baseline0,
                   baselines[0],
                   baselines[1],
                   baselines[2],
                   baselines[3] });

    // --- Buscar máximo en el pretrigger ---
    auto it_max = std::max_element(wf0.begin(), wf0.end());
    if (it_max == wf0.end()) return false;

    size_t max_pos = std::distance(wf0.begin(), it_max);

    // Ventana alrededor del máximo (con protección de bordes)
    constexpr int LEFT  = 2;
    constexpr int RIGHT = 9;

    size_t i0 = (max_pos > LEFT) ? max_pos - LEFT : 0;
    size_t i1 = std::min(max_pos + RIGHT, wf0.size());

    auto wf_max = std::vector<float>(wf0.begin() + i0,
                                           wf0.begin() + i1);

    float baseline_max = EstimateBaseline(wf_max);
    //const float std_max      = BaselineSTD(wf_max, baseline_max);

    float diff = std::abs(baseline_max - min_baseline);

    // --- Corte final ---
    constexpr float BASELINE_DIFF_THR = 1.0f;

    if (diff > BASELINE_DIFF_THR)
        return false;   // reject waveform

    return true;        // accept waveform
}

void callos::CALLOS::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Increment the number of events analyzed
  number_of_events_analyzed++;
  // Load the waveforms
  art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
  e.getByLabel(fInputLabel, wfHandle);

  //Load

  if (fFull_AverageWaveform){
    for (int i=0; i<fNSelectedChannels; i++)  
    {
      Baseline_SelectedChannels[i]={}; // Clear the baseline container
      Peak_SelectedChannels[i]={}; // Clear the peak container
    }
  }
  else{
    for (int i=0; i<fNSelectedChannels; i++)  
    {
      Charge_SelectedChannels[i]={}; // Clear the charge container
      Peak_SelectedChannels[i]={}; // Clear the peak container
      PedestalCharge_SelectedChannels[i]={}; // Clear the pedestal charge container
      PeakTick_SelectedChannels[i]={}; // Clear the peakTick container
      Baseline_SelectedChannels[i]={}; // Clear the baseline container
      BaselineSTD_SelectedChannels[i]={}; // Clear the baselineSTD container
      Charge_highPEs_SelectedChannels[i]={}; // Clear the high PE charge container
      Event_number_SelectedChannels[i]={}; // Clear the event number container
      Amplitude_highPEs_SelectedChannels[i]={}; // Clear the high PE amplitude container
      TimeStamp_SelectedChannels[i]={}; // Clear the timestamp container
      if (fPEMode)
      {
        PE_1_Peak_SelectedChannels[i]={}; // Clear the 1 PE peak container
      }
    }
    waveform_wrong_baseline_channel.clear();
    waveform_no_signals_channel.clear();
  }

  if (!wfHandle.isValid()) 
  {
   mf::LogError("SBNDCALLOS")<<"Input waveforms with input label "<<fInputLabel<<" couldn't be loaded..."<<std::endl;
   throw cet::exception("SBNDCALLOS") << "Input waveforms with input label " << fInputLabel << " not found\n";
  }
  
  // Get the Raw waveforms
  // int iWvf=0;
    for(auto const& wf : *wfHandle)
    {
      // Get raw wvf info
      // const int wfChannel = wf.ChannelNumber();
      // if (wfChannel > 699 && wfChannel < 732) number_waveforms++;


      const int wfChannel = wf.ChannelNumber();

      if (std::find(fSelectedChannels.begin(),
                    fSelectedChannels.end(),
                    wfChannel) != fSelectedChannels.end()) {
          number_waveforms++;
      }
      else continue;

      // if not in selected channels, skip

      // std::cout<<"CALLOS: Event "<<e.id().event()<<" Channel "<<wfChannel<<std::endl;
      const double wfStartTime = wf.TimeStamp();
      // if (wfChannel == 701){
      //   std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << " TimeStamp " << wfStartTime << std::endl;
      // }

      // select only relevant channels (this allows to sepparate diferent XAs, PMTs, etc based on PDS map)
      // if (fPDSchannelStatus[wfChannel]) // faster than checking if the channel is in the list of selected channels
      // {

        // move Raw wvf to float
      const size_t wfsize=wf.Waveform().size();
      std::vector<float> wave;
      wave.reserve(wfsize);
      wave.assign(wf.Waveform().begin(), wf.Waveform().end());
      std::vector<float> wf1;
      wf1.assign(wave.begin(), wave.begin()+0.16 * wave.size());

      if (wfsize < 700) {
        // if (fDEBUG) std::cout << "SimpleROIAlg::ProcessWaveform: waveform size is too small"<<std::endl;
        continue; // Skip waveforms that are too small
      }

      // std::vector<float> wf2;
      // wf2.assign(wave.begin()+38, wave.begin()+76);
      double baseline1 = std::accumulate(wf1.begin(), wf1.end(), 0.0f) / wf1.size();
      double variance1 = 0.0;
      for (const auto& value : wf1) {
        variance1 += (value - baseline1) * (value - baseline1);
      }
      double baselinestd1=std::sqrt(variance1/wf1.size());
      if (baselinestd1 > 1.0){
        // if (fDEBUG) std::cout << "SimpleROIAlg::ProcessWaveform: baseline std is too high"<<std::endl;
        waveform_wrong_baseline_channel.push_back(wfChannel);
        //std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << " Baseline STD " << baselinestd1 <<std::endl;
        // event_display_fft(wave, wfChannel);
        continue; // Skip noisy waveforms
      }
      bool single_pretrigger = SingleinPretrigger(wave);
      if (!single_pretrigger){
        // if (fDEBUG) std::cout << "SimpleROIAlg::ProcessWaveform: multiple signals in pretrigger"<<std::endl;
        waveform_wrong_baseline_channel.push_back(wfChannel);
        //std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << " Baseline STD " << baselinestd1 <<std::endl;
        // event_display_fft(wave, wfChannel);
        continue; // Skip noisy waveforms
      }
      // if (baselinestd1 > 0.7 && baselinestd1 <= 10.0){
      //   // if (fDEBUG) std::cout << "SimpleROIAlg::ProcessWaveform: baseline std is too high"<<std::endl;
      //   // waveform_wrong_baseline_channel.push_back(wfChannel);
      //   std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << " Baseline STD " << baselinestd1 <<std::endl;
      //   event_display_fft(wave, wfChannel);
      //   continue; // Skip noisy waveforms
      // }
      // float max_waveform = *std::max_element(wave.begin(), wave.end());
      std::vector<float> wave_wo_baseline(wave.size());

      std::transform(
          wave.begin(),
          wave.end(),
          wave_wo_baseline.begin(),
          [baseline1](float x) { return x - baseline1; }
      );

      //find how many samples are 20 bins over the baseline
      // int samples_over_baseline = 0;
      // for (size_t i = 0; i < wave_wo_baseline.size(); ++i) {
      //   if (wave_wo_baseline[i] > 2) {
      //     samples_over_baseline++;
      //   }
      // }

      // if (samples_over_baseline > 6)
      // {
      //   // number_of_signals++;
      // } 
      // Compute the charge of the full waveform
      float charge_full_waveform = std::accumulate(wave_wo_baseline.begin(), wave_wo_baseline.end(), 0.0f);

      
      if (charge_full_waveform < -1000 || charge_full_waveform > 1e7) continue; // Skip unphysical charges

      double max_waveform = *std::max_element(wave_wo_baseline.begin(), wave_wo_baseline.end());
      //double min_waveform = *std::min_element(wave_wo_baseline.begin(), wave_wo_baseline.end());
      double max_position = std::distance(wave_wo_baseline.begin(), std::max_element(wave_wo_baseline.begin(), wave_wo_baseline.end()));
      bool real_wvfm = false;
      int samples_over_bsln = 0;
      for (size_t i = (max_position - 5); i < (max_position + 20); ++i) {
        if (wave_wo_baseline[i] > 1.5) {
          samples_over_bsln++;
        }
      }
      if (samples_over_bsln > 15) real_wvfm = true;
      // int number_signals = 0;
      // if (wfChannel == 718 && max_waveform > 600) continue; // Skip channel 718 if max_waveform is too high (noisy channel) Replace for other condition
      if (max_waveform > 4 && max_waveform < 3500 && real_wvfm) {
          // Si max_waveform está en el rango válido y real_wvfm es verdadero
          // if (fDEBUG && wfChannel == 726 && min_waveform > -4) {
          //     // Mostrar la forma de onda en el dominio de la frecuencia (modo depuración)
          //     event_display_fft(wave_wo_baseline, wfChannel);
          // }
          // Incrementar el contador de señales detectadas
          number_of_signals++;
        // auto it = std::find(fSelectedChannels.begin(), fSelectedChannels.end(), wfChannel);
        // if (it != fSelectedChannels.end()) {
        //   int index = std::distance(fSelectedChannels.begin(), it);
        //   if (static_cast<std::vector<std::vector<float>>::size_type>(index) < Baseline_SelectedChannels.size()) {
        //     //Select only the events with charge around the 1P_charge_channel
        //     // if (charge > max_charge_channel[index] -10 && charge < max_charge_channel[index]+10) {
        //     AverageWaveform_highPEs_SelectedChannels[index].push_back(wave_wo_baseline);
        //     Charge_highPEs_SelectedChannels[index].push_back(charge_full_waveform);
        //     Amplitude_highPEs_SelectedChannels[index].push_back(max_waveform);
        //     TimeStamp_SelectedChannels[index].push_back(wfStartTime);
        //     // TimeStamp_SelectedChannels[index].push_back(e.id().event());
        //   }
        //   else {
        //     std::cerr << "Error: index " << index << " out of range for containers." << std::endl;
        //   }
        }

        auto it = std::find(fSelectedChannels.begin(), fSelectedChannels.end(), wfChannel);
        if (it != fSelectedChannels.end()) {
        int index = std::distance(fSelectedChannels.begin(), it);
        if (static_cast<std::vector<std::vector<float>>::size_type>(index) < Baseline_SelectedChannels.size()) {
          //Select only the events with charge around the 1P_charge_channel
          // if (charge > max_charge_channel[index] -10 && charge < max_charge_channel[index]+10) {
          AverageWaveform_highPEs_SelectedChannels[index].push_back(wave_wo_baseline);
          Charge_highPEs_SelectedChannels[index].push_back(charge_full_waveform);
          Event_number_SelectedChannels[index].push_back(e.id().event());
          Amplitude_highPEs_SelectedChannels[index].push_back(max_waveform);
          TimeStamp_SelectedChannels[index].push_back(wfStartTime);
          // TimeStamp_SelectedChannels[index].push_back(e.id().event());
        }
        else {
          std::cerr << "Error: index " << index << " out of range for containers." << std::endl;
        }
      }
      // if (fDEBUG) {
      //     std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << std::endl;

      //   if (max_waveform < 4000 && max_waveform - baseline1 > 20 && wfsize == 768 )// && samples_over_baseline > 20)
      //   {
      //     std::cout << "First half from 10% baseline"<< baseline1 <<"+-" << baselinestd1  << std::endl;
      //     std::cout << "Second half from 10% baseline"<< baseline2 <<"+-" << baselinestd2  << std::endl;
      //     event_display_fft(wave);
      //     event_display(wave_wo_baseline,1000,1000,baseline1,baselinestd1,1000,1000,wfChannel);
      //     event_display_fft(wave_wo_baseline);
      //   }
      // }

      //Prepare ROI container
      std::vector<SimpleROI> ROIs={};
      // Call the tool for selected channels
      fROIFinderAlgPtr->ProcessWaveform(wave, ROIs, wfChannel);
      // float max_waveform = *std::max_element(wave.begin(), wave.end());


      // if (fDEBUG && (wfChannel == 713 || wfChannel == 715))
      //   {
      //     event_display(wave,-1000,-1000,100000,0.5,10000,1000000);
      //   }
      // Loop over found ROIs
      // display waveform

      // if (fDEBUG) // && (wfChannel == 705 || wfChannel == 761)) {
      //     auto maxElementIter = *std::max_element(wave.begin(), wave.end());
      //     auto minElementIter = *std::min_element(wave.begin(), wave.end());
      //     if ((maxElementIter - minElementIter > 4) && (maxElementIter - minElementIter < 4000)) {
      //         // std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel 
      //         //           << " has a max over baseline of " << maxElementIter << std::endl;
      //         event_display_waveform(wave);
      //     }
      // }

      std::vector<float> roi;
      float charge;
      float Baseline;
      float BaselineSTD;
      int startTick;
      int endTick;
      double Peak;
      int PeakTick;
      float Pedestalcharge;
      int FirstBin;
      int LastBin;
      // std::vector<float> full_wvfm;
      // std::vector<float> derivate;
      if (ROIs.size() == 0){
        waveform_no_signals_channel.push_back(wfChannel);
      }
      for (unsigned int i=0; i<ROIs.size(); i++)
      {
        if (fFull_AverageWaveform){
          roi = ROIs[i].Waveform();
          Baseline = ROIs[i].Baseline();
          Peak = ROIs[i].Peak();
        }
        else{
          roi = ROIs[i].Waveform();
          
          //fft_roi = ROIs[i].Full_wvfm();
          // smoothed = ROIs[i].SmoothROI();
          // derivate = ROIs[i].DerivateROI();
          charge = ROIs[i].Charge();
          Baseline = ROIs[i].Baseline();
          BaselineSTD = ROIs[i].BaselineSTD();
          startTick = ROIs[i].StartTick();
          endTick = ROIs[i].EndTick();
          Peak = ROIs[i].Peak();
          PeakTick = ROIs[i].PeakTick();
          Pedestalcharge = ROIs[i].PedestalCharge();
          FirstBin = ROIs[i].FirstBin();
          LastBin = ROIs[i].LastBin();
          if (charge < 35 && fDEBUG && wfChannel < 764)
          {
            std::cout<<"CALLOS: Event "<<e.id().event()<<" Channel "<<wfChannel<<" ROI "<<i<<" Charge "<<charge<<std::endl;
            // event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick, wfChannel);
            std::cout<<"alignment_point "<<PeakTick<<std::endl;
            event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick, wfChannel, FirstBin, LastBin);
            
            //event_display_fft(wave, wfChannel);
            //event_display_fft(fft_roi);
          }
        }
        //Charge for 1PE selection OV = 6.5
        //std::vector<float> max_charge_channel_1PE = {63,59,63,64,0,0,75,69,61,77,58,59,62,61,57,55,83,78,0,80,90,80,70,58,85,85,0,70,90,90,78,76};
        //Charge for 1PE selection Gain HS equalized 65ADC x ticks
        //std::vector<float> max_charge_channel_1PE = {65,65,65,65,0,0,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65,65};
        // std::vector<float> max_charge_channel_1PE = {
        // 70,70,70,70,70,70,70,70,
        // 70,70,70,70,70,70,70,70,
        // 70,70,70,70,70,70,70,70,
        // 70,70,70,70,70,70,70,70
        // };
        std::vector<float> max_charge_channel_1PE(64, 65.0f);


        // if (fDEBUG)
        // {
        //   std::cout<<"CALLOS: Event "<<e.id().event()<<" Channel "<<wfChannel<<" ROI "<<i<<" Charge "<<charge<<std::endl;
        //   //event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick);
        //   event_display_fft(roi);
        //   //event_display_fft(fft_roi);
        // }

        // Verificar que wfChannel esté en fSelectedChannels
        auto it = std::find(fSelectedChannels.begin(), fSelectedChannels.end(), wfChannel);
        if (it != fSelectedChannels.end()) {
          int index = std::distance(fSelectedChannels.begin(), it);
          if (static_cast<std::vector<std::vector<float>>::size_type>(index) < Baseline_SelectedChannels.size()) {
            //Select only the events with charge around the 1P_charge_channel
            // if (charge > max_charge_channel[index] -10 && charge < max_charge_channel[index]+10) {
            if (fFull_AverageWaveform){
              Full_AverageWaveform_SelectedChannels[index].push_back(roi);
              Baseline_SelectedChannels[index].push_back(Baseline);
              Peak_SelectedChannels[index].push_back(Peak);
            }
            else{

              Charge_SelectedChannels[index].push_back(charge);
              Peak_SelectedChannels[index].push_back(Peak);
              PedestalCharge_SelectedChannels[index].push_back(Pedestalcharge);
              PeakTick_SelectedChannels[index].push_back(PeakTick);
              Baseline_SelectedChannels[index].push_back(Baseline);
              BaselineSTD_SelectedChannels[index].push_back(BaselineSTD);
              //std::cerr << "Debug: Before accessing AverageWaveform_SelectedChannels for index " << index << std::endl;
              AverageWaveform_fewPEs_SelectedChannels[index].push_back(roi);
              //std::cerr << "Debug: After accessing AverageWaveform_SelectedChannels" << std::endl;
            }
            if (fPEMode)
            {
              if (charge > max_charge_channel_1PE[index] -10 && charge < max_charge_channel_1PE[index]+10) {
                PE_1_AverageWaveform_SelectedChannels[index].push_back(roi);
                PE_1_Peak_SelectedChannels[index].push_back(Peak);
              }
            }
          }
          else {
            std::cerr << "Error: index " << index << " out of range for containers." << std::endl;
          }
        }
        roi.clear();
     
      }
    }

  // Deshabilitar ramas específicas
  // fTree->SetBranchStatus("EventNumber", 0); // 0 para deshabilitar
  // fTree->SetBranchStatus("SignalsNumber", 0);   
  // fTree->SetBranchStatus("WaveformsNumber", 0);
  // std::cout<<"CALLOS: Event "<<e.id().event()<<" analized "<<iWvf<<" waveforms \n";
  // Store info in containers.
  fTree->Fill();

  // Clean variables needed for next event.
  if (fFull_AverageWaveform){
    for (int i=0; i<fNSelectedChannels; i++)
    {
      Baseline_SelectedChannels[i].clear();
      Peak_SelectedChannels[i].clear();
    }
  }
  else{
    for (int i=0; i<fNSelectedChannels; i++)
    {
      Charge_SelectedChannels[i].clear();
      Peak_SelectedChannels[i].clear();
      PedestalCharge_SelectedChannels[i].clear();
      PeakTick_SelectedChannels[i].clear();
      Baseline_SelectedChannels[i].clear();
      BaselineSTD_SelectedChannels[i].clear();
      Charge_highPEs_SelectedChannels[i].clear();
      Amplitude_highPEs_SelectedChannels[i].clear();
      TimeStamp_SelectedChannels[i].clear();
      if (fPEMode)
      {
        PE_1_Peak_SelectedChannels[i].clear();
      }
    }
    waveform_wrong_baseline_channel.clear();
    waveform_no_signals_channel.clear();
  }
}
  // Clear the AverageWaveform_SelectedChannels
//   for (int i = 0; i < fNSelectedChannels; ++i) {
//     AverageWaveform_SelectedChannels[i].clear();
//   }
// }

void callos::CALLOS::endJob() 
{

  //double saturation = 2 * fMaxADC;
  // when all the roi are found, we need to average them
        //Rellenar los histogramas
  
  number_of_events_analyzed_total = number_of_events_analyzed;
  number_of_signals_total = number_of_signals;
  number_waveforms_total = number_waveforms;

  fSummaryTree->Fill(); // Guarda todas las ramas activas por defecto
  // hEventNumber->SetBinContent(number_of_events_analyzed_total,1);
  // hSignalsNumber->SetBinContent(number_of_signals_total,1);
  // hWaveformsNumber->SetBinContent(number_waveforms_total,1);

  for (int i = 0; i < fNSelectedChannels; ++i) 
  {
    std::vector<float> AverageWaveform_fewPEs(fROISamples);
    std::vector<float> Full_AverageWaveform(fROISamples);
    std::vector<float> AverageWaveform_highPEs(fROISamples);
    std::vector<float> AverageWaveform_1PE(fROISamples);
    // std::vector<float> Full_AverageWaveform(fROISamples);
    //std::cout << "Debug: Before Summing al the wavefroms" << std::endl;
    //std::cout << "Debug: AverageWaveform_SelectedChannels[i].size() = " << AverageWaveform_SelectedChannels[i].size() << std::endl;
    int wvf_number = 0;
    if (fFull_AverageWaveform)
    {
      for (size_t j = 0; j < Full_AverageWaveform_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < Full_AverageWaveform_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            Full_AverageWaveform[k] += Full_AverageWaveform_SelectedChannels[i][j][k];
            //std::cout << "Debug: AverageWaveform_SelectedChannels[i][j].size() = " << AverageWaveform_SelectedChannels[i][j].size() << std::endl;
          //}        
        }
        wvf_number++;
      }
    }
    else{
      int roi_number = 0;
      for (size_t j = 0; j < AverageWaveform_fewPEs_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < AverageWaveform_fewPEs_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            AverageWaveform_fewPEs[k] += AverageWaveform_fewPEs_SelectedChannels[i][j][k];
            //std::cout << "Debug: AverageWaveform_SelectedChannels[i][j].size() = " << AverageWaveform_SelectedChannels[i][j].size() << std::endl;
          //}        
        }
        roi_number++;
      }
      int highPEs_number = 0;
      for (size_t j = 0; j < AverageWaveform_highPEs_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < AverageWaveform_highPEs_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            AverageWaveform_highPEs[k] += AverageWaveform_highPEs_SelectedChannels[i][j][k];
            //std::cout << "Debug: AverageWaveform_SelectedChannels[i][j].size() = " << AverageWaveform_SelectedChannels[i][j].size() << std::endl;
          //}        
        }
        highPEs_number++;
      }
    }
    int roi_number_1PE = 0;
    if (fPEMode)
    {

      for (size_t j = 0; j < PE_1_AverageWaveform_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < PE_1_AverageWaveform_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            AverageWaveform_1PE[k] += PE_1_AverageWaveform_SelectedChannels[i][j][k];
            //std::cout << "Debug: AverageWaveform_SelectedChannels[i][j].size() = " << AverageWaveform_SelectedChannels[i][j].size() << std::endl;
          //}        
        }
        roi_number_1PE++;
      }
    }


   
    //Since data is in different files we can not average waveform per file, we have to sum all the waveforms and then divide later by the number of waveforms
    if (fFull_AverageWaveform){
      for (int k = 0; k < fROISamples; ++k)
      {
        if (abs(Full_AverageWaveform[k]) < 4000 ) //To avoid random numbers 
        {
          hFull_AverageWaveformHistograms[i]->SetBinContent(k, Full_AverageWaveform[k]);
        }
        else 
        {
          hFull_AverageWaveformHistograms[i]->SetBinContent(k, 0);
        }
        
      }
    }
    else{
      for (int k = 0; k < fROISamples; ++k)
      {
        if (abs(AverageWaveform_fewPEs[k]) < 4000) 
        {
          hAverageWaveformHistograms_fewPEs[i]->SetBinContent(k, AverageWaveform_fewPEs[k]);
        }
        else 
        {
          hAverageWaveformHistograms_fewPEs[i]->SetBinContent(k, 0);
        }
        
      }
      for (int k = 0; k < fROISamples; ++k)
      {
        // if (abs(AverageWaveform_highPEs[k]) < 4000) 
        // {
        hAverageWaveformHistograms_highPEs[i]->SetBinContent(k, AverageWaveform_highPEs[k]);
        // }
        // else 
        // {
        //   hAverageWaveformHistograms_highPEs[i]->SetBinContent(k, 0);
        // }
        
      }
    }
    //std::cout << "Number of waveforms summed for Channel " << i << ":" << roi_number  << std::endl;
    if (fPEMode)
    {
      for (int k = 0; k < fROISamples; ++k)
      {
        if (abs(AverageWaveform_1PE[k]) < 4000) 
        {
          hAverageWaveformHistograms_1PE[i]->SetBinContent(k, AverageWaveform_1PE[k]);
        }
        else 
        {
          hAverageWaveformHistograms_1PE[i]->SetBinContent(k, 0);
        }
        
      }
    }
    
  }  
  // Pending real implementation
  // Save the histograms and NTuples to the output file.
  // std::cout<<"CALLOS: Filling the tree..."<<std::endl;
  // std::cout<<"CALLOS: End of Job"<<std::endl;
  // Close the output file
  // free the memory?
}



DEFINE_ART_MODULE(callos::CALLOS)
