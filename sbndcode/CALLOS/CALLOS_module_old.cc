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

  bool fPEMode;

  bool fFull_AverageWaveform;
  // std::vector<std::string> fElectronics; //Not needed for now

  //PDS map
  //opdet::sbndPDMapAlg pdsmap;

  // Total number of channels
  //int fNPDSChannels= pdsmap.size();
  
  // Filter for the loop and mapping
  //std::vector<bool> fPDSchannelStatus;
  //std::vector<int>  fPDSchannelMap;
  // Channel number
  std::vector<int> fChNumbers;
  // Selected channels
  std::vector<int> fSelectedChannels;
  int number_of_events_analyzed = 0; // Variable to store the number of events analyzed
  // Number of selected channels
  int fNSelectedChannels;
  int fMaxADC;
  std::vector<int> fDictionary; //Need to assign a number to each channel
  // Size of the Region of Interest(ROI) of each signal, must be smaller than RawWaveforms size. rodrigoa: maybe should be a tool parameter
  int fROISamples;

  // Number of samples from start of the ROI to peak. rodrigoa: maybe should be a tool parameter
  int fStartToPeak;

  // Specific channels to be analyzed (only used if not empty)
  //std::vector<int> fSpecificChannels;

  std::vector<TH1F*> hAverageWaveformHistograms_fewPEs;

  std::vector<TH1F*> hAverageWaveformHistograms_1PE;
  std::vector<TH1F*> hAverageWaveformHistograms_2PE;
  std::vector<TH1F*> hAverageWaveformHistograms_3PE;

  std::vector<TH1F*> hFull_AverageWaveformHistograms;

  //Charge container
  std::vector<std::vector<float>> Charge_SelectedChannels;
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
  //PE Mode Avg Waveform and peak container
  std::vector<std::vector<std::vector<float>>> PE_1_AverageWaveform_SelectedChannels;
  std::vector<std::vector<float>> PE_1_Peak_SelectedChannels;
  std::vector<std::vector<std::vector<float>>> PE_2_AverageWaveform_SelectedChannels;
  std::vector<std::vector<float>> PE_2_Peak_SelectedChannels;
  std::vector<std::vector<std::vector<float>>> PE_3_AverageWaveform_SelectedChannels;
  std::vector<std::vector<float>> PE_3_Peak_SelectedChannels;


  // Tool pointers
  std::unique_ptr<callos::ROIFINDERALG> fROIFinderAlgPtr;

  // le  tree
  TTree *fTree;
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
  fPEMode = p.get<bool>("PEMode", false);
  fFull_AverageWaveform = p.get<bool>("Full_AverageWaveform", false);
  // get map info
  // for each pdtype in fPDTypes, get the channels and add them to the list of selected channels
  
  fSelectedChannels.reserve(64); //hardcoded for now
  // Only signal in 700 - 763 channels
  // for (int i = 1900; i <= 1963; i++) {
  //   fSelectedChannels.push_back(i);
  // }
    for (int i = 700; i <= 763; i++) {
    fSelectedChannels.push_back(i);
  }
  //   for (int i = 1600; i <= 1663; i++) {
  //   fSelectedChannels.push_back(i);
  // }
  //   for (int i = 1300; i <= 1363; i++) {
  //   fSelectedChannels.push_back(i);
  // }

  //for (auto pdtype : fPDTypes)
  //{
  //  std::vector<int> pdtype_channels = pdsmap.getChannelsOfType(pdtype);
  //  fSelectedChannels.insert(fSelectedChannels.end(), pdtype_channels.begin(), pdtype_channels.end());
  //}

  //if (fSpecificChannels.size()>0) fSelectedChannels = fSpecificChannels;
  //short the vector
  //std::sort(fSelectedChannels.begin(), fSelectedChannels.end());
  fNSelectedChannels = fSelectedChannels.size();
  // std::cout<<"CALLOS: Selected "<<fNSelectedChannels<<" channels of selected type"<<std::endl;

  // Prepare the filter and map
  //fPDSchannelStatus.resize(fNPDSChannels, false);
  //fPDSchannelMap   .resize(fNPDSChannels);
  //for (int i=0; i<fNSelectedChannels; i++)
  //{ 
  //  fPDSchannelStatus[fSelectedChannels[i]] = true; 
  //  fPDSchannelMap   [fSelectedChannels[i]] = i;
  //}

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
  }
  if (fPEMode){
    PE_1_AverageWaveform_SelectedChannels.resize(fNSelectedChannels);
    PE_1_Peak_SelectedChannels.resize(fNSelectedChannels);
    PE_2_AverageWaveform_SelectedChannels.resize(fNSelectedChannels);
    PE_2_Peak_SelectedChannels.resize(fNSelectedChannels);
    PE_3_AverageWaveform_SelectedChannels.resize(fNSelectedChannels);
    PE_3_Peak_SelectedChannels.resize(fNSelectedChannels);
  }
  // Initialize the ROIFinderAlg tool
  fROIFinderAlgPtr = art::make_tool<callos::ROIFINDERALG>(p.get<fhicl::ParameterSet>("ROIFinderAlg"));

  // Now that we know the number of selected channels, we can create the tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("callos","callos");

  fTree->Branch("EventNumber", &number_of_events_analyzed, "EventNumber/I");

  // art::TFileDirectory chargeDir = tfs->mkdir("Charge");
  // art::TFileDirectory baselineDir = tfs->mkdir("Baseline");
  // art::TFileDirectory baselineStdDir = tfs->mkdir("BaselineSTD");
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
      
      if (fPEMode)
      {
        hAverageWaveformHistograms_1PE.push_back(tfs->make<TH1F>(
            ("hAverageWaveform_1PE_" + std::to_string(i)).c_str(),
            ("Average Waveform Channel " + std::to_string(i)).c_str(),
            fROISamples, 0, fROISamples));
        hAverageWaveformHistograms_2PE.push_back(tfs->make<TH1F>(
            ("hAverageWaveform_2PE_" + std::to_string(i)).c_str(),
            ("Average Waveform Channel " + std::to_string(i)).c_str(),
            fROISamples, 0, fROISamples));
        hAverageWaveformHistograms_3PE.push_back(tfs->make<TH1F>(
            ("hAverageWaveform_3PE_" + std::to_string(i)).c_str(),
            ("Average Waveform Channel " + std::to_string(i)).c_str(),
            fROISamples, 0, fROISamples));
      }
    }
    for (int i=0; i<fNSelectedChannels; i++)
    {

      fTree->Branch(("Charge_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Charge_SelectedChannels[i]);
      fTree->Branch(("Peak_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Peak_SelectedChannels[i]);
      fTree->Branch(("PedestalCharge_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PedestalCharge_SelectedChannels[i]);
      if (fPEMode)
      {
        fTree->Branch(("Peak_1PE_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PE_1_Peak_SelectedChannels[i]);
        fTree->Branch(("Peak_2PE_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PE_2_Peak_SelectedChannels[i]);
        fTree->Branch(("Peak_3PE_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &PE_3_Peak_SelectedChannels[i]);
      }
      fTree->Branch(("Baseline_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &Baseline_SelectedChannels[i]);
      fTree->Branch(("BaselineSTD_ch_" + std::to_string(fSelectedChannels[i])).c_str(), &BaselineSTD_SelectedChannels[i]);

    }
  }
}

void callos::CALLOS::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Increment the number of events analyzed
  number_of_events_analyzed++;
  // Load the waveforms
  art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
  e.getByLabel(fInputLabel, wfHandle);

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
      if (fPEMode)
      {
        PE_1_Peak_SelectedChannels[i]={}; // Clear the peak container
        PE_2_Peak_SelectedChannels[i]={}; // Clear the peak container
        PE_3_Peak_SelectedChannels[i]={}; // Clear the peak container
      }

    }
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
      const int wfChannel = wf.ChannelNumber();

      //const double wfStartTime = wf.TimeStamp();
      //std::cout << "CALLOS: Event " << e.id().event() << " Channel " << wfChannel << " TimeStamp " << wfStartTime << std::endl;


      // select only relevant channels (this allows to sepparate diferent XAs, PMTs, etc based on PDS map)
      // if (fPDSchannelStatus[wfChannel]) // faster than checking if the channel is in the list of selected channels
      // {

        // move Raw wvf to float
      const size_t wfsize=wf.Waveform().size();
      std::vector<float> wave;
      wave.reserve(wfsize);
      wave.assign(wf.Waveform().begin(), wf.Waveform().end());
      // std::vector<float> wf1;
      // wf1.assign(wave.begin(), wave.begin()+120);
      // std::vector<float> wf2;
      // wf2.assign(wave.begin()+38, wave.begin()+76);
      // double baseline1 = std::accumulate(wf1.begin(), wf1.end(), 0.0f) / wf1.size();
      // double baseline2 = std::accumulate(wf2.begin(), wf2.end(), 0.0f) / wf2.size();
      // double variance1 = 0.0;
      // for (const auto& value : wf1) {
      //   variance1 += (value - baseline1) * (value - baseline1);
      // }
      // double baselinestd1=std::sqrt(variance1/wf1.size());
      // double variance2 = 0.0;
      // for (const auto& value : wf2) {
      //   variance2 += (value - baseline2) * (value - baseline2);
      // }
      // double baselinestd2=std::sqrt(variance2/wf2.size());
      // float max_waveform = *std::max_element(wave.begin(), wave.end());
      // std::vector<float> wave_wo_baseline;
      // wave_wo_baseline.reserve(wave.size());
      // // Subtract the baseline from the waveform
      // std::transform(wave.begin(), wave.end(), std::back_inserter(wave_wo_baseline), 
      //                [&](float x) { return x - baseline1; });
      // //find how many samples are 20 bins over the baseline
      // int samples_over_baseline = 0;
      // for (size_t i = 0; i < wave_wo_baseline.size(); ++i) {
      //   if (wave_wo_baseline[i] > 20) {
      //     samples_over_baseline++;
      //   }
      // }
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
      // std::vector<float> full_wvfm;
      // std::vector<float> derivate;
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
          if (fDEBUG && wfChannel < 764 )
          {
            std::cout<<"CALLOS: Event "<<e.id().event()<<" Channel "<<wfChannel<<" ROI "<<i<<" Charge "<<charge<<std::endl;
            // event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick, wfChannel);
            std::cout<<"alignment_point "<<PeakTick<<std::endl;
            event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick, wfChannel);
            
            event_display_fft(wave);
            //event_display_fft(fft_roi);
          }
        }

        // if (fDEBUG)
        // {
        //   std::cout<<"CALLOS: Event "<<e.id().event()<<" Channel "<<wfChannel<<" ROI "<<i<<" Charge "<<charge<<std::endl;
        //   //event_display(wave,startTick,endTick,Baseline,BaselineSTD, Peak, PeakTick);
        //   event_display_fft(roi);
        //   //event_display_fft(fft_roi);
        // }
        //Put max charge distribution manually to select the events with charge aroun this value
        // //Run 18320 1PE
        std::vector<float> max_charge_channel_1PE = {0, 60.0, 0, 50.0, 0, 60,
        0, 50, 0, 0, 0,
        0, 0, 70, 0, 70,
        0, 90, 0, 75, 0,
        50, 0, 50, 0, 50,
        0, 50, 0, 60, 0,
        60, 0, 80, 0, 80,
        0 , 0, 0, 60, 0,
        60, 0 , 60, 0, 60,
        0, 60, 0, 60, 0, 
        60, 0, 50, 0 , 60,
        0,60,0,60, 0, 75, 0, 75};

                //Run 18320 2PE
        std::vector<float> max_charge_channel_2PE = {0, 100, 0, 100, 0, 100,
        0, 100, 0, 0, 0,
        0, 0, 154, 0, 154,
        0, 185, 0, 162, 0,
        100, 0, 100, 0, 100,
        0, 100, 0, 100, 0,
        100, 0, 166, 0, 164,
        0 , 0, 0, 60, 0,
        100, 0 , 100, 0, 100,
        0, 100, 0, 100, 0, 
        100, 0, 100, 0 , 100,
        0,100,0,100, 0, 167, 0, 162};

                //Run 18320 3PE
        std::vector<float> max_charge_channel_3PE = {0, 200.0, 0, 200.0, 0, 200,
        0, 200, 0, 0, 0,
        0, 0, 236, 0, 229,
        0, 270, 0, 239, 0,
        200, 0, 200, 0, 200,
        0, 200, 0, 200, 0,
        200, 0, 252, 0, 251,
        0 , 0, 0, 200, 0,
        200, 0 , 200, 0, 200,
        0, 200, 0, 200, 0, 
        200, 0, 200, 0 , 200,
        0,200,0,200, 0, 250, 0, 248};


        //   //Run 18338/18337 1PE
        // std::vector<float> max_charge_channel_1PE = {0, 70.0, 0, 80.0, 0, 70,
        // 0, 70, 0, 0, 0,
        // 0, 0, 80, 0, 80,
        // 0, 70, 0, 120, 0,
        // 80, 0, 100, 0, 80,
        // 0, 100, 0, 80, 0,
        // 70, 0, 85, 0, 150,
        // 0 , 0, 0, 100, 0,
        // 110, 0 , 100, 0, 80,
        // 0, 90, 0, 80, 0, 
        // 100, 0, 70, 0 , 70,
        // 0,70,0,70, 0, 80, 0, 80};

        //         //Run 18338/18337 2PE
        // std::vector<float> max_charge_channel_2PE = {0, 70.0, 0, 80.0, 0, 70,
        // 0, 70, 0, 0, 0,
        // 0, 0, 160, 0, 160,
        // 0, 140, 0, 250, 0,
        // 80, 0, 100, 0, 80,
        // 0, 100, 0, 80, 0,
        // 70, 0, 170, 0, 280,
        // 0 , 0, 0, 100, 0,
        // 110, 0 , 100, 0, 80,
        // 0, 90, 0, 80, 0, 
        // 100, 0, 70, 0 , 70,
        // 0,70,0,70, 0, 160, 0, 160};

        //         //Run 18338/18337 3PE
        // std::vector<float> max_charge_channel_3PE = {0, 70.0, 0, 80.0, 0, 70,
        // 0, 70, 0, 0, 0,
        // 0, 0, 240, 0, 240,
        // 0, 210, 0, 400, 0,
        // 80, 0, 100, 0, 80,
        // 0, 100, 0, 80, 0,
        // 70, 0, 230, 0, 430,
        // 0 , 0, 0, 100, 0,
        // 110, 0 , 100, 0, 80,
        // 0, 90, 0, 80, 0, 
        // 100, 0, 70, 0 , 70,
        // 0,70,0,70, 0, 240, 0, 240};

        // //Run 18322/18324 1PE
        // std::vector<float> max_charge_channel_1PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 60, 0, 65,
        // 0, 40, 0, 75, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 73, 0, 73,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 80, 0, 80};



        // //Run 18322(18324) 2PE
        // std::vector<float> max_charge_channel_2PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 130, 0, 130,
        // 0, 83, 0, 150, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 142, 0, 146,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 162, 0, 152};


        // //Run 18322/18324 3PE
        // std::vector<float> max_charge_channel_3PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 200, 0, 200,
        // 0, 145, 0, 220, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 217, 0, 216,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 238, 0, 227};

        //Run 18317 1PE
        // std::vector<float> max_charge_channel_1PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 56, 0, 63,
        // 0, 40, 0, 71, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 70, 0, 71,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 77, 0, 75};


        // // Run 18317 2PE
        // std::vector<float> max_charge_channel_2PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 122, 0, 125,
        // 0, 80, 0, 135, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 130, 0, 130,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 148, 0, 138};


        // //Run 18317 3PE
        // std::vector<float> max_charge_channel_3PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 183, 0, 190,
        // 0, 120, 0, 205, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 195, 0, 195,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 215, 0, 198};

        //         //Run 18174 1PE
        // std::vector<float> max_charge_channel_1PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 57, 0, 70,
        // 0, 74, 0, 64, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 60, 0, 64,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 67, 0, 61};


        // //Run 18174 2PE
        // std::vector<float> max_charge_channel_2PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 106, 0, 117,
        // 0, 74, 0, 117, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 117, 0, 120,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 136, 0, 120};


        // //Run 18174 3PE
        // std::vector<float> max_charge_channel_3PE = {0, 60.0, 0, 50.0, 0, 60,
        // 0, 50, 0, 0, 0,
        // 0, 0, 171, 0, 164,
        // 0, 74, 0, 180, 0,
        // 50, 0, 50, 0, 50,
        // 0, 50, 0, 60, 0,
        // 60, 0, 180, 0, 182,
        // 0 , 0, 0, 60, 0,
        // 60, 0 , 60, 0, 60,
        // 0, 60, 0, 60, 0, 
        // 60, 0, 50, 0 , 60,
        // 0,60,0,60, 0, 200, 0, 186};

         // //Run 18336 1PE
        //  std::vector<float> max_charge_channel_1PE = {0, 60.0, 0, 80.0, 0, 60,
        // 0, 65, 0, 0, 0,
        // 0, 0, 80, 0, 80,
        // 0, 70, 0, 120, 0,
        // 70, 0, 80, 0, 70,
        // 0, 80, 0, 70, 0,
        // 60, 0, 90, 0, 120,
        // 0 , 0, 0, 75, 0,
        // 90, 0 , 90, 0, 90,
        // 0, 80, 0, 50, 0, 
        // 70, 0, 60, 0 , 55,
        // 0,60,0,65, 0, 75, 0, 75};

        //         //Run 18336 2PE
        // std::vector<float> max_charge_channel_2PE = {0, 60.0, 0, 80.0, 0, 60,
        // 0, 65, 0, 0, 0,
        // 0, 0, 160, 0, 160,
        // 0, 140, 0, 240, 0,
        // 80, 0, 100, 0, 80,
        // 0, 100, 0, 80, 0,
        // 70, 0, 170, 0, 245,
        // 0 , 0, 0, 100, 0,
        // 110, 0 , 100, 0, 80,
        // 0, 90, 0, 80, 0, 
        // 100, 0, 70, 0 , 70,
        // 0,70,0,70, 0, 140, 0, 140};

        //         //Run 18336 3PE
        // std::vector<float> max_charge_channel_3PE = {0, 60.0, 0, 80.0, 0, 60,
        // 0, 70, 0, 0, 0,
        // 0, 0, 240, 0, 240,
        // 0, 210, 0, 400, 0,
        // 80, 0, 100, 0, 80,
        // 0, 100, 0, 80, 0,
        // 70, 0, 230, 0, 420,
        // 0 , 0, 0, 100, 0,
        // 110, 0 , 100, 0, 80,
        // 0, 90, 0, 80, 0, 
        // 100, 0, 70, 0 , 70,
        // 0,70,0,70, 0, 240, 0, 240};

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
                PE_1_Peak_SelectedChannels[index].push_back(Peak);
                PE_1_AverageWaveform_SelectedChannels[index].push_back(roi);
              }
              if (charge > max_charge_channel_2PE[index] -10 && charge < max_charge_channel_2PE[index]+10) {
                PE_2_Peak_SelectedChannels[index].push_back(Peak);
                PE_2_AverageWaveform_SelectedChannels[index].push_back(roi);
              }
              if (charge > max_charge_channel_3PE[index] -10 && charge < max_charge_channel_3PE[index]+10) {
                PE_3_Peak_SelectedChannels[index].push_back(Peak);
                PE_3_AverageWaveform_SelectedChannels[index].push_back(roi);
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

      if (fPEMode)
      {
        PE_1_Peak_SelectedChannels[i].clear();
        PE_2_Peak_SelectedChannels[i].clear();
        PE_3_Peak_SelectedChannels[i].clear();
      }
    }
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
  for (int i = 0; i < fNSelectedChannels; ++i) 
  {
    std::vector<float> AverageWaveform_fewPEs(fROISamples);
    std::vector<float> Full_AverageWaveform(fROISamples);
    std::vector<float> AverageWaveform_1PE(fROISamples);
    std::vector<float> AverageWaveform_2PE(fROISamples);
    std::vector<float> AverageWaveform_3PE(fROISamples);
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
    }


      // int fft_roi_number = 0;
      // for (size_t j = 0; j < FFTAverageWaveform_SelectedChannels[i].size(); ++j) 
      // {
      //   for (size_t k = 0; k < FFTAverageWaveform_SelectedChannels[i][j].size(); ++k) 
      //   {
      //     //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
      //     //{
      //       Full_AverageWaveform[k] += FFTAverageWaveform_SelectedChannels[i][j][k];
      //       //std::cout << "Debug: FFTAverageWaveform_SelectedChannels[i][j].size() = " << FFTAverageWaveform_SelectedChannels[i][j].size() << std::endl;
      //     //}        
      //   }
      //   fft_roi_number++;
      // }

    int roi_number_1PE = 0;
    int roi_number_2PE = 0;
    int roi_number_3PE = 0;
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
      for (size_t j = 0; j < PE_2_AverageWaveform_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < PE_2_AverageWaveform_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            AverageWaveform_2PE[k] += PE_2_AverageWaveform_SelectedChannels[i][j][k];
            //std::cout << "Debug: AverageWaveform_SelectedChannels[i][j].size() = " << AverageWaveform_SelectedChannels[i][j].size() << std::endl;
          //}        
        }
        roi_number_2PE++;
      }
      for (size_t j = 0; j < PE_3_AverageWaveform_SelectedChannels[i].size(); ++j) 
      {
        for (size_t k = 0; k < PE_3_AverageWaveform_SelectedChannels[i][j].size(); ++k) 
        {
          //if (AverageWaveform_SelectedChannels[i][j][k] < saturation)
          //{
            AverageWaveform_3PE[k] += PE_3_AverageWaveform_SelectedChannels[i][j][k];
            
          //}
        }
        roi_number_3PE++;
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
      for (int k = 0; k < fROISamples; ++k)
      {
        if (abs(AverageWaveform_2PE[k]) < 4000) 
        {
          hAverageWaveformHistograms_2PE[i]->SetBinContent(k, AverageWaveform_2PE[k]);
        }
        else 
        {
          hAverageWaveformHistograms_2PE[i]->SetBinContent(k, 0);
        }
        
      }
      for (int k = 0; k < fROISamples; ++k)
      {
        if (abs(AverageWaveform_3PE[k]) < 4000) 
        {
          hAverageWaveformHistograms_3PE[i]->SetBinContent(k, AverageWaveform_3PE[k]);
        }
        else 
        {
          hAverageWaveformHistograms_3PE[i]->SetBinContent(k, 0);
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
