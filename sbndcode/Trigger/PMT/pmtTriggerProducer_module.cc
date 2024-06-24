/*
 File:    pmtTriggerProducer_module.cc
 Purpose: Simulates pmt hardware trigger
 Author:  Erin Yandel (eyandel@fnal.gov)
 Date:    October 26, 2021
 Version: 1.0

 This is a module to simulate the PMT hardware trigger on simulated data.
 It takes a largeant file with optical detectors waveforms and produces a
 vector whose index is time within the trigger (beam) window and content
 is number of PMT pairs passsing the trigger threshold.

 More information can be found at:
 https://sbnsoftware.github.io/sbndcode_wiki/SBND_Trigger
 and
 sbn-docdb: 23922

 The trigger logic implemented currently:
 - Turn raw OpDetWaveforms (in ADC) into binary waveforms (0=not above threshold, 1=above threshold)
 - Downsample binary waveform by 4 (simulate 2ns rate, hardware reads every 8ns)
 - Pair waveforms -> combine 2 downsampled waveforms based on parameters set in fhicl files
 - Extend wavelength when a rising edge is seen
 - Count the number of paired waveforms above threshold for various times within the trigger window

 Input:
 - output from OpDetSim module (in particular, OpDetWaveforms)

 Output:
 - std::vector<sbnd::comm::pmtTrigger> > : PMT trigger array

 fhicl: run_pmttriggerproducer.fcl
*/

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"

// SBN/SBND includes
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

class pmtTriggerProducer: public art::EDProducer {
public:
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    // Plugins should not be copied or assigned.
   explicit pmtTriggerProducer(fhicl::ParameterSet const & p);
   pmtTriggerProducer(pmtTriggerProducer const &) = delete;
   pmtTriggerProducer(pmtTriggerProducer &&) = delete;
   pmtTriggerProducer & operator = (pmtTriggerProducer const &) = delete;
   pmtTriggerProducer & operator = (pmtTriggerProducer &&) = delete;

   // Required functions.
   void produce(art::Event & evt) override;
   void reconfigure(fhicl::ParameterSet const & p);

   opdet::sbndPDMapAlg pdMap; //map for photon detector types

private:
   // Define producer-specific functions


   // Define global variables
   int run;
   int subrun;
   int event;
   size_t fEvNumber;
   size_t fChNumber; //channel number
   double fSampling; //sampling rate
   double fStartTime; //start time (in us) of raw waveform
   double fEndTime; //end time (in us) of raw waveform

   // double fBaseline = 8000.0; //baseline ADC (set in simulation)
   std::stringstream histname; //raw waveform hist name
   std::stringstream histname2; //other hists names
   std::string opdetType; //opdet wavform's opdet type (required to be pmt_coated or pmt_uncoated)

   std::vector<int> passed_trigger; //index =time (us, triger window only), content = number of pmt pairs passed threshold
   int max_passed = 0; //maximum number of pmt pairs passing threshold at the same time within trigger window
   std::vector<int> channel_numbers = {6,7,8,9,10,11,12,13,14,15,16,17,36,37,38,39,40,41,60,61,62,63,64,65,66,67,68,69,70,71,
     84,85,86,87,88,89,90,91,92,93,94,95,114,115,116,117,118,119,138,139,140,141,142,143,144,145,146,147,148,149,
     162,163,164,165,166,167,168,169,170,171,172,173,192,193,194,195,196,197,216,217,218,219,220,221,222,223,224,225,226,227,
     240,241,242,243,244,245,246,247,248,249,250,251,270,271,272,273,274,275,294,295,296,297,298,299,300,301,302,303,304,305};
   std::vector<std::vector<char>> channel_bin_wvfs;
   std::vector<char> wvf_bin_0;
   std::vector<char> paired;
   std::vector<std::vector<char>> unpaired_wvfs;

   // List parameters for the fcl file
   std::vector<double> fThreshold = {7960.0,7976.0}; //individual pmt threshold in ADC (set in fcl, passes if ADC is LESS THAN threshold), [coated, uncoated]
   bool fIndividualThresholds;
   int fOVTHRWidth;//over-threshold width, page 40 of manual (set in fcl)
   std::vector<int> fPair1 = {6,8,10,12,14,16,36,38,40,84,86,88,90,92,94,114,116,118,138,140,142,144,146,148,162,164,166,168,170,172,192,194,196,216,218,220,222,224,226,240,242,244,246,248,250,270,272,274,294,296}; //channel numbers for first set of paired pmts (set in fcl)
   std::vector<int> fPair2 = {7,9,11,13,15,17,37,39,41,85,87,89,91,93,95,115,117,119,139,141,143,145,147,149,163,165,167,169,171,173,193,195,197,217,219,221,223,225,227,241,243,245,247,249,251,271,273,275,295,297}; //channel numbers for second set of paired pmts (set in fcl)
   std::vector<int> fUnpaired = {298,299,300,301,302,303,304,305};//channel numbers for unpired pmts (set in fcl)
   std::string fPairLogic;
   double fWindowStart; //start time (in us) of trigger window (set in fcl, 0 for beam spill)
   double fWindowEnd; //end time (in us) of trigger window (set in fcl, 1.6 for beam spill)
   std::string fInputModuleName; //opdet waveform module name (set in fcl)
   std::vector<std::string> fOpDetsToPlot = {"pmt_coated", "pmt_uncoated"}; //types of optical detetcors (e.g. "pmt_coated", "xarapuca_vuv", etc.), should only be pmt_coated and pmt_uncoated (set in fcl)
   bool fSaveHists; //save raw, binary, etc. histograms (set in fcl)
   std::vector<int> fEvHists = {1,2,3}; //if fSaveHists=true, which event hists to save? (set in fcl)
   bool fVerbose; //true=output all cout statements, false=no non-error cout statements (set in fcl)

   // services
   art::ServiceHandle<art::TFileService> tfs;

}; // class pmtTriggerProducer

pmtTriggerProducer::pmtTriggerProducer(fhicl::ParameterSet const & p)
   : EDProducer(p)
{
   produces< std::vector<sbnd::comm::pmtTrigger> >();

   auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
   fSampling = clockData.OpticalClock().Frequency(); // MHz
   this->reconfigure(p);
}

// Constructor
void pmtTriggerProducer::reconfigure(fhicl::ParameterSet const & p)
{
   // Initialize member data here
   fInputModuleName = p.get< std::string >("InputModule", "opdaq");
   fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
   fIndividualThresholds = p.get<bool>("IndividualThresholds",false);
   fThreshold       = p.get<std::vector<double> >("Threshold");
   fOVTHRWidth        = p.get<int>("OVTHRWidth",11);
   fPair1    = p.get<std::vector<int> >("Pair1");
   fPair2    = p.get<std::vector<int> >("Pair2");
   fUnpaired    = p.get<std::vector<int> >("Unpaired");
   fPairLogic   = p.get<std::string>("PairLogic","OR");
   fWindowStart = p.get<double>("WindowStart",0.0);
   fWindowEnd = p.get<double>("WindowEnd",1.8);
   fSaveHists = p.get<bool>("SaveHists",true);
   fEvHists    = p.get<std::vector<int> >("EvHists");
   fVerbose = p.get<bool>("Verbose", true);

}

void pmtTriggerProducer::produce(art::Event & e)
{
   // load event info here:

   if (fVerbose){std::cout << "Processing event " << e.id().event() << std::endl;}

   std::unique_ptr< std::vector<sbnd::comm::pmtTrigger> > pmts_passed( new std::vector<sbnd::comm::pmtTrigger>);

   fEvNumber = e.id().event();
   run = e.run();
   subrun = e.subRun();
   event = e.id().event();

   art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
   e.getByLabel(fInputModuleName, waveHandle);

   if(!waveHandle.isValid()) {
      std::cout << Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
   }

  //  // example of usage for pdMap.getCollectionWithProperty()
  //  //
  //  // define a container
  //  auto inBoxTwo = pdMap.getCollectionWithProperty("pds_box", 2);
  //  // you can cout the whole json object
  //  std::cout << "inBoxTwo:\t" << inBoxTwo << "\n";
  //  // traverse its components in a loop
  //  for (auto const &e: inBoxTwo) {
  //    std::cout << e["pd_type"] << " " << e["channel"] << ' ' << "\n";
  //  }

   // // example of usage for pdMap.getCollectionFromCondition()
   // // define a lambda function with the conditions
   // auto subsetCondition = [](auto const& e)->bool
   //   // modify conditions as you want in the curly braces below
   //   {return e["pd_type"] == "pmt_uncoated" && e["tpc"] == 0;};
   // // get the container that satisfies the conditions
   // auto uncoatedsInTPC0 = pdMap.getCollectionFromCondition(subsetCondition);
   // std::cout << "uncoatedsInTPC0.size():\t" << uncoatedsInTPC0.size() << "\n";
   // for(auto const& e:uncoatedsInTPC0){
   //   std::cout << "e:\t" << e << "\n";
   // }

  // determine whether to this event should be plotted
   int i_ev = -1;
   auto iev = std::find(fEvHists.begin(), fEvHists.end(), fEvNumber);
   if (iev != fEvHists.end() && fSaveHists){
      i_ev = iev - fEvHists.begin();
   }

   if (i_ev!=-1 && i_ev<4){if (fVerbose){std::cout << "Outputting Hists" << std::endl;}}

   max_passed = 0;

   int num_pmt_wvf = 0;
   int num_pmt_ch = 0;

   if (fVerbose){
     std::cout << "fOpDetsToPlot:\t";
     for (auto const& opdet : fOpDetsToPlot){std::cout << opdet << " ";}
     std::cout << std::endl;
   }

   //find min and max start and end times to initialize all vectors to same full waveform length
   double fMinStartTime = -1510.0;//in us
   double fMaxEndTime = 1510.0;//in us

   for(auto const& wvf : (*waveHandle)) {
     fChNumber = wvf.ChannelNumber();
     opdetType = pdMap.pdType(fChNumber);
     if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
     if (wvf.TimeStamp() < fMinStartTime){ fMinStartTime = wvf.TimeStamp(); }
     if ((double(wvf.size()) / fSampling + wvf.TimeStamp()) > fMaxEndTime){ fMaxEndTime = double(wvf.size()) / fSampling + wvf.TimeStamp();}
  }
  if (fVerbose){std::cout<<"MinStartTime: "<<fMinStartTime<<" MaxEndTime: "<<fMaxEndTime<<std::endl;}

   wvf_bin_0.reserve(int((fMaxEndTime-fMinStartTime)/(1./fSampling)));
   channel_bin_wvfs.reserve(120);
   paired.reserve(fPair1.size());
   unpaired_wvfs.reserve(fPair1.size());
   passed_trigger.reserve(int((fWindowEnd-fWindowStart)/(4./fSampling)));

  // create a vector w/ the number of entries necessary for the sampling rate
  // e.g. if sampling rate is 500 MHz, each bin has width of 0.002 us or 2 ns, vector length of ~75000
   for (double i = fMinStartTime; i<fMaxEndTime+(1./fSampling); i+=(1./fSampling)){
      wvf_bin_0.push_back(0);
   }
   for (size_t i = 0; i<120; i++){
      channel_bin_wvfs.push_back(wvf_bin_0);
   }
   for (size_t i = 0; i<fPair1.size(); i++){
      paired.push_back(0);
   }
   for (size_t i = 0; i<fPair1.size(); i++){
      unpaired_wvfs.push_back(wvf_bin_0);
   }
  // window of the beam spill, 0.0 to 1.6 us
  // e.g. if sampling rate is 500 MHz, each bin has width of 0.008 us or 8 ns
   for (double i = fWindowStart; i<fWindowEnd+(4./fSampling); i+=(4./fSampling)){
      passed_trigger.push_back(0);
   }


   if (fPair2.size()!=paired.size()){std::cout<<"Pair lists mismatched sizes!"<<std::endl;}

   size_t wvf_id = -1;
   int hist_id = -1;
   for(auto const& wvf : (*waveHandle)) {
      wvf_id++;
      hist_id++;
      fChNumber = wvf.ChannelNumber();
      opdetType = pdMap.pdType(fChNumber);
      if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
      num_pmt_wvf++;

      fStartTime = wvf.TimeStamp(); //in us
      fEndTime = double(wvf.size()) / fSampling + fStartTime; //in us

      //find threshold in ADC
      double adc_threshold = 0;
      if(fThreshold.size()>0){
        adc_threshold = fThreshold.at(0); //if fIndividualThresholds=false, coated threshold
        if (!fIndividualThresholds){
          if (opdetType == "pmt_uncoated" && fThreshold.size()>1.){
          adc_threshold = fThreshold.at(1); //if fIndividualThresholds=false, uncoated threshold
        }
      }else{
        adc_threshold = fThreshold.at(num_pmt_wvf-1); //if fIndividualThresholds=true, pmt threshold
      }
    }else{
      std::cout<<"Threshold Array Empty!"<<std::endl;
    }
    // if(fVerbose){std::cout<<"Channel "<<fChNumber<<" is "<<opdetType<<" and is using theshold "<<adc_threshold<<" ADC."<<std::endl;}

      //create binary waveform
      std::vector<char> wvf_bin;
      wvf_bin.reserve(wvf.size());
      // start histo
      if (i_ev!=-1 && i_ev<3){
         histname.str(std::string());
         histname << "event_" << fEvNumber
                 << "_opchannel_" << fChNumber
                 << "_" << opdetType
                 << "_" << hist_id
                 << "_raw";
         //Create a new histogram for binary waveform
         TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), "Raw Waveform", wvf.size(), fStartTime, fEndTime);
         wvfHist->GetXaxis()->SetTitle("t (#mus)");
         for(unsigned int i = 0; i < wvf.size(); i++) {
            wvfHist->SetBinContent(i + 1, (double)wvf[i]);
         }
      } // end histo

      if (fStartTime > fMinStartTime){
         for (double i = fStartTime-fMinStartTime; i>0.; i-=(1./fSampling)){
            wvf_bin.push_back(0);
         }
      }

      for(unsigned int i = 0; i < wvf.size(); i++) {
         if((double)wvf[i]<adc_threshold){wvf_bin.push_back(1);}else{wvf_bin.push_back(0);}
      }

      if (fEndTime < fMaxEndTime){
         for (double i = fMaxEndTime-fEndTime; i>0.; i-=(1./fSampling)){
            wvf_bin.push_back(0);
         }
      }

       //combine wavform with any other waveforms from same channel
      int i_ch = -1.;
      auto ich = std::find(channel_numbers.begin(), channel_numbers.end(), fChNumber);
      if (ich != channel_numbers.end()){
         i_ch = ich - channel_numbers.begin(); // the index for the channel of interest
      }
      // if the number of bins is mismatched
      if (channel_bin_wvfs.at(i_ch).size() < wvf_bin.size()){
	       std::cout<<"Previous Channel" << fChNumber <<" Size: "<<channel_bin_wvfs.at(i_ch).size()<<"New Channel" << fChNumber <<" Size: "<<wvf_bin.size()<<std::endl;
         for(unsigned int i = channel_bin_wvfs.at(i_ch).size(); i < wvf_bin.size(); i++) {
            channel_bin_wvfs.at(i_ch).push_back(0);
         }
      }
      for(unsigned int i = 0; i < wvf_bin.size(); i++) {
         if(channel_bin_wvfs.at(i_ch).at(i)==1 || wvf_bin[i]==1){channel_bin_wvfs.at(i_ch)[i]=1;}else{channel_bin_wvfs.at(i_ch)[i]=0;}
      }
      wvf_bin.clear();
      wvf_bin.shrink_to_fit();
      waveHandle.clear();

   }//wave handle loop

     int wvf_num = -1;

     for (auto wvf_bin : channel_bin_wvfs){ // wvf_bin is a vector with entries making up the wvf, one for every channel
       wvf_num++;
       fChNumber = channel_numbers.at(wvf_num);
       fStartTime = fMinStartTime;
       fEndTime = fMaxEndTime;

       //downscale binary waveform by 4
       std::vector<char> wvf_bin_down;
       wvf_bin_down.reserve(int(wvf_bin.size()/4));
       for(unsigned int i = 0; i < wvf_bin.size(); i++) {
         if(i%4==0){wvf_bin_down.push_back(wvf_bin[i]);}
       }

       num_pmt_ch++;

      if (i_ev!=-1 && i_ev<3){
       histname2.str(std::string());
       histname2 << "event_" << fEvNumber
                << "_opchannel_" << fChNumber
                << "_binary";
       TH1D *wvfbHist = tfs->make< TH1D >(histname2.str().c_str(), "Binary Waveform", wvf_bin.size(), fStartTime, fEndTime);
       wvfbHist->GetXaxis()->SetTitle("t (#mus)");
       for(unsigned int i = 0; i < wvf_bin.size(); i++) {
         wvfbHist->SetBinContent(i + 1, wvf_bin[i]);
       }
     }

     if (i_ev!=-1 && i_ev<3){
       histname2.str(std::string());
       histname2 << "event_" << fEvNumber
                << "_opchannel_" << fChNumber
                << "_binary_down";

       TH1D *wvfbdHist = tfs->make< TH1D >(histname2.str().c_str(), "Downsampled Binary Waveform", wvf_bin_down.size(), fStartTime, fEndTime);
       wvfbdHist->GetXaxis()->SetTitle("t (#mus)");
       for(unsigned int i = 0; i < wvf_bin_down.size(); i++) {
         wvfbdHist->SetBinContent(i + 1, wvf_bin_down[i]);
       }
     }

       bool combine = false;
       bool found = false;
       bool unpaired = false;
       size_t pair_num = -1;

       for (size_t i = 0; i < fUnpaired.size(); i++){
         if (fUnpaired.at(i) == (int)fChNumber){found=true; unpaired=true;}
       }

       if (!found){
         // if a pmt is found in pair1 vector and has passed threshold (for some reason)
         // set found == true, set pair_num = index of location in Pair1
         for (size_t i = 0; i < fPair1.size(); i++){
           if (fPair1.at(i) == (int)fChNumber && paired.at(i)==1){found=true; pair_num=i; combine=true; break;}
            // else if pmt is found in pair1 vector and has not yet passed threshold,
            // set found == true, set the corresponding channel's wvf equal to downscaled wvf in the unpaired vector
            // should be paired now
           else if (fPair1.at(i) == (int)fChNumber && paired.at(i)==0){found=true; unpaired_wvfs.at(i)=wvf_bin_down; paired.at(i)=1; break;}
         }
         if (!found){
         // if still not found in pair1, check pair2 channels
           for (size_t i = 0; i < fPair2.size(); i++){
             if (fPair2.at(i) == (int)fChNumber && paired.at(i)==1){found=true; pair_num=i; combine=true; break;}
             else if (fPair2.at(i) == (int)fChNumber && paired.at(i)==0){found=true; unpaired_wvfs.at(i)=wvf_bin_down; paired.at(i)=1; break;}
           }
         }
       }

       //pair waveforms
       if (combine || unpaired){
         std::vector<char> wvf_combine;
         wvf_combine.reserve(wvf_bin_down.size());
         if (combine){
           if (unpaired_wvfs.at(pair_num).size()!=wvf_bin_down.size()){std::cout<<"Mismatched paired waveform size"<<std::endl;}
           for(unsigned int i = 0; i < wvf_bin_down.size(); i++) {
            if (fPairLogic=="OR"){
              if(unpaired_wvfs.at(pair_num)[i]==1 || wvf_bin_down[i]==1){wvf_combine.push_back(1);}else{wvf_combine.push_back(0);}
            }else if (fPairLogic=="AND"){
              if(unpaired_wvfs.at(pair_num)[i]==1 && wvf_bin_down[i]==1){wvf_combine.push_back(1);}else{wvf_combine.push_back(0);}
            }
           }
         }else if(unpaired){
           wvf_combine = wvf_bin_down;
         }

      if (i_ev!=-1 && i_ev<3){
       histname2.str(std::string());
       if (unpaired){
         histname2 << "event_" << fEvNumber
                  << "_opchannels_" << fChNumber
                  << "_unpaired"
                  << "_combined";
       }else{
         histname2 << "event_" << fEvNumber
                  << "_opchannels_" << fPair1.at(pair_num)
                  << "_" << fPair2.at(pair_num)
                  << "_combined";
       }

       TH1D *wvfcHist = tfs->make< TH1D >(histname2.str().c_str(), "Paired Waveform", wvf_combine.size(), fStartTime, fEndTime);
       wvfcHist->GetXaxis()->SetTitle("t (#mus)");
       for(unsigned int i = 0; i < wvf_combine.size(); i++) {
         wvfcHist->SetBinContent(i + 1, wvf_combine[i]);
       }
     }

       //implement over threshold trigger signal width
       //(Every time the combined waveform transitions from 0 to 1, change the next fOVTHRWidth values to 1 (ex: fOVTHRWidth=11 -> 12 high -> 12*8=96 ns true) )
       for(unsigned int i = 1; i < wvf_combine.size()-fOVTHRWidth; i++) {
         if(wvf_combine[i]==1 && wvf_combine[i-1]==0){
           for(unsigned int j = i+1; j < i+fOVTHRWidth+1; j++){
             wvf_combine[j] = 1;
           }
           }
       }

      if (i_ev!=-1 && i_ev<3){
       histname2.str(std::string());
       if (unpaired){
         histname2 << "event_" << fEvNumber
                  << "_opchannels_" << fChNumber
                  << "_unpaired"
                  << "_combined_width";
       }else{
         histname2 << "event_" << fEvNumber
                  << "_opchannels_" << fPair1.at(pair_num)
                  << "_" << fPair2.at(pair_num)
                  << "_combined_width";
       }

       TH1D *wvfcwHist = tfs->make< TH1D >(histname2.str().c_str(), "Over Threshold Paired Waveform", wvf_combine.size(), fStartTime, fEndTime);
       wvfcwHist->GetXaxis()->SetTitle("t (#mus)");
       for(unsigned int i = 0; i < wvf_combine.size(); i++) {
         wvfcwHist->SetBinContent(i + 1, wvf_combine[i]);
       }
     }

       //Combine the waveforms to get a 1D array of integers where the value corresponds to the number of pairs ON and the
       //index corresponds to the tick in the waveform
       double binspermus = wvf_combine.size()/(fEndTime-fStartTime);
       unsigned int startbin = std::floor(binspermus*(fWindowStart - fStartTime));
       unsigned int endbin = std::ceil(binspermus*(fWindowEnd - fStartTime));
       if (startbin<0.){startbin=0;}
       if (endbin > wvf_combine.size() - 1){endbin = wvf_combine.size() - 1;}
       if (passed_trigger.size() < endbin-startbin){
         for (unsigned int i = passed_trigger.size(); i<endbin; i++){
           passed_trigger.push_back(0);
         }
       }
       unsigned int i_p = 0;
       for(unsigned int i = startbin; i<endbin; i++){
         if (wvf_combine.at(i)==1){passed_trigger.at(i_p)++;}
         i_p++;
       }

       wvf_bin_down.clear();
       wvf_bin_down.shrink_to_fit();
       wvf_combine.clear();
       wvf_combine.shrink_to_fit();

     }

   }


  if (i_ev!=-1 && i_ev<3){
   histname.str(std::string());
   histname << "event_" << fEvNumber
            << "_passed_trigger";

   TH1D *passedHist = tfs->make< TH1D >(histname.str().c_str(), "Number of PMTs Passing Trigger During Beam", passed_trigger.size(), fWindowStart, fWindowEnd+(4./fSampling));
   passedHist->GetXaxis()->SetTitle("t (#mus)");
   for(unsigned int i = 0; i < passed_trigger.size(); i++) {
     passedHist->SetBinContent(i + 1, passed_trigger[i]);
   }
 }

  sbnd::comm::pmtTrigger pmt_time;

   for (int pmts: passed_trigger){
     pmt_time.numPassed.push_back(pmts);
     if (pmts > max_passed) max_passed = pmts;
   }
   pmt_time.maxPMTs = max_passed;
   pmts_passed->push_back(pmt_time);

   if (fVerbose){std::cout << "Length of passed trigger: "  << pmt_time.numPassed.size() << std::endl;}
   if (fVerbose){std::cout << "Max number of PMTs passed: " << pmt_time.maxPMTs << std::endl;}

   // the following lines "push" the relevant products you want to produce
   // EXAMPLE:
   // evt.put(std::move(muon_tracks));
   // evt.put(std::move(muon_tracks_assn));

   e.put(std::move(pmts_passed));

   //clear variables
   passed_trigger.clear();
   passed_trigger.shrink_to_fit();
   max_passed = 0;

   if (fVerbose){std::cout << "Number of PMT waveforms: " << num_pmt_wvf << std::endl;}
   if (fVerbose){std::cout << "Number of PMT channels: " << num_pmt_ch << std::endl;}

   channel_bin_wvfs.clear();
   channel_bin_wvfs.shrink_to_fit();
   paired.clear();
   paired.shrink_to_fit();
   unpaired_wvfs.clear();
   unpaired_wvfs.shrink_to_fit();

   wvf_bin_0.clear();
   wvf_bin_0.shrink_to_fit();

} // pmtTriggerProducer::produce()

// A macro required for a JobControl module.
DEFINE_ART_MODULE(pmtTriggerProducer)