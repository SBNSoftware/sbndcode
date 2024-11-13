////////////////////////////////////////////////////////////////////////
// Class:       BeamRatesCalib
// Module Type: analyzer
// File:        BeamRatesCalib.cc
//
// Analyzer taking in every-beam-spill data and determines rate of flash coincidences with scanned MON thresholds and MTC/A threshold
//
// Authors: J. McLaughlin; Based off wvfAna_module.cc
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet { //OpDet means optical detector 

  class BeamRateCalib;

  class BeamRateCalib : public art::EDAnalyzer {
  public:
    explicit BeamRateCalib(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    //idk why things are doubled here
    BeamRateCalib(BeamRateCalib const &) = delete;
    BeamRateCalib(BeamRateCalib &&) = delete;
    BeamRateCalib & operator = (BeamRateCalib const &) = delete;
    BeamRateCalib & operator = (BeamRateCalib &&) = delete;

    //Standard ART loop fucntions
    // Required functions.
    void analyze(art::Event const & e) override;
    //Selected optional functions
    void beginJob() override;
    void endJob() override;

    void SaveChannelWaveforms(const raw::OpDetWaveform &wvf, int EventCounter, int Threshold);
    void SaveMonWaveforms(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, int EventCounter, int WaveformSize);
    void ConstructMonPulse(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, std::vector<int> *MonPulse, bool Saving=false);
    void analyzeTrigger(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int WaveformSize);
    void Validate_MonSim_Outputs(art::Event const & e);
    double GetWaveformMedian(raw::OpDetWaveform wvf);
    opdet::sbndPDMapAlg pdMap; //map for photon detector types
  private:
    //Data members for analysis/output
    //Start with fcl params
    int fMonWidth; //Value*10ns MON pulse width when treshold is crossed
    //TTree *fWaveformTree;
    // Declare member data here.
    std::string fInputModuleName;
    std::string fInputProcessName;
    std::string fInputInstanceName;
    std::vector<std::string> fOpDetsToPlot; //keep for now but probably don't need
    std::string opdetType;
    std::string opdetElectronics;
    TH2D* hist_PropTriggers; //Initializes as nullptr
    TH2D* hist_PropTriggers_OffBeam; //Initializes as nullptr
    TH2D* hist_PropTriggers_NoRestriction;
    TH1D* hist_TriggerThreshold; //Initializes as nullptr
    //Helpful members to have for the analsis 
    std::vector<int> MONThresholds; //ADC counts to output a MON pulse
    std::vector<int> MTCA_Thresholds; //Number of pairs needed to trigger MTCA output
    std::vector<TH1D*> hist_TopHatPlots;
    TH1D* hist_MSUM_DutyCycle;
    TH1D* hist_MSUM_Energy;
    TH1D* hist_MSUM_Trigger_Eff;
    int EventCounter;
    int fChNumber;
    int fEvNumber;
    int fMonStart;
    int fMonStop;
    int fMonStep;
    int fMTCAStart;
    int fMTCAStop;
    int fMTCAStep;
    int fFCLthreshold;
    int fBeamWindowStart;
    int fBeamWindowEnd;
    bool fSaveAllMON;
    bool fCheckTriggers;
    int hist_id;

  };

  BeamRateCalib::BeamRateCalib(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p), 
    hist_PropTriggers(nullptr), 
    hist_PropTriggers_OffBeam(nullptr),
    hist_PropTriggers_NoRestriction(nullptr),
    hist_TriggerThreshold(nullptr), 
    hist_MSUM_DutyCycle(nullptr),
    hist_MSUM_Energy(nullptr), 
    hist_MSUM_Trigger_Eff(nullptr)
    // More initializers here.
  {
    //Read in assorted fcl parameters
    fInputModuleName = p.get< std::string >("InputModule" );
    fInputProcessName = p.get< std::string >("InputProcess" );
    fInputInstanceName = p.get< std::string >("InputInstance" );
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
    fMonWidth        = p.get<int>("MonWidth");
    fMonStart        = p.get<int>("MonStart", 20);
    fMonStop        = p.get<int>("MonStop", 250);
    fMonStep        = p.get<int>("MonStep", 20);
    fMTCAStart        = p.get<int>("MTCAStart", 3);
    fMTCAStop        = p.get<int>("MTCAStop", 64);
    fMTCAStep        = p.get<int>("MTCAStep", 1);
    fBeamWindowStart = p.get<int>("BeamWindowStart", 828+680);
    fBeamWindowEnd = p.get<int>("BeamWindowEnd", 1688+680);
    fSaveAllMON = p.get<bool>("SaveAllMon", false);
    fCheckTriggers = p.get<bool>("CheckHardwareTriggers", false); //Needs MTCA LLT to be digitized which is unusual (run 15670)
    fFCLthreshold        = p.get<int>("FCLthreshold", 15);
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
    fChNumber = 0;
    fEvNumber = 0;
    hist_id=0;
  }

  void BeamRateCalib::beginJob()
  {
    //Initilize our set of thresholds
    int MonStep=fMonStep;
    int MonStart = fMonStart;
    int MonEnd = fMonStop;
    for(int i=0; (MonStart+i*MonStep)<=MonEnd; i++)
    {
        MONThresholds.push_back( MonStart+i*MonStep);
    }
    int MTCAStep = fMTCAStep;
    int MTCAStart = fMTCAStart;
    int MTCAEnd = fMTCAStop;
    for(int i=0; MTCAStart+i*MTCAStep<=MTCAEnd; i++)
    {
        MTCA_Thresholds.push_back( MTCAStart+i*MTCAStep);
    }
    EventCounter=0;
    art::ServiceHandle<art::TFileService> tfs;
    hist_PropTriggers = tfs->make<TH2D>("hist_PropTriggers", "Proportion of Beam Spill Making Light Trigger", int((MonEnd-MonStart)/MonStep), MonStart, MonEnd, int((MTCAEnd-MTCAStart)/MTCAStep), MTCAStart, MTCAEnd); //Give it a real address + bins
    hist_PropTriggers_NoRestriction = tfs->make<TH2D>("hist_PropTriggers_LLT", "Proportion of Beam Spills with light LLT", int((MonEnd-MonStart)/MonStep), MonStart, MonEnd, int((MTCAEnd-MTCAStart)/MTCAStep), MTCAStart, MTCAEnd);
    hist_PropTriggers_OffBeam = tfs->make<TH2D>("hist_PropTriggers_OffBeam", "Proportion of Shifted Beam Spill Making Light Trigger", int((MonEnd-MonStart)/MonStep), MonStart, MonEnd, int((MTCAEnd-MTCAStart)/MTCAStep), MTCAStart, MTCAEnd); //Give it a real address + bins
    for(int i=0; (MonStart+i*MonStep)<=MonEnd; i++ )
    {
      for(int j=0; MTCAStart+j*MTCAStep<=MTCAEnd; j++)
      {
        std::string histName = "hist_TopHat_Mon_" + std::to_string(MonStart+i*MonStep) + "_MTCA_" + std::to_string(MTCAStart+j*MTCAStep);
        std::string histTitle = "hist_TopHat_Mon:" + std::to_string(MonStart+i*MonStep) + "_MTCA:" + std::to_string(MTCAStart+j*MTCAStep);
        auto tempHist = tfs->make<TH1D>(histName.c_str(), histTitle.c_str(), 5000, 0, 4999);
        hist_TopHatPlots.push_back(tempHist);
        int Index = j + i*int((MTCAEnd-MTCAStart)/MTCAStep);
        hist_TopHatPlots[Index]->GetXaxis()->SetTitle("Virtual Flash Trigger Index");
      }
    }
    if(fCheckTriggers)
    {
      //Make histgoram of MaxPassed for events with a hardware trigger. See if we can pull out the threshold from it
      std::string histname = "PMT_Passed_Per_Trigger";
      hist_TriggerThreshold = tfs->make<TH1D>(histname.c_str(), histname.c_str(), 121, -0.5, 120.5);
      std::string histname2 = "MSUM_Dutycyle";
      hist_MSUM_DutyCycle = tfs->make<TH1D>(histname2.c_str(),histname2.c_str(), 201, -0.25, 100.25 );
      std::string histname3 = "MSUM_Energy_Joulse";
      hist_MSUM_Energy = tfs->make<TH1D>(histname3.c_str(),histname3.c_str(), 100, 0, 0.000005 );
      std::string histname4 = "MSUM_Trigger_Eff_";
      histname4 = histname4 + std::to_string(fFCLthreshold);
      hist_MSUM_Trigger_Eff = tfs->make<TH1D>(histname4.c_str(),histname4.c_str(), int((MTCAEnd-MTCAStart)/MTCAStep), MTCAStart, MTCAEnd );
    }
  }

  void BeamRateCalib::analyze(art::Event const & e)
  {
    // Implementation of required member function here.
    art::ServiceHandle<art::TFileService> tfs; //Common art service should read about
    fEvNumber = e.id().event();

    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
    e.getByLabel(fInputModuleName, fInputInstanceName, fInputProcessName, waveHandle);
    if(!waveHandle.isValid() || waveHandle->size() == 0){
    std::cout << "Bad wavehandle found on event " << fEvNumber << "  " << waveHandle.isValid()  << "  " << waveHandle->size() << std::endl;
    return;
    }
    
    int WaveformSize=0;
    for(auto const& wvf : (*waveHandle)) 
    { //Iterate through each raw:OpDetWaveforms
      fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
      if(fChNumber>899) continue; //copy of trigger in signal doesn't need to be analyzed. 899 is for commissioning test
      opdetType = pdMap.pdType(fChNumber); //opdet::sbndPDMapAlg pdMap; //map for photon detector types. Lists arapuca vs PMT by channel num I guess
      opdetElectronics = pdMap.electronicsType(fChNumber); //? Maybe what kind of CAEN it is for digitization rate
      if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
      WaveformSize = wvf.size(); //Get the channel number for this waveform
      break;
    }
    //Seems to be a vector with some extra stuff like isValid()
    //raw::OpDetWaveform is a class with start time, vector of samples, and channel ID
    analyzeTrigger(waveHandle, WaveformSize); //Update the TH2D in other loop...Although maybe we dont need to 
    if(EventCounter<20 || fSaveAllMON)
    {
        //SaveChannelWaveforms(waveHandle, fEvNumber);
        SaveMonWaveforms(waveHandle, fFCLthreshold, fEvNumber, WaveformSize);
    }
    if(fCheckTriggers)
    {
      Validate_MonSim_Outputs(e);
    }
    EventCounter=EventCounter+1;
  }

  void BeamRateCalib::SaveMonWaveforms(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, int EventCounter, int WaveformSize)
  {
    std::vector<int> *MonPulse = new std::vector<int>(WaveformSize); //Maybe works? Should really add a waveform size getter 
    if(!fSaveAllMON) ConstructMonPulse(waveHandle, MonThreshold, MonPulse, true);
    else ConstructMonPulse(waveHandle, MonThreshold, MonPulse, false); //Dont save individual channels for now
    std::stringstream histname;
    histname << "event_" << EventCounter <<"_Mon"<<"_"<<hist_id;
    hist_id=hist_id+1;
    art::ServiceHandle<art::TFileService> tfs;
    TH1D *MonHist = tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(), 
                                    MonPulse->size(), 0.0, MonPulse->size()-1); //so this just breaks
    for(unsigned int i = 0; i < MonPulse->size(); i++) {
        MonHist->SetBinContent(i + 1, (double)(*MonPulse)[i]); //Loop over waveform and set bin content
      }
    delete MonPulse;
  }

  void BeamRateCalib::SaveChannelWaveforms(const raw::OpDetWaveform &wvf, int EventCounter, int Threshold)
  {
    std::stringstream histname;
    art::ServiceHandle<art::TFileService> tfs;
      int fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
      histname.str(std::string()); //Resets string stream to nothing
      histname << "event_" << EventCounter
               << "_opchannel_" << fChNumber
               << "_" << opdetType;
      //Create a new histogram
      //Make TH1D to hold waveform. String name is gross, right size and start and end time match timestamps
      TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), wvf.size(), 0, wvf.size()-1);
      for(unsigned int i = 0; i < wvf.size(); i++) {
        wvfHist->SetBinContent(i + 1, (double)wvf[i]); //Loop over waveform and set bin content
      }
      histname.str(std::string()); //reset name to name graph
      histname << "event_" << EventCounter
               << "_opchannel_" << fChNumber
               << "_" << opdetType << "_g";
      TH1D *ThreshHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), 3, 0, wvf.size()-1);
      for(unsigned int i = 0; i < 3; i++) {
        ThreshHist->SetBinContent(i + 1, (double)Threshold); //Loop over waveform and set bin content
      }
        
  }

  void BeamRateCalib::ConstructMonPulse(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, std::vector<int> *MonPulse, bool Saving)
  {
    std::vector<int> Crossings;
    std::fill(MonPulse->begin(), MonPulse->end(), 0);
    //Loop over the entries in our waveform vector
    //We care about getting the pairing correct
    int PMTIndex=0;
    std::vector<int> Pair2= { 6,   8,  10,  12,  14,  16,  36,  38,  40,  60,  62,  66,  68, 70,  84,  86,  88,  90,  92,  94, 114, 116, 118, 138, 140, 144, 146, 148, 162, 164, 168, 170, 172, 192, 194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304};
    std::vector<int> Pair1= { 7,   9,  11,  13,  15,  17,  37,  39,  41,  61,  63,  67,  69, 71,  85,  87,  89,  91,  93,  95, 115, 117, 119, 139, 141, 145, 147, 149, 163, 165, 169, 171, 173, 193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 247, 249, 251, 271, 273, 275, 295, 297, 299, 301, 303, 305};
    std::vector<int> Unpaired= {65,  64, 143, 142, 167, 166, 245, 244};
    for(auto const& wvf : (*waveHandle))
    {
        PMTIndex=PMTIndex+1; // Just counting wires in CAEN boards. Every 16th waveform is skippable, but I think the decoder handles that for me
        bool FirstChannelCross=true;
        bool SecondInPair=false;
        //Second pairs are one index below the previous value
        //Always even 
        //Not the 15th pair 
        //Add in a catch for the next CAEN start too
        
        SecondInPair = std::any_of(Pair2.begin(), Pair2.end(), [wvf](int x){return ( int(wvf.ChannelNumber()) == x); } );
        if(!SecondInPair)
        {
            Crossings.clear(); //Reset list of crossings
        }
        fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
        //if((PMTIndex)%16==0 || fChNumber>899) continue; //copy of trigger in signal doesn't need to be analyzed. 899 is for commissioning test
        if(fChNumber>899) continue; //copy of trigger in signal doesn't need to be analyzed. 899 is for commissioning test
        opdetType = pdMap.pdType(fChNumber); //opdet::sbndPDMapAlg pdMap; //map for photon detector types. Lists arapuca vs PMT by channel num I guess
        opdetElectronics = pdMap.electronicsType(fChNumber); //? Maybe what kind of CAEN it is for digitization rate
        if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
        int k=1;
        int Baseline=0;
        //int NumPoints = 10;
        //for(int Index=0; Index<NumPoints; Index++)
        //{
        //  Baseline = Baseline + wvf[Index];
        //}
        Baseline=14250; //Confirmed its been this value since Oct 4 but probably much earlier too   //Baseline/NumPoints;
        while(k < int(MonPulse->size())-1 )
        {
            //Check for crossings
            bool CrossedThreshold = (((wvf[k-1]-Baseline)>-MonThreshold) && ((wvf[k]-Baseline)<-MonThreshold));
            if(CrossedThreshold)
            {
              if(Saving && FirstChannelCross){
                SaveChannelWaveforms(wvf, fEvNumber, Baseline-MonThreshold);
                FirstChannelCross=false; // avoid overwrite of hist
              }
              bool Overlap=false;
              bool CrossLeft = std::any_of(Crossings.begin(), Crossings.end(), [k, this](int x){return (k>(x-4*fMonWidth))&&(k<x) ; } );
              bool CrossRight = std::any_of(Crossings.begin(), Crossings.end(), [k, this](int x){return (k<(x+4*fMonWidth)&&(k>x)) ;} );
              Overlap = (CrossLeft || CrossRight);
              //Need some logic for adjacent channels and for non-overlapping pulses
              //Increment the MON Pulse
              if(Overlap)
              {
                  int StartIndex=0; 
                  int EndIndex=0;
                  if(CrossLeft) //adjust end index to be first crossing bigger than k
                  {
                    StartIndex = k;
                    EndIndex = (*std::lower_bound(Crossings.begin(), Crossings.end(), k)) - 1; //First value that is not less than given val
                    if(CrossRight)
                    {
                      StartIndex = (*std::lower_bound(Crossings.begin(), Crossings.end(), k-4*fMonWidth))+4*fMonWidth+1; // adjust start index
                    }
                  }
                  else if(CrossRight)
                  {
                    EndIndex = k+4*fMonWidth;
                    StartIndex = (*std::lower_bound(Crossings.begin(), Crossings.end(), k-4*fMonWidth))+4*fMonWidth+1;
                    if(CrossLeft)
                    {
                      EndIndex = (*std::lower_bound(Crossings.begin(), Crossings.end(), k)) - 1; // adjust end index
                    }
                  }
                  if(StartIndex>=int(MonPulse->size())) StartIndex=MonPulse->size()-1; //should never happen
                  if(EndIndex>=int(MonPulse->size())) EndIndex=MonPulse->size()-1; //May happen
                  if(StartIndex<EndIndex) std::for_each(MonPulse->begin()+StartIndex, MonPulse->begin()+EndIndex, [](int &X){return X=X+1;}); //extra check for safety                  
                  if(!SecondInPair) Crossings.push_back(k); // only need to track other channel
                  k=EndIndex+1;
              }
              else
              {
                int StartIndex = k;
                int EndIndex = k+4*fMonWidth;
                if(StartIndex>=int(MonPulse->size())) StartIndex=MonPulse->size()-1; //should never happen
                if(EndIndex>=int(MonPulse->size())) EndIndex=MonPulse->size()-1; //May happen
                if(StartIndex<EndIndex) std::for_each(MonPulse->begin()+StartIndex, MonPulse->begin()+EndIndex, [](int &X){return X=X+1;});
                if(!SecondInPair) Crossings.push_back(k);
                k=EndIndex+1;
              }
            }
            else 
            {
              k=k+1;
            }
        }
    }
  }

  void BeamRateCalib::analyzeTrigger(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int WaveformSize)
  { //Use fill and not set entries for the th2d interface
    //Define some MON pulse to do trigger logics on
    //Loop over MON thresholds
    for(int i=0; i<int(MONThresholds.size()); i++)
    {
        //Loop over the entries in our waveform vector
        //We care about getting the pairing correct
        std::vector<int> *MonPulse = new std::vector<int>(WaveformSize); 
        ConstructMonPulse(waveHandle, MONThresholds[i], MonPulse, false);
        //Constructed MON pulse now lets compare it to the MTCA requirement
        //Reduce this to range where beam acceptance is valid too
        //20% of waveform is trigger
        //Accept 1.8 - 3 us 
        // 2000 - 3800 ns  with 2 ns samples
        // 1000 - 1900 samples
        int BeamAcceptanceStart = fBeamWindowStart;
        int BeamAcceptanceEnd = fBeamWindowEnd;
        int PeakMon = *std::max_element(MonPulse->begin()+BeamAcceptanceStart, MonPulse->begin()+BeamAcceptanceEnd);
        int PeakMon_NoRange = *std::max_element(MonPulse->begin(), MonPulse->end());
        int PeakOffBeam = *std::max_element(MonPulse->end()-1-(BeamAcceptanceEnd-BeamAcceptanceStart), MonPulse->end()-1);
        //Loop over MTCA thresholds
        for(int j=0; j<int(MTCA_Thresholds.size()); j++)
        {
            int PairVal = MTCA_Thresholds[j];
            int Index = j + i*int(MTCA_Thresholds.size());
            for(int Q=1; Q<int(MonPulse->size()); Q++)
            {
              if( ((*MonPulse)[Q] >=PairVal) && ((*MonPulse)[Q-1] <PairVal)   )
              {
                hist_TopHatPlots[Index]->Fill(Q); //Remake tophat plot to look at every crossing
              }
            }
            if(PeakMon>=PairVal)
            {
                //fill histogram
                double FillX = MONThresholds[i]+float((MONThresholds[1]-MONThresholds[0])/2.);
                double FillY = MTCA_Thresholds[j]+float((MTCA_Thresholds[1]-MTCA_Thresholds[0])/2.);
                hist_PropTriggers->Fill(FillX, FillY); //Fill histogram if we trigger
            }
            if(PeakOffBeam>=PairVal)
            {
              double FillX = MONThresholds[i]+float((MONThresholds[1]-MONThresholds[0])/2.);
              double FillY = MTCA_Thresholds[j]+float((MTCA_Thresholds[1]-MTCA_Thresholds[0])/2.);
              hist_PropTriggers_OffBeam->Fill(FillX, FillY);
            }
            if(PeakMon_NoRange>=PairVal)
            {
              double FillX = MONThresholds[i]+float((MONThresholds[1]-MONThresholds[0])/2.);
              double FillY = MTCA_Thresholds[j]+float((MTCA_Thresholds[1]-MTCA_Thresholds[0])/2.);
              hist_PropTriggers_NoRestriction->Fill(FillX, FillY); //Fill histogram if we trigger
            }
            //else // out of things to fill in break for time
            //{
            //  break;
            //}
        }
        delete MonPulse;
    }

  }

//Function made to run over run 15670 where we digitized the MTCA outputs in the timing CAEN 
//Assumes we have run the pmttriggerproducer as well
//Validates MON construction in producer, and extracts MTCA threshold settings. Records into histogram
  void BeamRateCalib::Validate_MonSim_Outputs(art::Event const & e)
  {
   /* //Load in waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
    e.getByLabel(fInputModuleName, fInputInstanceName, fInputProcessName, waveHandle);
    //Load in trigger information
    art::Handle< std::vector<sbnd::comm::pmtTrigger> > TriggerHandle; //User handle for vector of OpDetWaveforms
    e.getByLabel(std::string("pmttriggerproducer"), std::string(""), std::string("pmtTriggerproducer"), waveHandle);
    //Build my MON pulse equivalent
    std::vector<int> *MonPulse = new std::vector<int>(WaveformSize); 
    ConstructMonPulse(waveHandle, MonThreshold, MonPulse, false);
    //Load up the relevant waveforms from the timing CAEN
    //Have an annoying sample offset and subsampling to handle
    //80% of samples come post-trigger and 20% pre-trigger
  */
  //sbnd::comm::pmtTrigger has a vector numPassed that says how many PMT passed MON threshold for this (MON) sample
  //Also has an a number for the max number of PMT simultaneously passing (basically max of numPassed vector)
  art::Handle< std::vector<sbnd::comm::pmtTrigger> > TriggerHandle;
  e.getByLabel(std::string("pmttriggerproducer"), std::string(""), std::string("PMTDecodeAndTriggerProd"), TriggerHandle);
  //Get channel 903 (I labeled timing CAEN to be 900+ channel ID)
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fInputModuleName, fInputInstanceName, fInputProcessName, waveHandle);
  bool GoodTrigger=false; //could be an vector when we have many flash
  bool OneChannel[] = {false, false, false};
  int count=0;
  double Energy=0;
  for(auto const& wvf : (*waveHandle))
  {
    if(wvf.ChannelNumber()<904 || wvf.ChannelNumber()>907) continue;
    //Determine if waveform has a pulse in it
    std::stringstream histname;
    art::ServiceHandle<art::TFileService> tfs;
    int k=0;
    int fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
    histname.str(std::string()); //Resets string stream to nothing
    histname << "event_" << fEvNumber
              << "_opchannel_" << fChNumber;
    //Create a new histogram
    //Make TH1D to hold waveform. String name is gross, right size and start and end time match timestamps
    TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), wvf.size(), 0, wvf.size()-1);
    for(unsigned int i = 0; i < wvf.size(); i++) {
      wvfHist->SetBinContent(i + 1, (double)wvf[i]); //Loop over waveform and set bin content
    }
    std::vector<int> MonPulse;
    while(k < int(wvf.size()))
    {
      if(wvf[k]-wvf[0] > 100 && TMath::Abs(wvf[wvf.size()-1]-wvf[0])<30) //Make sure no baseline shift 
      {
        if(wvf.ChannelNumber()==905) OneChannel[0] = true;
        if(wvf.ChannelNumber()==906) OneChannel[1] = true;
        if(wvf.ChannelNumber()==907) OneChannel[2] = true;
        //break;
      }
      if(wvf.ChannelNumber()==904 && wvf[k]-wvf[0] > 17)
      {
        count=count+1; //We also care about coincident crossings though
      }
      if(wvf.ChannelNumber()==904)
      {
        //They disappate more power
        //Approximate power through the resistor is ( int(wvf[k]-wvf[0]/30)*0.125 V )**2/50 ohm
        //Total energy would be power * 10 us
        //Could do both energy and average power
        double ADC_Per_Pair = 17;
        int wvfMedian = GetWaveformMedian( wvf);
        int NPairs = std::max(int(-1*(wvf[k]-wvfMedian)/ADC_Per_Pair), 0);
        MonPulse.push_back(NPairs);
        double Volts_Per_Pair = 0.125; // 125 mV from CAEN. Knocked down by 10x by attenuators
        double Voltage = NPairs*Volts_Per_Pair;
        double Resistance = 50; //Ohms 
        double Power = Voltage*Voltage/Resistance;
        Energy =Energy +  Power*2*TMath::Power(10, -9); //Power integrated over 2ns. Joules
      }
      k=k+1;
    }
    if(wvf.ChannelNumber()==904)
    {
      double DutyCycle = double(count)/double(wvf.size());
      hist_MSUM_DutyCycle->Fill(DutyCycle);
      hist_MSUM_Energy->Fill(Energy); //Can convert to average power by dividing by readout time
      std::cout << Energy << std::endl;
      //Make th1d of trigger rate plot of column enabled in fcl
      int BeamAcceptanceStart = 828; //380; //   //1520? 380*4
      int BeamAcceptanceEnd = 1688; //2380; //    //2380? 595*4
      int PeakMon = *std::max_element(MonPulse.begin()+BeamAcceptanceStart, MonPulse.begin()+BeamAcceptanceEnd);
      for(int j=0; j<int(MTCA_Thresholds.size()); j++)
        {
            int PairVal = MTCA_Thresholds[j];
            if(PeakMon>=PairVal)
            {
                //fill histogram
                double FillY = MTCA_Thresholds[j]+float((MTCA_Thresholds[1]-MTCA_Thresholds[0])/2.);
                hist_MSUM_Trigger_Eff->Fill(FillY); //Fill histogram if we trigger
            }
        }
    }


  }
  GoodTrigger = OneChannel[0] && OneChannel[1] && OneChannel[2];
  std::stringstream histname;
  art::ServiceHandle<art::TFileService> tfs;
  histname.str(std::string()); //Resets string stream to nothing
  histname << "event_" << fEvNumber << "pmtriggerprod";
  for(auto& Trig: (*TriggerHandle))
  {
    std::vector< int > 	numPassed = Trig.numPassed;
    TH1D *pmtTrigHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), numPassed.size(), 0, numPassed.size()-1);
    for(int i=0; i<int(numPassed.size()); i++)
    {
      pmtTrigHist->SetBinContent(i+1, (double)numPassed[i]);
    }
  }
  if(GoodTrigger)
  {
    for(auto& Trig: (*TriggerHandle))
    {
      hist_TriggerThreshold->Fill(Trig.maxPMTs);
    }
  }
  //Record how proporition of readout where MSUM is over threshold by a fixed amount for a given run
  }

  void BeamRateCalib::endJob()
  {
    float Scale = 1/float(EventCounter);
    hist_PropTriggers->Scale(Scale);
    hist_PropTriggers->GetXaxis()->SetTitle("MON-Sigma Threshold (adc)");
    hist_PropTriggers->GetYaxis()->SetTitle("MTC/A Threshold (Pairs)");
    hist_PropTriggers->SetEntries(double(EventCounter));
    hist_PropTriggers_OffBeam->Scale(Scale);
    hist_PropTriggers_OffBeam->GetXaxis()->SetTitle("MON-Sigma Threshold (adc)");
    hist_PropTriggers_OffBeam->GetYaxis()->SetTitle("MTC/A Threshold (Pairs)");
    hist_PropTriggers_OffBeam->SetEntries(double(EventCounter));
    hist_PropTriggers_NoRestriction->Scale(Scale);
    hist_PropTriggers_NoRestriction->GetXaxis()->SetTitle("MON-Sigma Threshold (adc)");
    hist_PropTriggers_NoRestriction->GetYaxis()->SetTitle("MTC/A Threshold (Pairs)");
    hist_PropTriggers_NoRestriction->SetEntries(double(EventCounter));
    if(fCheckTriggers)
    {
    hist_MSUM_Trigger_Eff->Scale(Scale);
    hist_MSUM_Trigger_Eff->GetXaxis()->SetTitle("MTC/A Threshold");
    hist_MSUM_Trigger_Eff->GetYaxis()->SetTitle("Proportion Events with HLT");
    hist_MSUM_Trigger_Eff->SetEntries(double(EventCounter));
    }
  }

  double BeamRateCalib::GetWaveformMedian(raw::OpDetWaveform wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }

    DEFINE_ART_MODULE(opdet::BeamRateCalib) // Magic line that has to be here

}