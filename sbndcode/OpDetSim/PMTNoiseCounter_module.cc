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
#include "lardata/Utilities/LArFFT.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet { //OpDet means optical detector 

  class PMTNoiseCounter;

  class PMTNoiseCounter : public art::EDAnalyzer {
  public:
    explicit PMTNoiseCounter(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    //idk why things are doubled here
    PMTNoiseCounter(PMTNoiseCounter const &) = delete;
    PMTNoiseCounter(PMTNoiseCounter &&) = delete;
    PMTNoiseCounter & operator = (PMTNoiseCounter const &) = delete;
    PMTNoiseCounter & operator = (PMTNoiseCounter &&) = delete;

    //Standard ART loop fucntions
    // Required functions.
    void analyze(art::Event const & e) override;
    //Selected optional functions
    void beginJob() override;
    void endJob() override;
    double GetWaveformMedian(raw::OpDetWaveform wvf);
    void SaveChannelWaveforms(const raw::OpDetWaveform &wvf, int EventCounter);
    double GetWaveformMedian(std::vector<short int> wvf);
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
    std::vector<int> fPMTChannels;
    //We will want to record channels and events with the noise of interest
    TH1D* BadChannelHist;
    TH1D* BadChannelsPerEvent;
    //Get Per channel RMS
    TH1D* SumChannelRMS;
    std::vector<TH1D*> PMT_ChannelRMS;
    //Give some frequency for the noise of interest
    double fInterestingFrequency;
    bool fFunnyBusiness;
    int EventCounter;
    int hist_id;
    int fNominalSize;
    int fTotalCAENBoards;
    std::string fBaselineSubName;
  
  };

  PMTNoiseCounter::PMTNoiseCounter(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p), 
    BadChannelHist(nullptr), 
    BadChannelsPerEvent(nullptr),  // ,
    SumChannelRMS(nullptr)
    // More initializers here.
  {
    //Read in assorted fcl parameters
    fInputModuleName = p.get< std::string >("InputModule" );
    fInputProcessName = p.get< std::string >("InputProcess" );
    fInputInstanceName = p.get< std::string >("InputInstance" );
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
    fInterestingFrequency        = p.get<double>("InterestingFrequency");
    fTotalCAENBoards = p.get<int>("totalCAENBaord", 9);
    fPMTChannels = p.get<std::vector<int>>("ChannelsToRMS", {17,16,71,70,41,40,15,14,69,68,13,12,67,66,65, // digitzer 0 
                                                              39,38,11,10,9,8,63,62,37,36,7,6,61,60,64, // digitizer 1 
                                                              95,94,149,148,119,118,93,92,147,146,91,90,145,144,143, // digitizer 2
                                                              117,116,89,88,87,86,141,140,115,114,85,84,139,138,142, // digitizer 3 
                                                              173,172,227,226,197,196,171,170,225,224,169,168,223,222,167, // digitizer 4
                                                              195,194,221,220,165,164,219,218,193,192,163,162,217,216,166, // digitizer 5
                                                              251,250,305,304,275,274,249,248,303,302,247,246,301,300,245, // digitizer 6
                                                              273,272,299,298,243,242,297,296,271,270,241,240,295,294,244, // digitizer 7
                                                              });
    fFunnyBusiness = p.get<bool>("FunnyBusiness", false); //Use baseline to detect noise rather than the waveform itself
    fBaselineSubName = p.get<std::string>("BaselineSubName", "");
    fNominalSize = p.get<int>("NominalSize", 5000);
    EventCounter=0;
    hist_id=0;
  }

  void PMTNoiseCounter::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    //Is there clustering of the noise?
    BadChannelHist = tfs->make<TH1D>("BadChannelHist", "Channels showing noise of interest", 501, -0.5, 500.5);
    //How often does this come up as a problem?
    BadChannelsPerEvent = tfs->make<TH1D>("BadChannelsPerEvent", "Channels per event showing noise of interest", 121, -0.5, 120.5);
    for(int i=0; i<int(fPMTChannels.size()); i++)
    {
      std::stringstream histname;
      histname << "Channel_" << fPMTChannels[i] << "_RMS_NoPulses";
      PMT_ChannelRMS.push_back(tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(), 401, 0, 40.05));
    }
    SumChannelRMS = tfs->make<TH1D>("SummedWaveformRMS", "SummedWaveformRMS", 4001, -0.05, 400.05);
  }


  void PMTNoiseCounter::analyze(art::Event const & e)
  {
    //bool GoodSumWaveform = true;
    //if(!fFunnyBusiness) GoodSumWaveform = CheckSumWave(e);
    //if(!GoodSumWaveform && !fFunnyBusiness) return;
    //Else we can extract all of our fourier peaks
    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
    art::Handle< std::vector< raw::OpDetWaveform > > BaselinSub_WaveHandle; //User handle for vector of OpDetWaveforms
    e.getByLabel(fInputModuleName, fInputInstanceName, fInputProcessName, waveHandle);
    if(fFunnyBusiness) e.getByLabel(fBaselineSubName, BaselinSub_WaveHandle);
    std::vector<int> BadChannels(0);
    double FourierAmpCut = 2000; //2000; //Updated on 9/23/24
    //Restructure to be organized by flash
    int PMTPerBoard = 15;
    int NumPMT = fTotalCAENBoards*PMTPerBoard; //15 baords per CAEN is fixed
    int NumFlash = (*waveHandle).size()/NumPMT;
    
    for(int FlashCounter=0; FlashCounter<NumFlash; FlashCounter++)
    {
      for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
        {
            //Loop over each PMT in a board
            for(int CAENChannel=0; CAENChannel<PMTPerBoard; CAENChannel++)
            {
              int WaveIndex = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NumFlash;
              auto const& wvf = (*waveHandle)[WaveIndex];
              //all the channels should be good but just throw out too high numbers for extra check
              std::string opdetType = pdMap.pdType(wvf.ChannelNumber());
              if(std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
              if(wvf.ChannelNumber()>350) continue;
              //use larsoft fft
              size_t n = wvf.size();
              if(int(n)>fNominalSize)
              {
                continue; //Breaks lar fft...........
              }
              size_t cont=0;
              while(n>0){
                  cont++; //count for every factor of 2
                  n=(n>>1);
              }
              size_t FFT_Size = 1<<cont;
              art::ServiceHandle<util::LArFFT> fft;
              std::vector<TComplex> EmptyVec(FFT_Size);
              std::vector<short int> TempWaveform = wvf.Waveform();
              int wvfMedian = GetWaveformMedian( wvf);
              if(!fFunnyBusiness)
              {
                for(int i=0; i<int(TempWaveform.size()); i++)
                {
                  TempWaveform[i] = TempWaveform[i] - wvfMedian;
                }
              }
              else //Just do FFT on baseline
              {
                int minVal = 9999;
                auto const& BaselineSub_ch_wvfm = (*BaselinSub_WaveHandle)[WaveIndex];
                for(int i=0; i<int(TempWaveform.size()); i++)
                {
                  TempWaveform[i] = TempWaveform[i] - BaselineSub_ch_wvfm[i];
                  if(BaselineSub_ch_wvfm[i]<minVal) minVal=BaselineSub_ch_wvfm[i];
                }
                //if(minVal<-100) continue; //We just look at baseline here so IDK if this is needed
                wvfMedian = GetWaveformMedian( TempWaveform);
                for(int i=0; i<int(TempWaveform.size()); i++) //Center the waveform back on zero
                {
                  TempWaveform[i] = TempWaveform[i] - wvfMedian;
                }
              }
              fft->DoFFT(TempWaveform, EmptyVec);
              double Freq = fInterestingFrequency/TMath::Power(10, 6); //1/us
              double WindowTime = fNominalSize*0.002; //us
              double NumCycles = WindowTime*Freq;
              int Index = fmod(NumCycles,1)==0?(NumCycles):(NumCycles)+1; 
              double Peak = TComplex::Abs(EmptyVec[Index]);
              bool PeakIsHighest=true;
              PeakIsHighest = Peak>((10.0)*(TComplex::Abs(EmptyVec[0])));
              if((Peak>FourierAmpCut) && PeakIsHighest)
              {
                  std::cout << "Got peak val " << Peak << "  zero val " << (TComplex::Abs(EmptyVec[0])) << std::endl;
                  BadChannels.push_back(wvf.ChannelNumber());
                  SaveChannelWaveforms(wvf, int(e.id().event()));
              }
              //Do RMS analysis
              double sum = std::accumulate(TempWaveform.begin(), TempWaveform.end(), 0.0);
              double mean = sum / TempWaveform.size();
              double sq_sum = std::inner_product(TempWaveform.begin(), TempWaveform.end(), TempWaveform.begin(), 0.0);
              double stdev = std::sqrt(sq_sum / TempWaveform.size() - mean * mean);
              //Get right index for fill
              int RMS_Index= std::distance(fPMTChannels.begin(), std::find(fPMTChannels.begin(), fPMTChannels.end(), wvf.ChannelNumber()) );
              PMT_ChannelRMS[RMS_Index]->Fill(stdev);
            }
        }
      BadChannelsPerEvent->Fill(BadChannels.size());
      for(int i=0; i<int(BadChannels.size()); i++)
      {
          BadChannelHist->Fill(BadChannels[i]);
          //Plot a particular channel waveform as an additional check
      }
      BadChannels.clear();
    }

    //for(auto const& wvf : (*waveHandle))
    //{
    //    wvfCounter=wvfCounter+1;
    //    //all the channels should be good but just throw out too high numbers for extra check
    //    std::string opdetType = pdMap.pdType(wvf.ChannelNumber());
    //    if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
    //    if(wvf.ChannelNumber()>350) continue;
    //    //use larsoft fft
    //    size_t n = wvf.size();
    //    if(int(n)>fNominalSize) continue; //Breaks lar fft...........
    //    size_t cont=0;
    //    while(n>0){
    //        cont++; //count for every factor of 2
    //        n=(n>>1);
    //    }
    //    size_t FFT_Size = 1<<cont;
    //    art::ServiceHandle<util::LArFFT> fft;
    //    std::vector<TComplex> EmptyVec(FFT_Size);
    //    std::vector<short int> TempWaveform = wvf.Waveform();
    //    int wvfMedian = GetWaveformMedian( wvf);
    //    if(!fFunnyBusiness)
    //    {
    //      for(int i=0; i<int(TempWaveform.size()); i++)
    //      {
    //        TempWaveform[i] = TempWaveform[i] - wvfMedian;
    //      }
    //    }
    //    else
    //    {
    //      int minVal = 9999;
    //      auto const& BaselineSub_ch_wvfm = (*BaselinSub_WaveHandle)[wvfCounter-1];
    //      for(int i=0; i<int(TempWaveform.size()); i++)
    //      {
    //        TempWaveform[i] = TempWaveform[i] - BaselineSub_ch_wvfm[i];
    //        if(BaselineSub_ch_wvfm[i]<minVal) minVal=BaselineSub_ch_wvfm[i];
    //      }
    //      if(minVal<-100) continue;
    //      wvfMedian = GetWaveformMedian( TempWaveform);
    //      
    //      for(int i=0; i<int(TempWaveform.size()); i++) //Center the waveform back on zero
    //      {
    //        TempWaveform[i] = TempWaveform[i] - wvfMedian;
    //      }
    //    }
    //    fft->DoFFT(TempWaveform, EmptyVec);
    //    //Read out index of interest
    //    //f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (n)   if n is even should give frequencies under normal ordering
    //    //So we would need to convert from sample frequencies to real units
    //    //n = 5000 
    //    //Take in frequency in Hz
    //    double Freq = fInterestingFrequency/TMath::Power(10, 6); //1/us
    //    double WindowTime = fNominalSize*0.002; //us
    //    double NumCycles = WindowTime*Freq;
    //    int Index = fmod(NumCycles,1)==0?(NumCycles):(NumCycles)+1; 
    //    //double Freq = fInterestingFrequency/(2*3.14159); //Probably need to do this in a smarter way
    //    //int Index = fmod(1/Freq,1)==0?(1/Freq):(1/Freq)+1;
    //    double Peak = TComplex::Abs(EmptyVec[Index]);
    //    bool PeakIsHighest=true;
    //   //for(int FreqIndex=Index-1; FreqIndex<=Index+1; FreqIndex++) //May want to narrow range a little bit
    //    //{
    //    //  if(TComplex::Abs(EmptyVec[FreqIndex]) > Peak) 
    //    //  {
    //    //    Peak = TComplex::Abs(EmptyVec[FreqIndex]);
    //    //  }
    //    //}
    //    PeakIsHighest = Peak>((3.0)*(TComplex::Abs(EmptyVec[0])));
    //    //Peak = TComplex::Abs(EmptyVec[Index]); //reset this for test
    //    if((Peak>FourierAmpCut) && PeakIsHighest)
    //    {
    //        std::cout << "Got peak val " << Peak << "  zero val " << (TComplex::Abs(EmptyVec[0])) << std::endl;
    //         BadChannels.push_back(wvf.ChannelNumber());
    //        //Save the channel waveform for later inspection
    //       SaveChannelWaveforms(wvf, int(e.id().event()));
    //    }
    //    //Do RMS analysis
    //    double sum = std::accumulate(TempWaveform.begin(), TempWaveform.end(), 0.0);
    //    double mean = sum / TempWaveform.size();
    //    double sq_sum = std::inner_product(TempWaveform.begin(), TempWaveform.end(), TempWaveform.begin(), 0.0);
    //    double stdev = std::sqrt(sq_sum / TempWaveform.size() - mean * mean);
    //    //Get right index for fill
    //    int RMS_Index= std::distance(fPMTChannels.begin(), std::find(fPMTChannels.begin(), fPMTChannels.end(), wvf.ChannelNumber()) );
    //    PMT_ChannelRMS[RMS_Index]->Fill(stdev);
    //}
    //
    ////Do same RMS for summed waveform
    //if(!fFunnyBusiness)
    // {
    //  std::vector<int> SummedWave = ConstructSumWave(e);
    //  double sum = std::accumulate(SummedWave.begin(), SummedWave.end(), 0.0);
    //  double mean = sum / SummedWave.size();
    //  double sq_sum = std::inner_product(SummedWave.begin(), SummedWave.end(), SummedWave.begin(), 0.0);
    //  double stdev = std::sqrt(sq_sum / SummedWave.size() - mean * mean);
    //  SumChannelRMS->Fill(stdev);
   // }
    ////Fill histograms of interest
    //BadChannelsPerEvent->Fill(BadChannels.size());
    //for(int i=0; i<int(BadChannels.size()); i++)
   // {
   //     BadChannelHist->Fill(BadChannels[i]);
   //     //Plot a particular channel waveform as an additional check
   // }
  }

  double PMTNoiseCounter::GetWaveformMedian(raw::OpDetWaveform wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }
  double PMTNoiseCounter::GetWaveformMedian(std::vector<short int> wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }

  void PMTNoiseCounter::SaveChannelWaveforms(const raw::OpDetWaveform &wvf, int EventCounter)
  {
    std::stringstream histname;
    art::ServiceHandle<art::TFileService> tfs;
      int fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
      histname.str(std::string()); //Resets string stream to nothing
      histname << "event_" << EventCounter
               << "_opchannel_" << fChNumber << "_" <<  pdMap.pdType(fChNumber)+"_"+hist_id;
      hist_id=hist_id+1;
      std::cout << histname.str() << std::endl;
      //Create a new histogram
      //Make TH1D to hold waveform. String name is gross, right size and start and end time match timestamps
      TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), wvf.size(), 0, wvf.size()-1);
      for(unsigned int i = 0; i < wvf.size(); i++) {
        wvfHist->SetBinContent(i + 1, (double)wvf[i]); //Loop over waveform and set bin content
      } 
  }

   void PMTNoiseCounter::endJob()
  {
  }
    DEFINE_ART_MODULE(opdet::PMTNoiseCounter) // Magic line that has to be here


}