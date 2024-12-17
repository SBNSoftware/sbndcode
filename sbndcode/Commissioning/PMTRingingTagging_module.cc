////////////////////////////////////////////////////////////////////////
// Class:       PMTRingingTagging
// Plugin Type: analyzer (Unknown Unknown)
// File:        PMTRingingTagging.cc
//
// Generated at Tue Nov 19 12:48:41 2024 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
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
#include "lardata/Utilities/LArFFT.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

class PMTRingingTagging;


class PMTRingingTagging : public art::EDAnalyzer {
public:
  explicit PMTRingingTagging(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTRingingTagging(PMTRingingTagging const&) = delete;
  PMTRingingTagging(PMTRingingTagging&&) = delete;
  PMTRingingTagging& operator=(PMTRingingTagging const&) = delete;
  PMTRingingTagging& operator=(PMTRingingTagging&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void SaveChannelWaveforms(const raw::OpDetWaveform &wvf);
  void ResetTree();
  double GetWaveformMedian(raw::OpDetWaveform wvf);
  double GetWaveformMedian(std::vector<short int> wvf);
  void NoRingingFound(double MaxElement, int TreeVecCounter);
  void FFTTest();

private:
  //FCL params
  int fTaggedHistsToSave;
  int fNominalWaveformSize;
  std::string fPMTWaveformLabel;
  int histCounter;
  int fTotalCAENBoards;
  int fMaxRings;
  int PeakGrabWindow;
  double fRingingAmpCut;
  //Building up the output tree with useful information
  TTree* evtTree;
  int           tree_event;                // event number
  int           tree_run;
  int           tree_NumBadChannels; 
  std::vector<double> tree_ChBiggestPulse;
  std::vector<int> tree_ringingChannelIDs; //1 per channel
  std::vector<int> tree_ChannelNumberOfRings; //1 per channel
  std::vector<int> tree_ChID; //Record to make sure I know the ID associated with each entry
  std::vector<double> tree_PeakFFT; //1 per channel
  std::vector<std::vector<double>> tree_RingingPeakAmplitudes; //Potentially many per channel
  std::vector<std::vector<int>> tree_RingingPeakSamples;
  std::vector<double> tree_RecoDampingConstants; //1 per channel
  double fFFTAmpCut;
};


void PMTRingingTagging::ResetTree()
  {
    tree_event = -1;
    tree_run = -1;
    tree_NumBadChannels = -1;
    tree_ChBiggestPulse.clear();
    tree_ChBiggestPulse.resize(fTotalCAENBoards*15);
    tree_ringingChannelIDs.clear();
    tree_ringingChannelIDs.resize(fTotalCAENBoards*15);
    tree_ChannelNumberOfRings.clear();
    tree_ChannelNumberOfRings.resize(fTotalCAENBoards*15);
    tree_ChID.clear();
    tree_ChID.resize(fTotalCAENBoards*15);
    tree_PeakFFT.clear();
    tree_PeakFFT.resize(fTotalCAENBoards*15);
    tree_RecoDampingConstants.clear();
    tree_RecoDampingConstants.resize(fTotalCAENBoards*15);
    tree_RingingPeakAmplitudes.clear();
    tree_RingingPeakSamples.clear();
    for(int i =0; i<fTotalCAENBoards*15; i++)
    {
      std::vector<double> tempRingingAmpRecording(int(fNominalWaveformSize/1000)); 
      std::vector<int> tempRingingSampleRecording(int(fNominalWaveformSize/1000)); 
      tree_RingingPeakAmplitudes.push_back(tempRingingAmpRecording);
      tree_RingingPeakSamples.push_back(tempRingingSampleRecording);
    }
  }
PMTRingingTagging::PMTRingingTagging(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fPMTWaveformLabel = p.get< std::string >("PMTWaveformLabel" );
  fTaggedHistsToSave = p.get<int>("TaggedHitsToSave");
  fNominalWaveformSize = p.get<int>("NominalWaveformSize", 5000);
  fTotalCAENBoards = p.get<int>("totalCAENBaord", 9);
  fRingingAmpCut = p.get<double>("RingingAmpCut", 70);
  fFFTAmpCut = p.get<double>("FFTminAmp", 5000);
  PeakGrabWindow=p.get<int>("PeakGrabWindow", 200);
  histCounter=0;
  int TotalPMT=fTotalCAENBoards*15;
  int MaxRings = int(fNominalWaveformSize/1000);
  fMaxRings = MaxRings;
  //Tree of objects to save
  ResetTree(); //call at start of every event
  art::ServiceHandle<art::TFileService> tfs;
  evtTree = tfs->make<TTree>("RingingMetrics","Ringing Metrics Tree");
  evtTree->Branch("event", &tree_event);
  evtTree->Branch("run", &tree_run);
  evtTree->Branch("ChBiggestPulse",&tree_ChBiggestPulse, TotalPMT, 0);
  evtTree->Branch("RingingFoundChannelIDs",&tree_ringingChannelIDs, TotalPMT, 0);
  evtTree->Branch("ChannelNumberOfRings",&tree_ChannelNumberOfRings, TotalPMT, 0);
  evtTree->Branch("ChannelID",&tree_ChID, TotalPMT, 0);
  evtTree->Branch("ChannelPeakFFTFreq",&tree_PeakFFT, TotalPMT, 0);
  evtTree->Branch("ReconstructedDampingConstant",&tree_RecoDampingConstants, TotalPMT, 0);
  //maybe change below to a size_of * TotalPMT
  evtTree->Branch("RingingPeakAmplitudes", &tree_RingingPeakAmplitudes, TotalPMT*MaxRings, 0); //Not sure if this will actually work
  evtTree->Branch("RingingPeakSamples", &  tree_RingingPeakSamples, TotalPMT*MaxRings, 0); //Not sure if this will actually work

}


void PMTRingingTagging::FFTTest()
{
  int NumberPoints=200;
  int PeriodToTest = 100;
  double FreqToTest = 1.0/PeriodToTest;
  std::vector<int> X(NumberPoints);
  std::iota(X.begin(), X.end(), 0);
  std::vector<double> Y(NumberPoints);
  for(int i=0; i<NumberPoints; i++)
  {
    Y[i] = TMath::Sin(2*TMath::Pi()*FreqToTest*X[i]);
    std::cout << "X:"<<X[i] << " Y:"<<Y[i] << std::endl;
  }
  size_t FFT_Size = 65536*5;//1<<cont;
  art::ServiceHandle<util::LArFFT> fft;
  std::vector<TComplex> EmptyVec(FFT_Size); //May need to make this uncessarily large to prevent seg fault
  //std::vector<double> EmptyVecSizeCopy(fFreqSize);
  fft->DoFFT(Y, EmptyVec);
  for(int i=0; i<500; i++)
  {
    std::cout << "Fourier Mode " << i << " Has Amp " << TComplex::Abs(EmptyVec[i]) << std::endl;
  }

  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &NumberPoints, "R2C");
  int j;
  for(j = 1; j < NumberPoints; j *= 2) {}
    int fSize = j;
    int fFreqSize = fSize / 2 + 1;

   fftr2c->SetPoints(&(Y[0]));
   fftr2c->Transform();
   double real; 
   double imaginary;
    for(int i=0; i<fFreqSize; ++i)
    {
        fftr2c->GetPointComplex(i, real, imaginary);
        EmptyVec[i]=TComplex(real, imaginary);
        std::cout << "Second go  mode " << i << " amp " << TComplex::Abs(EmptyVec[i]) << " Real " << real << " imaginary " << imaginary << std::endl;
    }


  std::cout << " did that go well?" << std::endl;
  std::cin >> NumberPoints;
}

void PMTRingingTagging::analyze(art::Event const& e)
{
  //FFTTest();
  //grab initial inputs for tree book keeping
  art::ServiceHandle<art::TFileService> tfs; //Common art service should read about
  ResetTree();
  tree_event = e.id().event();
  tree_run = e.run();
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fPMTWaveformLabel, waveHandle);
  int PMTPerBoard = 15;
  int NumPMT = fTotalCAENBoards*PMTPerBoard; //15 baords per CAEN is fixed
  int NumFlash = (*waveHandle).size()/NumPMT; 
  //if(NumFlash>1) NumFlash=1; //assume 1 flash per event
  int TreeVecCounter=0;
  double Baseline=14250;
  for(int FlashCounter=0; FlashCounter<NumFlash; FlashCounter++)
    {
      for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
        {
            //Loop over each PMT in a board
            for(int CAENChannel=0; CAENChannel<PMTPerBoard; CAENChannel++)
            {
              int WaveIndex = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NumFlash;
              auto const& wvf = (*waveHandle)[WaveIndex];
              tree_ChID[TreeVecCounter]=wvf.ChannelNumber();
              //Grab min element of waveform and if its too small assume no ringing was fully caught
              int MinIndex = std::distance(wvf.begin(), std::min_element(wvf.begin(), wvf.end())  );
              double MinADC = wvf[MinIndex];
              if( (TMath::Abs(MinADC - Baseline) < fRingingAmpCut) || ((wvf.size()-MinIndex)<5000) )
              {
                NoRingingFound(TMath::Abs(MinADC - Baseline), TreeVecCounter);
                TreeVecCounter=TreeVecCounter+1;
                continue;
              }
              //Else do our fft for the remaining hunk of waveform
              int ShiftWindow=2000;
              if(MinIndex+ShiftWindow<int(wvf.size())) MinIndex=MinIndex+ShiftWindow;
              //size_t n = wvf.size()-MinIndex;
              //size_t cont=0;
              //while(n>0){
              //    cont++; //count for every factor of 2
              //    n=(n>>1);
              //}
              //size_t FFT_Size = 65536*5;//1<<cont;
              //std::vector<TComplex> EmptyVec(FFT_Size); //May need to make this uncessarily large to prevent seg fault
              //std::vector<double> EmptyVecSizeCopy(FFT_Size);
              //std::vector<short int> TempWaveform(wvf.size()-MinIndex);
              int wvfMedian = GetWaveformMedian( wvf);
              int EndIndex = int(wvf.size());
              if(EndIndex-MinIndex > 14000) EndIndex=14000+MinIndex;
              int NumberPoints = EndIndex-MinIndex;
              std::vector<double> TempWaveform(NumberPoints);
              TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &NumberPoints, "R2C"); //1 dimensional, size of sample, Real entrys to complex?
              int fSize, fFreqSize;
              for(fSize = 1; fSize < NumberPoints; fSize *= 2) {}
              fFreqSize  = fSize/2 + 1;
              fFreqSize=fFreqSize+0; //get compiler to be happy
              std::vector<double> EmptyVec(fFreqSize); //used to be type TComplex but I don't touch the complex value directly now
              //art::ServiceHandle<util::LArFFT> fft;
              for(int i=MinIndex; i<EndIndex; i++)
              {
                int AccumBegin = i;
                int AverageWindow=10;
                int AccumEnd=i+AverageWindow;
                if(AccumEnd>=int(wvf.size())) AccumEnd=int(wvf.size())-1;
                TempWaveform[i-MinIndex] = std::accumulate(wvf.begin()+AccumBegin,wvf.begin()+AccumEnd, 0.0)/(AccumEnd-AccumBegin) - wvfMedian;
                //Add zeros to waveform to avoid later peaks
                if(TempWaveform[i-MinIndex]<-100)
                {
                  int ZeroWidth = 20;
                  int Start = i-MinIndex-ZeroWidth;
                  int End = i-MinIndex+ZeroWidth;
                  if(Start<0) Start=0;
                  if(End>=int(TempWaveform.size())) End=int(TempWaveform.size())-1;
                  for(int j=Start; j<=End; j++) TempWaveform[j] = 0;
                }
              }
              fftr2c->SetPoints(&(TempWaveform[0]));
              fftr2c->Transform();
              double real; 
              double imaginary;
              //fft->DoFFT(TempWaveform, EmptyVec);
              //FFT done so now we can explicitly check for ringing
              for(int i=0; i<int(EmptyVec.size()); i++)
              {
                fftr2c->GetPointComplex(i, real, imaginary);
                //EmptyVec[i]=TComplex(real, imaginary);
                EmptyVec[i] = TMath::Sqrt(real*real + imaginary*imaginary);
                //EmptyVecSizeCopy[i] = TComplex::Abs(EmptyVec[i]);
              }
              int MaxFFTIndex = std::distance(EmptyVec.begin(), std::max_element(EmptyVec.begin(), EmptyVec.end()) );
              double MaxFFTAmp = EmptyVec[MaxFFTIndex];
              double WindowSize =  (EndIndex-MinIndex)*(2e-9); //seconds
              double FundamentalFreq = 1/WindowSize; // Hz
              double PeakFreq = FundamentalFreq*MaxFFTIndex; //Hz
              int NumberPeaks = (wvf.size()-MinIndex)*(2e-9)*PeakFreq-1; //We want a full cycle back to peak so don't count 1
              if( (MaxFFTAmp < fFFTAmpCut) || (NumberPeaks>=fMaxRings) || (PeakFreq==0) || NumberPeaks<3)
              {
                NoRingingFound(TMath::Abs(MinADC - wvfMedian), TreeVecCounter);
                TreeVecCounter=TreeVecCounter+1;
                continue;
              }
              //We found ringing!
              //Extract peak frequency
              tree_ringingChannelIDs[TreeVecCounter] = wvf.ChannelNumber();
              tree_ChBiggestPulse[TreeVecCounter] = TMath::Abs(MinADC - wvfMedian);
              tree_ChannelNumberOfRings[TreeVecCounter] = NumberPeaks;
              tree_PeakFFT[TreeVecCounter] = PeakFreq;
              int StepSize = (EndIndex-MinIndex)/MaxFFTIndex; //Peak to peak distance
              std::vector<int> PeakTime(NumberPeaks);
              std::vector<int> PeakYFill(NumberPeaks);
              int StartIndex = std::distance(wvf.begin(), std::max_element(wvf.begin()+MinIndex+StepSize, wvf.begin()+MinIndex+int(StepSize*1.05)) );
              tree_RingingPeakSamples[TreeVecCounter][0] = StartIndex;
              tree_RingingPeakAmplitudes[TreeVecCounter][0] =wvf[tree_RingingPeakSamples[TreeVecCounter][0]]- wvfMedian;
              PeakYFill[0] =tree_RingingPeakAmplitudes[TreeVecCounter][0];
              PeakTime[0]=tree_RingingPeakSamples[TreeVecCounter][0]*2e-9;
              int FitIndexHelper=1;
              for(int i=1; i<NumberPeaks; i++)
              {
                PeakGrabWindow=StepSize*0.05;
                int InterestingIndex = StartIndex+StepSize*(i);
                if(InterestingIndex-PeakGrabWindow<0) InterestingIndex=PeakGrabWindow;
                int EndSearch = InterestingIndex+PeakGrabWindow;
                if(EndSearch>int(wvf.size())) EndSearch=int(wvf.size())-1;
                tree_RingingPeakAmplitudes[TreeVecCounter][i] = *std::max_element(wvf.begin()+InterestingIndex-PeakGrabWindow, wvf.begin()+EndSearch) - wvfMedian;;
                int IndexGrabbed = std::distance(wvf.begin(), std::max_element(wvf.begin()+InterestingIndex-PeakGrabWindow, wvf.begin()+EndSearch));
                tree_RingingPeakSamples[TreeVecCounter][i] = IndexGrabbed;
                if( tree_RingingPeakAmplitudes[TreeVecCounter][i] > 0)
                {
                  PeakTime[FitIndexHelper]= IndexGrabbed*2e-9;//StepSize*(i+0.25)*2e-9; //should change to actual time units
                  PeakYFill[FitIndexHelper] = TMath::Log( tree_RingingPeakAmplitudes[TreeVecCounter][i] );
                  FitIndexHelper=FitIndexHelper+1;
                }
                else NumberPeaks = NumberPeaks-1;
              }
              TGraph g = TGraph(NumberPeaks, &PeakTime[0], &PeakYFill[0]);
              TF1 f1 = TF1("f1","[0]*x+[1]",0., NumberPeaks+10);
              g.Fit("f1");
              tree_RecoDampingConstants[TreeVecCounter] = -1/f1.GetParameter(0); //time constant in seconds
              if((histCounter<fTaggedHistsToSave) || 
              (PeakFreq>=100e3 && PeakFreq<=400e3 && tree_ChBiggestPulse[TreeVecCounter]> 200 && NumberPeaks>3)){ 
                SaveChannelWaveforms(wvf);
              }
              TreeVecCounter=TreeVecCounter+1;
            } //CAEN channel loop
        } //CAEN board loop
    }//flash loop
    evtTree->Fill();
}//Analyze


void PMTRingingTagging::SaveChannelWaveforms(const raw::OpDetWaveform &wvf)
{
  std::stringstream histname;
  art::ServiceHandle<art::TFileService> tfs;
  int fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
  histname.str(std::string()); //Resets string stream to nothing
  histname << "event_" << tree_event
            << "_opchannel_" << fChNumber
            << "_hist_" << histCounter;
  //Create a new histogram
  //Make TH1D to hold waveform. String name is gross, right size and start and end time match timestamps
  TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), wvf.size(), 0, wvf.size()-1);
  for(unsigned int i = 0; i < wvf.size(); i++) {
    wvfHist->SetBinContent(i + 1, (double)wvf[i]); //Loop over waveform and set bin content
  }   
  histCounter++;
}



 double PMTRingingTagging::GetWaveformMedian(raw::OpDetWaveform wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }
  double PMTRingingTagging::GetWaveformMedian(std::vector<short int> wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }


void PMTRingingTagging::NoRingingFound(double MaxElement, int TreeVecCounter)
{
  //assume no ringing. Not strictly true as could start before waveform
  tree_ChBiggestPulse[TreeVecCounter] = MaxElement;
  tree_ringingChannelIDs[TreeVecCounter] = -1;
  tree_ChannelNumberOfRings[TreeVecCounter] = 0;
  tree_RecoDampingConstants[TreeVecCounter] = -1;
  tree_PeakFFT[TreeVecCounter] = -1;
  for(int j=0; j<fMaxRings; j++)
  {
    tree_RingingPeakAmplitudes[TreeVecCounter][j]=-1;
    tree_RingingPeakSamples[TreeVecCounter][j]=-1;
  }
}

DEFINE_ART_MODULE(PMTRingingTagging)
