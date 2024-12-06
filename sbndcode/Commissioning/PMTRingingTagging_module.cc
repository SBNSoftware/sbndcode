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

private:
  //FCL params
  int fTaggedHistsToSave;
  int fNominalWaveformSize;
  std::string fPMTWaveformLabel;
  int histCounter;
  int fTotalCAENBoards;
  int fMaxRings;
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
    for(int i =0; i<fTotalCAENBoards*15; i++)
    {
      std::vector<double> tempRingingAmpRecording(int(fNominalWaveformSize/1000)); 
      tree_RingingPeakAmplitudes.push_back(tempRingingAmpRecording);
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

}

void PMTRingingTagging::analyze(art::Event const& e)
{
  //grab initial inputs for tree book keeping
  art::ServiceHandle<art::TFileService> tfs; //Common art service should read about
  ResetTree();
  tree_event = e.id().event();
  tree_run = e.run();
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fPMTWaveformLabel, waveHandle);
  int PMTPerBoard = 15;
  //int NumPMT = fTotalCAENBoards*PMTPerBoard; //15 baords per CAEN is fixed
  int NumFlash = 1;//(*waveHandle).size()/NumPMT; //assuming exactly 1 flash for this analysis
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
              int ShiftWindow=30;
              if(MinIndex+ShiftWindow<int(wvf.size())) MinIndex=MinIndex+ShiftWindow;
              size_t n = wvf.size()-MinIndex;
              size_t cont=0;
              while(n>0){
                  cont++; //count for every factor of 2
                  n=(n>>1);
              }
              size_t FFT_Size = 65536*5;//1<<cont;
              art::ServiceHandle<util::LArFFT> fft;
              std::vector<TComplex> EmptyVec(FFT_Size); //May need to make this uncessarily large to prevent seg fault
              std::vector<double> EmptyVecSizeCopy(FFT_Size);
              std::vector<short int> TempWaveform(wvf.size()-MinIndex);
              int wvfMedian = GetWaveformMedian( wvf);
              for(int i=MinIndex; i<int(TempWaveform.size()); i++)
              {
                TempWaveform[i-MinIndex] = wvf[i] - wvfMedian;
              }
              fft->DoFFT(TempWaveform, EmptyVec);
              //FFT done so now we can explicitly check for ringing
              for(int i=0; i<int(EmptyVec.size()); i++)
              {
                EmptyVec[i] = TComplex::Abs(EmptyVec[i]);
                EmptyVecSizeCopy[i] = TComplex::Abs(EmptyVec[i]);
              }
              int MaxFFTIndex = std::distance(EmptyVecSizeCopy.begin(), std::max_element(EmptyVecSizeCopy.begin(), EmptyVecSizeCopy.end()) );
              double MaxFFTAmp = EmptyVecSizeCopy[MaxFFTIndex];
              double WindowSize =  (wvf.size()-MinIndex)*(2e-9); //seconds
              double FundamentalFreq = 1/WindowSize; // Hz
              double PeakFreq = FundamentalFreq*MaxFFTIndex; //Hz
              int NumberPeaks = WindowSize*PeakFreq + 1; //quarter cycle to get the first peak
              if( (MaxFFTAmp < fFFTAmpCut) || (NumberPeaks>=fMaxRings) || (PeakFreq==0))
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
              int StartIndex = MinIndex-ShiftWindow;
              int StepSize = 1/PeakFreq*2e9; //Number of PMT samples corresponding to desired period
              std::vector<int> PeakTime(NumberPeaks);
              std::vector<int> PeakYFill(NumberPeaks);
              int PeakGrabWindow = 30;
              tree_RingingPeakAmplitudes[TreeVecCounter][0] = *std::max_element(wvf.begin()+int(StartIndex+StepSize*0.25)-PeakGrabWindow, wvf.begin()+int(StartIndex+StepSize*0.25)+PeakGrabWindow) - wvfMedian;
              if( (wvf[int(StartIndex+StepSize*0.25)] - wvfMedian) > 0)
              {
                PeakYFill[0] = TMath::Log(wvf[int(StartIndex+StepSize*0.25)] - wvfMedian);
                PeakTime[0]=StepSize*0.25*2e-9; //seconds
              }
              int FitIndexHelper=1;
              for(int i=1; i<NumberPeaks; i++)
              {
                int InterestingIndex = StartIndex+StepSize*(i+0.25);
                tree_RingingPeakAmplitudes[TreeVecCounter][i] = *std::max_element(wvf.begin()+InterestingIndex-PeakGrabWindow, wvf.begin()+InterestingIndex+PeakGrabWindow) - wvfMedian;
                if( tree_RingingPeakAmplitudes[TreeVecCounter][i] > 0)
                {
                  PeakTime[FitIndexHelper]=StepSize*(i+0.25)*2e-9; //should change to actual time units
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
  }
}

DEFINE_ART_MODULE(PMTRingingTagging)
