////////////////////////////////////////////////////////////////////////
// Class:       michelETagger
// Plugin Type: filter (Unknown Unknown)
// File:        michelETagger_module.cc
//
// Generated at Mon Jan 27 14:29:41 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbnobj/SBND/CRT/CRTCluster.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "TH1D.h"
#include <memory>

class michelETagger;


class michelETagger : public art::EDFilter {
public:
  explicit michelETagger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  michelETagger(michelETagger const&) = delete;
  michelETagger(michelETagger&&) = delete;
  michelETagger& operator=(michelETagger const&) = delete;
  michelETagger& operator=(michelETagger&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;
  bool DoubleFlashCheck(std::vector<double> &SummedVector);
  void ConstructSummedWaveform(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, std::vector<double> &SummedVector_TPC1, std::vector<double> &SummedVector_TPC2, int &FlashCounter);
  void ConvolveWithAnyKernel(std::vector<double> &Waveform, std::vector<double> &Kernel, std::vector<double> &Out);
  void SaveVector(std::vector<double> &HistEntries, std::string Name);



private:
  opdet::sbndPDMapAlg pdMap; //map for photon detector types
  int fGaussianConvlSize;
  double fGaussianConvlWidth;
  double fCoincidentWindow;
  int fTotalCAENBoards;
  std::vector<int> PMTPerBoard;
  std::string fCRTclusterLabel;
  std::string fPMTLabel;
  double fMuonADCCutoff;
  double fMichelADCCutoff;
  int fNsPerSample;
  int fBaseline;
  double fMuonLifetimes;
  int EventNum;
  int FlashNumForName;
  int TPCNumForName;
  // Declare member data here.

};


michelETagger::michelETagger(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGaussianConvlSize = p.get< int >("GaussianConvlSize", 501 );
  fGaussianConvlWidth = p.get< double >("GaussianConvlWidth", 160. );
  fCoincidentWindow = p.get< double >("CoincidentWindow", 250.); //ns
  fTotalCAENBoards = p.get<int>("TotalCAENBoards", 8);
  PMTPerBoard = p.get<std::vector<int>>("PMTPerBoard", {10, 13, 12, 7, 9, 10, 11, 12});
  fCRTclusterLabel = p.get<std::string>("CRTclusterLabel");
  fPMTLabel = p.get<std::string>("PMTLabel");
  fMuonADCCutoff = p.get<double>("MuonADCCutoff", 1000.);
  fMichelADCCutoff = p.get<double>("MichelADCCutoff", 100);
  fNsPerSample = p.get<int>("NsPerSample", 2);
  fBaseline = p.get<int>("Baseline", 0);
  fMuonLifetimes = p.get<double>("MuonLifetimes", 2.5);
}

bool michelETagger::DoubleFlashCheck(std::vector<double> &SummedVector)
{
  //Baseline subtracted summed waveform for a give TPC is available
  //Check for at least two big excursions from baseline
  //We have a 4 PE* 3PMT threshold so maybe calling 12PE big enough is okay? 
  //Lets make it a fcl parameter and have it default to 8PE
  //Deconvolution does not necessarily preserve amplitude but its of similar scale

  //Basic procedure will be to pick out fast light, with an alg more designed to deconvolved waveforms
  //We will convolve the Summed waveform with a guassian smoothing kernel, then an edge detection kernel
  //Finding the 0-crossing points in the edge detection kernel output gives us a sense of peaks in the smoothed waveform
  //These peaks are ~coincident with the true waveform peaks
  //We will grab the two largest peaks and make a decision on Michel-ness based on delta T and overall flash sizes
  
  //Construct the gaussian kernel
  bool DoubleFlash=false;
  std::vector<int> X(fGaussianConvlSize);
  std::iota(X.begin(), X.end(), -int(fGaussianConvlSize/2)); // Indices
  std::vector<double> GaussianKernel(fGaussianConvlSize);
  for(int i=0; i<fGaussianConvlSize; i++)
  {
    GaussianKernel[i] = (1/TMath::Sqrt(2*TMath::Pi() * TMath::Power(double(fGaussianConvlWidth), 2.0) )*
    TMath::Exp( - TMath::Power(double(X[i]), 2.0) / (2*TMath::Power(double(fGaussianConvlWidth), 2.0)) ) ); 
  }
  //Do convolution to smooth the waveform
  //Make edge detection kernel
  std::vector<double> EdgeDetectionKernel = {0, 1, 1, -1, -1, 0};
  std::vector<double> SmoothedWaveform(SummedVector.size());
  std::vector<double> EdgeWaveform(SummedVector.size());
  ConvolveWithAnyKernel(SummedVector, GaussianKernel, SmoothedWaveform);
  SaveVector(SmoothedWaveform, std::string("Event_") + std::to_string(EventNum)+std::string("_Flash_")+
    std::to_string(FlashNumForName)+std::string("TPC_")+std::to_string(TPCNumForName)+std::string("_Smoothed"));
  //Do edge detection on waveform 
  ConvolveWithAnyKernel(SmoothedWaveform, EdgeDetectionKernel, EdgeWaveform); //Summed vector passed by reference and modified
  SaveVector(EdgeWaveform, std::string("Event_") + std::to_string(EventNum)+std::string("_Flash_")+
    std::to_string(FlashNumForName)+std::string("TPC_")+std::to_string(TPCNumForName)+std::string("_EdgeWaveform"));
  //Apply selection cuts to our edge detection waveform 
  std::vector<int> CrossingIndecies;
  for(int i=1; i<int(EdgeWaveform.size()); i++)
  {
    if(EdgeWaveform[i-1]>0 && EdgeWaveform[i]<0)
    {
      CrossingIndecies.push_back(i);
    }
  }
  //Add saving for each analysis step
  //Will format histogram names as event_N_Flash_Y_Step_Name_TPC_Z
  if(int(CrossingIndecies.size())<2) return DoubleFlash; // else we have enough indecies
  //May want to add some histogram saving of summed waveforms, smoothed, edge, etc to check on algorithm
  //Now take the two largest smoothed values with a crossing index
  //Lambda function notes for future reference
  //[&] captures used variables by reference
  //(int...) function paraemters that are entries in CrossingIndecies
  //-> bool lets me define the type of the return 
  //Body compares values and returns true when something is greater than other
  //So the first entry is the one that always returns greater
  std::sort(CrossingIndecies.begin(), CrossingIndecies.end(), 
  [&](int Index_1, int Index_2)->bool { return SmoothedWaveform[Index_1] > SmoothedWaveform[Index_2]; } ); //sorts index from max to min
  //Now we can apply the actual waveform selection we want
  //Michel e- has the second largest peak follow the largest
  bool MichelFollowsMuon = (CrossingIndecies[0] < CrossingIndecies[1] );
  //Require that both muon and michel are large enough
  bool BigEnoughMuonFlash = SmoothedWaveform[CrossingIndecies[0]] >= fMuonADCCutoff;
  bool BigEnoughMichelFlash = SmoothedWaveform[CrossingIndecies[1]] >= fMichelADCCutoff;
  //Finally require the Muon lifetime is approximately correct
  double TimeToMichel = (CrossingIndecies[1]-CrossingIndecies[0])*fNsPerSample;
  double MuonLifetime = 2197; // 2197 ns mu+ lifetime;   616 mu- ns lifetime; These are all timeconstants not t1/2
  bool GoodLifetime = TimeToMichel < (fMuonLifetimes*MuonLifetime);
  DoubleFlash = MichelFollowsMuon && BigEnoughMuonFlash && BigEnoughMichelFlash && GoodLifetime;
  return DoubleFlash;
}

//This should really be done by TPC
void michelETagger::ConstructSummedWaveform(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, std::vector<double> &SummedVector_TPC1, 
std::vector<double> &SummedVector_TPC2, int &FlashCounter)
{
  int NumFlash = waveHandle->size()/(std::accumulate(PMTPerBoard.begin(), PMTPerBoard.end(), 0));
  //Loop over each CAEN
  for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
      {
          //Loop over each PMT in a board
          for(int CAENChannel=0; CAENChannel<PMTPerBoard[CurrentBoard]; CAENChannel++)
          {
            //Extract the waveforms associated with this particular flash
            //Waveforms are ordered by CAEN, Flash, PMT channel
            //At least in the decoded files
            //Have to verify for the deconvolved files, but I assume order is preserved
            int WaveIndex = (CAENChannel + 
            FlashCounter*PMTPerBoard[CurrentBoard] + (std::accumulate(PMTPerBoard.begin(), PMTPerBoard.begin()+CurrentBoard, 0))*NumFlash);
            auto const& wvf = (*waveHandle)[WaveIndex];
            for(int Sample=0; Sample<int(wvf.size()); Sample++)
            {
              int Baseline = fBaseline; //Fixed fcl baseline for use with either decoded or deconvolved waveforms
              if( pdMap.pdTPC(wvf.ChannelNumber()) == 0 ) SummedVector_TPC1[Sample] = wvf[Sample] - Baseline; //Update summed vector
              else if( pdMap.pdTPC(wvf.ChannelNumber()) == 1 ) SummedVector_TPC2[Sample] = wvf[Sample] - Baseline; //Update summed vector
            }
          }
      } // Finish loop over all waveforms in this flash and Summed vector for each TPC is ready for use in Analyze

}


void michelETagger::ConvolveWithAnyKernel(std::vector<double> &Waveform, std::vector<double> &Kernel, std::vector<double> &Out) {
    int KernelSize = Kernel.size();
    if(KernelSize%2==0)
    {
      //Insert a value at the mid point that averages the closest two values
      double NewVal = (Kernel[KernelSize/2] + Kernel[KernelSize/2+1])/2;
      Kernel.insert(Kernel.begin()+KernelSize/2, NewVal);
      KernelSize=KernelSize+1;
    }
    std::vector<int> X_indices(KernelSize);
    std::iota(X_indices.begin(), X_indices.end(), -KernelSize/2); // Indices
    // Now do the convolution
    for (int i = KernelSize/2; i < int(Waveform.size()) - (KernelSize)/2; i++) 
    {
        double PointSum = 0;
        for (int Index : X_indices) 
        {
            PointSum += Waveform[i - Index] * Kernel[KernelSize/2 + Index];
        }
        Out[i] = PointSum;
    }
    // Handle edges properly
    //Front edge
    for (int i = 0; i < KernelSize/2; i++) 
    {
        double PointSum = 0;
        for (int Index : std::vector<int>(X_indices.begin(), X_indices.begin()+ KernelSize/2 - i) ) 
        {
            PointSum += Waveform[i - Index] * Kernel[KernelSize/2 + Index];
        }
        Out[i] = PointSum;
    }
    //Back edge
    for (int i = int(Waveform.size()) - (KernelSize/2); i < int(Waveform.size()); i++) 
    {
        double PointSum = 0;
        for (int Index : std::vector<int>(X_indices.begin() + KernelSize/2 - (Waveform.size() - i), X_indices.end()) ) 
        {
            PointSum += Waveform[i - Index] * Kernel[KernelSize/2 + Index];
        }
        Out[i] = PointSum;
    }
    return;
}


void michelETagger::SaveVector(std::vector<double> &HistEntries, std::string Name)
{
    std::stringstream histname;
    art::ServiceHandle<art::TFileService> tfs;
    histname.str(Name); //Resets string stream to nothing
    //Create a new histogram
    //Make TH1D to hold waveform. String name is gross, right size and start and end time match timestamps
    TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), HistEntries.size(), 0, HistEntries.size()-1);
    for(unsigned int i = 0; i < HistEntries.size(); i++) {
      wvfHist->SetBinContent(i + 1, HistEntries[i]); //Loop over waveform and set bin content
    }
}


bool michelETagger::filter(art::Event& e)
{
  // Implementation of required member function here.
  //Will need CRT hit clusters
  //OpDetWaveforms
  EventNum = e.id().event();
  art::Handle< std::vector<sbnd::crt::CRTCluster>> crtClusterHandle;
  e.getByLabel(fCRTclusterLabel, crtClusterHandle);
  //Could instead use deconvolved waveforms if they are ready enough. Might give better performance
  //Plus we have to go through reco1 already so no real processing overhead associated with that choice
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fPMTLabel,waveHandle);
  //Loop over OpDetWaveforms and check for double peak signature
  int TotalFlash =  waveHandle->size()/(std::accumulate(PMTPerBoard.begin(), PMTPerBoard.end(), 0));
  bool MichelFound=false;
  for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
  {
    FlashNumForName=FlashCounter;
    int WaveIndex = FlashCounter*PMTPerBoard[0];
    //Get timestamp for comparison with CRT clusters
    double currentTimeStamp = (*waveHandle)[WaveIndex].TimeStamp(); 
    //Get size of this waveform to build a summed waveform vector
    int FlashSamples = int((*waveHandle)[WaveIndex].size()); 
    std::vector<double> SummedWaveform_TPC1(FlashSamples);
    std::vector<double> SummedWaveform_TPC2(FlashSamples);
    ConstructSummedWaveform(waveHandle, SummedWaveform_TPC1, SummedWaveform_TPC2, FlashCounter);
    SaveVector(SummedWaveform_TPC1, std::string("Event_") + std::to_string(EventNum)+std::string("_Flash_")+
    std::to_string(FlashNumForName)+std::string("TPC_1_SummedWaveform"));
    SaveVector(SummedWaveform_TPC2, std::string("Event_") + std::to_string(EventNum)+std::string("_Flash_")+
    std::to_string(FlashNumForName)+std::string("TPC_2_SummedWaveform"));
    //Check for double flash in each TPC
    TPCNumForName=1;
    bool DoubleFlash_TPC1 = DoubleFlashCheck( SummedWaveform_TPC1 );
    TPCNumForName=2;
    bool DoubleFlash_TPC2 = DoubleFlashCheck( SummedWaveform_TPC2 );
    bool CoincidentCRT = false;
    if(DoubleFlash_TPC1 || DoubleFlash_TPC2)
    {
      //Loop over crt clusters to check for good time with this opDetWaveform
      for(int ClustID; ClustID<int(crtClusterHandle->size()); ClustID++)
      {
        double CRTTimeStamp = (*crtClusterHandle)[ClustID].Ts0();
        double TDiff = (currentTimeStamp-CRTTimeStamp); //CRT should be slightly early so this is a little positive
        if(TDiff>0 && TDiff<=fCoincidentWindow) CoincidentCRT=true;
        if(TDiff > fCoincidentWindow) break; // stop looping over later clusters than this flash
      }
    }
    //end of flash processing check if we had good result 
    MichelFound = CoincidentCRT & (DoubleFlash_TPC1 || DoubleFlash_TPC2); 
    if(MichelFound) 
    {
      std::cout << "Found a Michel Candidate in Event " << EventNum << "  Flash " << FlashNumForName << std::endl;
      break;
    }
  }
  return MichelFound;
}

DEFINE_ART_MODULE(michelETagger)
