////////////////////////////////////////////////////////////////////////
// Class:       RawHitFinder
// Plugin Type: producer (Unknown Unknown)
// File:        RawHitFinder_module.cc
//
// Generated at Tue Oct 22 08:58:30 2024 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
//Takes in a baseline-subtracted pmt OpDetWaveform then performs an integrator
//based hit finding algorithm
//Has 3 options: 
// 1. perform hit finding on all channels, 
// 2. perform hit finding on summed channel then copy hit bounds to each pmt
// 3. perform hit finding on summed channel then look for hits in those bounds on each pmt

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"



#include <memory>
#include <numeric>

namespace opdet {
  class RawHitFinder;
}


class opdet::RawHitFinder : public art::EDProducer {
public:
  explicit RawHitFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RawHitFinder(RawHitFinder const&) = delete;
  RawHitFinder(RawHitFinder&&) = delete;
  RawHitFinder& operator=(RawHitFinder const&) = delete;
  RawHitFinder& operator=(RawHitFinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  void IntegratorHitFinder(std::vector<int> &HitStarts, std::vector<int> &HitEnds, 
                                    double StartThreshold, double EndThreshold, int MinWidth, int IntegrationWidth, int PreSample, int PostSample,
                                    int ScanStartSample, int ScanEndSample, raw::OpDetWaveform wvf);
  void MakeOpHits(std::vector<recob::OpHit> &Hit_Vec, std::vector<int> HitStarts, std::vector<int> HitEnds, raw::OpDetWaveform wvf);
  // Declare member data here.
  std::string fBaselineSubtractedChName;
  std::string fSummedSubtractedName;
  bool fAllAlgs;
  bool fRunCopy;
  bool fRunRestrictedSearch;
  bool fRunPerChannelSearch;
  double fSummedAreaThreshold;
  double fSummedAreaEndThreshold;
  int fSummedMinWidth;
  int fSummedIntegrationWidth;
  int fSummedPreSample;
  int fSummedPostSample;
  double fPerChannelAreaThreshold;
  double fPerChannelAreaEndThreshold;
  int fPerChannelMinWidth;
  int fPerChannelIntegrationWidth;
  int fPerChannelPreSample;
  int fPerChannelPostSample;
  int fTotalCAENBoards;
  std::string CopyInstanceName;
  std::string RestrictInstanceName;
  std::string FreeInstanceName;
  std::vector<std::string> fOpDetsToPlot;
  opdet::sbndPDMapAlg pdMap;
  
};


opdet::RawHitFinder::RawHitFinder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  
  fBaselineSubtractedChName =  p.get< std::string >("BaselineSubtractedChannelNames" );
  fSummedSubtractedName = p.get< std::string >("BaselineSubtractedSummedNames" );
  fAllAlgs = p.get<bool>("RunAllAlg", true);
  fRunCopy = p.get<bool>("RunCopyAlg", false);
  fRunRestrictedSearch = p.get<bool>("RunRestrictedSearch", false);
  fRunPerChannelSearch = p.get<bool>("RunPerChannleSearch", false);
  fTotalCAENBoards = p.get<int>("totalCAENBaord", 9);
  //Summed search fcl params
  fSummedAreaThreshold = p.get<double>("SummedAreaThreshold", 400);
  fSummedAreaEndThreshold = p.get<double>("SummedAreaEndThreshold", 80);
  fSummedMinWidth = p.get<int>("SummedHitMinWidth", 10);
  fSummedIntegrationWidth = p.get<int>("SummedIntegrationWidth", 5);
  fSummedPreSample = p.get<int>("SummedHitPreSamples", 10);
  fSummedPostSample = p.get<int>("SummedHitPostSamples", 10);
  //Perchannel fcl params
  fPerChannelAreaThreshold = p.get<double>("ChannelAreaThreshold", 120);
  fPerChannelAreaEndThreshold = p.get<double>("ChannelAreaEndThreshold", 40);
  fPerChannelMinWidth = p.get<int>("ChannelHitMinWidth", 10);
  fPerChannelIntegrationWidth = p.get<int>("ChannelIntegrationWidth", 5);
  fPerChannelPreSample = p.get<int>("ChannelHitPreSamples", 2);
  fPerChannelPostSample = p.get<int>("ChannelHitPostSamples", 2);
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  CopyInstanceName = std::string("RawHitsSummedCopy");
  RestrictInstanceName = std::string("RawHitsSummedRestricted") ;
  FreeInstanceName = std::string("RawHitsFreeBounds");
  fOpDetsToPlot = p.get<std::vector<std::string>>("OpDetsToPlot", {"pmt_coated", "pmt_uncoated"});
  if(fAllAlgs || fRunCopy) produces< std::vector< recob::OpHit > >(CopyInstanceName);
  if(fAllAlgs || fRunRestrictedSearch) produces< std::vector< recob::OpHit > >(RestrictInstanceName);
  if(fAllAlgs || fRunPerChannelSearch) produces< std::vector< recob::OpHit > >(FreeInstanceName);
}

void opdet::RawHitFinder::produce(art::Event& e)
{
  // Implementation of required member function here.
  //Get channel and summed waveform handles
  std::cout << "My Hit Finding module on event #" << e.id().event() << std::endl;
  art::Handle< std::vector< raw::OpDetWaveform > > channel_waveHandle;
  e.getByLabel(art::InputTag(fBaselineSubtractedChName), channel_waveHandle);
  art::Handle< std::vector< raw::OpDetWaveform > > summed_waveHandle;
  e.getByLabel(art::InputTag(fSummedSubtractedName), summed_waveHandle);

  std::unique_ptr< std::vector< recob::OpHit > > PerChannelHitPtr_Copy(new std::vector< recob::OpHit >);
  std::unique_ptr< std::vector< recob::OpHit > > PerChannelHitPtr_Restricted(new std::vector< recob::OpHit >);
  std::unique_ptr< std::vector< recob::OpHit > > PerChannelHitPtr_Free(new std::vector< recob::OpHit >);

  int NFlashes = (*summed_waveHandle).size();
  int PMTPerBoard = 15;
  if(fAllAlgs || fRunCopy || fRunRestrictedSearch)
  {
    for(int FlashCounter=0; FlashCounter<NFlashes; FlashCounter++)
    {
      std::vector<int> SummedHitStarts; 
      std::vector<int> SummedHitEnds;
      auto const& wvf = (*summed_waveHandle)[FlashCounter]; // raw::OpDetWaveform 
      IntegratorHitFinder(SummedHitStarts, SummedHitEnds, fSummedAreaThreshold, fSummedAreaEndThreshold, fSummedMinWidth, fSummedIntegrationWidth, 
      fSummedPreSample, fSummedPostSample, 0, int(wvf.size())-1, wvf); //Get all summed hit starts
      for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
        {
            //Loop over each PMT in a board
            for(int CAENChannel=0; CAENChannel<PMTPerBoard; CAENChannel++)
            {
                //To move to the next flash in a given CAEN we need to shift forward by 15 PMT
                //To move to the next board we have to shift forward by 15 PMT*Number flash
                //To move to next PMT we have to shift forward by 1
                //I think ARAPUCA live in their own object? I hope so
                int Index = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NFlashes;
                auto const& ch_wvfm = (*channel_waveHandle)[Index];
                int fChNumber = ch_wvfm.ChannelNumber();
                if(fChNumber>899) continue;
                std::string opdetType = pdMap.pdType(fChNumber);
                if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
                if(fRunCopy || fAllAlgs)
                {
                  MakeOpHits(*PerChannelHitPtr_Copy, SummedHitStarts, SummedHitEnds, ch_wvfm);
                }
                if(fRunRestrictedSearch || fAllAlgs)
                {
                  std::vector<int> ChannelHitStarts; 
                  std::vector<int> ChannelHitEnds;
                  for(int k=0; k<int(SummedHitStarts.size()); k++)
                  { //SummedHitStarts is refreshed for every flash so this should work fine
                    if(SummedHitEnds[k]>=int(ch_wvfm.size())) continue; // sometimes fragements are missing
                    IntegratorHitFinder(ChannelHitStarts, ChannelHitEnds, fPerChannelAreaThreshold, fPerChannelAreaEndThreshold, fPerChannelMinWidth, 
                    fPerChannelIntegrationWidth, fPerChannelPreSample, fPerChannelPostSample, SummedHitStarts[k], SummedHitEnds[k], ch_wvfm); //Look in tiny window for this channel
                  }
                  MakeOpHits(*PerChannelHitPtr_Restricted, ChannelHitStarts, ChannelHitEnds, ch_wvfm);
                  //Think we just want to put object into event once, so have to finish loop over flashes
                }
            }
        }
    }
  }
  //I'll make all sets of hits nicely time-ordered
  if(fAllAlgs || fRunPerChannelSearch)
  {
    for(int FlashCounter=0; FlashCounter<NFlashes; FlashCounter++)
    {
      for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
        {
          for(int CAENChannel=0; CAENChannel<PMTPerBoard; CAENChannel++)
            {
              int Index = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NFlashes;
              auto const& ch_wvfm = (*channel_waveHandle)[Index];
              int fChNumber = ch_wvfm.ChannelNumber();
              if(fChNumber>899) continue;
              std::string opdetType = pdMap.pdType(fChNumber);
              if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
              std::vector<int> ChannelHitStarts; 
              std::vector<int> ChannelHitEnds;
              IntegratorHitFinder(ChannelHitStarts, ChannelHitEnds, fPerChannelAreaThreshold, fPerChannelAreaEndThreshold, fPerChannelMinWidth, 
                        fPerChannelIntegrationWidth, fPerChannelPreSample, fPerChannelPostSample, 0, int(ch_wvfm.size())-1, ch_wvfm); //Get hits for this channel
              MakeOpHits(*PerChannelHitPtr_Free, ChannelHitStarts, ChannelHitEnds, ch_wvfm);
            }
        }
    }
  }
  if(fAllAlgs || fRunCopy) e.put(std::move(PerChannelHitPtr_Copy), CopyInstanceName);
  if(fAllAlgs || fRunRestrictedSearch) e.put(std::move(PerChannelHitPtr_Restricted), RestrictInstanceName);
  if(fAllAlgs || fRunPerChannelSearch) e.put(std::move(PerChannelHitPtr_Free), FreeInstanceName);
}


void opdet::RawHitFinder::MakeOpHits(std::vector<recob::OpHit> &Hit_Vec, std::vector<int> HitStarts, std::vector<int> HitEnds, raw::OpDetWaveform wvf)
{
  //Takes in list of hit starts and hit ends to make a vector of recob ophits
  //recob:Ophit construction needs 
 	//OpHit (int opchannel, double peaktime, double peaktimeabs, unsigned short frame, double width, double area, 
  //        double peakheight, double pe, double fasttototal)
  //with these inputs we can get timing information from index + wvf time
  //Amplitude and area come from vf
  //width from index
  //whats fasttotal and frame?
  //SBNDOpHitFinder_module has fasttotal=0;
  //Im going to just make frame the start index....
  //I guess peak time is relative to trigger? While peaktimeabs is unix? I'll make them equal
  int opChannel = wvf.ChannelNumber();
  double PMTSampleSize=2; //ns
  double SinglePEArea = 150; //adc ns
  for(int i=0; i<int(HitStarts.size()); i++)
  {
    unsigned short frame = HitStarts[i];
    double width = HitEnds[i] - HitStarts[i];
    double Area = std::accumulate(wvf.begin()+HitStarts[i], wvf.begin()+HitEnds[i], 0)*PMTSampleSize;
    int PeakIndex = HitStarts[i];
    for(int j=HitStarts[i]; j<=HitEnds[i]; j++)
    {
      if(wvf[PeakIndex]>wvf[j]) PeakIndex=j;
    }
    double peaktime = wvf.TimeStamp() + PMTSampleSize*(PeakIndex);
    double peaktimeabs = peaktime; //Not sure what this is supposed to be
    double fasttototal = 0;
    double pe = Area/SinglePEArea;
    double PeakHeight = wvf[PeakIndex];
    recob::OpHit TempHit(opChannel, peaktime, peaktimeabs, frame, width, Area, PeakHeight, pe, fasttototal);
    Hit_Vec.push_back(TempHit);
  }
}

void opdet::RawHitFinder::IntegratorHitFinder(std::vector<int> &HitStarts, std::vector<int> &HitEnds, 
                                    double StartThreshold, double EndThreshold, int MinWidth, int IntegrationWidth, int PreSample, int PostSample,
                                    int ScanStartSample, int ScanEndSample, raw::OpDetWaveform wvf)
{
  //Generate area waveform and look for threshold crossings appropriately
  //Fills in hit start index and hit end index
  //Area waveform will only get filled in between Min sample and max sample

  //Start with making the area vector
  std::vector<int> Area(wvf.size());
  double PMTSampleSize = 2; //ns
  bool Fire=false; //bool for tracking if we are in a pulse
  int HitStart=0; //Holder for current hit start
  int HitEnd=0;  //Holder for current hit end
  int i =ScanStartSample;
  while(i<=ScanEndSample)
  {
    bool PulseAdded=false;
    int End = i+IntegrationWidth;
    if(End>int(wvf.size())-1 ) End = wvf.size()-1;
    //Could add a postive requirement
    Area[i] = -1*std::accumulate(wvf.begin()+i, wvf.begin()+End, 0)*PMTSampleSize;
    if(!Fire && Area[i] >= StartThreshold)
    {
      Fire=true;
      HitStart = i - PreSample;
      if(HitStart<0) HitStart=0;
    }
    if(Fire && ((Area[i] <= EndThreshold) || i==(ScanEndSample-1)))
    {
      Fire=false;
      HitEnd = i+PostSample;
      if( (HitEnd - HitStart) >= MinWidth)
      {
        if(HitEnd>=int(wvf.size())) HitEnd=int(wvf.size())-1;
        HitStarts.push_back(HitStart);
        HitEnds.push_back(HitEnd);
        PulseAdded=true; 
      }
    }
    if(PulseAdded) i = int(HitEnd+1);
    else i=i+1;
  }
}

DEFINE_ART_MODULE(opdet::RawHitFinder)
