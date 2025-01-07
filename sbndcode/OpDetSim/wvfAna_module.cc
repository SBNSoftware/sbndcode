////////////////////////////////////////////////////////////////////////
// Class:       wvfAna
// Module Type: analyzer
// File:        wvfAna_module.cc
//
// Analyzer to read optical waveforms
//
// Authors: L. Paulucci and F. Marinho
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
#include "sbndcode/Decoders/PTB/sbndptb.h"


#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet {

  class wvfAna;

  class wvfAna : public art::EDAnalyzer {
  public:
    explicit wvfAna(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    wvfAna(wvfAna const &) = delete;
    wvfAna(wvfAna &&) = delete;
    wvfAna & operator = (wvfAna const &) = delete;
    wvfAna & operator = (wvfAna &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    //Selected optional functions
    void beginJob() override;
    void endJob() override;

    opdet::sbndPDMapAlg pdMap; //map for photon detector types

  private:

    size_t fEvNumber;
    size_t fChNumber;
    double fSampling;
    double fSampling_Daphne;
    double fStartTime;
    double fEndTime;
    bool fSaveSummedWaveform;
    bool fJustOne;
    int fPMTperCAEN;
    int fChToPlot;
    bool fCheckTiming;
    std::string fSummedInputModuleName;
    std::string fPTBLabel;
    int fTotalCAENBoards;
    bool fTriggerOnly;
    //TTree *fWaveformTree;

    // Declare member data here.
    std::string fInputModuleName;
    std::vector<std::string> fOpDetsToPlot;
    std::stringstream histname;
    std::string opdetType;
    std::string opdetElectronics;
  };


  wvfAna::wvfAna(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p)  // ,
    // More initializers here.
  {
    fInputModuleName = p.get< std::string >("InputModule" );
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
    fSaveSummedWaveform = p.get<bool>("SaveSummedWaveform", false);
    fSummedInputModuleName = p.get< std::string >("SummedInputModule", ":::");

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampling = clockData.OpticalClock().Frequency(); // MHz
    fSampling_Daphne = p.get<double>("DaphneFrequency" );
    fJustOne = p.get<bool>("JustOne", false );
    fChToPlot = p.get<int>("ChToPlot", -1 );
    fCheckTiming = p.get<bool>("CheckTiming", false);
    fPTBLabel = p.get< std::string >("PTBLabel",  "ptbdecoder::DECODE");
    fTotalCAENBoards = p.get<int>("TotalCAENBoards", 8);
    fTriggerOnly = p.get<bool>("TriggerOnly", false);
    fPMTperCAEN = p.get<int>("PMTPerCAEN", 15);

  }

  void wvfAna::beginJob()
  {

  }

  void wvfAna::analyze(art::Event const & e)
  {
    // Implementation of required member function here.
    std::cout << "My module on run # " << e.id().run() << "event #" << e.id().event() << std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEvNumber = e.id().event();

    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
    e.getByLabel(art::InputTag(fInputModuleName), waveHandle);
    //delete from here Nov 27
    std::cout << "Printing timing information for event + HLT + LLT + Waveforms" << std::endl;
    std::cout << "Run #" <<  e.id().run() << " event #" << e.id().event() << " has timestamp " << e.time().timeHigh() << "."<<e.time().timeLow() <<" s" << std::endl;
    art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle_2;
    e.getByLabel(fPTBLabel,ptbHandle_2);
    std::cout << "\t It has " << ptbHandle_2->size() << " PTB fragments" << std::endl;
    for(int index=0; index<int(ptbHandle_2->size()); index++)
    {
      auto ptb = (*ptbHandle_2)[index];
      auto hltrigs = ptb.GetHLTriggers();
      auto LLtrigs = ptb.GetLLTriggers();
      std::cout << "\t \t ptb fragment " << index << " with " << hltrigs.size() << " HLTs" <<std::endl;
      for(int HLT=0; HLT<int(hltrigs.size()); HLT++)
      {
        int Power=0;
        while(Power<64)
        {
          if(hltrigs[HLT].trigger_word & (0x1 << Power)) break;
          else Power=Power+1;
        }
        std::cout << "\t \t \t HLT " << HLT << " has trigger word " <<  hltrigs[HLT].trigger_word  << " HLT word of type " << Power << " at time " << hltrigs[HLT].timestamp*20 << " ns UNIX" << std::endl;
      }
      std::cout << "\t \t ptb fragment " << index << " with " << LLtrigs.size() << " LLTs" <<std::endl;
      for(int LLT=0; LLT<int(LLtrigs.size()); LLT++)
      {
        int Power=0;
        while(Power<64)
        {
          if(LLtrigs[LLT].trigger_word & (0x1 << Power)) break;
          else Power=Power+1;
        }
        std::cout << "\t \t \t LLT " << LLT << " has trigger word " <<  LLtrigs[LLT].trigger_word  << " LLT word of type " << Power << " at time " << LLtrigs[LLT].timestamp*20 << " ns UNIX" << std::endl;
      }
    }
    int PMTPerCAEN=fPMTperCAEN;
    int TotalFlash = waveHandle->size()/(fTotalCAENBoards*PMTPerCAEN);
    std::cout << "\t It also has " << TotalFlash << " opDetWaveforms " << std::endl;
    for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
    {
      int WaveIndex = FlashCounter*PMTPerCAEN;
      double currentTimeStamp = (*waveHandle)[WaveIndex].TimeStamp(); //Getting channel 17 every time!
      std::cout << "\t\t\t Flash " << FlashCounter << " is at " << currentTimeStamp << " us relative to the event?" << std::endl;
      
    }


    int GoodFlashIndex=0;
    double SmallestTimestamp=0;
    for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
    {
      int WaveIndex = FlashCounter*PMTPerCAEN;
      double currentTimeStamp = (*waveHandle)[WaveIndex].TimeStamp(); //Getting channel 17 every time!
      if(FlashCounter==0) SmallestTimestamp=TMath::Abs(currentTimeStamp);
      if(TMath::Abs(currentTimeStamp)<SmallestTimestamp) 
      {
        GoodFlashIndex=FlashCounter;
        SmallestTimestamp=TMath::Abs(currentTimeStamp);
      }
    }
    //delete to here Nov 27



    if(!waveHandle.isValid()) {
      std::cout << Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
    }

    // // example of usage for pdMap.getCollectionWithProperty()
    // //
    // // define a container
    // auto inBoxTwo = pdMap.getCollectionWithProperty("pds_box", 2);
    // // you can cout the whole json object
    // std::cout << "inBoxTwo:\t" << inBoxTwo << "\n";
    // // traverse its components in a loop
    // for (auto const &e: inBoxTwo) {
    //   std::cout << e["pd_type"] << " " << e["channel"] << ' ' << "\n";
    // }

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

    int hist_id = 0;

    for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
    {
    for(int CurrentBoard=0; CurrentBoard<fTotalCAENBoards; CurrentBoard++)
        {
            //Loop over each PMT in a board
            for(int CAENChannel=0; CAENChannel<PMTPerCAEN; CAENChannel++)
            {
              int WaveIndex = CAENChannel + FlashCounter*PMTPerCAEN + CurrentBoard*PMTPerCAEN*TotalFlash;
              auto const& wvf = (*waveHandle)[WaveIndex];
              fChNumber = wvf.ChannelNumber();
              if(fJustOne && fChToPlot!=int(fChNumber)) continue;
              if(fCheckTiming)
              {
                art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
                e.getByLabel(fPTBLabel,ptbHandle);
                for(int index=0; index<int(ptbHandle->size()); index++)
                {
                  auto ptb = (*ptbHandle)[index];
                  auto hltrigs = ptb.GetHLTriggers();
                  std::cout << "ptb fragment " << index << " with " << hltrigs.size() << " HLTs" <<std::endl;
                  for(int HLT=0; HLT<int(hltrigs.size()); HLT++)
                  {
                    std::cout << "This guy has " <<  hltrigs[HLT].trigger_word  << " HLT word of type " << hltrigs[HLT].word_type << " at time " << hltrigs[HLT].timestamp*20 << " ns UNIX" << std::endl;
                    std::cout << "The event has timestamp " << e.time().timeHigh() << " and my wvfm has " << wvf.TimeStamp() << std::endl;
                    std::cout << " event time low " <<  e.time().timeLow() << " event time high " << e.time().timeHigh() <<std::endl;
                  }
                }
              }
              if(fChNumber>=900)
              {
                opdetType = "pmt_coated";
                opdetElectronics = "";
              }
              else
              {
                opdetType = pdMap.pdType(fChNumber);
                opdetElectronics = pdMap.electronicsType(fChNumber);
              }
              if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
              histname.str(std::string());
              if(GoodFlashIndex!=FlashCounter)
              {
              histname <<"run_"<< e.id().run()  <<"_event_" << fEvNumber
                      << "_opchannel_" << fChNumber
                      << "_" << opdetType
                      << "_" << FlashCounter;
              if(fTriggerOnly) continue;
              }
              else
              {
              histname <<"run_"<< e.id().run()  <<"_event_" << fEvNumber
                      << "_opchannel_" << fChNumber
                      << "_" << opdetType
                      << "_" << FlashCounter << "_TriggerPulse";
              }
              fStartTime = wvf.TimeStamp(); //in us
              if (opdetElectronics == "daphne"){
                fEndTime = double(wvf.size()) / fSampling_Daphne + fStartTime;
              } //in us
              else{
                fEndTime = double(wvf.size()) / fSampling + fStartTime;
              } //in us

              //Create a new histogram
              TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", fStartTime), wvf.size(), fStartTime, fEndTime);
              for(unsigned int i = 0; i < wvf.size(); i++) {
                wvfHist->SetBinContent(i + 1, (double)wvf[i]);
              }
              hist_id++;
            }
        }
    }
    if(fSaveSummedWaveform)
    {
      //Load approrpirate 
      art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
      e.getByLabel(art::InputTag(fSummedInputModuleName), waveHandle);
      if(!waveHandle.isValid()) return;
      for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
      {
        auto const& wvf = (*waveHandle)[FlashCounter];
        histname.str(std::string());
        if(GoodFlashIndex!=FlashCounter){ 
          histname <<"run_"<< e.id().run()  <<"_event_" << fEvNumber << "_SummedWaveform_" << FlashCounter;
          if(fTriggerOnly) continue;
        }
        else histname <<"run_"<< e.id().run()  <<"_event_" << fEvNumber << "_SummedWaveform_" << FlashCounter<< "_TriggerPulse";
        TH1D *wvfHist = tfs->make< TH1D >(histname.str().c_str(), TString::Format(";t - %f (#mus);", fStartTime), wvf.size(), 0, wvf.size());
        for(unsigned int i = 0; i < wvf.size(); i++) 
        {
          wvfHist->SetBinContent(i + 1, (double)wvf[i]);
        }
      }
    }
  }

  void wvfAna::endJob()
  {
  }

  DEFINE_ART_MODULE(opdet::wvfAna)

}