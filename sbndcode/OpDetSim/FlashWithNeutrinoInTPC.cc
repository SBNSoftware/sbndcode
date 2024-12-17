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
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet { //OpDet means optical detector 

  class FlashWithNeutrinoInTPC;

  class FlashWithNeutrinoInTPC : public art::EDAnalyzer {
  public:
    explicit FlashWithNeutrinoInTPC(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    //idk why things are doubled here
    FlashWithNeutrinoInTPC(FlashWithNeutrinoInTPC const &) = delete;
    FlashWithNeutrinoInTPC(FlashWithNeutrinoInTPC &&) = delete;
    FlashWithNeutrinoInTPC & operator = (FlashWithNeutrinoInTPC const &) = delete;
    FlashWithNeutrinoInTPC & operator = (FlashWithNeutrinoInTPC &&) = delete;

    //Standard ART loop fucntions
    // Required functions.
    void analyze(art::Event const & e) override;
    //Selected optional functions
    void beginJob() override;
    void endJob() override;

    void SaveMonWaveforms(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, int EventCounter, int WaveformSize);
    void ConstructMonPulse(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, std::vector<int> *MonPulse);
    
    opdet::sbndPDMapAlg pdMap; //map for photon detector types
  private:
    //Data members for analysis/output
    //Start with fcl params
    int fMonWidth; //Value*10ns MON pulse width when treshold is crossed
    //TTree *fWaveformTree;
    // Declare member data here.
    std::string fInputModuleName;
    std::vector<std::string> fOpDetsToPlot; //keep for now but probably don't need
    std::string opdetType;
    std::string opdetElectronics;
    //Helpful members to have for the analsis 
    std::vector<int> MONThresholds; //ADC counts to output a MON pulse
    std::vector<int> MTCA_Thresholds; //Number of pairs needed to trigger MTCA output
    int EventCounter;
    int fChNumber;

  };

  FlashWithNeutrinoInTPC::FlashWithNeutrinoInTPC(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p)    // More initializers here.
  {
    //Read in assorted fcl parameters
    fInputModuleName = p.get< std::string >("InputModule" );
    fOpDetsToPlot    = p.get<std::vector<std::string> >("OpDetsToPlot");
    fMonWidth        = p.get<int>("MonWidth");
    fChNumber = 0;
  }

  void FlashWithNeutrinoInTPC::beginJob()
  {
    //Initilize our set of thresholds
    int MonStep=20;
    int MonStart = 20;
    int MonEnd = 250;
    for(int i=0; (MonStart+i*MonStep)<=MonEnd; i++)
    {
        MONThresholds.push_back( MonStart+i*MonStep);
    }
    int MTCAStep = 1;
    int MTCAStart = 3;
    int MTCAEnd = 64;
    for(int i=0; MTCAStart+i*MTCAStep<=MTCAEnd; i++)
    {
        MTCA_Thresholds.push_back( MTCAStart+i*MTCAStep);
    }
    EventCounter=0;
    art::ServiceHandle<art::TFileService> tfs;
  }

  void FlashWithNeutrinoInTPC::analyze(art::Event const & e)
  {
    // Implementation of required member function here.

    art::ServiceHandle<art::TFileService> tfs; //Common art service should read about
    //int fEvNumber = e.id().event();

    art::Handle< std::vector< raw::OpDetWaveform > > waveHandle; //User handle for vector of OpDetWaveforms
    std::string Instance("PMTChannels");
    std::string Process("PMTDecoder");
    e.getByLabel(fInputModuleName, Instance, Process, waveHandle);
    auto PFPparticles = e.getValidHandle<std::vector<recob::PFParticle>>("Reco2");
    auto PFPMetaData = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>("Reco2");;
    for(auto const &MetaData : (*PFPMetaData))
    {
      auto MetaMap = MetaData.GetPropertiesMap(); //map from string to something
      for (const auto& [key, value] : MetaMap)
      {
          std::cout << '[' << key << "] = " << value << "; ";
      }
      double NuScore = MetaMap["NuScore"];
      std::cout << NuScore << std::endl;
      break;
    }
  }

  void FlashWithNeutrinoInTPC::SaveMonWaveforms(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, int EventCounter, int WaveformSize)
  {
    std::vector<int> *MonPulse = new std::vector<int>(WaveformSize); //Maybe works? Should really add a waveform size getter 
    ConstructMonPulse(waveHandle, MonThreshold, MonPulse);
    std::stringstream histname;
    histname << "event_" << EventCounter <<"_Mon";
    art::ServiceHandle<art::TFileService> tfs;
    TH1D *MonHist = tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(), 
                                    MonPulse->size(), 0.0, MonPulse->size()-1); //so this just breaks
    for(unsigned int i = 0; i < MonPulse->size(); i++) {
        MonHist->SetBinContent(i + 1, (double)(*MonPulse)[i]); //Loop over waveform and set bin content
      }
    delete MonPulse;
  }

  void FlashWithNeutrinoInTPC::ConstructMonPulse(art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, int MonThreshold, std::vector<int> *MonPulse)
  {
    std::vector<int> Crossings;
    std::fill(MonPulse->begin(), MonPulse->end(), 0);
    //Loop over the entries in our waveform vector
    //We care about getting the pairing correct
    for(auto const& wvf : (*waveHandle))
    {
        bool SecondInPair=false;
        if(fChNumber==int(wvf.ChannelNumber()+1) && fChNumber%2==0) //Channels seem to be in reverse order
        {
            SecondInPair=true;
        }
        if(!SecondInPair)
        {
            Crossings.clear(); //Reset list of crossings
        }
        fChNumber = wvf.ChannelNumber(); //Get the channel number for this waveform
        if(fChNumber==15) continue; //copy of trigger in signal doesn't need to be analyzed
        opdetType = pdMap.pdType(fChNumber); //opdet::sbndPDMapAlg pdMap; //map for photon detector types. Lists arapuca vs PMT by channel num I guess
        opdetElectronics = pdMap.electronicsType(fChNumber); //? Maybe what kind of CAEN it is for digitization rate
        if (std::find(fOpDetsToPlot.begin(), fOpDetsToPlot.end(), opdetType) == fOpDetsToPlot.end()) {continue;}
        int k=1;
        int Baseline=0;
        int NumPoints = 10;
        for(int Index=0; Index<NumPoints; Index++)
        {
          Baseline = Baseline + wvf[Index];
        }
        Baseline=Baseline/NumPoints;
        while(k < int(MonPulse->size())-1 )
        {
            //Check for crossings
            bool CrossedThreshold = (((wvf[k-1]-Baseline)>-MonThreshold) && ((wvf[k]-Baseline)<-MonThreshold));
            if(CrossedThreshold)
            {
              bool Overlap=false;
              bool CrossLeft = std::any_of(Crossings.begin(), Crossings.end(), [k, this](int x){return (k>(x-4*fMonWidth))&&(k<x) ; } );
              bool CrossRight = std::any_of(Crossings.begin(), Crossings.end(), [k, this](int x){return (k<(x+4*fMonWidth)&&(k>x)) ;} );
              Overlap = (CrossLeft || CrossRight);
              Crossings.push_back(k);
              //Need some logic for adjacent channels and for non-overlapping pulses
              //Increment the MON Pulse
              if(Overlap)
              {
                  int StartIndex=0; 
                  int EndIndex=0;
                  int LowBound = *std::lower_bound(Crossings.begin(), Crossings.end(), k); //First value that is not less than k
                  if(CrossLeft)
                  {
                    StartIndex = k;
                    EndIndex = LowBound;
                  }
                  else if(CrossRight)
                  {
                    EndIndex = k+4*fMonWidth;
                    StartIndex = *std::lower_bound(Crossings.begin(), Crossings.end(), k-4*fMonWidth); //Maybe has bug when two crossings overlap in one channel
                  }
                  if(StartIndex>=int(MonPulse->size())) StartIndex=MonPulse->size()-1; //should never happen
                  if(EndIndex>=int(MonPulse->size())) EndIndex=MonPulse->size()-1; //May happen
                  if(StartIndex<EndIndex) std::for_each(MonPulse->begin()+StartIndex, MonPulse->begin()+EndIndex, [](int &X){return X=X+1;}); //extra check for safety
                  k=EndIndex+1;
              }
              else
              {
                int StartIndex = k;
                int EndIndex = k+4*fMonWidth;
                if(StartIndex>=int(MonPulse->size())) StartIndex=MonPulse->size()-1; //should never happen
                if(EndIndex>=int(MonPulse->size())) EndIndex=MonPulse->size()-1; //May happen
                if(StartIndex<EndIndex) std::for_each(MonPulse->begin()+StartIndex, MonPulse->begin()+EndIndex, [](int &X){return X=X+1;});
                k=k+4*fMonWidth+1;
              }
            }
            else k=k+1;
        }
    }
  }

  void FlashWithNeutrinoInTPC::endJob()
  {
  }
    DEFINE_ART_MODULE(opdet::FlashWithNeutrinoInTPC) // Magic line that has to be here

}