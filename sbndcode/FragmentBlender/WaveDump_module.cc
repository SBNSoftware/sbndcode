////////////////////////////////////////////////////////////////////////
// Class:       WaveDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        WaveDump_module.cc
//
// Generated at Thu Nov  5 13:45:51 2026 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "lardataobj/RawData/RawDigit.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

namespace sbnd {
  class WaveDump;
}


class sbnd::WaveDump : public art::EDAnalyzer {
public:
  explicit WaveDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveDump(WaveDump const&) = delete;
  WaveDump(WaveDump&&) = delete;
  WaveDump& operator=(WaveDump const&) = delete;
  WaveDump& operator=(WaveDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  TTree* evtTree;
  std::vector<int> WirePlane;
  std::vector<int> WireChannel;
  std::vector<int> WireWireNum;
  std::vector<int> WireTPC;
  std::vector<std::vector<float>> WireADC;
  std::string treeName      = "waveDump";

  void Clear(){
    WireADC.clear();
    WireTPC.clear();
    WireChannel.clear();
    WireWireNum.clear();
    WirePlane.clear();
    }

  void MakeTree(){
    art::ServiceHandle<art::TFileService> tfs;
    evtTree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    evtTree->Branch("Wire_Plane","std::vector<int>",&WirePlane); 
    evtTree->Branch("Wire_Channel","std::vector<int>",&WireChannel); 
    evtTree->Branch("Wire_Channel","std::vector<int>",&WireWireNum); 
    evtTree->Branch("Wire_TPC","std::vector<int>",&WireTPC); 
    evtTree->Branch("Wire_ADC","std::vector<std::vector<float>>",&WireADC);
    }
private:
    std::string fDataLabel;

};


sbnd::WaveDump::WaveDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
    fDataLabel= p.get<std::string>   ("DataLabel");
    MakeTree();
}

void sbnd::WaveDump::analyze(art::Event const& e)
{
    Clear();
    //Read in the waveforms (assuming raw::Digits)
    art::Handle< std::vector<raw::RawDigit> > digitHandle;
    //Loop over every raw digit and record appropriate info
    for(int i=0; i<int((*digitHandle)->size()); i++)
    {
        auto ThisDigit = (*digitHandle)[i];
        int Channel = ThisDigit.Channel();
        WireChannel.push_back(Channel);
        int InductionWires=1984;
        int CollectionWires=1664;
        //Construct map to plane/tpc 
        int TPC=0;
        if(Channel>= 2*InductionWires + CollectionWires) TPC=1;
        int TempChannel=Channel;
        if(TPC==1) TempChannel = TempChannel - (2*InductionWires + CollectionWires);
        int Plane=0;
        if(TempChannel<InductionWires) Plane=0;
        else if(TempChannel<2*InductionWires) Plane=1;
        else Plane=2;
        if(Plane==1) TempChannel=TempChannel-InductionWires;
        if(Plane==2) TempChannel=TempChannel-2*InductionWires;
        WirePlane.push_back(Plane);
        WireTPC.push_back(TPC);
        WireWireNum.push_back(TempChannel);
        //loop over ADCs
        auto TempADCs = ThisDigit.ADCs();
        WireADC.push_back({});
        for(int j=0; j<int(TempADCs.size()); j++)
        {
            WireADC[WireADC.size()-1].push_back(TempADCs[j]);
        }
    }
    evtTree->Fill();
}

DEFINE_ART_MODULE(sbnd::WaveDump)
