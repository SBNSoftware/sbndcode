////////////////////////////////////////////////////////////////////////
// Class:       RecoBWireDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoBWireDump_module.cc
//
// Generated at Tue Mar 18 10:43:26 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Wire.h"
#include "art_root_io/TFileService.h"


#include "TTree.h"

class RecoBWireDump;


class RecoBWireDump : public art::EDAnalyzer {
public:
  explicit RecoBWireDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoBWireDump(RecoBWireDump const&) = delete;
  RecoBWireDump(RecoBWireDump&&) = delete;
  RecoBWireDump& operator=(RecoBWireDump const&) = delete;
  RecoBWireDump& operator=(RecoBWireDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  void ResetVars();
  std::string fWireLabel;
  // Declare member data here.
  std::vector<int> ChID;
  std::vector<int> view;
  std::vector<std::vector<float>> ADC;
  int _run;        ///< The run number
  int _event;      ///< The event number
  TTree* evtTree;
};


void RecoBWireDump::ResetVars()
{
  ChID.clear();
  view.clear();
  ADC.clear();
  _run=-1;
  _event=-1;
}

RecoBWireDump::RecoBWireDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{
  fWireLabel = p.get< std::string >("InputInstance", "sptpc2d:gauss:Reco1" );
  art::ServiceHandle<art::TFileService> tfs;
  evtTree = tfs->make<TTree>("SomeWires","Wire Dump Tree");
  evtTree->Branch("event",&_event);
  evtTree->Branch("run",&_run);
  evtTree->Branch("ChID",&ChID);
  evtTree->Branch("View",&view);
  evtTree->Branch("Signal",&ADC);
}

void RecoBWireDump::analyze(art::Event const& e)
{
  //Load recob wires
  ResetVars();
  _event = e.id().event();
  _run = e.run();
  art::Handle<std::vector<recob::Wire>> wireHandle;
  e.getByLabel(fWireLabel,wireHandle);
  for(int i_wire=0; i_wire<int(wireHandle->size()); i_wire++)
  {
    recob::Wire ThisWire = (*wireHandle)[i_wire];
    ChID.push_back(ThisWire.Channel());
    view.push_back(ThisWire.View());
    ADC.push_back(ThisWire.Signal());
  }
  evtTree->Fill();
}

DEFINE_ART_MODULE(RecoBWireDump)
