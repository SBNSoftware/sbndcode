////////////////////////////////////////////////////////////////////////
// Class:       EventIDFilterCAF
// Plugin Type: filter (Unknown Unknown)
// File:        EventIDFilterCAF_module.cc
//
// Generated at Thu Oct 21 03:39:48 2021 by Francisco Nicolas-Arnaldos using cetskelgen
// from  version .
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

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <fstream>

namespace hypana {
  class EventIDFilterCAF;
}


class hypana::EventIDFilterCAF : public art::EDFilter {
public:
  explicit EventIDFilterCAF(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EventIDFilterCAF(EventIDFilterCAF const&) = delete;
  EventIDFilterCAF(EventIDFilterCAF&&) = delete;
  EventIDFilterCAF& operator=(EventIDFilterCAF const&) = delete;
  EventIDFilterCAF& operator=(EventIDFilterCAF&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  std::string GetEventLabel(int r, int sr, int e);

  std::string fEventList;
  std::string fMCTruthLabel;

  std::ifstream fInputFile;

  int fEventID, fRunID, fSubrunID;
  std::string fEventLabel;

  std::vector<std::string> fEventLabelList;
  std::vector<int> fEventIntMode;
  std::vector<int> fEventIntEnergy;

  bool fKeepEvent;
};


hypana::EventIDFilterCAF::EventIDFilterCAF(fhicl::ParameterSet const& p)
  : EDFilter{p},
  fEventList( p.get<std::string>("EventList") ),
  fMCTruthLabel ( p.get<std::string>("MCTruthLabel",  "generator") )
{
  fEventLabelList.clear();
  fEventIntMode.clear();
  fEventIntEnergy.clear();
  fInputFile.open(fEventList);
  
  if(!fInputFile.is_open()){
    std::cout<<"Could not open input event list\n";
  }
  else{
    int e, sr, r, intmode, intenergy;
    while ( fInputFile ) {
      fInputFile>>r;
      fInputFile>>sr;
      fInputFile>>e;
      fInputFile>>intmode;
      fInputFile>>intenergy;
      fEventLabelList.push_back(GetEventLabel(r, sr, e));
      fEventIntMode.push_back(intmode);
      fEventIntEnergy.push_back(intenergy);
    }
  }

  std::cout<<"Event IDs to filter:\n";
  for(size_t i=0; i<fEventLabelList.size(); i++){
    std::cout<<fEventLabelList[i]<<std::endl;
  }
  
}

bool hypana::EventIDFilterCAF::filter(art::Event& e)
{
  //EventID
  fRunID = e.id().run();
  fSubrunID = e.id().subRun();
  fEventID = e.id().event();
  fEventLabel = GetEventLabel(fRunID, fSubrunID, fEventID);
  
  fKeepEvent=false;
  std::vector<std::string>::iterator it = std::find(fEventLabelList.begin(),fEventLabelList.end(), fEventLabel) ;

  if( it != fEventLabelList.end() ){

    art::Handle< std::vector<simb::MCTruth> > MCTruthListHandle;
    e.getByLabel(fMCTruthLabel, MCTruthListHandle);
    for (size_t n = 0; n < MCTruthListHandle->size(); n++) {
      simb::MCTruth const&  evtTruth = (*MCTruthListHandle)[n];
      int intMode = evtTruth.GetNeutrino().Mode() ;
      int intEnergy = 1000*evtTruth.GetNeutrino().Nu().E() ;

      size_t ix = it - fEventLabelList.begin();

      if(intMode == fEventIntMode[ix] && intEnergy == fEventIntEnergy[ix] ) {

        std::cout<<"Keeping Event \n";
        fKeepEvent = true;
      }
    }
  }

  return fKeepEvent;
}

void hypana::EventIDFilterCAF::beginJob()
{
  // Implementation of optional member function here.
}

void hypana::EventIDFilterCAF::endJob()
{
  // Implementation of optional member function here.
}

std::string hypana::EventIDFilterCAF::GetEventLabel(int r, int sr, int e){
  return std::to_string(r)+"_"+std::to_string(sr)+"_"+std::to_string(e);
}

DEFINE_ART_MODULE(hypana::EventIDFilterCAF)
