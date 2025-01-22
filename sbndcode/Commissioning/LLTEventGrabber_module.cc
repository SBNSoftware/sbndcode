////////////////////////////////////////////////////////////////////////
// Class:       LLTEventGrabber
// Plugin Type: filter (Unknown Unknown)
// File:        LLTEventGrabber_module.cc
//
// Generated at Wed Jan 22 15:33:13 2024 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
#ifndef LLTEventGrabber_H
#define LLTEventGrabber_H
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

#include <memory>

namespace TriggerWork {
  class LLTEventGrabber;
}


class TriggerWork::LLTEventGrabber : public art::EDFilter {
public:
  explicit LLTEventGrabber(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LLTEventGrabber(LLTEventGrabber const&) = delete;
  LLTEventGrabber(LLTEventGrabber&&) = delete;
  LLTEventGrabber& operator=(LLTEventGrabber const&) = delete;
  LLTEventGrabber& operator=(LLTEventGrabber&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  std::vector<int> fLLT;
  int fEventHLT;
  int fTrigWindow_ns;
  std::string fPTBLabel;

private:

  // Declare member data here.

};


TriggerWork::LLTEventGrabber::LLTEventGrabber(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fLLT = p.get<std::vector<int>>("LLT");
  fEventHLT = p.get<int>("eventHLT");
  fPTBLabel = p.get<std::string>("PTBLabel");
  fTrigWindow_ns = p.get<int>("TrigWindowNs");
}

bool TriggerWork::LLTEventGrabber::filter(art::Event& e)
{
  bool GrabThisEvent=false;
  //Get timestamp of HLT that triggered event
  long eventTimestamp=0;
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  e.getByLabel(fPTBLabel,ptbHandle);
  for(int index=0; index<int(ptbHandle->size()); index++)
  {
    auto ptb = (*ptbHandle)[index];
    auto hltrigs = ptb.GetHLTriggers();
    for(int HLT=0; HLT<int(hltrigs.size()); HLT++)
    {
      int Power=0;
      std::vector<int> GrabbedHLTs;
      while(Power<32)
      {
        if(hltrigs[HLT].trigger_word & (0x1 << Power))
        {
          GrabbedHLTs.push_back(Power);
        }
        Power=Power+1;
      }
      //Do a loop over vector elements of HLT types
      for(int i=0; i<int(GrabbedHLTs.size()); i++ )
      {
        if(GrabbedHLTs[i] == fEventHLT) //Takes first instance of a given HLT
        {
            eventTimestamp = hltrigs[HLT].timestamp*20;
            break; 
        }
      }
      if(eventTimestamp>0) break;
    }
    if(eventTimestamp>0) break;
  }
   //Loop over LLT
  for(int index=0; index<int(ptbHandle->size()); index++)
  {
    auto ptb = (*ptbHandle)[index];
    auto lltrigs = ptb.GetLLTriggers();
    for(int LLT=0; LLT<int(lltrigs.size()); LLT++)
    {
      int Power=0;
      std::vector<int> GrabbedLLTs;
      while(Power<32)
      {
        if(lltrigs[LLT].trigger_word & (0x1 << Power))
        {
          GrabbedLLTs.push_back(Power);
        }
        Power=Power+1;
      }
      for(int i=0; i<int(GrabbedLLTs.size()); i++ )
      {
        long temp_LLT_Time = lltrigs[LLT].timestamp*20;
        if((temp_LLT_Time >= eventTimestamp) && int(temp_LLT_Time-eventTimestamp)<=fTrigWindow_ns)
        {
            //Loop over LLT we wanted to grab and check if type is correct
            if( std::find(fLLT.begin(), fLLT.end(), GrabbedLLTs[i]) != fLLT.end() )
            {
                GrabThisEvent=true;
                break;
            }
        }
      }
      if(GrabThisEvent) break;
    }
    if(GrabThisEvent) break;
  }


  return GrabThisEvent;
}

void TriggerWork::LLTEventGrabber::beginJob()
{
  // Implementation of optional member function here.
}

void TriggerWork::LLTEventGrabber::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerWork::LLTEventGrabber)
#endif