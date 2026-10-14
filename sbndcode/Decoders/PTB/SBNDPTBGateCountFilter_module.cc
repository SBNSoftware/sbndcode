////////////////////////////////////////////////////////////////////////
// Class:       SBNDPTBGateCountFilter
// Plugin Type: art::EDFilter
// File:        SBNDPTBGateCountFilter_module.cc
// Purpose:     Filter events from the BNBLight stream with
//              corrupted POT normalization data
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <iostream>
#include <string>
#include <vector>

// art includes
//#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"  
#include <bitset>


#include "sbndcode/Decoders/PTB/sbndptb.h"



class SBNDPTBGateCountFilter : public art::EDFilter {
  public:

    explicit SBNDPTBGateCountFilter(fhicl::ParameterSet const & pset);
    virtual bool filter(art::Event& evt) override;

  private:

    // Fhicl parameters
    std::string fPTBInputLabel;
    int fVerbose;

    const int bnbLightHLT_id = 2; // HLT2
    const int bnbLowLightLLT_id = 6; // HLT6

    int run;
    int subrun;
    int event;
  
};


SBNDPTBGateCountFilter::SBNDPTBGateCountFilter(fhicl::ParameterSet const & pset)
: EDFilter(pset),
fPTBInputLabel(pset.get<std::string>("PTBInputLabel"))
{
  std::cout << "SBNDPTBGateCountFilter configured with PTBInputLabel: " << fPTBInputLabel << std::endl;
}


bool SBNDPTBGateCountFilter::filter(art::Event & evt){

  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  if(fVerbose > 1){
    std::cout << "SBNDPTBGateCountFilter::Processing event: " << run << ":" << subrun << ":" << event << std::endl;
  }

  // Get PTB decoded data
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  evt.getByLabel(fPTBInputLabel, ptbHandle);
  if(!ptbHandle.isValid()){
    // If the handle is not valid, thow an exception
    throw art::Exception(art::errors::ProductNotFound)
      << "Could not find PTB data with label: " << fPTBInputLabel << std::endl;
  }
  std::vector<art::Ptr<raw::ptb::sbndptb>> ptbVec;
  art::fill_ptr_vector(ptbVec, ptbHandle);

  // Loop over PTB data and check for HLT2 and HLT6 triggers
  for (auto const& ptb : ptbVec) {

      // Loop over HLT trigger words in the PTB data
      for (unsigned i = 0; i < ptb->GetNHLTriggers(); ++i) {
          auto const& hlt = ptb->GetHLTrigger(i);
          int hlt_word = hlt.trigger_word;
          uint64_t timestamp = hlt.timestamp;

          if(fVerbose > 1){
            std::cout << "  HLT word: " << std::bitset<32>(hlt_word)
                    << ", timestamp: " << timestamp << std::endl;
          }

          // Look for the first HLT2 (BNBLight) instance in the event
          if (hlt_word & (1 << bnbLightHLT_id)) {

              // If first HLT2 instance in the event — decision point
              // Check if HLT6 (BNBLowLight) is also present in the same trigger word
              bool hasHLT6 = hlt_word & (1 << bnbLowLightLLT_id);

              if( fVerbose > 0 && hasHLT6){
                std::cout << "  FILTERING OUT Event " << run << ":" << subrun << ":" << event
                          << " - first HLT2 instance coincident with HLT6." << std::endl;
              }
              
              return !hasHLT6;  // false filters out, true keeps event
          }
      }
  }

  return true; // return true for other data streams (i.e. not BNBLight HLT2)
  
}

DEFINE_ART_MODULE(SBNDPTBGateCountFilter)
