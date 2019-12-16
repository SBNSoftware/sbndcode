#ifndef OnlinePMTAna_h
#define OnlinePMTAna_h
////////////////////////////////////////////////////////////////////////
// Class:       OnlinePMTAna
// Plugin Type: analyzer (art v3_03_01)
// File:        OnlinePMTAna.h
//
// Generated at Sun Dec 15 20:17:27 2019 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <TTree.h>

namespace sbnom {
  class OnlinePMTAna;
}

namespace sbnom { 

class OnlinePMTAna : public art::EDAnalyzer {
public:
  explicit OnlinePMTAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OnlinePMTAna(OnlinePMTAna const&) = delete;
  OnlinePMTAna(OnlinePMTAna&&) = delete;
  OnlinePMTAna& operator=(OnlinePMTAna const&) = delete;
  OnlinePMTAna& operator=(OnlinePMTAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  std::string fWaveformTag;
  TTree *fOutput;
  unsigned fNWaveforms;

};

} // end namespace sbnom

#endif /* OnlinePMTAna_h */
