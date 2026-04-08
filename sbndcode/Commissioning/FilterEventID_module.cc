// #include <iostream>
// #include <stdlib.h>
// #include <string>
// #include <vector>
// #include <numeric>
// #include <getopt.h>
// #include <chrono>
// #include <float.h>


//some ROOT includes
// #include "TInterpreter.h"
// #include "TROOT.h"
// #include "TH1F.h"
// #include "TH2D.h"
// #include "TCanvas.h"
// #include "TTree.h"
// #include "TFile.h"
// #include "TStyle.h"
// #include "TSystem.h"
// #include "TGraph.h"
// #include "TFFTReal.h"
// #include "TTimeStamp.h"
// #include "TLatex.h"
// #include "TLine.h"


// art includes
#include <algorithm>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "canvas/Utilities/InputTag.h"
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
// #include "messagefacility/MessageLogger/MessageLogger.h"  
// #include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RawData/RawDigit.h"
// #include "lardataobj/RawData/RDTimeStamp.h"
// #include "lardataobj/RawData/raw.h"
// #include "larcore/Geometry/Geometry.h"
// #include "canvas/Persistency/Provenance/Timestamp.h"

// #include "art_root_io/TFileService.h"

// #include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"


namespace filt{

  class FilterEventID : public art::EDFilter {
    public:

      explicit FilterEventID(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& evt) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:

      int run;
      int subrun;
      int event;
      std::vector<int> filter_evts;
      std::vector<int> filter_runs;
    
  };


  FilterEventID::FilterEventID(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
    filter_evts = pset.get<std::vector<int> >("filterevts");
    filter_runs = pset.get<std::vector<int> >("filterruns", std::vector<int>{});

    if (filter_runs.size() != filter_evts.size()) {
      throw std::runtime_error("FilterEventID requires filterruns and filterevts to have the same length");
    }
  }

  void FilterEventID::reconfigure(fhicl::ParameterSet const& pset){
  }


  bool FilterEventID::filter(art::Event & evt){

    run = evt.id().run();
    subrun = evt.id().subRun();
    event = evt.id().event();

    bool is_in_list = false;
    for (std::size_t index = 0; index < filter_evts.size(); ++index) {
      if (filter_runs[index] == run && filter_evts[index] == event) {
        is_in_list = true;
        break;
      }
    }
    std::cout << "run: " << run << " event: " << event << " is_in_list: " << is_in_list << std::endl;
    return is_in_list;
    
  }

  void FilterEventID::beginJob() {
  }

  DEFINE_ART_MODULE(FilterEventID)

}
