#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <getopt.h>
#include <chrono>
#include <float.h>


//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TFFTReal.h"
#include "TTimeStamp.h"
#include "TLatex.h"
#include "TLine.h"


// art includes
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
#include "messagefacility/MessageLogger/MessageLogger.h"  
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

#include "art_root_io/TFileService.h"

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"


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
    
  };


  FilterEventID::FilterEventID(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
    filter_evts = pset.get<std::vector<int> >("filterevts");
  }

  void FilterEventID::reconfigure(fhicl::ParameterSet const& pset){
  }


  bool FilterEventID::filter(art::Event & evt){

    run = evt.id().run();
    subrun = evt.id().subRun();
    event = evt.id().event();

    bool is_in_list = std::find(filter_evts.begin(), filter_evts.end(), event) != filter_evts.end();
    std::cout << "event: " << event << " is_in_list: " << is_in_list << std::endl;
    return is_in_list;
    
  }

  void FilterEventID::beginJob() {
  }

  DEFINE_ART_MODULE(FilterEventID)

}
