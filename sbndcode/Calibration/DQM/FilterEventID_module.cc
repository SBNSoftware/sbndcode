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
      virtual void endJob() override;
    private:

      int run;
      int subrun;
      int event;
      std::vector<int> filter_run;
      std::vector<int> filter_subrun;
      std::vector<int> filter_evts;
      std::vector<int> filter_slcid;
      std::vector<int> filter_pfpid;
      std::ofstream outfile;
  };


  FilterEventID::FilterEventID(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
    filter_evts = pset.get<std::vector<int> >("filterevts");
    filter_run = pset.get<std::vector<int> >("filterrun");
    filter_subrun = pset.get<std::vector<int> >("filtersubrun");
    filter_slcid = pset.get<std::vector<int> >("filterslcid");
    filter_pfpid = pset.get<std::vector<int> >("filterpfpid");
  }

  void FilterEventID::reconfigure(fhicl::ParameterSet const& pset){
  }


  bool FilterEventID::filter(art::Event & evt){

    // check that the filter lists are the same length
    if (filter_run.size() != filter_subrun.size() || filter_run.size() != filter_evts.size()){
      throw cet::exception("FilterEventID") << "Filter lists are not the same length";
    }

    run = evt.id().run();
    subrun = evt.id().subRun();
    event = evt.id().event();

    // bool is_in_list = std::find(filter_evts.begin(), filter_evts.end(), event) != filter_evts.end();
    // std::cout << "event: " << event << " is_in_list: " << is_in_list << std::endl;
    // return is_in_list;

    for (size_t i = 0; i < filter_run.size(); i++){
      // std::cout << "run: " << run << " subrun: " << subrun << " event: " << event << " filter_run: " << filter_run[i] << " filter_subrun: " << filter_subrun[i] << " filter_evts: " << filter_evts[i] << std::endl;
      if (run == filter_run[i] && subrun == filter_subrun[i] && event == filter_evts[i]){
        std::cout << "FOUND A MATCH" << std::endl;
        // write out to a txt file
        outfile << run << ", " << subrun << ", " << event << ", " << filter_slcid[i] << ", " << filter_pfpid[i] << std::endl;
        return true;
      }
    }
    return false;
    
  }

  void FilterEventID::beginJob() {
    outfile.open("filter_evts.txt");
    outfile << "run, subrun, event, slcID, id" << std::endl;
  }

  void FilterEventID::endJob() {
    outfile.close();
  }

  DEFINE_ART_MODULE(FilterEventID)

}
