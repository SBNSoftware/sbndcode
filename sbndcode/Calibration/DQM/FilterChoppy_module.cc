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
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

#include "art_root_io/TFileService.h"

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"


namespace filt{

  class FilterChoppy : public art::EDFilter {
    public:

      explicit FilterChoppy(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& evt) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

      std::vector<art::Ptr<raw::RDTimeStamp>> raw_timestamps_handle;

    private:

      int run;
      int subrun;
      int event;
      long delta_tmaxmin;
      std::vector<long> ts_vec;
      std::vector<long> ts_hist_vec;
      int npeak;
    
  };


  FilterChoppy::FilterChoppy(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
  }

  void FilterChoppy::reconfigure(fhicl::ParameterSet const& pset){
  }


  bool FilterChoppy::filter(art::Event & evt){

    raw_timestamps_handle.clear();
    run = evt.id().run();
    subrun = evt.id().subRun();
    event = evt.id().event();

    ts_vec.clear();
    ts_hist_vec.clear();
    npeak = 0;
    delta_tmaxmin = 0;

    // get the timestamps
    art::Handle<std::vector<raw::RDTimeStamp>> timestamp_handle;
    evt.getByLabel("daq", timestamp_handle);

    // exit if the data isn't present
    if (!timestamp_handle.isValid()) {
      std::cerr << "Error: missing timestamps with producer (" << "daq" << ")" << std::endl;
      return false;
    }
    art::fill_ptr_vector(raw_timestamps_handle, timestamp_handle);

    // collect TS
    unsigned timeindex = 0;
    for (auto const& timestamps: raw_timestamps_handle) {
      long this_ts = timestamps->GetTimeStamp();
      ts_vec.push_back(this_ts);
      timeindex++;
    }
    // find min, max TS values
    long ts_min = ts_vec[0];
    long ts_max = ts_vec[0];
    for (long ts : ts_vec) {
      if (ts < ts_min) {
        ts_min = ts;
      }
      if (ts > ts_max) {
        ts_max = ts;
      }
    }
    delta_tmaxmin = std::abs(ts_max - ts_min);

    // if (delta_tmaxmin > 100000){
    //   std::cout << "Event " << event << " has a delta_tmaxmin of " << delta_tmaxmin << std::endl;
    // }
    return delta_tmaxmin < 100000;


    // // Identify choppy events by looking at the # of peaks in the TS distribution
    // make binned hist
    // float nbin = 100.;
    // float bin_width = (ts_max - ts_min)/nbin;
    // for (int i = 0; i < 100; ++i) {
    //   float bin_low = ts_min + i*bin_width;
    //   float bin_high = ts_min + (i+1)*bin_width;
    //   int bin_count = 0;
    //   for (long ts : ts_vec) {
    //     if ((bin_low < ts) & (ts <= bin_high)) {
    //       bin_count += 1;
    //     }
    //   }
    //   ts_hist_vec.push_back(bin_count);
    // }
    // // find peak bins
    // for (int bc : ts_hist_vec) {
    //   //TODO: make threshold fcl parameter
    //   if (bc > 1000) {
    //     npeak += 1;
    //   }
    // }

    // // If the energy deposit within the beam time is greater than some limit then trigger the event
    // if (npeak > 1){
    //   std::cout << "Event " << event << " has " << npeak << " peaks" << std::endl;
    // }
    // return npeak < 2;
    
  }

  void FilterChoppy::beginJob() {
  }

  DEFINE_ART_MODULE(FilterChoppy)

}
