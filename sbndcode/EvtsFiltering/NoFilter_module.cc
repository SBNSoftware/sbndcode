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

  class NoFilter : public art::EDFilter {
    public:

      explicit NoFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& evt) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:
    
  };


  NoFilter::NoFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
  }

  void NoFilter::reconfigure(fhicl::ParameterSet const& pset){
  }


  bool NoFilter::filter(art::Event & evt){
    return true;
  }

  void NoFilter::beginJob() {
  }

  DEFINE_ART_MODULE(NoFilter)

}
