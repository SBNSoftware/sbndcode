// Framework includes 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// ROOT includes 

// C++ includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <bitset>

const int DEFAULT_VALUE = -9999;

class testFilter : public art::EDFilter {
  public:
    explicit testFilter(fhicl::ParameterSet const& p);
    virtual bool filter(art::Event& e) override;
    void    reconfigure(fhicl::ParameterSet const& p);
    virtual ~testFilter() { }

  private:
    void ResetWireHitsVars(int n);
    art::ServiceHandle<art::TFileService> tfs;
    int max_hits;
    std::string fHitsModuleLabel; 

    // Wire hits variables
    int                 nhits;            
    std::vector<int>    hit_tpc;      
    std::vector<int>    hit_plane;      
    std::vector<int>    hit_wire;      
    std::vector<double> hit_peakT;    
  };

  testFilter::testFilter(fhicl::ParameterSet const& p): EDFilter{p} 
  {
    // read in the parameters from the .fcl file 
    this->reconfigure(p);
  }

  void testFilter::reconfigure(fhicl::ParameterSet const& p)
  {
    max_hits = p.get<int>("MaxHits", 50000);
    fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
  }

  bool testFilter::filter(art::Event& evt)
  { 
    bool KeepMe = true;

    // stuff for inside the event
    //int run = evt.run();
    int event = evt.id().event();

    // get nhits 
    art::Handle<std::vector<recob::Hit>> hitListHandle;
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle)) {
      art::fill_ptr_vector(hitlist, hitListHandle);
      nhits = hitlist.size();
      std::cout << "got data product,check" << std::endl;
    }
    else {
      std::cout << "Failed to get recob::Hit data product." << std::endl;
      nhits = 0;
    }
    if (nhits > max_hits) {
      std::cout << "Available hits are " << nhits 
                << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << max_hits << "hits." << std::endl;
      nhits = max_hits;
    }

    ResetWireHitsVars(nhits);

    for (int i = 0; i < nhits; ++i) {
      //std::cout << "pre-wire ID" << std::endl; 
      geo::WireID wireid = hitlist[i]->WireID();
      //std::cout << "post-wire ID" << std::endl;
      hit_plane[i] = wireid.Plane;
      //std::cout << "post-hitplane ID" << std::endl;
      hit_tpc[i] = wireid.TPC;
      hit_wire[i] = wireid.Wire;
      hit_peakT[i] = hitlist[i]->PeakTime();
      //std::cout << "post-hit_peakT" << std::endl;
      if (i<5){
        std::cout << "event #" << event << std::endl;
        std:: cout << "plane, tpc, wire, peakT: " << hit_plane[i] << ", " << hit_tpc[i] << ", " 
                   << hit_wire[i] << ", " <<  hit_peakT[i] << std::endl;
      }
    }
    return KeepMe; 
  }

  void testFilter::ResetWireHitsVars(int n) {
  hit_tpc.assign(n, DEFAULT_VALUE);
  hit_plane.assign(n, DEFAULT_VALUE);
  hit_wire.assign(n, DEFAULT_VALUE);
  hit_peakT.assign(n, DEFAULT_VALUE);
}
  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(testFilter)