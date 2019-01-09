#ifndef CRTTRACKMATCHALG_H_SEEN
#define CRTTRACKMATCHALG_H_SEEN


///////////////////////////////////////////////
// CRTTrackMatchAlg.h
//
// Functions for CRT TPC track matching
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"


namespace sbnd{

  struct RecoCRTTrack{
    int crtID;
    int tpc;
    TVector3 start; //[cm]
    TVector3 end; //[cm]
    double trueTime; // [us]
    bool complete;
  };

  class CRTTrackMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
    };

    CRTTrackMatchAlg(const Config& config);

    CRTTrackMatchAlg(const fhicl::ParameterSet& pset) :
      CRTTrackMatchAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrackMatchAlg();

    ~CRTTrackMatchAlg();

    void reconfigure(const Config& config);

    // Function to transform a CRTTrack into an expected reconstructed track
    std::vector<RecoCRTTrack> CrtToRecoTrack(crt::CRTTrack, int id);

    // Function to shift CRTTrack in X and work out how much is reconstructed
    std::vector<RecoCRTTrack> CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                                 int tpc, int id, double time, bool complete);

    // Function to calculate if a CRTTrack crosses the TPC volume
    bool CrossesTPC(crt::CRTTrack track);

  private:

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorProperties const* fDetectorProperties;
    detinfo::DetectorClocks const* fDetectorClocks;

  };

}

#endif
