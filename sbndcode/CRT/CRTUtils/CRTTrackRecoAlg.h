#ifndef CRTTRACKRECOALG_H_SEEN
#define CRTTRACKRECOALG_H_SEEN


///////////////////////////////////////////////
// CRTTrackRecoAlg.h
//
// Functions for CRT track reconstruction
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

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"

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

  class CRTTrackRecoAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> AverageHitDistance {
        Name("AverageHitDistance"),
        Comment("Distance to average hits over on same plane")
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("Distance to combine CRT hits into track")
      };

    };

    CRTTrackRecoAlg(const Config& config);

    CRTTrackRecoAlg(const fhicl::ParameterSet& pset) :
      CRTTrackRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrackRecoAlg(double aveHitDist, double distLim);

    ~CRTTrackRecoAlg();

    void reconfigure(const Config& config);

    // Function to make creating CRTTracks easier
    crt::CRTTrack FillCrtTrack(crt::CRTHit hit1, crt::CRTHit hit2, bool complete);

    // Function to average hits within a certain distance of each other
    std::vector<crt::CRTHit> AverageHits(std::vector<art::Ptr<crt::CRTHit>> hits);

    // Take a list of hits and find average parameters
    crt::CRTHit DoAverage(std::vector<art::Ptr<crt::CRTHit>> hits);

    // Create CRTTracks from list of hits
    std::vector<crt::CRTTrack> CreateTracks(std::vector<crt::CRTHit> hits);

    // Calculate the tagger crossing point of CRTTrack candidate
    TVector3 CrossPoint(crt::CRTHit hit, TVector3 start, TVector3 diff);


  private:

    geo::GeometryCore const* fGeometryService;
    art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    const geo::AuxDetGeometry* fAuxDetGeo;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

    double fAverageHitDistance;
    double fDistanceLimit;

    CRTHitRecoAlg hitAlg;

  };

}

#endif
