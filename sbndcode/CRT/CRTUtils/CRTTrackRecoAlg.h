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
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

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

      fhicl::Atom<double> TimeLimit {
        Name("TimeLimit"),
        Comment("")
      };

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

    std::vector<std::vector<art::Ptr<sbn::crt::CRTHit>>> CreateCRTTzeros(std::vector<art::Ptr<sbn::crt::CRTHit>>);

    // Function to make creating CRTTracks easier
    sbn::crt::CRTTrack FillCrtTrack(sbn::crt::CRTHit hit1, sbn::crt::CRTHit hit2, bool complete);
    sbn::crt::CRTTrack FillCrtTrack(sbn::crt::CRTHit hit1, sbn::crt::CRTHit hit2, size_t nhits);

    // Function to average hits within a certain distance of each other
    std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> AverageHits(std::vector<art::Ptr<sbn::crt::CRTHit>> hits, std::map<art::Ptr<sbn::crt::CRTHit>, int> hitIds);
    std::vector<sbn::crt::CRTHit> AverageHits(std::vector<art::Ptr<sbn::crt::CRTHit>> hits);

    // Take a list of hits and find average parameters
    sbn::crt::CRTHit DoAverage(std::vector<art::Ptr<sbn::crt::CRTHit>> hits);

    // Create CRTTracks from list of hits
    std::vector<std::pair<sbn::crt::CRTTrack, std::vector<int>>> CreateTracks(std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> hits);
    std::vector<sbn::crt::CRTTrack> CreateTracks(std::vector<sbn::crt::CRTHit> hits);


  private:

    double fTimeLimit;
    double fAverageHitDistance;
    double fDistanceLimit;

    CRTHitRecoAlg hitAlg;
    CRTGeoAlg fCrtGeo;

  };

}

#endif
