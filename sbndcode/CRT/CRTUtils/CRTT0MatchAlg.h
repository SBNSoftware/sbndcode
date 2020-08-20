#ifndef CRTT0MATCHALG_H_SEEN
#define CRTT0MATCHALG_H_SEEN


///////////////////////////////////////////////
// CRTT0MatchAlg.h
//
// Functions for CRT t0 matching
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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"

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

  class CRTT0MatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> MinTrackLength {
        Name("MinTrackLength"),
        Comment("")
      };

      fhicl::Atom<double> TrackDirectionFrac {
        Name("TrackDirectionFrac"),
        Comment("")
      };

      fhicl::Atom<int> TSMode {
        Name("TSMode"),
        Comment(""),
        1 /* Value for SBND */
      }; 

      fhicl::Atom<double> TimeCorrection {
        Name("TimeCorrection"),
        Comment(""),
        0.
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("")
      };

    };

    CRTT0MatchAlg(const Config& config);
    CRTT0MatchAlg(const Config& config, geo::GeometryCore const *GeometryService);

    CRTT0MatchAlg(const fhicl::ParameterSet& pset) :
      CRTT0MatchAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTT0MatchAlg(const fhicl::ParameterSet& pset, geo::GeometryCore const *GeometryService) :
      CRTT0MatchAlg(fhicl::Table<Config>(pset, {})(), GeometryService) {}

    CRTT0MatchAlg();

    void reconfigure(const Config& config);

    // Utility function that determines the possible x range of a track
    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
                                           double startX, double endX, int driftDirection, std::pair<double, double> xLimits);

    // Calculate the distance of closest approach (DCA) between the end of a track and a crt hit
    double DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                 TVector3 trackPos, TVector3 trackDir, crt::CRTHit crtHit, int driftDirection, double t0);

    std::pair<TVector3, TVector3> TrackDirectionAverage(recob::Track track, double frac);
    std::pair<TVector3, TVector3> TrackDirectionAverageFromPoints(recob::Track track, double frac);

    // Return the closest CRT hit to a TPC track and the DCA
    std::pair<crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
                                                 recob::Track tpcTrack, std::pair<double, double> t0MinMax, std::vector<sbnd::crt::CRTHit> crtHits, int driftDirection);
    std::pair<crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
                                                 recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event);
    std::pair<crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
                                                 recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTHit> crtHits);

    // Match track to T0 from CRT hits
    double T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                         recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event);
    double T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                         recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTHit> crtHits);

    // Match track to T0 from CRT hits, also return the DCA
    std::pair<double, double> T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                                                  recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event);
    std::pair<double, double> T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                                                  recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTHit> crtHits);


  private:

    geo::GeometryCore const* fGeometryService;

    double fMinTrackLength;
    double fTrackDirectionFrac;
    double fDistanceLimit;
    int fTSMode;
    double fTimeCorrection;

    art::InputTag fTPCTrackLabel;

  };

}

#endif
