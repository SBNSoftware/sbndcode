#ifndef CRTSPACEPOINTMATCHALG_H_SEEN
#define CRTSPACEPOINTMATCHALG_H_SEEN


///////////////////////////////////////////////
// CRTSpacePointMatchAlg.h
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
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
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

#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
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


  struct  matchCand {
    sbnd::crt::CRTSpacePoint thisSP;
    double t0;
    double dca;
    double extrapLen;
  };


  class CRTSpacePointMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> MinTrackLength {
        Name("MinTrackLength"),
          Comment(""),
          20.0
          };

      fhicl::Atom<double> TrackDirectionFrac {
        Name("TrackDirectionFrac"),
          Comment(""),
          0.5
          };

      fhicl::Atom<double> TimeCorrection {
        Name("TimeCorrection"),
          Comment(""),
          0.
          };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
          Comment(""),
          100.
          };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
          Comment("")
          };


      fhicl::Atom<int> DirMethod {
        Name("DirMethod"),
          Comment("1=endpoints (default), 2=average;  must use endpoints if applying SCE position corrections"),
          1
          };

      fhicl::Atom<bool> DCAuseBox {
        Name("DCAuseBox"),
          Comment("false = distance to point (default), true = distance to box edge"),
          false
          };

      fhicl::Atom<bool> DCAoverLength {
        Name("DCAoverLength"),
          Comment("false = use DCA to select closest CRTSpacePoint (default), true = use DCA/extrapolation_length"),
          false
          };

      fhicl::Atom<double> DoverLLimit {
        Name("DoverLLimit"),
          Comment(""),
          1.
          };

      fhicl::Atom<double> PEcut {
        Name("PEcut"),
          Comment("Only consider CRTSpacePoints with PE values larger than this"),
          0.0
          };

      fhicl::Atom<double> MaxUncert {
        Name("MaxUncert"),
          Comment("Only consider CRTSpacePoints with position uncertainties below this value (cm)"),
          1000.0
          };
    };


    CRTSpacePointMatchAlg(const Config& config);
    CRTSpacePointMatchAlg(const Config& config, geo::GeometryCore const *GeometryService);

  CRTSpacePointMatchAlg(const fhicl::ParameterSet& pset) :
    CRTSpacePointMatchAlg(fhicl::Table<Config>(pset, {})()) {}
    
  CRTSpacePointMatchAlg(const fhicl::ParameterSet& pset, geo::GeometryCore const *GeometryService) :
    CRTSpacePointMatchAlg(fhicl::Table<Config>(pset, {})(), GeometryService) {}
    
    CRTSpacePointMatchAlg();
    
    void reconfigure(const Config& config);
    


    // Utility function that determines the possible x range of a track
    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
                                           double startX, double endX, int driftDirection, std::pair<double, double> xLimits);

    // Calculate the distance of closest approach (DCA) between the end of a track and a crt hit
    double DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                 geo::Point_t trackPos, geo::Vector_t trackDir, sbnd::crt::CRTSpacePoint crtSP, int driftDirection, double t0);

    std::pair<geo::Vector_t, geo::Vector_t> TrackDirectionAverage(recob::Track track, double frac);
    std::pair<geo::Vector_t, geo::Vector_t> TrackDirection(detinfo::DetectorPropertiesData const& detProp,recob::Track track, double frac, double CRTtime, int driftDirection);
    std::pair<geo::Vector_t, geo::Vector_t> TrackDirectionAverageFromPoints(recob::Track track, double frac);

    // Return the closest CRTSpacePoint to a TPC track and the DCA
    matchCand GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                      recob::Track tpcTrack, std::pair<double, double> t0MinMax, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints, int driftDirection);

    matchCand GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                      recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints, const art::Event& event);

    matchCand GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                      recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints);

    // Match track to T0 from CRTSpacePoints
    double T0FromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints, const art::Event& event);
    double T0FromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints);
    
    // Match track to T0 from CRTSpacePoints, also return the DCA
    std::pair<double, double>  T0AndDCAFromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                          recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints, const art::Event& event);
    std::pair<double, double>  T0AndDCAFromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                          recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSpacePoints);
 

  private:

    geo::GeometryCore const* fGeometryService;

    double fMinTrackLength;
    double fTrackDirectionFrac;
    double fDistanceLimit;
    double fTimeCorrection;
    int fDirMethod;
    bool fDCAuseBox;
    bool fDCAoverLength;
    double fDoverLLimit;
    double fPEcut;
    double fMaxUncert;

    art::InputTag fTPCTrackLabel;

  };


}
#endif
