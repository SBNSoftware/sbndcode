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
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
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

  class CRTTrackMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> MaxAngleDiff {
        Name("MaxAngleDiff"),
        Comment("")
      };

      fhicl::Atom<double> MaxDistance {
        Name("MaxDistance"),
        Comment("")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("")
      };

      fhicl::Atom<bool> MinimizeAngle {
        Name("MinimizeAngle"),
        Comment("")
      };

    };

    CRTTrackMatchAlg(const Config& config);
    CRTTrackMatchAlg(const Config& config, geo::GeometryCore const* GeometryService,  detinfo::DetectorProperties const* DetectorProperties);

    CRTTrackMatchAlg(const fhicl::ParameterSet& pset) :
      CRTTrackMatchAlg(fhicl::Table<Config>(pset, {})()) {}
    CRTTrackMatchAlg(const fhicl::ParameterSet& pset, geo::GeometryCore const* GeometryService,  detinfo::DetectorProperties const* DetectorProperties) :
      CRTTrackMatchAlg(fhicl::Table<Config>(pset, {})(), GeometryService, DetectorProperties) {}

    CRTTrackMatchAlg();

    ~CRTTrackMatchAlg();

    void reconfigure(const Config& config);

    // Calculate intersection between CRT track and TPC
    std::pair<TVector3, TVector3> TpcIntersection(const geo::TPCGeo& tpcGeo, crt::CRTTrack track);

    // Function to calculate if a CRTTrack crosses the TPC volume
    bool CrossesTPC(crt::CRTTrack track);

    // Function to calculate if a CRTTrack crosses the TPC volume
    bool CrossesAPA(crt::CRTTrack track);

    double T0FromCRTTracks(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event);
    double T0FromCRTTracks(recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<crt::CRTTrack> crtTracks);

    // Find the closest valid matching CRT track ID
    int GetMatchedCRTTrackId(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event);
    int GetMatchedCRTTrackId(recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<crt::CRTTrack> crtTracks);

    // Get all CRT tracks that cross the right TPC within an allowed time
    std::vector<crt::CRTTrack> AllPossibleCRTTracks(recob::Track tpcTrack, 
                                                    std::vector<crt::CRTTrack> crtTracks, 
                                                    const art::Event& event); 

    std::vector<crt::CRTTrack> AllPossibleCRTTracks(recob::Track tpcTrack, 
                                                    std::vector<art::Ptr<recob::Hit>> hits,
                                                    std::vector<crt::CRTTrack> crtTracks);

    // Find the closest matching crt track by angle between tracks within angle and DCA limits
    std::pair<crt::CRTTrack, double> ClosestCRTTrackByAngle(recob::Track tpcTrack, 
                                                            std::vector<crt::CRTTrack> crtTracks, 
                                                            const art::Event& event,
                                                            double minDCA = 0.); 
    std::pair<crt::CRTTrack, double> ClosestCRTTrackByAngle(recob::Track tpcTrack, 
                                                            std::vector<art::Ptr<recob::Hit>> hits, 
                                                            std::vector<crt::CRTTrack> crtTracks, 
                                                            double minDCA = 0.); 
    // Find the closest matching crt track by average DCA between tracks within angle and DCA limits
    std::pair<crt::CRTTrack, double> ClosestCRTTrackByDCA(recob::Track tpcTrack, 
                                                          std::vector<crt::CRTTrack> crtTracks, 
							                                            const art::Event& event,
                                                          double minAngle = 0.); 
    std::pair<crt::CRTTrack, double> ClosestCRTTrackByDCA(recob::Track tpcTrack, 
                                                          std::vector<art::Ptr<recob::Hit>> hits, 
                                                          std::vector<crt::CRTTrack> crtTracks, 
                                                          double minAngle = 0.); 

    // Calculate the angle between tracks assuming start is at the largest Y
    double AngleBetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack);

    // Calculate the average DCA between tracks
    double AveDCABetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack, double shift);
    double AveDCABetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack, const art::Event& event);
    double AveDCABetweenTracks(recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, crt::CRTTrack crtTrack);

  private:

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorProperties const* fDetectorProperties;

    CRTBackTracker fCrtBackTrack;

    double fMaxAngleDiff;
    double fMaxDistance;
    bool fMinimizeAngle;

    art::InputTag fTPCTrackLabel;

  };

}

#endif
