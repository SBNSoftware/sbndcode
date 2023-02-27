#ifndef CRTTRACKMATCHALG_H_SEEN
#define CRTTRACKMATCHALG_H_SEEN

///////////////////////////////////////////////
// CRTTrackMatchAlg.h
//
// Functions for CRT TPC track matching
// T Brooks (tbrooks@fnal.gov), November 2018
// H Lay (h.lay@lancaster.ac.uk) February 2023
///////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"

namespace sbnd::crt {

  struct MatchCandidate {
    art::Ptr<CRTTrack> thisTrack;
    double             time;
    double             score;
    bool               valid;

    MatchCandidate(const art::Ptr<CRTTrack> _thisTrack = art::Ptr<CRTTrack>(), const double _time = -std::numeric_limits<double>::max(),
                   const double _score = -std::numeric_limits<double>::max(), const bool _valid = false)
    {
      thisTrack = _thisTrack;
      time      = _time;
      score     = _score;
      valid     = _valid;
    }
  };

  class CRTTrackMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> MaxAngleDiff {
        Name("MaxAngleDiff"),
          Comment("")
          };

      fhicl::Atom<double> MaxDCA {
        Name("MaxDCA"),
          Comment("")
          };

      fhicl::Atom<double> MaxScore {
        Name("MaxScore"),
          Comment("")
          };

      fhicl::Atom<std::string> SelectionMetric {
        Name("SelectionMetric"),
          Comment("")
          };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
          Comment("")
          };
    };

    CRTTrackMatchAlg(const Config& config);

    CRTTrackMatchAlg(const Config& config, geo::GeometryCore const* GeometryService);

    CRTTrackMatchAlg(const fhicl::ParameterSet& pset) :
    CRTTrackMatchAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrackMatchAlg(const fhicl::ParameterSet& pset, geo::GeometryCore const* GeometryService) :
    CRTTrackMatchAlg(fhicl::Table<Config>(pset, {})(), GeometryService) {}

    CRTTrackMatchAlg();

    ~CRTTrackMatchAlg();

    void reconfigure(const Config& config);

    bool TPCIntersection(const geo::TPCGeo &tpcGeo, const art::Ptr<CRTTrack> &track, geo::Point_t &entry, geo::Point_t &exit);

    MatchCandidate GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e);

    MatchCandidate GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks);

    std::vector<art::Ptr<CRTTrack>> AllPossibleCRTTracks(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                         const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e);

    std::vector<art::Ptr<CRTTrack>> AllPossibleCRTTracks(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                         const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks);

    MatchCandidate ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxDCA);

    MatchCandidate ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks, const double maxDCA);

    MatchCandidate ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                        const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxAngle);

    MatchCandidate ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                        const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks, const double maxAngle);

    MatchCandidate ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e);

    MatchCandidate ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks);

    double AngleBetweenTracks(const art::Ptr<recob::Track> &tpcTrack, const art::Ptr<CRTTrack> &crtTrack);

    double AveDCABetweenTracks(const art::Ptr<recob::Track> &tpcTrack, const art::Ptr<CRTTrack> &crtTrack, const double shift);

  private:

    geo::GeometryCore const* fGeometryService;

    double      fMaxAngleDiff;
    double      fMaxDCA;
    double      fMaxScore;
    std::string fSelectionMetric;

    art::InputTag fTPCTrackLabel;
  };
}

#endif
