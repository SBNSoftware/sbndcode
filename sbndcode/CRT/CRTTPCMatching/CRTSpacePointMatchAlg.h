#ifndef CRTSPACEPOINTMATCHALG_H_SEEN
#define CRTSPACEPOINTMATCHALG_H_SEEN

///////////////////////////////////////////////
// CRTSpacePointMatchAlg.h
//
// Functions for CRT space point to TPC
// track matching
// T Brooks (tbrooks@fnal.gov), November 2018
// H Lay (h.lay@lancaster.ac.uk) February 2023
///////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"

namespace sbnd::crt {

  struct SPMatchCandidate
  {
    art::Ptr<CRTSpacePoint> thisSP;
    art::Ptr<recob::Track>  thisTrack;
    double                  time;
    double                  score;
    bool                    valid;

    SPMatchCandidate(const art::Ptr<CRTSpacePoint> _thisSP = art::Ptr<CRTSpacePoint>(), const art::Ptr<recob::Track> _thisTrack = art::Ptr<recob::Track>(),
                     const double _time = -std::numeric_limits<double>::max(), const double _score = -std::numeric_limits<double>::max(),
                     const bool _valid = false)
    {
      thisSP    = _thisSP;
      thisTrack = _thisTrack;
      time      = _time;
      score     = _score;
      valid     = _valid;
    }
  };


  class CRTSpacePointMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

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

      fhicl::Atom<double> DCALimit {
        Name("DCALimit"),
          Comment(""),
          200.
          };

      fhicl::Atom<double> MinTPCTrackLength {
        Name("MinTPCTrackLength"),
          Comment("Only consider TPC tracks with a length longer than this (cm)"),
          0.
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

      fhicl::Atom<double> PECut {
        Name("PECut"),
          Comment("Only consider CRTSpacePoints with PE values larger than this"),
          0.0
          };

      fhicl::Atom<double> MaxUncert {
        Name("MaxUncert"),
          Comment("Only consider CRTSpacePoints with position uncertainties below this value (cm)"),
          1000.0
          };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
          Comment("")
          };

      fhicl::Atom<art::InputTag> CRTSpacePointLabel {
        Name("CRTSpacePointLabel"),
          Comment("")
          };
    };

    CRTSpacePointMatchAlg(const Config& config);

    CRTSpacePointMatchAlg(const Config& config, geo::GeometryCore const *GeometryService);

    CRTSpacePointMatchAlg(const fhicl::ParameterSet& pset) :
    CRTSpacePointMatchAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTSpacePointMatchAlg(const fhicl::ParameterSet& pset, geo::GeometryCore const *GeometryService) :
    CRTSpacePointMatchAlg(fhicl::Table<Config>(pset, {})(), GeometryService) {}

    CRTSpacePointMatchAlg();

    ~CRTSpacePointMatchAlg();

    void reconfigure(const Config& config);

    SPMatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                             const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const art::Event &e);

    SPMatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                             const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs,
                                             const art::Event &e);

    SPMatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                             const std::pair<double, double> t0MinMax, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const int driftDirection,
                                             const art::Event &e);

    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const &detProp, const double startX,
                                           const double endX, const int driftDirection, const std::pair<double, double> xLimits);

    double DistOfClosestApproach(detinfo::DetectorPropertiesData const &detProp, geo::Point_t trackStart,
                                 const geo::Vector_t &trackDir, const art::Ptr<CRTSpacePoint> &crtSP,
                                 const int driftDirection, const double t0, const art::Event &e);

    std::pair<geo::Vector_t, geo::Vector_t> AverageTrackDirections(const art::Ptr<recob::Track> &track, const double frac);
    std::pair<geo::Vector_t, geo::Vector_t> TrackDirections(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                            const double frac, const double CRTtime, const int driftDirection);

  private:

    geo::GeometryCore const* fGeometryService;

    double fTrackDirectionFrac;
    double fDCALimit;
    double fTimeCorrection;
    double fMinTPCTrackLength;
    int    fDirMethod;
    bool   fDCAuseBox;
    bool   fDCAoverLength;
    double fPECut;
    double fMaxUncert;

    art::InputTag fTPCTrackLabel;
    art::InputTag fCRTSpacePointLabel;
  };
}
#endif
