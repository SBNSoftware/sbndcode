#ifndef CRTSPACEPOINTMATCHALG_H_SEEN
#define CRTSPACEPOINTMATCHALG_H_SEEN

///////////////////////////////////////////////
// CRTSpacePointMatchAlg.h
//
// Functions for CRT t0 matching
// T Brooks (tbrooks@fnal.gov), November 2018
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

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"

namespace sbnd::crt {

  struct  MatchCandidate {
    art::Ptr<CRTSpacePoint> thisSP;
    double                  t0;
    double                  dca;
    double                  extrapLen;

    MatchCandidate(const art::Ptr<CRTSpacePoint> _thisSP, const double _t0,
                   const double _dca, const double _extrapLen)
    {
      thisSP    = _thisSP;
      t0        = _t0;
      dca       = _dca;
      extrapLen = _extrapLen;
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

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
          Comment(""),
          100.
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

    MatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                           const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const art::Event &e);

    MatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                           const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs);

    MatchCandidate GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                           const std::pair<double, double> t0MinMax, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const int driftDirection);

    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const &detProp, const double startX,
                                           const double endX, const int driftDirection, const std::pair<double, double> xLimits);

    double DistOfClosestApproach(detinfo::DetectorPropertiesData const &detProp, geo::Point_t trackStart,
                                 const geo::Vector_t &trackDir, const art::Ptr<CRTSpacePoint> &crtSP,
                                 const int driftDirection, const double t0);

    std::pair<geo::Vector_t, geo::Vector_t> AverageTrackDirections(const art::Ptr<recob::Track> &track, const double frac);
    std::pair<geo::Vector_t, geo::Vector_t> TrackDirections(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                            const double frac, const double CRTtime, const int driftDirection);

  private:

    geo::GeometryCore const* fGeometryService;

    double fTrackDirectionFrac;
    double fDistanceLimit;
    double fTimeCorrection;
    int    fDirMethod;
    bool   fDCAuseBox;
    bool   fDCAoverLength;
    double fPECut;
    double fMaxUncert;

    art::InputTag fTPCTrackLabel;
  };
}
#endif
