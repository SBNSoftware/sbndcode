#ifndef CRTCOMMONUTILS_H_SEEN
#define CRTCOMMONUTILS_H_SEEN

///////////////////////////////////////////////
// CRTCommonUtils.h
//
// Common functions for CRT reconstruction
///////////////////////////////////////////////

#include "canvas/Persistency/Common/Ptr.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

#include <string>

namespace sbnd::crt {

  namespace CRTCommonUtils{

    // Returns the CRT tagger index given the tagger name
    enum CRTTagger GetTaggerEnum(const std::string tagger);

    // Returns the CRT tagger name given the tagger index
    std::string GetTaggerName(const CRTTagger tagger);

    // Returns the coordinate constrained by virtue of the tagger's position
    enum CoordSet GetTaggerDefinedCoordinate(const CRTTagger tagger);

    // Returns the coordinate direction of the strip's width given the tagger & local orientation
    enum CoordSet GetStripWidthGlobalCoordinate(const CRTTagger tagger, const uint16_t orientation);

    // Returns whether the tagger is one of the top taggers, useful for track building
    bool IsTopTagger(const CRTTagger tagger);

    // Returns whether a set of three taggers includes both of the top taggers, useful for track building
    bool CoverTopTaggers(const CRTTagger tagger1, const CRTTagger tagger2, const CRTTagger tagger3);

    // Returns the intersection between a cuboid and an infinite line (TPC & CRTTrack)
    bool CuboidIntersection(const geo::Point_t &min, const geo::Point_t &max, const geo::Point_t &start, const geo::Point_t &end,
                            geo::Point_t &entry, geo::Point_t &exit);

    // Returns the intersection between an infinte line and an axis-aligned plane
    double LinePlaneIntersection(const geo::Point_t &start ,const geo::Vector_t &dir, const CoordSet plane, const double value);

    // Returns whether a point lies inside a square plane of a cuboid defined by its minimum and maximum points
    bool IsInsideRectangle(const geo::Point_t &min, const geo::Point_t &max, const geo::Point_t &p, const CoordSet plane);

    // Returns a simple distance of closest approach between an infinite track and CRTSpacePoint
    double SimpleDCA(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Vector_t &direction);

    // Returns the minimum distance from an infinite tracks to a CRTSpacePoint assuming its a 2D rectangle
    double DistToCRTSpacePoint(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Point_t &end, const CRTTagger tagger);

    // Returns the distance between an infinite line and a segment
    double LineSegmentDistance(const geo::Point_t &start1, const geo::Point_t &end1, const geo::Point_t &start2, const geo::Point_t &end2);
  }
}

#endif
