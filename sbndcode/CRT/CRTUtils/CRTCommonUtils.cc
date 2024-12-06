#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

namespace sbnd::crt {

  enum CRTTagger CRTCommonUtils::GetTaggerEnum(const std::string tagger)
  {
    if      (tagger == "volTaggerBot_0"     ) return kBottomTagger;
    else if (tagger == "volTaggerSouth_0"   ) return kSouthTagger;
    else if (tagger == "volTaggerNorth_0"   ) return kNorthTagger;
    else if (tagger == "volTaggerWest_0"    ) return kWestTagger;
    else if (tagger == "volTaggerEast_0"    ) return kEastTagger;
    else if (tagger == "volTaggerTopLow_0"  ) return kTopLowTagger;
    else if (tagger == "volTaggerTopHigh_0" ) return kTopHighTagger;
    else {
      mf::LogWarning("CRTCommonUtils") << "CRT tagger unknown: " << tagger << std::endl;
      return kUndefinedTagger;
    }
  }

  std::string CRTCommonUtils::GetTaggerName(const CRTTagger tagger)
  {
    if      (tagger == kBottomTagger)  return "volTaggerBot_0";
    else if (tagger == kSouthTagger)   return "volTaggerSouth_0";
    else if (tagger == kNorthTagger)   return "volTaggerNorth_0";
    else if (tagger == kWestTagger)    return "volTaggerWest_0";
    else if (tagger == kEastTagger)    return "volTaggerEast_0";
    else if (tagger == kTopLowTagger)  return "volTaggerTopLow_0";
    else if (tagger == kTopHighTagger) return "volTaggerTopHigh_0";
    else {
      mf::LogWarning("CRTCommonUtils") << "CRT tagger unknown: " << tagger << std::endl;
      return "";
    }
  }

  enum CoordSet CRTCommonUtils::GetTaggerDefinedCoordinate(const CRTTagger tagger)
  {
    switch(tagger)
      {
      case kWestTagger:
      case kEastTagger:
        return kX;
      case kBottomTagger:
      case kTopLowTagger:
      case kTopHighTagger:
        return kY;
      case kSouthTagger:
      case kNorthTagger:
        return kZ;
      case kUndefinedTagger:
        return kUndefinedSet;
      }

    return kUndefinedSet;
  }

  enum CoordSet CRTCommonUtils::GetStripWidthGlobalCoordinate(const CRTTagger tagger, const uint16_t orientation)
  {
    if((tagger == kBottomTagger && orientation == 1) ||
       ((tagger == kTopLowTagger || tagger == kTopHighTagger) && orientation == 0) ||
       ((tagger == kSouthTagger || tagger == kNorthTagger) && orientation == 0))
      return kX;
    else if(((tagger == kSouthTagger || tagger == kNorthTagger) && orientation == 1) ||
            ((tagger == kWestTagger || tagger == kEastTagger) && orientation == 1))
      return kY;
    else if((tagger == kBottomTagger && orientation == 0) ||
            ((tagger == kTopLowTagger || tagger == kTopHighTagger) && orientation == 1) ||
            ((tagger == kWestTagger || tagger == kEastTagger) && orientation == 0))
      return kZ;

    return kUndefinedSet;
  }

  bool CRTCommonUtils::IsTopTagger(const CRTTagger tagger)
  {
    return tagger == kTopLowTagger || tagger == kTopHighTagger;
  }

  bool CRTCommonUtils::CoverTopTaggers(const CRTTagger tagger1, const CRTTagger tagger2, const CRTTagger tagger3)
  {
    return ( tagger1 == kTopLowTagger || tagger2 == kTopLowTagger || tagger3 == kTopLowTagger)
      && ( tagger1 == kTopHighTagger || tagger2 == kTopHighTagger || tagger3 == kTopHighTagger);
  }

  bool CRTCommonUtils::CuboidIntersection(const geo::Point_t &min, const geo::Point_t &max, const geo::Point_t &start, const geo::Point_t &end,
                                          geo::Point_t &entry, geo::Point_t &exit)
  {
    const geo::Vector_t dir = (end - start);

    std::vector<std::pair<double, CoordSet>> distances;

    distances.emplace_back(LinePlaneIntersection(start, dir, kX, min.X()), kX);
    distances.emplace_back(LinePlaneIntersection(start, dir, kX, max.X()), kX);
    distances.emplace_back(LinePlaneIntersection(start, dir, kY, min.Y()), kY);
    distances.emplace_back(LinePlaneIntersection(start, dir, kY, max.Y()), kY);
    distances.emplace_back(LinePlaneIntersection(start, dir, kZ, min.Z()), kZ);
    distances.emplace_back(LinePlaneIntersection(start, dir, kZ, max.Z()), kZ);

    std::vector<double> chosen_distances;

    for(auto const& [k, plane] : distances)
      {
        const geo::Point_t intersection = start + k * dir;

        if(IsInsideRectangle(min, max, intersection, plane))
          chosen_distances.push_back(k);
      }

    if(chosen_distances.size() == 0)
      return false;
    else if(chosen_distances.size() == 2)
      {
        entry = start + chosen_distances[0] * dir;
        exit  = start + chosen_distances[1] * dir;

        if(chosen_distances[1] < chosen_distances[0])
          std::swap(entry, exit);

        return true;
      }
    else
      return false;
  }

  double CRTCommonUtils::LinePlaneIntersection(const geo::Point_t &start ,const geo::Vector_t &dir, const CoordSet plane, const double value)
  {
    if(plane == kX)
      return (value - start.X()) / dir.X();
    else if(plane == kY)
      return (value - start.Y()) / dir.Y();
    else if(plane == kZ)
      return (value - start.Z()) / dir.Z();

    std::cout << "Well this is disconcerting..." << std::endl;
    return -std::numeric_limits<double>::max();
  }

  bool CRTCommonUtils::IsInsideRectangle(const geo::Point_t &min, const geo::Point_t &max, const geo::Point_t &p, const CoordSet plane)
  {
    if(plane == kX)
      return (p.Y() >= min.Y() && p.Y() <= max.Y())
        && (p.Z() >= min.Z() && p.Z() <= max.Z());
    else if(plane == kY)
      return (p.X() >= min.X() && p.X() <= max.X())
        && (p.Z() >= min.Z() && p.Z() <= max.Z());
    else if(plane == kZ)
      return (p.X() >= min.X() && p.X() <= max.X())
        && (p.Y() >= min.Y() && p.Y() <= max.Y());

    std::cout << "Even more disconcerting..." << std::endl;
    return false;
  }

  double CRTCommonUtils::SimpleDCA(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Vector_t &direction)
  {
    const geo::Point_t pos = sp->Pos();
    const geo::Point_t end = start + direction;

    const double denominator = direction.R();
    const double numerator   = (pos - start).Cross(pos - end).R();

    return numerator / denominator;
  }

  double CRTCommonUtils::DistToCRTSpacePoint(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Point_t &end, const CRTTagger tagger)
  {
    const geo::Point_t min = sp->Pos() - geo::Vector_t(sp->Err());
    const geo::Point_t max = sp->Pos() + geo::Vector_t(sp->Err());

    geo::Point_t entry, exit;

    if(CuboidIntersection(min, max, start, end, entry, exit))
      return 0;

    const CoordSet constrainedPlane = GetTaggerDefinedCoordinate(tagger);

    geo::Point_t vertex1, vertex2, vertex3, vertex4;

    if(constrainedPlane == kX)
      {
        vertex1 = {sp->X(), sp->Y() - sp->YErr(), sp->Z() - sp->ZErr()};
        vertex2 = {sp->X(), sp->Y() + sp->YErr(), sp->Z() - sp->ZErr()};
        vertex3 = {sp->X(), sp->Y() - sp->YErr(), sp->Z() + sp->ZErr()};
        vertex4 = {sp->X(), sp->Y() + sp->YErr(), sp->Z() + sp->ZErr()};
      }
    else if(constrainedPlane == kY)
      {
        vertex1 = {sp->X() - sp->XErr(), sp->Y(), sp->Z() - sp->ZErr()};
        vertex2 = {sp->X() + sp->XErr(), sp->Y(), sp->Z() - sp->ZErr()};
        vertex3 = {sp->X() - sp->XErr(), sp->Y(), sp->Z() + sp->ZErr()};
        vertex4 = {sp->X() + sp->XErr(), sp->Y(), sp->Z() + sp->ZErr()};
      }
    else if(constrainedPlane == kZ)
      {
        vertex1 = {sp->X() - sp->XErr(), sp->Y() - sp->YErr(), sp->Z()};
        vertex2 = {sp->X() + sp->XErr(), sp->Y() - sp->YErr(), sp->Z()};
        vertex3 = {sp->X() - sp->XErr(), sp->Y() + sp->YErr(), sp->Z()};
        vertex4 = {sp->X() + sp->XErr(), sp->Y() + sp->YErr(), sp->Z()};
      }

    double dist1 = LineSegmentDistance(vertex1, vertex2, start, end);
    double dist2 = LineSegmentDistance(vertex1, vertex3, start, end);
    double dist3 = LineSegmentDistance(vertex4, vertex2, start, end);
    double dist4 = LineSegmentDistance(vertex4, vertex3, start, end);

    return std::min({dist1, dist2, dist3, dist4});
  }

  double CRTCommonUtils::LineSegmentDistance(const geo::Point_t &start1, const geo::Point_t &end1, const geo::Point_t &start2, const geo::Point_t &end2)
  {
    const double smallNum = std::numeric_limits<double>::epsilon(); //0.00001;

    const geo::Vector_t direction1 = end1 - start1;
    const geo::Vector_t direction2 = end2 - start2;

    const geo::Vector_t w = start1 - start2;

    const double a = direction1.Mag2();
    const double b = direction1.Dot(direction2);
    const double c = direction2.Mag2();
    const double d = direction1.Dot(w);
    const double e = direction2.Dot(w);
    const double D = a * c - b * b;

    double sc, sN, sD = D;
    double tc, tN, tD = D;

    if(D < smallNum)
      {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
      }
    else
      {
        sN = (b * e - c * d)/D;
        tN = (a * e - b * d)/D;

        if(sN < 0.)
          {
            sN = 0.;
            tN = e;
            tD = c;
          }
        else if(sN > sD)
          {
            sN = sD;
            tN = e + b;
            tD = c;
          }
      }

    sc = (std::abs(sN) < smallNum ? 0.0 : sN / sD);
    tc = (std::abs(tN) < smallNum ? 0.0 : tN / tD);

    const geo::Vector_t dP = w + (sc * direction1) - (tc * direction2);

    return dP.R();
  }
}
