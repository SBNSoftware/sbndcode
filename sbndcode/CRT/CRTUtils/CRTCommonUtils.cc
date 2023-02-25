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

  std::pair<geo::Point_t, geo::Point_t> CRTCommonUtils::CubeIntersection(const geo::Point_t &min, const geo::Point_t &max,
                                                                         const geo::Point_t &start, const geo::Point_t &end)
  {
    geo::Vector_t dir = (end - start);
    geo::Vector_t invDir(1./dir.X(), 1./dir.Y(), 1./dir.Z());

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    geo::Point_t enter (-99999, -99999, -99999);
    geo::Point_t exit (-99999, -99999, -99999);

    // Find the intersections with the X plane
    if(invDir.X() >= 0){
      tmin = (min.X() - start.X()) * invDir.X();
      tmax = (max.X() - start.X()) * invDir.X();
    }
    else{
      tmin = (max.X() - start.X()) * invDir.X();
      tmax = (min.X() - start.X()) * invDir.X();
    }

    // Find the intersections with the Y plane
    if(invDir.Y() >= 0){
      tymin = (min.Y() - start.Y()) * invDir.Y();
      tymax = (max.Y() - start.Y()) * invDir.Y();
    }
    else{
      tymin = (max.Y() - start.Y()) * invDir.Y();
      tymax = (min.Y() - start.Y()) * invDir.Y();
    }

    // Check that it actually intersects
    if((tmin > tymax) || (tymin > tmax)) return std::make_pair(enter, exit);

    // Max of the min points is the actual intersection
    if(tymin > tmin) tmin = tymin;

    // Min of the max points is the actual intersection
    if(tymax < tmax) tmax = tymax;

    // Find the intersection with the Z plane
    if(invDir.Z() >= 0){
      tzmin = (min.Z() - start.Z()) * invDir.Z();
      tzmax = (max.Z() - start.Z()) * invDir.Z();
    }
    else{
      tzmin = (max.Z() - start.Z()) * invDir.Z();
      tzmax = (min.Z() - start.Z()) * invDir.Z();
    }

    // Check for intersection
    if((tmin > tzmax) || (tzmin > tmax)) return std::make_pair(enter, exit);

    // Find final intersection points
    if(tzmin > tmin) tmin = tzmin;

    // Find final intersection points
    if(tzmax < tmax) tmax = tzmax;

    // Calculate the actual crossing points
    double xmin = start.X() + tmin * dir.X();
    double xmax = start.X() + tmax * dir.X();
    double ymin = start.Y() + tmin * dir.Y();
    double ymax = start.Y() + tmax * dir.Y();
    double zmin = start.Z() + tmin * dir.Z();
    double zmax = start.Z() + tmax * dir.Z();

    // Return pair of entry and exit points
    enter.SetXYZ(xmin, ymin, zmin);
    exit.SetXYZ(xmax, ymax, zmax);
    return std::make_pair(enter, exit);
  }

  double CRTCommonUtils::SimpleDCA(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Vector_t &direction)
  {
    geo::Point_t pos = sp->Pos();
    geo::Point_t end = start + direction;
    double denominator = direction.R();
    double numerator = (pos - start).Cross(pos - end).R();
    return numerator/denominator;
  }

  double CRTCommonUtils::DistToCRTSpacePoint(const art::Ptr<CRTSpacePoint> &sp, const geo::Point_t &start, const geo::Point_t &end)
  {
    // Check if track goes inside hit
    geo::Point_t min = sp->Pos() - geo::Vector_t(sp->Err());
    geo::Point_t max = sp->Pos() + geo::Vector_t(sp->Err());

    if(CubeIntersection(min, max, start, end).first.X() != -99999) return 0;

    // Calculate the closest distance to each edge of the CRT hit
    // Assume min error is the fixed position of tagger
    geo::Point_t vertex1 (sp->X(), sp->Y() - sp->YErr(), sp->Z() - sp->ZErr());
    geo::Point_t vertex2 (sp->X(), sp->Y() + sp->YErr(), sp->Z() - sp->ZErr());
    geo::Point_t vertex3 (sp->X(), sp->Y() - sp->YErr(), sp->Z() + sp->ZErr());
    geo::Point_t vertex4 (sp->X(), sp->Y() + sp->YErr(), sp->Z() + sp->ZErr());
    if(sp->YErr() < sp->XErr() && sp->YErr() < sp->ZErr()){
      vertex1.SetXYZ(sp->X() - sp->XErr(), sp->Y(), sp->Z() - sp->ZErr());
      vertex2.SetXYZ(sp->X() + sp->XErr(), sp->Y(), sp->Z() - sp->ZErr());
      vertex3.SetXYZ(sp->X() - sp->XErr(), sp->Y(), sp->Z() + sp->ZErr());
      vertex4.SetXYZ(sp->X() + sp->XErr(), sp->Y(), sp->Z() + sp->ZErr());
    }
    if(sp->ZErr() < sp->XErr() && sp->ZErr() < sp->YErr()){
      vertex1.SetXYZ(sp->X() - sp->XErr(), sp->Y() - sp->YErr(), sp->Z());
      vertex2.SetXYZ(sp->X() + sp->XErr(), sp->Y() - sp->YErr(), sp->Z());
      vertex3.SetXYZ(sp->X() - sp->XErr(), sp->Y() + sp->YErr(), sp->Z());
      vertex4.SetXYZ(sp->X() + sp->XErr(), sp->Y() + sp->YErr(), sp->Z());
    }

    double dist1 = LineSegmentDistance(vertex1, vertex2, start, end);
    double dist2 = LineSegmentDistance(vertex1, vertex3, start, end);
    double dist3 = LineSegmentDistance(vertex4, vertex2, start, end);
    double dist4 = LineSegmentDistance(vertex4, vertex3, start, end);

    return std::min(std::min(dist1, dist2), std::min(dist3, dist4));
  }

  double CRTCommonUtils::LineSegmentDistance(const geo::Point_t &start1, const geo::Point_t &end1, const geo::Point_t &start2, const geo::Point_t &end2)
  {
    double smallNum = 0.00001;

    // 1 is segment
    geo::Vector_t direction1 = end1 - start1;
    // 2 is infinite line
    geo::Vector_t direction2 = end2 - start2;

    geo::Vector_t u = direction1;
    geo::Vector_t v = direction2;
    geo::Vector_t w = start1 - start2;

    double a = u.Dot(u);
    double b = u.Dot(v);
    double c = v.Dot(v);
    double d = u.Dot(w);
    double e = v.Dot(w);
    double D = a * c - b * b;
    double sc, sN, sD = D; // sc = sN/sD
    double tc, tN, tD = D; // sc = sN/sD

    // Compute the line parameters of the two closest points
    if(D < smallNum){ // Lines are almost parallel
      sN = 0.0;
      sD = 1.0;
      tN = e;
      tD = c;
    }
    else{
      sN = (b * e - c * d)/D;
      tN = (a * e - b * d)/D;
      if(sN < 0.){ // sc < 0, the s = 0 edge is visible
        sN = 0.;
        tN = e;
        tD = c;
      }
      else if(sN > sD){ // sc > 1, the s = 1 edge is visible
        sN = sD;
        tN = e + b;
        tD = c;
      }
    }

    sc = (std::abs(sN) < smallNum ? 0.0 : sN / sD);
    tc = (std::abs(tN) < smallNum ? 0.0 : tN / tD);
    // Get the difference of the two closest points
    geo::Vector_t dP = w + (sc * u) - (tc * v);

    return dP.R();
  }
}
