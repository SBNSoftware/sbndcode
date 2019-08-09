#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

namespace sbnd{

// Simple distance of closest approach between infinite track and centre of hit
double CRTCommonUtils::SimpleDCA(crt::CRTHit hit, TVector3 start, TVector3 direction){

  TVector3 pos (hit.x_pos, hit.y_pos, hit.z_pos);
  TVector3 end = start + direction;
  double denominator = direction.Mag();
  double numerator = (pos - start).Cross(pos - end).Mag();
  return numerator/denominator;

}

// Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
double CRTCommonUtils::DistToCrtHit(crt::CRTHit hit, TVector3 start, TVector3 end){

  // Check if track goes inside hit
  TVector3 min (hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
  TVector3 max (hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
  if(CubeIntersection(min, max, start, end).first.X() != -99999) return 0;

  // Calculate the closest distance to each edge of the CRT hit
  // Assume min error is the fixed position of tagger
  TVector3 vertex1 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
  TVector3 vertex2 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos - hit.z_err);
  TVector3 vertex3 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos + hit.z_err);
  TVector3 vertex4 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
  if(hit.y_err < hit.x_err && hit.y_err < hit.z_err){
    vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
    vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
    vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
    vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
  }
  if(hit.z_err < hit.x_err && hit.z_err < hit.y_err){
    vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
    vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
    vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
    vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
  }

  double dist1 = LineSegmentDistance(vertex1, vertex2, start, end);
  double dist2 = LineSegmentDistance(vertex1, vertex3, start, end);
  double dist3 = LineSegmentDistance(vertex4, vertex2, start, end);
  double dist4 = LineSegmentDistance(vertex4, vertex3, start, end);

  return std::min(std::min(dist1, dist2), std::min(dist3, dist4));

}


// Distance between infinite line (2) and segment (1)
// http://geomalgorithms.com/a07-_distance.html
double CRTCommonUtils::LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2){

  double smallNum = 0.00001;

  // 1 is segment
  TVector3 direction1 = end1 - start1;
  // 2 is infinite line
  TVector3 direction2 = end2 - start2;

  TVector3 u = direction1;
  TVector3 v = direction2;
  TVector3 w = start1 - start2;

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
  TVector3 dP = w + (sc * u) - (tc * v);

  return dP.Mag();

}

// Intersection between axis-aligned cube and infinite line
// (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
std::pair<TVector3, TVector3> CRTCommonUtils::CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end){

  TVector3 dir = (end - start);
  TVector3 invDir (1./dir.X(), 1./dir.Y(), 1/dir.Z());

  double tmin, tmax, tymin, tymax, tzmin, tzmax;

  TVector3 enter (-99999, -99999, -99999);
  TVector3 exit (-99999, -99999, -99999);

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

}

