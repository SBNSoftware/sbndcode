#ifndef CRTCOMMONUTILS_H_SEEN
#define CRTCOMMONUTILS_H_SEEN


///////////////////////////////////////////////
// CRTCommonUtils.h
//
// Common functions for CRT reconstruction
///////////////////////////////////////////////

// sbndcode includes
#include "sbnobj/Common/CRT/CRTHit.hh"

// c++
#include <vector>
#include <utility>

// ROOT
#include "TVector3.h"

namespace sbnd{
namespace CRTCommonUtils{
  
  // Simple distance of closest approach between infinite track and centre of hit
  double SimpleDCA(sbn::crt::CRTHit hit, TVector3 start, TVector3 direction);

  // Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
  double DistToCrtHit(sbn::crt::CRTHit hit, TVector3 start, TVector3 end);

  // Distance between infinite line (2) and segment (1)
  // http://geomalgorithms.com/a07-_distance.html
  double LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2);

  // Intersection between axis-aligned cube and infinite line
  // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
  std::pair<TVector3, TVector3> CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end);

}
}

#endif
