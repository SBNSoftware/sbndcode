/**
 * @file   GeoObjectSorterSBND.cxx
 * @brief  Algorithm class for sorting standard geo::XXXGeo objects for SBND.
 * @date   April 6, 2017
 * @author petrillo@fnal.gov
 */

#include "sbndcode/Geometry/GeoObjectSorterSBND.h"

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// framework libraries
#include <cmath> // std::abs()

namespace {
  constexpr double tolerance = 1e-4; /// Comparison tolerance, in centimeters

/// Returns whether the two values are equal (including some tolerance).
inline bool equal(double a, double b)
  { return std::abs(a - b) <= tolerance; }
}

//----------------------------------------------------------------------------
// Define sort order for planes in APA configuration
//   same as standard, but implemented differently
bool PlaneSorter(geo::PlaneGeo const& p1, geo::PlaneGeo const& p2) {

  /*
   * Sort the wire planes so that the first faces the center of the TPC.
   *
   * This might be moved to geo::TPCGeo, which has all the information to
   * enforce such a policy.
   * In fact, we don't have as many here.
   *
   * We rely on a trick, that assumes the cathode to be located at x = 0,
   * and we sort after |x| of the wire planes.
   *
   */
  decltype(auto) p1c = p1.GetBoxCenter();
  decltype(auto) p2c = p2.GetBoxCenter();

  /*
   *   #2  #1  #0                   cathode                   #0  #1  #2
   *   |   |   |                       |                       |   |   |
   *   |   |   |                       |                       |   |   |
   *   |   |   |                       |                       |   |   |
   *   |   |   |                       |                       |   |   |
   *
   */
  return std::abs(p1c.X()) < std::abs(p2c.X());

} // PlaneSorter()


//----------------------------------------------------------------------------
//---  geo::GeoObjectSorterSBND implementation
//---
geo::GeoObjectSorterSBND::GeoObjectSorterSBND(fhicl::ParameterSet const&)
  {}

//----------------------------------------------------------------------------
bool geo::GeoObjectSorterSBND::compareCryostats(CryostatGeo const& c1, CryostatGeo const& c2) const
  {
  //
  // sort order for cryostats: by x
  // (not that we have that many in SBND...)
  //
  return (c1.CenterX()) < (c2.CenterX());
} // CryostatSorter()


//----------------------------------------------------------------------------
bool geo::GeoObjectSorterSBND::compareOpDets(OpDetGeo const& t1, OpDetGeo const& t2) const
{
  auto const xyz1 = t1.GetCenter();
  auto const xyz2 = t2.GetCenter();

  if(xyz1.Z() != xyz2.Z())
    return xyz1.Z() < xyz2.Z();
  if(xyz1.Y() != xyz2.Y())
    return xyz1.Y() < xyz2.Y();
  return xyz1.X() < xyz2.X();
} // OpDetsSorter

//----------------------------------------------------------------------------
bool geo::GeoObjectSorterSBND::compareTPCs(TPCGeo const& t1, TPCGeo const& t2) const
{
  //
  // Define sort order for TPCs (in SBND case, by x).
  //

  // The goal is to number TPCs first in the x direction so that,
  // in the case of APA configuration, TPCs 2c and 2c+1 make up APA c.
  // then numbering will go in y then in z direction.

  // First sort all TPCs belonging to different "z groups"
  if (!equal(t1.CenterZ(), t2.CenterZ()))
    return t1.CenterZ() < t2.CenterZ();

  // Within the same-z groups, sort TPCs belonging to different "y groups"
  if (!equal(t1.CenterY(), t2.CenterY()))
    return t1.CenterY() < t2.CenterY();

  // Within the same z and y groups, sort TPCs belonging to different
  // "x groups";
  // if the x is also the same, then t1 and t2 are the same TPC and strict
  // ordering requires us to return false.
  return t1.CenterX() < t2.CenterX();

} // TPCSorter()
