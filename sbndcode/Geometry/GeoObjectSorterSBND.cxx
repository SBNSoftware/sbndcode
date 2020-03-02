/**
 * @file   GeoObjectSorterSBND.cxx
 * @brief  Algorithm class for sorting standard geo::XXXGeo objects for SBND.
 * @date   April 6, 2017
 * @author petrillo@fnal.gov
 */

#include "sbndcode/Geometry/GeoObjectSorterSBND.h"

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// framework libraries
#include "cmath" // std::abs()


//----------------------------------------------------------------------------
//---  local helper functions
//---
static const double tol = 1e-4; /// Comparison tolerance, in centimeters


/// Returns whether the two values are equal (including some tolerance).
inline bool equal(double a, double b)
  { return std::abs(a - b) <= tol; }


//----------------------------------------------------------------------------
bool CryostatSorter(geo::CryostatGeo const& c1, geo::CryostatGeo const& c2) {
  //
  // sort order for cryostats: by x
  // (not that we have that many in SBND...)
  //
  return (c1.CenterX()) < (c2.CenterX());

} // CryostatSorter()


//----------------------------------------------------------------------------
bool TPCSorter(geo::TPCGeo const& t1, geo::TPCGeo const& t2) {
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
bool WireSorter(geo::WireGeo const& w1, geo::WireGeo const& w2) {

  /*
   * Wire comparison algorithm: compare wire centers:
   *
   * 1. if they have different z, sort by increasing z
   * 2. if they have different y, sort by y:
   *    1. increasing y if thetaZ is larger than 0
   *    2. decreasing y if thetaZ is smaller than 0
   * 3. if they have different x, sort by increasing x
   * 4. otherwise, they are the same wire: return false
   *
   * This is strict weak ordering, as required by std::sort().
   *
   * The definition of "different" is that the difference is
   * larger than the chosen tolerance.
   *
   * The wire plane dimensions are roughly 5 x 4 meters.
   * The angle from z axis is +/- pi/6 for the inclined wires.
   * The same angle for the diagonal of the plane is roughly pi/5.
   * Therefore, the TPC is so "narrow" that there are wires connected to both
   * sides, whose z is the same:
   *        ,------,-------,-----.          ,--------------------.
   *        |  ,.-'    ,.-'    ,.|       uD |.,                  |
   *     vA |-'    ,.-'    ,.-'  |          |  `-.,              |
   *        |  ,.-'    ,.-'    ,.|       uC |.,    '-.,          |     y
   *     vB |-'    ,.-'    ,.-'  |          |  '-.,    '-.,      |    ^
   *        |  ,.-'    ,.-'      |       uB |.,    '-.,    '-.,  |    |
   *     vC |-'    ,.-'          |          |  '-.,    '-.,    '-|    |
   *        |  ,.-'              |       uA |.,    '-.,    '-.,  | ---+----> z
   *     vD |-'                  |          |  '-.,    '-.,    '-|    |
   *        `--------------------'          `------^-------^-----'
   *
   * This cartoon shows two views on the same APA, facing the same drift
   * volume.
   * We ask for wires to be sorted by increasing z coordinate.
   * In the left view (V) the order should be vA, vB, vC, vD, that is left to
   * right. That is also top to bottom.
   * In the right view (U) the order should be uA, uB, uC, uD, that is left to
   * right. In this case that is instead bottom to top.
   * Wires can inform us of their start, end and middle. Start and end might
   * be arbitrarily swapped and are not a reliable starting point. In the
   * order we seek, the middle point ("center") monotonically increases with
   * z, but because the plane is "narrow", some wires can cross it completely
   * (vC and vD). In that case, the z coordinate of the center will be the
   * same and it's not resolving the ambiguity. The y coordinate, on the other
   * end, is always decreasing, and it may be relied upon for all wires.
   * On the other view, the situation is similar, with the ambiguity on wires
   * uC and uD resolved by y coordinate, but in this case the condition is for
   * increasing y coordinate.
   * The algorithm needs to figure out which is the right y direction from
   * just the two wires being compared. The observation is that wires with
   * angle larger than 0 with z are in planes with the leftmost disposition
   * and need a decreasing y coordinate, while the others need an increasing
   * y coordinate.
   *
   * In the case of a "short" TPC, the role of z and y would be swapped.
   *
   * This is stubbornly pretending there is no discontinuity in the wire plane
   * (that is, it is ignoring the junction half way along z.
   */
  decltype(auto) c1 = w1.GetCenter(), c2 = w2.GetCenter();

  //
  // we do z first, which easily resolves the vertical wires:
  //
  if (!equal(c1.Z(), c2.Z())) return c1.Z() < c2.Z();

  //
  // here, z is the same: sort by y
  //
  if (!equal(c1.Y(), c2.Y())) {
    // need to figure out the angle of the wires (we assume both share the same)
    decltype(auto) e1 = w1.GetEnd();

    //
    // We work with end - center (e1 - c1):
    // * if its delta y and delta z have the same sign, thetaZ is positive and
    //   we need a decreasing y
    // * otherwise, we want increasing y
    //
    bool const decreasingY = ((e1.Y() - c1.Y()) > 0) == ((e1.Z() - c1.Z()) > 0);
    if (decreasingY) return c1.Y() > c2.Y(); // decreasing => first upper wires
    else             return c1.Y() < c2.Y(); // increasing => first lower wires
  } // if same y

  //
  // also y is the same... go check x
  //
  if (!equal(c1.X(), c2.X())) {
    // mmh... here we are well beyond SBND realm.
    throw cet::exception("GeoObjectSorterSBND")
      << "Wires differ only for x coordinate... this is not SBND any more!\n";
  //  return c1.X() < c2.X();
  }

  //
  // same center, same wire; strict ordering requires us to return false
  //
  return false;

} // WireSorter()


//----------------------------------------------------------------------------
//---  geo::GeoObjectSorterSBND implementation
//---
geo::GeoObjectSorterSBND::GeoObjectSorterSBND(fhicl::ParameterSet const& p)
  : fDetVersion(p.get<std::string>("DetectorVersion", "SBND"))
  {}


//----------------------------------------------------------------------------
void geo::GeoObjectSorterSBND::SortCryostats
  (std::vector<geo::CryostatGeo>& cgeo) const
  { std::sort(cgeo.begin(), cgeo.end(), CryostatSorter); }

//----------------------------------------------------------------------------
void geo::GeoObjectSorterSBND::SortTPCs(std::vector<geo::TPCGeo>& tgeo) const
  { std::sort(tgeo.begin(), tgeo.end(), TPCSorter); }

//----------------------------------------------------------------------------
void geo::GeoObjectSorterSBND::SortPlanes
  (std::vector<geo::PlaneGeo>& pgeo, geo::DriftDirection_t const driftDir)
  const
{
  // The drift direction has to be set before this method is called.
  // Using the drift direction would render the trick of the sorter unnecessary.

  std::sort(pgeo.begin(), pgeo.end(), PlaneSorter);

  /*
  switch (driftDir) {
    case geo::kPosX:
      std::sort(pgeo.rbegin(), pgeo.rend(), SortPlanes);
      break;
    case geo::kNegX:
      std::sort(pgeo.begin(),  pgeo.end(),  SortPlanes);
      break;
    case geo::kUnknownDrift:
      throw cet::exception("TPCGeo")
        << "Drift direction is unknown, can't sort the planes\n";
  } // driftDir
  */

} // geo::GeoObjectSorterSBND::SortPlanes()

//----------------------------------------------------------------------------
void geo::GeoObjectSorterSBND::SortWires
  (std::vector<geo::WireGeo>& wgeo) const
  { std::sort(wgeo.begin(), wgeo.end(), WireSorter); }



//----------------------------------------------------------------------------

void geo::GeoObjectSorterSBND::SortAuxDets
  (std::vector<geo::AuxDetGeo>& adgeo) const
{
//  std::sort(adgeo.begin(), adgeo.end(), sortAuxDetSBND);
}

void geo::GeoObjectSorterSBND::SortAuxDetSensitive
  (std::vector<geo::AuxDetSensitiveGeo>& adsgeo) const
  {}
