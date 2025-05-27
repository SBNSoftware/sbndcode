#include "sbndcode/Geometry/WireReadoutSorterSBND.h"

// LArSoft libraries
#include "larcorealg/Geometry/WireGeo.h"

// framework libraries
#include "cetlib_except/exception.h"

// C++ standard libraries
#include <cmath> // std::abs()

namespace {
  constexpr double tolerance = 1e-4; /// Comparison tolerance, in centimeters

  inline bool equal(double a, double b)
    { return std::abs(a - b) <= tolerance; }
}

//----------------------------------------------------------------------------
//---  geo::WireReadoutSorterSBND implementation
//---
geo::WireReadoutSorterSBND::WireReadoutSorterSBND() = default;
geo::WireReadoutSorterSBND::WireReadoutSorterSBND(fhicl::ParameterSet const&)
  {}

//----------------------------------------------------------------------------
bool geo::WireReadoutSorterSBND::compareWires(WireGeo const& w1,
                                              WireGeo const& w2) const
{
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
    throw cet::exception("WireReadoutSorterSBND")
      << "Wires differ only for x coordinate... this is not SBND any more!\n";
  //  return c1.X() < c2.X();
  }

  //
  // same center, same wire; strict ordering requires us to return false
  //
  return false;

} // compareWires
