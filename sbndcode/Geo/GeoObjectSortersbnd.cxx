////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSortersbnd.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "sbndcode/Geo/GeoObjectSortersbnd.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"

// #include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {

//----------------------------------------------------------------------------
// Define sort order for cryostats in APA configuration
//   same as standard
static bool sortCryosbnd(const CryostatGeo* c1, const CryostatGeo* c2)
{
  double xyz1[3] = {0.}, xyz2[3] = {0.};
  double local[3] = {0.};
  c1->LocalToWorld(local, xyz1);
  c2->LocalToWorld(local, xyz2);

  return xyz1[0] < xyz2[0];
}


//----------------------------------------------------------------------------
// Define sort order for tpcs in APA configuration.
static bool sortTPCsbnd(const TPCGeo* t1, const TPCGeo* t2)
{
  double xyz1[3] = {0.};
  double xyz2[3] = {0.};
  double local[3] = {0.};
  t1->LocalToWorld(local, xyz1);
  t2->LocalToWorld(local, xyz2);

  // very useful for aligning volume sorting with GDML bounds
  //  --> looking at this output is one way to find the z-borders between APAs,
  //      which tells soreWiresbnd to sort depending on z position via "InVertSplitRegion"
  //mf::LogVerbatim("sortTPCsbnd") << "tpx z range = " << xyz1[2] - t1->Length()/2
  //         << " to " << xyz1[2] + t1->Length()/2;

  // The goal is to number TPCs first in the x direction so that,
  // in the case of APA configuration, TPCs 2c and 2c+1 make up APA c.
  // then numbering will go in y then in z direction.

  // First sort all TPCs into same-z groups
  if (xyz1[2] < xyz2[2]) return true;

  // Within a same-z group, sort TPCs into same-y groups
  if (xyz1[2] == xyz2[2] && xyz1[1] < xyz2[1]) return true;

  // Within a same-z, same-y group, sort TPCs according to x
  if (xyz1[2] == xyz2[2] && xyz1[1] == xyz2[1] && xyz1[0] < xyz2[0]) return true;

  // none of those are true, so return false
  return false;
}


//----------------------------------------------------------------------------
// Define sort order for planes in APA configuration
//   same as standard, but implemented differently
static bool sortPlanesbnd(const PlaneGeo* p1, const PlaneGeo* p2)
{
  double xyz1[3] = {0.};
  double xyz2[3] = {0.};
  double local[3] = {0.};
  p1->LocalToWorld(local, xyz1);
  p2->LocalToWorld(local, xyz2);

  //mf::LogVerbatim("sortPlanessbnd") << "Sorting planes: ("
  //            << xyz1[0] <<","<< xyz1[1] <<","<< xyz1[2] << ") and ("
  //            << xyz2[0] <<","<< xyz2[1] <<","<< xyz2[2] << ")";

  return xyz1[0] > xyz2[0];
}

//----------------------------------------------------------------------------
// we want the wires to be sorted such that the smallest corner wire
// on the readout end of a plane is wire zero, with wire number
// increasing away from that wire.

// Since sbndt has an APA which is both above and below the world origin,
// we cannot use the APA trick. If we could ask where wire 0 was, we could
// still do this in a single implimentation, but we aren't sure what wire
// center we will be getting, so this reversed sorting must be handled
// at the plane level where there is one vertical center.
// If the plane center is above, count from top down (the top stacked and
// largest APAs) If the plane is below (bottom stacked APA) count bottom up
struct sortWiresbnd {

  std::string detVersion;

  sortWiresbnd(std::string detv)
    : detVersion(detv)
  {}

  bool operator()(WireGeo* w1, WireGeo* w2) {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    w1->GetCenter(xyz1); w2->GetCenter(xyz2);


    //mf::LogVerbatim("sortWiresbnd") << "Sorting wires: ("
    //            << xyz1[0] <<","<< xyz1[1] <<","<< xyz1[2] << ") and ("
    //            << xyz2[0] <<","<< xyz2[1] <<","<< xyz2[2] << ")";


    // immedieately take care of vertical wires regardless of which TPC
    // vertical wires should always have same y, and always increase in z direction
    if ( xyz1[1] == xyz2[1] && xyz1[2] < xyz2[2] ) return true;



    ///////////////////////////////////////////////////////////

    // we want the wires to be sorted such that the smallest corner wire
    // on the readout end of a plane is wire zero, with wire number
    // increasing away from that wire.

    if ( xyz1[1] > xyz2[1] ) return true;



    return false;
  }
};

//----------------------------------------------------------------------------
GeoObjectSortersbnd::GeoObjectSortersbnd(fhicl::ParameterSet const& p)
  : fDetVersion(p.get< std::string >("DetectorVersion", "sbnd"))
{
}

//----------------------------------------------------------------------------
GeoObjectSortersbnd::~GeoObjectSortersbnd()
{
}

//----------------------------------------------------------------------------
void GeoObjectSortersbnd::SortCryostats(std::vector<geo::CryostatGeo*> & cgeo) const
{
  std::sort(cgeo.begin(), cgeo.end(), sortCryosbnd);

  return;
}

//----------------------------------------------------------------------------
void GeoObjectSortersbnd::SortTPCs(std::vector<geo::TPCGeo*>  & tgeo) const
{

  std::sort(tgeo.begin(), tgeo.end(), sortTPCsbnd);

  return;
}

//----------------------------------------------------------------------------
void GeoObjectSortersbnd::SortPlanes(std::vector<geo::PlaneGeo*> & pgeo,
                                       geo::DriftDirection_t  const& driftDir) const
{
  // sort the planes to increase in drift direction
  // The drift direction has to be set before this method is called.  It is set when
  // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
  if     (driftDir == geo::kPosX) std::sort(pgeo.rbegin(), pgeo.rend(), sortPlanesbnd);
  else if (driftDir == geo::kNegX) std::sort(pgeo.begin(),  pgeo.end(),  sortPlanesbnd);
  else if (driftDir == geo::kUnknownDrift)
    throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes\n";

  return;
}

//----------------------------------------------------------------------------
void GeoObjectSortersbnd::SortWires(std::vector<geo::WireGeo*> & wgeo) const
{
  sortWiresbnd swsbnd(fDetVersion);

  std::sort(wgeo.begin(), wgeo.end(), swsbnd);

  return;
}



//----------------------------------------------------------------------------

void GeoObjectSortersbnd::SortAuxDets(std::vector<geo::AuxDetGeo*> & adgeo) const
{
//        std::sort(adgeo.begin(), adgeo.end(), sortAuxDetSBND);

  return;
}

void GeoObjectSortersbnd::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const
{
  return;
}


}
