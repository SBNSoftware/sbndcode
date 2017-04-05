////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSortersbnd.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov  // updated to sbnd by A. Szelc
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERSBND_H
#define GEO_GEOOBJECTSORTERSBND_H

#include <vector>
#include <string>
#include "larcore/Geometry/GeoObjectSorter.h"



namespace geo {

class GeoObjectSortersbnd : public GeoObjectSorter {

public:

  GeoObjectSortersbnd(fhicl::ParameterSet const& p);
  ~GeoObjectSortersbnd();

  void SortCryostats(std::vector<geo::CryostatGeo*> & cgeo)     const;
  void SortTPCs     (std::vector<geo::TPCGeo*>      & tgeo)     const;
  void SortPlanes   (std::vector<geo::PlaneGeo*>    & pgeo,
                     geo::DriftDirection_t     const& driftDir) const;
  void SortWires    (std::vector<geo::WireGeo*>     & wgeo)     const;
  void SortAuxDets  (std::vector<geo::AuxDetGeo*>   & adgeo)    const;

  void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const;
  
private:

  std::string fDetVersion;  ///< string of the detector version
};

}

#endif // GEO_GEOOBJECTSORTERsbnd_H
