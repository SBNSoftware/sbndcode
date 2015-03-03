////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterlar1nd.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov  // updated to lar1nd by A. Szelc
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERLAR1ND_H
#define GEO_GEOOBJECTSORTERLAR1ND_H

#include <vector>
#include <string>
#include "Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorterlar1nd : public GeoObjectSorter {

  public:

    GeoObjectSorterlar1nd(fhicl::ParameterSet const& p);
    ~GeoObjectSorterlar1nd();

    void SortCryostats(std::vector<geo::CryostatGeo*> & cgeo)     const;
    void SortTPCs     (std::vector<geo::TPCGeo*>      & tgeo)     const;
    void SortPlanes   (std::vector<geo::PlaneGeo*>    & pgeo,
		       geo::DriftDirection_t     const& driftDir) const;
    void SortWires    (std::vector<geo::WireGeo*>     & wgeo)     const;
    void SortAuxDets  (std::vector<geo::AuxDetGeo*>   & adgeo)    const;    
  private:
    
    std::string fDetVersion;  ///< string of the detector version
  };

}

#endif // GEO_GEOOBJECTSORTERlar1nd_H
