/**
 * @file   GeoObjectSortersbnd.h
 * @brief  Algorithm class for sorting standard geo::XXXGeo objects for SBND.
 * @date   April 6, 2017
 * @author petrillo@fnal.gov
 */

#ifndef GEO_GEOOBJECTSORTERSBND_H
#define GEO_GEOOBJECTSORTERSBND_H


// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorter.h"

// C/C++ standard libraries
#include <vector>
#include <string>


namespace geo {

  class GeoObjectSorterSBND: public GeoObjectSorter {
    
      public:
     
    GeoObjectSorterSBND(fhicl::ParameterSet const& p);
  
    virtual void SortCryostats(std::vector<geo::CryostatGeo*> & cgeo)     const override;
    virtual void SortTPCs     (std::vector<geo::TPCGeo*>      & tgeo)     const override;
    virtual void SortPlanes   (std::vector<geo::PlaneGeo*>    & pgeo,
                               geo::DriftDirection_t     const& driftDir) const override;
    virtual void SortWires    (std::vector<geo::WireGeo*>     & wgeo)     const override;
    virtual void SortAuxDets  (std::vector<geo::AuxDetGeo*>   & adgeo)    const override;
    virtual void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const override;
    
  private:
    
    std::string fDetVersion;  ///< String of the detector version.
    
  }; // class GeoObjectSorterSBND

} // namespace geo

#endif // GEO_GEOOBJECTSORTERSBND_H
