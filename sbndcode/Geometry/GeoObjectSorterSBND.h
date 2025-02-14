/**
 * @file   GeoObjectSortersbnd.h
 * @brief  Algorithm class for sorting standard geo::XXXGeo objects for SBND.
 * @date   April 6, 2017
 * @author petrillo@fnal.gov
 */

#ifndef GEO_GEOOBJECTSORTERSBND_H
#define GEO_GEOOBJECTSORTERSBND_H

// Framework libraries
#include "fhiclcpp/fwd.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorter.h"

namespace geo {

  class GeoObjectSorterSBND: public GeoObjectSorter {
  public:

    explicit GeoObjectSorterSBND(fhicl::ParameterSet const&);

  private:

    bool compareCryostats(CryostatGeo const& c1, CryostatGeo const& c2) const override;
    bool compareOpDets(OpDetGeo const& t1, OpDetGeo const& t2) const override;
    bool compareTPCs(TPCGeo const& t1, TPCGeo const& t2) const override;

  }; // class GeoObjectSorterSBND

} // namespace geo

#endif // GEO_GEOOBJECTSORTERSBND_H
