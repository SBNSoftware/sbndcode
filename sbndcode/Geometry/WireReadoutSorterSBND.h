/**
 * @file   GeoObjectSortersbnd.h
 * @brief  Algorithm class for sorting standard geo::XXXGeo objects for SBND.
 * @date   April 6, 2017
 * @author petrillo@fnal.gov
 */

#ifndef GEO_WIREREADOUTSORTERSBND_H
#define GEO_WIREREADOUTSORTERSBND_H

// Framework libraries
#include "fhiclcpp/fwd.h"

// LArSoft libraries
#include "larcorealg/Geometry/WireReadoutSorter.h"
#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class WireReadoutSorterSBND: public WireReadoutSorter {
  public:
    WireReadoutSorterSBND();
    explicit WireReadoutSorterSBND(fhicl::ParameterSet const&);

  private:
    bool compareWires(WireGeo const& w1, WireGeo const& w2) const override;
  }; // class WireReadoutSorterSBND

} // namespace geo

#endif // GEO_WIREREADOUTSORTERSBND_H
