/**
 * @file   ChannelMapSBNDAlg.h
 * @brief  Channel mapping for SBND.
 * @date   April 5, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * The channel mapping is derived from the standard one, with the
 * exception of the sorting algorithm which is still custom.
 *
 */

#ifndef SBNDCODE_GEOMETRY_CHANNELMAPSBNDALG_H
#define SBNDCODE_GEOMETRY_CHANNELMAPSBNDALG_H

// SBND libraries
#include "sbndcode/Geometry/GeoObjectSorterSBND.h"

// LArSoft libraries
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"


namespace geo {
  
  /**
   * @brief Custom channel mapping algorithm for SBND.
   *
   * This uses the standard channel mapping, and the custom sorter
   * `GeoObjectSortersbnd`.
   *
   * 
   */
  class ChannelMapSBNDAlg : public ChannelMapStandardAlg {
    
    geo::GeoObjectSorterSBND fSBNDsorter; ///< Sorts geo::XXXGeo objects.
    
      public:
    
    ChannelMapSBNDAlg(fhicl::ParameterSet const& p)
      : ChannelMapStandardAlg(p)
      , fSBNDsorter(p)
      {}
    
    /// Returns a custom SBND sorter.
    virtual geo::GeoObjectSorter const& Sorter() const override 
      { return fSBNDsorter; }
    
    
  }; // class ChannelMapSBNDAlg

} // namespace geo

#endif // SBNDCODE_GEOMETYR_CHANNELMAPSBNDALG_H

