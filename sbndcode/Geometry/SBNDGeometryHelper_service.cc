/**
 * @file   SBNDGeometryHelper_service.cc
 * @brief  Geometry helper service for SBND geometries (implementation).
 * @author Erica Snider (erica@fnal.gov)
 * 
 * Handles SBND-specific information for the generic Geometry service
 * within LArSoft. Derived from the `geo::ExptGeoHelperInterface` class.
 */


#include "sbndcode/Geometry/SBNDGeometryHelper.h"

// SBND code
#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

// framework libraries
#include "canvas/Utilities/Exception.h"


namespace sbnd
{

  SBNDGeometryHelper::SBNDGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & )
    :  fPset( pset ),
       fChannelMap()
  {}


  void  SBNDGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
  {
    fChannelMap.reset();


    fChannelMap = std::make_shared<geo::ChannelMapSBNDAlg>( sortingParameters );
//    fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>( sortingParameters );
    if (!fChannelMap) {
      throw art::Exception(art::errors::NullPointerError)
        << "Failed to create a channel map for SBND geometry!\n";
    }
    geom->ApplyChannelMap(fChannelMap);
  }


  std::shared_ptr<const geo::ChannelMapAlg> SBNDGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

} // namespace sbnd

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper, geo::ExptGeoHelperInterface)
