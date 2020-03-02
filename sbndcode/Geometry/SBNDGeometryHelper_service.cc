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
#include "larcorealg/Geometry/ChannelMapAlg.h"

namespace sbnd
{
  SBNDGeometryHelper::SBNDGeometryHelper(fhicl::ParameterSet const&)
  {}

  std::unique_ptr<geo::ChannelMapAlg>
  SBNDGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                                               std::string const& /*detectorName*/) const
  {
    return std::make_unique<geo::ChannelMapSBNDAlg>(sortingParameters);
  }

} // namespace sbnd

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper, geo::ExptGeoHelperInterface)
