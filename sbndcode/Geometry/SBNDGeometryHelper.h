/**
 * @file   SBNDGeometryHelper.h
 * @brief  Geometry helper service for SBND geometries.
 * @author Erica Snider (erica@fnal.gov)
 *
 * Handles SBND-specific information for the generic Geometry service
 * within LArSoft. Derived from the `geo::ExptGeoHelperInterface` class.
 */

#ifndef SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H
#define SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H


// LArSoft libraries
#include "larcore/Geometry/ExptGeoHelperInterface.h"

// C++ libraries
#include <memory> // std::unique_ptr<>

namespace geo{
  class ChannelMapAlg;
}

namespace sbnd
{
  class SBNDGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
    SBNDGeometryHelper(fhicl::ParameterSet const& pset);

  private:
    std::unique_ptr<geo::ChannelMapAlg>
    doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                             std::string const& detectorName) const override;
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper,
                                   geo::ExptGeoHelperInterface,
                                   SHARED)

#endif // SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H
