///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeometryHelper.h
/// \brief Auxiliary detector geometry helper service for SBND geometries.
///
/// Handles SBND-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class.
///
/// Ported from LArIATAuxDetGeometryHelper.h (Author: brebel@fnal.gov)
///
/// \verion $Id
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef SBND_CRTExptGeoHelperInterface_h
#define SBND_CRTExptGeoHelperInterface_h

#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"
#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"
#include "sbndcode/CRT/CRTChannelMapAlg.h"

namespace sbnd {

  class CRTGeometryHelper : public geo::AuxDetExptGeoHelperInterface {
  public:

    explicit CRTGeometryHelper(fhicl::ParameterSet const & pset);

  private:

    AuxDetChannelMapAlgPtr_t
    doConfigureAuxDetChannelMapAlg(fhicl::ParameterSet const& sortingParameters) const override;

    fhicl::ParameterSet fPset; ///< Copy of configuration parameter set
  };

}  // namespace sbnd

DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::CRTGeometryHelper,
                                   geo::AuxDetExptGeoHelperInterface,
                                   SHARED)

#endif  // SBND_CRTExptGeoHelperInterface_h
