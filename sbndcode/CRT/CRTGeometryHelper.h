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
#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "sbndcode/CRT/CRTChannelMapAlg.h"
#include <memory>

namespace sbnd {

  class CRTGeometryHelper : public geo::AuxDetExptGeoHelperInterface {
  public:

    CRTGeometryHelper(fhicl::ParameterSet const & pset,
                      art::ActivityRegistry &);

  private:

    virtual void doConfigureAuxDetChannelMapAlg(
        fhicl::ParameterSet const& sortingParameters,
        geo::AuxDetGeometryCore* geom) override;

    virtual AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;

    fhicl::ParameterSet fPset; ///< Copy of configuration parameter set
    std::shared_ptr<geo::CRTChannelMapAlg> fChannelMap; ///< Channel map

  };

}  // namespace sbnd

DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface, LEGACY)

#endif  // SBND_CRTExptGeoHelperInterface_h

