///////////////////////////////////////////////////////////////////////////////
/// \file SBNDAuxDetGeometryHelper.h
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

#ifndef SBND_AuxDetExptGeoHelperInterface_h
#define SBND_AuxDetExptGeoHelperInterface_h

#include <memory>
#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"

namespace sbnd {

  class SBNDAuxDetGeometryHelper : public geo::AuxDetExptGeoHelperInterface {
  public:

    SBNDAuxDetGeometryHelper(fhicl::ParameterSet const & pset,
			     art::ActivityRegistry &);

  private:

    virtual void doConfigureAuxDetChannelMapAlg(
        fhicl::ParameterSet const& sortingParameters,
        geo::AuxDetGeometryCore* geom) override;

    virtual AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;

    fhicl::ParameterSet fPset; ///< Copy of configuration parameter set
    std::shared_ptr<geo::AuxDetChannelMapAlg> fChannelMap; ///< Channel map

  };

}  // namespace sbnd

DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDAuxDetGeometryHelper, geo::AuxDetExptGeoHelperInterface, LEGACY)

#endif  // SBND_AuxDetExptGeoHelperInterface_h

