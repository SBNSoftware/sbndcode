///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeometryHelper_service.cc
///
/// Ported from LArIATAuxDetGeometryHelper_service.cc (Author: rs@fnal.gov)
///
/// \version $Id
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "sbndcode/CRT/CRTGeometryHelper.h"
#include "sbndcode/CRT/CRTChannelMapAlg.h"
#include <memory>

namespace sbnd {

  //------------------------------------------------------------------------
  CRTGeometryHelper::CRTGeometryHelper(
      fhicl::ParameterSet const& pset, art::ActivityRegistry &)
  : fPset(pset), fChannelMap() {}

  //------------------------------------------------------------------------
  void CRTGeometryHelper::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const& sortingParameters,
      geo::AuxDetGeometryCore* geom) {
    fChannelMap = \
      std::make_shared<geo::CRTChannelMapAlg>(sortingParameters);

    if (fChannelMap) {
      geom->ApplyChannelMap(fChannelMap);
    }
  }

  //------------------------------------------------------------------------
  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t
  CRTGeometryHelper::doGetAuxDetChannelMapAlg() const {
    return fChannelMap;
  }

}  // namespace sbnd

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface)

