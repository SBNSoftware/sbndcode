///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeometryHelper_service.cc
///
/// Ported from LArIATAuxDetGeometryHelper_service.cc (Author: rs@fnal.gov)
///
/// \version $Id
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "sbndcode/CRT/CRTGeometryHelper.h"
#include "sbndcode/CRT/CRTChannelMapAlg.h"

#include <memory>

namespace sbnd {

  //------------------------------------------------------------------------
  CRTGeometryHelper::CRTGeometryHelper(fhicl::ParameterSet const& pset)
  : fPset(pset) {}

  //------------------------------------------------------------------------
  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t
  CRTGeometryHelper::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const& sortingParameters) const
  {
    return std::make_unique<geo::CRTChannelMapAlg>(sortingParameters);
  }

}  // namespace sbnd

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::CRTGeometryHelper,
                                  geo::AuxDetExptGeoHelperInterface)
