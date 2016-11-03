///////////////////////////////////////////////////////////////////////////////
/// \file SBNDAuxDetGeometryHelper_service.cc
///
/// Ported from LArIATAuxDetGeometryHelper_service.cc (Author: rs@fnal.gov)
///
/// \version $Id
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "larcore/Geometry/AuxDetGeometryCore.h"
#include "sbndcode/Geo/SBNDAuxDetGeometryHelper.h"
#include "sbndcode/Geo/AuxDetChannelMapSBNDAlg.h"
#include <memory>

namespace sbnd {

  //------------------------------------------------------------------------
  SBNDAuxDetGeometryHelper::SBNDAuxDetGeometryHelper(
      fhicl::ParameterSet const& pset, art::ActivityRegistry &)
  : fPset(pset), fChannelMap() {}

  //------------------------------------------------------------------------
  void SBNDAuxDetGeometryHelper::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const& sortingParameters,
      geo::AuxDetGeometryCore* geom) {
    fChannelMap = \
      std::make_shared<geo::AuxDetChannelMapSBNDAlg>(sortingParameters);

    if (fChannelMap) {
      geom->ApplyChannelMap(fChannelMap);
    }
  }

  //------------------------------------------------------------------------
  SBNDAuxDetGeometryHelper::AuxDetChannelMapAlgPtr_t
  SBNDAuxDetGeometryHelper::doGetAuxDetChannelMapAlg() const {
    return fChannelMap;
  }

}  // namespace sbnd

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDAuxDetGeometryHelper, geo::AuxDetExptGeoHelperInterface)

