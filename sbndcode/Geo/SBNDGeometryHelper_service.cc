////////////////////////////////////////////////////////////////////////////////
/// \file SBNDGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> sbnd/Geometry
#include "larcore/Geometry/GeometryCore.h"
#include "sbndcode/Geo/SBNDGeometryHelper.h"

#include "larcore/Geometry/ChannelMapAlg.h"

// Migration note:
// #include "sbndcode/Geo/ChannelMapsbndAlg.h"
#include "larcore/Geometry/ChannelMapStandardAlg.h"


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


//    fChannelMap = std::make_shared<geo::ChannelMapsbndAlg>( sortingParameters );
    fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>( sortingParameters );
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

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper, geo::ExptGeoHelperInterface)
