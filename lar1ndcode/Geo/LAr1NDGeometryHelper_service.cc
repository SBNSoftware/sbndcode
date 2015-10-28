////////////////////////////////////////////////////////////////////////////////
/// \file LAr1NDGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> lar1nd/Geometry
#include "Geometry/GeometryCore.h"
#include "lar1ndcode/Geo/LAr1NDGeometryHelper.h"

#include "Geometry/ChannelMapAlg.h"

// Migration note:
// Geometry --> lar1nd/Geometry for the two below
#include "lar1ndcode/Geo/ChannelMaplar1ndAlg.h"


#include "TString.h"


namespace lar1nd
{

  LAR1NDGeometryHelper::LAR1NDGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & )
    :  fPset( pset ),
       fChannelMap()
  {}


  void  LAR1NDGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
  {
    fChannelMap.reset();


    fChannelMap = std::make_shared<geo::ChannelMaplar1ndAlg>( sortingParameters );
    if ( fChannelMap )
    {
      geom->ApplyChannelMap(fChannelMap);
    }
  }


  std::shared_ptr<const geo::ChannelMapAlg> LAR1NDGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lar1nd::LAR1NDGeometryHelper, geo::ExptGeoHelperInterface)
