////////////////////////////////////////////////////////////////////////////////
/// \file LAr1NDGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> lar1nd/Geometry
#include "lar1ndcode/Geo/LAr1NDGeometryHelper.h"

#include "Geometry/ChannelMapAlg.h"

// Migration note:
// Geometry --> lar1nd/Geometry for the two below
#include "lar1ndcode/Geo/ChannelMaplar1ndAlg.h"


#include "TString.h"


namespace lar1nd
{

  LAR1NDGeometryHelper::LAR1NDGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg )
  :  fPset( pset ),
     fReg( reg ),
     fChannelMap()
  {}

  LAR1NDGeometryHelper::~LAR1NDGeometryHelper() throw()
  {}  
  
  void LAR1NDGeometryHelper::doConfigureChannelMapAlg( const TString & detectorName,
                                                     fhicl::ParameterSet const & sortingParam,
                                                     std::vector<geo::CryostatGeo*> & c )
  {
    fChannelMap = nullptr;
    
  
      fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMaplar1ndAlg( sortingParam ) );
  
    if ( fChannelMap )
    {
      fChannelMap->Initialize( c );
    }
  }
  
  
  std::shared_ptr<const geo::ChannelMapAlg> LAR1NDGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lar1nd::LAR1NDGeometryHelper, geo::ExptGeoHelperInterface)
