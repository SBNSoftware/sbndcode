////////////////////////////////////////////////////////////////////////////////
/// \file LArIATGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "lar1ndcode/Geo/LArIATGeometryHelper.h"
#include "Geometry/ChannelMapAlg.h"

#include "Geometry/ChannelMapStandardAlg.h"

#include "TString.h"


namespace lariatgeo
{

  //------------------------------------------------------------------------
  LArIATGeometryHelper::LArIATGeometryHelper(fhicl::ParameterSet const& pset, 
					     art::ActivityRegistry    & reg )
  : fPset(pset)
  , fReg(reg)
  , fChannelMap()
  {}

  //------------------------------------------------------------------------
  LArIATGeometryHelper::~LArIATGeometryHelper() throw()
  {}  
  
  //------------------------------------------------------------------------
  void LArIATGeometryHelper::doConfigureChannelMapAlg(const TString                  & detectorName,
						      fhicl::ParameterSet      const & sortingParam,
						      std::vector<geo::CryostatGeo*> & c )
  {
    fChannelMap = nullptr;
    
    fChannelMap = std::shared_ptr<geo::ChannelMapAlg>(new geo::ChannelMapStandardAlg(sortingParam));
    if(fChannelMap) fChannelMap->Initialize(c);

    return;
  }
  
  //------------------------------------------------------------------------
  std::shared_ptr<const geo::ChannelMapAlg> LArIATGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATGeometryHelper, geo::ExptGeoHelperInterface)
