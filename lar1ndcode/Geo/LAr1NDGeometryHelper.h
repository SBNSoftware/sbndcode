////////////////////////////////////////////////////////////////////////////////
/// \file LAR1NDGeometryHelper.h
/// \brief Geometry helper service for LAR1ND geometries. 
/// 
/// Handles LAR1ND-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef LAR1ND_ExptGeoHelperInterface_h
#define LAR1ND_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"

#include <memory>


namespace geo{
  class ChannelMapAlg;
  class GeometryCore;
}

namespace lar1nd
{
  class LAR1NDGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    LAR1NDGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );

    // Public interface for ExptGeoHelperInterface (for reference purposes)
    //
    // Configure and initialize the channel map.
    //
    // void  ConfigureChannelMapAlg( const TString & detectorName, 
    //                               fhicl::ParameterSet const & sortingParam,
    //                               std::vector<geo::CryostatGeo*> & c );
    //
    // Returns null pointer if the initialization failed
    // NOTE:  the sub-class owns the ChannelMapAlg object
    //
    // std::shared_ptr<const geo::ChannelMapAlg> & GetChannelMapAlg() const;
  
  private:
    

    void  doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom) override;
    virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
        
    fhicl::ParameterSet const & fPset;
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lar1nd::LAR1NDGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // LAR1ND_ExptGeoHelperInterface_h
