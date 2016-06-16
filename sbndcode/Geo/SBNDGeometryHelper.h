////////////////////////////////////////////////////////////////////////////////
/// \file SBNDGeometryHelper.h
/// \brief Geometry helper service for SBND geometries. 
/// 
/// Handles SBND-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef SBND_ExptGeoHelperInterface_h
#define SBND_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"

#include <memory>


namespace geo{
  class ChannelMapAlg;
  class GeometryCore;
}

namespace sbnd
{
  class SBNDGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    SBNDGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );

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
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // SBND_ExptGeoHelperInterface_h
