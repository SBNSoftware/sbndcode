/**
 * @file   SBNDGeometryHelper.h
 * @brief  Geometry helper service for SBND geometries.
 * @author Erica Snider (erica@fnal.gov)
 * 
 * Handles SBND-specific information for the generic Geometry service
 * within LArSoft. Derived from the `geo::ExptGeoHelperInterface` class.
 */

#ifndef SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H
#define SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H


// LArSoft libraries
#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

// C++ libraries
#include <memory> // std::shared_ptr<>



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
        
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbnd::SBNDGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // SBNDCODE_GEOMETY_SBNDGEOMETRYHELPER_H
