////////////////////////////////////////////////////////////////////////////////
/// \file LArIATGeometryHelper.h
/// \brief Geometry helper service for LArIAT geometries. 
/// 
/// Handles LArIAT-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef LArIAT_ExptGeoHelperInterface_h
#define LArIAT_ExptGeoHelperInterface_h

#include "Geometry/ExptGeoHelperInterface.h"

#include <memory>
#include <vector>

// Forward declarations
//
class TString;

namespace geo
{
  class ChannelMapAlg;
  class CryostaGeo;
  class ExptGeoHelperInterface;
}

// Declaration
//
namespace lariatgeo
{
  class LArIATGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    LArIATGeometryHelper(fhicl::ParameterSet const & pset, 
			 art::ActivityRegistry &reg );
    ~LArIATGeometryHelper() throw();

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
    
    void  doConfigureChannelMapAlg(const TString & detectorName,
				   fhicl::ParameterSet const & sortingParam,
				   std::vector<geo::CryostatGeo*> & c ) override;
    std::shared_ptr<const geo::ChannelMapAlg> doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet   const&        fPset;       ///< configuration parameter set
    art::ActivityRegistry      &        fReg;        ///< activity registry
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap; ///< channel map
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // LArIAT_ExptGeoHelperInterface_h
