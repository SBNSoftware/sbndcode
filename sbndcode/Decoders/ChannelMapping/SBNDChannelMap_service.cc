/**
*  @file   sbndcode/Decoders/ChannelMapping/SBNDChannelMap_service.cc
 * @brief  Wrapper service for `sbndDB::SBNDChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @author Vu Chi Lan Nguyen (vclnguyen1@sheffield.ac.uk) adapted for SBND
 */

#include "sbndcode/Decoders/ChannelMapping/SBNDChannelMapProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace sbndDB { class SBNDChannelMap; }
class sbndDB::SBNDChannelMap: public SBNDChannelMapProvider {
    
    public:
  
  SBNDChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class sbndDB::SBNDChannelMap


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbndDB::SBNDChannelMap::SBNDChannelMap
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& /* reg */)
  : SBNDChannelMapProvider(pset)
  {}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbndDB::SBNDChannelMap, sbndDB::ISBNDChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(sbndDB::SBNDChannelMap, sbndDB::ISBNDChannelMap)


// -----------------------------------------------------------------------------
