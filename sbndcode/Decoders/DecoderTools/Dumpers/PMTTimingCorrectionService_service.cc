/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrectionService_service.cc
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 */

#include "sbndcode/Decoders/DecoderTools/Dumpers/IPMTTimingCorrectionService.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrectionsProvider.h"

// framework libraries
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace sbndDB { class PMTTimingCorrectionService; }
class sbndDB::PMTTimingCorrectionService
  : public IPMTTimingCorrectionService, private PMTTimingCorrectionsProvider {
    
      void preBeginRun(const art::Run& run);
    
      /// Returns a constant pointer to the service provider
      virtual PMTTimingCorrectionsProvider const* do_provider() const override
         { return this; }
      
    public:
  
      PMTTimingCorrectionService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class sbndDB::PMTTimingCorrectionService


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbndDB::PMTTimingCorrectionService::PMTTimingCorrectionService
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : PMTTimingCorrectionsProvider(pset)
{
  reg.sPreBeginRun.watch(this, &PMTTimingCorrectionService::preBeginRun);
}


// -----------------------------------------------------------------------------
void sbndDB::PMTTimingCorrectionService::preBeginRun(const art::Run& run)
{
  readTimeCorrectionDatabase(run);
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbndDB::PMTTimingCorrectionService, sbndDB::IPMTTimingCorrectionService, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(sbndDB::PMTTimingCorrectionService, sbndDB::IPMTTimingCorrectionService)


// -----------------------------------------------------------------------------
