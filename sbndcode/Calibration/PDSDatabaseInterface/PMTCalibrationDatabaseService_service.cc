/**
 *  Wrapper class for 'PMTCalibrationDatabaseProvider.h'
 *  Andrea Scarpelli (ascarpell@bnl.gov)
 */
// Ported to SBND by Alejandro Sanchez-Castillo, Jan. 2025

#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabaseProvider.h"

// framework libraries
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

// -----------------------------------------------------------------------------
namespace sbndDB {
  class PMTCalibrationDatabaseService;
}
class sbndDB::PMTCalibrationDatabaseService : public IPMTCalibrationDatabaseService,
                                              private PMTCalibrationDatabaseProvider {

  void preBeginRun(const art::Run& run);

  /// Returns a constant pointer to the service provider
  virtual PMTCalibrationDatabaseProvider const* do_provider() const override { return this; }

public:
  PMTCalibrationDatabaseService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

}; // class sbndDB::PMTCalibrationDatabaseService

// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbndDB::PMTCalibrationDatabaseService::PMTCalibrationDatabaseService(
  const fhicl::ParameterSet& pset,
  art::ActivityRegistry& reg)
  : PMTCalibrationDatabaseProvider(pset)
{
  reg.sPreBeginRun.watch(this, &PMTCalibrationDatabaseService::preBeginRun);
}

// -----------------------------------------------------------------------------
void sbndDB::PMTCalibrationDatabaseService::preBeginRun(const art::Run& run)
{
  readPMTCalibrationDatabase(run);
}

// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(sbndDB::PMTCalibrationDatabaseService,
                                   sbndDB::IPMTCalibrationDatabaseService,
                                   SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(sbndDB::PMTCalibrationDatabaseService,
                                  sbndDB::IPMTCalibrationDatabaseService)

// -----------------------------------------------------------------------------