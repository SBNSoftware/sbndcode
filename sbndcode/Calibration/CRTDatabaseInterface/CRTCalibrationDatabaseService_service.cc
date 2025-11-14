#include "sbndcode/Calibration/CRTDatabaseInterface/ICRTCalibrationDatabaseService.h"
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabaseProvider.h"

// framework libraries
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

// -----------------------------------------------------------------------------
namespace sbndDB
{
  class CRTCalibrationDatabaseService;
}

class sbndDB::CRTCalibrationDatabaseService : public ICRTCalibrationDatabaseService,
                                              private CRTCalibrationDatabaseProvider
{
  void preBeginRun(const art::Run& run);

  virtual CRTCalibrationDatabaseProvider const* do_provider() const override { return this; }

public:
  CRTCalibrationDatabaseService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

};

sbndDB::CRTCalibrationDatabaseService::CRTCalibrationDatabaseService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : CRTCalibrationDatabaseProvider(pset)
{
  reg.sPreBeginRun.watch(this, &CRTCalibrationDatabaseService::preBeginRun);
}

void sbndDB::CRTCalibrationDatabaseService::preBeginRun(const art::Run& run)
{
  readCRTCalibrationDatabase(run);
}

DECLARE_ART_SERVICE_INTERFACE_IMPL(sbndDB::CRTCalibrationDatabaseService,
                                   sbndDB::ICRTCalibrationDatabaseService,
                                   SHARED)

DEFINE_ART_SERVICE_INTERFACE_IMPL(sbndDB::CRTCalibrationDatabaseService,
                                  sbndDB::ICRTCalibrationDatabaseService)
