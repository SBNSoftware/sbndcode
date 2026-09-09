#ifndef SBNDCODE_DATABASEINTERFACE_ICRTCALIBRATIONDATABASESERVICE_H
#define SBNDCODE_DATABASEINTERFACE_ICRTCALIBRATIONDATABASESERVICE_H

#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabase.h"
#include "larcore/CoreUtils/ServiceProviderWrappers.h"

namespace sbndDB
{
  using ICRTCalibrationDatabaseService = lar::ServiceProviderInterfaceWrapper<CRTCalibrationDatabase>;
}

DECLARE_ART_SERVICE_INTERFACE(sbndDB::ICRTCalibrationDatabaseService, SHARED)

#endif
