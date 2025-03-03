/**
 * @file   icaruscode/Timing/IPMTTimingCorrectionService.h
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Andrea Scarpelli (ascarpell@bnl.gov), Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef SBNDCCODE_DATABASEINTERFACE_IPMTCALIBRATIONDATABASESERVICE_H
#define SBNDCCODE_DATABASEINTERFACE_IPMTCALIBRATIONDATABASESERVICE_H

// ICARUS libraries
#include "sbndcode/DatabaseInterface/PMTCalibrationDatabase.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceProviderWrappers.h"


// -----------------------------------------------------------------------------
namespace sbndDB {
  /// The only thing this service does is to return its service provider of type
  /// `sbndDB::PMTCalibrationDatabase`.
  using IPMTCalibrationDatabaseService
    = lar::ServiceProviderInterfaceWrapper<PMTCalibrationDatabase>;
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE(sbndDB::IPMTCalibrationDatabaseService, SHARED)


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_DATABASEINTERFACE_IPMTCALIBRATIONDATABASESERVICE_H