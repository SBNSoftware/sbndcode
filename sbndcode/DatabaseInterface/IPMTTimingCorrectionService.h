/**
 * @file   icaruscode/Timing/IPMTTimingCorrectionService.h
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Andrea Scarpelli (ascarpell@bnl.gov), Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H
#define ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H

// ICARUS libraries
#include "sbndcode/DatabaseInterface/PMTTimingCorrections.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceProviderWrappers.h"


// -----------------------------------------------------------------------------
namespace sbndDB {
  /// The only thing this service does is to return its service provider of type
  /// `icarusDB::PMTTimingCorrections`.
  using IPMTTimingCorrectionService
    = lar::ServiceProviderInterfaceWrapper<PMTTimingCorrections>;
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE(sbndDB::IPMTTimingCorrectionService, SHARED)


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H