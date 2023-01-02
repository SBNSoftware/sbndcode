/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/IPMTTimingCorrectionService.h
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_IPMTTIMINGCORRECTIONSERVICE_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_IPMTTIMINGCORRECTIONSERVICE_H

// SBND libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrections.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceProviderWrappers.h"


// -----------------------------------------------------------------------------
namespace sbndDB {
  /// The only thing this service does is to return its service provider of type
  /// `sbndDB::PMTTimingCorrections`.
  using IPMTTimingCorrectionService
    = lar::ServiceProviderInterfaceWrapper<PMTTimingCorrections>;
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE(sbndDB::IPMTTimingCorrectionService, SHARED)


// -----------------------------------------------------------------------------


#endif // SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_IPMTTIMINGCORRECTIONSERVICE_H
