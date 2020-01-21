// AdcNoiseRemovalService.h
#ifndef AdcNoiseRemovalService_H
#define AdcNoiseRemovalService_H
// David Adams
// May 2016
//
// Interface for a service that removes noise from a collection of ADC
// channel samples. A map of samples are passed so that coherent noise
// may be evaluated and removed.
#include <iostream>
#include "sbndcode/Calibration/AdcChannelData.h"
class AdcNoiseRemovalService {
 public:
  // Remove the noise for selected channels from datamap[chan].samples and store results
  // in that same vector.
  // Return 0 for success.
  virtual int update(AdcChannelDataMap& datamap) const =0;
  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;
};
#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(AdcNoiseRemovalService, LEGACY)
#endif
#endif
