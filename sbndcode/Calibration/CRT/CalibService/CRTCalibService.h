///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTCalibService
// Module type: service
// File:        CRTCalibService.h
// Author:      Henry Lay, May 2024
//
// Service for inputing CRT calibration values to reconstruction
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SBNDCRTCalibService_H
#define SBNDCRTCalibService_H

#include <unordered_map>
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace SBND {
  class CRTCalibService;
}

  
class SBND::CRTCalibService {

public:

  CRTCalibService(fhicl::ParameterSet const& pset);
  CRTCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  double GetTimingOffsetFromFEBMAC5(unsigned int feb_mac5) const;

  double GetPedestalFromFEBMAC5AndChannel(unsigned int feb_mac5, unsigned int ch) const;

private:

  std::unordered_map<unsigned int, double> fTimingOffsetFromFEBMAC5;

  std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> fPedestalFromFEBMAC5AndChannel;

};

DECLARE_ART_SERVICE(SBND::CRTCalibService, LEGACY)

#endif
