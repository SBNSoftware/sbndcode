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
#include "sbnobj/SBND/CRT/CRTEnums.hh"

namespace SBND {
  class CRTCalibService;
}

  
class SBND::CRTCalibService {

public:

  CRTCalibService(fhicl::ParameterSet const& pset);
  CRTCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  double GetT0CableOffsetFromFEBMAC5(unsigned int feb_mac5) const;

  double GetT1CableOffsetFromFEBMAC5(unsigned int feb_mac5) const;

  double GetT0CalibOffsetFromFEBMAC5(unsigned int feb_mac5) const;

  double GetT1CalibOffsetFromFEBMAC5(unsigned int feb_mac5) const;

  double GetPedestalFromFEBMAC5AndChannel(unsigned int feb_mac5, unsigned int ch) const;

  double GetGainFromFEBMAC5AndChannel(unsigned int feb_mac5, unsigned int ch) const;

  enum sbnd::crt::CRTChannelStatus GetChannelStatusFromFEBMAC5AndChannel(unsigned int feb_mac5,
                                                                         unsigned int ch) const;

private:

  std::unordered_map<unsigned int, double> fT0CableOffsetFromFEBMAC5, fT1CableOffsetFromFEBMAC5,
    fT0CalibOffsetFromFEBMAC5, fT1CalibOffsetFromFEBMAC5;

  std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> fPedestalFromFEBMAC5AndChannel;

  std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> fGainFromFEBMAC5AndChannel;

  std::unordered_map<unsigned int, std::unordered_map<unsigned int, sbnd::crt::CRTChannelStatus>> fChannelStatusFromFEBMAC5AndChannel;

};

DECLARE_ART_SERVICE(SBND::CRTCalibService, LEGACY)

#endif
