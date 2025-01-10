///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTChannelMapService
// Module type: service
// File:        CRTChannelMapService.h
// Author:      Henry Lay, May 2024
//
// Implementation of hardware-offline channel mapping reading from a file or from the hardware database 
// SBND CRT
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SBNDCRTChannelMapService_H
#define SBNDCRTChannelMapService_H

#include <unordered_map>
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace SBND {
  class CRTChannelMapService;
}

  
class SBND::CRTChannelMapService {

public:

  CRTChannelMapService(fhicl::ParameterSet const& pset);
  CRTChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  typedef struct ModuleInfo {
    unsigned int feb_mac5;              // The MAC5 address of the hardware front end board
    bool         channel_order_swapped; // Are the order of the channels inverted with respect to simulation (FEB direction)
    unsigned int offline_module_id;     // The module number in simulation
    bool         valid;
  } ModuleInfo_t;

  ModuleInfo_t GetModuleInfoFromFEBMAC5(unsigned int feb_mac5) const;

  ModuleInfo_t GetModuleInfoFromOfflineID(unsigned int offline_module_id) const;

private:

  // look up channel info by offline module number
  
  std::unordered_map<unsigned int, ModuleInfo_t> fModuleInfoFromOfflineID;

  // look up channel info by FEB MAC5
  
  std::unordered_map<unsigned int, ModuleInfo_t> fModuleInfoFromFEBMAC5;

};

DECLARE_ART_SERVICE(SBND::CRTChannelMapService, LEGACY)

#endif
