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

  unsigned int ConstructOfflineChannelIDFromOfflineModuleIDAndStripNumber(const unsigned int offline_module_id, const unsigned int strip_number);

  unsigned int ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(const unsigned int offline_module_id, const unsigned int local_offline_channel);

  unsigned int ConstructOnlineChannelIDFromMAC5AndLocalOnlineChannel(const unsigned int mac5, const unsigned int local_online_channel);

  unsigned int GetOfflineModuleIDFromOfflineChannelID(const unsigned int offline_channel_id);

  unsigned int GetLocalOfflineChannelFromOfflineChannelID(const unsigned int offline_channel_id);

  unsigned int GetMAC5FromOnlineChannelID(const unsigned int online_channel_id);

  unsigned int GetLocalOnlineChannelFromOnlineChannelID(const unsigned int online_channel_id);

  bool MAC5IsValid(const unsigned int mac5);

  bool OfflineModuleIDIsValid(const unsigned int offline_module_id);

  unsigned int GetOfflineModuleIDFromMAC5(const unsigned int mac5);

  unsigned int GetMAC5FromOfflineModuleID(const unsigned int offline_module_id);

  unsigned int GetOfflineChannelIDFromOnlineChannelID(const unsigned int online_channel_id);

  unsigned int GetOnlineChannelIDFromOfflineChannelID(const unsigned int offline_channel_id);

  bool GetInversionFromOfflineModuleID(const unsigned int offline_module_id);

  unsigned int GetMAC5FromOfflineChannelID(const unsigned int offline_channel_id);

private:

  typedef struct ModuleInfo {
    unsigned int mac5;                  // The MAC5 address of the hardware front end board
    bool         channel_order_swapped; // Are the order of the channels inverted with respect to simulation (FEB direction)
    unsigned int offline_module_id;     // The module number in simulation
    bool         valid;
  } ModuleInfo_t;

  // look up channel info by offline module number
  std::unordered_map<unsigned int, ModuleInfo_t> fModuleInfoFromOfflineID;
  ModuleInfo_t GetModuleInfoFromOfflineModuleID(unsigned int offline_module_id) const;

  // look up channel info by MAC5
  std::unordered_map<unsigned int, ModuleInfo_t> fModuleInfoFromMAC5;
  ModuleInfo_t GetModuleInfoFromMAC5(unsigned int mac5) const;
};

DECLARE_ART_SERVICE(SBND::CRTChannelMapService, LEGACY)

#endif
