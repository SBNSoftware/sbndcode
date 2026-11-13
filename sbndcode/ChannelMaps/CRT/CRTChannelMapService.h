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

  unsigned int ConstructOfflineChannelIDFromOfflineModuleIDAndStripNumber(const unsigned int offline_module_id, const unsigned int strip_number) const;

  unsigned int ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(const unsigned int offline_module_id, const unsigned int local_offline_channel) const;

  unsigned int ConstructOnlineChannelIDFromMAC5AndLocalOnlineChannel(const unsigned int mac5, const unsigned int local_online_channel) const;

  unsigned int GetOfflineModuleIDFromOfflineChannelID(const unsigned int offline_channel_id) const;

  unsigned int GetLocalOfflineChannelFromOfflineChannelID(const unsigned int offline_channel_id) const;

  unsigned int GetMAC5FromOnlineChannelID(const unsigned int online_channel_id) const;

  unsigned int GetLocalOnlineChannelFromOnlineChannelID(const unsigned int online_channel_id) const;

  bool MAC5IsValid(const unsigned int mac5) const;

  bool OfflineModuleIDIsValid(const unsigned int offline_module_id) const;

  unsigned int GetOfflineModuleIDFromMAC5(const unsigned int mac5) const;

  unsigned int GetMAC5FromOfflineModuleID(const unsigned int offline_module_id) const;

  unsigned int GetOfflineChannelIDFromOnlineChannelID(const unsigned int online_channel_id) const;

  unsigned int GetOnlineChannelIDFromOfflineChannelID(const unsigned int offline_channel_id) const;

  bool GetInversionFromOfflineModuleID(const unsigned int offline_module_id) const;

  unsigned int GetMAC5FromOfflineChannelID(const unsigned int offline_channel_id) const;

  unsigned int GetTimingChainFromOfflineModuleID(unsigned int offline_module_id) const;

  unsigned int GetTimingChainFromOfflineChannelID(unsigned int offline_channel_id) const;

  unsigned int GetTimingChainFromMAC5(unsigned int mac5) const;

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

  // look up timing chain by offline module number
  std::unordered_map<unsigned int, unsigned int> fTimingChainFromOfflineID;

  // look up channel info by MAC5
  std::unordered_map<unsigned int, unsigned int> fTimingChainFromMAC5;
};

// The maps are filled once in the constructor and only read afterwards, so the service can be shared
// between threads: every lookup is const and none of them modifies the maps.
DECLARE_ART_SERVICE(SBND::CRTChannelMapService, SHARED)

#endif
