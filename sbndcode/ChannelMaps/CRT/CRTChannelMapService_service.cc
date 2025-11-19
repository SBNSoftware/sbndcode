//////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTChannelMapService
// Module type: service
// File:        SBND::CRTChannelMapService_service.cc
// Author:      Henry Lay, May 2024.
//
// Implementation of hardware-offline channel mapping reading from a file.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "CRTChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
SBND::CRTChannelMapService::CRTChannelMapService(fhicl::ParameterSet const& pset)
{
  const std::string channelMapFile = pset.get<std::string>("FileName");

  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);

  if(fullname.empty())
    {
      std::cout << "SBND::CRTChannelMapService Input file " << channelMapFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  std::cout << "SBND CRT Channel Map: Building map from file " << channelMapFile << std::endl;
  std::ifstream inFile(fullname, std::ios::in);
  std::string line;

  while(std::getline(inFile,line))
    {
      std::stringstream linestream(line);
    
      SBND::CRTChannelMapService::ModuleInfo_t m;
      linestream 
        >> m.offline_module_id
        >> m.mac5
        >> m.channel_order_swapped;

      m.valid = true;

      fModuleInfoFromOfflineID[m.offline_module_id] = m;
      fModuleInfoFromMAC5[m.mac5]                   = m;
    }
  
  inFile.close();
}

SBND::CRTChannelMapService::ModuleInfo_t SBND::CRTChannelMapService::GetModuleInfoFromMAC5(unsigned int mac5) const
{
  SBND::CRTChannelMapService::ModuleInfo_t bad;
  bad.valid = false;

  auto moduleIter = fModuleInfoFromMAC5.find(mac5);

  if(moduleIter == fModuleInfoFromMAC5.end())
    {
      mf::LogInfo("SBND CRT Channel Map") << "Asked for FEB with MAC5: " << mac5 << '\n'
                                          << "This FEB does not appear in the channel map." << std::endl;

      return bad;
    }

  return moduleIter->second;
}

SBND::CRTChannelMapService::ModuleInfo_t SBND::CRTChannelMapService::GetModuleInfoFromOfflineModuleID(unsigned int offline_module_id) const
{
  SBND::CRTChannelMapService::ModuleInfo_t bad;
  bad.valid = false;

  auto moduleIter = fModuleInfoFromOfflineID.find(offline_module_id);

  if(moduleIter == fModuleInfoFromOfflineID.end())
    {
      mf::LogInfo("SBND CRT Channel Map") << "Asked for module with offline ID: " << offline_module_id << '\n'
                                          << "This module does not appear in the channel map." << std::endl;

      return bad;
    }

  return moduleIter->second;
}

unsigned int SBND::CRTChannelMapService::ConstructOfflineChannelIDFromOfflineModuleIDAndStripNumber(const unsigned int offline_module_id, const unsigned int strip_number)
{
  return 32 * offline_module_id + 2 * strip_number;
}

unsigned int SBND::CRTChannelMapService::ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(const unsigned int offline_module_id, const unsigned int local_offline_channel)
{
  return 32 * offline_module_id + local_offline_channel;
}

unsigned int SBND::CRTChannelMapService::ConstructOnlineChannelIDFromMAC5AndLocalOnlineChannel(const unsigned int mac5, const unsigned int local_online_channel)
{
  return 100 * mac5 + local_online_channel;
}

unsigned int SBND::CRTChannelMapService::GetOfflineModuleIDFromOfflineChannelID(const unsigned int offline_channel_id)
{
  return offline_channel_id / 32;
}

unsigned int SBND::CRTChannelMapService::GetLocalOfflineChannelFromOfflineChannelID(const unsigned int offline_channel_id)
{
  return offline_channel_id % 32;
}

unsigned int SBND::CRTChannelMapService::GetMAC5FromOnlineChannelID(const unsigned int online_channel_id)
{
  return online_channel_id / 100;
}

unsigned int SBND::CRTChannelMapService::GetLocalOnlineChannelFromOnlineChannelID(const unsigned int online_channel_id)
{
  return online_channel_id % 100;
}

bool SBND::CRTChannelMapService::MAC5IsValid(const unsigned int mac5)
{
  ModuleInfo_t moduleinfo = GetModuleInfoFromMAC5(mac5);

  return moduleinfo.valid;
}

bool SBND::CRTChannelMapService::OfflineModuleIDIsValid(const unsigned int offline_module_id)
{
  ModuleInfo_t moduleinfo = GetModuleInfoFromOfflineModuleID(offline_module_id);

  return moduleinfo.valid;
}

unsigned int SBND::CRTChannelMapService::GetOfflineModuleIDFromMAC5(const unsigned int mac5)
{
  ModuleInfo_t moduleinfo = GetModuleInfoFromMAC5(mac5);

  if(!moduleinfo.valid)
    return 0;

  return moduleinfo.offline_module_id;
}

unsigned int SBND::CRTChannelMapService::GetMAC5FromOfflineModuleID(const unsigned int offline_module_id)
{
  ModuleInfo_t moduleinfo = GetModuleInfoFromOfflineModuleID(offline_module_id);

  if(!moduleinfo.valid)
    return 0;

  return moduleinfo.mac5;
}

unsigned int SBND::CRTChannelMapService::GetOfflineChannelIDFromOnlineChannelID(const unsigned int online_channel_id)
{
  const unsigned int mac5                 = GetMAC5FromOnlineChannelID(online_channel_id);
  const unsigned int local_online_channel = GetLocalOnlineChannelFromOnlineChannelID(online_channel_id);

  ModuleInfo_t moduleinfo = GetModuleInfoFromMAC5(mac5);

  if(!moduleinfo.valid)
    return 0;

  const unsigned int local_offline_channel = moduleinfo.channel_order_swapped ? 31 - local_online_channel : local_online_channel;

  return ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(moduleinfo.offline_module_id, local_offline_channel);
}

unsigned int SBND::CRTChannelMapService::GetOnlineChannelIDFromOfflineChannelID(const unsigned int offline_channel_id)
{
  const unsigned int offline_module_id     = GetOfflineModuleIDFromOfflineChannelID(offline_channel_id);
  const unsigned int local_offline_channel = GetLocalOfflineChannelFromOfflineChannelID(offline_channel_id);

  ModuleInfo_t moduleinfo = GetModuleInfoFromOfflineModuleID(offline_module_id);

  if(!moduleinfo.valid)
    return 0;

  const unsigned int local_online_channel = moduleinfo.channel_order_swapped ? 31 - local_offline_channel : local_offline_channel;

  return ConstructOnlineChannelIDFromMAC5AndLocalOnlineChannel(moduleinfo.mac5, local_online_channel);
}

bool SBND::CRTChannelMapService::GetInversionFromOfflineModuleID(const unsigned int offline_module_id)
{
  ModuleInfo_t moduleinfo = GetModuleInfoFromOfflineModuleID(offline_module_id);

  if(!moduleinfo.valid)
    return false;

  return moduleinfo.channel_order_swapped;
}

unsigned int SBND::CRTChannelMapService::GetMAC5FromOfflineChannelID(const unsigned int offline_channel_id)
{
  return GetMAC5FromOfflineModuleID(GetOfflineModuleIDFromOfflineChannelID(offline_channel_id));
}

DEFINE_ART_SERVICE(SBND::CRTChannelMapService)
