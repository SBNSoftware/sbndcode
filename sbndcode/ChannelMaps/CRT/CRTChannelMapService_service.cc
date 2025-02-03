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
        >> m.feb_mac5
        >> m.channel_order_swapped;

      m.valid = true;

      fModuleInfoFromOfflineID[m.offline_module_id] = m;
      fModuleInfoFromFEBMAC5[m.feb_mac5]            = m;
    }
  
  inFile.close();
}

SBND::CRTChannelMapService::CRTChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : SBND::CRTChannelMapService(pset)
{
}

SBND::CRTChannelMapService::ModuleInfo_t SBND::CRTChannelMapService::GetModuleInfoFromFEBMAC5(unsigned int feb_mac5) const
{
  SBND::CRTChannelMapService::ModuleInfo_t bad;
  bad.valid = false;

  auto moduleIter = fModuleInfoFromFEBMAC5.find(feb_mac5);

  if(moduleIter == fModuleInfoFromFEBMAC5.end())
    {
      mf::LogInfo("SBND CRT Channel Map") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                          << "This FEB does not appear in the channel map." << std::endl;

      return bad;
    }

  return moduleIter->second;
}

SBND::CRTChannelMapService::ModuleInfo_t SBND::CRTChannelMapService::GetModuleInfoFromOfflineID(unsigned int offline_module_id) const
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

DEFINE_ART_SERVICE(SBND::CRTChannelMapService)
