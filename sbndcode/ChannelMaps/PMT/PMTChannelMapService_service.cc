//////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::PMTChannelMapService
// Module type: service
// File:        SBND::PMTChannelMapService_service.cc
// Author:      Lan Nguyen, Dec 2024.
//
// Implementation of hardware-offline channel mapping reading from a file.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "PMTChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
SBND::PMTChannelMapService::PMTChannelMapService(fhicl::ParameterSet const& pset)
{
  const std::string channelMapFile = pset.get<std::string>("FileName");

  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);

  if(fullname.empty())
    {
      std::cout << "SBND::PMTChannelMapService Input file " << channelMapFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  std::cout << "SBND PMT Channel Map: Building map from file " << channelMapFile << std::endl;
  std::ifstream inFile(fullname, std::ios::in);
  std::string line;

  //Read the first header line
  std::getline(inFile,line);

  while(std::getline(inFile,line))
    {
      std::stringstream linestream(line);
    
      SBND::PMTChannelMapService::PMTInfo_t m;
      linestream 
        >> m.channelID
        >> m.installationID
        >> m.breakoutBoardID
        >> m.slotID
        >> m.digitiserBoardID
        >> m.digitiserChannel
        >> m.warmCableTransit
        >> m.PMTtoBBDelay
        >> m.PMTtoBBTransit
        >> m.BBtoDigitiserDelay
        >> m.BBtoDigitiserTransit
        >> m.TotalTransit
        >> m.dqmID;

      m.valid = true;

      fPMTInfoFromChannelID[m.channelID] = m;
    }
  
  inFile.close();
}

SBND::PMTChannelMapService::PMTChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : SBND::PMTChannelMapService(pset)
{
}

SBND::PMTChannelMapService::PMTInfo_t SBND::PMTChannelMapService::GetPMTInfoFromChannelID(unsigned int pmt_ch_id) const
{
  SBND::PMTChannelMapService::PMTInfo_t bad;
  bad.valid = false;

  auto moduleIter = fPMTInfoFromChannelID.find(pmt_ch_id);

  if(moduleIter == fPMTInfoFromChannelID.end())
    {
      mf::LogInfo("SBND PMT Channel Map") << "Asked for pmt with channel ID: " << pmt_ch_id << '\n'
                                          << "This pmt does not appear in the channel map." << std::endl;

      return bad;
    }

  return moduleIter->second;
}

DEFINE_ART_SERVICE(SBND::PMTChannelMapService)
