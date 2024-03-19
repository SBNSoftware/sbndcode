//////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTChannelMapService
// Module type: service
// File:        SBND::CRTChannelMapService_service.cc
// Author:      Linyan Wan, Mar 2024.  
//
// Implementation of hardware-offline channel mapping reading from a file.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "CRTChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

  
SBND::CRTChannelMapService::CRTChannelMapService(fhicl::ParameterSet const& pset) {

  bool useDB = pset.get<bool>("UseHWDB",false);
  bool useFile = pset.get<bool>("ReadMapFromFile",true);
  if (useDB && useFile) {
    throw cet::exception("SBND::CRTChannelMapService: UseHWDB and ReadMapFromFile are both true");
  }
  if (!useDB && !useFile) {
    throw cet::exception("SBND::CRTChannelMapService: UseHWDB and ReadMapFromFile are both false");
  }
  if (useFile) {
    std::string channelMapFile = pset.get<std::string>("FileName");

    std::string fullname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(channelMapFile, fullname);

    if (fullname.empty()) {
      std::cout << "SBND::CRTChannelMapService Input file " << channelMapFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }
    std::cout << "SBND CRT Channel Map: Building CRT wiremap from file " << channelMapFile << std::endl;
    std::ifstream inFile(fullname, std::ios::in);
    std::string line;

    while (std::getline(inFile,line)) {
      std::stringstream linestream(line);
    
      SBND::CRTChannelMapService::ChanInfo_t c;
      linestream 
	>> c.offlchan
	>> c.Wall
	>> c.MAC5
	>> c.Swap;

      c.valid = true;
		if (c.offlchan < 0 || c.MAC5 < 0) c.valid = false;

      fChanInfoFromMAC5[c.MAC5] = c;
      fChanInfoFromOfflChan[c.offlchan] = c;
    }
    inFile.close();
  }
  else
    {
      throw cet::exception("SBND:CRTChannelMapService: Database access to be implemented.");
    }
}

SBND::CRTChannelMapService::CRTChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) : SBND::CRTChannelMapService
													    (pset) {
}

SBND::CRTChannelMapService::ChanInfo_t SBND::CRTChannelMapService::GetChanInfoFromMAC5(unsigned int mac5) const {

  SBND::CRTChannelMapService::ChanInfo_t badinfo{};
  badinfo.valid = false;

  // look up one map at a time in order to handle cases where the item is not found
  // without throwing and catching exception which can make debugging hard
  
  auto fm1 = fChanInfoFromMAC5.find(mac5);
  if (fm1 == fChanInfoFromMAC5.end()) return badinfo;
  return fm1->second;
}


SBND::CRTChannelMapService::ChanInfo_t SBND::CRTChannelMapService::GetChanInfoFromOfflChan(unsigned int offlineChannel) const {

  SBND::CRTChannelMapService::ChanInfo_t badinfo{};
  badinfo.valid = false;
  auto fm = fChanInfoFromOfflChan.find(offlineChannel);
  if (fm == fChanInfoFromOfflChan.end()) return badinfo;
  return fm->second;
}

DEFINE_ART_SERVICE(SBND::CRTChannelMapService)
