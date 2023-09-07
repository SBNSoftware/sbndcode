//////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::TPCChannelMapService
// Module type: service
// File:        SBND::TPCChannelMapService_service.cc
// Author:      Tom Junk, August 2023.  Input file from Nupur Oza
//
// Implementation of hardware-offline channel mapping reading from a file.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "TPCChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

  
SBND::TPCChannelMapService::TPCChannelMapService(fhicl::ParameterSet const& pset) {

  bool useDB = pset.get<bool>("UseHWDB",false);
  bool useFile = pset.get<bool>("ReadMapFromFile",true);
  if (useDB && useFile) {
    throw cet::exception("SBND::TPCChannelMapService: UseHWDB and ReadMapFromFile are both true");
  }
  if (!useDB && !useFile) {
    throw cet::exception("SBND::TPCChannelMapService: UseHWDB and ReadMapFromFile are both false");
  }
  if (useFile) {
    std::string channelMapFile = pset.get<std::string>("FileName");

    std::string fullname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(channelMapFile, fullname);

    if (fullname.empty()) {
      std::cout << "SBND::TPCChannelMapService Input file " << channelMapFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }
    std::cout << "SBND TPC Channel Map: Building TPC wiremap from file " << channelMapFile << std::endl;
    std::ifstream inFile(fullname, std::ios::in);
    std::string line;

    while (std::getline(inFile,line)) {
      std::stringstream linestream(line);
      std::string planestr;
      std::string qfspstr;
    
      SBND::TPCChannelMapService::ChanInfo_t c;
      linestream 
	>> c.wireno
	>> planestr
	>> c.EastWest
	>> c.NorthSouth
	>> c.SideTop
	>> c.FEMBPosition
	>> c.FEMBSerialNum
	>> c.FEMBOnWIB
	>> c.FEMBCh
	>> c.asic
	>> c.WIBCrate
	>> c.WIB
	>> c.WIBCh 
	>> qfspstr 
	>> c.QFSPFiber 
	>> c.FEMCrate
	>> c.FEM
	>> c.FEMCh
	>> c.offlchan;

      c.valid = true;
      c.plane = 10;
      if (planestr == "U") c.plane = 0;
      if (planestr == "V") c.plane = 1;
      if (planestr == "Y") c.plane = 2;
      if (c.plane == 10) c.valid = false;
      c.WIBQFSP = atoi(qfspstr.substr(3,1).c_str());

      fChanInfoFromFEMInfo[c.FEMCrate][c.FEM][c.FEMCh] = c;
      fChanInfoFromOfflChan[c.offlchan] = c;
    }
    inFile.close();
  }
  else
    {
      throw cet::exception("SBND:TPCChannelMapService: Database access to be implemented.");
    }
}

SBND::TPCChannelMapService::TPCChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) : SBND::TPCChannelMapService
													    (pset) {
}

SBND::TPCChannelMapService::ChanInfo_t SBND::TPCChannelMapService::GetChanInfoFromFEMElements(unsigned int femcrate,
											      unsigned int fem,
											      unsigned int femchan) const {

  SBND::TPCChannelMapService::ChanInfo_t badinfo{};
  badinfo.valid = false;

  // look up one map at a time in order to handle cases where the item is not found
  // without throwing and catching exception which can make debugging hard
  
  auto fm1 = fChanInfoFromFEMInfo.find(femcrate);
  if (fm1 == fChanInfoFromFEMInfo.end()) {
    unsigned int substituteCrate = 1;  // a hack -- ununderstood crates get mapped to crate 1
    fm1 = fChanInfoFromFEMInfo.find(substituteCrate);
    if (fm1 == fChanInfoFromFEMInfo.end()) {
      return badinfo;
    }
  }
  auto& m1 = fm1->second;
  auto fm2 = m1.find(fem);
  if (fm2 == m1.end()) return badinfo;
  auto& m2 = fm2->second;
  auto fm3 = m2.find(femchan);
  if (fm3 == m2.end()) return badinfo;    
  return fm3->second;
}


SBND::TPCChannelMapService::ChanInfo_t SBND::TPCChannelMapService::GetChanInfoFromOfflChan(unsigned int offlineChannel) const {

  SBND::TPCChannelMapService::ChanInfo_t badinfo{};
  badinfo.valid = false;
  auto fm = fChanInfoFromOfflChan.find(offlineChannel);
  if (fm == fChanInfoFromOfflChan.end()) return badinfo;
  return fm->second;
}

DEFINE_ART_SERVICE(SBND::TPCChannelMapService)
