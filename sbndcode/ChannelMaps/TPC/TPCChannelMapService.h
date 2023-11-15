///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::TPCChannelMapService
// Module type: service
// File:        TPCChannelMapService.h
// Author:      Tom Junk and Nupur Oza, August 2023
//
// Implementation of hardware-offline channel mapping reading from a file or from the hardware database 
// SBND TPC
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SBNDTPCChannelMapService_H
#define SBNDTPCChannelMapService_H

#include <unordered_map>
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace SBND {
  class TPCChannelMapService;
}

  
class SBND::TPCChannelMapService {

public:

  TPCChannelMapService(fhicl::ParameterSet const& pset);
  TPCChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  typedef struct ChanInfo {
    unsigned int wireno;          // wire number
    unsigned int plane;           // 0: U,  1: V,  2: Y
    std::string  EastWest;        // "East" or "West"
    std::string  NorthSouth;      // "North" or "South"
    std::string  SideTop;         // "T" or "S"
    std::string  FEMBPosition;    // e.g. "A10"
    std::string  FEMBSerialNum;   // Serial number.  Some contain letters
    unsigned int FEMBOnWIB;       // 0:3
    unsigned int FEMBCh;          // channel on FEMB -- 0:127
    unsigned int asic;            // 0:7
    unsigned int asicchan;        // ASIC channel:  0 to 15
    unsigned int WIBCrate;        // 1:4
    unsigned int WIB;             // 1:6
    unsigned int WIBCh;           // 0:895   (7 FEMBs)
    unsigned int WIBQFSP;         // 1:2
    unsigned int QFSPFiber;       // 1:4
    unsigned int FEMCrate;        // 1:11
    unsigned int FEM;             // 1:16
    unsigned int FEMCh;           // 0:63 channel in a FEM
    unsigned int offlchan;        // in gdml and channel sorting convention
    bool valid;          // true if valid, false if not
  } ChanInfo_t;

  ChanInfo_t GetChanInfoFromFEMElements(
					unsigned int femcrate,
					unsigned int fem,
					unsigned int femchan) const;

  ChanInfo_t GetChanInfoFromOfflChan(unsigned int offlchan) const;

private:

  // look up channel info by offline channel number
  
  std::unordered_map<unsigned int, ChanInfo_t> fChanInfoFromOfflChan;

  // look up channel info by FEMCrate, FEM, and FEMCh
  
  std::unordered_map<unsigned int,
		     std::unordered_map<unsigned int,
					std::unordered_map< unsigned int, ChanInfo_t > > > fChanInfoFromFEMInfo;

};


DECLARE_ART_SERVICE(SBND::TPCChannelMapService, LEGACY)

#endif
