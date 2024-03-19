///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTChannelMapService
// Module type: service
// File:        CRTChannelMapService.h
// Author:      Linyan Wan, Mar 2024
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

  typedef struct ChanInfo {
    std::string  Wall;        // "East", "North", "Bottom"
    unsigned int Swap;           // same=0, swapped=1
    unsigned int MAC5;           // MAC5
    unsigned int offlchan;        // in gdml and channel sorting convention
    bool valid;          // true if valid, false if not
  } ChanInfo_t;

  ChanInfo_t GetChanInfoFromMAC5(
					unsigned int mac5) const;

  ChanInfo_t GetChanInfoFromOfflChan(unsigned int offlchan) const;

private:

  // look up channel info by offline channel number
  
  std::unordered_map<unsigned int, ChanInfo_t> fChanInfoFromOfflChan;

  // look up channel info by FEMCrate, FEM, and FEMCh
  
  std::unordered_map< unsigned int, ChanInfo_t > fChanInfoFromMAC5;

};


DECLARE_ART_SERVICE(SBND::CRTChannelMapService, LEGACY)

#endif
