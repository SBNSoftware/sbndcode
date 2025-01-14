///////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::PMTChannelMapService
// Module type: service
// File:        PMTChannelMapService.h
// Author:      Henry Lay, May 2024
//
// Implementation of hardware-offline channel mapping reading from a file or from the hardware database 
// SBND CRT
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SBNDPMTChannelMapService_H
#define SBNDPMTChannelMapService_H

#include <unordered_map>
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace SBND {
  class PMTChannelMapService;
}

  
class SBND::PMTChannelMapService {

public:

  PMTChannelMapService(fhicl::ParameterSet const& pset);
  PMTChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  typedef struct PMTInfo {
    unsigned int channelID;
    std::string installationID;
    unsigned int breakoutBoardID;
    unsigned int slotID;
    unsigned int digitiserBoardID;
    unsigned int digitiserChannel;
    double warmCableTransit; //ns
    double PMTtoBBDelay; //ns
    double PMTtoBBTransit; //ns
    double BBtoDigitiserDelay; //ns
    double BBtoDigitiserTransit; //ns
    double TotalTransit; //ns
    unsigned int dqmID;
    bool valid;
  } PMTInfo_t;

  PMTInfo_t GetPMTInfoFromChannelID(unsigned int pmt_ch_id) const;

private:

  // look up channel info by offline channel ID
  
  std::unordered_map<unsigned int, PMTInfo_t> fPMTInfoFromChannelID;

};

DECLARE_ART_SERVICE(SBND::PMTChannelMapService, LEGACY)

#endif
