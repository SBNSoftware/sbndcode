/// \file    SBNDPTBRawUtils.h
/// \brief   Standalone C++ methods to interact with the sbndptb data product, separated from
///          the data product defintion itself so ROOT can persist it
/// \author  trj@fnal.gov

#ifndef  SBNDPTBRawUtils_H
#define  SBNDPTBRawUtils_H

#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndaq-artdaq-core/Overlays/SBND/PTB_content.h"

namespace raw {
  namespace ptb {
    const std::vector<raw::ptb::ChStatus>  GetChStatusBeforeHLTs(const raw::ptb::sbndptb &pdata);
  }
}

#endif
