////////////////////////////////////////////////////////////////////////
// Class:       DigitalNoiseChannelStatus
// Plugin Type: service (Unknown Unknown)
// File:        DigitalNoiseChannelStatus_service.cc
//
// Generated at Thu May  2 17:24:05 2024 by Thomas Junk using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#ifndef DigitalNoiseChannelStatusService_H
#define DigitalNoiseChannelStatusService_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"

#include <unordered_set>

namespace sbnd {
  class DigitalNoiseChannelStatus;
}


class sbnd::DigitalNoiseChannelStatus {
public:
  explicit DigitalNoiseChannelStatus(fhicl::ParameterSet const& p, art::ActivityRegistry& areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  bool IsBad(raw::ChannelID_t chan) const;

  const std::unordered_set<raw::ChannelID_t> & GetSetOfBadChannels() const;

  size_t NBadChannels() const;

  void pub_PrepEvent(const art::Event& evt, art::ScheduleContext);

private:

  std::unordered_set<raw::ChannelID_t> fDNChannels;   // set of channels with digital noise on them on this event
  float fRMSCutWire;
  float fRMSCutRawDigit;
  int fNBADCutRawDigit;
  std::string fRawDigitLabel;
  std::string fRecobWireLabel;

  int fNAwayFromPedestalRawDigit;
  int fDistFromPedestalRawDigit;
  int fNAwayFromPedestalRecobWire;
  float fDistFromPedestalRecobWire;
  
};

DECLARE_ART_SERVICE(sbnd::DigitalNoiseChannelStatus, LEGACY)

#endif
