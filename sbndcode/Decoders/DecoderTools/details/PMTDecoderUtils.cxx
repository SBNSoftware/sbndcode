/**
 * @file   sbndcode/Decoders/DecoderTools/details/PMTDecoderUtils.h
 * @brief  Some helpers for PMT decoder tool.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date  Jan 28, 2023
 */


// library header
#include "sbndcode/Decoders/DecoderTools/details/PMTDecoderUtils.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
#include <ostream>
#include <cassert>


// -----------------------------------------------------------------------------
std::ostream& daq::details::operator<<
  (std::ostream& out, BoardSetup_t::ChannelSetup_t const& chSetup)
{
  out << "channel ";
  if (chSetup.hasChannel()) {
    if (chSetup.category.empty()) out << chSetup.channelID;
    else out << "0x" << std::hex << chSetup.channelID << std::dec;
  }
  else out << "from channel DB";
  
  if (chSetup.mustSkip()) out << ", always skipped";
  else {
    if (chSetup.mustSave() || chSetup.hasChannel()) {
      if (chSetup.mustSave() || (chSetup.minSpan == 0))
        out << ", always saved";
      else
        out << ", saved if span larger than " << chSetup.minSpan << " ADC";
    }
    else {
      out << ", saved if in channel DB";
      if (chSetup.minSpan > 0)
        out << " and if span larger than " << chSetup.minSpan << " ADC";
    }
    
    if (chSetup.onGlobalOnly) out << " only if on global trigger";
    
    if (!chSetup.category.empty())
      out << " in special category '" << chSetup.category << "'";
    
  } // if not forced to skip
  
  return out;
} // daq::details::operator<< (daq::details::BoardSetup_t::ChannelSetup_t)


// -----------------------------------------------------------------------------
std::ostream& daq::details::operator<<
  (std::ostream& out, BoardInfoLookup const& db)
{
  
  out << "Information on " << db.nBoardInfo() << " boards recorded:";
  for (BoardInfoLookup::BoardInfo_t const& boardInfo: db.allBoardInfo()) {
    assert(boardInfo.setup);
    out << "\n  board \"" << boardInfo.setup->name
      << "\" (fragment ID " << std::hex << boardInfo.fragmentID << std::dec
      << "): trigger delay " << boardInfo.setup->triggerDelay
      << ", TTT reset delay " << boardInfo.setup->TTTresetDelay
      << ", pre-trigger buffer length " << boardInfo.facts.preTriggerTime;
    if (boardInfo.config) {
      out << ", buffer " << boardInfo.config->bufferLength
        << " tick long, board ID " << boardInfo.config->boardID
        << " with " << boardInfo.config->nChannels << " channels";
    }
    else {
      out << " (no PMT configuration)";
    }
    
    std::vector<std::size_t> specialChannels;
    for (auto const& [ iChannel, chSetup ]
      : util::enumerate(boardInfo.setup->channelSettings)
    ) {
      if (!chSetup.isDefault()) specialChannels.push_back(iChannel);
    }
    if (specialChannels.empty()) continue;
    out << "; " << specialChannels.size() << " special channel configurations:";
    for (std::size_t const iCh: specialChannels) {
      out << "\n    ch.number " << iCh << ": "
        << boardInfo.setup->channelSettings[iCh];
    }
    
  } // for board info
  out << '\n';
  return out;
} // daq::details::operator<< (daq::details::BoardInfoLookup)


// -----------------------------------------------------------------------------
