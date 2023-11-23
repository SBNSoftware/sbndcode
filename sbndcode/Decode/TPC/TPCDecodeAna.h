#ifndef _sbnddaq_analysis_HeaderData
#define _sbnddaq_analysis_HeaderData

#include <string>
#include <iostream>
#include <sstream> 
#include <stdlib.h>
#include <ctime>

namespace tpcAnalysis {

// HeaderData: assumes the each board is 
// generating one fragment, and the each fragment has one header
// which is associated with a seuqence of ADC counts 
class HeaderData {
  public:
  uint8_t crate; //!< Index of readout electronics crate
  uint8_t slot; //!< Index of "slot" of readout board 
  uint32_t event_number; //!< Event number for this header 
  uint32_t checksum; //!< checksum associated with header
  uint32_t timestamp; //!< timestamp for this header
  uint32_t adc_word_count; //!< word count of ADC counts associated with this header

  unsigned index; //!< Globally usable index for this header information. 
    // No two header objects in the same event shold have the same index

  // by default make words noticable
  // Nevis uses DEADBEEF as a default, so distinguish from
  // that use BEEFDEAD
  HeaderData():
    crate(0xFF),
    slot(0xFF),
    event_number(0xBEEFDEAD),
    checksum(0xBEEFDEAD),
    timestamp(0xBEEFDEAD),
    adc_word_count(0xBEEFDEAD)
  {}

  // print the data -- for debugging
  std::string Print() const {
    std::stringstream buffer;
    buffer << "crate: " << ((unsigned)crate) << std::endl;
    buffer << "slot: " << ((unsigned)slot) << std::endl;
    buffer << "event no: " << event_number << std::endl;
    buffer << "checksum: " << checksum << std::endl;
    buffer << "adc word count: " << adc_word_count << std::endl;
    buffer << "timestamp: " << timestamp << std::endl;
    buffer << std::hex << "HEX checksum: " << checksum << std::endl;

    return buffer.str();
  }

};
}

#endif
