#ifndef RAW_NTB_SBNDNTB_HEADER_GUARD
#define RAW_NTB_SBNDNTB_HEADER_GUARD

#include <stdlib.h>

namespace raw {
  namespace ntb {
    struct sbndntb {
      uint64_t boardreader_timestamp;  // from artdaq::fragment
      uint32_t event_number;  // from fragment metadata
      uint32_t frame_number;  // from fragment metadata
      uint16_t sample_number;  // from fragment metadata
      uint16_t busy;          // from NTB header words
      uint16_t sixteen_mhz_remainder;  // from NTB header words
      uint16_t two_mhz_sample;    // from NTB header words
      uint32_t frame;     // from NTB header words
      uint32_t ntrig;     // from NTB header words
      uint16_t pmt_trig_data;  // from NTB trigger data words
      uint16_t pc;             // from NTB trigger data words
      uint16_t external;       // from NTB trigger data words
      uint16_t active;         // from NTB trigger data words
      uint16_t gate2;          // from NTB trigger data words
      uint16_t gate1;          // from NTB trigger data words
      uint16_t veto;           // from NTB trigger data words
      uint16_t calib;          // from NTB trigger data words
      uint16_t phase;          // from NTB trigger data words
      uint16_t gatefake;       // from NTB trigger data words
      uint16_t beamfake;       // from NTB trigger data words
      uint16_t spare1;         // from NTB trigger data words
      uint32_t trailer;        // from NTB trailer words
    };
  }
}

#endif
