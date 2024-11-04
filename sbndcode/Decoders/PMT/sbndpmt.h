////////////////////////////////////////////////////////////////////////
//
// An offline object to associate original CAEN1730 timing information to 
// the decoded OpDetWaveforms to aid in downstream timing reconstruction
// author: Lynn Tung
// with inspiration taken from Tom Junk's offline ptb object
//
////////////////////////////////////////////////////////////////////////

#ifndef  sbndpmt_H
#define  sbndpmt_H

#include <stdint.h>
#include <vector>

namespace raw {
  namespace pmt {
    class eventTimingInfo {
      public:
        eventTimingInfo() {}; // constructor
        eventTimingInfo(uint16_t timingType, 
                        uint16_t timingChannel) :
        timingType(timingType), timingChannel(timingChannel) {};
    
        uint16_t timingType; // e.g. SPECTDC = 0; PTB HLT = 1; CAEN-only = 3
        uint16_t timingChannel; // e.g. TDC ETRIG = 4; PTB BNB Beam+Light = 2
    };

    class boardTimingInfo {
      // Associate one of these to every opdetwaveform in the board/digitizer, one per trigger
      public:
        boardTimingInfo() {}; // constructor
        boardTimingInfo(uint16_t postPercent, 
                        uint32_t triggerTimeTag) :
        postPercent(postPercent), triggerTimeTag(triggerTimeTag) {};

        uint16_t    postPercent; // # 0-100, represents a percentage
        uint32_t triggerTimeTag; // ns
    };
  }
}

#endif