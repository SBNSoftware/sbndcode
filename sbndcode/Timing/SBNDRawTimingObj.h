////////////////////////////////////////////////////////////////////////
//
// An offline object to associate original CAEN1730 timing information to 
// the decoded OpDetWaveforms to aid in downstream timing reconstruction
// author: Lynn Tung
// with inspiration taken from Tom Junk's offline ptb object
//
////////////////////////////////////////////////////////////////////////

#ifndef  sbndtiming_H
#define  sbndtiming_H

#include <stdint.h>
#include <vector>

namespace raw {
    class TimingReferenceInfo {
      public:
        TimingReferenceInfo() {}; // constructor
        TimingReferenceInfo(uint16_t timingType, 
                        uint16_t timingChannel) :
        timingType(timingType), timingChannel(timingChannel) {};
    
        uint16_t timingType; // e.g. SPECTDC = 0; PTB HLT = 1; CAEN-only = 3
        uint16_t timingChannel; // e.g. TDC ETRIG = 4; PTB BNB Beam+Light = 2
    };

    class FrameShiftInfo {
      public:
        FrameShiftInfo() {}; // constructor
        FrameShiftInfo(double _frameTdcCrtt1, double _frameTdcBes, double _frameTdcRwm, double _frameHltCrtt1, double _frameHltBeamGate) :
        frameTdcCrtt1(_frameTdcCrtt1), frameTdcBes(_frameTdcBes), frameTdcRwm(_frameTdcRwm), frameHltCrtt1(_frameHltCrtt1), frameHltBeamGate(_frameHltBeamGate) {};
    
        double frameTdcCrtt1;
        double frameTdcBes;
        double frameTdcRwm;
        double frameHltCrtt1;
        double frameHltBeamGate;
    };

  namespace pmt {
    class BoardTimingInfo {
      // Associate one of these to every opdetwaveform in the board/digitizer, one per trigger
      public:
        BoardTimingInfo() {}; // constructor
        BoardTimingInfo(uint16_t postPercent, 
                        std::vector<uint32_t> &triggerTimeTag) :
        postPercent(postPercent), triggerTimeTag(triggerTimeTag) {};

        uint16_t              postPercent; // # 0-100, represents a percentage
        std::vector<uint32_t> triggerTimeTag; // ns
    };
  }
}

#endif
