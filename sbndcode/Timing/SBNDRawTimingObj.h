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

    //To save important timing info in unix timestamp format
    //Products made by FrameShift module
    class TimingInfo {
      public:
        TimingInfo() {}; // constructor
        TimingInfo(uint64_t _rawDAQHeaderTimestamp, uint64_t _tdcCrtt1, uint64_t _tdcBes, uint64_t _tdcRwm, uint64_t _tdcEtrig, uint64_t _hltCrtt1, uint64_t _hltEtrig, uint64_t _hltBeamGate) :
        rawDAQHeaderTimestamp(_rawDAQHeaderTimestamp), tdcCrtt1(_tdcCrtt1), tdcBes(_tdcBes), tdcRwm(_tdcRwm), tdcEtrig(_tdcEtrig), hltCrtt1(_hltCrtt1), hltEtrig(_hltEtrig), hltBeamGate(_hltBeamGate) {};

        uint64_t rawDAQHeaderTimestamp;
        uint64_t tdcCrtt1;
        uint64_t tdcBes;
        uint64_t tdcRwm;
        uint64_t tdcEtrig; //Global frame for timingType = 0
        uint64_t hltCrtt1;
        uint64_t hltEtrig; //Global frame for timingType = 1
        uint64_t hltBeamGate;
    };

    //Object to shift timing to a specific frame after being decoded (by defaul: TDC ETRIG ch4)
    //Products made by FrameShift module
    class FrameShiftInfo {
      public:
        FrameShiftInfo() {}; // constructor
        FrameShiftInfo(uint16_t _timingType, double _frameTdcCrtt1, double _frameTdcBes, double _frameTdcRwm, double _frameHltCrtt1, double _frameHltBeamGate, double _frameApplyAtCaf) :
        timingType(_timingType), frameTdcCrtt1(_frameTdcCrtt1), frameTdcBes(_frameTdcBes), frameTdcRwm(_frameTdcRwm), frameHltCrtt1(_frameHltCrtt1), frameHltBeamGate(_frameHltBeamGate), frameApplyAtCaf(_frameApplyAtCaf) {};
    
        uint16_t timingType; // e.g. SPECTDC = 0; PTB HLT = 1; CAEN-only = 3
        double frameTdcCrtt1;
        double frameTdcBes;
        double frameTdcRwm;
        double frameHltCrtt1;
        double frameHltBeamGate;
        double frameApplyAtCaf;
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
