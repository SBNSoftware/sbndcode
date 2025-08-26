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

    class BoardAlignment {
      // Associate one of these to every opdetwaveform in the board/digitizer, one per board
      std::vector<double> fShift; //ns
      std::vector<int> fStatus;

      public:
        BoardAlignment() : fShift({}), fStatus({}) {}; // constructor
        virtual ~BoardAlignment() {}; // destructror
       
        BoardAlignment(std::vector<double> &_shift, std::vector<int> &_status) : fShift(_shift), fStatus(_status) {};

        std::vector<double> Shift() const {return fShift;}
        std::vector<int> Status() const {return fStatus;}
    };
  }
}

#endif
