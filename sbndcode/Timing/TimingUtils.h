#ifndef TIMINGUTILS_H_SEEN
#define TIMINGUTILS_H_SEEN

///////////////////////////////////////////////
// TimingUtils.h
//
// Common functions for Timing reconstruction
///////////////////////////////////////////////

#include "canvas/Persistency/Common/Ptr.h"

namespace sbnd {

  namespace TimingUtils{
    
    //Subtract two timestamps in UTC format
    //ts1 and ts2 are in nanoseconds in uint64_t format
    //Return the difference in nanoseconds as double
    double SubtractUTCTimestmap(const uint64_t& ts1, const uint64_t& ts2);

  }
}

#endif