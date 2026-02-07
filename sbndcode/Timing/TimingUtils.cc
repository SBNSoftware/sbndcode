#include "sbndcode/Timing/TimingUtils.h"

namespace sbnd {
  double TimingUtils::SubtractUTCTimestmap(const uint64_t& ts1, const uint64_t& ts2)
  {
  
    double ts1_s = ts1 / uint64_t(1e9);
    double ts1_ns = ts1 % uint64_t(1e9);
    double ts2_s = ts2 / uint64_t(1e9);
    double ts2_ns = ts2 % uint64_t(1e9);
    
    double diff_s = 0;
    double diff_ns = 0;
    
    //If the same PPS, just subtract the nanoseconds
    if(ts1_s == ts2_s){
      diff_ns = ts1_ns - ts2_ns;
    }
    //If ts1 is later than ts2, then subtract the seconds and add the nanoseconds
    else{
      diff_s = ts1_s - ts2_s;
      diff_ns = diff_s * (double)1e9 + ts1_ns - ts2_ns;
    }
  
    return diff_ns;
  }
}
