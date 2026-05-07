#include "sbndcode/Timing/TimingUtils.h"

namespace sbnd {
  double TimingUtils::SubtractUTCTimestmap(uint64_t ts1, uint64_t ts2)
  {
    return (ts1 > ts2)? (double)(ts1 - ts2): -(double)(ts2 - ts1);
  }
}