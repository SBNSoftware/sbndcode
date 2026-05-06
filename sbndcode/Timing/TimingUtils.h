#ifndef TIMINGUTILS_H_SEEN
#define TIMINGUTILS_H_SEEN

/**
 * @file TimingUtils.h
 * @brief Common utility functions for timing reconstruction.
 */

#include "canvas/Persistency/Common/Ptr.h"

namespace sbnd {

  /**
   * @namespace TimingUtils
   * @brief Collection of helper utilities for timing calculations.
   */
  namespace TimingUtils{

    /**
     * @brief Subtracts two UTC timestamps.
     * @param[in] ts1 First timestamp in nanoseconds.
     * @param[in] ts2 Second timestamp in nanoseconds.
     * @return Difference between `ts1 - ts2` in nanoseconds as double.
     */
    double SubtractUTCTimestmap(uint64_t ts1, uint64_t ts2);

  }
}

#endif