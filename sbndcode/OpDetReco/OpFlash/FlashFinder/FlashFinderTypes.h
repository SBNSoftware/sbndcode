#ifndef __FLASHFINDERTYPES_H__
#define __FLASHFINDERTYPES_H__

#include <vector>
#include <climits>
#include <limits>
#include <stdlib.h>

namespace lightana {

  static const double kINVALID_TIME = std::numeric_limits<double>::max();
  static const unsigned int kINVALID_CHANNEL = std::numeric_limits<unsigned int>::max();

  struct LiteOpHit_t {
    size_t channel;
    double peak_time;
    double pe;
    LiteOpHit_t() : channel(kINVALID_CHANNEL),
                    peak_time(kINVALID_TIME),
              pe(0)
    {}
  };

  typedef std::vector<lightana::LiteOpHit_t> LiteOpHitArray_t;

  struct LiteOpFlash_t {
    std::vector<double> channel_pe;
    std::vector<unsigned int> asshit_idx;
    double time;
    double time_err;
    
    LiteOpFlash_t() : channel_pe(),
                     time(kINVALID_TIME)
    {}

    LiteOpFlash_t(double flash_time,
                  double flash_time_err,
                  std::vector<double>&& pe_array,
                  std::vector<unsigned int>&& hit_idx)
      : channel_pe(std::move(pe_array))
      , asshit_idx(std::move(hit_idx))
    { time = flash_time; time_err = flash_time_err; }
    
    void Register(size_t channel, double pe)
    {
      if(channel >= channel_pe.size()) channel_pe.resize(channel+1,0);
      channel_pe[channel] = pe;
    }
  };

  typedef std::vector<lightana::LiteOpFlash_t> LiteOpFlashArray_t;
}
#endif
