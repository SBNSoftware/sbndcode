#ifndef CRTData_hh_
#define CRTData_hh_

#include <stdint.h> //uint32_t
#include <vector>
#include <utility>

namespace sbnd {
namespace crt {

  class CRTData {

    uint32_t fChannel;
    uint32_t fT0;
    uint32_t fT1;
    uint32_t fADC;
    int      fTrackID; //FIXME lazy truth matching, to be removed

   public:

    CRTData();
    CRTData(uint32_t channel, uint32_t t0, uint32_t t1, uint32_t adc, int trackID);

    virtual ~CRTData();
    uint32_t Channel() const;
    uint32_t T0() const;
    uint32_t T1() const;
    uint32_t ADC() const;
    int TrackID() const;

  };

} // namespace crt
} // namespace sbnd

#endif
