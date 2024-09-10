///////////////////////////////////////////////////////////////////////
/// File: FlashT0Base.h
///
/// Interfacce class for a tool to calculate the recob::OpFlash t0
/// from the associated recob::OpHits
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_FLASHT0BASE_H
#define SBND_FLASHT0BASE_H

#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderTypes.h"

namespace lightana
{
  class FlashT0Base{

  public:
    // Default destructor
    virtual ~FlashT0Base() noexcept = default;

    // Method to calculate the OpFlash t0
    virtual double GetFlashT0(double flash_peaktime, LiteOpHitArray_t ophit_list) = 0;

  private:

  };
}

#endif
