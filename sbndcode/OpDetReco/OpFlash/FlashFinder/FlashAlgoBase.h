#ifndef FLASHALGOBASE_H
#define FLASHALGOBASE_H

#include <iostream>
#include "FlashFinderTypes.h"
#include "FlashFinderFMWKInterface.h"
#include <vector>

namespace lightana
{
  class FlashAlgoBase {

  public:
    FlashAlgoBase(const std::string name);
    const std::string& Name() const { return _name; }
    virtual ~FlashAlgoBase();
    virtual void Configure(const Config_t &p) = 0;
    virtual LiteOpFlashArray_t RecoFlash(const LiteOpHitArray_t ophits) = 0;
    virtual void Reset();

  private:
    std::string _name;
    
  };
}

#endif
/** @} */ // end of doxygen group
