//really bad flash framework
//by vic

#ifndef FLASHFINDERMANAGER_H
#define FLASHFINDERMANAGER_H

#include "FlashAlgoBase.h"
#include "FlashAlgoFactory.h"

namespace lightana
{

  class FlashFinderManager {

  public:
    FlashFinderManager();
    
    ~FlashFinderManager();

    void SetFlashAlgo (FlashAlgoBase* algo);

    LiteOpFlashArray_t RecoFlash(const LiteOpHitArray_t& ophits) const;
    
  private:
    
    FlashAlgoBase* _flash_algo;
  };
  
}
#endif

