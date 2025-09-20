////////////////////////////////////////////////////////////////////////
////
////  FlasFinderManager source
////
//////////////////////////////////////////////////////////////////////////

#ifndef FLASHFINDERMANAGER_CXX
#define FLASHFINDERMANAGER_CXX

#include "FlashFinderManager.h"
#include <sstream>

namespace lightana{
  
  FlashFinderManager::FlashFinderManager() :
    _flash_algo(nullptr)  
  {}

  FlashFinderManager::~FlashFinderManager() {}
  
  void FlashFinderManager::SetFlashAlgo (FlashAlgoBase* algo)
  
  {
    if(!algo) {
      std::cerr << "\n\t Not a valid algo\n";
      throw std::exception();
    }
    _flash_algo = algo;
  }

  LiteOpFlashArray_t FlashFinderManager::RecoFlash(const LiteOpHitArray_t& ophits) const
                                             
  {
    if(!_flash_algo) {
      std::cerr << "\n\t No flashing to run!\n";
      throw std::exception();
    }
    return _flash_algo->RecoFlash(ophits);
  }

}
#endif
