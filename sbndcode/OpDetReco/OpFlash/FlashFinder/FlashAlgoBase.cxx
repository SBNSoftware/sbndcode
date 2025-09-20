#ifndef FLASHALGOBASE_CXX
#define FLASHALGOBASE_CXX

#include "FlashAlgoBase.h"

namespace lightana{

  FlashAlgoBase::FlashAlgoBase(const std::string name)
  {
    _name = name;
    Reset();
  }

  FlashAlgoBase::~FlashAlgoBase() {}

  void FlashAlgoBase::Reset() {}
  
}

#endif
