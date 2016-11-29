#include "art/Persistency/Common/Wrapper.h"
#include "sbndcode/CRT/CRTData.hh"

namespace {
  struct dictionary {
    art::Wrapper<crt::CRTData> w;
    std::vector<crt::CRTData> v;
    art::Wrapper<std::vector<crt::CRTData> > wv;
  };
}


