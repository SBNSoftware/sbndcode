#include "canvas/Persistency/Common/Wrapper.h"
#include "sbndcode/CRT/CRTData.hh"
#include <vector>

namespace {
  struct dictionary {
    art::Wrapper<sbnd::crt::CRTData> w;
    std::vector<sbnd::crt::CRTData> v;
    art::Wrapper<std::vector<sbnd::crt::CRTData> > wv;
  };
}


