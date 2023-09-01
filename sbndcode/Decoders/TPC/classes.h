#include "canvas/Persistency/Common/Wrapper.h"
#include <vector>
#include "sbndcode/Decoders/TPC/TPCDecodeAna.h"


namespace {
  struct dictionary {
    tpcAnalysis::TPCDecodeAna h;
    std::vector<tpcAnalysis::TPCDecodeAna> h_v;
    art::Wrapper<tpcAnalysis::TPCDecodeAna> h_w;
    art::Wrapper<std::vector<tpcAnalysis::TPCDecodeAna>> h_v_w;
  };
}


