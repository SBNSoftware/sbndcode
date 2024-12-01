#include "sbndcode/SBNDCVN/module_helpers/SBNDICVNMapper.cxx"
#include "sbndcode/SBNDCVN/module_helpers/SBNDICVNMapper.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMapProducer.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"

namespace lcvn {

  typedef SBNDICVNMapper<lcvn::SBNDPixelMapHitProducer, recob::Hit> SBNDCVNMapper;
  template class SBNDICVNMapper<lcvn::SBNDPixelMapHitProducer, recob::Hit>;

DEFINE_ART_MODULE(lcvn::SBNDCVNMapper)
}
