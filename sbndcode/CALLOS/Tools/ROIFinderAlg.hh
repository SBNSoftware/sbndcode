///////////////////////////////////////////////////////////////////////
///
/// Interface class for ROIFindetAlg tool
///
////////////////////////////////////////////////////////////////////////

#ifndef SBND_ROIFINDERALG_H
#define SBND_ROIFINDERALG_H

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include <cmath>
#include <vector>

namespace callos {
  class ROIFINDERALG;
}

//Base class
class callos::ROIFINDERALG {
public:
  //Constructor
  virtual ~ROIFINDERALG() noexcept = default;

  // Required functions.
  virtual bool ProcessWaveform(std::vector<float> const& wvf ,std::vector<float> & ROI, std::vector<float> & charges)=0;
};

#endif
