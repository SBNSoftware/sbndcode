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
#include "sbndcode/CALLOS/AverageWaveform.h"

namespace callos {
  class ROIFINDERALG;
}

//Base class
class callos::ROIFINDERALG {
public:
  //Constructor
  virtual ~ROIFINDERALG() noexcept = default;

  // Required functions.
  virtual bool ProcessWaveform(std::vector<float> const& wvf ,std::vector<SimpleROI> & ROI)=0;
};

#endif
