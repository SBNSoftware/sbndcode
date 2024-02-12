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


namespace callos {
  class ROIFINDERALG;
}

//Base class
class callos::ROIFINDERALG {
public:
  //Constructor
  virtual ~ROIFINDERALG() noexcept = default;

  // Required functions.
  // virtual std::vector<std::vector<std::float>> FindROIS(std::vector<short> const& wfHandle);
};

#endif
