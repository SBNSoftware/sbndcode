///////////////////////////////////////////////////////////////////////
///
/// Interface class for ROIFindetAlg tool
///
////////////////////////////////////////////////////////////////////////

#ifndef SBND_ROIFINDERALG_H
#define SBND_ROIFINDERALG_H

namespace opdet {
  class ROIFINDERALG;
}

//Base class
class opdet::ROIFINDERALG {
public:
  //Constructor
  virtual ~ROIFINDERALG() noexcept = default;

  // Required functions.
  virtual std::vector<std::vector<std::float>> FindROIS(std::vector<short> const& wfHandle) = 0;
};

#endif
