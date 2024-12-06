///////////////////////////////////////////////////////////////////////
///
/// Interface class for OpDeconvolution tool
///
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDECONVOLUTIONALG_H
#define SBND_OPDECONVOLUTIONALG_H

#include "lardataobj/RawData/OpDetWaveform.h"

namespace opdet {
  class OpDeconvolutionAlg;
}

//Base class
class opdet::OpDeconvolutionAlg {
public:
  //Constructor
  virtual ~OpDeconvolutionAlg() noexcept = default;

  // Required functions.
  virtual std::vector<raw::OpDetWaveform> RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfHandle) = 0;
};

#endif
