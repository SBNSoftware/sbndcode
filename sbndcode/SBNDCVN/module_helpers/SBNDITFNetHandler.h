/////////////////////////////////////////////////////////////////////////////////////////////////
/// \file    SBNDTFNetHandler.h
/// \brief   SBNDTFNetHandler for CVN
/// \author  Varuna Meddage by looking at the format of larrecodnn/CVN/interfaces/ITFNetHandler.h
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LCVN_SBNDTFNETHANDLER_H
#define LCVN_SBNDTFNETHANDLER_H

#include <memory>
#include <vector>

#include "larrecodnn/CVN/func/InteractionType.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"

namespace lcvn {

  /// Wrapper for caffe::Net which handles construction and prediction
  class SBNDITFNetHandler {
  public:
    virtual ~SBNDITFNetHandler() noexcept = default;
    /// Return prediction arrays for SBNDPixelMap
    virtual std::vector<std::vector<float>> Predict(const SBNDPixelMap& pm) const = 0;
  };

}

#endif // LCVN_SBNDTFNETHANDLER_H
