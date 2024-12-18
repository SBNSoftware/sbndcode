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
#include "larrecodnn/CVN/func/PixelMap.h"

namespace lcvn {

  /// Wrapper for caffe::Net which handles construction and prediction
  class SBNDITFNetHandler {
  public:
    virtual ~SBNDITFNetHandler() noexcept = default;
    /// Return prediction arrays for PixelMap
    virtual std::vector<std::vector<float>> Predict(const PixelMap& pm) const = 0;
  };

}

#endif // LCVN_SBNDTFNETHANDLER_H
