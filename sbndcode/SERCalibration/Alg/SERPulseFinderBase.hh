#ifndef SBND_SERPULSEFINDERBASE_H
#define SBND_SERPULSEFINDERBASE_H

#include "lardataobj/RawData/OpDetWaveform.h"
#include "TH1D.h"

namespace opdet {
  class SERPulseFinderBase;
}

class opdet::SERPulseFinderBase{
public:

    virtual ~SERPulseFinderBase() noexcept = default;

    virtual void RunSERCalibration(std::vector<raw::OpDetWaveform> const&  , std::vector<TH1D>& ) = 0;
    virtual std::vector<std::vector<int>> GetPeakAmplitudeVector()=0;
    virtual std::vector<int> GetNumberOfPeaksVector()=0;
};

#endif