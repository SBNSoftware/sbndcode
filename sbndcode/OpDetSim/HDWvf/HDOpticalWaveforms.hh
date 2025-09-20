///////////////////////////////////////////////////////////////////////
///
/// Interface class for HDOpticalWaveforms tool
///
//////////////////////////IncludeHDOpticalWaveforms//////////////////////////////////////////////

#ifndef SBND_HDOpticalWaveform_H
#define SBND_HDOpticalWaveform_H


namespace opdet {
  class HDOpticalWaveform;
}

//Base class
class opdet::HDOpticalWaveform {
public:
  //Constructor
  virtual ~HDOpticalWaveform() noexcept = default;

  //Returns fluctuated factor for SPR
  virtual size_t TimeBinShift(double TimeBin_HD) = 0;

  //Returns fluctuated factor for SPR
  virtual int produceSER_HD(std::vector<std::vector<double>> &SER_HD, std::vector<double>& SER) = 0;
};

#endif
