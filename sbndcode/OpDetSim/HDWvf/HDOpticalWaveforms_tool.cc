////////////////////////////////////////////////////////////////////////
// Specific class tool for HDOpticalWaveform
// File: HDOpticalWaveforms_tool.hh
// Base class:        HDOpticalWaveform.hh
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <vector>
#include <cmath>

#include "sbndcode/OpDetSim/HDWvf/HDOpticalWaveforms.hh"


namespace opdet {
  class HDOpticalWaveforms;
}


class opdet::HDOpticalWaveforms : opdet::HDOpticalWaveform {
public:

  //Configuration parameters
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<double> Resolution {
      Name("ResolutionIncrease"),
      Comment("Increase in resolution of the optical wvf, i.e: real sampling is 12.5ns/tick but HD waveform is binned in 1.25ns/tick-> ResolutionIncrease = 10 ")
    };

  };

  explicit HDOpticalWaveforms(art::ToolConfigTable<Config> const& config);

  int produceSER_HD(std::vector<std::vector<double>> &SER_HD, std::vector<double>& SER) override;
  
  size_t TimeBinShift(double TimeBin_HD) override;
  

private:
  //Configuration parameters
  double fResolution;
};


opdet::HDOpticalWaveforms::HDOpticalWaveforms(art::ToolConfigTable<Config> const& config)
  : fResolution { config().Resolution() }
  // , fGain { config().gain() }
{
}

int opdet::HDOpticalWaveforms::produceSER_HD(std::vector<std::vector<double>> &SER_HD, std::vector<double>& SER)
{
    int N=fResolution;
    SER_HD.resize(N);

    for(int i=1; i<N;i++)SER_HD[i%N].push_back(0);

    for(int i=0; i < (int)SER.size(); i++)
    {
      if (i%N == 0)SER_HD[i%N].push_back(SER.at(i)); //prepare xRF SPEs from higher sampling rate estimation
      else SER_HD[N-i%N].push_back(SER.at(i));
    }

  return 0;
}

size_t opdet::HDOpticalWaveforms::TimeBinShift(double TimeBin_HD){
  
  // Get only the decimal part
  float whole, fractional;

  fractional = std::modf(TimeBin_HD, &whole);
  
  size_t shift=std::floor(fractional*fResolution); // 1/fResolution it's the HD tick width: Res=10-> HDtick =1/10
  
  shift=shift% ((int)fResolution);

  return shift;
}


DEFINE_ART_CLASS_TOOL(opdet::HDOpticalWaveforms)
