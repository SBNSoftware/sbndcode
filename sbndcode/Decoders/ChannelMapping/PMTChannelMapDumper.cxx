/**
 * @file   sbndcode/Decoders/ChannelMapping/PMTChannelMapDumper.cxx
 * @brief  Utility dumping the content of PMT channel mapping on screen.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @author Vu Chi Lan Nguyen (vclnguyen1@sheffield.ac.uk) adapted for SBND
 * 
 * This utility can be run with any configuration file including a configuration
 * for `SBNDChannel` service.
 * 
 * It is using _art_ facilities for tool loading, but it does not run in _art_
 * environment. So it may break without warning and without solution.
 * 
 */


// SBND libraries
#include "sbndcode/Decoders/ChannelMapping/SBNDChannelMapProvider.h"

// LArSoft and framework libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/TestUtils/unit_test_base.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <iostream>
#include <algorithm>
#include <numeric> // std::iota()
#include <array>


// -----------------------------------------------------------------------------
template <std::size_t KeyNo = 0U>
struct SortByElement {
  
  template <typename TupleA, typename TupleB>
  bool operator() (TupleA const& A, TupleB const& B) const
    { return less(A, B); }
  
  template <typename Tuple>
  static decltype(auto) key(Tuple const& t) { return std::get<KeyNo>(t); }
  
  template <typename TupleA, typename TupleB>
  static bool less(TupleA const& A, TupleB const& B)
    { using std::less; return less{}(key(A), key(B)); }
  
}; // SortByElement<>


// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
  
  using Environment
    = testing::TesterEnvironment<testing::BasicEnvironmentConfiguration>;
  
  testing::BasicEnvironmentConfiguration config("PMTchannelMappingDumper");

  //
  // parameter parsing
  //
  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc)
    config.SetConfigurationPath(argv[iParam]);
  else {
    std::cerr << "FHiCL configuration file path required as first argument!"
      << std::endl;
    return 1;
  }

  Environment const Env { config };
  
  sbndDB::SBNDChannelMapProvider const channelMapping
    { Env.ServiceParameters("SBNDChannelMap") };
  
  // hard-coded list of fragment ID; don't like it?
  // ask for an extension of the channel mapping features.
  std::array<unsigned int, 24U> FragmentIDs;
  std::iota(FragmentIDs.begin(), FragmentIDs.end(), 0x2000);
  
  mf::LogVerbatim("PMTchannelMappingDumper") << "Fragment IDs:";
  for (auto const [ iFragment, fragmentID]: util::enumerate(FragmentIDs)) {
    unsigned int const effFragmentID = fragmentID & 0xFF;
    
    if (!channelMapping.hasPMTDigitizerID(effFragmentID)) {
      mf::LogVerbatim("PMTchannelMappingDumper")
        << "[" << iFragment << "] " << std::hex << fragmentID << std::dec
        << " not found in the database (as "
        << std::hex << effFragmentID << std::dec << ")";
      continue;
    }
    
    sbndDB::DigitizerChannelChannelIDPairVec digitizerChannels
      = channelMapping.getChannelIDPairVec(effFragmentID);
    
    std::sort
      (digitizerChannels.begin(), digitizerChannels.end(), SortByElement<1U>{});
    
    
    mf::LogVerbatim log("PMTchannelMappingDumper");
    log
      << "[" << iFragment << "] " << std::hex << fragmentID << std::dec
      << " includes " << digitizerChannels.size()
      << " LArSoft channels between " << std::get<1U>( digitizerChannels.front() )
      << " and " << std::get<1U>( digitizerChannels.back() )
      << " [board channel index in brackets]:";
    constexpr unsigned int Cols = 8U;
    unsigned int n = 0;
    for(auto const [ digitizerChannel, channelID, laserChannel ]: digitizerChannels) {
      if (n-- == 0) { log << "\n     "; n = Cols - 1U; }
      log << " "  << std::setw(3) << channelID
        << " [" << std::setw(3) << digitizerChannel << "]";
    } // for channel
    
  } // for fragment
  
  return 0;
} // main()


