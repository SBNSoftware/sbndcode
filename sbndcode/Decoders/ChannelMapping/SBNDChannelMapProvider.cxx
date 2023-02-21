////////////////////////////////////////////////////////////////////////
/// \file   sbndcode/Decoders/ChannelMapping/SBNDChannelMapProvider.cxx
/// \author T. Usher (factorized by Gianluca Petrillo, petrillo@slac.stanford.edu)
//  \author L. Nguyen (vclnguyen1@sheffield.ac.uk) adapted for SBND
////////////////////////////////////////////////////////////////////////

// library header
#include "sbndcode/Decoders/ChannelMapping/SBNDChannelMapProvider.h"

#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

#include <string>
#include <iostream>
#include <cassert>

namespace sbndDB
{


//----------------------------------------------------------------------
// Constructor.
SBNDChannelMapProvider::SBNDChannelMapProvider(const fhicl::ParameterSet& pset) {

    mf::LogInfo("SBNDChannelMapProvider") << "Building the channel mapping" ;

    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet&channelMappingParams = pset.get<fhicl::ParameterSet>("ChannelMappingTool");

    // Get instance of the mapping tool (allowing switch between database instances)
    fChannelMappingTool = art::make_tool<IChannelMapping>(channelMappingParams);

    cet::cpu_timer theClockReadoutIDs;

    theClockReadoutIDs.start();


    // Do the channel mapping initialization
    if (fChannelMappingTool->BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
      {
	throw cet::exception("SBNDChannelMapProvider") << "Cannot recover the Fragment ID channel map from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size() << " Fragment IDs";

// 	for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::hex << pair.first << ", # pairs: " 
//								   << std::dec << pair.second.size() << std::endl;
	 
         for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::dec << pair.first << ", # pairs: " 
								   << std::dec << pair.second.size() << std::endl;
      }
    
    double readoutIDsTime = theClockReadoutIDs.accumulated_real_time();


    mf::LogInfo("SBNDChannelMapProvider") << "==> Readout IDs time: " << readoutIDsTime << std::endl;
    
    return;
}

bool SBNDChannelMapProvider::hasPMTDigitizerID(const unsigned int fragmentID)   const
{
    return findPMTfragmentEntry(fragmentID) != nullptr;
}

unsigned int SBNDChannelMapProvider::nPMTfragmentIDs() const {
  return fFragmentToDigitizerMap.size();
}

const DigitizerChannelChannelIDPairVec& SBNDChannelMapProvider::getChannelIDPairVec(const unsigned int fragmentID) const
{
    DigitizerChannelChannelIDPairVec const* digitizerPair = findPMTfragmentEntry(fragmentID);
    
    if (digitizerPair) return *digitizerPair;
    throw cet::exception("SBNDChannelMapProvider") << "Fragment ID " << fragmentID << " not found in lookup map when looking for PMT channel info \n";
    
}

auto SBNDChannelMapProvider::findPMTfragmentEntry(unsigned int fragmentID) const
  -> DigitizerChannelChannelIDPairVec const*
{
  auto it = fFragmentToDigitizerMap.find(PMTfragmentIDtoDBkey(fragmentID));
  return (it == fFragmentToDigitizerMap.end())? nullptr: &(it->second);
}


constexpr unsigned int SBNDChannelMapProvider::PMTfragmentIDtoDBkey
  (unsigned int fragmentID)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert(((fragmentID & ~0xFF) == 0x00) || ((fragmentID & ~0xFF) == 0x20));
  
  return fragmentID & 0xFF;
  
} // SBNDChannelMapProvider::PMTfragmentIDtoDBkey()


constexpr unsigned int SBNDChannelMapProvider::DBkeyToPMTfragmentID
  (unsigned int DBkey)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert((DBkey & 0xFF) < 24);
  
  return (DBkey & 0xFF) | 0x2000;
  
} // SBNDChannelMapProvider::PMTfragmentIDtoDBkey()


} // end namespace
