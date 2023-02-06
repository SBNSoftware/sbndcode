/**
 *  @file   IChannelMapping.h
 *
 *  @brief  This provides an art tool interface definition for tools handle the channel mapping
 *          The idea is to be able to switch between postgres and sqlite implementations
 *
 *  @author usher@slac.stanford.edu
 *  @author vclnguyen1@sheffield.ac.uk adapted for SBND
 *
 */
#ifndef IChannelMapping_h
#define IChannelMapping_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace sbndDB
{
  /**
   *  @brief  IChannelMapping interface class definiton
   */
  class IChannelMapping
  {
  public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IChannelMapping() noexcept = default;
    
    /**
     *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
     *         and the related crate and readout information. 
     *         Then define the function interface to fill these data structures 
     */
    using DigitizerChannelChannelIDPair    = std::tuple<size_t,size_t,size_t>; // std::tuple<DigitizerChannel, ChannelID, LaserChannel>
    using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;
    using FragmentToDigitizerChannelMap    = std::map<size_t, DigitizerChannelChannelIDPairVec>;

    virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const = 0;

};

} // namespace sbndDB 
#endif
