///////////////////////////////////////////////////////////////////////
///
/// \file   SBNDChannelMap
///
/// \brief  Interface class for hardware/software channel mapping 
///         for SBNDS
///
/// \author Afroditi Papadopoulou (apapadopoulou@anl.gov)
///
////////////////////////////////////////////////////////////////////////

#ifndef SBNDChannelMap_H
#define SBNDChannelMap_H

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include <vector>
#include <string>

namespace sbndDB {

using ReadoutIDVec                     = std::vector<unsigned int>;
using ChannelPlanePair                 = std::pair<unsigned int, unsigned int>;
using ChannelPlanePairVec              = std::vector<ChannelPlanePair>;

using DigitizerChannelChannelIDPair    = std::tuple<size_t,size_t, size_t>;
using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;

using ChannelPlanePair                 = std::pair<unsigned int, unsigned int>;
using ChannelPlanePairVec              = std::vector<ChannelPlanePair>;
using SlotChannelVecPair               = std::pair<unsigned int, ChannelPlanePairVec>;
using TPCReadoutBoardToChannelMap      = std::map<unsigned int, SlotChannelVecPair>;

class SBNDChannelMap //: private lar::EnsureOnlyOneSchedule
{
public:
    virtual ~SBNDChannelMap() noexcept = default;

    // Section to access fragment to board mapping
    virtual bool                                    hasFragmentID(const unsigned int)       const = 0;
    /// Returns the number of TPC fragment IDs known to the service.
    virtual unsigned int                            nTPCfragmentIDs() const = 0;
    virtual const std::string&                      getCrateName(const unsigned int)        const = 0;
    virtual const ReadoutIDVec&                     getReadoutBoardVec(const unsigned int)  const = 0;
    virtual const TPCReadoutBoardToChannelMap&      getReadoutBoardToChannelMap()           const = 0;

    // Section to access channel information for a given board
    virtual bool                                    hasBoardID(const unsigned int)          const = 0;
    /// Returns the number of TPC board IDs known to the service.
    virtual unsigned int                            nTPCboardIDs() const = 0;
    virtual unsigned int                            getBoardSlot(const unsigned int)        const = 0;
    virtual const ChannelPlanePairVec&              getChannelPlanePair(const unsigned int) const = 0;

    // Section for recovering PMT information
    virtual bool                                    hasPMTDigitizerID(const unsigned int)   const = 0;
    /// Returns the number of PMT fragment IDs known to the service.
    virtual unsigned int                            nPMTfragmentIDs() const = 0;
    virtual const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const = 0;

    virtual unsigned int                            getSimMacAddress(const unsigned int)    const = 0;    
    virtual unsigned int                         gettopSimMacAddress(const unsigned int)    const = 0;    

    virtual std::pair<double, double>          getSideCRTCalibrationMap(int mac5, int chan) const = 0;
};

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::SBNDChannelMap, SHARED)

#endif
