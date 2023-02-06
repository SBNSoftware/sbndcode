///////////////////////////////////////////////////////////////////////
///
/// \file   ISBNDChannelMap
///
/// \brief  Interface class for hardware/software channel mapping 
///         for SBND
///
/// \author L. Nguyen (vclnguyen1@sheffield.ac.uk) adapted for SBN
///
////////////////////////////////////////////////////////////////////////

#ifndef ISBNDChannelMap_H
#define ISBNDChannelMap_H

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include <vector>
#include <string>

namespace sbndDB {

using DigitizerChannelChannelIDPair    = std::tuple<size_t,size_t, size_t>;
using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;

class ISBNDChannelMap //: private lar::EnsureOnlyOneSchedule
{
public:
    virtual ~ISBNDChannelMap() noexcept = default;

    // Section for recovering PMT information
    virtual bool                                    hasPMTDigitizerID(const unsigned int)   const = 0;
    /// Returns the number of PMT fragment IDs known to the service.
    virtual unsigned int                            nPMTfragmentIDs() const = 0;
    virtual const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const = 0;

};

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::ISBNDChannelMap, SHARED)

#endif
