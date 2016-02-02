/**
 *  @file  lar1ndcode/SBNDPandora/SBNDGeometryHelper.h
 *
 *  @brief helper function for SBND geometry
 *
 */
#ifndef SBND_PANDORA_HELPER_H
#define SBND_PANDORA_HELPER_H

namespace lar_pandora 
{

class SBNDGeometryHelper 
{
public:

    enum SBNDVolume
    {
        kLeftVolume = 0,
        kRightVolume = 1,
        kUnknownVolume = 2
    };

    /**
     *  @brief Assign a drift volume ID based on cryostate and TPC
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc 
     */
     static SBNDGeometryHelper::SBNDVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_SBND_PANDORA_HELPER_H
