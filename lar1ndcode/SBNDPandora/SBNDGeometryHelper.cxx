/**
 *  @file  lar1ndcode/SBNDPandora/SBNDGeometryHelper.cxx
 *
 *  @brief helper function for SBND geometry
 *
 */

#include "lar1ndcode/SBNDPandora/SBNDGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

SBNDGeometryHelper::SBNDVolume SBNDGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Left drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return SBNDGeometryHelper::kLeftVolume;
    }

    // Right drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return SBNDGeometryHelper::kRightVolume;
    }

    throw cet::exception("LArPandora") << " SBNDGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
