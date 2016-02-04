/**
 *  @file   lar1ndcode/SBNDPandora/SBNDTransformationPlugin.cxx
 *
 *  @brief  Implementation of the SBND transformation plugin class.
 *
 *  $Log: $
 */

#include "lar1ndcode/SBNDPandora/SBNDTransformationPlugin.h"

#include <cmath>

namespace lar_pandora
{

SBNDTransformationPlugin::SBNDTransformationPlugin(const bool isForward) :
    lar_content::LArRotationalTransformationPlugin( (isForward ? M_PI / 3.f : -M_PI / 3.f),
        (isForward ? M_PI / 3.f : -M_PI / 3.f), 1.f)
{
}

} // namespace lar_pandora
