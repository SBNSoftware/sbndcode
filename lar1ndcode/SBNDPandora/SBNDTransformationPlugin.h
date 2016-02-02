/**
 *  @file   lar1ndcode/SBNDPandora/SBNDTransformationPlugin.h
 *
 *  @brief  Header file for the SBND transformation plugin class.
 *
 *  $Log: $
 */
#ifndef SBND_TRANSFORMATION_PLUGIN_H
#define SBND_TRANSFORMATION_PLUGIN_H 1

#include "LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  SBNDTransformationPlugin class
 */
class SBNDTransformationPlugin : public lar_content::LArRotationalTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    SBNDTransformationPlugin(const bool isForward);
};

} // namespace lar_pandora

#endif // #ifndef SBND_TRANSFORMATION_PLUGIN_H
