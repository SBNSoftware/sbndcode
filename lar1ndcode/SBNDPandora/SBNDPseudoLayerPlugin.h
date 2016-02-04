/**
 *  @file   lar1ndcode/SBNDPandora/SBNDPseudoLayerPlugin.h
 *
 *  @brief  Header file for the SBND pseudo layer plugin class.
 *
 *  $Log: $
 */
#ifndef SBND_PSEUDO_LAYER_PLUGIN_H
#define SBND_PSEUDO_LAYER_PLUGIN_H 1

#include "LArPlugins/LArPseudoLayerPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  SBNDPseudoLayerPlugin class
 */
class SBNDPseudoLayerPlugin : public lar_content::LArPseudoLayerPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    SBNDPseudoLayerPlugin();
};

} // namespace lar

#endif // #ifndef SBND_PSEUDO_LAYER_PLUGIN_H
