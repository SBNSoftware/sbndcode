/**
 *  @file   lar1ndcode/SBNDPandora/SBNDPandora_module.cc
 *
 *  @brief  Producer module for SBND 4APA detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "LArPandoraInterface/LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  SBNDPandora class
 */
class SBNDPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    SBNDPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~SBNDPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useLeftVolume;      ///<
    bool            m_useRightVolume;       ///<
};

DEFINE_ART_MODULE(SBNDPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArContent.h"

// Local includes (LArPandora)
#include "lar1ndcode/SBNDPandora/SBNDPseudoLayerPlugin.h"
#include "lar1ndcode/SBNDPandora/SBNDTransformationPlugin.h"
#include "lar1ndcode/SBNDPandora/SBNDGeometryHelper.h"

namespace lar_pandora {

SBNDPandora::SBNDPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useLeftVolume = pset.get<bool>("UseLeftVolume",true);
    m_useRightVolume = pset.get<bool>("UseRightVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

SBNDPandora::~SBNDPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SBNDPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** SBNDPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("lar1nd"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " SBNDPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new SBNDPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new SBNDTransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int SBNDPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const SBNDGeometryHelper::SBNDVolume volumeID(SBNDGeometryHelper::GetVolumeID(cstat, tpc));

    if (SBNDGeometryHelper::kLeftVolume == volumeID) 
    {
        if (m_useLeftVolume) 
            return 0;
    }

    if (SBNDGeometryHelper::kRightVolume == volumeID) 
    {
        if (m_useRightVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
