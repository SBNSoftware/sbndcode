#include "CRTGeoAlg.h"

namespace sbnd{

  CRTGeoAlg::CRTGeoAlg(fhicl::ParameterSet const &p) :
    CRTGeoAlg(p, lar::providerFrom<geo::Geometry>(), 
              ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr())
  {}

  CRTGeoAlg::CRTGeoAlg(fhicl::ParameterSet const &p, geo::GeometryCore const *geometry, 
                       geo::AuxDetGeometryCore const *auxdet_geometry)
    : fCableLengthCorrectionsVector(p.get<std::vector<std::pair<unsigned, double>>>("CableLengthCorrections", std::vector<std::pair<unsigned, double>>()))
    , fSiPMPedestalsVector(p.get<std::vector<std::pair<unsigned, double>>>("SiPMPedestals", std::vector<std::pair<unsigned, double>>()))
    , fChannelInversionVector(p.get<std::vector<std::pair<unsigned, bool>>>("InvertedChannelOrder", std::vector<std::pair<unsigned, bool>>()))
  {
    fGeometryService = geometry;
    fAuxDetGeoCore   = auxdet_geometry;
    TGeoManager* manager = fGeometryService->ROOTGeoManager();

    fCableLengthCorrections = std::map<unsigned, double>(fCableLengthCorrectionsVector.begin(), 
                                                         fCableLengthCorrectionsVector.end());
    fSiPMPedestals          = std::map<unsigned, double>(fSiPMPedestalsVector.begin(), 
                                                         fSiPMPedestalsVector.end());
    fChannelInversion       = std::map<unsigned, bool>(fChannelInversionVector.begin(), 
                                                       fChannelInversionVector.end());

    // Record used objects
    std::vector<std::string> usedTaggers;
    std::vector<std::string> usedModules;
    std::vector<std::string> usedStrips;

    // Loop through aux dets
    const std::vector<geo::AuxDetGeo> &auxDets = fAuxDetGeoCore->AuxDetGeoVec();

    for(unsigned ad_i = 0; ad_i < auxDets.size(); ++ad_i)
      {
        const geo::AuxDetGeo &auxDet = auxDets[ad_i];

        // Get the geometry object for the auxDet
        std::set<std::string> volNames = {auxDet.TotalVolume()->GetName()};
        std::vector<std::vector<TGeoNode const*>> paths = fGeometryService->FindAllVolumePaths(volNames);

        // Build up a path that ROOT understands
        std::string path = "";
        for (size_t inode = 0; inode < paths.at(0).size(); inode++){
          path += paths.at(0).at(inode)->GetName();
          if(inode < paths.at(0).size() - 1){
            path += "/";
          }
        }

        // Loop through strips
        for(unsigned ads_i = 0; ads_i < auxDet.NSensitiveVolume(); ads_i++)
          {
            const geo::AuxDetSensitiveGeo &auxDetSensitive = auxDet.SensitiveVolume(ads_i);

            // Access the path
            manager->cd(path.c_str());

            TGeoNode* nodeArray = manager->GetCurrentNode();
            TGeoNode* nodeStrip = nodeArray->GetDaughter(ads_i);
            TGeoNode* nodeModule = manager->GetMother(1);
            TGeoNode* nodeTagger = manager->GetMother(2);
            TGeoNode* nodeDet = manager->GetMother(3);

            // Fill the tagger information
            const std::string taggerName = nodeTagger->GetName();
            if(std::find(usedTaggers.begin(), usedTaggers.end(), taggerName) == usedTaggers.end())
              {
                usedTaggers.push_back(taggerName);
                CRTTaggerGeo tagger  = CRTTaggerGeo(nodeTagger, nodeDet);
                fTaggers.insert(std::pair<std::string, CRTTaggerGeo>(taggerName, tagger));
              }

            // Fill the module information
            const std::string moduleName = nodeModule->GetName();
	    const bool invert = fChannelInversion.size() ? fChannelInversion.at(ad_i) : false;
            if(std::find(usedModules.begin(), usedModules.end(), moduleName) == usedModules.end())
              {
                const uint32_t cableDelayCorrection = fCableLengthCorrections.size() ? 
                  fCableLengthCorrections.at(ad_i) : 0;
                usedModules.push_back(moduleName);
                CRTModuleGeo module  = CRTModuleGeo(nodeModule, auxDet, ad_i, taggerName,
                                                    cableDelayCorrection, invert);
                fModules.insert(std::pair<std::string, CRTModuleGeo>(moduleName, module));
              }

            // Fill the strip information
            const std::string stripName = nodeStrip->GetName();
	    // Some modules need their channel numbers counted in reverse as they're inverted relative to the geometry
            const uint32_t channel0 = invert ? 32 * ad_i + (31 - 2 * ads_i) : 32 * ad_i + 2 * ads_i;
            const uint32_t channel1 = invert ? 32 * ad_i + (31 - 2 * ads_i -1) : 32 * ad_i + 2 * ads_i + 1;
            if(std::find(usedStrips.begin(), usedStrips.end(), stripName) == usedStrips.end())
              {
                usedStrips.push_back(stripName);
                CRTStripGeo strip  = CRTStripGeo(nodeStrip, auxDetSensitive, ads_i, moduleName,
                                                 channel0, channel1);
                fStrips.insert(std::pair<std::string, CRTStripGeo>(stripName, strip));
              }

            double halfWidth  = auxDetSensitive.HalfWidth1();
            double halfHeight = auxDetSensitive.HalfHeight();

            // Calcualte SiPM locations
            // SiPM0 is on the left in local coordinates
            const double sipm0Y = -halfHeight;
            const double sipm1Y = halfHeight;
            const double sipmX  = fModules.at(moduleName).top ? halfWidth : -halfWidth;

            // Find world coordinates
            const double sipm0XYZ[3] = {sipmX, sipm0Y, 0};
            double sipm0XYZWorld[3];
            auxDetSensitive.LocalToWorld(sipm0XYZ, sipm0XYZWorld);

            const double sipm1XYZ[3] = {sipmX, sipm1Y, 0};
            double sipm1XYZWorld[3];
            auxDetSensitive.LocalToWorld(sipm1XYZ, sipm1XYZWorld);

            const uint32_t pedestal0 = fSiPMPedestals.size() ? fSiPMPedestals.at(channel0) : 0;
            const uint32_t pedestal1 = fSiPMPedestals.size() ? fSiPMPedestals.at(channel1) : 0;

            // Fill SiPM information
            CRTSiPMGeo sipm0 = CRTSiPMGeo(stripName, channel0, sipm0XYZWorld, pedestal0);
            CRTSiPMGeo sipm1 = CRTSiPMGeo(stripName, channel1, sipm1XYZWorld, pedestal1);
            fSiPMs.insert(std::pair<uint16_t, CRTSiPMGeo>(channel0, sipm0));
            fSiPMs.insert(std::pair<uint16_t, CRTSiPMGeo>(channel1, sipm1));
          }
      }
  }

  CRTGeoAlg::~CRTGeoAlg() {}

  std::vector<double> CRTGeoAlg::CRTLimits() const {
    std::vector<double> limits;

    std::vector<double> minXs;
    std::vector<double> minYs;
    std::vector<double> minZs;
    std::vector<double> maxXs;
    std::vector<double> maxYs;
    std::vector<double> maxZs;
    for(auto const& tagger : fTaggers){
      minXs.push_back(tagger.second.minX);
      minYs.push_back(tagger.second.minY);
      minZs.push_back(tagger.second.minZ);
      maxXs.push_back(tagger.second.maxX);
      maxYs.push_back(tagger.second.maxY);
      maxZs.push_back(tagger.second.maxZ);
    }
    limits.push_back(*std::min_element(minXs.begin(), minXs.end()));
    limits.push_back(*std::min_element(minYs.begin(), minYs.end()));
    limits.push_back(*std::min_element(minZs.begin(), minZs.end()));
    limits.push_back(*std::max_element(maxXs.begin(), maxXs.end()));
    limits.push_back(*std::max_element(maxYs.begin(), maxYs.end()));
    limits.push_back(*std::max_element(maxZs.begin(), maxZs.end()));

    return limits;
  }

  size_t CRTGeoAlg::NumTaggers() const
  {
    return fTaggers.size();
  }

  size_t CRTGeoAlg::NumModules() const
  {
    return fModules.size();
  }

  size_t CRTGeoAlg::NumStrips() const
  {
    return fStrips.size();
  }

  std::map<std::string, CRTTaggerGeo> CRTGeoAlg::GetTaggers() const
  {
    return fTaggers;
  }

  std::map<std::string, CRTModuleGeo> CRTGeoAlg::GetModules() const
  {
    return fModules;
  }

  std::map<std::string, CRTStripGeo> CRTGeoAlg::GetStrips() const
  {
    return fStrips;
  }

  CRTTaggerGeo CRTGeoAlg::GetTagger(const std::string taggerName) const
  {
    return fTaggers.at(taggerName);
  }

  CRTModuleGeo CRTGeoAlg::GetModule(const std::string moduleName) const
  {
    return fModules.at(moduleName);
  }

  CRTModuleGeo CRTGeoAlg::GetModule(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName);
  }

  CRTStripGeo CRTGeoAlg::GetStrip(const std::string stripName) const
  {
    return fStrips.at(stripName);
  }

  CRTStripGeo CRTGeoAlg::GetStrip(const uint16_t channel) const
  {
    return fStrips.at(fSiPMs.at(channel).stripName);
  }

  CRTSiPMGeo CRTGeoAlg::GetSiPM(const uint16_t channel) const
  {
    return fSiPMs.at(channel);
  }

  std::string CRTGeoAlg::GetTaggerName(const std::string name) const
  {
    if(fStrips.find(name) != fStrips.end())
      return fModules.at(fStrips.at(name).moduleName).taggerName;
    else if(fModules.find(name) != fModules.end())
      return fModules.at(name).taggerName;
    
    return "";
  }

  std::string CRTGeoAlg::ChannelToStripName(const uint16_t channel) const
  {
    return fSiPMs.at(channel).stripName;
  }

  std::string CRTGeoAlg::ChannelToTaggerName(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).taggerName;
  }

  size_t CRTGeoAlg::ChannelToOrientation(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).orientation;
  }

  std::vector<double> CRTGeoAlg::StripHit3DPos(const std::string stripName, const double x, 
                                               const double ex)
  {
    const CRTStripGeo &strip = fStrips.at(stripName);

    const uint16_t adsID = strip.adsID;
    const uint16_t adID  = fModules.at(strip.moduleName).adID;
  
    const geo::AuxDetSensitiveGeo &auxDetSensitive = fAuxDetGeoCore->AuxDetGeoVec()[adID].SensitiveVolume(adsID);

    double halfWidth  = auxDetSensitive.HalfWidth1();
    double halfHeight = auxDetSensitive.HalfHeight();
    double halfLength = auxDetSensitive.Length()/2.;
  
    double l1[3] = {halfWidth, -halfHeight + x + ex, halfLength};
    double w1[3];
    auxDetSensitive.LocalToWorld(l1, w1);

    double l2[3] = {-halfWidth, -halfHeight + x - ex, -halfLength};
    double w2[3];
    auxDetSensitive.LocalToWorld(l2, w2);

    std::vector<double> limits = {std::min(w1[0],w2[0]), std::max(w1[0],w2[0]),
                                  std::min(w1[1],w2[1]), std::max(w1[1],w2[1]),
                                  std::min(w1[2],w2[2]), std::max(w1[2],w2[2])};
    return limits;
  }

  TVector3 CRTGeoAlg::ChannelToSipmPosition(const uint16_t channel) const
  {
    const CRTSiPMGeo &sipm = fSiPMs.at(channel);
    return {sipm.x, sipm.y, sipm.z};
  }

  std::pair<int, int> CRTGeoAlg::GetStripSipmChannels(const std::string stripName) const
  {
    const CRTStripGeo &strip = fStrips.at(stripName);
    return std::make_pair(strip.channel0, strip.channel1);
  }

  double CRTGeoAlg::DistanceDownStrip(const TVector3 position, const std::string stripName) const
  {
    // === TO-DO ===
    // This assumes that the CRT is arranged such that its three axes map onto the
    // three world axes. This would potentially not always be the case... e.g. using
    // an A-frame. We would need to amend this in that scenario.

    const CRTStripGeo &strip = fStrips.at(stripName);
    double distance = std::numeric_limits<double>::max();

    const TVector3 pos = ChannelToSipmPosition(strip.channel0);

    const double xdiff = std::abs(strip.maxX-strip.minX);
    const double ydiff = std::abs(strip.maxY-strip.minY);
    const double zdiff = std::abs(strip.maxZ-strip.minZ);
  
    if(xdiff > ydiff && xdiff > zdiff) distance = position.X() - pos.X();
    else if(ydiff > xdiff && ydiff > zdiff) distance = position.Y() - pos.Y();
    else if(zdiff > xdiff && zdiff > ydiff) distance = position.Z() - pos.Z();
  
    return std::abs(distance);
  }

  double CRTGeoAlg::DistanceDownStrip(const TVector3 position, const uint16_t channel) const
  {
    const CRTSiPMGeo sipm = fSiPMs.at(channel);
   
    return DistanceDownStrip(position, sipm.stripName);
  }

  bool CRTGeoAlg::CheckOverlap(const CRTStripGeo &strip1, const CRTStripGeo &strip2)
  {
    const double minX = std::max(strip1.minX, strip2.minX);
    const double maxX = std::min(strip1.maxX, strip2.maxX);
    const double minY = std::max(strip1.minY, strip2.minY);
    const double maxY = std::min(strip1.maxY, strip2.maxY);
    const double minZ = std::max(strip1.minZ, strip2.minZ);
    const double maxZ = std::min(strip1.maxZ, strip2.maxZ);

    // If the two strips overlap in 2 dimensions then return true
    return (minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ);
  }
}
