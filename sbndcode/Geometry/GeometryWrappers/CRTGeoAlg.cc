#include "CRTGeoAlg.h"

namespace sbnd{

  // Constructor - get values from the auxdet geometry service
  CRTGeoAlg::CRTGeoAlg() :
    CRTGeoAlg::CRTGeoAlg(lar::providerFrom<geo::Geometry>(), 
                         ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr())
  {}

  CRTGeoAlg::CRTGeoAlg(fhicl::ParameterSet const &p)
  {
    fCableLengthCorrectionsVector = p.get<std::vector<std::pair<unsigned, double>>>("CableLengthCorrectionsMap", std::vector<std::pair<unsigned, double>>());
    fCableLengthCorrections       = std::map<unsigned, double>(fCableLengthCorrectionsVector.begin(), fCableLengthCorrectionsVector.end());

    fSiPMPedestalsVector          = p.get<std::vector<std::pair<unsigned, double>>>("SiPMPedestalsVector", std::vector<std::pair<unsigned, double>>());
    fSiPMPedestals                = std::map<unsigned, double>(fSiPMPedestalsVector.begin(), fSiPMPedestalsVector.end());

    CRTGeoAlg(lar::providerFrom<geo::Geometry>(), 
	      ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr());
  }

  CRTGeoAlg::CRTGeoAlg(geo::GeometryCore const *geometry, 
                       geo::AuxDetGeometryCore const *auxdet_geometry) 
  {
    fGeometryService = geometry;
    fAuxDetGeoCore   = auxdet_geometry;
    TGeoManager* manager = fGeometryService->ROOTGeoManager();

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
                fTaggers.at(taggerName) = CRTTaggerGeo(nodeTagger, nodeDet);
              }

            // Fill the module information
            const std::string moduleName = nodeModule->GetName();
            if(std::find(usedModules.begin(), usedModules.end(), moduleName) == usedModules.end())
              {
		const uint32_t cableDelayCorrection = fCableLengthCorrections.size() ? 
		  fCableLengthCorrections.at(32 * ad_i) : 0;
                usedModules.push_back(moduleName);
                fModules.at(moduleName) = CRTModuleGeo(nodeModule, auxDet, ad_i, taggerName,
						       cableDelayCorrection);
              }

            // Fill the strip information
            const std::string stripName = nodeStrip->GetName();
            const uint32_t channel0 = 32 * ad_i + 2 * ads_i + 0;
            const uint32_t channel1 = 32 * ad_i + 2 * ads_i + 1;          
            if(std::find(usedStrips.begin(), usedStrips.end(), stripName) == usedStrips.end())
              {
                usedStrips.push_back(stripName);
                fStrips.at(stripName) = CRTStripGeo(nodeStrip, auxDetSensitive, ads_i, moduleName,
                                                    channel0, channel1);
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
            fSiPMs.at(channel0) = CRTSiPMGeo(stripName, channel0, sipm0XYZWorld, pedestal0);
            fSiPMs.at(channel1) = CRTSiPMGeo(stripName, channel1, sipm0XYZWorld, pedestal1);
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

  size_t CRTGeoAlg::ChannelToPlaneID(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).planeID;
  }

  std::vector<double> CRTGeoAlg::StripLimitsWithChargeSharing(const std::string stripName, const double x, 
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

  geo::Point_t CRTGeoAlg::ChannelToSipmPosition(const uint16_t channel) const
  {
    const CRTSiPMGeo &sipm = fSiPMs.at(channel);
    return {sipm.x, sipm.y, sipm.z};
  }

  std::pair<int, int> CRTGeoAlg::GetStripSipmChannels(const std::string stripName) const
  {
    const CRTStripGeo &strip = fStrips.at(stripName);
    return std::make_pair(strip.channel0, strip.channel1);
  }

  double CRTGeoAlg::DistanceDownStrip(const geo::Point_t position, const std::string stripName) const
  {
    // === TO-DO ===
    // This assumes that the CRT is arranged such that its three axes map onto the
    // three world axes. This would potentially not always be the case... e.g. using
    // an A-frame. We would need to amend this in that scenario.

    const CRTStripGeo &strip = fStrips.at(stripName);
    double distance = std::numeric_limits<double>::max();

    const geo::Point_t pos = ChannelToSipmPosition(strip.channel0);

    const double xdiff = std::abs(strip.maxX-strip.minX);
    const double ydiff = std::abs(strip.maxY-strip.minY);
    const double zdiff = std::abs(strip.maxZ-strip.minZ);
  
    if(xdiff > ydiff && xdiff > zdiff) distance = position.X() - pos.X();
    else if(ydiff > xdiff && ydiff > zdiff) distance = position.Y() - pos.Y();
    else if(zdiff > xdiff && zdiff > ydiff) distance = position.Z() - pos.Z();
  
    return std::abs(distance);
  }

  bool CRTGeoAlg::CheckOverlap(const CRTModuleGeo &module1, const CRTModuleGeo &module2)
  {
    const double minX = std::max(module1.minX, module2.minX);
    const double maxX = std::min(module1.maxX, module2.maxX);
    const double minY = std::max(module1.minY, module2.minY);
    const double maxY = std::min(module1.maxY, module2.maxY);
    const double minZ = std::max(module1.minZ, module2.minZ);
    const double maxZ = std::min(module1.maxZ, module2.maxZ);

    // If the two strips overlap in 2 dimensions then return true
    return (minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ);
  }

  bool CRTGeoAlg::HasOverlap(const CRTModuleGeo &module)
  {
    const size_t planeID         = module.planeID;
    const std::string taggerName = module.taggerName;

    for(auto const &[moduleName, module2] : fModules)
      {
	
        if(module.taggerName != taggerName || module2.planeID == planeID) 
          continue;

        if(CheckOverlap(module, module2)) 
          return true;
      }
    return false;
  }

  bool CRTGeoAlg::StripHasOverlap(const std::string stripName)
  {
    return HasOverlap(fModules.at(fStrips.at(stripName).moduleName));
  }
}
