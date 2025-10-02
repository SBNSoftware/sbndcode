#include "CRTGeoService.h"

namespace sbnd::crt {

  CRTGeoService::CRTGeoService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
    : fDefaultGain(pset.get<double>("DefaultGain"))
    , fMC(pset.get<bool>("MC"))
  {
    reg.sPreBeginRun.watch(this, &CRTGeoService::preBeginRun);

    fGeometryService               = lar::providerFrom<geo::Geometry>();
    fAuxDetGeoCore                 = ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr();
    fCRTCalibrationDatabaseService = lar::providerFrom<sbndDB::ICRTCalibrationDatabaseService const>();

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

        art::ServiceHandle<SBND::CRTChannelMapService> ChannelMapService;
        const bool invert = ChannelMapService->GetInversionFromOfflineModuleID(ad_i);

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

            if(std::find(usedModules.begin(), usedModules.end(), moduleName) == usedModules.end())
              {
                const std::string stripName = nodeStrip->GetVolume()->GetName();
                const bool minos = stripName.find("MINOS") != std::string::npos ? true : false;

                usedModules.push_back(moduleName);
                CRTModuleGeo module  = CRTModuleGeo(nodeModule, auxDet, ad_i, taggerName, 0, 0, invert, minos);
                fModules.insert(std::pair<std::string, CRTModuleGeo>(moduleName, module));
              }

            // Fill the strip information
            const std::string stripName = nodeStrip->GetName();
            const uint32_t channel0     = ChannelMapService->ConstructOfflineChannelIDFromOfflineModuleIDAndStripNumber(ad_i, ads_i);
            const uint32_t channel1     = channel0 + 1;

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
            geo::AuxDetSensitiveGeo::LocalPoint_t const sipm0XYZ{sipmX, sipm0Y, 0};
            auto const sipm0XYZWorld = auxDetSensitive.toWorldCoords(sipm0XYZ);

            geo::AuxDetSensitiveGeo::LocalPoint_t const sipm1XYZ{sipmX, sipm1Y, 0};
            auto const sipm1XYZWorld = auxDetSensitive.toWorldCoords(sipm1XYZ);

            // Fill SiPM information
            CRTSiPMGeo sipm0 = CRTSiPMGeo(stripName, channel0, sipm0XYZWorld, 0, fDefaultGain, CRTChannelStatus::kGoodChannel);
            CRTSiPMGeo sipm1 = CRTSiPMGeo(stripName, channel1, sipm1XYZWorld, 0, fDefaultGain, CRTChannelStatus::kGoodChannel);
            fSiPMs.insert(std::pair<uint16_t, CRTSiPMGeo>(channel0, sipm0));
            fSiPMs.insert(std::pair<uint16_t, CRTSiPMGeo>(channel1, sipm1));
          }
      }
  }

  CRTGeoService::~CRTGeoService() {}

  void CRTGeoService::preBeginRun(art::Run const& run)
  {
    if(fMC)
      return;

    art::ServiceHandle<SBND::CRTChannelMapService> ChannelMapService;

    for(auto&& [ name, module ] : fModules)
      {
        const unsigned int mac5 = ChannelMapService->GetMAC5FromOfflineModuleID(module.adID);

        module.t0DelayCorrection = fCRTCalibrationDatabaseService->getT0CableLengthOffset(mac5) + fCRTCalibrationDatabaseService->getT0CalibratedOffset(mac5);
        module.t1DelayCorrection = fCRTCalibrationDatabaseService->getT1CableLengthOffset(mac5) + fCRTCalibrationDatabaseService->getT1CalibratedOffset(mac5);
      }

    for(auto&& [ channel, sipm ] : fSiPMs)
      {
        const unsigned int online_channel_id = ChannelMapService->GetOnlineChannelIDFromOfflineChannelID(sipm.channel);

        sipm.pedestal = fCRTCalibrationDatabaseService->getPedestal(online_channel_id);
        sipm.status   = fCRTCalibrationDatabaseService->getChannelStatus(online_channel_id);
        sipm.gain     = fCRTCalibrationDatabaseService->getGainFactor(online_channel_id);
      }
  }

  std::vector<double> CRTGeoService::CRTLimits() const
  {
    std::vector<double> limits;

    std::vector<double> minXs, minYs, minZs, maxXs, maxYs, maxZs;

    for(auto const& tagger : fTaggers)
      {
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

  size_t CRTGeoService::NumTaggers() const
  {
    return fTaggers.size();
  }

  size_t CRTGeoService::NumModules() const
  {
    return fModules.size();
  }

  size_t CRTGeoService::NumStrips() const
  {
    return fStrips.size();
  }

  size_t CRTGeoService::NumSiPMs() const
  {
    return fSiPMs.size();
  }

  std::map<std::string, CRTTaggerGeo> CRTGeoService::GetTaggers() const
  {
    return fTaggers;
  }

  std::map<std::string, CRTModuleGeo> CRTGeoService::GetModules() const
  {
    return fModules;
  }

  std::map<std::string, CRTStripGeo> CRTGeoService::GetStrips() const
  {
    return fStrips;
  }

  std::map<uint16_t, CRTSiPMGeo> CRTGeoService::GetSiPMs() const
  {
    return fSiPMs;
  }

  CRTTaggerGeo CRTGeoService::GetTagger(const std::string taggerName) const
  {
    return fTaggers.at(taggerName);
  }

  CRTModuleGeo CRTGeoService::GetModule(const std::string moduleName) const
  {
    return fModules.at(moduleName);
  }

  CRTModuleGeo CRTGeoService::GetModule(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName);
  }

  CRTModuleGeo CRTGeoService::GetModuleByAuxDetIndex(const unsigned ad_i) const
  {
    for(auto const &[name, module] : fModules)
      {
        if(module.adID == ad_i)
          return module;
      }

    CRTModuleGeo void_return;
    return void_return;
  }

  CRTStripGeo CRTGeoService::GetStrip(const std::string stripName) const
  {
    return fStrips.at(stripName);
  }

  CRTStripGeo CRTGeoService::GetStrip(const uint16_t channel) const
  {
    return fStrips.at(fSiPMs.at(channel).stripName);
  }

  CRTStripGeo CRTGeoService::GetStripByAuxDetIndices(const unsigned ad_i, const unsigned ads_i) const
  {
    const CRTModuleGeo module = GetModule(ad_i);
    const uint16_t channel = 32 * ad_i + 2 * ads_i;

    return GetStrip(channel);
  }

  CRTSiPMGeo CRTGeoService::GetSiPM(const uint16_t channel) const
  {
    return fSiPMs.at(channel);
  }

  std::string CRTGeoService::GetTaggerName(const std::string name) const
  {
    if(fStrips.find(name) != fStrips.end())
      return fModules.at(fStrips.at(name).moduleName).taggerName;
    else if(fModules.find(name) != fModules.end())
      return fModules.at(name).taggerName;

    return "";
  }

  std::string CRTGeoService::ChannelToStripName(const uint16_t channel) const
  {
    return fSiPMs.at(channel).stripName;
  }

  std::string CRTGeoService::ChannelToTaggerName(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).taggerName;
  }

  enum CRTTagger CRTGeoService::ChannelToTaggerEnum(const uint16_t channel) const
  {
    return CRTCommonUtils::GetTaggerEnum(ChannelToTaggerName(channel));
  }

  enum CRTTagger CRTGeoService::AuxDetIndexToTaggerEnum(const unsigned ad_i) const
  {
    return CRTCommonUtils::GetTaggerEnum(GetModuleByAuxDetIndex(ad_i).taggerName);
  }

  size_t CRTGeoService::ChannelToOrientation(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).orientation;
  }

  std::array<double, 6> CRTGeoService::StripHit3DPos(const uint16_t channel, const double x,
                                                     const double ex)
  {
    const CRTStripGeo &strip = GetStrip(channel);

    const uint16_t adsID = strip.adsID;
    const uint16_t adID  = fModules.at(strip.moduleName).adID;

    const geo::AuxDetSensitiveGeo &auxDetSensitive = fAuxDetGeoCore->AuxDetGeoVec()[adID].SensitiveVolume(adsID);

    double halfWidth  = auxDetSensitive.HalfWidth1();
    double halfHeight = auxDetSensitive.HalfHeight();
    double halfLength = auxDetSensitive.Length()/2.;

    geo::AuxDetSensitiveGeo::LocalPoint_t l1{halfWidth, -halfHeight + x + ex, halfLength};
    auto const w1 = auxDetSensitive.toWorldCoords(l1);

    geo::AuxDetSensitiveGeo::LocalPoint_t l2{-halfWidth, -halfHeight + x - ex, -halfLength};
    auto const w2 = auxDetSensitive.toWorldCoords(l2);

    std::array<double, 6> limits = {std::min(w1.X(),w2.X()), std::max(w1.X(),w2.X()),
                                    std::min(w1.Y(),w2.Y()), std::max(w1.Y(),w2.Y()),
                                    std::min(w1.Z(),w2.Z()), std::max(w1.Z(),w2.Z())};
    return limits;
  }

  std::vector<double> CRTGeoService::StripWorldToLocalPos(const CRTStripGeo &strip, const double x,
                                                          const double y, const double z)
  {
    const uint16_t adsID = strip.adsID;
    const uint16_t adID  = fModules.at(strip.moduleName).adID;

    const geo::AuxDetSensitiveGeo &auxDetSensitive = fAuxDetGeoCore->AuxDetGeoVec()[adID].SensitiveVolume(adsID);

    geo::Point_t const worldpos{x, y, z};
    auto const localpos = auxDetSensitive.toLocalCoords(worldpos);

    std::vector<double> localvec = {localpos.X(), localpos.Y(), localpos.Z()};
    return localvec;
  }

  std::vector<double> CRTGeoService::StripWorldToLocalPos(const uint16_t channel, const double x,
                                                          const double y, const double z)
  {
    const CRTStripGeo strip = GetStrip(channel);
    return StripWorldToLocalPos(strip, x, y, z);
  }

  std::array<double, 6> CRTGeoService::FEBWorldPos(const CRTModuleGeo &module)
  {
    const geo::AuxDetGeo &auxDet = fAuxDetGeoCore->AuxDetGeoVec()[module.adID];

    const double halfWidth  = auxDet.HalfWidth1();
    const double halfLength = auxDet.Length() / 2.;

    geo::AuxDetGeo::LocalPoint_t limits { halfWidth + 5, -15, -halfLength - 5};
    geo::AuxDetGeo::LocalPoint_t limits2{ halfWidth,     15,  -halfLength + 5};

    if(!module.top)
      {
        limits.SetX(-limits.X());
        limits2.SetX(-limits2.X());
      }

    auto const limitsWorld  = auxDet.toWorldCoords(limits);
    auto const limitsWorld2 = auxDet.toWorldCoords(limits2);

    const double minX = std::min(limitsWorld.X(), limitsWorld2.X());
    const double maxX = std::max(limitsWorld.X(), limitsWorld2.X());
    const double minY = std::min(limitsWorld.Y(), limitsWorld2.Y());
    const double maxY = std::max(limitsWorld.Y(), limitsWorld2.Y());
    const double minZ = std::min(limitsWorld.Z(), limitsWorld2.Z());
    const double maxZ = std::max(limitsWorld.Z(), limitsWorld2.Z());

    return {minX, maxX, minY, maxY, minZ, maxZ};
  }

  std::array<double, 6> CRTGeoService::FEBChannel0WorldPos(const CRTModuleGeo &module)
  {
    const geo::AuxDetGeo &auxDet = fAuxDetGeoCore->AuxDetGeoVec()[module.adID];

    const double halfWidth  = auxDet.HalfWidth1();
    const double halfLength = auxDet.Length() / 2.;

    geo::AuxDetGeo::LocalPoint_t limits { halfWidth + 5, -15, -halfLength - 5};
    geo::AuxDetGeo::LocalPoint_t limits2{ halfWidth,     -12, -halfLength + 5};

    if(!module.top)
      {
        limits.SetX(-limits.X());
        limits2.SetX(-limits2.X());
      }

    if(module.invertedOrdering)
      {
        limits.SetY(-limits.Y());
        limits2.SetY(-limits2.Y());
      }

    auto const limitsWorld  = auxDet.toWorldCoords(limits);
    auto const limitsWorld2 = auxDet.toWorldCoords(limits2);

    const double minX = std::min(limitsWorld.X(), limitsWorld2.X());
    const double maxX = std::max(limitsWorld.X(), limitsWorld2.X());
    const double minY = std::min(limitsWorld.Y(), limitsWorld2.Y());
    const double maxY = std::max(limitsWorld.Y(), limitsWorld2.Y());
    const double minZ = std::min(limitsWorld.Z(), limitsWorld2.Z());
    const double maxZ = std::max(limitsWorld.Z(), limitsWorld2.Z());

    return {minX, maxX, minY, maxY, minZ, maxZ};
  }

  geo::Point_t CRTGeoService::ChannelToSipmPosition(const uint16_t channel) const
  {
    const CRTSiPMGeo &sipm = fSiPMs.at(channel);
    return {sipm.x, sipm.y, sipm.z};
  }

  std::pair<int, int> CRTGeoService::GetStripSipmChannels(const std::string stripName) const
  {
    const CRTStripGeo &strip = fStrips.at(stripName);
    return std::make_pair(strip.channel0, strip.channel1);
  }

  double CRTGeoService::DistanceDownStrip(const geo::Point_t position, const std::string stripName) const
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

  double CRTGeoService::DistanceDownStrip(const geo::Point_t position, const uint16_t channel) const
  {
    const CRTSiPMGeo sipm = fSiPMs.at(channel);

    return DistanceDownStrip(position, sipm.stripName);
  }

  bool CRTGeoService::CheckOverlap(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer)
  {
    const CRTTagger tagger1 = CRTCommonUtils::GetTaggerEnum(ChannelToTaggerName(strip1.channel0));
    const CRTTagger tagger2 = CRTCommonUtils::GetTaggerEnum(ChannelToTaggerName(strip2.channel0));

    if(tagger1 != tagger2)
      return false;

    const double minX = std::max(strip1.minX, strip2.minX) - overlap_buffer / 2.;
    const double maxX = std::min(strip1.maxX, strip2.maxX) + overlap_buffer / 2.;
    const double minY = std::max(strip1.minY, strip2.minY) - overlap_buffer / 2.;
    const double maxY = std::min(strip1.maxY, strip2.maxY) + overlap_buffer / 2.;
    const double minZ = std::max(strip1.minZ, strip2.minZ) - overlap_buffer / 2.;
    const double maxZ = std::min(strip1.maxZ, strip2.maxZ) + overlap_buffer / 2.;

    const CoordSet  taggerCoord = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger1);

    if(taggerCoord == kX)
      return minY<maxY && minZ<maxZ;
    else if(taggerCoord == kY)
      return minX<maxX && minZ<maxZ;
    else if(taggerCoord == kZ)
      return minX<maxX && minY<maxY;
    else
      {
        std::cout << "Tagger coordinate ill-defined: " << taggerCoord << std::endl;
        return false;
      }
  }

  bool CRTGeoService::CheckOverlap(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer)
  {
    CRTStripGeo strip1 = GetStrip(channel1);
    CRTStripGeo strip2 = GetStrip(channel2);

    return CheckOverlap(strip1, strip2, overlap_buffer);
  }

  bool CRTGeoService::AdjacentStrips(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer)
  {
    CRTModuleGeo module1 = GetModule(strip1.channel0);
    CRTModuleGeo module2 = GetModule(strip2.channel0);

    if(module1.taggerName != module2.taggerName || module1.orientation != module2.orientation)
      return false;

    const double minX = std::max(strip1.minX, strip2.minX);
    const double maxX = std::min(strip1.maxX, strip2.maxX);
    const double minY = std::max(strip1.minY, strip2.minY);
    const double maxY = std::min(strip1.maxY, strip2.maxY);
    const double minZ = std::max(strip1.minZ, strip2.minZ);
    const double maxZ = std::min(strip1.maxZ, strip2.maxZ);

    if(strip1.minX == strip2.minX && strip1.maxX == strip2.maxX &&
       strip1.minY == strip2.minY && strip1.maxY == strip2.maxY)
      return minZ < (maxZ + overlap_buffer);
    if(strip1.minX == strip2.minX && strip1.maxX == strip2.maxX &&
       strip1.minZ == strip2.minZ && strip1.maxZ == strip2.maxZ)
      return minY < (maxY + overlap_buffer);
    if(strip1.minY == strip2.minY && strip1.maxY == strip2.maxY &&
       strip1.minZ == strip2.minZ && strip1.maxZ == strip2.maxZ)
      return minX < (maxX + overlap_buffer);
    else
      return false;
  }

  bool CRTGeoService::AdjacentStrips(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer)
  {
    CRTStripGeo strip1 = GetStrip(channel1);
    CRTStripGeo strip2 = GetStrip(channel2);

    return AdjacentStrips(strip1, strip2, overlap_buffer);
  }

  bool CRTGeoService::DifferentOrientations(const CRTStripGeo &strip1, const CRTStripGeo &strip2)
  {
    const CRTModuleGeo module1 = GetModule(strip1.moduleName);
    const CRTModuleGeo module2 = GetModule(strip2.moduleName);

    return module1.orientation != module2.orientation;
  }

  enum CRTTagger CRTGeoService::WhichTagger(const double &x, const double &y, const double &z, const double &buffer)
  {
    for(auto const& [name, tagger] : fTaggers)
      {
        if(x > tagger.minX - buffer &&
           x < tagger.maxX + buffer &&
           y > tagger.minY - buffer &&
           y < tagger.maxY + buffer &&
           z > tagger.minZ - buffer &&
           z < tagger.maxZ + buffer)
          return CRTCommonUtils::GetTaggerEnum(name);
      }
    return kUndefinedTagger;
  }

  enum CoordSet CRTGeoService::GlobalConstrainedCoordinates(const uint16_t channel)
  {
    const std::string taggerName = ChannelToTaggerName(channel);
    const CRTTagger tagger       = CRTCommonUtils::GetTaggerEnum(taggerName);
    const uint16_t orientation   = GetModule(channel).orientation;

    const CoordSet widthdir    = CRTCommonUtils::GetStripWidthGlobalCoordinate(tagger, orientation);
    const CoordSet taggercoord = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger);

    return widthdir | taggercoord;
  }

  bool CRTGeoService::IsPointInsideCRTLimits(const geo::Point_t &point)
  {
    const std::vector<double> lims = CRTLimits();

    return (point.X() > lims[0] && point.X() < lims[3]) &&
           (point.Y() > lims[1] && point.Y() < lims[4]) &&
           (point.Z() > lims[2] && point.Z() < lims[5]);
  }
}

DEFINE_ART_SERVICE(sbnd::crt::CRTGeoService)
