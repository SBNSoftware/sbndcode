#include "CRTGeoAlg.h"

namespace sbnd::crt {

  CRTGeoAlg::CRTGeoAlg(fhicl::ParameterSet const &p) :
    CRTGeoAlg(p, lar::providerFrom<geo::Geometry>(),
              ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr())
  {}

  CRTGeoAlg::CRTGeoAlg(fhicl::ParameterSet const &p, geo::GeometryCore const *geometry,
                       geo::AuxDetGeometryCore const *auxdet_geometry)
    : fT0CableLengthCorrectionsVector(p.get<std::vector<std::pair<unsigned, double>>>("T0CableLengthCorrections", std::vector<std::pair<unsigned, double>>()))
    , fT1CableLengthCorrectionsVector(p.get<std::vector<std::pair<unsigned, double>>>("T1CableLengthCorrections", std::vector<std::pair<unsigned, double>>()))
    , fDefaultPedestal(p.get<double>("DefaultPedestal", 0.))
    , fSiPMPedestalsVector(p.get<std::vector<std::pair<unsigned, double>>>("SiPMPedestals", std::vector<std::pair<unsigned, double>>()))
    , fDefaultGain(p.get<double>("DefaultGain", 0.025))
    , fSiPMGainsVector(p.get<std::vector<std::pair<unsigned, double>>>("SiPMGains", std::vector<std::pair<unsigned, double>>()))
    , fChannelInversionVector(p.get<std::vector<std::pair<unsigned, bool>>>("InvertedChannelOrder", std::vector<std::pair<unsigned, bool>>()))
  {
    fGeometryService = geometry;
    fAuxDetGeoCore   = auxdet_geometry;
    TGeoManager* manager = fGeometryService->ROOTGeoManager();

    fT0CableLengthCorrections = std::map<unsigned, double>(fT0CableLengthCorrectionsVector.begin(),
                                                           fT0CableLengthCorrectionsVector.end());
    fT1CableLengthCorrections = std::map<unsigned, double>(fT1CableLengthCorrectionsVector.begin(),
                                                           fT1CableLengthCorrectionsVector.end());
    fSiPMPedestals            = std::map<unsigned, double>(fSiPMPedestalsVector.begin(),
                                                           fSiPMPedestalsVector.end());
    fSiPMGains                = std::map<unsigned, double>(fSiPMGainsVector.begin(),
                                                           fSiPMGainsVector.end());
    fChannelInversion         = std::map<unsigned, bool>(fChannelInversionVector.begin(),
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
                const int32_t t0CableDelayCorrection = fT0CableLengthCorrections.size() ?
                  fT0CableLengthCorrections.at(ad_i) : 0;
                const int32_t t1CableDelayCorrection = fT1CableLengthCorrections.size() ?
                  fT1CableLengthCorrections.at(ad_i) : 0;

                const std::string stripName = nodeStrip->GetVolume()->GetName();
                const bool minos = stripName.find("MINOS") != std::string::npos ? true : false;

                usedModules.push_back(moduleName);
                CRTModuleGeo module  = CRTModuleGeo(nodeModule, auxDet, ad_i, taggerName,
                                                    t0CableDelayCorrection, t1CableDelayCorrection,
                                                    invert, minos);
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
            geo::AuxDetSensitiveGeo::LocalPoint_t const sipm0XYZ{sipmX, sipm0Y, 0};
            auto const sipm0XYZWorld = auxDetSensitive.toWorldCoords(sipm0XYZ);

            geo::AuxDetSensitiveGeo::LocalPoint_t const sipm1XYZ{sipmX, sipm1Y, 0};
            auto const sipm1XYZWorld = auxDetSensitive.toWorldCoords(sipm1XYZ);

            const uint32_t pedestal0 = fSiPMPedestals.size() ? fSiPMPedestals.at(channel0) : fDefaultPedestal;
            const uint32_t pedestal1 = fSiPMPedestals.size() ? fSiPMPedestals.at(channel1) : fDefaultPedestal;

            const double gain0 = fSiPMGains.size() ? fSiPMGains.at(channel0) : fDefaultGain;
            const double gain1 = fSiPMGains.size() ? fSiPMGains.at(channel1) : fDefaultGain;

            // Fill SiPM information
            CRTSiPMGeo sipm0 = CRTSiPMGeo(stripName, channel0, sipm0XYZWorld, pedestal0, gain0);
            CRTSiPMGeo sipm1 = CRTSiPMGeo(stripName, channel1, sipm1XYZWorld, pedestal1, gain1);
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

  size_t CRTGeoAlg::NumSiPMs() const
  {
    return fSiPMs.size();
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

  std::map<uint16_t, CRTSiPMGeo> CRTGeoAlg::GetSiPMs() const
  {
    return fSiPMs;
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

  CRTModuleGeo CRTGeoAlg::GetModuleByAuxDetIndex(const unsigned ad_i) const
  {
    for(auto const &[name, module] : fModules)
      {
        if(module.adID == ad_i)
          return module;
      }

    CRTModuleGeo void_return;
    return void_return;
  }

  CRTStripGeo CRTGeoAlg::GetStrip(const std::string stripName) const
  {
    return fStrips.at(stripName);
  }

  CRTStripGeo CRTGeoAlg::GetStrip(const uint16_t channel) const
  {
    return fStrips.at(fSiPMs.at(channel).stripName);
  }

  CRTStripGeo CRTGeoAlg::GetStripByAuxDetIndices(const unsigned ad_i, const unsigned ads_i) const
  {
    const CRTModuleGeo module = GetModule(ad_i);
    const uint16_t channel =
      module.invertedOrdering ? 32 * ad_i + (31 -2 *ads_i) : 32 * ad_i + 2 * ads_i;

    return GetStrip(channel);
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

  enum CRTTagger CRTGeoAlg::ChannelToTaggerEnum(const uint16_t channel) const
  {
    return CRTCommonUtils::GetTaggerEnum(ChannelToTaggerName(channel));
  }

  size_t CRTGeoAlg::ChannelToOrientation(const uint16_t channel) const
  {
    return fModules.at(fStrips.at(fSiPMs.at(channel).stripName).moduleName).orientation;
  }

  std::array<double, 6> CRTGeoAlg::StripHit3DPos(const uint16_t channel, const double x,
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

  std::vector<double> CRTGeoAlg::StripWorldToLocalPos(const CRTStripGeo &strip, const double x, 
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

  std::vector<double> CRTGeoAlg::StripWorldToLocalPos(const uint16_t channel, const double x,
                                                      const double y, const double z)
  {
    const CRTStripGeo strip = GetStrip(channel);
    return StripWorldToLocalPos(strip, x, y, z);
  }

  std::array<double, 6> CRTGeoAlg::FEBWorldPos(const CRTModuleGeo &module)
  {
    const geo::AuxDetGeo &auxDet = fAuxDetGeoCore->AuxDetGeoVec()[module.adID];

    const double halfWidth  = auxDet.HalfWidth1();
    const double halfLength = auxDet.Length() / 2.;

    geo::AuxDetGeo::LocalPoint_t limits { - halfWidth - 5, -15, -halfLength - 5};
    geo::AuxDetGeo::LocalPoint_t limits2{ - halfWidth,     15,  -halfLength + 5};

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

  std::array<double, 6> CRTGeoAlg::FEBChannel0WorldPos(const CRTModuleGeo &module)
  {
    const geo::AuxDetGeo &auxDet = fAuxDetGeoCore->AuxDetGeoVec()[module.adID];

    const double halfWidth  = auxDet.HalfWidth1();
    const double halfLength = auxDet.Length() / 2.;

    geo::AuxDetGeo::LocalPoint_t limits { - halfWidth - 5, -15, -halfLength - 5};
    geo::AuxDetGeo::LocalPoint_t limits2{ - halfWidth,     -12, -halfLength + 5};

    if(!module.top)
      {
        limits.SetX(-limits.X());
        limits.SetX(-limits2.X());
      }

    if(module.invertedOrdering)
      {
        limits.SetY(-limits.Y());
        limits.SetY(-limits2.Y());
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

  double CRTGeoAlg::DistanceDownStrip(const geo::Point_t position, const uint16_t channel) const
  {
    const CRTSiPMGeo sipm = fSiPMs.at(channel);

    return DistanceDownStrip(position, sipm.stripName);
  }

  bool CRTGeoAlg::CheckOverlap(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer)
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

  bool CRTGeoAlg::CheckOverlap(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer)
  {
    CRTStripGeo strip1 = GetStrip(channel1);
    CRTStripGeo strip2 = GetStrip(channel2);

    return CheckOverlap(strip1, strip2, overlap_buffer);
  }

  bool CRTGeoAlg::AdjacentStrips(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer)
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

  bool CRTGeoAlg::AdjacentStrips(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer)
  {
    CRTStripGeo strip1 = GetStrip(channel1);
    CRTStripGeo strip2 = GetStrip(channel2);

    return AdjacentStrips(strip1, strip2, overlap_buffer);
  }

  bool CRTGeoAlg::DifferentOrientations(const CRTStripGeo &strip1, const CRTStripGeo &strip2)
  {
    const CRTModuleGeo module1 = GetModule(strip1.moduleName);
    const CRTModuleGeo module2 = GetModule(strip2.moduleName);

    return module1.orientation != module2.orientation;
  }

  enum CRTTagger CRTGeoAlg::WhichTagger(const double &x, const double &y, const double &z, const double &buffer)
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

  enum CoordSet CRTGeoAlg::GlobalConstrainedCoordinates(const uint16_t channel)
  {
    const std::string taggerName = ChannelToTaggerName(channel);
    const CRTTagger tagger       = CRTCommonUtils::GetTaggerEnum(taggerName);
    const uint16_t orientation   = GetModule(channel).orientation;

    const CoordSet widthdir    = CRTCommonUtils::GetStripWidthGlobalCoordinate(tagger, orientation);
    const CoordSet taggercoord = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger);

    return widthdir | taggercoord;
  }

  bool CRTGeoAlg::IsPointInsideCRTLimits(const geo::Point_t &point)
  {
    const std::vector<double> lims = CRTLimits();

    return point.X() > lims[0] &&
           point.X() < lims[3] &&
           point.Y() > lims[1] &&
           point.Y() < lims[4] &&
           point.Z() > lims[2] &&
           point.Z() < lims[5];
  }
}
