#include "CRTGeoAlg.h"

namespace sbnd{

// Constructor - get values from the auxdet geometry service
CRTGeoAlg::CRTGeoAlg():
  CRTGeoAlg::CRTGeoAlg(lar::providerFrom<geo::Geometry>(), ((const geo::AuxDetGeometry*)&(*art::ServiceHandle<geo::AuxDetGeometry>()))->GetProviderPtr())
{}

CRTGeoAlg::CRTGeoAlg(geo::GeometryCore const *geometry, geo::AuxDetGeometryCore const *auxdet_geometry) {
  fGeometryService = geometry;
  fAuxDetGeoCore = auxdet_geometry;

  // Keep track of which objects we've recorded
  std::vector<std::string> usedTaggers;
  std::vector<std::string> usedModules;
  std::vector<std::string> usedStrips;

  // Get the auxdets (strip arrays for some reason)
  const std::vector<geo::AuxDetGeo>& auxDets = fAuxDetGeoCore->AuxDetGeoVec();

  int ad_i = 0;
  // Loop over them
  for(auto const& auxDet : auxDets){

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

    int sv_i = 0;
    // Loop over the strips in the arrays
    for(size_t i = 0; i < auxDet.NSensitiveVolume(); i++){

      geo::AuxDetSensitiveGeo const& auxDetSensitive = auxDet.SensitiveVolume(i);

      // Access the path
      TGeoManager* manager = fGeometryService->ROOTGeoManager();
      manager->cd(path.c_str());

      // We get the array of strips first, which is the AuxDet,
      // then from the AuxDet, we get the strip by picking the
      // daughter with the ID of the AuxDetSensitive, and finally
      // from the AuxDet, we go up and pick the module, tagger, and det
      TGeoNode* nodeArray = manager->GetCurrentNode();
      TGeoNode* nodeStrip = nodeArray->GetDaughter(i);
      TGeoNode* nodeModule = manager->GetMother(1);
      TGeoNode* nodeTagger = manager->GetMother(2);
      TGeoNode* nodeDet = manager->GetMother(3);

      // Fill the tagger information
      std::string taggerName = nodeTagger->GetName();
      if(std::find(usedTaggers.begin(), usedTaggers.end(), taggerName) == usedTaggers.end()){
        usedTaggers.push_back(taggerName);

        // Get the limits in local coords
        double halfWidth = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDX();
        double halfHeight = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDY();
        double halfLength = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDZ()/2;

        // Transform those limits to world coords
        double limits[3] = {halfWidth, halfHeight, halfLength};
        double limitsDet[3];
        nodeTagger->LocalToMaster(limits, limitsDet);
        double limitsWorld[3];
        nodeDet->LocalToMaster(limitsDet, limitsWorld);

        double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
        double limitsDet2[3];
        nodeTagger->LocalToMaster(limits2, limitsDet2);
        double limitsWorld2[3];
        nodeDet->LocalToMaster(limitsDet2, limitsWorld2);

        // Create a tagger geometry object
        CRTTaggerGeo tagger = {};
        tagger.name = taggerName;
        tagger.minX = std::min(limitsWorld[0], limitsWorld2[0]);
        tagger.maxX = std::max(limitsWorld[0], limitsWorld2[0]);
        tagger.minY = std::min(limitsWorld[1], limitsWorld2[1]);
        tagger.maxY = std::max(limitsWorld[1], limitsWorld2[1]);
        tagger.minZ = std::min(limitsWorld[2], limitsWorld2[2]);
        tagger.maxZ = std::max(limitsWorld[2], limitsWorld2[2]);
        tagger.null = false;
        fTaggers[taggerName] = tagger;
      }

      // Fill the module information
      std::string moduleName = nodeModule->GetName();
      if(std::find(usedModules.begin(), usedModules.end(), moduleName) == usedModules.end()){
        usedModules.push_back(moduleName);

        // Technically the auxdet is the strip array but this is basically the same as the module
        // Get the limits in local coordinates
        double halfWidth = auxDet.HalfWidth1();
        double halfHeight = auxDet.HalfHeight();
        double halfLength = auxDet.Length()/2;

        // Transform to world coordinates
        double limits[3] = {halfWidth, halfHeight, halfLength};
        double limitsWorld[3];
        auxDet.LocalToWorld(limits, limitsWorld);

        double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
        double limitsWorld2[3];
        auxDet.LocalToWorld(limits2, limitsWorld2);

        // Determine which plane the module is in in the tagger (XY configuration)
        double origin[3] = {0, 0, 0};
        double modulePosMother[3];
        nodeModule->LocalToMaster(origin, modulePosMother);
        size_t planeID = (modulePosMother[2] > 0);
        // Are the sipms at the top or bottom
        bool top = (planeID == 1) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

        // Create a module geometry object
        CRTModuleGeo module = {};
        module.name = moduleName;
        module.auxDetID = ad_i;
        module.minX = std::min(limitsWorld[0], limitsWorld2[0]);
        module.maxX = std::max(limitsWorld[0], limitsWorld2[0]);
        module.minY = std::min(limitsWorld[1], limitsWorld2[1]);
        module.maxY = std::max(limitsWorld[1], limitsWorld2[1]);
        module.minZ = std::min(limitsWorld[2], limitsWorld2[2]);
        module.maxZ = std::max(limitsWorld[2], limitsWorld2[2]);
        module.normal = auxDet.GetNormalVector();
        module.null = false;
        module.planeID = planeID;
        module.top = top;
        module.tagger = taggerName;
        fModules[moduleName] = module;
      }

      // Fill the strip information
      std::string stripName = nodeStrip->GetName();
      if(std::find(usedStrips.begin(), usedStrips.end(), stripName) == usedStrips.end()){
        usedStrips.push_back(stripName);

        // Get the limits in local coordinates
	// Note that these dimensions DO NOT conform to the expected mapping
	// width  = conventional strip length (one end to another)
	// height = conventional strip width (distance between the two SiPMs)
	// length = conventional thickness
        double halfWidth = auxDetSensitive.HalfWidth1();
        double halfHeight = auxDetSensitive.HalfHeight();
        double halfLength = auxDetSensitive.Length()/2;

        // Transform to world coordinates
        double limits[3] = {halfWidth, halfHeight, halfLength};
        double limitsWorld[3];
        auxDetSensitive.LocalToWorld(limits, limitsWorld);

        double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
        double limitsWorld2[3];
        auxDetSensitive.LocalToWorld(limits2, limitsWorld2);

        // Create a strip geometry object
        CRTStripGeo strip = {};
        strip.name = stripName;
        strip.sensitiveVolumeID = sv_i;
        strip.minX = std::min(limitsWorld[0], limitsWorld2[0]);
        strip.maxX = std::max(limitsWorld[0], limitsWorld2[0]);
        strip.minY = std::min(limitsWorld[1], limitsWorld2[1]);
        strip.maxY = std::max(limitsWorld[1], limitsWorld2[1]);
        strip.minZ = std::min(limitsWorld[2], limitsWorld2[2]);
        strip.maxZ = std::max(limitsWorld[2], limitsWorld2[2]);
        strip.normal = auxDetSensitive.GetNormalVector();
        strip.width = halfHeight * 2.;
        strip.null = false;
        strip.module = moduleName;

        // Create a sipm geometry object
        uint32_t channel0 = 32 * ad_i + 2 * sv_i + 0;
        uint32_t channel1 = 32 * ad_i + 2 * sv_i + 1;
        // Sipm0 is on the left in local coords
        double sipm0Y = -halfHeight;
        double sipm1Y = halfHeight;
        // In local coordinates the X position is at half width (top if top) (bottom if not)
        double sipmX = halfWidth;
        if(!fModules[moduleName].top) sipmX = - halfWidth;
        double sipm0XYZ[3] = {sipmX, sipm0Y, 0};
        double sipm0XYZWorld[3];
        auxDetSensitive.LocalToWorld(sipm0XYZ, sipm0XYZWorld);

        CRTSipmGeo sipm0 = {};
        sipm0.channel = channel0;
        sipm0.x = sipm0XYZWorld[0];
        sipm0.y = sipm0XYZWorld[1];
        sipm0.z = sipm0XYZWorld[2];
        sipm0.strip = stripName;
        sipm0.null = false;
        fSipms[channel0] = sipm0;

        double sipm1XYZ[3] = {sipmX, sipm1Y, 0};
        double sipm1XYZWorld[3];
        auxDetSensitive.LocalToWorld(sipm1XYZ, sipm1XYZWorld);

        CRTSipmGeo sipm1 = {};
        sipm1.channel = channel1;
        sipm1.x = sipm1XYZWorld[0];
        sipm1.y = sipm1XYZWorld[1];
        sipm1.z = sipm1XYZWorld[2];
        sipm1.strip = stripName;
        sipm1.null = false;
        fSipms[channel1] = sipm1;

        strip.sipms = std::make_pair(channel0, channel1);
        fStrips[stripName] = strip;
        // Add the strip to the relevant module
        fModules[moduleName].strips[stripName] = strip;
      }
      sv_i++;
    }
    ad_i++;
  }

  // Need to fill the tagger modules after all the strip associations have been made
  // Loop over the modules
  for(auto const& module : fModules){
    // Fill the tagger map
    std::string taggerName = module.second.tagger;
    std::string moduleName = module.second.name;
    fTaggers[taggerName].modules[moduleName] = module.second;
  }

}


CRTGeoAlg::~CRTGeoAlg(){

}

// ----------------------------------------------------------------------------------
// Return the volume enclosed by the whole CRT system
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

// ----------------------------------------------------------------------------------
// Get the number of taggers in the geometry
size_t CRTGeoAlg::NumTaggers() const{
  return fTaggers.size();
}


// ----------------------------------------------------------------------------------
// Get the total number of modules in the geometry
size_t CRTGeoAlg::NumModules() const{
  return fModules.size();
}


// ----------------------------------------------------------------------------------
// Get the total number of strips in the geometry
size_t CRTGeoAlg::NumStrips() const{
  return fStrips.size();
}


// ----------------------------------------------------------------------------------
// Get the tagger geometry object by name
CRTTaggerGeo CRTGeoAlg::GetTagger(std::string taggerName) const{
  for(auto const& tagger : fTaggers){
    if(taggerName == tagger.first) return tagger.second;
  }
  CRTTaggerGeo nullTagger = {};
  nullTagger.null = true;
  return nullTagger;
}

// Get the tagger geometry object by index
CRTTaggerGeo CRTGeoAlg::GetTagger(size_t tagger_i) const{
  size_t index = 0;
  for(auto const& tagger : fTaggers){
    if(tagger_i == index) return tagger.second;
    index++;
  }
  CRTTaggerGeo nullTagger = {};
  nullTagger.null = true;
  return nullTagger;
}


// ----------------------------------------------------------------------------------
// Get the module geometry object by name
CRTModuleGeo CRTGeoAlg::GetModule(std::string moduleName) const{
  for(auto const& module : fModules){
    if(moduleName == module.first) return module.second;
  }
  CRTModuleGeo nullModule = {};
  nullModule.null = true;
  return nullModule;
}

// Get the module geometry object by global index
CRTModuleGeo CRTGeoAlg::GetModule(size_t module_i) const{
  size_t index = 0;
  for(auto const& module : fModules){
    if(module_i == index) return module.second;
    index++;
  }
  CRTModuleGeo nullModule = {};
  nullModule.null = true;
  return nullModule;
}



// ----------------------------------------------------------------------------------
// Get the strip geometry object by name
CRTStripGeo CRTGeoAlg::GetStrip(std::string stripName) const{
  for(auto const& strip : fStrips){
    if(stripName == strip.first) return strip.second;
  }
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;
  return nullStrip;
}

// Get the strip geometry object by global index
CRTStripGeo CRTGeoAlg::GetStrip(size_t strip_i) const{
  size_t index = 0;
  for(auto const& strip : fStrips){
    if(strip_i == index) return strip.second;
    index++;
  }
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;
  return nullStrip;
}


// Get the tagger name from strip or module name
std::string CRTGeoAlg::GetTaggerName(std::string name) const{
  if(fModules.find(name) != fModules.end()){
    return fModules.at(name).tagger;
  }
  if(fStrips.find(name) != fStrips.end()){
    return fModules.at(fStrips.at(name).module).tagger;
  }
  return "";
}

// Get the name of the strip from the SiPM channel ID
std::string CRTGeoAlg::ChannelToStripName(size_t channel) const{
  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      return sipm.second.strip;
    }
  }
  return "";
}


// Recalculate strip limits including charge sharing
std::vector<double> CRTGeoAlg::StripLimitsWithChargeSharing(std::string stripName, double x, double ex){
  int module = fModules.at(fStrips.at(stripName).module).auxDetID;
  std::string moduleName = fGeometryService->AuxDet(module).TotalVolume()->GetName();
  auto const& sensitiveGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(moduleName,
                                                                      2*fStrips.at(stripName).sensitiveVolumeID);

  double halfWidth = sensitiveGeo.HalfWidth1();
  double halfHeight = sensitiveGeo.HalfHeight();
  double halfLength = sensitiveGeo.HalfLength();

  // Note that for the purpose of this reconstruction the "height" coordinates (in y)
  // is the dimension we would conventionally think of as width (distance between the SiPMs)

  // Get the maximum strip limits in world coordinates
  double l1[3] = {halfWidth, -halfHeight + x + ex, halfLength};
  double w1[3];
  sensitiveGeo.LocalToWorld(l1, w1);

  // Get the minimum strip limits in world coordinates
  double l2[3] = {-halfWidth, -halfHeight + x - ex, -halfLength};
  double w2[3];
  sensitiveGeo.LocalToWorld(l2, w2);

  // Use this to get the limits in the two variable directions
  std::vector<double> limits = {std::min(w1[0],w2[0]), std::max(w1[0],w2[0]),
                                std::min(w1[1],w2[1]), std::max(w1[1],w2[1]),
                                std::min(w1[2],w2[2]), std::max(w1[2],w2[2])};
  return limits;
}

// Get the world position of Sipm from the channel ID
geo::Point_t CRTGeoAlg::ChannelToSipmPosition(size_t channel) const{
  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      geo::Point_t position {sipm.second.x, sipm.second.y, sipm.second.z};
      return position;
    }
  }
  geo::Point_t null {-99999, -99999, -99999};
  return null;
}

// Get the sipm channels on a strip
std::pair<int, int> CRTGeoAlg::GetStripSipmChannels(std::string stripName) const{
  for(auto const& strip : fStrips){
    if(stripName == strip.first) return strip.second.sipms;
  }
  return std::make_pair(-99999, -99999);
}

// Return the distance to a sipm in the plane of the sipms
double CRTGeoAlg::DistanceBetweenSipms(geo::Point_t position, size_t channel) const{
  double distance = -99999;

  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      geo::Point_t pos {sipm.second.x, sipm.second.y, sipm.second.z};
      // Get the other sipm
      int otherChannel = channel + 1;
      if(channel % 2) otherChannel = channel - 1;
      // Work out which coordinate is different
      if(fSipms.at(otherChannel).x != pos.X()) distance = position.X() - pos.X();
      if(fSipms.at(otherChannel).y != pos.Y()) distance = position.Y() - pos.Y();
      if(fSipms.at(otherChannel).z != pos.Z()) distance = position.Z() - pos.Z();
      // Return distance in that coordinate
      return distance;
    }
  }
  return distance;
}

// Returns max distance from sipms in strip
double CRTGeoAlg::DistanceBetweenSipms(geo::Point_t position, std::string stripName) const{
  std::pair<int, int> sipms = GetStripSipmChannels(stripName);
  double sipmDist = std::max(DistanceBetweenSipms(position, sipms.first), DistanceBetweenSipms(position, sipms.second));
  return sipmDist;
}

// Return the distance along the strip (from sipm end)
double CRTGeoAlg::DistanceDownStrip(geo::Point_t position, std::string stripName) const{
  double distance = -99999;
  for(auto const& strip : fStrips){
    if(stripName == strip.first){
      geo::Point_t pos = ChannelToSipmPosition(strip.second.sipms.first);
      // Work out the longest dimension of strip
      double xdiff = std::abs(strip.second.maxX-strip.second.minX);
      double ydiff = std::abs(strip.second.maxY-strip.second.minY);
      double zdiff = std::abs(strip.second.maxZ-strip.second.minZ);
      if(xdiff > ydiff && xdiff > zdiff) distance = position.X() - pos.X();
      if(ydiff > xdiff && ydiff > zdiff) distance = position.Y() - pos.Y();
      if(zdiff > xdiff && zdiff > ydiff) distance = position.Z() - pos.Z();
      return std::abs(distance);
    }
  }
  return distance;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside CRT volume
bool CRTGeoAlg::IsInsideCRT(TVector3 point){
  geo::Point_t pt {point.X(), point.Y(), point.Z()};
  return IsInsideCRT(pt);
}

bool CRTGeoAlg::IsInsideCRT(geo::Point_t point){
  std::vector<double> limits = CRTLimits();
  return point.X() > limits[0] && point.Y() > limits[1] && point.Z() > limits[2]
    && point.X() < limits[3] && point.Y() < limits[4] && point.Z() < limits[5];
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a tagger
bool CRTGeoAlg::IsInsideTagger(const CRTTaggerGeo& tagger, geo::Point_t point){
  if(tagger.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = tagger.minX;
  double ymin = tagger.minY;
  double zmin = tagger.minZ;
  double xmax = tagger.maxX;
  double ymax = tagger.maxY;
  double zmax = tagger.maxZ;

  return x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a module
bool CRTGeoAlg::IsInsideModule(const CRTModuleGeo& module, geo::Point_t point){
  if(module.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = module.minX;
  double ymin = module.minY;
  double zmin = module.minZ;
  double xmax = module.maxX;
  double ymax = module.maxY;
  double zmax = module.maxZ;

  return x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a strip
bool CRTGeoAlg::IsInsideStrip(const CRTStripGeo& strip, geo::Point_t point){
  if(strip.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = strip.minX;
  double ymin = strip.minY;
  double zmin = strip.minZ;
  double xmax = strip.maxX;
  double ymax = strip.maxY;
  double zmax = strip.maxZ;

  return x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax;
}

// ----------------------------------------------------------------------------------
// Check if two modules overlap in 2D
bool CRTGeoAlg::CheckOverlap(const CRTModuleGeo& module1, const CRTModuleGeo& module2){
  // Get the minimum and maximum X, Y, Z coordinates
  double minX = std::max(module1.minX, module2.minX);
  double maxX = std::min(module1.maxX, module2.maxX);
  double minY = std::max(module1.minY, module2.minY);
  double maxY = std::min(module1.maxY, module2.maxY);
  double minZ = std::max(module1.minZ, module2.minZ);
  double maxZ = std::min(module1.maxZ, module2.maxZ);

  // If the two strips overlap in 2 dimensions then return true
  return (minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ);
}

// ----------------------------------------------------------------------------------
// Check is a module overlaps with a perpendicual module in the same tagger
bool CRTGeoAlg::HasOverlap(const CRTModuleGeo& module){
  // Record plane of mother module
  size_t planeID = module.planeID;
  // Get mother tagger of module
  std::string taggerName = module.tagger;
  // Loop over other modules in tagger
  for(auto const& module2 : fTaggers[taggerName].modules){
    // If in other plane loop over strips
    if(module2.second.planeID == planeID) continue;
    // Check for overlaps
    if(CheckOverlap(module, module2.second)) return true;
  }
  return false;
}

bool CRTGeoAlg::StripHasOverlap(std::string stripName){
  return HasOverlap(fModules.at(fStrips.at(stripName).module));
}

// ----------------------------------------------------------------------------------
// Find the average of the tagger entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::TaggerCrossingPoint(std::string taggerName, const simb::MCParticle& particle){
  const CRTTaggerGeo& tagger = fTaggers.at(taggerName);
  return TaggerCrossingPoint(tagger, particle);
}

geo::Point_t CRTGeoAlg::TaggerCrossingPoint(const CRTTaggerGeo& tagger, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideTagger(tagger, point)) continue;
    if(first){
      entry = point;
      first = false;
    }
    exit = point;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesTagger(const CRTTaggerGeo& tagger, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideTagger(tagger, point)) return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Find the average of the module entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::ModuleCrossingPoint(std::string moduleName, const simb::MCParticle& particle){
  const CRTModuleGeo& module = fModules.at(moduleName);
  return ModuleCrossingPoint(module, particle);
}

geo::Point_t CRTGeoAlg::ModuleCrossingPoint(const CRTModuleGeo& module, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideModule(module, point)){
      // Look at the mid point
      if(i == particle.NumberTrajectoryPoints()-1) continue;
      geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
      geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
      if(!IsInsideModule(module, mid)) continue;
      point = mid;
    }
    if(first){
      entry = point;
      first = false;
    }
    exit = point;

    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(!IsInsideModule(module, mid)) continue;
    exit = mid;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesModule(const CRTModuleGeo& module, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideModule(module, point)) return true;

    // Also look at midpoint with next point
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(IsInsideModule(module, mid)) return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Find the average of the strip entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::StripCrossingPoint(std::string stripName, const simb::MCParticle& particle){
  const CRTStripGeo& strip = fStrips.at(stripName);
  return StripCrossingPoint(strip, particle);
}

geo::Point_t CRTGeoAlg::StripCrossingPoint(const CRTStripGeo& strip, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideStrip(strip, point)){
      if(i == particle.NumberTrajectoryPoints()-1) continue;
      geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
      geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
      if(!IsInsideStrip(strip, mid)) continue;
      point = mid;
    }
    if(first){
      entry = point;
      first = false;
    }
    exit = point;
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(!IsInsideStrip(strip, mid)) continue;
    exit = mid;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesStrip(const CRTStripGeo& strip, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideStrip(strip, point)) return true;

    // Also look at midpoint with next point
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(IsInsideStrip(strip, mid)) return true;
  }
  return false;
}


// ----------------------------------------------------------------------------------
// Work out which strips the true particle crosses
std::vector<std::string> CRTGeoAlg::CrossesStrips(const simb::MCParticle& particle){
  std::vector<std::string> stripNames;
  for(auto const& tagger : fTaggers){
    if(!CrossesTagger(tagger.second, particle)) continue;
    for(auto const& module : tagger.second.modules){
      if(!CrossesModule(module.second, particle)) continue;
      for(auto const& strip : module.second.strips){
        if(!CrossesStrip(strip.second, particle)) continue;
        if(std::find(stripNames.begin(), stripNames.end(), strip.first) != stripNames.end()) continue;
        stripNames.push_back(strip.first);
      }
    }
  }
  return stripNames;
}


// ----------------------------------------------------------------------------------
// Find the angle of true particle trajectory to tagger
double CRTGeoAlg::AngleToTagger(std::string taggerName, const simb::MCParticle& particle){
  // Get normal to tagger using the top modules
  TVector3 normal (0,0,0);
  for(auto const& module : GetTagger(taggerName).modules){
    normal.SetXYZ(module.second.normal.X(), module.second.normal.Y(), module.second.normal.Z());
    break;
  }
  //FIXME this is pretty horrible
  if(normal.X()<0.5 && normal.X()>-0.5) normal.SetX(0);
  if(normal.Y()<0.5 && normal.Y()>-0.5) normal.SetY(0);
  if(normal.Z()<0.5 && normal.Z()>-0.5) normal.SetZ(0);

  if(std::abs(normal.X())==1 && GetTagger(taggerName).minX < 0) normal.SetX(-1);
  if(std::abs(normal.X())==1 && GetTagger(taggerName).minX > 0) normal.SetX(1);
  if(std::abs(normal.Y())==1 && GetTagger(taggerName).minY < 0) normal.SetY(-1);
  if(std::abs(normal.Y())==1 && GetTagger(taggerName).minY > 0) normal.SetY(1);
  if(std::abs(normal.Z())==1 && GetTagger(taggerName).minZ < 0) normal.SetZ(-1);
  if(std::abs(normal.Z())==1 && GetTagger(taggerName).minZ > 0) normal.SetZ(1);

  TVector3 start (particle.Vx(), particle.Vy(), particle.Vz());
  TVector3 end (particle.EndX(), particle.EndY(), particle.EndZ());
  TVector3 diff = end - start;

  if(normal.X() == 1 && start.X() > end.X()) diff = start - end;
  if(normal.X() == -1 && start.X() < end.X()) diff = start - end;
  if(normal.Y() == 1 && start.Y() > end.Y()) diff = start - end;
  if(normal.Y() == -1 && start.Y() < end.Y()) diff = start - end;
  if(normal.Z() == 1 && start.Z() > end.Z()) diff = start - end;
  if(normal.Z() == -1 && start.Z() < end.Z()) diff = start - end;

  return normal.Angle(diff);
}


// ----------------------------------------------------------------------------------
// Check if a particle enters the CRT volume
bool CRTGeoAlg::EntersVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  std::vector<double> limits = CRTLimits();
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideCRT(point)){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  return enters && (startOutside || endOutside);
}

// ----------------------------------------------------------------------------------
// Check if a particle crosses the CRT volume
bool CRTGeoAlg::CrossesVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  std::vector<double> limits = CRTLimits();
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideCRT(point)){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  return startOutside && enters && endOutside;
}


// ----------------------------------------------------------------------------------
// Determine if a particle would be able to produce a hit in a tagger
bool CRTGeoAlg::ValidCrossingPoint(std::string taggerName, const simb::MCParticle& particle){

  // Get all the crossed strips in the tagger
  std::vector<std::string> crossedModules;
  for(auto const& module : fTaggers[taggerName].modules){
    geo::Point_t crossPoint = ModuleCrossingPoint(module.second.name, particle);
    if(crossPoint.X() != -99999) crossedModules.push_back(module.second.name);
  }

  // Check if the strip has a possible overlap, return true if not
  for(size_t i = 0; i < crossedModules.size(); i++){
    if(!HasOverlap(fModules[crossedModules[i]])) return true;
    // Check if any of the crossed strips overlap, return true if they do
    for(size_t j = i; j < crossedModules.size(); j++){
      if(CheckOverlap(fModules[crossedModules[i]], fModules[crossedModules[j]])) return true;
    }
  }
  return false;
}

}
