#include "CRTHitRecoAlg.h"

namespace sbnd{

CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){

  this->reconfigure(config);
  
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

}


CRTHitRecoAlg::CRTHitRecoAlg(){

  fGeometryService = lar::providerFrom<geo::Geometry>();
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

}


CRTHitRecoAlg::~CRTHitRecoAlg(){

}


void CRTHitRecoAlg::reconfigure(const Config& config){

  return;
}
  

// Function to calculate the strip position limits in real space from channel
std::vector<double> CRTHitRecoAlg::ChannelToLimits(CRTStrip stripHit){

  // Get strip geometry from the channel ID
  int strip = (stripHit.channel >> 1) & 15;
  int module = (stripHit.channel >> 5);
  std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
  const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);

  double halfWidth = stripGeo.HalfWidth1();
  double halfHeight = stripGeo.HalfHeight();
  double halfLength = stripGeo.HalfLength();

  // Get the maximum strip limits in world coordinates
  double l1[3] = {-halfWidth+stripHit.x+stripHit.ex, halfHeight, halfLength};
  double w1[3] = {0,0,0};
  stripGeo.LocalToWorld(l1, w1);

  // Get the minimum strip limits in world coordinates
  double l2[3] = {-halfWidth+stripHit.x-stripHit.ex, -halfHeight, -halfLength};
  double w2[3] = {0,0,0};
  stripGeo.LocalToWorld(l2, w2);

  // Use this to get the limits in the two variable directions
  std::vector<double> limits = {std::min(w1[0],w2[0]), std::max(w1[0],w2[0]), 
                                std::min(w1[1],w2[1]), std::max(w1[1],w2[1]), 
                                std::min(w1[2],w2[2]), std::max(w1[2],w2[2])};
  return limits;

} // CRTHitRecoAlg::ChannelToLimits()


// Function to calculate the overlap between two crt strips
std::vector<double> CRTHitRecoAlg::CrtOverlap(std::vector<double> strip1, std::vector<double> strip2){

  // Get the minimum and maximum X, Y, Z coordinates
  double minX = std::max(strip1[0], strip2[0]);
  double maxX = std::min(strip1[1], strip2[1]);
  double minY = std::max(strip1[2], strip2[2]);
  double maxY = std::min(strip1[3], strip2[3]);
  double minZ = std::max(strip1[4], strip2[4]);
  double maxZ = std::min(strip1[5], strip2[5]);

  std::vector<double> null = {-99999, -99999, -99999, -99999, -99999, -99999};
  std::vector<double> overlap = {minX, maxX, minY, maxY, minZ, maxZ};

  // If the two strips overlap in 2 dimensions then return the overlap
  if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return overlap;
  // Otherwise return a "null" value
  return null;

} // CRTHitRecoAlg::CRTOverlap()


// Function to return the CRT tagger name and module position from the channel ID
std::pair<std::string,unsigned> CRTHitRecoAlg::ChannelToTagger(uint32_t channel){

  // Get the strip geometry from the channel ID
  int strip = (channel >> 1) & 15;
  int module = (channel >> 5);
  std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
  TVector3 center = fAuxDetGeoCore->AuxDetChannelToPosition(2*strip, name);
  const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);

  // Get the full volume path string
  std::set<std::string> volNames = {stripGeo.TotalVolume()->GetName()};
  std::vector<std::vector<TGeoNode const*> > paths = fGeometryService->FindAllVolumePaths(volNames);
  std::string path = "";
  for (size_t inode=0; inode<paths.at(0).size(); inode++) {
    path += paths.at(0).at(inode)->GetName();
    if (inode < paths.at(0).size() - 1) {
      path += "/";
    }
  }

  // Retrive the geometry manager from the path
  TGeoManager* manager = fGeometryService->ROOTGeoManager();
  manager->cd(path.c_str());

  // Get the parent module and tagger
  TGeoNode* nodeModule = manager->GetMother(2);
  TGeoNode* nodeTagger = manager->GetMother(3);

  // Module position in parent (tagger) frame
  double origin[3] = {0, 0, 0};
  double modulePosMother[3];
  nodeModule->LocalToMaster(origin, modulePosMother);
  unsigned planeID = (modulePosMother[2] > 0);
  // Get the name of the tagger
  std::string tagName = nodeTagger->GetName();
  std::pair<std::string, unsigned> output = std::make_pair(tagName, planeID);

  return output;

} // CRTHitRecoAlg::ChannelToTagger()


// Function to check if a CRT strip overlaps with a perpendicular module
bool CRTHitRecoAlg::CheckModuleOverlap(uint32_t channel){

  // FIXME: Would be better to check overlap of all individual strips rather than modules
  bool hasOverlap = false;

  // Get the strip geometry from the channel ID
  int strip = (channel >> 1) & 15;
  int module = (channel >> 5);
  std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
  const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);

  // Get the parent module and tagger from the the geometry manager
  std::set<std::string> volNames = {stripGeo.TotalVolume()->GetName()};
  std::vector<std::vector<TGeoNode const*> > paths = fGeometryService->FindAllVolumePaths(volNames);
  std::string path = "";
  for (size_t inode=0; inode<paths.at(0).size(); inode++) {
    path += paths.at(0).at(inode)->GetName();
    if (inode < paths.at(0).size() - 1) {
      path += "/";
    }
  }
  TGeoManager* manager = fGeometryService->ROOTGeoManager();
  manager->cd(path.c_str());
  TGeoNode* nodeModule = manager->GetMother(2);
  TGeoNode* nodeTagger = manager->GetMother(3);
  std::string modName = nodeModule->GetName();

  // Get the limits of the module in the tagger frame
  double height = fGeometryService->AuxDet(module).HalfHeight();
  double width = fGeometryService->AuxDet(module).HalfWidth1();
  double length = fGeometryService->AuxDet(module).Length()/2.;
  double pos1[3] = {width, height, length};
  double tagp1[3];
  nodeModule->LocalToMaster(pos1, tagp1);
  double pos2[3] = {-width, -height, -length};
  double tagp2[3];
  nodeModule->LocalToMaster(pos2, tagp2);
  std::vector<double> limits = {std::min(tagp1[0], tagp2[0]),
                                std::max(tagp1[0], tagp2[0]),
                                std::min(tagp1[1], tagp2[1]),
                                std::max(tagp1[1], tagp2[1]),
                                std::min(tagp1[2], tagp2[2]),
                                std::max(tagp1[2], tagp2[2])};

  // Get which layer the module is in the tagger
  double origin[3] = {0, 0, 0};
  double modulePosMother[3];
  nodeModule->LocalToMaster(origin, modulePosMother);

  unsigned planeID = (modulePosMother[2] > 0);

  // Get the number of daughters from the tagger
  int nDaughters = nodeTagger->GetNdaughters();

  // Loop over the daughters
  for(int mod_i = 0; mod_i < nDaughters; mod_i++){
    // Check the name not the same as the current module
    TGeoNode* nodeDaughter = nodeTagger->GetDaughter(mod_i);
    std::string d_name = nodeDaughter->GetName();
    // Remove last two characters from name to match the AuxDet name
    if(d_name == modName) continue;

    // Get the limits of the module in the tagger frame
    double d_tagp1[3];
    nodeDaughter->LocalToMaster(pos1, d_tagp1);
    double d_tagp2[3];
    nodeDaughter->LocalToMaster(pos2, d_tagp2);
    std::vector<double> d_limits = {std::min(d_tagp1[0], d_tagp2[0]),
                                    std::max(d_tagp1[0], d_tagp2[0]),
                                    std::min(d_tagp1[1], d_tagp2[1]),
                                    std::max(d_tagp1[1], d_tagp2[1]),
                                    std::min(d_tagp1[2], d_tagp2[2]),
                                    std::max(d_tagp1[2], d_tagp2[2])};

    // Get which layer the module is in the tagger
    double d_modulePosMother[3];
    nodeDaughter->LocalToMaster(origin, d_modulePosMother);
    unsigned d_planeID = (d_modulePosMother[2] > 0);

    // Check the overlap of the two modules
    std::vector<double> overlap = CrtOverlap(limits, d_limits);
    // If there is an overlap set to true and the modules are in different layers
    if(overlap[0]!=-99999 && d_planeID!=planeID) hasOverlap = true;
  }

  return hasOverlap;

} // CRTHitRecoAlg::CheckModuleOverlap


// Function to make filling a CRTHit a bit faster
sbnd::crt::CRTHit CRTHitRecoAlg::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                              std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, int plane, 
                              double x, double ex, double y, double ey, double z, double ez, std::string tagger){

  sbnd::crt::CRTHit crtHit;

  crtHit.feb_id      = tfeb_id;
  crtHit.pesmap      = tpesmap;
  crtHit.peshit      = peshit;
  crtHit.ts0_s_corr  = 0; 
  crtHit.ts0_ns      = time * 1e3;
  crtHit.ts0_ns_corr = 0; 
  crtHit.ts1_ns      = time * 1e3;
  crtHit.ts0_s       = time * 1e-6;
  crtHit.x_pos       = x;
  crtHit.x_err       = ex;
  crtHit.y_pos       = y; 
  crtHit.y_err       = ey;
  crtHit.z_pos       = z;
  crtHit.z_err       = ez;
  crtHit.tagger      = tagger;

  return crtHit;

} // CRTHitRecoAlg::FillCrtHit()

}
