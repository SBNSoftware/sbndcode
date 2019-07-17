#include "GeometryCosmicIdAlg.h"

namespace sbnd{

GeometryCosmicIdAlg::GeometryCosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


GeometryCosmicIdAlg::GeometryCosmicIdAlg(){

}


GeometryCosmicIdAlg::~GeometryCosmicIdAlg(){

}


void GeometryCosmicIdAlg::reconfigure(const Config& config){

  return;
}

bool GeometryCosmicIdAlg::GeometryCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, bool tpc0Flash, bool tpc1Flash){

  // Remove any tracks that are detected in one TPC and reconstructed in another
  int tpc = fTpcGeo.DetectedInTPC(hits);
  double startX = track.Start().X();
  double endX = track.End().X();

  // Check if track is stitched
  // If it isn't check the start/end points are in same TPC
  if(tpc == 0 && (startX>0 || endX>0)) return true;
  else if(tpc == 1 && (startX<0 || endX<0)) return true;

  if(tpc0Flash && !tpc1Flash && tpc != 0) return true;
  if(tpc1Flash && !tpc0Flash && tpc != 1) return true;
  
  return false;
  
}


}
