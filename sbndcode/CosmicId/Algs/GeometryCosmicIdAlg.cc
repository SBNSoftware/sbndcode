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

  fPeLimit = config.PeLimit();
  return;
}

// Remove any tracks in different TPC to beam activity
bool GeometryCosmicIdAlg::GeometryCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, bool tpc0Flash, bool tpc1Flash){

  // Remove any tracks that are detected in one TPC and reconstructed in another
  int tpc = fTpcGeo.DetectedInTPC(hits);
  double startX = track.Start().X();
  double endX = track.End().X();

  // Check if track is stitched
  if(tpc == -1) return false;

  // Check the start/end points are in same TPC (track shifted into other TPC because time outside of beam)
  if(tpc == 0 && (startX>0 || endX>0)) return true;
  else if(tpc == 1 && (startX<0 || endX<0)) return true;

  // Check if track was detected in TPC where there was no beam activity
  if(tpc0Flash && !tpc1Flash && tpc == 1) return true;
  if(tpc1Flash && !tpc0Flash && tpc == 0) return true;
  
  return false;
  
}

// Remove any tracks in different TPC to beam activity with cut on PE
bool GeometryCosmicIdAlg::GeometryCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, double tpc0_pe, double tpc1_pe){

  // Remove any tracks that are detected in one TPC and reconstructed in another
  int tpc = fTpcGeo.DetectedInTPC(hits);
  double startX = track.Start().X();
  double endX = track.End().X();

  // Check if track is stitched
  if(tpc == -1) return false;

  // Check the start/end points are in same TPC (track shifted into other TPC because time outside of beam)
  if(tpc == 0 && (startX>0 || endX>0)) return true;
  else if(tpc == 1 && (startX<0 || endX<0)) return true;

  // Check if track was detected in TPC where there was no beam activity
  if(tpc == 1 && tpc1_pe < fPeLimit) return true;
  if(tpc == 0 && tpc0_pe < fPeLimit) return true;
  
  return false;
  
}


}
