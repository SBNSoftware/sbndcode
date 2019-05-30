#include "FiducialVolumeCosmicTagAlg.h"

namespace sbnd{

FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTagAlg(){

}


FiducialVolumeCosmicTagAlg::~FiducialVolumeCosmicTagAlg(){

}


void FiducialVolumeCosmicTagAlg::reconfigure(const Config& config){

  fMinX = config.FiducialCuts().MinX(); 
  fMinY = config.FiducialCuts().MinY(); 
  fMinZ = config.FiducialCuts().MinZ(); 
  fMaxX = config.FiducialCuts().MaxX(); 
  fMaxY = config.FiducialCuts().MaxY(); 
  fMaxZ = config.FiducialCuts().MaxZ();

  return;
}

bool FiducialVolumeCosmicTagAlg::InFiducial(geo::Point_t point){

  return CosmicRemovalUtils::InFiducial(point, fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

}

bool FiducialVolumeCosmicTagAlg::FiducialVolumeCosmicTag(recob::Track track){
  
  bool startInFiducial = InFiducial(track.Vertex());

  bool endInFiducial = InFiducial(track.End());

  if(!startInFiducial && !endInFiducial)  return true;
  
  return false;
  
}


}
