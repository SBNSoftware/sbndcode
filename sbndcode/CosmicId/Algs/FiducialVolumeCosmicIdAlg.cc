#include "FiducialVolumeCosmicIdAlg.h"

namespace sbnd{

FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicIdAlg(){

}


FiducialVolumeCosmicIdAlg::~FiducialVolumeCosmicIdAlg(){

}


void FiducialVolumeCosmicIdAlg::reconfigure(const Config& config){

  fMinX = config.FiducialCuts().MinX(); 
  fMinY = config.FiducialCuts().MinY(); 
  fMinZ = config.FiducialCuts().MinZ(); 
  fMaxX = config.FiducialCuts().MaxX(); 
  fMaxY = config.FiducialCuts().MaxY(); 
  fMaxZ = config.FiducialCuts().MaxZ();

  return;
}

bool FiducialVolumeCosmicIdAlg::InFiducial(geo::Point_t point){

  return fTpcGeo.InFiducial(point, fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

}

bool FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicId(recob::Track track){
  
  bool startInFiducial = InFiducial(track.Vertex());

  bool endInFiducial = InFiducial(track.End());

  if(!startInFiducial && !endInFiducial)  return true;
  
  return false;
  
}


}
