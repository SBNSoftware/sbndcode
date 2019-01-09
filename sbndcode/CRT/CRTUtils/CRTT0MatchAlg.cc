#include "CRTT0MatchAlg.h"

namespace sbnd{

CRTT0MatchAlg::CRTT0MatchAlg(const Config& config){

  this->reconfigure(config);
  
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

}


CRTT0MatchAlg::CRTT0MatchAlg(){

  fGeometryService = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

}


CRTT0MatchAlg::~CRTT0MatchAlg(){

}


void CRTT0MatchAlg::reconfigure(const Config& config){

  return;

}
 

// Utility function that determines the possible t0 range of a track
std::pair<double, double> CRTT0MatchAlg::TrackT0Range(double startX, double endX, int tpc){

  std::pair<double, double> result;
  double Vd = fDetectorProperties->DriftVelocity();

  // Whole track must be within tpc
  // Find which TPC the track hits are detected in
  if(tpc == 0){
    // Lowest |X| is furthest from APA
    double lowX = std::max(startX, endX);
    // xmin is shift from furthest to 0 (the CPA)
    double xmax = 0 - lowX;
    // Highest |X| is closest to APA
    double highX = std::min(startX, endX);
    // xmax is shift from closest to APA
    double xmin = -(2.0*fGeometryService->DetHalfWidth()+3.) - highX;
    double t0max = -(xmin/Vd);
    double t0min = -(xmax/Vd);
    result = std::make_pair(t0min, t0max);
  }

  else{
    // Lowest |X| is furthest from APA
    double lowX = std::min(startX, endX);
    // xmin is shift from furthest to 0 (the CPA)
    double xmin = 0 - lowX;
    // Highest |X| is closest to APA
    double highX = std::max(startX, endX);
    // xmax is shift from closest to APA
    double xmax = (2.0*fGeometryService->DetHalfWidth()+3.) - highX;
    double t0min = xmin/Vd; 
    double t0max = xmax/Vd; 
    result = std::make_pair(t0min, t0max);
  }

  return result;

} // CRTT0MatchAlg::TrackT0Range()


double CRTT0MatchAlg::DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, crt::CRTHit crtHit, int tpc, double t0){

  double minDist = 99999;

  // Convert the t0 into an x shift
  double shift = t0 * fDetectorProperties->DriftVelocity();
  // Apply the shift depending on which TPC the track is in
  if (tpc == 1) trackPos[0] += shift;
  if (tpc == 0) trackPos[0] -= shift;

  TVector3 endPos = trackPos + trackDir;
  double denominator = (endPos - trackPos).Mag();
  // 1D hits should only have a lot of variance in one direction
  if(crtHit.x_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double xpos = crtHit.x_pos + ((i+1.)/10. - 1.)*crtHit.x_err;
      TVector3 crtPoint(xpos, crtHit.y_pos, crtHit.z_pos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else if(crtHit.y_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double ypos = crtHit.y_pos + ((i+1.)/10. - 1.)*crtHit.y_err;
      TVector3 crtPoint(crtHit.x_pos, ypos, crtHit.z_pos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else if(crtHit.y_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double zpos = crtHit.z_pos + ((i+1.)/10. - 1.)*crtHit.z_err;
      TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, zpos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else{
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
    double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
    double dca = numerator/denominator;
    if(dca < minDist) minDist = dca;
  }

  return minDist;

} // CRTT0MatchAlg::DistToOfClosestApproach()


}
