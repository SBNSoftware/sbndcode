#include "ApaCrossCosmicIdAlg.h"

namespace sbnd{

ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg(const Config& config){

  this->reconfigure(config);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

}


ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg(){

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

}


ApaCrossCosmicIdAlg::~ApaCrossCosmicIdAlg(){

}


void ApaCrossCosmicIdAlg::reconfigure(const Config& config){

  fDistanceLimit = config.DistanceLimit(); 
  fMaxApaDistance = config.MaxApaDistance();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}

std::pair<double, double> ApaCrossCosmicIdAlg::MinApaDistance(recob::Track track, std::vector<double> t0List, int tpc){

  double crossTime = -99999;
  double xmax = fTpcGeo.MaxX();

  double minDist = 99999;
  double startX = track.Vertex().X();
  double endX = track.End().X();
  geo::Point_t point = track.Vertex();

  // Don't try to shift tracks near the Apa
  if(std::abs(startX) > xmax - fMaxApaDistance 
     || std::abs(endX) > xmax - fMaxApaDistance) return std::make_pair(minDist, crossTime);

  // If in tpc 0 use start/end with lowest X
  if(tpc == 0 && endX < startX) point = track.End();

  // If in tpc 1 use start/end with highest X
  if(tpc == 1 && endX > startX) point = track.End();

  //Shift track by all t0's
  for(auto const& t0 : t0List){
    // If particle crosses the APA before t = 0 the crossing point won't be reconstructed
    if(t0 < 0) continue;
    double shiftedX = point.X();
    double shift = t0 * fDetectorProperties->DriftVelocity();
    if(tpc == 0) shiftedX = point.X() - shift;
    if(tpc == 1) shiftedX = point.X() + shift;

    //Check track still in TPC
    if(std::abs(shiftedX) > (xmax + fDistanceLimit)) continue;
    //Calculate distance between start/end and APA
    double dist = std::abs(std::abs(shiftedX) - xmax);
    if(dist < minDist) {
      minDist = dist;
      crossTime = t0;
    }
  }

  return std::make_pair(minDist, crossTime);

}

double ApaCrossCosmicIdAlg::T0FromApaCross(recob::Track track, std::vector<double> t0List, int tpc){

  std::pair<double, double> min = MinApaDistance(track, t0List, tpc);
  if(min.first < fDistanceLimit) return min.second;
  return -99999;

}


double ApaCrossCosmicIdAlg::ApaDistance(recob::Track track, double t0, std::vector<art::Ptr<recob::Hit>> hits){

  std::vector<double> t0List {t0};
  int tpc = fTpcGeo.DetectedInTPC(hits);
  std::pair<double, double> min = MinApaDistance(track, t0List, tpc);
  return min.first;

}

std::pair<double, double> ApaCrossCosmicIdAlg::MinApaDistance(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  int tpc = fTpcGeo.DetectedInTPC(hits);
  if(tpc == 0){
    return MinApaDistance(track, t0Tpc0, tpc);
  }
  if(tpc == 1){
    return MinApaDistance(track, t0Tpc1, tpc);
  }
  return std::make_pair(-99999, -99999);
} 


bool ApaCrossCosmicIdAlg::ApaCrossCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  int tpc = fTpcGeo.DetectedInTPC(hits);

  if(tpc == 0){
    double crossTimeThrough = T0FromApaCross(track, t0Tpc0, tpc);
    if(crossTimeThrough != -99999 && (crossTimeThrough < fBeamTimeMin || crossTimeThrough > fBeamTimeMax)) return true;
  }
  if(tpc == 1){
    double crossTimeThrough = T0FromApaCross(track, t0Tpc1, tpc);
    if(crossTimeThrough != -99999 && (crossTimeThrough < fBeamTimeMin || crossTimeThrough > fBeamTimeMax)) return true;
  }

  return false;

}

}
