#include "ApaCrossCosmicIdAlg.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

namespace sbnd{

ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg(){

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


// Get the minimum distance from track to APA for different times
  std::pair<double, double> ApaCrossCosmicIdAlg::MinApaDistance(detinfo::DetectorPropertiesData const& detProp,
                                                                recob::Track track, std::vector<double> t0List, int tpc){

  double crossTime = -99999;
  double xmax = fTpcGeo.MaxX();

  double minDist = 99999;
  double startX = track.Vertex().X();
  double endX = track.End().X();
  geo::Point_t point = track.Vertex();

  // Don't try to shift tracks near the APA (may give artificially small distances)
  if(std::abs(startX) > xmax - fMaxApaDistance 
     || std::abs(endX) > xmax - fMaxApaDistance) return std::make_pair(minDist, crossTime);

  // If in TPC 0 use start/end with lowest X
  if(tpc == 0 && endX < startX) point = track.End();

  // If in TPC 1 use start/end with highest X
  if(tpc == 1 && endX > startX) point = track.End();

  // If in both TPCs (stitched) return null values
  if(tpc == -1) return std::make_pair(minDist, crossTime);

  // Shift track by all t0's
  for(auto const& t0 : t0List){
    // If particle crosses the APA before t = 0 the crossing point won't be reconstructed
    if(t0 < 0) continue;
    double shiftedX = point.X();
    double shift = t0 * detProp.DriftVelocity();
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


// Get time by matching tracks which cross the APA
double ApaCrossCosmicIdAlg::T0FromApaCross(detinfo::DetectorPropertiesData const& detProp,
                                           recob::Track track, std::vector<double> t0List, int tpc){

  // Get the minimum distance to the APA and corresponding time
  std::pair<double, double> min = MinApaDistance(detProp, track, t0List, tpc);
  // Check the distance is within allowed limit
  if(min.first < fDistanceLimit) return min.second;
  return -99999;

}


// Get the distance from track to APA at fixed time
double ApaCrossCosmicIdAlg::ApaDistance(detinfo::DetectorPropertiesData const& detProp,
                                        recob::Track track, double t0, std::vector<art::Ptr<recob::Hit>> hits){

  std::vector<double> t0List {t0};
  // Determine the TPC from hit collection
  int tpc = fTpcGeo.DetectedInTPC(hits);
  // Get the distance to the APA at the given time
  std::pair<double, double> min = MinApaDistance(detProp, track, t0List, tpc);
  return min.first;

}

// Work out what TPC track is in and get the minimum distance from track to APA for different times
std::pair<double, double> ApaCrossCosmicIdAlg::MinApaDistance(detinfo::DetectorPropertiesData const& detProp,
                                                              recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  // Determine the TPC from hit collection
  int tpc = fTpcGeo.DetectedInTPC(hits);
  // Get the minimum distance from the APA in the corresponding TPC (with corresponding flashes)
  if(tpc == 0){
    return MinApaDistance(detProp, track, t0Tpc0, tpc);
  }
  if(tpc == 1){
    return MinApaDistance(detProp, track, t0Tpc1, tpc);
  }
  return std::make_pair(-99999, -99999);
} 


// Tag tracks with times outside the beam
bool ApaCrossCosmicIdAlg::ApaCrossCosmicId(detinfo::DetectorPropertiesData const& detProp,
                                           recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  // Determine the TPC from hit collection
  int tpc = fTpcGeo.DetectedInTPC(hits);

  // Get the minimum distance from the APA in the corresponding TPC (with corresponding flashes)
  if(tpc == 0){
    // Get the closest matched time with distance below limit
    double crossTimeThrough = T0FromApaCross(detProp, track, t0Tpc0, tpc);
    // If the matched time is outside of the beam time then tag as a cosmic
    if(crossTimeThrough != -99999 && (crossTimeThrough < fBeamTimeMin || crossTimeThrough > fBeamTimeMax)) return true;
  }
  if(tpc == 1){
    double crossTimeThrough = T0FromApaCross(detProp, track, t0Tpc1, tpc);
    if(crossTimeThrough != -99999 && (crossTimeThrough < fBeamTimeMin || crossTimeThrough > fBeamTimeMax)) return true;
  }

  return false;

}

}
