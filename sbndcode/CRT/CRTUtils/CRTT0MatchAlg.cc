#include "CRTT0MatchAlg.h"

namespace sbnd{

CRTT0MatchAlg::CRTT0MatchAlg(const Config& config){

  this->reconfigure(config);
  
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

}


CRTT0MatchAlg::CRTT0MatchAlg(){

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

}


CRTT0MatchAlg::~CRTT0MatchAlg(){

}


void CRTT0MatchAlg::reconfigure(const Config& config){

  fMinTrackLength = config.MinTrackLength();
  fTrackDirectionFrac = config.TrackDirectionFrac();
  fDistanceLimit = config.DistanceLimit();

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
    // xmin is shift from furthest to 0 (the CPA) //FIXME CPA width
    double xmax = 0 - lowX;
    // Highest |X| is closest to APA
    double highX = std::min(startX, endX);
    // xmax is shift from closest to APA
    double xmin = fTpcGeo.MinX() - highX;
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
    double xmax = fTpcGeo.MaxX() - highX;
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

std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirectionAverage(recob::Track track, double frac){

  // Calculate direction as an average over directions
  size_t nTrackPoints = track.NumberTrajectoryPoints();
  recob::TrackTrajectory trajectory  = track.Trajectory();
  std::vector<geo::Vector_t> validDirections;
  for(size_t i = 0; i < nTrackPoints; i++){
    if(trajectory.FlagsAtPoint(i)!=recob::TrajectoryPointFlags::InvalidHitIndex) continue;
    validDirections.push_back(track.DirectionAtPoint(i));
  }

  size_t nValidPoints = validDirections.size();
  int endPoint = (int)floor(nValidPoints*frac);
  double xTotStart = 0; double yTotStart = 0; double zTotStart = 0;
  double xTotEnd = 0; double yTotEnd = 0; double zTotEnd = 0;
  for(int i = 0; i < endPoint; i++){
    geo::Vector_t dirStart = validDirections.at(i);
    geo::Vector_t dirEnd = validDirections.at(nValidPoints - (i+1));
    xTotStart += dirStart.X();
    yTotStart += dirStart.Y();
    zTotStart += dirStart.Z();
    xTotEnd += dirEnd.X();
    yTotEnd += dirEnd.Y();
    zTotEnd += dirEnd.Z();
  }
  TVector3 startDir = {-xTotStart/endPoint, -yTotStart/endPoint, -zTotStart/endPoint};
  TVector3 endDir = {xTotEnd/endPoint, yTotEnd/endPoint, zTotEnd/endPoint};

  return std::make_pair(startDir, endDir);

} // CRTT0MatchAlg::TrackDirectionAverage()

std::pair<crt::CRTHit, double> CRTT0MatchAlg::ClosestCRTHit(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc){

  // Calculate direction as an average over directions
  std::pair<TVector3, TVector3> startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
  TVector3 startDir = startEndDir.first;
  TVector3 endDir = startEndDir.second;

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();

  // ====================== Matching Algorithm ========================== //
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = TrackT0Range(start.X(), end.X(), tpc);
  std::vector<std::pair<crt::CRTHit, double>> t0Candidates;

  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3;
    if (!(crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.)) continue;
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    // Calculate the distance between the crossing point and the CRT hit
    double startDist = DistOfClosestApproach(start, startDir, crtHit, tpc, crtTime);
    double endDist = DistOfClosestApproach(end, endDir, crtHit, tpc, crtTime);
    // If the distance is less than some limit record the time
    if ((crtPoint-start).Mag() < (crtPoint-end).Mag()){ 
      t0Candidates.push_back(std::make_pair(crtHit, startDist));
    }
    else{
      t0Candidates.push_back(std::make_pair(crtHit, endDist));
    }
  
  }

  // Sort the candidates by distance
  std::sort(t0Candidates.begin(), t0Candidates.end(), [](auto& left, auto& right){
            return left.second < right.second;});

  if(t0Candidates.size() > 0){
    return t0Candidates[0];
  }
  crt::CRTHit hit;
  return std::make_pair(hit, -99999);

}

double CRTT0MatchAlg::T0FromCRTHits(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc){

  if (tpcTrack.Length() < fMinTrackLength) return -99999; 

  std::pair<crt::CRTHit, double> closestHit = ClosestCRTHit(tpcTrack, crtHits, tpc);
  if(closestHit.second == -99999) return -99999;

  double crtTime = ((double)(int)closestHit.first.ts1_ns) * 1e-3;
  if(closestHit.second < fDistanceLimit) return crtTime;

  return -99999;

}

}
