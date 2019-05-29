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

  fTPCTrackLabel = config.TPCTrackLabel();

  return;

}
 

// Utility function that determines the possible t0 range of a track
std::pair<double, double> CRTT0MatchAlg::TrackT0Range(double startX, double endX, int driftDirection, std::pair<double, double> xLimits){

  // If track is stitched return zeros
  if(driftDirection == 0) return std::make_pair(0, 0);

  std::pair<double, double> result;
  double Vd = driftDirection * fDetectorProperties->DriftVelocity();

  // Shift the most postive end to the most positive limit
  double maxX = std::max(startX, endX);
  double maxLimit = std::max(xLimits.first, xLimits.second);
  double maxShift = maxLimit - maxX;
  // Shift the most negative end to the most negative limit
  double minX = std::min(startX, endX);
  double minLimit = std::min(xLimits.first, xLimits.second);
  double minShift = minLimit - minX;
  // Convert to time
  double t0max = maxShift/Vd;
  double t0min = minShift/Vd;

  return std::make_pair(std::min(t0min, t0max), std::max(t0min, t0max));

} // CRTT0MatchAlg::TrackT0Range()


double CRTT0MatchAlg::DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, crt::CRTHit crtHit, int driftDirection, double t0){

  //double minDist = 99999;

  // Convert the t0 into an x shift
  double shift = t0 * fDetectorProperties->DriftVelocity();
  // Apply the shift depending on which TPC the track is in
  trackPos[0] += driftDirection * shift;

  TVector3 end = trackPos + trackDir;

  return CRTCommonUtils::DistToCrtHit(crtHit, trackPos, end);

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


std::pair<crt::CRTHit, double> CRTT0MatchAlg::ClosestCRTHit(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event){

  // Get the drift direction from the TPC
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
  int driftDirection = fTpcGeo.DriftDirectionFromHits(hits);
  std::pair<double, double> xLimits = fTpcGeo.XLimitsFromHits(hits);

  // Calculate direction as an average over directions
  std::pair<TVector3, TVector3> startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
  TVector3 startDir = startEndDir.first;
  TVector3 endDir = startEndDir.second;

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();

  // ====================== Matching Algorithm ========================== //
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = TrackT0Range(start.X(), end.X(), driftDirection, xLimits);
  std::vector<std::pair<crt::CRTHit, double>> t0Candidates;

  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3;
    // If track is stitched then try all hits
    if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) continue;
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    // Calculate the distance between the crossing point and the CRT hit
    double startDist = DistOfClosestApproach(start, startDir, crtHit, driftDirection, crtTime);
    double endDist = DistOfClosestApproach(end, endDir, crtHit, driftDirection, crtTime);
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


double CRTT0MatchAlg::T0FromCRTHits(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event){

  if (tpcTrack.Length() < fMinTrackLength) return -99999; 

  std::pair<crt::CRTHit, double> closestHit = ClosestCRTHit(tpcTrack, crtHits, event);
  if(closestHit.second == -99999) return -99999;

  double crtTime = ((double)(int)closestHit.first.ts1_ns) * 1e-3;
  if(closestHit.second < fDistanceLimit) return crtTime;

  return -99999;

}

}
