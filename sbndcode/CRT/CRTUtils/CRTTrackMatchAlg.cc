#include "CRTTrackMatchAlg.h"

namespace sbnd{

CRTTrackMatchAlg::CRTTrackMatchAlg(const Config& config){

  this->reconfigure(config);
  
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

}


CRTTrackMatchAlg::CRTTrackMatchAlg(){

  fGeometryService = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

}


CRTTrackMatchAlg::~CRTTrackMatchAlg(){

}


void CRTTrackMatchAlg::reconfigure(const Config& config){

  fMaxAngleDiff = config.MaxAngleDiff();
  fMaxDistance = config.MaxDistance();
  fTPCTrackLabel = config.TPCTrackLabel();

  return;

}
 

// Calculate intersection between CRT track and TPC (AABB Ray-Box intersection)
// (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
std::pair<TVector3, TVector3> CRTTrackMatchAlg::TpcIntersection(const geo::TPCGeo& tpcGeo, crt::CRTTrack track){

  // Find the intersection between the track and the TPC
  TVector3 start (track.x1_pos, track.y1_pos, track.z1_pos);
  TVector3 end (track.x2_pos, track.y2_pos, track.z2_pos);
  TVector3 min (tpcGeo.MinX(), tpcGeo.MinY(), tpcGeo.MinZ());
  TVector3 max (tpcGeo.MaxX(), tpcGeo.MaxY(), tpcGeo.MaxZ());

  std::pair<TVector3, TVector3> intersection = CRTCommonUtils::CubeIntersection(min, max, start, end);
  if(intersection.first.X() != -99999) return intersection;

  // Allow variations in track start/end points
  // Try the corners of the allowed region
  TVector3 min1 (track.x1_pos - track.x1_err, track.y1_pos - track.y1_err, track.z1_pos - track.z1_err);
  TVector3 min2 (track.x2_pos - track.x2_err, track.y2_pos - track.y2_err, track.z2_pos - track.z2_err);
  intersection = CRTCommonUtils::CubeIntersection(min, max, min1, min2);
  if(intersection.first.X() != -99999) return intersection;

  TVector3 max1 (track.x1_pos + track.x1_err, track.y1_pos + track.y1_err, track.z1_pos + track.z1_err);
  TVector3 max2 (track.x2_pos + track.x2_err, track.y2_pos + track.y2_err, track.z2_pos + track.z2_err);
  intersection = CRTCommonUtils::CubeIntersection(min, max, max1, max2);
  return intersection;

}


// Function to calculate if a CRTTrack crosses the TPC volume
bool CRTTrackMatchAlg::CrossesTPC(crt::CRTTrack track){

  for(size_t c = 0; c < fGeometryService->Ncryostats(); c++){
    const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(c);
    for(size_t t = 0; t < cryostat.NTPC(); t++){
      const geo::TPCGeo& tpcGeo = cryostat.TPC(t);
      std::pair<TVector3, TVector3> intersection = TpcIntersection(tpcGeo, track);
      if(intersection.first.X() != -99999) return true;
    }
  }
  return false;

} // CRTTrackMatchAlg::CrossesTPC()


// Function to calculate if a CRTTrack crosses the wire planes
bool CRTTrackMatchAlg::CrossesAPA(crt::CRTTrack track){

  for(size_t c = 0; c < fGeometryService->Ncryostats(); c++){
    const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(c);
    for(size_t t = 0; t < cryostat.NTPC(); t++){
      const geo::TPCGeo& tpcGeo = cryostat.TPC(t);
      std::pair<TVector3, TVector3> intersection = TpcIntersection(tpcGeo, track);
      if(std::abs(intersection.first.X()) == fTpcGeo.CpaWidth()) return true;
      if(std::abs(intersection.second.X()) == fTpcGeo.CpaWidth()) return true;
    }
  }
  return false;

} // CRTTrackMatchAlg::CrossesAPA()


double CRTTrackMatchAlg::T0FromCRTTracks(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event){

  std::pair<crt::CRTTrack, double> closestDCA = ClosestCRTTrackByDCA(tpcTrack, crtTracks, event);

  if(closestDCA.second == -99999 || closestDCA > fMaxDistance) return -99999;

  double crtTime = ((double)(int)closestDCA.first.ts1_ns) * 1e-3; // [us]

  return crtTime;

}


// Find the closest valid matching CRT track ID
int CRTTrackMatchAlg::GetMatchedCRTTrackId(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event){

  std::pair<crt::CRTTrack, double> closestDCA = ClosestCRTTrackByDCA(tpcTrack, crtTracks, event);

  if(closestDCA.second == -99999 || closestDCA.second > fMaxAngleDiff) return -99999;

  int crt_i = 0;
  for(auto const& track : crtTracks){
    if(fCrtBackTrack.TrackCompare(closestDCA.first, track)) return crt_i;
    crt_i++;
  }

  return -99999;

}


// Get all CRT tracks that cross the right TPC within an allowed time
std::vector<crt::CRTTrack> CRTTrackMatchAlg::AllPossibleCRTTracks(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event){

   std::vector<crt::CRTTrack> trackCandidates;

  // Get the hits associated with the tpc track
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

  // Get the drift direction (0 for stitched tracks)
  int driftDirection = fTpcGeo.DriftDirectionFromHits(hits);

  // Get the TPC Geo object from the tpc track
  geo::TPCID tpcID = hits[0]->WireID().asTPCID();
  const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);

  // Loop over the crt tracks
  for(auto const& crtTrack : crtTracks){
    // Calculate the intersection points for that TPC
    std::pair<TVector3, TVector3> intersection = TpcIntersection(tpcGeo, crtTrack);

    // Skip if it doesn't intersect
    if(intersection.first.X() == -99999) continue;

    // Shift the track to the CRT track
    double crtTime = ((double)(int)crtTrack.ts1_ns) * 1e-3; // [us]
    double shift = driftDirection * crtTime * fDetectorProperties->DriftVelocity();
    geo::Point_t start = tpcTrack.Vertex();
    geo::Point_t end = tpcTrack.End();
    start.SetX(start.X() + shift);
    end.SetX(end.X() + shift);

    // Check the track is fully contained in the TPC
    if(!fTpcGeo.InsideTPC(start, tpcGeo, 2.) && shift != 0) continue;
    if(!fTpcGeo.InsideTPC(end, tpcGeo, 2.) && shift != 0) continue;

    trackCandidates.push_back(crtTrack);
    
  }

  return trackCandidates;
}


// Find the closest matching crt track by angle between tracks within angle and DCA limits
std::pair<crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByAngle(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event, double minDCA){

  // Get the hits associated with the tpc track
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

  // Get the drift direction (0 for stitched tracks)
  int driftDirection = fTpcGeo.DriftDirectionFromHits(hits);

  std::vector<crt::CRTTrack> possTracks = AllPossibleCRTTracks(tpcTrack, crtTracks, event);

  std::vector<std::pair<crt::CRTTrack, double>> candidates;
  for(auto const& possTrack : possTracks){
    double angle = AngleBetweenTracks(tpcTrack, possTrack);

    if(minDCA != -1){
      if(minDCA == 0) minDCA = fMaxDistance;
      double crtTime = ((double)(int)possTrack.ts1_ns) * 1e-3; // [us]
      double shift = driftDirection * crtTime * fDetectorProperties->DriftVelocity();
      double DCA = AveDCABetweenTracks(tpcTrack, possTrack, shift);
      if(DCA > minDCA) continue;
    }
      
    candidates.push_back(std::make_pair(possTrack, angle));
  }

  std::sort(candidates.begin(), candidates.end(), [](auto& left, auto& right){
            return left.second < right.second;});

  if(candidates.size() > 0){
    return candidates[0];
  }
  crt::CRTTrack track;
  return std::make_pair(track, -99999);
}


// Find the closest matching crt track by average DCA between tracks within angle and DCA limits
std::pair<crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByDCA(recob::Track tpcTrack, std::vector<crt::CRTTrack> crtTracks, const art::Event& event, double minAngle){

  // Get the hits associated with the tpc track
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

  // Get the drift direction (0 for stitched tracks)
  int driftDirection = fTpcGeo.DriftDirectionFromHits(hits);

  std::vector<crt::CRTTrack> possTracks = AllPossibleCRTTracks(tpcTrack, crtTracks, event);

  std::vector<std::pair<crt::CRTTrack, double>> candidates;
  for(auto const& possTrack : possTracks){

    double crtTime = ((double)(int)possTrack.ts1_ns) * 1e-3; // [us]
    double shift = driftDirection * crtTime * fDetectorProperties->DriftVelocity();

    double DCA = AveDCABetweenTracks(tpcTrack, possTrack, shift);

    if(minAngle != -1){
      if(minAngle == 0) minAngle = fMaxAngleDiff;
      double angle = AngleBetweenTracks(tpcTrack, possTrack);
      if(angle > minAngle) continue;
    }
    candidates.push_back(std::make_pair(possTrack, DCA));
  }

  std::sort(candidates.begin(), candidates.end(), [](auto& left, auto& right){
            return left.second < right.second;});

  if(candidates.size() > 0){
    return candidates[0];
  }
  crt::CRTTrack track;
  return std::make_pair(track, -99999);

}


// Calculate the angle between tracks assuming start is at the largest Y
double CRTTrackMatchAlg::AngleBetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack){

  // Calculate the angle between the tracks
  TVector3 crtStart (crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
  TVector3 crtEnd (crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
  if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);

  TVector3 tpcStart = tpcTrack.Vertex<TVector3>();
  TVector3 tpcEnd  = tpcTrack.End<TVector3>();
  if(tpcStart.Y() < tpcEnd.Y()) std::swap(tpcStart, tpcEnd);

  return (tpcStart - tpcEnd).Angle(crtStart - crtEnd);

}


// Calculate the average DCA between tracks
double CRTTrackMatchAlg::AveDCABetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack, double shift){

  TVector3 crtStart (crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
  TVector3 crtEnd (crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
  if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
  double denominator = (crtEnd - crtStart).Mag();

  size_t npts = tpcTrack.NumberTrajectoryPoints();

  double aveDCA = 0;
  int usedPts = 0;
  for(size_t i = 0; i < npts; i++){
    TVector3 point = tpcTrack.LocationAtPoint<TVector3>(i);

    // Pandora produces dummy points
    if(!tpcTrack.HasValidPoint(i)) continue;

    point.SetX(point.X() + shift);
    aveDCA += (point - crtStart).Cross(point - crtEnd).Mag()/denominator;
    usedPts++;
  }

  return aveDCA/usedPts;

}


// Calculate the average DCA between tracks
double CRTTrackMatchAlg::AveDCABetweenTracks(recob::Track tpcTrack, crt::CRTTrack crtTrack, const art::Event& event){

  // Get the hits associated with the tpc track
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

  // Get the drift direction (0 for stitched tracks)
  int driftDirection = fTpcGeo.DriftDirectionFromHits(hits);
  double crtTime = ((double)(int)crtTrack.ts1_ns) * 1e-3; // [us]
  double shift = driftDirection * crtTime * fDetectorProperties->DriftVelocity();

  TVector3 crtStart (crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
  TVector3 crtEnd (crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
  if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
  double denominator = (crtEnd - crtStart).Mag();

  size_t npts = tpcTrack.NumberTrajectoryPoints();

  double aveDCA = 0;
  int usedPts = 0;
  for(size_t i = 0; i < npts; i++){
    TVector3 point = tpcTrack.LocationAtPoint<TVector3>(i);

    // Pandora produces dummy points
    if(!tpcTrack.HasValidPoint(i)) continue;

    point.SetX(point.X() + shift);
    aveDCA += (point - crtStart).Cross(point - crtEnd).Mag()/denominator;
    usedPts++;
  }

  return aveDCA/usedPts;

}


}
