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

  return;

}
 

// Function to transform a CRTTrack into an expected reconstructed track
std::vector<RecoCRTTrack> CRTTrackMatchAlg::CrtToRecoTrack(crt::CRTTrack track, int id){

  // TODO: could be more detector agnositic (pos of cpa, drift directions, max length of track)
  std::vector<RecoCRTTrack> recoCrtTracks;
  // Get the time of the track
  double crtTime = ((double)(int)track.ts1_ns) * 1e-3; // [us]
  // Convert time into a x shift
  double xShift = crtTime * fDetectorProperties->DriftVelocity();

  // Shift track, remembering to take into account the tpc, if the track crosses the cpa and 
  //the size of the readout window
  TVector3 crtStart(track.x1_pos, track.y1_pos, track.z1_pos);
  TVector3 crtEnd(track.x2_pos, track.y2_pos, track.z2_pos);
  if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
  if(!track.complete){
    // Find point where track crosses bottom plane
    TVector3 diff = (crtEnd - crtStart).Unit();
    TVector3 newEnd = crtStart + 5000*diff; 
    crtEnd = newEnd;
  }

  TVector3 cpaCrossStart, cpaCrossEnd;
  // Calculate the expected reconstructed length, angle and start and end position
  std::vector<int> crtTpc;
  if(crtStart.X() < 0. && crtEnd.X() < 0.){

    // Track in TPC 0
    std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 0, 
                                                              id, crtTime, track.complete);
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
  }
  else if(crtStart.X() > 0. && crtEnd.X() > 0.){

    // Track in TPC 1
    std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 1, 
                                                              id, crtTime, track.complete);
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
  }
  else {
    // Track in both TPCs and will be split
    TVector3 direction = crtStart - crtEnd;
    double step = (0. - crtStart.X())/direction.X(); 
    TVector3 cpaCross(0., crtStart.Y() + step*direction.Y(), crtStart.Z() + step*direction.Z());
    cpaCrossStart = cpaCross;
    cpaCrossEnd = cpaCross;

    if(crtStart.X() < 0.){ 
      std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 0, 
                                                                 id, crtTime, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

      std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 1, 
                                                                 id, crtTime, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
    }
    else {
      std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 0, 
                                                                 id, crtTime, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

      std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 1, 
                                                                 id, crtTime, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
    }

  }

  return recoCrtTracks;

} // CRTTrackMatchAlg::CRTToRecoTrack()


// Function to shift CRTTrack in X and work out how much is reconstructed
std::vector<RecoCRTTrack> CRTTrackMatchAlg::CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                                               int tpc, int id, double time, bool complete){

  std::vector<RecoCRTTrack> recoCrtTracks;

  // Get the true entry and exit points in the TPC
  double xmin = -2.0 * fGeometryService->DetHalfWidth();
  double xmax = 2.0 * fGeometryService->DetHalfWidth();
  double ymin = -fGeometryService->DetHalfHeight();
  double ymax = fGeometryService->DetHalfHeight();
  double zmin = 0.;
  double zmax = fGeometryService->DetLength();

  // Get track info
  TVector3 diff = end - start;
  TVector3 startTPC, endTPC;
  bool first = true;

  // Loop over trajectory points
  int npts = 1000;
  for (int traj_i = 0; traj_i <= npts; traj_i++){
    TVector3 trajPoint = start + ((traj_i)/((double)npts))*diff;
    // Check if point is within reconstructable volume
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
      first = false;
      startTPC = trajPoint;
    }
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
      endTPC = trajPoint;
    }
  }

  // Don't shift if not inside TPC
  if(startTPC == endTPC){
    return recoCrtTracks;
  }

  // Shift in x depending on TPC
  if(tpc == 1){
    // Track in TPC 1
    startTPC[0] -= shift;
    endTPC[0] -= shift;
  }
  else if(tpc == 0){
    // Track in TPC 0
    startTPC[0] += shift;
    endTPC[0] += shift;
  }
  
  double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
  double driftTimeMuS = (2.*fGeometryService->DetHalfWidth())/fDetectorProperties->DriftVelocity(); // [us]
  double deltaX = (readoutWindowMuS - driftTimeMuS) * fDetectorProperties->DriftVelocity(); // [cm]

  if(tpc == 0) xmax = deltaX;
  if(tpc == 1) xmin = -deltaX;

  // Get track info
  TVector3 diffTPC = endTPC - startTPC;
  TVector3 startCut, endCut;
  first = true;
  // Loop over trajectory points
  for (int traj_i = 0; traj_i <= npts; traj_i++){
    TVector3 trajPoint = startTPC + ((traj_i)/((double)npts))*diffTPC;
    // Check if point is within reconstructable volume
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
      first = false;
      startCut = trajPoint;
    }
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
      endCut = trajPoint;
    }
  }

  if(!(startCut.X() == endCut.X())){
    RecoCRTTrack recoCrtTrack = {id, tpc, startCut, endCut, time, complete};
    recoCrtTracks.push_back(recoCrtTrack);
  }

  return recoCrtTracks;

} // CRTTrackMatchAlg::CreateRecoCRTTrack()


// Function to calculate if a CRTTrack crosses the TPC volume
bool CRTTrackMatchAlg::CrossesTPC(crt::CRTTrack track){

  // Check if particle enters the TPC
  bool enters = false;
  double xmin = -2.0 * fGeometryService->DetHalfWidth();
  double xmax = 2.0 * fGeometryService->DetHalfWidth();
  double ymin = -fGeometryService->DetHalfHeight();
  double ymax = fGeometryService->DetHalfHeight();
  double zmin = 0.;
  double zmax = fGeometryService->DetLength();

  if(track.complete){
    // Get track info
    TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
    TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
    TVector3 diff = end - start;
    // Loop over trajectory points
    int npts = 100;
    for (int traj_i = 0; traj_i < npts; traj_i++){
      TVector3 trajPoint = start + ((traj_i+1)/100.)*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
        enters = true;
      }
    }
  }

  // If track just between top two planes
  else{
    //
    TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
    TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
    if(start.Y() < end.Y()) std::swap(start, end);
    TVector3 diff = (end - start).Unit();
    int npts = 100;
    for (int traj_i = 0; traj_i < npts; traj_i++){
      // TODO: length of track needs to be more detector agnostic
      TVector3 trajPoint = start + (100.*(traj_i+1))*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
        enters = true;
      }
    }
  }

  return enters;

} // CRTTrackMatchAlg::CrossesTPC()


}
