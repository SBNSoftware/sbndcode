#include "CRTT0MatchAlg.h"

namespace sbnd{

CRTT0MatchAlg::CRTT0MatchAlg(const Config& config) : CRTT0MatchAlg(config, lar::providerFrom<geo::Geometry>(), 
lar::providerFrom<spacecharge::SpaceChargeService>()) {}

CRTT0MatchAlg::CRTT0MatchAlg(const Config& config, geo::GeometryCore const *GeometryService,  spacecharge::SpaceCharge const *SCE){
  
  
  this->reconfigure(config);
  fGeometryService = GeometryService;
  fSCE = SCE;
}
CRTT0MatchAlg::CRTT0MatchAlg() = default;
  

 
  
void CRTT0MatchAlg::reconfigure(const Config& config){
    
  fMinTrackLength = config.MinTrackLength();
  fTrackDirectionFrac = config.TrackDirectionFrac();
  fDistanceLimit = config.DistanceLimit();
  fTSMode = config.TSMode();
  fTimeCorrection = config.TimeCorrection();
  fTPCTrackLabel = config.TPCTrackLabel();
  fSCEposCorr = config.SCEposCorr();
  fDirMethod = config.DirMethod();
  fDCAuseBox = config.DCAuseBox();
  fDCAoverLength = config.DCAoverLength();
  fDoverLLimit = config.DoverLLimit();
  fPEcut = config.PEcut();
  fMaxUncert = config.MaxUncert();
  //  fDistEndpointAVedge = config.DistEndpointAVedge();

  return;

}

matchCand makeNULLmc (){
    sbn::crt::CRTHit hit;
    matchCand null;
    null.thishit = hit;
    null.t0 = -99999;
    null.dca = -99999;
    null.extrapLen = -99999;
    return null;
}

 

// Utility function that determines the possible t0 range of a track
std::pair<double, double> CRTT0MatchAlg::TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
                                                      double startX, double endX, int driftDirection, std::pair<double, double> xLimits){

  // If track is stitched return zeros
  if(driftDirection == 0) return std::make_pair(0, 0);

  //std::pair<double, double> result; // unused
  double Vd = driftDirection * detProp.DriftVelocity();

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

  //  if (t0min>2500)  std::cout << " t0 min " << t0min << " t0max " << t0max << std::endl;
  return std::make_pair(std::min(t0min, t0max), std::max(t0min, t0max));


} // CRTT0MatchAlg::TrackT0Range()


double CRTT0MatchAlg::DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                            TVector3 trackPos, TVector3 trackDir, sbn::crt::CRTHit crtHit, int driftDirection, double t0){

  //double minDist = 99999;

  // Convert the t0 into an x shift
  double xshift = driftDirection* t0 * detProp.DriftVelocity();
  trackPos[0] += xshift;

  if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
    geo::Point_t temppt = {trackPos.X(),trackPos.Y(),trackPos.Z()};
    geo::TPCID tpcid = fGeometryService->PositionToTPCID(temppt);
    geo::Vector_t  fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
    trackPos[0] += fPosOffsets.X();
    trackPos[1] += fPosOffsets.Y();
    trackPos[2] += fPosOffsets.Z();
  }

  TVector3 end = trackPos + trackDir;

  // calculate distance of closest approach (DCA)
  //  default is the distance to the point specified by the CRT hit (Simple DCA)
  //    useBox is the distance to the closest edge of the rectangle with the CRT hit at the center and the sides defined
  //   the position uncertainties on the CRT hits.
  double thisdca;

  if (fDCAuseBox) thisdca =   CRTCommonUtils::DistToCrtHit(crtHit, trackPos, end);
  else thisdca =  CRTCommonUtils::SimpleDCA(crtHit, trackPos, trackDir);
  return thisdca;

} // CRTT0MatchAlg::DistToOfClosestApproach()


std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirectionAverage(recob::Track track, double frac)
  {
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


  std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirection(detinfo::DetectorPropertiesData const& detProp,recob::Track track, double frac, double CRTtime, int driftDirection){
							      
    size_t nTrackPoints = track.NPoints();
    int midPt = (int)floor(nTrackPoints*frac);
    geo::Point_t startP = track.Start();
    geo::Point_t endP = track.End();
    geo::Point_t midP = track.LocationAtPoint(midPt);

    double xshift = driftDirection * CRTtime * detProp.DriftVelocity();
    TVector3  startPoint = {startP.X()+xshift,startP.Y(),startP.Z()};
    TVector3  endPoint = {endP.X()+xshift,endP.Y(),endP.Z()};
    TVector3  midPoint = {midP.X()+xshift,midP.Y(),midP.Z()};
    if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {

      // Apply the shift depending on which TPC the track is in                                 
      geo::Point_t fTrackPos = startP;
      fTrackPos.SetX(startPoint.X());
      geo::TPCID tpcid = fGeometryService->PositionToTPCID(fTrackPos);                        
      geo::Vector_t fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      startPoint.SetX(fTrackPos.X() + fPosOffsets.X());                                       
      startPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());                                       
      startPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());                                       
      fTrackPos = endP;
      fTrackPos.SetX(endPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //      fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      endPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      endPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      endPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
      fTrackPos = midP;
      fTrackPos.SetX(midPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      midPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      midPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      midPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
    }
    
    TVector3 startDir = {midPoint.X()-startPoint.X(),midPoint.Y()-startPoint.Y(),midPoint.Z()-startPoint.Z()};
    float norm = startDir.Mag();
    if (norm>0)  startDir *=(1.0/norm);
    TVector3 endDir = {midPoint.X()-endPoint.X(),midPoint.Y()-endPoint.Y(),midPoint.Z()-endPoint.Z()};    
    norm = endDir.Mag();
    if (norm>0)  endDir *=(1.0/norm);
    
    return std::make_pair(startDir, endDir);
    
  } // CRTT0MatchAlg::TrackDirection()                                                                  

std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirectionAverageFromPoints(recob::Track track, double frac){

  // Calculate direction as an average over directions
  size_t nTrackPoints = track.NumberTrajectoryPoints();
  recob::TrackTrajectory trajectory  = track.Trajectory();
  std::vector<TVector3> validPoints;
  for(size_t i = 0; i < nTrackPoints; i++){
    if(trajectory.FlagsAtPoint(i) != recob::TrajectoryPointFlags::InvalidHitIndex) continue;
    validPoints.push_back(track.LocationAtPoint<TVector3>(i));
  }

  size_t nValidPoints = validPoints.size();
  int endPoint = (int)floor(nValidPoints*frac);
  TVector3 startDir = validPoints.at(0) - validPoints.at(endPoint-1);
  TVector3 endDir = validPoints.at(nValidPoints - 1) - validPoints.at(nValidPoints - (endPoint));

  return std::make_pair(startDir.Unit(), endDir.Unit());

} // CRTT0MatchAlg::TrackDirectionAverageFromPoints()


// Keeping ClosestCRTHit function for backward compatibility only
// *** use GetClosestCRTHit instead

std::pair<sbn::crt::CRTHit, double> CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								 recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return ClosestCRTHit(detProp, tpcTrack, hits, crtHits);
}


std::pair<sbn::crt::CRTHit, double>  CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
				       recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbn::crt::CRTHit> crtHits) {

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();
  // Get the drift direction from the TPC
  int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
  std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

  return ClosestCRTHit(detProp, tpcTrack, t0MinMax, crtHits, driftDirection);
}

std::pair<sbn::crt::CRTHit, double> CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
				       recob::Track tpcTrack, std::pair<double, double> t0MinMax, std::vector<sbn::crt::CRTHit> crtHits, int driftDirection) {
  matchCand bestmatch = GetClosestCRTHit(detProp, tpcTrack,t0MinMax,crtHits,driftDirection);
  return std::make_pair(bestmatch.thishit,bestmatch.dca);

}


matchCand CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
				       recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbn::crt::CRTHit> crtHits) {

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();
  // Get the drift direction from the TPC
  int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
  std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

  return GetClosestCRTHit(detProp, tpcTrack, t0MinMax, crtHits, driftDirection);
}

matchCand CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								 recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);
}


matchCand CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
				       recob::Track tpcTrack, std::pair<double, double> t0MinMax, std::vector<sbn::crt::CRTHit> crtHits, int driftDirection) {
  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();


  // ====================== Matching Algorithm ========================== //
  //  std::vector<std::pair<sbn::crt::CRTHit, double>> t0Candidates;
  std::vector<matchCand> t0Candidates;



  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = -99999.;  // units are us
    if (fTSMode == 1) {
      crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3 + fTimeCorrection;
    }
    else {
      crtTime = ((double)(int)crtHit.ts0_ns) * 1e-3 + fTimeCorrection;
    }
    //    if (crtTime>3000) std::cout << "crt hit times " << crtTime << std::endl;
    // If track is stitched then try all hits
    if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) continue;
    // cut on CRT hit PE value
    if (crtHit.peshit<fPEcut) continue;
    if (crtHit.x_err>fMaxUncert) continue;
    if (crtHit.y_err>fMaxUncert) continue;
    if (crtHit.z_err>fMaxUncert) continue;

    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    //Calculate Track direction
    std::pair<TVector3, TVector3> startEndDir;
    // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
    if (fDirMethod==2)  startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
    else startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
    TVector3 startDir = startEndDir.first;
    TVector3 endDir = startEndDir.second;
    
    // Calculate the distance between the crossing point and the CRT hit, SCE corrections are done inside but dropped
    double startDist = DistOfClosestApproach(detProp, start, startDir, crtHit, driftDirection, crtTime);
    double endDist = DistOfClosestApproach(detProp, end, endDir, crtHit, driftDirection, crtTime);

    
    double xshift = driftDirection * crtTime * detProp.DriftVelocity();
    auto thisstart = start; 
    thisstart.SetX(start.X()+xshift);
    auto thisend = end; 
    thisend.SetX(end.X()+xshift);

    // repeat SCE correction for endpoints
    if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
      geo::Point_t temppt = {thisstart.X(),thisstart.Y(),thisstart.Z()};
      geo::TPCID tpcid = fGeometryService->PositionToTPCID(temppt);
      geo::Vector_t  fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
      thisstart[0] += fPosOffsets.X();
      thisstart[1] += fPosOffsets.Y();
      thisstart[2] += fPosOffsets.Z();
      temppt.SetX(thisend.X());
      temppt.SetY(thisend.Y());
      temppt.SetZ(thisend.Z());
      tpcid = fGeometryService->PositionToTPCID(temppt);
      fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
      thisend[0] += fPosOffsets.X();
      thisend[1] += fPosOffsets.Y();
      thisend[2] += fPosOffsets.Z();
    }


    matchCand newmc = makeNULLmc();
    if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
    double distS = (crtPoint-thisstart).Mag();
    double distE =  (crtPoint-thisend).Mag();
    // std::cout << " distS " << distS << " distE " << distE << std::endl;
    // std::cout << "startdis " << startDist << " endDist " << endDist << " dca "  << std::endl;
    // std::cout << " doL start " << startDist/distS << " doL end " << endDist/distE << std::endl;

    if (distS < distE){ 
      newmc.thishit = crtHit;
      newmc.t0= crtTime;
      newmc.dca = startDist;
      newmc.extrapLen = distS;
      t0Candidates.push_back(newmc);
    }
    else{
      newmc.thishit = crtHit;
      newmc.t0= crtTime;
      newmc.dca = endDist;
      newmc.extrapLen = distE;
      t0Candidates.push_back(newmc);
    }
    }
  }


  //  std::cout << " found " << t0Candidates.size() << " candidates" << std::endl;
  matchCand bestmatch = makeNULLmc();
  if(t0Candidates.size() > 0){
  // Find candidate with shortest DCA or DCA/L value
    bestmatch=t0Candidates[0];
    double sin_angle = bestmatch.dca/bestmatch.extrapLen;
    if (fDCAoverLength) { // Use dca/extrapLen to judge best
      for(auto &thisCand : t0Candidates){
	double this_sin_angle = thisCand.dca/thisCand.extrapLen;
	if (bestmatch.dca<0 )	bestmatch=thisCand;
	else if (this_sin_angle<sin_angle && thisCand.dca>=0)	bestmatch=thisCand;
      }
    }
    else { // use Dca to judge best
      for(auto &thisCand : t0Candidates){
	if (bestmatch.dca<0 )	bestmatch=thisCand;
	else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0)	bestmatch=thisCand;
      }
    }
  }

  //  std::cout << "best match has dca of " << bestmatch.dca << std::endl;
  return bestmatch;

}


double CRTT0MatchAlg::T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                                    recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, const art::Event& event){
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
  return T0FromCRTHits(detProp, tpcTrack, hits, crtHits);
}

double CRTT0MatchAlg::T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
                                    recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbn::crt::CRTHit> crtHits) {

  if (tpcTrack.Length() < fMinTrackLength) return -99999; 

  matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);
  if(closestHit.dca <0) return -99999;

  double crtTime;
  if (fTSMode == 1) {
    crtTime = ((double)(int)closestHit.thishit.ts1_ns) * 1e-3 + fTimeCorrection;
  }
  else {
    crtTime = ((double)(int)closestHit.thishit.ts0_ns) * 1e-3 + fTimeCorrection;
  }
  if (closestHit.dca < fDistanceLimit && (closestHit.dca/closestHit.extrapLen) < fDoverLLimit) return crtTime;

  return -99999;

}

  std::pair<double, double>  CRTT0MatchAlg::T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
					     recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, const art::Event& event){
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
  return T0AndDCAFromCRTHits(detProp, tpcTrack, hits, crtHits);
}

std::pair<double, double> CRTT0MatchAlg::T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
					     recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbn::crt::CRTHit> crtHits) {

  if (tpcTrack.Length() < fMinTrackLength) return std::make_pair(-9999., -9999.);

  matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);

  if(closestHit.dca < 0 ) return std::make_pair(-9999., -9999.);
  if (closestHit.dca < fDistanceLimit && (closestHit.dca/closestHit.extrapLen) < fDoverLLimit) return std::make_pair(closestHit.t0, closestHit.dca);

  return std::make_pair(-9999., -9999.);


}

}
