#include "CRTSpacePointMatchAlg.h"

namespace sbnd{

  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg(const Config& config) : CRTSpacePointMatchAlg(config, lar::providerFrom<geo::Geometry>()) {}

  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg(const Config& config, geo::GeometryCore const *GeometryService) {
  
  
    this->reconfigure(config);
    fGeometryService = GeometryService;
  }
  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg() = default;
  

 
  
  void CRTSpacePointMatchAlg::reconfigure(const Config& config){
    
    fMinTrackLength = config.MinTrackLength();
    fTrackDirectionFrac = config.TrackDirectionFrac();
    fDistanceLimit = config.DistanceLimit();
    fTimeCorrection = config.TimeCorrection();
    fTPCTrackLabel = config.TPCTrackLabel();
    fDirMethod = config.DirMethod();
    fDCAuseBox = config.DCAuseBox();
    fDCAoverLength = config.DCAoverLength();
    fDoverLLimit = config.DoverLLimit();
    fPEcut = config.PEcut();
    fMaxUncert = config.MaxUncert();

    return;

  }

  matchCand makeNULLmc (){
    sbnd::crt::CRTSpacePoint sp;
    matchCand null;
    null.thisSP = sp;
    null.t0 = -99999;
    null.dca = -99999;
    null.extrapLen = -99999;
    return null;
  }

 

  // Utility function that determines the possible t0 range of a track
  std::pair<double, double> CRTSpacePointMatchAlg::TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
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


  } // CRTSpacePointMatchAlg::TrackT0Range()


  double CRTSpacePointMatchAlg::DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                                      geo::Point_t trackPos, geo::Vector_t trackDir, sbnd::crt::CRTSpacePoint crtSP, int driftDirection, double t0){

    //double minDist = 99999;

    // Convert the t0 into an x shift
    double xshift = driftDirection* t0 * detProp.DriftVelocity();
    trackPos.SetX(trackPos.X() + xshift);

    geo::Point_t end = trackPos + trackDir;

    // calculate distance of closest approach (DCA)
    //  default is the distance to the point specified by the CRTSpacePoint (Simple DCA)
    //    useBox is the distance to the closest edge of the rectangle with the CRTSpacePoint at the center and the sides defined
    //   the position uncertainties on the CRTSpacePoints.
    double thisdca;

    if (fDCAuseBox) thisdca = crt::CRTCommonUtils::DistToCRTSpacePoint(crtSP, trackPos, end);
    else thisdca = crt::CRTCommonUtils::SimpleDCA(crtSP, trackPos, trackDir);
    return thisdca;

  } // CRTSpacePointMatchAlg::DistToOfClosestApproach()


  std::pair<geo::Vector_t, geo::Vector_t> CRTSpacePointMatchAlg::TrackDirectionAverage(recob::Track track, double frac)
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
    geo::Vector_t startDir = {-xTotStart/endPoint, -yTotStart/endPoint, -zTotStart/endPoint};
    geo::Vector_t endDir = {xTotEnd/endPoint, yTotEnd/endPoint, zTotEnd/endPoint};

    return std::make_pair(startDir, endDir);

  } // CRTSpacePointMatchAlg::TrackDirectionAverage()


  std::pair<geo::Vector_t, geo::Vector_t> CRTSpacePointMatchAlg::TrackDirection(detinfo::DetectorPropertiesData const& detProp,recob::Track track, double frac, double CRTtime, int driftDirection){

    size_t nTrackPoints = track.NPoints();
    int midPt = (int)floor(nTrackPoints*frac);
    geo::Point_t startP = track.Start();
    geo::Point_t endP = track.End();
    geo::Point_t midP = track.LocationAtPoint(midPt);

    double xshift = driftDirection * CRTtime * detProp.DriftVelocity();
    geo::Point_t  startPoint = {startP.X()+xshift,startP.Y(),startP.Z()};
    geo::Point_t  endPoint = {endP.X()+xshift,endP.Y(),endP.Z()};
    geo::Point_t  midPoint = {midP.X()+xshift,midP.Y(),midP.Z()};
    geo::Vector_t startDir = {midPoint.X()-startPoint.X(),midPoint.Y()-startPoint.Y(),midPoint.Z()-startPoint.Z()};
    float norm = startDir.R();
    if (norm>0)  startDir *=(1.0/norm);
    geo::Vector_t endDir = {midPoint.X()-endPoint.X(),midPoint.Y()-endPoint.Y(),midPoint.Z()-endPoint.Z()};    
    norm = endDir.R();
    if (norm>0)  endDir *=(1.0/norm);
    
    return std::make_pair(startDir, endDir);
    
  } // CRTSpacePointMatchAlg::TrackDirection()                                                                  

  std::pair<geo::Vector_t, geo::Vector_t> CRTSpacePointMatchAlg::TrackDirectionAverageFromPoints(recob::Track track, double frac){

    // Calculate direction as an average over directions
    size_t nTrackPoints = track.NumberTrajectoryPoints();
    recob::TrackTrajectory trajectory  = track.Trajectory();
    std::vector<geo::Point_t> validPoints;
    for(size_t i = 0; i < nTrackPoints; i++){
      if(trajectory.FlagsAtPoint(i) != recob::TrajectoryPointFlags::InvalidHitIndex) continue;
      validPoints.push_back(track.LocationAtPoint(i));
    }

    size_t nValidPoints = validPoints.size();
    int endPoint = (int)floor(nValidPoints*frac);
    geo::Vector_t startDir = validPoints.at(0) - validPoints.at(endPoint-1);
    geo::Vector_t endDir = validPoints.at(nValidPoints - 1) - validPoints.at(nValidPoints - (endPoint));

    return std::make_pair(startDir.Unit(), endDir.Unit());

  } // CRTSpacePointMatchAlg::TrackDirectionAverageFromPoints()



  matchCand CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                                           recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSPs) {

    geo::Point_t start = tpcTrack.Vertex();
    geo::Point_t end = tpcTrack.End();
    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTSpacePoint(detProp, tpcTrack, t0MinMax, crtSPs, driftDirection);
  }

  matchCand CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                                           recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSPs, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return GetClosestCRTSpacePoint(detProp, tpcTrack, hits, crtSPs);
  }


  matchCand CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                                           recob::Track tpcTrack, std::pair<double, double> t0MinMax, std::vector<sbnd::crt::CRTSpacePoint> crtSPs, int driftDirection) {
    geo::Point_t start = tpcTrack.Vertex();
    geo::Point_t end = tpcTrack.End();


    // ====================== Matching Algorithm ========================== //
    //  std::vector<std::pair<sbnd::crt::CRTSpacePoint, double>> t0Candidates;
    std::vector<matchCand> t0Candidates;



    // Loop over all the CRTSpacePoints
    for(auto &crtSP : crtSPs){
      // Check if space point is within the allowed t0 range
      double crtTime = crtSP.Time() * 1e-3 + fTimeCorrection;

      // If track is stitched then try all hits
      if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.)
            || t0MinMax.first == t0MinMax.second)) continue;
      // cut on CRTSpacePoint PE value
      if (crtSP.PE()<fPEcut) continue;
      if (crtSP.XErr()>fMaxUncert) continue;
      if (crtSP.YErr()>fMaxUncert) continue;
      if (crtSP.ZErr()>fMaxUncert) continue;

      geo::Point_t crtPoint = crtSP.Pos();
  
      //Calculate Track direction
      std::pair<geo::Vector_t, geo::Vector_t> startEndDir;
      // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
      if (fDirMethod==2)  startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
      else startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
      geo::Vector_t startDir = startEndDir.first;
      geo::Vector_t endDir = startEndDir.second;
    
      // Calculate the distance between the crossing point and the CRTSpacePoint
      double startDist = DistOfClosestApproach(detProp, start, startDir, crtSP, driftDirection, crtTime);
      double endDist = DistOfClosestApproach(detProp, end, endDir, crtSP, driftDirection, crtTime);

    
      double xshift = driftDirection * crtTime * detProp.DriftVelocity();
      auto thisstart = start;
      thisstart.SetX(start.X()+xshift);
      auto thisend = end;
      thisend.SetX(end.X()+xshift);

      matchCand newmc = makeNULLmc();
      if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
        double distS = (crtPoint-thisstart).R();
        double distE =  (crtPoint-thisend).R();
        // std::cout << " distS " << distS << " distE " << distE << std::endl;
        // std::cout << "startdis " << startDist << " endDist " << endDist << " dca "  << std::endl;
        // std::cout << " doL start " << startDist/distS << " doL end " << endDist/distE << std::endl;

        if (distS < distE){
          newmc.thisSP = crtSP;
          newmc.t0= crtTime;
          newmc.dca = startDist;
          newmc.extrapLen = distS;
          t0Candidates.push_back(newmc);
        }
        else{
          newmc.thisSP = crtSP;
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
          if (bestmatch.dca<0 )   bestmatch=thisCand;
          else if (this_sin_angle<sin_angle && thisCand.dca>=0)   bestmatch=thisCand;
        }
      }
      else { // use Dca to judge best
        for(auto &thisCand : t0Candidates){
          if (bestmatch.dca<0 )   bestmatch=thisCand;
          else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0) bestmatch=thisCand;
        }
      }
    }

    //  std::cout << "best match has dca of " << bestmatch.dca << std::endl;
    return bestmatch;

  }


  double CRTSpacePointMatchAlg::T0FromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                     recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSPs, const art::Event& event){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return T0FromCRTSpacePoints(detProp, tpcTrack, hits, crtSPs);
  }

  double CRTSpacePointMatchAlg::T0FromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                     recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSPs) {

    if (tpcTrack.Length() < fMinTrackLength) return -99999;

    matchCand closestSP = GetClosestCRTSpacePoint(detProp, tpcTrack, hits, crtSPs);
    if(closestSP.dca <0) return -99999;

    double crtTime = closestSP.thisSP.Time() * 1e-3 + fTimeCorrection;

    if (closestSP.dca < fDistanceLimit && (closestSP.dca/closestSP.extrapLen) < fDoverLLimit) return crtTime;

    return -99999;

  }

  std::pair<double, double>  CRTSpacePointMatchAlg::T0AndDCAFromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                                               recob::Track tpcTrack, std::vector<sbnd::crt::CRTSpacePoint> crtSPs, const art::Event& event){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return T0AndDCAFromCRTSpacePoints(detProp, tpcTrack, hits, crtSPs);
  }

  std::pair<double, double> CRTSpacePointMatchAlg::T0AndDCAFromCRTSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                                                              recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTSpacePoint> crtSPs) {

    if (tpcTrack.Length() < fMinTrackLength) return std::make_pair(-9999., -9999.);

    matchCand closestSP = GetClosestCRTSpacePoint(detProp, tpcTrack, hits, crtSPs);

    if(closestSP.dca < 0 ) return std::make_pair(-9999., -9999.);
    if (closestSP.dca < fDistanceLimit && (closestSP.dca/closestSP.extrapLen) < fDoverLLimit) return std::make_pair(closestSP.t0, closestSP.dca);

    return std::make_pair(-9999., -9999.);


  }

}
