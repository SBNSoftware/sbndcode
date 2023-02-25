#include "CRTTrackMatchAlg.h"

namespace sbnd{

  CRTTrackMatchAlg::CRTTrackMatchAlg(const Config& config) : CRTTrackMatchAlg(config, lar::providerFrom<geo::Geometry>())
  {}

  CRTTrackMatchAlg::CRTTrackMatchAlg(const Config& config, geo::GeometryCore const* GeometryService){

    this->reconfigure(config);

    fGeometryService = GeometryService;
  
  }


  CRTTrackMatchAlg::CRTTrackMatchAlg(){


  }


  CRTTrackMatchAlg::~CRTTrackMatchAlg(){

  }


  void CRTTrackMatchAlg::reconfigure(const Config& config){

    fMaxAngleDiff = config.MaxAngleDiff();
    fMaxDistance = config.MaxDistance();
    fMaxScore = config.MaxScore();
    fTPCTrackLabel = config.TPCTrackLabel();
    fSelectionMetric = config.SelectionMetric();

    return;

  }
 

  // Calculate intersection between CRT track and TPC (AABB Ray-Box intersection)
  // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
  std::pair<geo::Point_t, geo::Point_t> CRTTrackMatchAlg::TpcIntersection(const geo::TPCGeo& tpcGeo, sbnd::crt::CRTTrack track){

    // Find the intersection between the track and the TPC
    geo::Point_t start = track.Start();
    geo::Point_t end   = track.End();
    geo::Point_t min(tpcGeo.MinX(), tpcGeo.MinY(), tpcGeo.MinZ());
    geo::Point_t max(tpcGeo.MaxX(), tpcGeo.MaxY(), tpcGeo.MaxZ());

    std::pair<geo::Point_t, geo::Point_t> intersection = crt::CRTCommonUtils::CubeIntersection(min, max, start, end);
    return intersection;
  }


  // Function to calculate if a CRTTrack crosses the TPC volume
  bool CRTTrackMatchAlg::CrossesTPC(sbnd::crt::CRTTrack track){

    for(size_t c = 0; c < fGeometryService->Ncryostats(); c++){
      const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(geo::CryostatID(c));
      for(size_t t = 0; t < cryostat.NTPC(); t++){
        const geo::TPCGeo& tpcGeo = cryostat.TPC(geo::TPCID(cryostat.ID(),t));
        std::pair<geo::Point_t, geo::Point_t> intersection = TpcIntersection(tpcGeo, track);
        if(intersection.first.X() != -99999) return true;
      }
    }
    return false;

  } // CRTTrackMatchAlg::CrossesTPC()

  double CRTTrackMatchAlg::T0FromCRTTracks(detinfo::DetectorPropertiesData const& detProp,
                                           recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return T0FromCRTTracks(detProp, tpcTrack, hits, crtTracks);
  }

  double CRTTrackMatchAlg::T0FromCRTTracks(detinfo::DetectorPropertiesData const& detProp,
                                           recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks) {

    std::pair<sbnd::crt::CRTTrack, double> closest;
    if(fSelectionMetric == "angle"){
      closest = ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxAngleDiff) return -99999;
    }
    else if(fSelectionMetric == "dca"){
      closest = ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxDistance) return -99999;
    }
    else{
      closest = ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxScore) return -99999;
    }

    double crtTime = closest.first.Time() * 1e-3; // [us]

    return crtTime;

  }

  int CRTTrackMatchAlg::GetMatchedCRTTrackId(detinfo::DetectorPropertiesData const& detProp,
                                             recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event){
    std::pair<int, double> result = GetMatchedCRTTrackIdAndScore(detProp, tpcTrack, crtTracks, event);
    return result.first;
  }

  // Find the closest valid matching CRT track ID
  int CRTTrackMatchAlg::GetMatchedCRTTrackId(detinfo::DetectorPropertiesData const& detProp,
                                             recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks) {
    std::pair<int, double> result = GetMatchedCRTTrackIdAndScore(detProp, tpcTrack, hits, crtTracks);
    return result.first;
  }

  std::pair<int,double> CRTTrackMatchAlg::GetMatchedCRTTrackIdAndScore(detinfo::DetectorPropertiesData const& detProp,
                                                                       recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return GetMatchedCRTTrackIdAndScore(detProp, tpcTrack, hits, crtTracks);
  }

  // Find the closest valid matching CRT track ID
  std::pair<int,double> CRTTrackMatchAlg::GetMatchedCRTTrackIdAndScore(detinfo::DetectorPropertiesData const& detProp,
                                                                       recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks) {

    std::pair<int, double> null = std::make_pair(-99999, -99999);

    std::pair<sbnd::crt::CRTTrack, double> closest;
    if(fSelectionMetric == "angle"){
      closest = ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxAngleDiff) return null;
    }
    else if(fSelectionMetric == "dca"){
      closest = ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxDistance) return null;
    }
    else{
      closest = ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);
      if(closest.second == -99999 || closest.second > fMaxScore) return null;
    }

    int crt_i = 0;
    for(auto const& track : crtTracks){
      if(track == closest.first) return std::make_pair(crt_i, closest.second);
      crt_i++;
    }

    return null;

  }

  std::vector<sbnd::crt::CRTTrack> CRTTrackMatchAlg::AllPossibleCRTTracks(detinfo::DetectorPropertiesData const& detProp,
                                                                          recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);
  }


  // Get all CRT tracks that cross the right TPC within an allowed time
  std::vector<sbnd::crt::CRTTrack> CRTTrackMatchAlg::AllPossibleCRTTracks(detinfo::DetectorPropertiesData const& detProp,
                                                                          recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks) {

    std::vector<sbnd::crt::CRTTrack> trackCandidates;

    // Get the hits associated with the tpc track

    // Get the drift direction (0 for stitched tracks)
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    // Get the TPC Geo object from the tpc track
    geo::TPCID tpcID = hits[0]->WireID().asTPCID();
    const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);

    // Loop over the crt tracks
    for(auto const& crtTrack : crtTracks){
      // Calculate the intersection points for that TPC
      std::pair<geo::Point_t, geo::Point_t> intersection = TpcIntersection(tpcGeo, crtTrack);

      // Skip if it doesn't intersect
      if(intersection.first.X() == -99999) continue;

      // Shift the track to the CRT track
      double crtTime = crtTrack.Time() * 1e-3; // [us]
      double shift = driftDirection * crtTime * detProp.DriftVelocity();
      geo::Point_t start = tpcTrack.Vertex();
      geo::Point_t end = tpcTrack.End();
      start.SetX(start.X() + shift);
      end.SetX(end.X() + shift);

      // Check the track is fully contained in the TPC
      if(!TPCGeoUtil::InsideTPC(start, tpcGeo, 2.) && shift != 0) continue;
      if(!TPCGeoUtil::InsideTPC(end, tpcGeo, 2.) && shift != 0) continue;

      trackCandidates.push_back(crtTrack);
    
    }

    return trackCandidates;
  }

  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const& detProp,
                                                                                  recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event, double minDCA){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks, minDCA);
  }

  // Find the closest matching crt track by angle between tracks within angle and DCA limits
  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const& detProp,
                                                                                  recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks, double minDCA){

    // Get the drift direction (0 for stitched tracks)
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    std::vector<sbnd::crt::CRTTrack> possTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<std::pair<sbnd::crt::CRTTrack, double>> candidates;
    for(auto const& possTrack : possTracks){
      double angle = AngleBetweenTracks(tpcTrack, possTrack);

      if(minDCA != -1){
        if(minDCA == 0) minDCA = fMaxDistance;
        double crtTime = possTrack.Time() * 1e-3; // [us]
        double shift = driftDirection * crtTime * detProp.DriftVelocity();
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
    sbnd::crt::CRTTrack track;
    return std::make_pair(track, -99999);
  }

  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const& detProp,
                                                                                recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event, double minAngle) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks, minAngle);
  }

  // Find the closest matching crt track by average DCA between tracks within angle and DCA limits
  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const& detProp,
                                                                                recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks,  double minAngle){

    // Get the drift direction (0 for stitched tracks)
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    std::vector<sbnd::crt::CRTTrack> possTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<std::pair<sbnd::crt::CRTTrack, double>> candidates;
    for(auto const& possTrack : possTracks){

      double crtTime = possTrack.Time() * 1e-3; // [us]
      double shift = driftDirection * crtTime * detProp.DriftVelocity();

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
    sbnd::crt::CRTTrack track;
    return std::make_pair(track, -99999);

  }


  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const& detProp,
                                                                                  recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);
  }

  // Find the closest matching crt track by average DCA between tracks within angle and DCA limits
  std::pair<sbnd::crt::CRTTrack, double> CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const& detProp,
                                                                                  recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks){

    // Get the drift direction (0 for stitched tracks)
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    std::vector<sbnd::crt::CRTTrack> possTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<std::pair<sbnd::crt::CRTTrack, double>> candidates;
    for(auto const& possTrack : possTracks){

      double crtTime = possTrack.Time() * 1e-3; // [us]
      double shift = driftDirection * crtTime * detProp.DriftVelocity();

      double DCA = AveDCABetweenTracks(tpcTrack, possTrack, shift);
      double angle = AngleBetweenTracks(tpcTrack, possTrack);
      double score = DCA + 4*180/TMath::Pi()*angle;

      candidates.push_back(std::make_pair(possTrack, score));
    }

    std::sort(candidates.begin(), candidates.end(), [](auto& left, auto& right){
        return left.second < right.second;});

    if(candidates.size() > 0){
      return candidates[0];
    }
    sbnd::crt::CRTTrack track;
    return std::make_pair(track, -99999);

  }


  // Calculate the angle between tracks assuming start is at the largest Y
  double CRTTrackMatchAlg::AngleBetweenTracks(recob::Track tpcTrack, sbnd::crt::CRTTrack crtTrack){

    // Calculate the angle between the tracks
    geo::Point_t crtStart = crtTrack.Start();
    geo::Point_t crtEnd   = crtTrack.End();
    if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);

    geo::Point_t tpcStart = tpcTrack.Vertex();
    geo::Point_t tpcEnd   = tpcTrack.End();
    if(tpcStart.Y() < tpcEnd.Y()) std::swap(tpcStart, tpcEnd);

    return TMath::ACos((tpcStart - tpcEnd).Dot(crtStart - crtEnd) / ((tpcStart - tpcEnd).R() * (crtStart - crtEnd).R()));

  }


  // Calculate the average DCA between tracks
  double CRTTrackMatchAlg::AveDCABetweenTracks(recob::Track tpcTrack, sbnd::crt::CRTTrack crtTrack, double shift){

    geo::Point_t crtStart = crtTrack.Start();
    geo::Point_t crtEnd   = crtTrack.End();
    if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
    double denominator = (crtEnd - crtStart).R();

    size_t npts = tpcTrack.NumberTrajectoryPoints();

    double aveDCA = 0;
    int usedPts = 0;
    for(size_t i = 0; i < npts; i++){
      geo::Point_t point = tpcTrack.LocationAtPoint(i);

      // Pandora produces dummy points
      if(!tpcTrack.HasValidPoint(i)) continue;

      point.SetX(point.X() + shift);
      aveDCA += (point - crtStart).Cross(point - crtEnd).R()/denominator;
      usedPts++;
    }

    return aveDCA/usedPts;

  }

  double CRTTrackMatchAlg::AveDCABetweenTracks(detinfo::DetectorPropertiesData const& detProp,
                                               recob::Track tpcTrack, sbnd::crt::CRTTrack crtTrack, const art::Event& event) {
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    return AveDCABetweenTracks(detProp, tpcTrack, hits, crtTrack);
  }


  // Calculate the average DCA between tracks
  double CRTTrackMatchAlg::AveDCABetweenTracks(detinfo::DetectorPropertiesData const& detProp,
                                               recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, sbnd::crt::CRTTrack crtTrack) {

    // Get the drift direction (0 for stitched tracks)
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    double crtTime     = crtTrack.Time() * 1e-3; // [us]
    double shift       = driftDirection * crtTime * detProp.DriftVelocity();

    geo::Point_t crtStart = crtTrack.Start();
    geo::Point_t crtEnd   = crtTrack.End();
    if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
    double denominator = (crtEnd - crtStart).R();

    size_t npts = tpcTrack.NumberTrajectoryPoints();

    double aveDCA = 0;
    int usedPts = 0;
    for(size_t i = 0; i < npts; i++){
      geo::Point_t point = tpcTrack.LocationAtPoint(i);

      // Pandora produces dummy points
      if(!tpcTrack.HasValidPoint(i)) continue;

      point.SetX(point.X() + shift);
      aveDCA += (point - crtStart).Cross(point - crtEnd).R()/denominator;
      usedPts++;
    }

    return aveDCA/usedPts;

  }


}
