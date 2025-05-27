#include "CRTTrackMatchAlg.h"

#include "larcore/Geometry/Geometry.h"

namespace sbnd::crt {

  CRTTrackMatchAlg::CRTTrackMatchAlg(const Config& config)
    : CRTTrackMatchAlg(config, lar::providerFrom<geo::Geometry>())
  {}

  CRTTrackMatchAlg::CRTTrackMatchAlg(const Config& config, geo::GeometryCore const *GeometryService)
  {
    this->reconfigure(config);
    fGeometryService = GeometryService;
  }

  CRTTrackMatchAlg::CRTTrackMatchAlg(){}

  CRTTrackMatchAlg::~CRTTrackMatchAlg(){}

  void CRTTrackMatchAlg::reconfigure(const Config& config)
  {
    fMaxAngleDiff      = config.MaxAngleDiff();
    fMaxDCA            = config.MaxDCA();
    fMaxScore          = config.MaxScore();
    fMinTPCTrackLength = config.MinTPCTrackLength();
    fTPCTrackLabel     = config.TPCTrackLabel();
    fSelectionMetric   = config.SelectionMetric();
    fUseTs0            = config.UseTs0();

    return;
  }
 
  TrackMatchCandidate CRTTrackMatchAlg::GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return GetBestMatchedCRTTrack(detProp, tpcTrack, hits, crtTracks);
  }

  TrackMatchCandidate CRTTrackMatchAlg::GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks)
  {
    if(tpcTrack->Length() < fMinTPCTrackLength)
      return TrackMatchCandidate();

    if(fSelectionMetric == "angle")
      {
        TrackMatchCandidate candidate = ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks, fMaxDCA);

        if(candidate.score > fMaxAngleDiff)
          return TrackMatchCandidate();
        else
          return candidate;
      }
    else if(fSelectionMetric == "dca")
      {
        TrackMatchCandidate candidate = ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks, fMaxAngleDiff);

        if(candidate.score > fMaxDCA)
          return TrackMatchCandidate();
        else
          return candidate;
      }
    else
      {
        TrackMatchCandidate candidate = ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);

        if(candidate.score > fMaxScore)
          return TrackMatchCandidate();
        else
          return candidate;
      }
  }

  bool CRTTrackMatchAlg::TPCIntersection(const geo::TPCGeo &tpcGeo, const art::Ptr<CRTTrack> &track, geo::Point_t &entry, geo::Point_t &exit)
  {
    const geo::Point_t start = track->Start();
    const geo::Point_t end   = track->End();

    const geo::Point_t min(tpcGeo.MinX(), tpcGeo.MinY(), tpcGeo.MinZ());
    const geo::Point_t max(tpcGeo.MaxX(), tpcGeo.MaxY(), tpcGeo.MaxZ());

    return CRTCommonUtils::CuboidIntersection(min, max, start, end, entry, exit);
  }

  std::vector<art::Ptr<CRTTrack>> CRTTrackMatchAlg::AllPossibleCRTTracks(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                                         const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);
  }

  std::vector<art::Ptr<CRTTrack>> CRTTrackMatchAlg::AllPossibleCRTTracks(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                                         const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks)
  {
    std::vector<art::Ptr<CRTTrack>> candidates;

    const int driftDirection                 = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    const std::pair<double, double> xLimits  = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    const std::pair<double, double> t0MinMax = TrackT0Range(detProp, tpcTrack->Vertex().X(), tpcTrack->End().X(), driftDirection, xLimits);

    const geo::TPCID tpcID    = hits[0]->WireID().asTPCID();
    const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);

    for(auto const &crtTrack : crtTracks)
      {
        geo::Point_t entry, exit;

        if(!TPCIntersection(tpcGeo, crtTrack, entry, exit))
          continue;

        const double crtTime = fUseTs0 ? crtTrack->Ts0() * 1e-3 : crtTrack->Ts1() * 1e-3;

        if(!(crtTime > t0MinMax.first - 10. && crtTime < t0MinMax.second + 10.))
          continue;

        const double shift = driftDirection * crtTime * detProp.DriftVelocity();

        geo::Point_t start = tpcTrack->Vertex();
        geo::Point_t end   = tpcTrack->End();
        start.SetX(start.X() + shift);
        end.SetX(end.X() + shift);

        if(!TPCGeoUtil::InsideTPC(start, tpcGeo, 2.) && shift != 0)
          continue;

        if(!TPCGeoUtil::InsideTPC(end, tpcGeo, 2.) && shift != 0)
          continue;

        candidates.push_back(crtTrack);
      }

    return candidates;
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxDCA)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks, maxDCA);
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks, const double maxDCA)
  {
    if(maxDCA == - 1.)
      return TrackMatchCandidate();

    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<TrackMatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = fUseTs0 ? possCRTTrack->Ts0() * 1e-3 : possCRTTrack->Ts1() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        if(DCA > maxDCA)
          continue;

        candidates.emplace_back(possCRTTrack, tpcTrack, crtTime, angle, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const TrackMatchCandidate &a, const TrackMatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return TrackMatchCandidate();
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                             const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxAngle)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks, maxAngle);
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                             const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks,  const double maxAngle)
  {
    if(maxAngle == -1.)
      return TrackMatchCandidate();

    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<TrackMatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = fUseTs0 ? possCRTTrack->Ts0() * 1e-3 : possCRTTrack->Ts1() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        if(angle > maxAngle)
          continue;
        candidates.emplace_back(possCRTTrack, tpcTrack, crtTime, DCA, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const TrackMatchCandidate &a, const TrackMatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return TrackMatchCandidate();
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);
  }

  TrackMatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                               const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks)
  {
    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<TrackMatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = fUseTs0 ? possCRTTrack->Ts0() * 1e-3 : possCRTTrack->Ts1() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        const double score   = DCA + 4 * 180 / TMath::Pi() * angle;

        candidates.emplace_back(possCRTTrack, tpcTrack, crtTime, score, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const TrackMatchCandidate &a, const TrackMatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return TrackMatchCandidate();
  }

  double CRTTrackMatchAlg::AngleBetweenTracks(const art::Ptr<recob::Track> &tpcTrack, const art::Ptr<CRTTrack> &crtTrack)
  {
    geo::Point_t crtStart = crtTrack->Start();
    geo::Point_t crtEnd   = crtTrack->End();
    if(crtStart.Y() < crtEnd.Y())
      std::swap(crtStart, crtEnd);

    geo::Point_t tpcStart = tpcTrack->Vertex();
    geo::Point_t tpcEnd   = tpcTrack->End();
    if(tpcStart.Y() < tpcEnd.Y())
      std::swap(tpcStart, tpcEnd);

    double angle = TMath::ACos((tpcStart - tpcEnd).Dot(crtStart - crtEnd) / ((tpcStart - tpcEnd).R() * (crtStart - crtEnd).R()));

    if(angle > TMath::Pi()/2. && angle < TMath::Pi())
      angle = TMath::Pi() - angle;

    return angle;
  }

  double CRTTrackMatchAlg::AveDCABetweenTracks(const art::Ptr<recob::Track> &tpcTrack, const art::Ptr<CRTTrack> &crtTrack, const double shift)
  {
    geo::Point_t crtStart = crtTrack->Start();
    geo::Point_t crtEnd   = crtTrack->End();
    if(crtStart.Y() < crtEnd.Y())
      std::swap(crtStart, crtEnd);

    const double denominator = (crtEnd - crtStart).R();
    const unsigned N         = tpcTrack->NumberTrajectoryPoints();

    double aveDCA = 0;
    int usedPts = 0;

    for(unsigned i = 0; i < N; ++i)
      {
        geo::Point_t point = tpcTrack->LocationAtPoint(i);

        if(!tpcTrack->HasValidPoint(i))
          continue;

        point.SetX(point.X() + shift);
        aveDCA += (point - crtStart).Cross(point - crtEnd).R() / denominator;
        ++usedPts;
      }

    return aveDCA / usedPts;
  }

  std::pair<double, double> CRTTrackMatchAlg::TrackT0Range(detinfo::DetectorPropertiesData const &detProp, const double startX,
                                                           const double endX, const int driftDirection, const std::pair<double, double> xLimits)
  {
    if(driftDirection == 0)
      return std::make_pair(0, 0);

    const double driftVel = driftDirection * detProp.DriftVelocity();

    const double maxX     = std::max(startX, endX);
    const double maxLimit = std::max(xLimits.first, xLimits.second);
    const double maxShift = maxLimit - maxX;

    const double minX     = std::min(startX, endX);
    const double minLimit = std::min(xLimits.first, xLimits.second);
    const double minShift = minLimit - minX;

    const double t0max = maxShift / driftVel;
    const double t0min = minShift / driftVel;

    return std::make_pair(std::min(t0min, t0max), std::max(t0min, t0max));
  }
}
