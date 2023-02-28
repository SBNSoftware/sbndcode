#include "CRTTrackMatchAlg.h"

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
    fMaxAngleDiff    = config.MaxAngleDiff();
    fMaxDCA          = config.MaxDCA();
    fMaxScore        = config.MaxScore();
    fTPCTrackLabel   = config.TPCTrackLabel();
    fSelectionMetric = config.SelectionMetric();

    return;
  }
 
  MatchCandidate CRTTrackMatchAlg::GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return GetBestMatchedCRTTrack(detProp, tpcTrack, hits, crtTracks);
  }

  MatchCandidate CRTTrackMatchAlg::GetBestMatchedCRTTrack(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks)
  {
    if(fSelectionMetric == "angle")
      {
        MatchCandidate candidate = ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks, fMaxDCA);

        if(candidate.score > fMaxAngleDiff)
          return MatchCandidate();
        else
          return candidate;
      }
    else if(fSelectionMetric == "dca")
      {
        MatchCandidate candidate = ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks, fMaxAngleDiff);

        if(candidate.score > fMaxDCA)
          return MatchCandidate();
        else
          return candidate;
      }
    else
      {
        MatchCandidate candidate = ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);

        if(candidate.score > fMaxScore)
          return MatchCandidate();
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

    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const geo::TPCID tpcID    = hits[0]->WireID().asTPCID();
    const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);

    for(auto const &crtTrack : crtTracks)
      {
        geo::Point_t entry, exit;

        if(!TPCIntersection(tpcGeo, crtTrack, entry, exit))
          continue;

        const double crtTime = crtTrack->Time() * 1e-3;
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

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxDCA)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByAngle(detProp, tpcTrack, hits, crtTracks, maxDCA);
  }

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByAngle(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks, const double maxDCA)
  {
    if(maxDCA == - 1.)
      return MatchCandidate();

    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<MatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = possCRTTrack->Time() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        if(DCA > maxDCA)
          continue;

        candidates.emplace_back(possCRTTrack, crtTime, angle, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const MatchCandidate &a, const MatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return MatchCandidate();
  }

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack, 
                                                        const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e, const double maxAngle)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByDCA(detProp, tpcTrack, hits, crtTracks, maxAngle);
  }

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByDCA(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                        const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks,  const double maxAngle)
  {
    if(maxAngle == -1.)
      return MatchCandidate();

    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<MatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = possCRTTrack->Time() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        if(angle > maxAngle)
          continue;
        candidates.emplace_back(possCRTTrack, crtTime, DCA, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const MatchCandidate &a, const MatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return MatchCandidate();
  }

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<CRTTrack>> &crtTracks, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> tpcTrackHandle;
    e.getByLabel(fTPCTrackLabel, tpcTrackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(tpcTrackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(tpcTrack.key());

    return ClosestCRTTrackByScore(detProp, tpcTrack, hits, crtTracks);
  }

  MatchCandidate CRTTrackMatchAlg::ClosestCRTTrackByScore(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &tpcTrack,
                                                          const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTTrack>> &crtTracks)
  {
    const int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);

    const std::vector<art::Ptr<CRTTrack>> possCRTTracks = AllPossibleCRTTracks(detProp, tpcTrack, hits, crtTracks);

    std::vector<MatchCandidate> candidates;

    for(auto const &possCRTTrack : possCRTTracks)
      {
        const double crtTime = possCRTTrack->Time() * 1e-3;
        const double shift   = driftDirection * crtTime * detProp.DriftVelocity();
        const double DCA     = AveDCABetweenTracks(tpcTrack, possCRTTrack, shift);
        const double angle   = AngleBetweenTracks(tpcTrack, possCRTTrack);

        const double score   = DCA + 4 * 180 / TMath::Pi() * angle;

        candidates.emplace_back(possCRTTrack, crtTime, score, true);
      }

    std::sort(candidates.begin(), candidates.end(),
              [](const MatchCandidate &a, const MatchCandidate &b)
              { return a.score < b.score; });

    if(candidates.size() > 0)
      return candidates[0];
    else
      return MatchCandidate();
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
}
