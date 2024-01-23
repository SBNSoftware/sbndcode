#include "CRTSpacePointMatchAlg.h"

namespace sbnd::crt {

  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg(const Config& config)
    : CRTSpacePointMatchAlg(config, lar::providerFrom<geo::Geometry>())
  {}

  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg(const Config& config, geo::GeometryCore const *GeometryService)
  {
    this->reconfigure(config);
    fGeometryService = GeometryService;
  }

  CRTSpacePointMatchAlg::CRTSpacePointMatchAlg(){}

  CRTSpacePointMatchAlg::~CRTSpacePointMatchAlg(){}

  void CRTSpacePointMatchAlg::reconfigure(const Config& config)
  {
    fTrackDirectionFrac = config.TrackDirectionFrac();
    fDCALimit           = config.DCALimit();
    fTimeCorrection     = config.TimeCorrection();
    fMinTPCTrackLength  = config.MinTPCTrackLength();
    fDirMethod          = config.DirMethod();
    fDCAuseBox          = config.DCAuseBox();
    fDCAoverLength      = config.DCAoverLength();
    fPECut              = config.PECut();
    fMaxUncert          = config.MaxUncert();
    fTPCTrackLabel      = config.TPCTrackLabel();
    fCRTSpacePointLabel = config.CRTSpacePointLabel();

    return;
  }

  SPMatchCandidate CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                                  const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const art::Event &e)
  {
    art::Handle<std::vector<recob::Track>> trackHandle;
    e.getByLabel(fTPCTrackLabel, trackHandle);

    const art::FindManyP<recob::Hit> tracksToHits(trackHandle, e, fTPCTrackLabel);
    const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(track.key());

    return GetClosestCRTSpacePoint(detProp, track, hits, crtSPs, e);
  }

  SPMatchCandidate CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                                  const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs,
                                                                  const art::Event &e)
  {
    const geo::Point_t start = track->Vertex();
    const geo::Point_t end   = track->End();

    const int driftDirection                = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    const std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);

    const std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTSpacePoint(detProp, track, t0MinMax, crtSPs, driftDirection, e);
  }

  SPMatchCandidate CRTSpacePointMatchAlg::GetClosestCRTSpacePoint(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                                  const std::pair<double, double> t0MinMax, const std::vector<art::Ptr<CRTSpacePoint>> &crtSPs, const int driftDirection,
                                                                  const art::Event &e)
  {
    if(track->Length() < fMinTPCTrackLength)
      return SPMatchCandidate();

    const geo::Point_t start = track->Vertex();
    const geo::Point_t end   = track->End();

    std::vector<SPMatchCandidate> candidates;

    for(auto &crtSP : crtSPs){

      const double crtTime = crtSP->Time() * 1e-3 + fTimeCorrection;

      if(!(crtTime > t0MinMax.first - 10. && crtTime < t0MinMax.second + 10.))
        continue;

      if(crtSP->PE() < fPECut || crtSP->XErr() > fMaxUncert || crtSP->YErr() > fMaxUncert || crtSP->ZErr() > fMaxUncert)
        continue;

      const geo::Point_t crtPoint = crtSP->Pos();
  
      std::pair<geo::Vector_t, geo::Vector_t> startEndDir;

      if(fDirMethod==2)
        startEndDir = AverageTrackDirections(track, fTrackDirectionFrac);
      else
        startEndDir = TrackDirections(detProp, track, fTrackDirectionFrac, crtTime, driftDirection);

      const geo::Vector_t startDir = startEndDir.first;
      const geo::Vector_t endDir   = startEndDir.second;
    
      const double startDCA = DistOfClosestApproach(detProp, start, startDir, crtSP, driftDirection, crtTime, e);
      const double endDCA   = DistOfClosestApproach(detProp, end, endDir, crtSP, driftDirection, crtTime, e);

      const double xshift = driftDirection * crtTime * detProp.DriftVelocity();

      geo::Point_t thisstart = start;
      thisstart.SetX(start.X()+xshift);
      geo::Point_t thisend = end;
      thisend.SetX(end.X()+xshift);

      if(startDCA < fDCALimit || endDCA < fDCALimit)
        {
          const double distS = (crtPoint - thisstart).R();
          const double distE = (crtPoint - thisend).R();

          const double scoreS = fDCAoverLength ? startDCA / distS : startDCA;
          const double scoreE = fDCAoverLength ? endDCA / distE : endDCA;

          if(distS < distE && startDCA < fDCALimit)
            candidates.emplace_back(crtSP, track, crtTime, scoreS, true);
          else if(endDCA < fDCALimit)
            candidates.emplace_back(crtSP, track, crtTime, scoreE, true);
        }
    }

    if(candidates.size() == 0)
      return SPMatchCandidate();

    std::sort(candidates.begin(), candidates.end(),
              [](const SPMatchCandidate &a, const SPMatchCandidate &b)
              {
                if(a.score < 0 && b.score > 0)
                  return false;
                else if(a.score > 0 && b.score < 0)
                  return true;
                else
                  return a.score < b.score;
              });

    return candidates[0];
  }

  std::pair<double, double> CRTSpacePointMatchAlg::TrackT0Range(detinfo::DetectorPropertiesData const &detProp, const double startX,
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

  double CRTSpacePointMatchAlg::DistOfClosestApproach(detinfo::DetectorPropertiesData const &detProp, geo::Point_t trackStart,
                                                      const geo::Vector_t &trackDir, const art::Ptr<CRTSpacePoint> &crtSP,
                                                      const int driftDirection, const double t0, const art::Event &e)
  {
    const double xshift = driftDirection* t0 * detProp.DriftVelocity();
    trackStart.SetX(trackStart.X() + xshift);

    const geo::Point_t end = trackStart + trackDir;

    art::Handle<std::vector<CRTSpacePoint>> spacePointHandle;
    e.getByLabel(fCRTSpacePointLabel, spacePointHandle);

    const art::FindOneP<CRTCluster> spacePointsToClusters(spacePointHandle, e, fCRTSpacePointLabel);
    const art::Ptr<CRTCluster> cluster = spacePointsToClusters.at(crtSP.key());

    if(fDCAuseBox)
      return CRTCommonUtils::DistToCRTSpacePoint(crtSP, trackStart, end, cluster->Tagger());
    else
      return CRTCommonUtils::SimpleDCA(crtSP, trackStart, trackDir);
  }

  std::pair<geo::Vector_t, geo::Vector_t> CRTSpacePointMatchAlg::AverageTrackDirections(const art::Ptr<recob::Track> &track, const double frac)
  {
    const unsigned N          = track->NumberTrajectoryPoints();
    const unsigned NValid     = track->CountValidPoints();
    const unsigned NValidFrac = std::floor(NValid * frac);

    geo::Vector_t forwardDir(0, 0, 0), backwardDir(0, 0, 0);

    unsigned iValidForward = 0, iValidBackward = 0;

    for(unsigned i = 0; i < N; ++i)
      {
        if(track->HasValidPoint(i) && iValidForward < NValidFrac)
          {
            forwardDir += track->DirectionAtPoint(i);
            ++iValidForward;
          }

        if(track->HasValidPoint(N - (i + 1)) && iValidBackward < NValidFrac)
          {
            backwardDir += track->DirectionAtPoint(N - (i + 1));
            ++iValidBackward;
          }
      }

    forwardDir  /= NValidFrac;
    backwardDir /= NValidFrac;

    return std::make_pair(forwardDir, backwardDir);

  }

  std::pair<geo::Vector_t, geo::Vector_t> CRTSpacePointMatchAlg::TrackDirections(detinfo::DetectorPropertiesData const &detProp, const art::Ptr<recob::Track> &track,
                                                                                 const double frac, const double CRTtime, const int driftDirection)
  {
    const unsigned N    = track->NumberTrajectoryPoints();
    const unsigned NMid = std::floor(N * frac);

    geo::Point_t start = track->Start();
    geo::Point_t end   = track->End();
    geo::Point_t mid   = track->LocationAtPoint(NMid);

    const double xshift = driftDirection * CRTtime * detProp.DriftVelocity();

    start.SetX(start.X() + xshift);
    end.SetX(end.X() + xshift);
    mid.SetX(mid.X() + xshift);

    const geo::Vector_t startDir = (mid - start).Unit();
    const geo::Vector_t endDir   = (mid - end).Unit();

    return std::make_pair(startDir, endDir);
  }
}
