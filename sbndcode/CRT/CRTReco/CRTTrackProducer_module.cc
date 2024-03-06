////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackProducer
// Plugin Type: producer
// File:        CRTTrackProducer_module.cc
//
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "Eigen/Dense"

#include "TMath.h"

#include <memory>


namespace geo {
  using Point2D_t = ROOT::Math::PositionVector2D<ROOT::Math::Cartesian2D<double>,
                                                 ROOT::Math::GlobalCoordinateSystemTag>;

  using Vector2D_t = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>,
                                                      ROOT::Math::GlobalCoordinateSystemTag>;
}

namespace sbnd::crt {
  class CRTTrackProducer;
}


class sbnd::crt::CRTTrackProducer : public art::EDProducer {
public:
  explicit CRTTrackProducer(fhicl::ParameterSet const& p);

  CRTTrackProducer(CRTTrackProducer const&) = delete;
  CRTTrackProducer(CRTTrackProducer&&) = delete;
  CRTTrackProducer& operator=(CRTTrackProducer const&) = delete;
  CRTTrackProducer& operator=(CRTTrackProducer&&) = delete;

  void produce(art::Event& e) override;

  void OrderSpacePoints(std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec);

  std::vector<std::pair<CRTTrack, std::set<unsigned>>> CreateTrackCandidates(const std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec,
                                                                             const art::FindOneP<CRTCluster> &spacePointsToCluster);

  void TimeErrorCalculator(const std::vector<double> &times, double &mean, double &err);

  void OrderTrackCandidates(std::vector<std::pair<CRTTrack, std::set<unsigned>>> &trackCandidates);

  std::vector<std::pair<CRTTrack, std::set<unsigned>>> ChoseTracks(std::vector<std::pair<CRTTrack, std::set<unsigned>>> &trackCandidates);

  double DistanceOfClosestApproach(const CRTTagger tagger, const art::Ptr<CRTSpacePoint> &spacePoint,
                                   const geo::Point_t &start, const geo::Vector_t &dir);

  bool IsPointInsideBox(const geo::Point_t &point, const geo::Point_t &centre, const geo::Point_t &widths);

  double DistanceOfClosestApproach(const CoordSet &constrainedPlane, const geo::Point_t &point,
                                   const geo::Point_t &centre, const geo::Point_t &widths);

  double MinimumApproach(const double &x, const double &dx, const double &y, const double &dy, const geo::Point2D_t &p);

  double DistanceOfClosestApproach(const geo::Point2D_t &v1, const geo::Point2D_t &v2, const geo::Point2D_t &p);

  void BestFitLine(const geo::Point_t &a, const geo::Point_t &b, const geo::Point_t &c, const CRTTagger &primary_tagger,
                   const CRTTagger &secondary_tagger, const CRTTagger &tertiary_tagger, geo::Point_t &start,
                   geo::Point_t &mid, geo::Point_t &end, double &gof);

  geo::Point_t LineTaggerIntersectionPoint(const geo::Point_t &start, const geo::Vector_t &dir, const CRTTagger &tagger);

  double TripleTrackToF(std::vector<double> times);

private:

  CRTGeoAlg   fCRTGeoAlg;
  std::string fCRTSpacePointModuleLabel;
  double      fCoincidenceTimeRequirement;
  double      fThirdSpacePointMaximumDCA;
};


sbnd::crt::CRTTrackProducer::CRTTrackProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<double>("CoincidenceTimeRequirement"))
  , fThirdSpacePointMaximumDCA(p.get<double>("ThirdSpacePointMaximumDCA"))
  {
    produces<std::vector<CRTTrack>>();
    produces<art::Assns<CRTSpacePoint, CRTTrack>>();
  }

void sbnd::crt::CRTTrackProducer::produce(art::Event& e)
{
  auto trackVec            = std::make_unique<std::vector<CRTTrack>>();
  auto trackSpacePointAssn = std::make_unique<art::Assns<CRTSpacePoint, CRTTrack>>();
  
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::FindOneP<CRTCluster> spacePointsToCluster(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  OrderSpacePoints(CRTSpacePointVec);

  std::vector<std::pair<CRTTrack, std::set<unsigned>>> trackCandidates = CreateTrackCandidates(CRTSpacePointVec, spacePointsToCluster);

  OrderTrackCandidates(trackCandidates);

  std::vector<std::pair<CRTTrack, std::set<unsigned>>> chosenTracks = ChoseTracks(trackCandidates);

  for(auto const& [track, spIDs] : chosenTracks)
    {
      trackVec->push_back(track);

      for(auto const& spID : spIDs)
        util::CreateAssn(*this, e, *trackVec, CRTSpacePointVec[spID], *trackSpacePointAssn);
    }

  e.put(std::move(trackVec));
  e.put(std::move(trackSpacePointAssn));
}

void sbnd::crt::CRTTrackProducer::OrderSpacePoints(std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec)
{
  std::sort(spacePointVec.begin(), spacePointVec.end(), 
            [](const art::Ptr<CRTSpacePoint> &a, const art::Ptr<CRTSpacePoint> &b) -> bool {
              return a->Time() < b->Time();
            });
}

std::vector<std::pair<sbnd::crt::CRTTrack, std::set<unsigned>>> sbnd::crt::CRTTrackProducer::CreateTrackCandidates(const std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec,
                                                                                                                   const art::FindOneP<CRTCluster> &spacePointsToCluster)
{
  std::vector<std::pair<CRTTrack, std::set<unsigned>>> candidates;

  for(unsigned i = 0; i < spacePointVec.size(); ++i)
    {
      const art::Ptr<CRTSpacePoint> primarySpacePoint = spacePointVec[i];
      const art::Ptr<CRTCluster> primaryCluster       = spacePointsToCluster.at(primarySpacePoint.key());

      for(unsigned ii = i+1; ii < spacePointVec.size(); ++ii)
        {
          const art::Ptr<CRTSpacePoint> secondarySpacePoint = spacePointVec[ii];
          const art::Ptr<CRTCluster> secondaryCluster       = spacePointsToCluster.at(secondarySpacePoint.key());

          if(secondarySpacePoint->Time() - primarySpacePoint->Time() > fCoincidenceTimeRequirement ||
             secondaryCluster->Tagger() == primaryCluster->Tagger())
            continue;

          const geo::Point_t &start = primarySpacePoint->Pos();
          const geo::Point_t &end   = secondarySpacePoint->Pos();
          const geo::Vector_t &dir  = (end - start).Unit();

          if(CRTCommonUtils::IsTopTagger(primaryCluster->Tagger()) || CRTCommonUtils::IsTopTagger(secondaryCluster->Tagger()))
            {
              for(unsigned iii = ii + 1; iii < spacePointVec.size(); ++iii)
                {
                  const art::Ptr<CRTSpacePoint> tertiarySpacePoint = spacePointVec[iii];
                  const art::Ptr<CRTCluster> tertiaryCluster       = spacePointsToCluster.at(tertiarySpacePoint.key());

                  if(!CRTCommonUtils::CoverTopTaggers(primaryCluster->Tagger(), secondaryCluster->Tagger(), tertiaryCluster->Tagger()))
                    continue;

                  if(tertiarySpacePoint->Time() - primarySpacePoint->Time() > fCoincidenceTimeRequirement ||
                     tertiaryCluster->Tagger() == primaryCluster->Tagger())
                    continue;

                  if(tertiarySpacePoint->Time() - secondarySpacePoint->Time() > fCoincidenceTimeRequirement ||
                     tertiaryCluster->Tagger() == secondaryCluster->Tagger())
                    continue;

                  const CRTTagger &tertiaryTagger = tertiaryCluster->Tagger();

                  const double dca = DistanceOfClosestApproach(tertiaryTagger, tertiarySpacePoint, start, dir);

                  if(dca < fThirdSpacePointMaximumDCA)
                    {
                      double time, etime;
                      const std::vector<double> times = {primarySpacePoint->Time(), secondarySpacePoint->Time(), tertiarySpacePoint->Time()};
                      TimeErrorCalculator(times, time, etime);
                      const double tof = TripleTrackToF(times);

                      const double pe = primarySpacePoint->PE() + secondarySpacePoint->PE() + tertiarySpacePoint->PE();

                      const std::set<CRTTagger> used_taggers = {primaryCluster->Tagger(), secondaryCluster->Tagger(), tertiaryCluster->Tagger()};
 
                      geo::Point_t fitStart, fitMid, fitEnd;
                      double gof;
                      
                      BestFitLine(primarySpacePoint->Pos(), secondarySpacePoint->Pos(), tertiarySpacePoint->Pos(), primaryCluster->Tagger(), 
                                  secondaryCluster->Tagger(), tertiaryCluster->Tagger(), fitStart, fitMid, fitEnd, gof);

                      const CRTTrack track({fitStart, fitMid, fitEnd}, time, etime, pe, tof, used_taggers);
                      const std::set<unsigned> used_spacepoints = {i, ii, iii};

                      candidates.emplace_back(track, used_spacepoints);
                    }
                }
            }

          double time, etime;
          TimeErrorCalculator({primarySpacePoint->Time(), secondarySpacePoint->Time()}, time, etime);
          const double tof = secondarySpacePoint->Time() - primarySpacePoint->Time();

          const double pe = primarySpacePoint->PE() + secondarySpacePoint->PE();

          const std::set<CRTTagger> used_taggers = {primaryCluster->Tagger(), secondaryCluster->Tagger()};

          const CRTTrack track(start, end, time, etime, pe, tof, used_taggers);
          const std::set<unsigned> used_spacepoints = {i, ii};

          candidates.emplace_back(track, used_spacepoints);
        }
    }
  return candidates;
}

void sbnd::crt::CRTTrackProducer::TimeErrorCalculator(const std::vector<double> &times, double &mean, double &err)
{
  double sum = 0.;
  for(auto const &time : times)
    sum += time;

  mean = sum / times.size();

  double summed_var = 0.;
  for(auto const &time : times)
    summed_var += std::pow((time - mean), 2);

  err = std::sqrt(summed_var / times.size());
}

void sbnd::crt::CRTTrackProducer::OrderTrackCandidates(std::vector<std::pair<CRTTrack, std::set<unsigned>>> &trackCandidates)
{
  std::sort(trackCandidates.begin(), trackCandidates.end(),
            [](const std::pair<CRTTrack, std::set<unsigned>> &a, const std::pair<CRTTrack, std::set<unsigned>> &b) -> bool {
              if(a.second.size() != b.second.size())
                return a.second.size() > b.second.size();
              else
                return a.first.TimeErr() < b.first.TimeErr();
            });
}

std::vector<std::pair<sbnd::crt::CRTTrack, std::set<unsigned>>> sbnd::crt::CRTTrackProducer::ChoseTracks(std::vector<std::pair<CRTTrack, std::set<unsigned>>> &trackCandidates)
{
  std::vector<std::pair<sbnd::crt::CRTTrack, std::set<unsigned>>> chosenTracks;

  std::set<unsigned> used;

  for(auto const& [track, spIDs] : trackCandidates)
    {
      bool keep = true;
      for(auto const& spID : spIDs)
        {
          if(used.count(spID) != 0)
            keep = false;
        }
      
      if(keep)
        {
          chosenTracks.emplace_back(track, spIDs);

          for(auto const& spID : spIDs)
            used.insert(spID);
        }
    }

  return chosenTracks;
}

double sbnd::crt::CRTTrackProducer::DistanceOfClosestApproach(const CRTTagger tagger, const art::Ptr<CRTSpacePoint> &spacePoint,
                                                              const geo::Point_t &start, const geo::Vector_t &dir)
{
  const CoordSet constrainedPlane = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger);
  const geo::Point_t spacePointCentre = spacePoint->Pos();
  const geo::Point_t spacePointWidths = spacePoint->Err();

  double k;

  switch(constrainedPlane)
    {
    case kX:
      k = (spacePointCentre.X() - start.X()) / dir.X();
      break;
    case kY:
      k = (spacePointCentre.Y() - start.Y()) / dir.Y();
      break;
    case kZ:
      k = (spacePointCentre.Z() - start.Z()) / dir.Z();
      break;
    default:
      std::cout << "Tagger not defined in one plane" << std::endl;
      k = 999999.;
      break;
    }

  const geo::Point_t planePoint = start + k * dir;

  if(IsPointInsideBox(planePoint, spacePointCentre, spacePointWidths))
    return 0.;

  return DistanceOfClosestApproach(constrainedPlane, planePoint, spacePointCentre, spacePointWidths);
}

bool sbnd::crt::CRTTrackProducer::IsPointInsideBox(const geo::Point_t &point, const geo::Point_t &centre, const geo::Point_t &widths)
{
  return (point.X() < centre.X() + widths.X())
    && (point.X() > centre.X() - widths.X())
    && (point.Y() < centre.Y() + widths.Y())
    && (point.Y() > centre.Y() - widths.Y())
    && (point.Z() < centre.Z() + widths.Z())
    && (point.Z() > centre.Z() - widths.Z());
}

double sbnd::crt::CRTTrackProducer::DistanceOfClosestApproach(const CoordSet &constrainedPlane, const geo::Point_t &point,
                                                              const geo::Point_t &centre, const geo::Point_t &widths)
{
  switch(constrainedPlane)
    {
    case kX:
      return MinimumApproach(centre.Y(), widths.Y(), centre.Z(), widths.Z(), {point.Y(), point.Z()});
    case kY:
      return MinimumApproach(centre.X(), widths.X(), centre.Z(), widths.Z(), {point.X(), point.Z()});
    case kZ:
      return MinimumApproach(centre.X(), widths.X(), centre.Y(), widths.Y(), {point.X(), point.Y()});
    default:
      std::cout << "Tagger not defined in one plane" << std::endl;
      return 999999.;
    }
}

double sbnd::crt::CRTTrackProducer::MinimumApproach(const double &x, const double &dx, const double &y, const double &dy, const geo::Point2D_t &p)
{
  const geo::Point2D_t v1(x - dx, y - dy);
  const geo::Point2D_t v2(x - dx, y + dy);
  const geo::Point2D_t v3(x + dx, y - dy);
  const geo::Point2D_t v4(x + dx, y + dy);

  return std::min({DistanceOfClosestApproach(v1, v2, p),
                   DistanceOfClosestApproach(v4, v2, p),
                   DistanceOfClosestApproach(v3, v4, p),
                   DistanceOfClosestApproach(v1, v3, p)
        });
}

double sbnd::crt::CRTTrackProducer::DistanceOfClosestApproach(const geo::Point2D_t &v1, const geo::Point2D_t &v2, const geo::Point2D_t &p)
{
  const geo::Vector2D_t line = v2 - v1;
  const geo::Vector2D_t v1p  = p - v1;
  const geo::Vector2D_t v2p  = p - v2;

  const double angle1 = TMath::RadToDeg() * TMath::ACos(line.Dot(v1p) / line.R() * v1p.R());
  const double angle2 = TMath::RadToDeg() * TMath::ACos(line.Dot(v2p) / line.R() * v2p.R());

  if((angle1 > 90 && angle2 > 90) || (angle1 < 90 && angle2 < 90))
    return std::min(v1p.R(), v2p.R());
  else
    return std::abs((line.Y() * (p.X() - v1.X())) + (line.X() * (p.Y() - v1.Y()))) / line.R();
}

void sbnd::crt::CRTTrackProducer::BestFitLine(const geo::Point_t &a, const geo::Point_t &b, const geo::Point_t &c, const CRTTagger &primary_tagger, 
                                              const CRTTagger &secondary_tagger, const CRTTagger &tertiary_tagger, geo::Point_t &start, 
                                              geo::Point_t &mid, geo::Point_t &end, double &gof)
{
  Eigen::Matrix3d X {
    {a.X(), a.Y(), a.Z()},
    {b.X(), b.Y(), b.Z()},
    {c.X(), c.Y(), c.Z()}
  };

  Eigen::Matrix3d P {
    {2/3., -1/3., -1/3.},
    {-1/3., 2/3., -1/3.},
    {-1/3., -1/3., 2/3.}
  };

  Eigen::Matrix3d sol = X.transpose() * P * X;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenMat(sol);

  if(eigenMat.info() != Eigen::ComputationInfo::Success)
    {
      std::cout << "Decomposition Issue" << std::endl;
      return;
    }

  const auto &eigenValues(eigenMat.eigenvalues());

  std::vector<std::pair<float, unsigned>> eigenValuesAndColumns = {
    {eigenValues[0], 0},
    {eigenValues[1], 1},
    {eigenValues[2], 2}
  };

  std::sort(eigenValuesAndColumns.begin(), eigenValuesAndColumns.end(),
            [](const auto &a, const auto& b) { return a.first > b.first; });

  const Eigen::Matrix3d eigenVectors = eigenMat.eigenvectors();
  const unsigned column = eigenValuesAndColumns[0].second;

  const Eigen::Vector3d eigenVector = { eigenVectors(0, column), eigenVectors(1, column), eigenVectors(2, column) };

  geo::Point_t mean = {a.X() + b.X() + c.X(), a.Y() + b.Y() + c.Y(), a.Z() + b.Z() + c.Z()};
  mean /= 3.;

  const geo::Vector_t dir = {eigenVector(0), eigenVector(1), eigenVector(2)};

  start = LineTaggerIntersectionPoint(mean, dir, primary_tagger);
  mid   = LineTaggerIntersectionPoint(mean, dir, secondary_tagger);
  end   = LineTaggerIntersectionPoint(mean, dir, tertiary_tagger);
  gof   = eigenValuesAndColumns[0].first;

}

geo::Point_t sbnd::crt::CRTTrackProducer::LineTaggerIntersectionPoint(const geo::Point_t &start, const geo::Vector_t &dir, const CRTTagger &tagger)
{
  const CoordSet constrainedPlane = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger);
  const CRTTaggerGeo taggerGeo    = fCRTGeoAlg.GetTagger(CRTCommonUtils::GetTaggerName(tagger));
  double k;

  switch(constrainedPlane)
    {
    case kX:
      {
        const double x = (taggerGeo.maxX + taggerGeo.minX) / 2.;
        k =  (x - start.X()) / dir.X();
      }
      break;
    case kY:
      {
        const double y = (taggerGeo.maxY + taggerGeo.minY) / 2.;
        k = (y - start.Y()) / dir.Y();
      }
      break;
    case kZ:
      {
        const double z = (taggerGeo.maxZ + taggerGeo.minZ) / 2.;
        k = (z - start.Z()) / dir.Z();
      }
      break;
    default:
      std::cout << "Tagger not defined in one plane" << std::endl;
      k = 999999.;
      break;
    }

  return start + k * dir;
}

double sbnd::crt::CRTTrackProducer::TripleTrackToF(std::vector<double> times)
{
  if(times.size() != 3)
    {
      std::cout << "There are not 3 space points in this track" << std::endl;
      return -std::numeric_limits<double>::max();
    }
  
  std::sort(times.begin(), times.end());

  return times[2] - times[0];
}

DEFINE_ART_MODULE(sbnd::crt::CRTTrackProducer)
