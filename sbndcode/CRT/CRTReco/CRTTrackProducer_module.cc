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

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include <memory>

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

private:

  std::string fCRTSpacePointModuleLabel;
  double      fCoincidenceTimeRequirement;
};


sbnd::crt::CRTTrackProducer::CRTTrackProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<double>("CoincidenceTimeRequirement"))
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

          if(CRTCommonUtils::IsHorizontalTagger(primaryCluster->Tagger()) && CRTCommonUtils::IsHorizontalTagger(secondaryCluster->Tagger()))
            {
              for(unsigned iii = ii + 1; iii < spacePointVec.size(); ++iii)
                {
                  const art::Ptr<CRTSpacePoint> tertiarySpacePoint = spacePointVec[ii];
                  const art::Ptr<CRTCluster> tertiaryCluster       = spacePointsToCluster.at(tertiarySpacePoint.key());

                  if(!CRTCommonUtils::IsHorizontalTagger(tertiaryCluster->Tagger()))
                    continue;

                  if(tertiarySpacePoint->Time() - primarySpacePoint->Time() > fCoincidenceTimeRequirement ||
                     tertiaryCluster->Tagger() == primaryCluster->Tagger())
                    continue;

                  // const geo::Point_t &tertiaryPos = tertiarySpacePoint->Pos();

                }
            }
          else
            {
              double time, etime;
              TimeErrorCalculator({primarySpacePoint->Time(), secondarySpacePoint->Time()}, time, etime);
              const double pe = primarySpacePoint->PE() + secondarySpacePoint->PE();
              const std::vector<CRTTagger> used_taggers = {primaryCluster->Tagger(), secondaryCluster->Tagger()};

              const CRTTrack track(start, dir, time, etime, pe, false, used_taggers);
              const std::set<unsigned> used_spacepoints = {i, ii};

              candidates.emplace_back(track, used_spacepoints);
            }
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

DEFINE_ART_MODULE(sbnd::crt::CRTTrackProducer)
