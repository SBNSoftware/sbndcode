////////////////////////////////////////////////////////////////////////
// Class:       CRTBlobProducer
// Plugin Type: producer
// File:        CRTBlobProducer_module.cc
//
// Author:      Henry Lay (h.lay@sheffield.ac.uk)
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
#include "sbnobj/SBND/CRT/CRTBlob.hh"

#include <numeric>

namespace sbnd::crt {
  class CRTBlobProducer;
}

class sbnd::crt::CRTBlobProducer : public art::EDProducer {
public:
  explicit CRTBlobProducer(fhicl::ParameterSet const& p);

  CRTBlobProducer(CRTBlobProducer const&) = delete;
  CRTBlobProducer(CRTBlobProducer&&) = delete;
  CRTBlobProducer& operator=(CRTBlobProducer const&) = delete;
  CRTBlobProducer& operator=(CRTBlobProducer&&) = delete;

  void produce(art::Event& e) override;

  void OrderSpacePoints(std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec, const art::FindOneP<CRTCluster> &spacePointsToCluster);

  std::vector<std::pair<CRTBlob, std::set<unsigned>>> CreateBlobCandidates(const std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec,
                                                                           const art::FindOneP<CRTCluster> &spacePointsToCluster);

  void TimeErrorCalculator(const std::vector<double> &times, double &mean, double &err);

  void OrderBlobCandidates(std::vector<std::pair<CRTBlob, std::set<unsigned>>> &blobCandidates);

  std::vector<std::pair<sbnd::crt::CRTBlob, std::set<unsigned>>> ChoseBlobs(std::vector<std::pair<CRTBlob, std::set<unsigned>>> &blobCandidates);

private:

  std::string      fCRTSpacePointModuleLabel;
  double           fCoincidenceTimeRequirement;
  bool             fUseTs0;
};


sbnd::crt::CRTBlobProducer::CRTBlobProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<double>("CoincidenceTimeRequirement"))
  , fUseTs0(p.get<bool>("UseTs0"))
{
  produces<std::vector<CRTBlob>>();
  produces<art::Assns<CRTSpacePoint, CRTBlob>>();
}

void sbnd::crt::CRTBlobProducer::produce(art::Event& e)
{
  auto blobVec            = std::make_unique<std::vector<CRTBlob>>();
  auto blobSpacePointAssn = std::make_unique<art::Assns<CRTSpacePoint, CRTBlob>>();

  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::FindOneP<CRTCluster> spacePointsToCluster(CRTSpacePointHandle, e, fCRTSpacePointModuleLabel);

  OrderSpacePoints(CRTSpacePointVec, spacePointsToCluster);

  std::vector<std::pair<CRTBlob, std::set<unsigned>>> blobCandidates = CreateBlobCandidates(CRTSpacePointVec, spacePointsToCluster);

  OrderBlobCandidates(blobCandidates);

  std::vector<std::pair<CRTBlob, std::set<unsigned>>> chosenBlobs = ChoseBlobs(blobCandidates);

  for(auto const& [blob, spIDs] : chosenBlobs)
    {
      blobVec->push_back(blob);

      for(auto const& spID : spIDs)
        util::CreateAssn(*this, e, *blobVec, CRTSpacePointVec[spID], *blobSpacePointAssn);
    }

  e.put(std::move(blobVec));
  e.put(std::move(blobSpacePointAssn));
}

void sbnd::crt::CRTBlobProducer::OrderSpacePoints(std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec, const art::FindOneP<CRTCluster> &spacePointsToCluster)
{
  std::sort(spacePointVec.begin(), spacePointVec.end(), 
            [&](const art::Ptr<CRTSpacePoint> &a, const art::Ptr<CRTSpacePoint> &b) -> bool {
              if(fUseTs0)
                return a->Ts0() < b->Ts0();
              else
                return a->Ts1() < b->Ts1();
            });
}

std::vector<std::pair<sbnd::crt::CRTBlob, std::set<unsigned>>> sbnd::crt::CRTBlobProducer::CreateBlobCandidates(const std::vector<art::Ptr<CRTSpacePoint>> &spacePointVec,
                                                                                                                const art::FindOneP<CRTCluster> &spacePointsToCluster)
{
  std::vector<std::pair<CRTBlob, std::set<unsigned>>> candidates;

  for(unsigned i = 0; i < spacePointVec.size(); ++i)
    {
      const art::Ptr<CRTSpacePoint> primarySpacePoint = spacePointVec[i];
      const art::Ptr<CRTCluster> primaryCluster       = spacePointsToCluster.at(primarySpacePoint.key());

      std::set<unsigned> used_spacepoints = { i };
      std::vector<double> t0s             = { primarySpacePoint->Ts0() };
      std::vector<double> t1s             = { primarySpacePoint->Ts1() };
      std::vector<double> pes             = { primarySpacePoint->PE() };

      std::map<CRTTagger, int> taggers;
      for(int i = 0; i < 7; ++i)
        taggers[(CRTTagger)i] = 0;
      ++taggers[primaryCluster->Tagger()];

      for(unsigned ii = i+1; ii < spacePointVec.size(); ++ii)
        {
          const art::Ptr<CRTSpacePoint> secondarySpacePoint = spacePointVec[ii];
          const art::Ptr<CRTCluster> secondaryCluster       = spacePointsToCluster.at(secondarySpacePoint.key());

          const double tdiff_prim_sec = fUseTs0 ? secondarySpacePoint->Ts0() - primarySpacePoint->Ts0() : secondarySpacePoint->Ts1() - primarySpacePoint->Ts1();

          if(tdiff_prim_sec > fCoincidenceTimeRequirement)
            break;

          used_spacepoints.insert(ii);
          t0s.push_back(secondarySpacePoint->Ts0());
          t1s.push_back(secondarySpacePoint->Ts1());
          pes.push_back(secondarySpacePoint->PE());
          ++taggers[secondaryCluster->Tagger()];
        }
      
      double t0, et0;
      TimeErrorCalculator(t0s, t0, et0);

      double t1, et1;
      TimeErrorCalculator(t1s, t1, et1);

      const double pe = std::accumulate(pes.begin(), pes.end(), 0.);

      const CRTBlob blob(t0, et0, t1, et1, pe, taggers);

      candidates.emplace_back(blob, used_spacepoints);
    }
  
  return candidates;
}

void sbnd::crt::CRTBlobProducer::TimeErrorCalculator(const std::vector<double> &times, double &mean, double &err)
{
  double sum = 0.;
  for(auto const &time : times)
    sum += time;

  mean = sum / times.size();

  double summed_var = 0.;
  for(auto const &time : times)
    summed_var += std::pow((time - mean), 2);

  err = std::sqrt(summed_var / (times.size() - 1));
}

void sbnd::crt::CRTBlobProducer::OrderBlobCandidates(std::vector<std::pair<CRTBlob, std::set<unsigned>>> &blobCandidates)
{
  std::sort(blobCandidates.begin(), blobCandidates.end(),
            [&](const std::pair<CRTBlob, std::set<unsigned>> &a, const std::pair<CRTBlob, std::set<unsigned>> &b) -> bool {
              return fUseTs0 ? a.first.Ts0Err() < b.first.Ts0Err() : a.first.Ts1Err() < b.first.Ts1Err();
            });
}

std::vector<std::pair<sbnd::crt::CRTBlob, std::set<unsigned>>> sbnd::crt::CRTBlobProducer::ChoseBlobs(std::vector<std::pair<CRTBlob, std::set<unsigned>>> &blobCandidates)
{
  std::vector<std::pair<sbnd::crt::CRTBlob, std::set<unsigned>>> chosenBlobs;

  std::set<unsigned> used;

  for(auto const& [blob, spIDs] : blobCandidates)
    {
      bool keep = true;
      for(auto const& spID : spIDs)
        {
          if(used.count(spID) != 0)
            keep = false;
        }
      
      if(keep)
        {
          chosenBlobs.emplace_back(blob, spIDs);

          for(auto const& spID : spIDs)
            used.insert(spID);
        }
    }

  return chosenBlobs;
}

DEFINE_ART_MODULE(sbnd::crt::CRTBlobProducer)
