////////////////////////////////////////////////////////////////////////
// Class:       CRTClusterProducer
// Plugin Type: producer
// File:        CRTClusterProducer_module.cc
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include <memory>

namespace sbnd::crt {
  class CRTClusterProducer;
}


class sbnd::crt::CRTClusterProducer : public art::EDProducer {
public:
  explicit CRTClusterProducer(fhicl::ParameterSet const& p);

  CRTClusterProducer(CRTClusterProducer const&) = delete;
  CRTClusterProducer(CRTClusterProducer&&) = delete;
  CRTClusterProducer& operator=(CRTClusterProducer const&) = delete;
  CRTClusterProducer& operator=(CRTClusterProducer&&) = delete;

  void produce(art::Event& e) override;

  std::map<sbnd::crt::CRTTagger, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> GroupStripHits(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec);

  std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits);

  sbnd::crt::CRTCluster CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits);

private:

  sbnd::crt::CRTGeoAlg fCRTGeoAlg;
  std::string          fCRTStripHitModuleLabel;
  uint32_t             fCoincidenceTimeRequirement;
};


sbnd::crt::CRTClusterProducer::CRTClusterProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fCRTStripHitModuleLabel(p.get<std::string>("CRTStripHitModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<uint32_t>("CoincidenceTimeRequirement"))
  {
    produces<std::vector<sbnd::crt::CRTCluster>>();
    produces<art::Assns<sbnd::crt::CRTCluster, sbnd::crt::CRTStripHit>>();
  }

void sbnd::crt::CRTClusterProducer::produce(art::Event& e)
{
  auto clusterVec          = std::make_unique<std::vector<sbnd::crt::CRTCluster>>();
  auto clusterStripHitAssn = std::make_unique<art::Assns<sbnd::crt::CRTCluster, sbnd::crt::CRTStripHit>>();
  
  art::Handle<std::vector<sbnd::crt::CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  
  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  std::map<sbnd::crt::CRTTagger, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap = GroupStripHits(CRTStripHitVec);

  for(auto& [tagger, stripHits] : taggerStripHitsMap)
    {
      std::sort(stripHits.begin(), stripHits.end(), [](art::Ptr<sbnd::crt::CRTStripHit> &a, art::Ptr<sbnd::crt::CRTStripHit> &b)->bool{
          return a->Ts1() < b->Ts1();});

      std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> clustersAndHits = CreateClusters(stripHits);

      for(auto const& [cluster, clusteredHits] : clustersAndHits)
        {
          clusterVec->push_back(cluster);
          util::CreateAssn(*this, e, *clusterVec, clusteredHits, *clusterStripHitAssn);
        }
    }

  e.put(std::move(clusterVec));
  e.put(std::move(clusterStripHitAssn));
}

std::map<sbnd::crt::CRTTagger, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> sbnd::crt::CRTClusterProducer::GroupStripHits(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec)
{
  std::map<sbnd::crt::CRTTagger, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap;

  for(const art::Ptr<sbnd::crt::CRTStripHit> &stripHit : CRTStripHitVec)
    {
      sbnd::crt::CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());

      if(taggerStripHitsMap.find(tagger) != taggerStripHitsMap.end())
        taggerStripHitsMap[tagger].push_back(stripHit);
      else
        {
          taggerStripHitsMap[tagger] = std::vector<art::Ptr<sbnd::crt::CRTStripHit>>();
          taggerStripHitsMap[tagger].push_back(stripHit);
        }
    }

  return taggerStripHitsMap;
}

std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> sbnd::crt::CRTClusterProducer::CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits)
{
  std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> clustersAndHits;

  std::vector<bool> used(stripHits.size(), false);

  for(uint16_t i = 0; i < stripHits.size(); ++i)
    {
      const art::Ptr<sbnd::crt::CRTStripHit> &initialStripHit = stripHits[i];

      if(!used[i])
        {
          std::vector<art::Ptr<sbnd::crt::CRTStripHit>> clusteredHits;
          clusteredHits.push_back(initialStripHit);
          used[i] = true;

          for(uint16_t ii = i+1; ii < stripHits.size(); ++ii)
            {
              const art::Ptr<sbnd::crt::CRTStripHit> &stripHit = stripHits[ii];
    
              if(!used[ii])
                {
                  if(stripHit->Ts1() - initialStripHit->Ts1() < fCoincidenceTimeRequirement)
                    {
                      clusteredHits.push_back(stripHit);
                      used[ii] = true;
                    }
                }
            }
          const sbnd::crt::CRTCluster &cluster = CharacteriseCluster(clusteredHits);
          clustersAndHits.emplace_back(cluster, clusteredHits);
        }      
    }
  return clustersAndHits;
}

sbnd::crt::CRTCluster sbnd::crt::CRTClusterProducer::CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits)
{
  const uint16_t nHits = clusteredHits.size();

  const CRTStripGeo strip0 = fCRTGeoAlg.GetStrip(clusteredHits.at(0)->Channel());
  const sbnd::crt::CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(clusteredHits.at(0)->Channel());

  uint32_t ts0 = 0, ts1 = 0, s = 0;
  bool threeD = false;

  for(uint16_t i = 0; i < clusteredHits.size(); ++i)
    {
      const art::Ptr<sbnd::crt::CRTStripHit> hit = clusteredHits[i];

      ts0 += hit->Ts0();
      ts1 += hit->Ts1();
      s   += hit->UnixS();

      const CRTStripGeo strip = fCRTGeoAlg.GetStrip(hit->Channel());
      threeD |= fCRTGeoAlg.DifferentOrientations(strip0, strip);
    }

  s   /= nHits;
  ts0 /= nHits;
  ts1 /= nHits;

  return sbnd::crt::CRTCluster(ts0, ts1, s, nHits, tagger, threeD);
}
  
DEFINE_ART_MODULE(sbnd::crt::CRTClusterProducer)
