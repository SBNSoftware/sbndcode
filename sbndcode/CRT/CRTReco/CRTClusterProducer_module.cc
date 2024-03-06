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

  std::map<CRTTagger, std::vector<art::Ptr<CRTStripHit>>> GroupStripHits(const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec);

  std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> CreateClusters(const std::vector<art::Ptr<CRTStripHit>> &stripHits);

  std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> SplitClusters(const std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> &initialClusters);

  CRTCluster CharacteriseCluster(const std::vector<art::Ptr<CRTStripHit>> &clusteredHits);

private:

  CRTGeoAlg   fCRTGeoAlg;
  std::string fCRTStripHitModuleLabel;
  uint32_t    fCoincidenceTimeRequirement;
  double      fOverlapBuffer;
};


sbnd::crt::CRTClusterProducer::CRTClusterProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fCRTStripHitModuleLabel(p.get<std::string>("CRTStripHitModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<uint32_t>("CoincidenceTimeRequirement"))
  , fOverlapBuffer(p.get<double>("OverlapBuffer"))
  {
    produces<std::vector<CRTCluster>>();
    produces<art::Assns<CRTCluster, CRTStripHit>>();
  }

void sbnd::crt::CRTClusterProducer::produce(art::Event& e)
{
  auto clusterVec          = std::make_unique<std::vector<CRTCluster>>();
  auto clusterStripHitAssn = std::make_unique<art::Assns<CRTCluster, CRTStripHit>>();
  
  art::Handle<std::vector<CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  
  std::vector<art::Ptr<CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  std::map<CRTTagger, std::vector<art::Ptr<CRTStripHit>>> taggerStripHitsMap = GroupStripHits(CRTStripHitVec);

  for(auto& [tagger, stripHits] : taggerStripHitsMap)
    {
      std::sort(stripHits.begin(), stripHits.end(), [](art::Ptr<CRTStripHit> &a, art::Ptr<CRTStripHit> &b)->bool{
          return a->Ts1() < b->Ts1();});

      std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> clustersAndHits = CreateClusters(stripHits);

      for(auto const& [cluster, clusteredHits] : clustersAndHits)
        {
          clusterVec->push_back(cluster);
          util::CreateAssn(*this, e, *clusterVec, clusteredHits, *clusterStripHitAssn);
        }
    }

  e.put(std::move(clusterVec));
  e.put(std::move(clusterStripHitAssn));
}

std::map<sbnd::crt::CRTTagger, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> sbnd::crt::CRTClusterProducer::GroupStripHits(const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec)
{
  std::map<CRTTagger, std::vector<art::Ptr<CRTStripHit>>> taggerStripHitsMap;

  for(const art::Ptr<CRTStripHit> &stripHit : CRTStripHitVec)
    {
      CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());

      if(taggerStripHitsMap.find(tagger) != taggerStripHitsMap.end())
        taggerStripHitsMap[tagger].push_back(stripHit);
      else
        {
          taggerStripHitsMap[tagger] = std::vector<art::Ptr<CRTStripHit>>();
          taggerStripHitsMap[tagger].push_back(stripHit);
        }
    }

  return taggerStripHitsMap;
}

std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> sbnd::crt::CRTClusterProducer::CreateClusters(const std::vector<art::Ptr<CRTStripHit>> &stripHits)
{
  std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> clustersAndHits;

  std::vector<bool> used(stripHits.size(), false);

  for(uint16_t i = 0; i < stripHits.size(); ++i)
    {
      const art::Ptr<CRTStripHit> &initialStripHit = stripHits[i];

      if(!used[i])
        {
          std::vector<art::Ptr<CRTStripHit>> clusteredHits;
          clusteredHits.push_back(initialStripHit);
          used[i] = true;

          for(uint16_t ii = i+1; ii < stripHits.size(); ++ii)
            {
              const art::Ptr<CRTStripHit> &stripHit = stripHits[ii];
    
              if(!used[ii])
                {
                  if(stripHit->Ts1() - initialStripHit->Ts1() < fCoincidenceTimeRequirement)
                    {
                      clusteredHits.push_back(stripHit);
                      used[ii] = true;
                    }
                }
            }
          
          const CRTCluster &cluster = CharacteriseCluster(clusteredHits);
          clustersAndHits.emplace_back(cluster, clusteredHits);
        }
    }
  return SplitClusters(clustersAndHits);
}

std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>>
  sbnd::crt::CRTClusterProducer::SplitClusters(const std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> &initialClusters)
{
  std::vector<std::pair<CRTCluster, std::vector<art::Ptr<CRTStripHit>>>> clustersAndHits;

  for(auto const& [cluster, hits] : initialClusters)
    {
      std::map<uint16_t, std::set<uint16_t>> overlaps;
      
      std::vector<bool> used(hits.size(), false);

      for(uint16_t j = 0; j < hits.size(); ++j)
        {
          const art::Ptr<CRTStripHit> &hit1 = hits[j];
          overlaps[j] = {j};
          
          for(uint16_t jj = 0; jj < hits.size(); ++jj)
            {
              const art::Ptr<CRTStripHit> &hit2 = hits[jj];
              if(fCRTGeoAlg.CheckOverlap(hit1->Channel(), hit2->Channel(), 10.))
                overlaps[j].insert(jj);
            }
        }

      std::set<uint16_t> leftovers;
      
      for(auto const& [id, overlap_set] : overlaps)
        {
          if(used[id])
            continue;

          bool exclusive = true;
          
          for(auto const& id2 : overlap_set)
            {
              if(overlaps[id2] != overlap_set)
                exclusive = false;
            }
          
          if(exclusive)
            {
              std::vector<art::Ptr<CRTStripHit>> newClusteredHits;
              
              for(auto const& id2 : overlap_set)
                {
                  newClusteredHits.push_back(hits[id2]);
                  used[id2] = true;
                }
              
              const CRTCluster &cluster = CharacteriseCluster(newClusteredHits);
              clustersAndHits.emplace_back(cluster, newClusteredHits);
            }
          else
            leftovers.insert(id);
        }
      
      std::vector<art::Ptr<CRTStripHit>> leftoverClusteredHits;
      
      for(auto const& id : leftovers)
        leftoverClusteredHits.push_back(hits[id]);
      
      if(leftoverClusteredHits.size() != 0)
        {
          const CRTCluster &cluster = CharacteriseCluster(leftoverClusteredHits);
          clustersAndHits.emplace_back(cluster, leftoverClusteredHits);
        }
    }

  return clustersAndHits;
}     

sbnd::crt::CRTCluster sbnd::crt::CRTClusterProducer::CharacteriseCluster(const std::vector<art::Ptr<CRTStripHit>> &clusteredHits)
{
  const uint16_t nHits = clusteredHits.size();

  const CRTStripGeo strip0 = fCRTGeoAlg.GetStrip(clusteredHits.at(0)->Channel());
  const CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(clusteredHits.at(0)->Channel());

  uint32_t ts0 = 0, ts1 = 0, s = 0;
  CoordSet composition = kUndefinedSet;

  for(uint16_t i = 0; i < clusteredHits.size(); ++i)
    {
      const art::Ptr<CRTStripHit> hit = clusteredHits[i];

      ts0 += hit->Ts0();
      ts1 += hit->Ts1();
      s   += hit->UnixS();

      const CRTStripGeo strip = fCRTGeoAlg.GetStrip(hit->Channel());
      if(fCRTGeoAlg.DifferentOrientations(strip0, strip))
        composition = kXYZ;
    }

  if(composition == kUndefinedSet)
    composition = fCRTGeoAlg.GlobalConstrainedCoordinates(strip0.channel0);

  s   /= nHits;
  ts0 /= nHits;
  ts1 /= nHits;

  return CRTCluster(ts0, ts1, s, nHits, tagger, composition);
}
  
DEFINE_ART_MODULE(sbnd::crt::CRTClusterProducer)
