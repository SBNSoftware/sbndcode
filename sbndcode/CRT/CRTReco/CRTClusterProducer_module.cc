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

#include <memory>

namespace sbnd {
  class CRTClusterProducer;
}


class sbnd::CRTClusterProducer : public art::EDProducer {
public:
  explicit CRTClusterProducer(fhicl::ParameterSet const& p);

  CRTClusterProducer(CRTClusterProducer const&) = delete;
  CRTClusterProducer(CRTClusterProducer&&) = delete;
  CRTClusterProducer& operator=(CRTClusterProducer const&) = delete;
  CRTClusterProducer& operator=(CRTClusterProducer&&) = delete;

  void produce(art::Event& e) override;

  std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> GroupStripHits(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec);

  std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits);

  sbnd::crt::CRTCluster CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits);

  std::array<double, 6> FindOverlap(const art::Ptr<sbnd::crt::CRTStripHit> &hit0, const art::Ptr<sbnd::crt::CRTStripHit> &hit1, 
                                    const CRTStripGeo &strip0, const CRTStripGeo &strip1);

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;
  std::string     fCRTStripHitModuleLabel;
  uint32_t        fCoincidenceTimeRequirement;
};


sbnd::CRTClusterProducer::CRTClusterProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fCRTStripHitModuleLabel(p.get<std::string>("CRTStripHitModuleLabel"))
  , fCoincidenceTimeRequirement(p.get<uint32_t>("CoincidenceTimeRequirement"))
  {
    produces<std::vector<sbnd::crt::CRTCluster>>();
    produces<art::Assns<sbnd::crt::CRTCluster, sbnd::crt::CRTStripHit>>();
  }

void sbnd::CRTClusterProducer::produce(art::Event& e)
{
  auto clusterVec          = std::make_unique<std::vector<sbnd::crt::CRTCluster>>();
  auto clusterStripHitAssn = std::make_unique<art::Assns<sbnd::crt::CRTCluster, sbnd::crt::CRTStripHit>>();
  
  art::Handle<std::vector<sbnd::crt::CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  
  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap = GroupStripHits(CRTStripHitVec);

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

std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> sbnd::CRTClusterProducer::GroupStripHits(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec)
{
  std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap;

  for(const art::Ptr<sbnd::crt::CRTStripHit> &stripHit : CRTStripHitVec)
    {
      std::string taggerName = fCRTGeoAlg.ChannelToTaggerName(stripHit->Channel());

      if(taggerStripHitsMap.find(taggerName) != taggerStripHitsMap.end())
        taggerStripHitsMap[taggerName].push_back(stripHit);
      else
        {
          taggerStripHitsMap[taggerName] = std::vector<art::Ptr<sbnd::crt::CRTStripHit>>();
          taggerStripHitsMap[taggerName].push_back(stripHit);
        }
    }

  return taggerStripHitsMap;
}

std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>>> sbnd::CRTClusterProducer::CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits)
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

sbnd::crt::CRTCluster sbnd::CRTClusterProducer::CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits)
{
  const uint16_t nHits = clusteredHits.size();
  uint32_t s = 0, ts0 = 0, ts1 = 0;
  uint16_t adc = 0, adcCorr = 0;

  std::string taggerName = fCRTGeoAlg.ChannelToTaggerName(clusteredHits.at(0)->Channel());
  
  std::array<double, 6> edges({-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                               -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                               -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()});

  if(clusteredHits.size() == 1)
    edges = fCRTGeoAlg.StripHit3DPos(fCRTGeoAlg.GetStrip(clusteredHits.at(0)->Channel()).name, clusteredHits.at(0)->Pos(), clusteredHits.at(0)->Error());
  
  for(uint16_t i = 0; i < clusteredHits.size(); ++i)
    {
      const art::Ptr<sbnd::crt::CRTStripHit> hit = clusteredHits[i];
      const CRTStripGeo strip                    = fCRTGeoAlg.GetStrip(hit->Channel());

      s       += hit->UnixS();
      ts0     += hit->Ts0();
      ts1     += hit->Ts1();
      adc     += hit->ADC1();
      adc     += hit->ADC2();
      adcCorr += hit->ADC1();
      adcCorr += hit->ADC2();

      for(uint16_t ii = i+1; ii < clusteredHits.size(); ++ii)
        {
          const art::Ptr<sbnd::crt::CRTStripHit> overlappinghit = clusteredHits[ii];
          const CRTStripGeo overlappingstrip                    = fCRTGeoAlg.GetStrip(overlappinghit->Channel());
          
          if(!fCRTGeoAlg.CheckOverlap(strip, overlappingstrip))
            continue;

          std::array<double, 6> overlap = FindOverlap(hit, overlappinghit, strip, overlappingstrip);

          for(uint16_t iii = 0; iii < 6; iii+=2)
            {
              if(overlap[iii] > edges[iii]) 
                edges[iii] = overlap[iii];
              
              if(overlap[iii+1] < edges[iii+1])
                edges[iii+1] = overlap[iii+1];
            }
        } 
    }

  s   /= nHits;
  ts0 /= nHits;
  ts1 /= nHits;

  return sbnd::crt::CRTCluster(ts0, ts1, s, edges, nHits, adc, adcCorr, taggerName);
}
  
std::array<double, 6> sbnd::CRTClusterProducer::FindOverlap(const art::Ptr<sbnd::crt::CRTStripHit> &hit0, const art::Ptr<sbnd::crt::CRTStripHit> &hit1,
                                                            const CRTStripGeo &strip0, const CRTStripGeo &strip1)
{
  const std::array<double, 6> hit0pos = fCRTGeoAlg.StripHit3DPos(strip0.name, hit0->Pos(), hit0->Error());
  const std::array<double, 6> hit1pos = fCRTGeoAlg.StripHit3DPos(strip1.name, hit1->Pos(), hit1->Error());
  
  std::array<double, 6> overlap({std::max(hit0pos[0], hit1pos[0]),
                                 std::min(hit0pos[1], hit1pos[1]),
                                 std::max(hit0pos[2], hit1pos[2]),
                                 std::min(hit0pos[3], hit1pos[3]),
                                 std::max(hit0pos[4], hit1pos[4]),
                                 std::min(hit0pos[5], hit1pos[5])});

  return overlap;
}

DEFINE_ART_MODULE(sbnd::CRTClusterProducer)
