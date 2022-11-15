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
  std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<uint16_t>>> CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits);
  sbnd::crt::CRTCluster CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits);

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
    produces<art::Assns<sbnd::crt::CRTStripHit, sbnd::crt::CRTCluster>>();
  }

void sbnd::CRTClusterProducer::produce(art::Event& e)
{
  auto clusterVec          = std::make_unique<std::vector<sbnd::crt::CRTCluster>>();
  auto clusterStripHitAssn = std::make_unique<art::Assns<sbnd::crt::CRTStripHit, sbnd::crt::CRTCluster>>();
  
  art::Handle<std::vector<sbnd::crt::CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  
  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap = GroupStripHits(CRTStripHitVec);

  for(auto& [tagger, stripHits] : taggerStripHitsMap)
    {
      std::sort(stripHits.begin(), stripHits.end(), [](art::Ptr<sbnd::crt::CRTStripHit> &a, art::Ptr<sbnd::crt::CRTStripHit> &b)->bool{
	  return a->Ts1() < b->Ts1();});

      std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<uint16_t>>> ClustersAndHitIDs = CreateClusters(stripHits);

    }

  e.put(std::move(clusterVec));
  e.put(std::move(clusterStripHitAssn));
}

std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> sbnd::CRTClusterProducer::GroupStripHits(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec)
{
  std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTStripHit>>> taggerStripHitsMap;

  for(const art::Ptr<sbnd::crt::CRTStripHit> &stripHit : CRTStripHitVec)
    {
      CRTStripGeo strip   = fCRTGeoAlg.GetStrip(stripHit->Channel());
      CRTTaggerGeo tagger = fCRTGeoAlg.GetTagger(strip.name);

      if(taggerStripHitsMap.find(tagger.name) != taggerStripHitsMap.end())
	taggerStripHitsMap[tagger.name].push_back(stripHit);
      else
	{
	  taggerStripHitsMap[tagger.name] = std::vector<art::Ptr<sbnd::crt::CRTStripHit>>();
	  taggerStripHitsMap[tagger.name].push_back(stripHit);
	}
    }

  return taggerStripHitsMap;
}

std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<uint16_t>>> sbnd::CRTClusterProducer::CreateClusters(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &stripHits)
{
  std::vector<std::pair<sbnd::crt::CRTCluster, std::vector<uint16_t>>> clustersAndHitIDs;

  std::vector<bool> used(stripHits.size(), false);

  for(uint16_t i = 0; i < stripHits.size(); ++i)
    {
      const art::Ptr<sbnd::crt::CRTStripHit> &initialStripHit = stripHits[i];

      if(!used[i])
	{
	  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> clusteredHits;
	  std::vector<uint16_t> hitIDs;
	  clusteredHits.push_back(initialStripHit);
	  hitIDs.push_back(i);
	  used[i] = true;

	  for(uint16_t ii = i+1; ii < stripHits.size(); ++ii)
	    {
	      const art::Ptr<sbnd::crt::CRTStripHit> &stripHit = stripHits[i];
	      
	      if(!used[ii])
		{
		  if(stripHit->Ts1() - initialStripHit->Ts1() < fCoincidenceTimeRequirement)
		    {
		      clusteredHits.push_back(stripHit);
		      hitIDs.push_back(ii);
		      used[ii] = true;
		    }
		}
	      const sbnd::crt::CRTCluster &cluster = CharacteriseCluster(clusteredHits);
	      clustersAndHitIDs.push_back(std::pair<sbnd::crt::CRTCluster, std::vector<uint16_t>>(cluster, hitIDs));
	    }
	}      
    }
  return clustersAndHitIDs;
}

sbnd::crt::CRTCluster sbnd::CRTClusterProducer::CharacteriseCluster(const std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &clusteredHits)
{
  const uint16_t nHits = clusteredHits.size();
  uint16_t s = 0, ts0 = 0, ts1 = 0;
  uint32_t adc = 0, adcCorr = 0;

  CRTStripGeo strip   = fCRTGeoAlg.GetStrip(clusteredHits.at(0)->Channel());
  CRTTaggerGeo tagger = fCRTGeoAlg.GetTagger(strip.name);
  
  for(auto const& hit : clusteredHits)
    {
      s       += hit->UnixS();
      ts0     += hit->Ts0();
      ts1     += hit->Ts1();
      adc     += hit->ADC1();
      adc     += hit->ADC2();
      adcCorr += hit->ADC1();
      adcCorr += hit->ADC2();
    }

  s   /= nHits;
  ts0 /= nHits;
  ts1 /= nHits;

  return sbnd::crt::CRTCluster(s, ts0, ts1, 0., 0., 0., 0., 0., 0., nHits, adc, adcCorr, tagger.name);
}

DEFINE_ART_MODULE(sbnd::CRTClusterProducer)
