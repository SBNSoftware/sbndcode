////////////////////////////////////////////////////////////////////////
// Class:       CRTSpacePointProducer
// Plugin Type: producer
// File:        CRTSpacePointProducer_module.cc
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

#include "sbndcode/CRT/CRTReco/CRTClusterCharacterisationAlg.h"

namespace sbnd::crt {
  class CRTSpacePointProducer;
}

class sbnd::crt::CRTSpacePointProducer : public art::EDProducer {
public:
  explicit CRTSpacePointProducer(fhicl::ParameterSet const& p);

  CRTSpacePointProducer(CRTSpacePointProducer const&) = delete;
  CRTSpacePointProducer(CRTSpacePointProducer&&) = delete;
  CRTSpacePointProducer& operator=(CRTSpacePointProducer const&) = delete;
  CRTSpacePointProducer& operator=(CRTSpacePointProducer&&) = delete;

  void produce(art::Event& e) override;

private:

  CRTClusterCharacterisationAlg fClusterCharacAlg;
  std::string fClusterModuleLabel;
};


sbnd::crt::CRTSpacePointProducer::CRTSpacePointProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fClusterCharacAlg(p.get<fhicl::ParameterSet>("ClusterCharacterisationAlg", fhicl::ParameterSet()))
  , fClusterModuleLabel(p.get<std::string>("ClusterModuleLabel"))
  {
    produces<std::vector<CRTSpacePoint>>();
    produces<art::Assns<CRTCluster, CRTSpacePoint>>();
  }

void sbnd::crt::CRTSpacePointProducer::produce(art::Event& e)
{
  auto spacePointVec         = std::make_unique<std::vector<CRTSpacePoint>>();
  auto spacePointClusterAssn = std::make_unique<art::Assns<CRTCluster, CRTSpacePoint>>();

  art::Handle<std::vector<CRTCluster>> clusterHandle;
  e.getByLabel(fClusterModuleLabel, clusterHandle);

  std::vector<art::Ptr<CRTCluster>> clusterVec;
  art::fill_ptr_vector(clusterVec, clusterHandle);

  art::FindManyP<CRTStripHit> clusterToStripHits(clusterHandle, e, fClusterModuleLabel);

  for(const art::Ptr<CRTCluster> &cluster : clusterVec)
    {
      const uint nhits = cluster->NHits();
      const std::vector<art::Ptr<CRTStripHit>> stripHits = clusterToStripHits.at(cluster.key());

      if(nhits == 1 && cluster->Tagger() == kBottomTagger)
        {
          CRTSpacePoint spacepoint = fClusterCharacAlg.CharacteriseSingleHitCluster(cluster, stripHits[0]);
            {
              spacePointVec->push_back(spacepoint);
              util::CreateAssn(*this, e, *spacePointVec, cluster, *spacePointClusterAssn);
            }
        }
      else if(nhits > 1)
        {
          CRTSpacePoint spacepoint;
          if(fClusterCharacAlg.CharacteriseMultiHitCluster(cluster, stripHits, spacepoint))
            {
              spacePointVec->push_back(spacepoint);
              util::CreateAssn(*this, e, *spacePointVec, cluster, *spacePointClusterAssn);
            }
        }
    }

  e.put(std::move(spacePointVec));
  e.put(std::move(spacePointClusterAssn));
}

DEFINE_ART_MODULE(sbnd::crt::CRTSpacePointProducer)
