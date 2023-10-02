////////////////////////////////////////////////////////////////////////
// Class:       CRTHitProducer
// Plugin Type: producer
// File:        CRTHitProducer_module.cc
//
// Generated at Fri Sep 30 10:25:35 2022 by Henry Lay using cetskelgen
// from  version .
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

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/SBND/CRT/FEBData.hh"

#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"

#include <memory>

namespace sbnd {
  class CRTHitProducer;
}


class sbnd::CRTHitProducer : public art::EDProducer {
public:
  explicit CRTHitProducer(fhicl::ParameterSet const& p);

  CRTHitProducer(CRTHitProducer const&) = delete;
  CRTHitProducer(CRTHitProducer&&) = delete;
  CRTHitProducer& operator=(CRTHitProducer const&) = delete;
  CRTHitProducer& operator=(CRTHitProducer&&) = delete;

  void produce(art::Event& e) override;

private:

  CRTHitRecoAlg fHitRecoAlg;
  std::string fFEBDataModuleLabel;

};


sbnd::CRTHitProducer::CRTHitProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fHitRecoAlg(p.get<fhicl::ParameterSet>("HitRecoAlg"))
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  {
    produces<std::vector<sbn::crt::CRTHit>>();
    produces<art::Assns<sbnd::crt::FEBData, sbn::crt::CRTHit>>();

  }

void sbnd::CRTHitProducer::produce(art::Event& e)
{ 
  auto crtHitVec      = std::make_unique<std::vector<sbn::crt::CRTHit>>();
  auto crtHitDataAssn = std::make_unique<art::Assns<sbnd::crt::FEBData, sbn::crt::CRTHit>>();
  
  art::Handle<std::vector<sbnd::crt::FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  
  std::vector<art::Ptr<sbnd::crt::FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Performs pedestal subtraction and cable length corrections
  // Only keeps strips passing threshold requirement
  const std::map<std::string, std::vector<std::vector<sbnd::CRTStripHit>>> stripHits = fHitRecoAlg.ProduceStripHits(FEBDataVec);

  // Finds coincidences between overlapping perpendicular strips
  const std::vector<sbn::crt::CRTHit> crtHits = fHitRecoAlg.ProduceCRTHits(stripHits);
  // Save hits and associations to FEBDatas to the event record
  for(auto const &hit : crtHits)
    {
      crtHitVec->push_back(hit);
      for(auto const &feb_id : hit.feb_id)
        util::CreateAssn(*this, e, *crtHitVec, FEBDataVec[feb_id], *crtHitDataAssn);
    }

  e.put(std::move(crtHitVec));
  e.put(std::move(crtHitDataAssn));
}

DEFINE_ART_MODULE(sbnd::CRTHitProducer)
