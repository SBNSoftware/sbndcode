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

  CRTHitRecoAlg *hitRecoAlg;

  std::string fFEBDataModuleLabel;

};


sbnd::CRTHitProducer::CRTHitProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  {
    produces<std::vector<sbn::crt::CRTHit> >();
    produces<art::Assns<sbn::crt::CRTHit, sbnd::crt::FEBData> >();
  }

void sbnd::CRTHitProducer::produce(art::Event& e)
{
  std::unique_ptr< std::vector<sbn::crt::CRTHit> > crtHitVec(new std::vector<sbn::crt::CRTHit>);
  std::unique_ptr< art::Assns<sbn::crt::CRTHit, sbnd::crt::FEBData> > crtHitDataAssn(new art::Assns<sbn::crt::CRTHit, sbnd::crt::FEBData>);
  
  art::Handle<std::vector<sbnd::crt::FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  
  std::vector<art::Ptr<sbnd::crt::FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Performs pedestal subtraction and cable length corrections
  // Only keeps strips passing threshold requirement
  std::vector<sbnd::CRTStripHit> stripHits = hitRecoAlg->ProduceStripHits(FEBDataVec);

  e.put(std::move(crtHitVec));
  e.put(std::move(crtHitDataAssn));
}

DEFINE_ART_MODULE(sbnd::CRTHitProducer)
