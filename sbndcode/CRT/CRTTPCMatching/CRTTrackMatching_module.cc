/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTTrackMatching
/// Module Type: producer
/// File:        CRTTrackMatching_module.cc
///
/// Author:         Thomas Brooks
/// E-mail address: tbrooks@fnal.gov
/////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "sbndcode/CRT/CRTTPCMatching/CRTTrackMatchAlg.h"

#include <memory>

namespace sbnd::crt {
  class CRTTrackMatching;
}

class sbnd::crt::CRTTrackMatching : public art::EDProducer {
public:

  explicit CRTTrackMatching(fhicl::ParameterSet const& p);

  CRTTrackMatching(CRTTrackMatching const&) = delete;
  CRTTrackMatching(CRTTrackMatching&&) = delete;
  CRTTrackMatching& operator=(CRTTrackMatching const&) = delete;
  CRTTrackMatching& operator=(CRTTrackMatching&&) = delete;

  void produce(art::Event& e) override;

private:

    CRTTrackMatchAlg fMatchingAlg;
    art::InputTag    fTPCTrackModuleLabel;
    art::InputTag    fCRTTrackModuleLabel;
};


sbnd::crt::CRTTrackMatching::CRTTrackMatching(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fMatchingAlg(p.get<fhicl::ParameterSet>("MatchingAlg"))
  , fTPCTrackModuleLabel(p.get<art::InputTag>("TPCTrackModuleLabel"))
  , fCRTTrackModuleLabel(p.get<art::InputTag>("CRTTrackModuleLabel"))
  {
    produces<std::vector<anab::T0>>();
    produces<art::Assns<recob::Track, anab::T0>>();
    produces<art::Assns<CRTTrack, anab::T0>>();
  }

void sbnd::crt::CRTTrackMatching::produce(art::Event& e)
{
  auto T0Vec          = std::make_unique<std::vector<anab::T0>>();
  auto tpcTrackT0Assn = std::make_unique<art::Assns<recob::Track, anab::T0>>();
  auto crtTrackT0Assn = std::make_unique<art::Assns<CRTTrack, anab::T0>>();

  art::Handle<std::vector<CRTTrack>> crtTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, crtTrackHandle);

  std::vector<art::Ptr<CRTTrack>> crtTrackVec;
  art::fill_ptr_vector(crtTrackVec, crtTrackHandle);

  art::Handle<std::vector<recob::Track>> tpcTrackHandle;
  e.getByLabel(fTPCTrackModuleLabel, tpcTrackHandle);

  std::vector<art::Ptr<recob::Track>> tpcTrackVec;
  art::fill_ptr_vector(tpcTrackVec, tpcTrackHandle);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  for(auto const &tpcTrack : tpcTrackVec)
    {
      MatchCandidate best = fMatchingAlg.GetBestMatchedCRTTrack(detProp, tpcTrack, crtTrackVec, e);

      if(best.valid)
        {
          T0Vec->push_back(anab::T0(best.time * 1e3, 0, tpcTrack->ID(), T0Vec->size(), best.score));
          util::CreateAssn(*this, e, *T0Vec, tpcTrack, *tpcTrackT0Assn);
          util::CreateAssn(*this, e, *T0Vec, best.thisTrack, *crtTrackT0Assn);
        }
    }

  e.put(std::move(T0Vec));
  e.put(std::move(tpcTrackT0Assn));
  e.put(std::move(crtTrackT0Assn));
}

DEFINE_ART_MODULE(sbnd::crt::CRTTrackMatching)
