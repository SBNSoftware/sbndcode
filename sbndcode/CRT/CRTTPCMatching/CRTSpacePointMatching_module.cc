/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTSpacePointMatching
/// Module Type: producer
/// File:        CRTSpacePointMatching_module.cc
///
/// Author:         Thomas Brooks
/// E-mail address: tbrooks@fnal.gov
///
/// Modified from CRTT0Matching by Thomas Warburton.
/////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

#include "sbndcode/CRT/CRTTPCMatching/CRTSpacePointMatchAlg.h"

#include <memory>

namespace sbnd::crt {
  class CRTSpacePointMatching;
}

class sbnd::crt::CRTSpacePointMatching : public art::EDProducer {
public:

  explicit CRTSpacePointMatching(fhicl::ParameterSet const& p);

  CRTSpacePointMatching(CRTSpacePointMatching const&) = delete;
  CRTSpacePointMatching(CRTSpacePointMatching&&) = delete;
  CRTSpacePointMatching& operator=(CRTSpacePointMatching const&) = delete;
  CRTSpacePointMatching& operator=(CRTSpacePointMatching&&) = delete;

  void produce(art::Event& e) override;

private:

  CRTSpacePointMatchAlg fMatchingAlg;
  art::InputTag         fTPCTrackModuleLabel;
  art::InputTag         fCRTSpacePointModuleLabel;
};


sbnd::crt::CRTSpacePointMatching::CRTSpacePointMatching(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fMatchingAlg(p.get<fhicl::ParameterSet>("MatchingAlg"))
  , fTPCTrackModuleLabel(p.get<art::InputTag>("TPCTrackModuleLabel"))
  , fCRTSpacePointModuleLabel(p.get<art::InputTag>("CRTSpacePointModuleLabel"))
  {
    produces<std::vector<anab::T0>>();
    produces<art::Assns<recob::Track, anab::T0>>();
    produces<art::Assns<CRTSpacePoint, anab::T0>>();
  }

void sbnd::crt::CRTSpacePointMatching::produce(art::Event& e)
{
  auto T0Vec       = std::make_unique<std::vector<anab::T0>>();
  auto trackT0Assn = std::make_unique<art::Assns<recob::Track, anab::T0>>();
  auto crtSPT0Assn = std::make_unique<art::Assns<CRTSpacePoint, anab::T0>>();

  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);

  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTPCTrackModuleLabel, trackHandle);

  std::vector<art::Ptr<recob::Track>> trackVec;
  art::fill_ptr_vector(trackVec, trackHandle);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  for(auto const &track : trackVec)
    {
      MatchCandidate closest = fMatchingAlg.GetClosestCRTSpacePoint(detProp, track, CRTSpacePointVec, e);

      if(closest.score >= 0)
        {
          T0Vec->push_back(anab::T0(closest.time*1e3, 0, track->ID(),  T0Vec->size(), closest.score));
          util::CreateAssn(*this, e, *T0Vec, track, *trackT0Assn);
          util::CreateAssn(*this, e, *T0Vec, closest.thisSP, *crtSPT0Assn);
        }
    }

  e.put(std::move(T0Vec));
  e.put(std::move(trackT0Assn));
  e.put(std::move(crtSPT0Assn));
}

DEFINE_ART_MODULE(sbnd::crt::CRTSpacePointMatching)
