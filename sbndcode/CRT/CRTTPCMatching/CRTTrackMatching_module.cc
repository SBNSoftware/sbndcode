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
#include "lardataobj/RecoBase/PFParticle.h"

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
    produces<art::Assns<CRTTrack, recob::Track, anab::T0>>();
  }

void sbnd::crt::CRTTrackMatching::produce(art::Event& e)
{
  auto crtTrackTPCTrackAssn = std::make_unique<art::Assns<CRTTrack, recob::Track, anab::T0>>();

  art::Handle<std::vector<CRTTrack>> crtTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, crtTrackHandle);

  std::vector<art::Ptr<CRTTrack>> crtTrackVec;
  art::fill_ptr_vector(crtTrackVec, crtTrackHandle);

  art::Handle<std::vector<recob::Track>> tpcTrackHandle;
  e.getByLabel(fTPCTrackModuleLabel, tpcTrackHandle);

  std::vector<art::Ptr<recob::Track>> tpcTrackVec;
  art::fill_ptr_vector(tpcTrackVec, tpcTrackHandle);

  art::FindOneP<recob::PFParticle> tracksToPFPs(tpcTrackHandle, e, fTPCTrackModuleLabel);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  std::vector<TrackMatchCandidate> candidates;

  for(auto const &tpcTrack : tpcTrackVec)
    {
      const art::Ptr<recob::PFParticle> pfp = tracksToPFPs.at(tpcTrack.key());

      if(pfp->PdgCode() != 13)
        continue;

      TrackMatchCandidate best = fMatchingAlg.GetBestMatchedCRTTrack(detProp, tpcTrack, crtTrackVec, e);

      if(best.valid)
        candidates.push_back(best);
    }

  std::sort(candidates.begin(), candidates.end(),
            [](const TrackMatchCandidate &a, const TrackMatchCandidate &b)
            { return a.score < b.score; });

  std::set<int> used_crt_tracks;

  for(auto const &candidate : candidates)
    {
      if(used_crt_tracks.count(candidate.thisCRTTrack.key()) == 0)
        {
          const anab::T0 t0(candidate.time * 1e3, 0, 0, 0, candidate.score);
          crtTrackTPCTrackAssn->addSingle(candidate.thisCRTTrack, candidate.thisTPCTrack, t0);

          used_crt_tracks.insert(candidate.thisCRTTrack.key());
        }
    }

  e.put(std::move(crtTrackTPCTrackAssn));
}

DEFINE_ART_MODULE(sbnd::crt::CRTTrackMatching)
