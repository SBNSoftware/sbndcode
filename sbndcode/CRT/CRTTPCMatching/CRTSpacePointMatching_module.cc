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
#include "lardataobj/RecoBase/PFParticle.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

#include "sbndcode/CRT/CRTTPCMatching/CRTSpacePointMatchAlg.h"
#include "sbndcode/CRT/CRTEventDisplay/CRTEventDisplayAlg.h"

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

  CRTEventDisplayAlg    fCRTEventDisplayAlg;

  bool fMakeEventDisplays;
};


sbnd::crt::CRTSpacePointMatching::CRTSpacePointMatching(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fMatchingAlg(p.get<fhicl::ParameterSet>("MatchingAlg"))
  , fTPCTrackModuleLabel(p.get<art::InputTag>("TPCTrackModuleLabel"))
  , fCRTSpacePointModuleLabel(p.get<art::InputTag>("CRTSpacePointModuleLabel"))
  , fCRTEventDisplayAlg(p.get<fhicl::ParameterSet>("CRTEventDisplayAlg"))
  , fMakeEventDisplays(p.get<bool>("MakeEventDisplays"))
{
  if(!fMakeEventDisplays)
    produces<art::Assns<CRTSpacePoint, recob::Track, anab::T0>>();
}

void sbnd::crt::CRTSpacePointMatching::produce(art::Event& e)
{
  auto crtSPTPCTrackAssn = std::make_unique<art::Assns<CRTSpacePoint, recob::Track, anab::T0>>();

  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);

  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTPCTrackModuleLabel, trackHandle);

  std::vector<art::Ptr<recob::Track>> trackVec;
  art::fill_ptr_vector(trackVec, trackHandle);

  art::FindOneP<recob::PFParticle> tracksToPFPs(trackHandle, e, fTPCTrackModuleLabel);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  std::vector<SPMatchCandidate> candidates;

  for(auto const &track : trackVec)
    {
      const art::Ptr<recob::PFParticle> pfp = tracksToPFPs.at(track.key());

      if(pfp->PdgCode() != 13)
        continue;

      SPMatchCandidate closest = fMatchingAlg.GetClosestCRTSpacePoint(detProp, track, CRTSpacePointVec, e);

      if(closest.valid)
        candidates.push_back(closest);
    }

  std::sort(candidates.begin(), candidates.end(),
            [](const SPMatchCandidate &a, const SPMatchCandidate &b)
            { return a.score < b.score; });

  std::set<int> used_crt_sps;

  int nmatches = 0;

  for(auto const &candidate : candidates)
    {
      if(used_crt_sps.count(candidate.thisSP.key()) == 0)
        {
          const anab::T0 t0(candidate.time * 1e3, 0, 0, 0, candidate.score);
          crtSPTPCTrackAssn->addSingle(candidate.thisSP, candidate.thisTrack, t0);

          if(fMakeEventDisplays)
            {
              fCRTEventDisplayAlg.SetPrint(false);
              fCRTEventDisplayAlg.SetMinTime(candidate.time * 1e3 - 1e3);
              fCRTEventDisplayAlg.SetMaxTime(candidate.time * 1e3 + 1e3);
              auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
              fCRTEventDisplayAlg.Draw(clockData, e, Form("crtEventDisplayEvent%iSPMatch%i", e.event(), nmatches));
            }

          used_crt_sps.insert(candidate.thisSP.key());
          ++nmatches;
        }
    }

  if(!fMakeEventDisplays)
    e.put(std::move(crtSPTPCTrackAssn));
}

DEFINE_ART_MODULE(sbnd::crt::CRTSpacePointMatching)
