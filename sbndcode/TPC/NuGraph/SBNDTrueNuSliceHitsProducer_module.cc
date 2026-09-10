/**
 * @file   sbndcode/TPC/NuGraph/SBNDTrueNuSliceHitsProducer_module.cc
 * @brief  Implementation of `SBNDTrueNuSliceHitsProducer` module
 * @author Alexander Antonakis and Giuseppe Cerati (cerati@fnal.gov)
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <algorithm> // std::find()
#include <cmath> // std::abs()
#include <memory>
#include <unordered_map>
#include <utility> // std::move()
#include <vector>

#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

class SBNDTrueNuSliceHitsProducer : public art::EDProducer {

/**
 * @brief Produce separate collections of hits and space point from the slice with most true hits from the neutrino interaction.
 *
 * Produce separate collections of hits and space point from the slice with most true hits from the neutrino interaction.
 * It also produces a vector of ints, that maps the new hit collection to the original one.
 * Optionally, it can store the association to the trigger flash matched to the slice.
 *
 */

public:
  explicit SBNDTrueNuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDTrueNuSliceHitsProducer(SBNDTrueNuSliceHitsProducer const&) = delete;
  SBNDTrueNuSliceHitsProducer(SBNDTrueNuSliceHitsProducer&&) = delete;
  SBNDTrueNuSliceHitsProducer& operator=(SBNDTrueNuSliceHitsProducer const&) = delete;
  SBNDTrueNuSliceHitsProducer& operator=(SBNDTrueNuSliceHitsProducer&&) = delete;


private:

  // Declare member data here.
  art::InputTag const fBaryMatchLabel;
  art::InputTag const fSliceLabel;
  art::InputTag const fSpsLabel;
  bool const doFlashAssn;

  // Required functions.
  void produce(art::Event& e) override;

};


SBNDTrueNuSliceHitsProducer::SBNDTrueNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fBaryMatchLabel( p.get<art::InputTag>("BaryMatchLabel") ),
  fSliceLabel(p.get<art::InputTag>("SliceLabel")),
  fSpsLabel(p.get<art::InputTag>("SpsLabel")),
  doFlashAssn(p.get<bool>("DoFlashAssn",false))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  produces<std::vector<unsigned int>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  if (doFlashAssn) produces<art::Assns<recob::Slice, recob::OpFlash>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SBNDTrueNuSliceHitsProducer::produce(art::Event& e)
{

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto assocSliceHitKeys = std::make_unique<std::vector<unsigned int> >();

  art::PtrMaker<recob::Hit> hitPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spsPtrMaker(e);

  auto outputSpacepoints = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputSpsHitAssns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  auto outputSlcFlhAssns = std::make_unique<art::Assns<recob::Slice, recob::OpFlash>>();

  auto const* pi = lar::providerFrom<cheat::ParticleInventoryService>();
  auto const* bt = lar::providerFrom<cheat::BackTrackerService>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSliceLabel);
  art::FindManyP<recob::Hit> assocSliceHit(inputSlice, e, fSliceLabel);
  int slkey = -1;
  int maxCnt = 0;
  // find the slice with the largest number of hits contributed by neutrino interactions
  for (size_t sk=0; sk<inputSlice->size(); sk++) {
    int slhcnt = 0;
    for (auto hit : assocSliceHit.at(sk)) {
      std::vector<sim::TrackIDE> const& h_ides = bt->ChannelToTrackIDEs(clockData, hit->Channel(), hit->StartTick(), hit->EndTick());
      for (auto const& tide : h_ides) {
      	int tid = tide.trackID;
	auto mct = pi->TrackIdToMCTruth_P(abs(tid));
	bool fromNu = mct.isAvailable() && mct->NeutrinoSet();
	if (fromNu) {
	  slhcnt++;
	  break;
	}
      }
    }
    if (slhcnt>maxCnt) {
      slkey = sk;
      maxCnt = slhcnt;
    }
  }
  mf::LogTrace{ "SBNDTrueNuSliceHitsProducer" } << "keys: slice=" << slkey;

  if (slkey>=0) {
    auto const& hits = assocSliceHit.at(slkey);
    mf::LogTrace{ "SBNDTrueNuSliceHitsProducer" } << "slice hits=" << hits.size();
    std::unordered_map<size_t,size_t> keyIdxMap;
    for (size_t ih=0; ih<hits.size(); ih++) {
      assocSliceHitKeys->push_back(hits[ih].key());
      keyIdxMap.insert({hits[ih].key(), ih});
    }
    std::sort(begin(*assocSliceHitKeys), end(*assocSliceHitKeys));
    outputHits->reserve(assocSliceHitKeys->size());
    for (auto const& key: *assocSliceHitKeys) outputHits->push_back(*hits[keyIdxMap.at(key)]);
  }

  // Get spacepoints from the event record
  auto splist = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpsLabel);
  size_t countsps = 0;
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> spacePointToHits(splist, e, fSpsLabel);
  for (size_t spIdx = 0; spIdx < splist->size(); ++spIdx) {
    auto const& assochits = spacePointToHits.at(spIdx);
    // Consider only space points with hits associated on all planes. This is enough for NuGraph.
    if (assochits.size()<3) continue;
    //
    std::vector<size_t> hitpos;
    for (size_t j = 0; j < assochits.size(); ++j) {
      auto pos = std::lower_bound(assocSliceHitKeys->begin(),assocSliceHitKeys->end(),assochits[j].key());
      if ( (pos == assocSliceHitKeys->end()) || (*pos != assochits[j].key()) ) break;
      hitpos.push_back(pos-assocSliceHitKeys->begin());
    }
    if (hitpos.size()<3) continue;
    outputSpacepoints->emplace_back((*splist)[spIdx]);
    //
    const art::Ptr<recob::SpacePoint> sps = spsPtrMaker(outputSpacepoints->size()-1);
    for (size_t j = 0; j < hitpos.size(); ++j) {
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(hitpos[j]);
      outputSpsHitAssns->addSingle(ahp,sps);
    }
    countsps++;
  } // for spacepoint
  mf::LogTrace{ "SBNDTrueNuSliceHitsProducer" } << "sps count=" << countsps;

  art::ValidHandle< std::vector< sbn::TPCPMTBarycenterMatch > > baryMatchListHandle = e.getValidHandle<std::vector<sbn::TPCPMTBarycenterMatch> >(fBaryMatchLabel);
  art::FindManyP<recob::Slice> msl(baryMatchListHandle, e, fBaryMatchLabel);
  art::FindManyP<recob::OpFlash> mfl(baryMatchListHandle, e, fBaryMatchLabel);
  float minRad = 99999.;
  int mbkey = -1;
  for (size_t ibm=0; ibm<baryMatchListHandle->size(); ++ibm) {
    art::Ptr<sbn::TPCPMTBarycenterMatch> bm(baryMatchListHandle,ibm);
    if (int(msl.at(ibm).at(0).key())!=slkey) continue;
    if (mfl.at(ibm).size()>0) {
      mf::LogTrace{ "SBNDTrueNuSliceHitsProducer" } << "ibm=" << ibm << " radius=" << bm->radius << " radius_Trigger=" << bm->radius_Trigger << " flashTime=" << bm->flashTime
						  << " slkey=" << msl.at(ibm).at(0).key() << " center=" << bm->chargeCenter
						  << " flkey=" << mfl.at(ibm).at(0).key() << " center=" << bm->flashCenter;
      if (bm->radius<minRad) {
	//this is the triggering flash
	mbkey = ibm;
	minRad = bm->radius;
      }
    }
  }
  mf::LogTrace{ "SBNDTrueNuSliceHitsProducer" } << "mbkey=" << mbkey;
  if (doFlashAssn){
    if (slkey>=0 && mbkey>=0) {
      art::Ptr<recob::Slice> slp(inputSlice,slkey);
      outputSlcFlhAssns->addSingle(slp,mfl.at(mbkey).at(0));
    }
  }

  e.put(std::move(outputHits));
  e.put(std::move(assocSliceHitKeys));
  e.put(std::move(outputSpacepoints));
  e.put(std::move(outputSpsHitAssns));
  if (doFlashAssn) e.put(std::move(outputSlcFlhAssns));
}

DEFINE_ART_MODULE(SBNDTrueNuSliceHitsProducer)