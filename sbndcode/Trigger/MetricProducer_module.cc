////////////////////////////////////////////////////////////////////////
// Class:       MetricProducer
// Plugin Type: producer (Unknown Unknown)
// File:        MetricProducer_module.cc
//
// Michelle Stancari
// August 2022
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

//#include "sbndaq-artdaq-core/Obj/SBND/CRTmetric.h"
#include "sbnobj/SBND/Trigger/CRTmetric.hh"

#include <memory>
#include <iostream>

namespace sbndaq {
  class MetricProducer;
}

class sbndaq::MetricProducer : public art::EDProducer {
public:
  explicit MetricProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MetricProducer(MetricProducer const&) = delete;
  MetricProducer(MetricProducer&&) = delete;
  MetricProducer& operator=(MetricProducer const&) = delete;
  MetricProducer& operator=(MetricProducer&&) = delete;

  // Required functions.
  void produce(art::Event& evt) override;

private:

  // fhicl parameters
  art::Persistable is_persistable_;
  int fBeamWindowStart;
  int fBeamWindowEnd;
  bool fVerbose;

  // event information
  int fRun;
  art::EventNumber_t fEvent;

  //metric variables
  int hitsperplane[7];


  void analyze_crt_fragment(artdaq::Fragment & frag);

};


sbndaq::MetricProducer::MetricProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  is_persistable_(p.get<bool>("is_persistable", true) ? art::Persistable::Yes : art::Persistable::No),
  fBeamWindowStart(p.get<int>("BeamWindowStart",320000)),
  fBeamWindowEnd(p.get<int>("BeamWindowEnd",350000)),
  fVerbose(p.get<bool>("Verbose",false))
  {
  // Call appropriate produces<>() functions here.
  produces< sbndaq::CRTmetric >("", is_persistable_);

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


void sbndaq::MetricProducer::produce(art::Event& evt)
{

  // load event information
  int fRun = evt.run();
  art::EventNumber_t fEvent = evt.event();

  if (fVerbose) std::cout << "Processing: Run " << fRun << ", Event " << fEvent << std::endl;

  // object to store required trigger information in
  std::unique_ptr<sbndaq::CRTmetric> CRTMetricInfo = std::make_unique<sbndaq::CRTmetric>();

  // clear variables at the beginning of the event
  // move this to constructor??
  for (int ip=0;ip<7;++ip)  { CRTMetricInfo->hitsperplane[ip]=0; hitsperplane[ip]=0;}
  //std::fill(hitsperplane, hitsperplane + sizeof(hitsperplane), 0);

  // get fragment handles
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();

  // loop over fragment handles
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      // container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() != sbndaq::detail::FragmentType::BERNCRTV2) continue;
	if (fVerbose)     std::cout << "    Found " << contf.block_count() << " CRT Fragments in container " << std::endl;
	for (size_t ii = 0; ii < contf.block_count(); ++ii) analyze_crt_fragment(*contf[ii].get());
      }
    }
    else {
      // normal fragment
      if (handle->front().type()!=sbndaq::detail::FragmentType::BERNCRTV2) continue;
      if (fVerbose)   std::cout << "   found CRT fragments " << handle->size() << std::endl;
      for (auto frag : *handle)	analyze_crt_fragment(frag);
    }
  } // loop over frag handles

/*
  // determine relevant metrics
  fNAboveThreshold = checkCoincidence();

  // add to trigger object
  CRTMetricInfo->nAboveThreshold = fNAboveThreshold;
  CRTMetricInfo->PMT_chA = fPMT_chA;
  CRTMetricInfo->PMT_chB = fPMT_chB;
  CRTMetricInfo->PMT_chC = fPMT_chC;
  CRTMetricInfo->PMT_chD = fPMT_chD;
  CRTMetricInfo->PMT_chE = fPMT_chE;
*/

  for (int i=0;i<7;++i) {CRTMetricInfo->hitsperplane[i] = hitsperplane[i];}

  if (fVerbose) {
    std::cout << "CRT hit count during beam spill ";
    for (int i=0;i<7;++i) std::cout << i << " " << CRTMetricInfo->hitsperplane[i] ;
    std::cout << std::endl;
  }

  // add to event
  evt.put(std::move(CRTMetricInfo));

}



void sbndaq::MetricProducer::analyze_crt_fragment(artdaq::Fragment & frag)
{


  sbndaq::BernCRTFragmentV2 bern_fragment(frag);

  // use  fragment ID to get plane information
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  //frag.sequenceID()
  auto thisone = frag.fragmentID();  uint plane = (thisone & 0x0700) >> 8;
  if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}

  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
    // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
    auto thisflag = bevt->flags;
    if (thisflag & 0x2 && !(thisflag & 0xC) ) {
      // check ts1 for beam window
      auto thistime=bevt->ts1;
      if ((int)thistime>fBeamWindowStart && (int)thistime<fBeamWindowEnd) hitsperplane[plane]++;
      //CRTMetricInfo->hitsperplane[plane]++;
    }
  }



}

// -------------------------------------------------

// -------------------------------------------------

DEFINE_ART_MODULE(sbndaq::MetricProducer)
