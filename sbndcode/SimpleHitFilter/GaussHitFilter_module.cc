////////////////////////////////////////////////////////////////////////
// Class:       GaussHitFilter
// Plugin Type: producer (Unknown Unknown)
// File:        GaussHitFilter_module.cc
//
// Generated at Tue Jun  3 15:57:00 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <vector>


class GaussHitFilter;


class GaussHitFilter : public art::EDProducer {
public:
  explicit GaussHitFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GaussHitFilter(GaussHitFilter const&) = delete;
  GaussHitFilter(GaussHitFilter&&) = delete;
  GaussHitFilter& operator=(GaussHitFilter const&) = delete;
  GaussHitFilter& operator=(GaussHitFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string             fHitProducer;
  std::vector<double> fMinHitHeight;

};


GaussHitFilter::GaussHitFilter(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fHitProducer    = p.get<std::string>   ("HitProducer");
  fMinHitHeight    = p.get<std::vector<double>>   ("minHitHeight");
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>> ();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void GaussHitFilter::produce(art::Event& e)
{
  // Implementation of required member function here.
  //Make vector of hits to save
  std::unique_ptr< std::vector<recob::Hit> > hit_v(std::make_unique<std::vector<recob::Hit>>());
  //Load in hit vector
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (e.getByLabel(fHitProducer,hitHandle))
    art::fill_ptr_vector(hitlist, hitHandle);
  //loop over entries in hitlist
  for(size_t i=0; i<hitlist.size(); i++){
    auto const& thisHit = hitlist[i];
    int PlaneIndex = (thisHit->WireID().Plane)%3;
    if(thisHit->PeakAmplitude() > fMinHitHeight[PlaneIndex])
    {
      hit_v->push_back(*thisHit);
    }
  }
  //Save filtered hits
  std::cout << " on this event I eliminated " << hitlist.size()-hit_v->size() << " hits " 
  << double(hitlist.size()-hit_v->size())/hitlist.size() << " %" << std::endl;
  e.put(std::move(hit_v));
}

DEFINE_ART_MODULE(GaussHitFilter)
