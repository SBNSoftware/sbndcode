////////////////////////////////////////////////////////////////////////
// Class:       ArtDAQFragmentBlender
// Plugin Type: producer (Unknown Unknown)
// File:        ArtDAQFragmentBlender_module.cc
//
// Generated at Tue Sep  1 12:05:42 2026 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "artdaq-core/Data/Fragment.hh"
#include "TRandom.h"
#include "gallery/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

class ArtDAQFragmentBlender;


class ArtDAQFragmentBlender : public art::EDProducer {
public:
  explicit ArtDAQFragmentBlender(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ArtDAQFragmentBlender(ArtDAQFragmentBlender const&) = delete;
  ArtDAQFragmentBlender(ArtDAQFragmentBlender&&) = delete;
  ArtDAQFragmentBlender& operator=(ArtDAQFragmentBlender const&) = delete;
  ArtDAQFragmentBlender& operator=(ArtDAQFragmentBlender&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string fNoiseFileList;// Declare member data here.
  std::string fTPCDAQLabel;
  int fNumberNoiseFiles;
  TRandom randDraws;
  int TotalNoiseEvents;
  std::unique_ptr<gallery::Event> noiseGalleryEvent =new gallery::Event;
};


ArtDAQFragmentBlender::ArtDAQFragmentBlender(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fNoiseFileList = param.get<std::string>("NoiseFileList");
  fNumberNoiseFiles = param.get<int>("NumberNoiseFiles");
  fTPCDAQLabel = param.get<std::string>("TPCDAQLabel", "::");
  produces<std::vector<artdaq::Fragments>>();
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  int FileToGrab = fNumberNoiseFiles*randDraws.Uniform(1.0);
  noiseGalleryEvent* = gallery::Event event(allInputFiles);
  TotalNoiseEvents =  noiseGalleryEvent->numberOfEventsInFile();
}

void ArtDAQFragmentBlender::produce(art::Event& e)
{
  // Implementation of required member function here.
  int EventToScramble = TotalNoiseEvents*randDraws.Uniform(1.0);
  noiseGalleryEvent->goToEntry(EventToScramble);
  //Initialize fragment collection to put into the event
  std::unique_ptr<std::vector<artdaq::Fragment>> ScrambledFragments(new std::vector<artdaq::Fragment>);
  //Read out the TPC fragments from both events
  //nominal file
  art::Handle< std::vector<artdaq::Fragment> > NominalFragHandle;
  std::vector< art::Ptr<artdaq::Fragment> > NominalTPCfragmentList;
  if (e.getByLabel(fTPCDAQLabel,NominalFragHandle))
    art::fill_ptr_vector(NominalTPCfragmentList, NominalFragHandle);
  //Noise file
  art::Handle< std::vector<artdaq::Fragment> > NoiseFragHandle;
  std::vector< art::Ptr<artdaq::Fragment> > NoiseTPCfragmentList;
  NoiseFragHandle = noiseGalleryEvent->getValidHandle(fTPCDAQLabel);
  art::fill_ptr_vector(NominalTPCfragmentList, NoiseFragHandle);
  //Now mix up the entries in our scrambled vector
  //Will need to do a smarter thing in end to get the right mix
  for(int i=0; i<NoiseTPCfragmentList.size(); i++)
  {
    if(i%2==0) ScrambledFragments.push_back(NoiseTPCfragmentList[i]);
    else ScrambledFragments.push_back(NominalTPCfragmentList[i]);
  }
  //Add the new collection to the event
  event.put(std::move(ScrambledFragments));

}

DEFINE_ART_MODULE(ArtDAQFragmentBlender)
