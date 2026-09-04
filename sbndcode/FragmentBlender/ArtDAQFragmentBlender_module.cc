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
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "artdaq-core/Data/Fragment.hh"
#include "TRandom.h"
#include "gallery/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "sbndcode/Decoders/TPC/SBNDTPCDecoder.h"

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
  daq::SBNDTPCDecoder tpcDecoderBusiness;
  std::string fNoiseFileList;// Declare member data here.
  std::string fTPCDAQLabel;
  int fNumberNoiseFiles;
  TRandom randDraws;
  int TotalNoiseEvents;
  std::unique_ptr<gallery::Event> noiseGalleryEvent;
};


ArtDAQFragmentBlender::ArtDAQFragmentBlender(fhicl::ParameterSet const& p)
  : EDProducer{p} ,
    tpcDecoderBusiness{p}
  // More initializers here.
{
  fNoiseFileList = p.get<std::string>("NoiseFileList");
  fNumberNoiseFiles = p.get<int>("NumberNoiseFiles");
  fTPCDAQLabel = p.get<std::string>("TPCDAQLabel", "::");
  produces<std::vector<raw::RawDigit>>("test");
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  unsigned int FileToGrab = fNumberNoiseFiles*randDraws.Uniform(1.0);
  std::vector<std::string> allInputFiles;
  std::ifstream file(fNoiseFileList);
  std::string line;
  while (getline(file, line)) {
            allInputFiles.push_back(line);
        }
  file.close();
  std::vector<std::string> TempOneFile;
  TempOneFile.push_back(allInputFiles[FileToGrab]);
  std::cout << " grabbing "<< allInputFiles[FileToGrab] << " as noise file " << std::endl;
  noiseGalleryEvent.reset( new gallery::Event(TempOneFile) );
  TotalNoiseEvents =  noiseGalleryEvent->numberOfEventsInFile();
}

void ArtDAQFragmentBlender::produce(art::Event& e)
{
  // Implementation of required member function here.
  int EventToScramble = TotalNoiseEvents*randDraws.Uniform(1.0);
  noiseGalleryEvent->goToEntry(EventToScramble);
  //Initialize fragment collection to put into the event
  //Update to recob wire
  //There is a small chance everything breaks and we need to save all the items from the decoder
  std::unique_ptr<std::vector<raw::RawDigit>> ScrambledFragments(new std::vector<raw::RawDigit>);
  std::unique_ptr<std::vector<raw::RawDigit>> NominalFragments = tpcDecoderBusiness.produce2(e);
  std::unique_ptr<std::vector<raw::RawDigit>> NoiseFragments = tpcDecoderBusiness.produce2(noiseGalleryEvent);
  //Loop through fragments and find the right IDs to add in 
  //Should update this logic to use ranges of channel IDs
  for(int i=0; i<int(NominalFragments.size()); i++)
  {
    int ChannelID = NominalFragments[i].Channel();
    if(   (ChannelID>=3968 && ChannelID<5632) || (ChannelID>=9600 && ChannelID<11264)  ) //grab collection clusters
    {
      ScrambledFragments->push_back(NominalFragments[i]);
    }
    else //must be induction
    {
      ScrambledFragments->push_back(NoiseFragments[i]);
    }
  }
  std::cout << "Scrambled size " << ScrambledFragments->size() << std::endl;
  //Add the new collection to the event
  std::cout << " adding in scrambled fragments " << std::endl;
  e.put(std::move(ScrambledFragments), "test");
}

DEFINE_ART_MODULE(ArtDAQFragmentBlender)
