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

  std::string fNoiseFileList;// Declare member data here.
  std::string fTPCDAQLabel;
  int fNumberNoiseFiles;
  TRandom randDraws;
  int TotalNoiseEvents;
  std::unique_ptr<gallery::Event> noiseGalleryEvent;
  std::vector<int> NominalFragmentsToTake;
  std::vector<int> NoiseFragmentsToTake;
};


ArtDAQFragmentBlender::ArtDAQFragmentBlender(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
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
  NominalFragmentsToTake={};
  NoiseFragmentsToTake={37123, 37124, 37125, 37126, 37127, 37128, 37129, 37130, 37131,
       37132, 37133, 37134, 37135, 37136, 37137, 37138, //crate 01
       37379, 37380, 37381, 37382, 37383, 37384, 37385, 37386, 37387,
       37388, 37389, 37390, 37391, 37392, 37393, 37394, //crate 02
       37635, 37636, 37637, 37638, 37639, 37640, 37641, 37642, 37643,
       37644, 37645, 37646, 37647, 37648, 37649, 37650, //crate 03
       37891, 37892, 37893, 37894, 37895, 37896, 37897, 37898, 37899,
       37900, 37901, 37902, 37903, 37904, 37905, 37906, //crate 04
       38147, 38148, 38149, 38150, 38151, 38152, 38153, 38154, 38155,
       38156, 38157, 38158, 38159, 38160, 38161, 38162, //crate 05
       38403, 38404, 38405, 38406, 38407, 38408, 38409, 38410, 38411,
       38412, 38413, 38414, 38415, 38416, 38417, 38418, //crate 06
       38659, 38660, 38661, 38662, 38663, 38664, 38665, 38666, 38667,
       38668, 38669, 38670, 38671, 38672, 38673, 38674, //crate07
       38915, 38916, 38917, 38918, 38919, 38920, 38921, 38922, 38923,
       38924, 38925, 38926, 38927, 38928, 38929, //crate08
       39171, 39172, 39173, 39174, 39175, 39176, 39177, 39178, 39179,
       39180, 39181, 39182, 39183, 39184, 39185, 39186, //crate09
       39427, 39428, 39429, 39430, 39431, 39432, 39433, 39434, 39435,
       39436, 39437, 39438, 39439, 39440, 39441, 39442, //crate10
       39683, 39684, 39685, 39686, 39687, 39688, 39689, 39690, 39691,
       39692, 39693, 39694, 39695, 39696, 39697, 39698 //crate11
     }
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
  std::unique_ptr<std::vector<raw::RawDigit>> NominalFragments = daq::SBNDTPCDecoder::produce2(e);
  std::unique_ptr<std::vector<raw::RawDigit>> NoiseFragments = daq::SBNDTPCDecoder::produce2(noiseGalleryEvent);
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
