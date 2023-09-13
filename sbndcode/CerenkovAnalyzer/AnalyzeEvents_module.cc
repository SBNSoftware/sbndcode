////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Wed Oct 13 08:28:13 2021 by Edward Tyley using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"


#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"



#include "lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h"
#include "lardata/RecoBaseProxy/ProxyBase/withAssociated.h"
#include "lardata/RecoBaseProxy/ProxyBase/withParallelData.h"
#include "lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h"
#include "lardata/RecoBaseProxy/ProxyBase/getCollection.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/Simulation/SimPhotons.h"


//  artg4tk includes:
#include "artg4tk/pluginDetectors/gdml/PhotonHit.hh"

// Additional Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"


// ROOT includes
#include <TH1F.h>
#include <TTree.h>

// STL includes
#include <string>
#include <vector>


namespace test {
class AnalyzeEvents;
}

class test::AnalyzeEvents : public art::EDAnalyzer {
  public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  private:
  // Declare member data here.

  // Create out output tree
  TTree* fTree;
 
// Tree Variables

  std::vector<float> fEnergy;
  std::vector<float> fWlen;
  std::vector<float> fxPos;
  std::vector<float> fyPos;
  std::vector<float> fzPos;
  std::vector<float> fTime;
  std::vector<int> fID;
  std::vector<int> fProcess;
  bool fAnalyzePhotonHit;
};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
     : EDAnalyzer{p}  
    // Initialise out input labels by reading the fhicl parameters
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  typedef std::vector<art::Handle<artg4tk::PhotonHitCollection>> HandleVector;
  std::map<int, int> photonsperdet;
  auto allSims = e.getMany<artg4tk::PhotonHitCollection>();
  std::cout << "CheckPhotonHits Event:  " << e.event() << "  Is empty?: " << allSims.empty() << std::endl;

  for (HandleVector::const_iterator i = allSims.begin(); i != allSims.end(); ++i) {
    photonsperdet.clear();
    const artg4tk::PhotonHitCollection& sims(**i);
    if (sims.empty())
     {
      std::cout<<"-----> Sims is empty"<< std::endl;
     }
    for (artg4tk::PhotonHitCollection::const_iterator j = sims.begin(); j != sims.end(); ++j) {
      const artg4tk::PhotonHit& hit = *j;
      fEnergy.push_back(hit.GetEdep()*1000000);
      fxPos.push_back(hit.GetXpos());
      fyPos.push_back(hit.GetYpos());
      fzPos.push_back(hit.GetZpos());
      fID.push_back(hit.GetID());
      fTime.push_back(hit.GetTime());
      fWlen.push_back(0.001239847/hit.GetEdep());
      fProcess.push_back(hit.GetProcessID());
      }
    }

  fTree->Fill();
}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");
 
  fTree->Branch("True_Energy", &fEnergy);
  fTree->Branch("True_Wlen", &fWlen);
  fTree->Branch("True_xpos", &fxPos);
  fTree->Branch("True_ypos", &fyPos);
  fTree->Branch("True_zpos", &fzPos);
  fTree->Branch("ProcessID", &fProcess);
  fTree->Branch("True_time", &fTime);
  fTree->Branch("ID", &fID);
}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
