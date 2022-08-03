////////////////////////////////////////////////////////////////////////
// Class:       HandScanalyser
// Plugin Type: analyzer (Unknown Unknown)
// File:        HandScanalyser_module.cc
//
// Generated at Wed Aug  3 09:21:28 2022 by Rhiannon Jones using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art objects
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"

// C++ objects
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace evd {
  class HandScanalyser;
}


class evd::HandScanalyser : public art::EDAnalyzer {
public:
  explicit HandScanalyser(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HandScanalyser(HandScanalyser const&) = delete;
  HandScanalyser(HandScanalyser&&) = delete;
  HandScanalyser& operator=(HandScanalyser const&) = delete;
  HandScanalyser& operator=(HandScanalyser&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  
  // Input Tags
  std::string fGenieGenModuleLabel;
  std::string fG4ModuleLabel;
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fVertexModuleLabel;

  // Analysis configuration parameters
  std::vector<std::string> fTopologyLabels; 
  std::vector<int> fNumberMuons;  
  std::vector<int> fNumberElectrons;
  std::vector<int> fNumberProtons;
  std::vector<int> fNumberCPions;
  std::vector<int> fNumberNPions;

  // General configuration parameters
  std::string fOutputFileName;

  // Setup the file to write information to
  std::ofstream fOutputFile;
};


evd::HandScanalyser::HandScanalyser(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
  // More initializers here.
  fGenieGenModuleLabel   (p.get< std::string >("GenieGenModuleLabel")),
  fG4ModuleLabel         (p.get< std::string >("G4ModuleLabel")),
  fTrackModuleLabel      (p.get< std::string >("TrackModuleLabel")),
  fShowerModuleLabel     (p.get< std::string >("ShowerModuleLabel")),
  fPFParticleModuleLabel (p.get< std::string >("PFParticleModuleLabel")),
  fVertexModuleLabel     (p.get< std::string >("VertexModuleLabel")),
  fTopologyLabels        (p.get< std::vector< std::string > >("TopologyLabels")),
  fNumberMuons           (p.get< std::vector< int > >("NumberMuons")),
  fNumberElectrons       (p.get< std::vector< int > >("NumberElectrons")),
  fNumberProtons         (p.get< std::vector< int > >("NumberProtons")),
  fNumberCPions          (p.get< std::vector< int > >("NumberCPions")),
  fNumberNPions          (p.get< std::vector< int > >("NumberNPions")),
  fOutputFileName        (p.get< std::string >("OutputFileName"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void evd::HandScanalyser::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fOutputFile << " Hello there kind madame" << std::endl;
}

void evd::HandScanalyser::beginJob()
{
  // Implementation of optional member function here.
  //
  // Open a csv file to output the relevant information
  fOutputFile.open(fOutputFileName.c_str());
}

void evd::HandScanalyser::endJob()
{
  // Implementation of optional member function here.
  fOutputFile.close();
}

DEFINE_ART_MODULE(evd::HandScanalyser)
