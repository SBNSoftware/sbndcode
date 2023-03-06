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
  std::vector<int> fNumberPhotons;
  std::vector<int> fNumberCPions;
  std::vector<int> fNumberNPions;

  // General configuration parameters
  std::string fOutputFileName;
  std::string fAllFileName;
  bool fVerbose;

  // Setup the file to write information to
  std::ofstream fOutputFile;
  std::ofstream fAllFile;

  // Counters
  unsigned int nEv;
  unsigned int nFile;
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
  fNumberPhotons         (p.get< std::vector< int > >("NumberPhotons")),
  fNumberCPions          (p.get< std::vector< int > >("NumberCPions")),
  fNumberNPions          (p.get< std::vector< int > >("NumberNPions")),
  fOutputFileName        (p.get< std::string >("OutputFileName")),
  fAllFileName           (p.get< std::string >("AllFileName")),
  fVerbose               (p.get< bool >("Verbose"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void evd::HandScanalyser::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  nEv++;

  // Print file numbers to file
  // Get the file ID based on 100 events per file
  if((std::floor(nEv/100.)-nFile+1) == 1){
    fOutputFile << " ----------------------------------------------------------- " << std::endl;
    fOutputFile << " File: " << nFile << std::endl;
    fOutputFile << " ----------------------------------------------------------- " << std::endl;
    fAllFile << " ----------------------------------------------------------- " << std::endl;
    fAllFile << " File: " << nFile << std::endl;
    fAllFile << " ----------------------------------------------------------- " << std::endl;
    nFile++;
  }

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mctList;
  if (e.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mctList, mctruthListHandle);

  // * MC particle information
  art::Handle< std::vector<simb::MCParticle> > mcparticleListHandle;
  std::vector<art::Ptr<simb::MCParticle> > mcpList;
  if (e.getByLabel(fG4ModuleLabel,mcparticleListHandle))
    art::fill_ptr_vector(mcpList, mcparticleListHandle);

  // Counters for particle types in the event
  int nMu = 0;
  int nEl = 0;
  int nPr = 0;
  int nPh = 0;
  int nCp = 0;
  int nNp = 0;

  // Loop over MCTruth objects and get associated MCParticles
  for(art::Ptr<simb::MCTruth> mct: mctList){
    // Make sure the origin is neutrino
    if(mct->Origin() != simb::kBeamNeutrino){
      std::cout << " Origin: " << mct->Origin() << std::endl;
      continue;
    }

    // Now loop over the MCParticles and check how many primary particles of each type there are
    for(int p = 0; p < mct->NParticles(); ++p){
      simb::MCParticle mcp = mct->GetParticle(p);

      // Make sure the parent particle is primary
      for(art::Ptr<simb::MCParticle> mcp1: mcpList){
        if(mcp.Mother() != mcp1->TrackId()) continue;
        if(mcp1->Process() != "primary") continue;
      }
      int pdg = mcp.PdgCode();

      if(abs(pdg) == 13)   nMu++;
      if(abs(pdg) == 11)   nEl++;
      if(abs(pdg) == 2212) nPr++;
      if(abs(pdg) == 22)   nPh++;
      if(abs(pdg) == 211)  nCp++;
      if(abs(pdg) == 111)  nNp++;

    } // mcp
   
    // Print the contents of the event
    if((nMu+nEl+nPr+nPh+nCp+nNp) == 0) continue;
    fAllFile << " Event " << std::setw(4) 
             << nEv << " has " 
             << nMu << " muons, " 
             << nEl << " electrons, " 
             << nPr << " protons, " 
             << nPh << " photons, " 
             << nCp << " charged pions and " 
             << nNp << " neutral pions " << std::endl;

    // Now loop over the configuration vectors and see if it's a match
    bool topologyFound = false;
    for(unsigned int i = 0; i < fTopologyLabels.size(); ++i){
      if((nMu == fNumberMuons.at(i) || fNumberMuons.at(i) == -1) && 
         (nEl == fNumberElectrons.at(i) || fNumberElectrons.at(i) == -1) &&
         (nPr == fNumberProtons.at(i) || fNumberProtons.at(i) == -1) &&
         (nPh == fNumberPhotons.at(i) || fNumberPhotons.at(i) == -1) &&
         (nCp == fNumberCPions.at(i) || fNumberCPions.at(i) == -1) &&
         (nNp == fNumberNPions.at(i) || fNumberNPions.at(i) == -1)){

        fOutputFile << " Event " << std::setw(4) << nEv << " is of the topology: " << std::setw(10) << fTopologyLabels.at(i) << std::endl;
        topologyFound = true;
      } // If topology matches
    } // Topology labels

    for(int p = 0; p < mct->NParticles(); ++p){
      simb::MCParticle mcp = mct->GetParticle(p);
      for(art::Ptr<simb::MCParticle> mcp1: mcpList){
        if(mcp.Mother() != mcp1->TrackId()) continue;
        if(mcp1->Process() != "primary") continue;
      }
      int pdg = mcp.PdgCode();

      if(topologyFound){
        // Now get the particle's energies and print them too
        if(nMu != 0 && abs(pdg) == 13)
          fOutputFile << "   Muon energy     = " << mcp.E() << " GeV " << std::endl;
        if(nEl != 0 && abs(pdg) == 11)
          fOutputFile << "   Electron energy = " << mcp.E() << " GeV " << std::endl;
        if(nPr != 0 && abs(pdg) == 2212)
          fOutputFile << "   Proton energy   = " << mcp.E() << " GeV " << std::endl;
        if(nPh != 0 && abs(pdg) == 22)
          fOutputFile << "   Photon energy   = " << mcp.E() << " GeV " << std::endl;
        if(nCp != 0 && abs(pdg) == 211)
          fOutputFile << "   Pi^{+/-} energy = " << mcp.E() << " GeV " << std::endl;
        if(nNp != 0 && abs(pdg) == 111)
          fOutputFile << "   Pi^{0} energy   = " << mcp.E() << " GeV " << std::endl;
      } // Topology found
      // Now get the particle's energies and print them too
      if(nMu != 0 && abs(pdg) == 13)
        fAllFile << "   Muon energy     = " << mcp.E() << " GeV " << std::endl;
      if(nEl != 0 && abs(pdg) == 11)
        fAllFile << "   Electron energy = " << mcp.E() << " GeV " << std::endl;
      if(nPr != 0 && abs(pdg) == 2212)
        fAllFile << "   Proton energy   = " << mcp.E() << " GeV " << std::endl;
      if(nPh != 0 && abs(pdg) == 22)
        fAllFile << "   Photon energy   = " << mcp.E() << " GeV " << std::endl;
      if(nCp != 0 && abs(pdg) == 211)
        fAllFile << "   Pi^{+/-} energy = " << mcp.E() << " GeV " << std::endl;
      if(nNp != 0 && abs(pdg) == 111)
        fAllFile << "   Pi^{0} energy   = " << mcp.E() << " GeV " << std::endl;
    } // MCP
    if(topologyFound) fOutputFile << std::endl;
    fAllFile << std::endl;
  } // MCTruth
} // Analyser

void evd::HandScanalyser::beginJob()
{
  // Implementation of optional member function here.
  //
  // Open a csv file to output the relevant information
  fOutputFile.open(fOutputFileName.c_str());
  fAllFile.open(fAllFileName.c_str());

  // File and event counter
  nFile = 0;
  nEv = 0;
}

void evd::HandScanalyser::endJob()
{
  // Implementation of optional member function here.
  fOutputFile.close();
  fAllFile.close();
}

DEFINE_ART_MODULE(evd::HandScanalyser)
