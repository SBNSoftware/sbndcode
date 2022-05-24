////////////////////////////////////////////////////////////////////////
// Class:       NCReco
// Plugin Type: analyzer (art v3_06_03)
// File:        NCReco_module.cc
//
// Generated at Tue Jul 20 11:25:17 2021 by Ross MacFadyen using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"

//C++ Includes
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

namespace ana {
  class NCReco;
}

class ana::NCReco : public art::EDAnalyzer {
public:
  explicit NCReco(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCReco(NCReco const&) = delete;
  NCReco(NCReco&&) = delete;
  NCReco& operator=(NCReco const&) = delete;
  NCReco& operator=(NCReco&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void endJob() override;
  void beginJob() override;


private:
    //fcl parameters
    string fGenieGenModuleLabel;
    string fLArGeantModuleLabel;
    string fMCTruthLabel;
    bool   fVerbose;

    //TTree
    TTree* fTree;

    //Service handlesacktracker
    art::ServiceHandle<cheat::BackTrackerService> backtracker;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    int numevents;

    // tree variables
    int fRun;    ///< art run number
    int fSubRun; ///< art subrun number
    int fEvent;  ///< art event number

    vector<double> fE;
    vector<double> fP;
    vector<double> fPx;
    vector<double> fPy;
    vector<double> fPz;
    int fMode;
    int fCCNC;
    int fNuPdg;
    vector<int> fPdg;
    int fNParticles;
};
///////////////////////////////////////////////////////////////////////////////

ana::NCReco::NCReco(const fhicl::ParameterSet& pset) : 
  EDAnalyzer(pset),
  numevents(0)
{

  fGenieGenModuleLabel = pset.get<string>("GenieGenModuleLabel","generator");
  fLArGeantModuleLabel = pset.get<string>("LArGeantModuleLabel","largeant");
  fVerbose             = pset.get<bool>("Verbose",false); 

} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::NCReco::beginJob() {
  
  fTree = tfs->make<TTree>("ncrecoTree", "Tree with NC Reco information");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");
 
  fTree->Branch("E", &fE);
  fTree->Branch("P", &fP);
  fTree->Branch("Px", &fPx);
  fTree->Branch("Py", &fPy);
  fTree->Branch("Pz", &fPz);
  fTree->Branch("Mode", &fMode, "Mode/I");
  fTree->Branch("CCNC", &fCCNC, "CCNC/I");
  fTree->Branch("NuPdg", &fNuPdg, "NuPdg/I");  //for double-checking only. can get rid of later
  fTree->Branch("Pdg", &fPdg);
  fTree->Branch("NParticles", &fNParticles, "NParticles/I");
}

////////////////////////////////////////////////////////////////////
void ana::NCReco::analyze(const art::Event& evt) {

  fRun    = evt.run();
  fSubRun = evt.subRun();
  fEvent  = evt.event();

  if(fVerbose) {
      cout << "Analysing run " << fRun << ", subrun " << fSubRun << ", event: " << fEvent << endl;
  }
  numevents++;

  //Getting  MC truth information
  art::Handle< vector<simb::MCTruth> > mctruthListHandle;
  vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){
      art::fill_ptr_vector(mclist, mctruthListHandle);
  }

  //###############################################
  //### Get the Truth information for the event ###
  //###############################################

  //List the /particles in the event
//  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  map<int,const simb::MCParticle*> trueParticles;
  map<int,bool> mcparticlescontained;
  map<int,float> trueParticleEnergy;

 
  for(auto const& truth: mclist) {
    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    fCCNC=neutrino.CCNC();
    fMode=neutrino.Mode();
  
//  for (const auto& particleIt: particles) {
//
//    const simb::MCParticle* particle = particleIt.second;
//    if (particle->Process()!="primary"){
//      continue;
//    }
//    trueParticles[particle->TrackId()] = particle; //not sure what this line does. Didn't reproduce it w/MCptcl loop above

    int fNParticlesTot = truth->NParticles();
    for ( int i=0; i<fNParticlesTot; i++) {
      auto const& particle = truth->GetParticle(i);
      // get truth info from mcp  
      //think I need to put in all the stuff below this (from the if id=12 or 14 to before filling trees)
      //and also include all the stuff that was in the commented out lines below

      if (particle.Process()!="primary"){
        continue;
      }

      fNParticles++;

      int id = particle.PdgCode();
          
      if(id==12 || id==14){ //if NC neutrino, get the true E
        fNuPdg=id;
      }
      
      fE.push_back(particle.E());
      fP.push_back(particle.P());
      fPx.push_back(particle.Px());
      fPy.push_back(particle.Py());
      fPz.push_back(particle.Pz());
      fPdg.push_back(particle.PdgCode());
      
      if(fVerbose){
        cout << "True Particle with track ID: " << particle.TrackId() << " Has code of: " 
        << particle.PdgCode() << " and Energy of: " << particle.E() << " With Mother: " 
        << particle.Mother() << " Proccess: " << particle.Process() << " End Process: "  
        << particle.EndProcess() 
        << endl;
      }// if verbose
    }
  } // for MCParticles

  fTree->Fill();
  fE.clear();
  fP.clear();
  fPx.clear();
  fPy.clear();
  fPz.clear();
  fPdg.clear();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::NCReco::endJob() {

  cout << "Number of events processed: " <<  numevents  << endl;

}// end endJob


DEFINE_ART_MODULE(ana::NCReco)
