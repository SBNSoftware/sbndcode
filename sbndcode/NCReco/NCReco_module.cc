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

    vector<double> fETrue;
    vector<int> fMode;
    double fEh;
    double fPh;
    double fPzh;
//    double fEReco; 
    bool fProtonOnly;
    int fCCNC;
    int fNProtons = 0;
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
 
  fTree->Branch("ETrue", &fETrue);
  fTree->Branch("Mode", &fMode);
  fTree->Branch("EReco", &fEReco);

  fTree->Branch("ProtonOnly", &fProtonOnly);  //for double-checking only. can get rid of later
  fTree->Branch("CCNC", &fCCNC);
  fTree->Branch("NProtons", &fNProtons);
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

  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  map<int,const simb::MCParticle*> trueParticles;
  map<int,bool> mcparticlescontained;
  map<int,float> trueParticleEnergy;

  fProtonOnly = true;

  for(auto const& truth: mclist) {
    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    fCCNC=neutrino.CCNC();
    fMode=neutrino.Mode();
  }

  if(fCCNC==1){ 
  //Loops over particles, Makes a map of Track id and pdgcode
    for (const auto& particleIt: particles) {

      const simb::MCParticle* particle = particleIt.second;
      trueParticles[particle->TrackId()] = particle;
      int id = particle->TrackId();
      double Energy = particle.E();

      if(id==14){ //if NC neutrino, get the true E
	fETrue.push_back(Energy);

      }else if(id==2212 && Energy>0.025 && fNProtons==0){ //If proton, check if energy is above threshold, and if so, if this is the only one so far.
	fNProtons = 1;
	fEh = particle.E();
	fPh = particle.P();
	fPzh = particle.Pz();

      } else if(id==2212 && Energy>0.025 && fNProtons>0){
	fNProtons = 2;
	break

      } else if((abs(id)==211 && Energy>0.01) || (id=22 && Energy> 0.03)) { //check if there are pions or photons. If so, skip event
	fProtonOnly = false;
	break
      }

      if(fVerbose){
          cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " 
               << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " 
               << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  
               << particle->EndProcess() 
          << endl;
      }// if verbose

    } // for MCParticles
  }

  if(fCCNC==1 && fNProtons==1 && fProtonOnly==true){
    fTree.Fill();
  }  

  fETrue.clear();
  fMode.clear();
////  fEReco.clear();
  fNProtons.clear();
  fProtonOnly.clear();
  fCCNC.clear();
  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::NCReco::endJob() {

  cout << "Number of events processed: " <<  numevents  << endl;

}// end endJob


DEFINE_ART_MODULE(ana::NCReco)
