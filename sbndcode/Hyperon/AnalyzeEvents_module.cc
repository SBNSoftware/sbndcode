////////////////////////////////////////////////////////////////////////
// Class:       hyperon:AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        hyperon:AnalyzeEvents_module.cc
//
// Generated at Fri Oct  4 15:30:11 2024 by Jarek Nowak using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////



#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"



#include "art_root_io/TFileService.h"
// LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// // Root Includes
#include <iostream>
#include <vector>
#include <TTree.h>

namespace hyperon {
  class AnalyzeEvents;
}

class hyperon::AnalyzeEvents : public art::EDAnalyzer {
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

  int fVerbose;

  std::string fHitLabel, fGenieGenModuleLabel;
  std::vector<std::string> fPFParticleLabels;

  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

// Neutrino Interaction variables
   int intType, CCNC, neutrinoPDG, numProtons, numNeutrons, numPi, numPi0, numTrueHits;
   float W, X, Y, QSqr, Pt, Theta, neutrinoE, leptonP;
   float trueVertexX, trueVertexY, trueVertexZ;


  TTree *fTree;
  unsigned int fEventID;
};


hyperon::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  // ,
//  , fHitLabel(pset.get<std::string>("HitLabel"))
  , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::AnalyzeEvents::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
 fEventID = evt.id().event(); 
 std::cout<<"Event# "<<evt.id().event()<<std::endl;


// * MC truth information
   art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
   std::vector<art::Ptr<simb::MCTruth> > mclist;
   if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

   art::FindManyP< simb::MCParticle > fmpart( mctruthListHandle, evt, "largeant" );

   for (size_t i_truth = 0; i_truth < mclist.size(); i_truth++)
   {
       art::Ptr<simb::MCTruth> truth(mclist.at(i_truth));
       std::cout<<"Neutrino at index " << i_truth << " has pdg: " << truth->GetNeutrino().Nu().PdgCode() << std::endl;


       std::vector< art::Ptr<simb::MCParticle> > assocParticles(fmpart.at(i_truth));
       for (size_t i_mcpart = 0; i_mcpart < assocParticles.size(); i_mcpart++)
       {
           art::Ptr<simb::MCParticle> mcParticle(assocParticles.at(i_mcpart));
           std::cout<<"--Particle at index " << i_mcpart << " has pdg: " << mcParticle->PdgCode() << std::endl;
       }

   }

//  int nGeniePrimaries = 0, nGEANTparticles = 0, nMCNeutrinos = 0;




/*
 const std::vector<art::Ptr<simb::MCTruth>> truthVec = particleInventory->MCTruthVector_Ps();

  std::cout << std::setprecision(1) << std::fixed;
  if (fVerbose) {
    for (auto const &truth : truthVec) {
      std::cout << "Truth: " << truth << std::endl;
      if (truth->NeutrinoSet()) {
        const simb::MCNeutrino neutrino = truth->GetNeutrino();
        std::cout << "Neutrino: " << neutrino << std::endl;

        const simb::MCParticle nu = neutrino.Nu();
        std::cout << "X: " << nu.Vx() << " Y: " << nu.Vy() << " Z " << nu.Vz() << std::endl;
      } // truth->NeutrinoSet
    } // fVerbose
  } // truth: truthVec
  std::cout << std::setprecision(2) << std::fixed;
*/
/// std::cout<<nuSliceKey<<std::endl;
 fTree->Fill();
}


void hyperon::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree> ("tree", "Output TTree");

  //add branches here
  fTree ->Branch("eventID", &fEventID);
}

void hyperon::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(hyperon::AnalyzeEvents)
