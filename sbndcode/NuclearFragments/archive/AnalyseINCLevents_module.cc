////////////////////////////////////////////////////////////////////////
// Class:       AnalyseINCLevents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyseINCLevents_module.cc
//
// Generated at Thu Nov  7 11:04:12 2024 by Anna Beever using cetskelgen
// from cetlib version 3.18.02.
//
// MCTruth heavily based on Lan Nguyen's code <3
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
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"


//Additional framework includes
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//ROOT includes
#include <TTree.h>

//C++ includes
#include <iostream>

//Additional LArSoft includes
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace hepMCincl {
  class AnalyseINCLevents;
}


class hepMCincl::AnalyseINCLevents : public art::EDAnalyzer {
public:
  explicit AnalyseINCLevents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseINCLevents(AnalyseINCLevents const&) = delete;
  AnalyseINCLevents(AnalyseINCLevents&&) = delete;
  AnalyseINCLevents& operator=(AnalyseINCLevents const&) = delete;
  AnalyseINCLevents& operator=(AnalyseINCLevents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Output TTree
  TTree *fTree;

  //Tree variables
  unsigned int fEventID;
  int fNParticles;
  int fPDG;
  unsigned int fNPFParticles;
  unsigned int fNPrimaryChildren;

  //Define input labels
  std::string fSliceLabel;
  art::InputTag fPFParticleLabel;

};


hepMCincl::AnalyseINCLevents::AnalyseINCLevents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel", "pandora"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hepMCincl::AnalyseINCLevents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Set the event ID
  fEventID = e.id().event();
  std::cout << "===============================" << std::endl;
  std::cout << "The event ID is: " << fEventID << std::endl;

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  for (size_t i=0; i<ev_mct->size(); i++){
    
    simb::MCTruth const & mct = ev_mct->at(i);

    fNParticles = mct.NParticles();

    for(int j=0; j<mct.NParticles();j++){

      simb::MCParticle const & mcp = mct.GetParticle(j);
      
      if(mcp.TrackId() != j) {
        std::cout << "ERROR: \nTrackId does not match index\n" << std::endl;
      }
      
      if(abs(mcp.Vx())>210 ||  abs(mcp.Vy())>210||mcp.Vz()>510 || mcp.Vz()<-1){
        std::cout<<"OUTSIDE TPC x y z ="<<mcp.Vx()<<" "<<mcp.Vy()<<" "<<mcp.Vz()<<std::endl;
      }

      fPDG = mcp.PdgCode();
      std::cout << mcp.PdgCode() << std::endl;
      std::cout<<"vx vy vz ="<<mcp.Vx()<<" "<<mcp.Vy()<<" "<<mcp.Vz()<<std::endl;
      std::cout<<"px py pz ="<<mcp.Px()<<" "<<mcp.Py()<<" "<<mcp.Pz()<<std::endl;

      std::cout << mcp.Mother() << std::endl;

    }


  }

  /*fNPFParticles = 0;
  fNPrimaryChildren = 0;

  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if(sliceHandle.isValid()){
    art::fill_ptr_vector(sliceVector, sliceHandle);
  }

  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

  int nuID = -1;
  //int nuSliceKey = -1;

  for(const art::Ptr<recob::Slice> &slice: sliceVector){

    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

    // organize the PFPlist into a map
    std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
    for (unsigned i = 0; i < slicePFPs.size(); i++) {
      id_to_pfp[slicePFPs[i]->Self()] = slicePFPs[i];
    }

    for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs){

      const bool isPrimary(slicePFP->IsPrimary());
      const bool isNeutrino((std::abs(slicePFP->PdgCode()) == 12) || (std::abs(slicePFP->PdgCode()) == 14));

      if(!(isPrimary && isNeutrino)){
        continue;
      }

      //nuSliceKey = slice.key();
      nuID = slicePFP->Self();
      std::cout << "PDG: " << slicePFP->PdgCode() << std::endl;
      fNPFParticles = slicePFPs.size();
      std::cout << "NParticles: " << fNPFParticles << std::endl;
      fNPrimaryChildren = slicePFP->NumDaughters();
      std::cout << "NPrimaryChildren: " << fNPrimaryChildren << std::endl;
      for(unsigned i : slicePFP->Daughters()){
        std::cout << id_to_pfp.at(i)->PdgCode() << std::endl;
      }

      //std::cout here with particles and children

      break;
    }

    if(nuID >= 0){
      break;
    }

  }

  // get the PFParticle's and the associated data
  art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
  e.getByLabel(fPFParticleLabel, pfparticle_handle);

  std::vector<art::Ptr<recob::PFParticle>> pfparticles;
  art::fill_ptr_vector(pfparticles, pfparticle_handle);

  // organize the PFPlist into a map
  std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
  for (unsigned i = 0; i < pfparticles.size(); i++) {
    id_to_pfp[pfparticles[i]->Self()] = pfparticles[i];
    std::cout << pfparticles[i]->PdgCode() << std::endl;
  }*/

  const std::vector<simb::MCParticle> &mcparticle_list = *e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  for (unsigned i = 0; i < mcparticle_list.size(); i++){
    //if(mcparticle_list[i].Mother() == 1){
      std::cout << "PDG: " << mcparticle_list[i].PdgCode() << " | Self: " << mcparticle_list[i].TrackId() << " | Mother: " << mcparticle_list[i].Mother() << " | Ndaughters: " << mcparticle_list[i].NumberDaughters() << std::endl;
    //}
  }

  // Fill tree
  fTree->Fill();
}

void hepMCincl::AnalyseINCLevents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "");

  //Add branches to TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("NParticles", &fNParticles);
  fTree->Branch("pdg", &fPDG);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaryChildren", &fNPrimaryChildren);

}

void hepMCincl::AnalyseINCLevents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(hepMCincl::AnalyseINCLevents)
