////////////////////////////////////////////////////////////////////////
// Class:       ValidateTracks
// Plugin Type: analyzer (art v3_05_01)
// File:        ValidateTracks_module.cc
//
// Generated at Sat Oct  3 17:53:53 2020 by Diana Mendez mendez using cetskelgen
// from cetlib version v3_10_00.
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

// Additional framework includes
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/PFParticle.h"

// ROOT includes
#include <TTree.h>

// Other includes
#include <string>
#include <vector>


namespace testValidation {
  class ValidateTracks;
}


class testValidation::ValidateTracks : public art::EDAnalyzer {
public:
  explicit ValidateTracks(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ValidateTracks(ValidateTracks const&) = delete;
  ValidateTracks(ValidateTracks&&) = delete;
  ValidateTracks& operator=(ValidateTracks const&) = delete;
  ValidateTracks& operator=(ValidateTracks&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;

  // Variable to fill the output tree with
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNDaughthers;
  std::string fPFParticleLabel;

};


testValidation::ValidateTracks::ValidateTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
	fPFParticleLabel = p.get<std::string>("PFParticleLabel");
}

void testValidation::ValidateTracks::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Define out event ID variable
  fEventID = e.id().event();

  // Initialize the counters for this event
  fNPFParticles = 0;
  fNPrimaries 	= 0;
  fNDaughthers = 0;

  // Accessing the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;

  if(e.getByLabel(fPFParticleLabel, pfparticleHandle)) // Make sure the handle is valid
  	art::fill_ptr_vector(pfparticleVect, pfparticleHandle); // Fill the vector with the art::Ptr PFParticles

  if(!pfparticleVect.size()) return; // Skip event if there are no reconstructed particles
  fNPFParticles = pfparticleVect.size();

  size_t neutrinoID = 99999; // Initiate neutrino ID to be non-physical for checks

  // Find the neutrino ID
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
  	// Check that we are looking at a primary particle and that it has a neutrino pdg code
  	// Move on to the next pfparticle in the list if we arent
  	if(! (pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12))) continue;
  	neutrinoID = pfp->Self();
  	fNDaughthers = pfp->NumDaughters();
  	fNPrimaries++;
  }

  if(neutrinoID == 99999) return; // Skip the event if no neutrino was found

  // Fill the output tree with all the relevant variables
  fTree->Fill();
}

void testValidation::ValidateTracks::beginJob()
{
  // Implementation of optional member function here.
  // The TFileService is used to define the TTree and writing it to the output file
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analysis Output Tree");

  //Add branches to out tree
  fTree->Branch("eventID", &fEventID, "eventID/i");
  fTree->Branch("nPFParticles", &fNPFParticles, "nPFParticles/i");
  fTree->Branch("nPrimaries", &fNPrimaries, "nPrimaries/i");
  fTree->Branch("nDaughters", &fNDaughthers, "nDaughters/i");
}

void testValidation::ValidateTracks::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(testValidation::ValidateTracks)
