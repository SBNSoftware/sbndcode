////////////////////////////////////////////////////////////////////////
// Class:       AnalyseEvents
// Plugin Type: analyzer (art v3_02_06)
// File:        AnalyseEvents_module.cc
//
// Generated at Sun Sep 29 15:12:48 2019 by Rhiannon Jones using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

// Core framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h" // Find associations as pointers
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// ROOT includes
#include <TTree.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>

namespace sbnd {
  class AnalyseEvents;
}


class sbnd::AnalyseEvents : public art::EDAnalyzer {
public:
  explicit AnalyseEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseEvents(AnalyseEvents const&) = delete;
  AnalyseEvents(AnalyseEvents&&) = delete;
  AnalyseEvents& operator=(AnalyseEvents const&) = delete;
  AnalyseEvents& operator=(AnalyseEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;

  // Output histogram declaration
  TH1D  *fTrackLengthHist;

  // Variables to access and analyse
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;
  std::vector<float> *fDaughterTrackLengths;

  // Module labels 
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;

};


sbnd::AnalyseEvents::AnalyseEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTree(nullptr),
    fTrackLengthHist(nullptr),
    fEventID(99999),
    fNPFParticles(99999),
    fNPrimaries(99999),
    fNPrimaryDaughters(-99999),
    fDaughterTrackLengths(nullptr)
    // Initialising all variables to obviously unphysical values 
    // so that we know when something just hasn't been filled
{
  // Get the producer name for all reco objects from the configuration fhicl file 
  //    analysisConfig.fcl
  fPFParticleLabel  = p.get<std::string>("PFParticleLabel");
  fTrackLabel       = p.get<std::string>("TrackLabel");
  fCalorimetryLabel = p.get<std::string>("CalorimetryLabel");
}  

void sbnd::AnalyseEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Define our event ID variable
  fEventID = e.id().event();
  
  // Initialise our counters for this event
  fNPFParticles         = 0;
  fNPrimaries           = 0;
  fDaughterTrackLengths->clear(); // Make sure the vector is empty at the beginning of the event

  // Access the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;
  if(e.getByLabel(fPFParticleLabel, pfparticleHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect, pfparticleHandle); // Fill the vector with art::Ptr PFParticles

  // Access the Tracks from pandoraTrack
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > trackVect;
  if(e.getByLabel(fTrackLabel, trackHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(trackVect, trackHandle); // Fill the vector with art::Ptr Tracks
  
  // Initiate neutrinoID to be non-physical so we can check if we've found it later
  size_t neutrinoID = 99999;

  // Let's find the ID of our neutrino so we can find it's daughters later on
  //
  //    NOTE: For this example, we will only search for a single neutrino in each event.
  //
  //    In a real, data-like analysis this would not be the case.
  //    Then we should really look for a list of neutrinos in each event and analyse them all.
  //
  //    However, we have only asked Pandora to reconstruct a single neutrino.
  //
  if(!pfparticleVect.size()) {
    std::cerr << " Error: No reconstructed particles found in the event. Skipping. " << std::endl;
    return;
  }
  
  // Find the neutrino ID
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    // Check if we are looking at a primary PFParticle & that it has a neutrino PdgCode
    // If we aren't move onto the next particle in the list
    if(!(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12))) continue;
    
    // We've found the neutrino! So let's define the neutrinoID variable
    neutrinoID         = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;

    std::cout << " Found the neutrino! PFParticle ID: " << neutrinoID << std::endl;
    break;
  }

  // Make sure we've found at the neutrino, if we haven't skip the event
  if(neutrinoID == 99999){
    std::cerr << " Error: No reconstructed neutrino found in the event. Skipping." << std::endl;
    return;
  }

  // Get the vector of vectors of tracks for each PFParticle
  // The vector size of associated tracks to a single PFParticle should be 0 or 1
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect, e, fTrackLabel);

  // Let's loop over our particles again and analyse the neutrino daughters
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    fNPFParticles++;

    // Check that the PFParticle is a daughter of the neutrino
    // If it isn't, move onto the next particle in the list
    if(pfp->Parent() != neutrinoID) continue;
      
    // Check whether there is an associated track to this PFParticle
    std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());

    if(pfpTracks.size()){
      for(art::Ptr<recob::Track> &trk : pfpTracks){
        // Get some information about the tracks
        // Fill the TTree variable
        fDaughterTrackLengths->push_back(trk->Length());

        // Fill the histogram
        fTrackLengthHist->Fill(trk->Length());
      }
    } // nTracks > 0
  } // pfparticleVect

  // Fill the output TTree with all relevant variables
  fTree->Fill();
}

void sbnd::AnalyseEvents::beginJob()
{
  // Implementation of optional member function here.
  // The TFileService is used to define the TTree and writing it to the output file
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","Analyser Output Tree");

  // Add a histogram to hold our track lengths to the output file
  // Since we've seen that they are distributed in ~5cm bins between 0 & 340cm 
  // we should set these as the definitions in the histogram constructor
  fTrackLengthHist = tfs->make<TH1D>("trackLengthHist","Reconstructed track lengths",70,0,340);

  // Add branches to our TTree
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/I");
  fTree->Branch("daughterTrackLengths",&fDaughterTrackLengths);

}

void sbnd::AnalyseEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::AnalyseEvents)
