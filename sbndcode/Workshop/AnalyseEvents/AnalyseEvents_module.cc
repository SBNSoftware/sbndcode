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

  // Variables to access and analyse
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;
  float fDaughterTrackLength;

  // Module labels 
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;

};


sbnd::AnalyseEvents::AnalyseEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  // Get the producer name for the PFParticles from the configuration fhicl file 
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
  fNPFParticles = 0;
  fNPrimaries   = 0;

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

  // Get the vector of vectors of tracks for each PFParticle
  // The vector size of associated tracks to a single PFParticle should be 0 or 1
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect, e, fTrackLabel);
  
  // Define a vector of IDs for the location of the neutrino daughters
  std::vector<size_t> daughterIDs;
  daughterIDs.clear();

  // Now, let's loop over our particles!
  if(!pfparticleVect.size()) return;
  for(art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    fNPFParticles++;
    if(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12)){ // Check we are looking at a primary PFParticle
      fNPrimaries++;
      fNPrimaryDaughters = pfp->NumDaughters();

      // Get the list of PFParticle IDs which correspond to neutrino daughters
      daughterIDs = pfp->Daughters();

      // Loop over the daughters and break out once finished
      for(unsigned int i = 0; i < daughterIDs.size(); ++i){
        art::Ptr<recob::PFParticle> daughter = pfparticleVect.at(daughterIDs[i]);

        // Check whether there is an associated track to this PFParticle
        std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(daughter->Self());

        if(pfpTracks.size()){
          for(art::Ptr<recob::Track> &trk : pfpTracks){
            // Get some information about the tracks
            fDaughterTrackLength = trk->Length();
          }
        } // nTracks > 0
      } // Daughters
    } // IsPrimary
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

  // Add branches to our TTree
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/I");
  fTree->Branch("daughterTrackLength",&fDaughterTrackLength,"daughterTrackLength/F");

}

void sbnd::AnalyseEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::AnalyseEvents)
