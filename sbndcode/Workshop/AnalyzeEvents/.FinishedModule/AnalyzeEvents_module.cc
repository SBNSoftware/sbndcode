////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (art v3_05_01)
// File:        AnalyzeEvents_module.cc
//
// Generated at Mon Nov  2 14:19:17 2020 by Owen Goodwin <owen.goodwin@postgrad.manchester.ac.uk> using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////


//default includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//Addtional frame work include
#include "canvas/Persistency/Common/FindManyP.h" // Find associations as pointers
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardataobj/RecoBase/Hit.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

//Root include
#include <TTree.h>
#include <TH1D.h>

//c++ include
#include <vector>
#include <string>

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

  // Output tree declaration
  TTree *fTree;


  //Out put histogram declaration 
  TH1D *fTrackLengthHist;


  // Variables to go in output tree
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;


  std::vector<bool>                 *fDaughterLongestTrack;


  std::vector<int>                  *fDaughterTrackTruePDG;
  
  std::vector<float> *fDaughterTrackLengths;
  std::vector<std::vector<float> > *fDaughterTrackdEdx;
  std::vector<std::vector<float> > *fDaughterTrackResidualRange;

  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;


  

};


test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  //Get the producer name for the PFParticles from the config fhicl file [analysisConfig.fcl]
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fCalorimetryLabel = p.get<std::string>("CalorimetryLabel");
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{

  //we need the particle inventory service
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  //and the clock data for event
 auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // Implementation of required member function here.
  // Define our event ID variable
  fEventID = e.id().event();


  //check is MC

  bool isMC = !e.isRealData();

  // Intiallise our counters for each event
  fNPFParticles = 0;
  fNPrimaries = 0;

  // Get the ID of longest track to add to out tree
  // Initialise longest track variable
  float longest   = -99999.;
  int   longestID = -1;
 

  fDaughterLongestTrack->clear();
  fDaughterTrackLengths->clear();  //otherwise it will keep particles from previous event
  fDaughterTrackdEdx->clear();
  fDaughterTrackResidualRange->clear();

  fDaughterTrackTruePDG->clear();

  // Access the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;
  if(e.getByLabel(fPFParticleLabel, pfparticleHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect, pfparticleHandle); // Fill the vector with art::Ptr PFParticles

  if(pfparticleVect.empty()) return; //If no reco'd particles skip event
  fNPFParticles = pfparticleVect.size();
  // Initiate neutrinoID to be non-physical so we can check if we've found it later



    // Access the Tracks from pandoraTrack
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > trackVect;
  if(e.getByLabel(fTrackLabel, trackHandle)) // Make sure the handle is valid
    art::fill_ptr_vector(trackVect, trackHandle); // Fill the vector with art::Ptr Tracks



  size_t neutrinoID = 99999;
  // Find the neutrino ID
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    // Check if we are looking at a primary PFParticle & that it has a neutrino PdgCode
    // If we aren't move onto the next particle in the list
    if(!(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12))) continue;

    // We've found the neutrino! So let's define the neutrinoID variable
    neutrinoID = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;
  }
  if(neutrinoID == 99999) return; //if we haven't found neutrino, skip event
  // Fill the output tree with all relevant varibles


   // Get the vector of vectors of tracks for each PFParticle & calorimetry for each track
  // The vector size of associated tracks to a single PFParticle should be 0 or 1
  art::FindManyP<recob::Track>      trackAssoc      (pfparticleVect, e, fTrackLabel);

      // The track-calo association vector size should be 0-3, 1 for each of the available planes
  art::FindManyP<anab::Calorimetry> calorimetryAssoc(trackVect,      e, fCalorimetryLabel);

  //get hits
  //association made by track produceder
  art::FindManyP<recob::Hit> hitAssoc(trackVect,      e, fTrackLabel);

  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
      if(pfp->Parent() != neutrinoID) continue;

    // Check that the PFParticle is a daughter of the neutrino
    // If it isn't, move onto the next particle in the list
     // Check whether there is an associated track to this PFParticle
    std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
    if(pfpTracks.empty()) continue;

  // Find the ID of the longest track and iterate over their lengths
  for(const art::Ptr<recob::Track> &trk : pfpTracks){
    if(trk->Length() > longest){
      longest   = trk->Length();
      longestID = trk->ID();
    }
  }

    //loop through tracks associated to the pfp [should only be one per pfp, (or zero)]
    for(const art::Ptr<recob::Track> &trk : pfpTracks){
      fDaughterTrackLengths->push_back(trk->Length());
          // Fill the histogram
      fTrackLengthHist->Fill(trk->Length());

          // Fill the longest track boolean
      if(trk->ID() == longestID) fDaughterLongestTrack->push_back(true);
      else fDaughterLongestTrack->push_back(false);

      std::vector< art::Ptr<anab::Calorimetry> > trackCalo = calorimetryAssoc.at(trk.key());


      if(trackCalo.empty()) continue;

      for(const art::Ptr<anab::Calorimetry> &cal : trackCalo){
        if(!cal->PlaneID().isValid) continue;

        // Get the plane number
        int planenum = cal->PlaneID().Plane;

        // Only look at the collection plane, since this is where the dEdx
        // is acquired
        if (planenum!=2) continue;

        // Fill the dEdx & residual range entries for this track
        fDaughterTrackdEdx->push_back(cal->dEdx());
        fDaughterTrackResidualRange->push_back(cal->ResidualRange());
      } // Calorimetry



      if(isMC){
         // Get the vector of hits for this track
        std::vector< art::Ptr<recob::Hit> > trackHits = hitAssoc.at(trk.key());
        // Use utilitiy to get the ID code of the true particle which has contributed most energy to the track
        int trkidtruth = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trackHits, true);
        // Use Particle Invertory Service to get the mc particle contributing to this ID
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trkidtruth);
      // Save the PDG code of the track
        fDaughterTrackTruePDG->push_back(particle->PdgCode());
      }


    } //end track loop
  } //end pfp loop



  
  fTree->Fill();

}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  // The TFileService is used to define the TTree and writing it to the output file

  //define service
  art::ServiceHandle<art::TFileService> tfs;


  //define tree
  fTree = tfs->make<TTree>("tree","Analyser Output Tree");



  // Add a histogram to hold our track lengths to the output file
  // Since we've seen that they are distributed in ~5cm bins between 0 & 340cm 
  // we should set these as the definitions in the histogram constructor
  fTrackLengthHist        = tfs->make<TH1D>("trackLengthHist",         "Reconstructed track lengths",70,0,340);


  // Add branches to our TTree
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/i");
  fTree->Branch("daughterTrackLengths",&fDaughterTrackLengths);
  fTree->Branch("daughterLongestTrack",&fDaughterLongestTrack);
  fTree->Branch("daughterTrackdEdx",&fDaughterTrackdEdx);
  fTree->Branch("daughterTrackResidualRange",&fDaughterTrackResidualRange);
  fTree->Branch("daughterTrackTruePDG",&fDaughterTrackTruePDG);
}

void test::AnalyzeEvents::endJob()
{
 


  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
