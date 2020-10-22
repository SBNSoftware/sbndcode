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
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileService.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "larsim/Utils/TruthMatchUtils.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"

// ROOT includes
#include <TTree.h>

// Other includes
#include <string>
#include <vector>


namespace sbnd {
  class ValidateTracks;
}


class sbnd::ValidateTracks : public art::EDAnalyzer {
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
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;

  // Variables to fill the output tree with
  unsigned int fEventID;

  // Truth
  std::vector<int>    *fTruePDG;
  std::vector<double> *fPosition;
  std::vector<double> *fPositionT;
  std::vector<double> *fEndPosition;
  std::vector<double> *fEndPositionT;
  std::vector<double> *fMomentum;
  std::vector<double> *fMomentumE;
  std::vector<double> *fMomentumP;
  std::vector<double> *fMomentumPt;
  std::vector<double> *fMomentumMass;
  std::vector<double> *fEndMomentum;
  std::vector<double> *fEndMomentumE;
  std::vector<double> *fEndMomentumP;
  std::vector<double> *fEndMomentumPt;
  std::vector<double> *fEndMomentumMass;

  // Reco
  unsigned int fNPFParticles;
  std::vector<int> *fNDaughthers;
  std::vector<int> *fParticleID;
  std::vector<int> *fParticlePDG;
  std::vector<bool> *fIsPrimary;
  std::vector<float> *fLengths;
  std::vector<float> *fVPoints;
  std::vector<float> *fStartDirPhi;
  std::vector<float> *fStartDirZ;
  std::vector<float> *fStartDirMag;
  std::vector<float> *fStartX;
  std::vector<float> *fStartY;
  std::vector<float> *fStartZ;
  std::vector<float> *fEndX;
  std::vector<float> *fEndY;
  std::vector<float> *fEndZ;

  // Module labels
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fG4Label;
  std::string fGenLabel;

  // Additional member functions
};


sbnd::ValidateTracks::ValidateTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // All vectors must be initialized in the class constructer
  // True
  fTruePDG(nullptr),
  fPosition(nullptr),
  fPositionT(nullptr),
  fEndPosition(nullptr),
  fEndPositionT(nullptr),
  fMomentum(nullptr),
  fMomentumE(nullptr),
  fMomentumP(nullptr),
  fMomentumPt(nullptr),
  fMomentumMass(nullptr),
  fEndMomentum(nullptr),
  fEndMomentumE(nullptr),
  fEndMomentumP(nullptr),
  fEndMomentumPt(nullptr),
  fEndMomentumMass(nullptr),
  // Reco
  fNDaughthers(nullptr),
  fParticleID(nullptr),
  fParticlePDG(nullptr),
  fIsPrimary(nullptr),
  fLengths(nullptr),
  fVPoints(nullptr),
  fStartDirPhi(nullptr),
  fStartDirZ(nullptr),
  fStartDirMag(nullptr),
  fStartX(nullptr),
  fStartY(nullptr),
  fStartZ(nullptr),
  fEndX(nullptr),
  fEndY(nullptr),
  fEndZ(nullptr)
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fG4Label = p.get<std::string>("G4Label");
  fGenLabel = p.get<std::string>("GenLabel");
}

void sbnd::ValidateTracks::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  // Define out event ID variable
  fEventID = evt.id().event();
  std::cout << "\n====================";
  std::cout << "Event ID: " << fEventID;
  /*
  // Make sure the vectors are empty and counters set to zero.
  // Truth
  fTruePDG->clear();
  fPosition->clear();
  fPositionT->clear();
  fEndPosition->clear();
  fEndPositionT->clear();
  fMomentum->clear();
  fMomentumE->clear();
  fMomentumP->clear();
  fMomentumPt->clear();
  fMomentumMass->clear();
  fEndMomentum->clear();
  fEndMomentumE->clear();
  fEndMomentumP->clear();
  fEndMomentumPt->clear();
  fEndMomentumMass->clear();
 */ 
  fTruePDG->clear();
  fPositionT->clear();
  fMomentumE->clear();

  // Reco
  fNPFParticles = 0;

  fNDaughthers  ->clear();
  fParticleID   ->clear();
  fParticlePDG  ->clear();
  fIsPrimary    ->clear();
  fLengths      ->clear();
  fVPoints      ->clear();
  fStartDirPhi  ->clear();
  fStartDirZ    ->clear();
  fStartDirMag  ->clear();
  fStartX       ->clear();
  fStartY       ->clear();
  fStartZ       ->clear();
  fEndX         ->clear();
  fEndY         ->clear();
  fEndZ         ->clear();

  // =========================================================
  // Truth stuff
  // Handles
  art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
  art::Handle<std::vector<simb::MCParticle>> mcpartHandle;
  // Object vectors
  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  std::vector<art::Ptr<simb::MCParticle>> mcparts;

  if(evt.getByLabel(fGenLabel, mctruthHandle)){
    art::fill_ptr_vector(mctruths, mctruthHandle);
  }
  if(evt.getByLabel(fG4Label, mcpartHandle)){
    art::fill_ptr_vector(mcparts, mcpartHandle);
  }

  for(const art::Ptr<simb::MCParticle> &part : mcparts){
    fTruePDG->push_back(part->PdgCode());
    fPositionT->push_back(part->Position().T());
    fMomentumE->push_back(part->Momentum().E());
  }

  // =========================================================
  // Reco
  // Accessing the PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfps;
  if(evt.getByLabel(fPFParticleLabel, pfpHandle)){ // Make sure the handle is valid
    art::fill_ptr_vector(pfps, pfpHandle); // Fill the vector with the art::Ptr PFParticles
  }

  // Access the Tracks from pandoraTrack
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > tracks;
  if(evt.getByLabel(fTrackLabel, trackHandle)){ // Make sure the handle is valid
    art::fill_ptr_vector(tracks, trackHandle); // Fill the vector with art::Ptr Tracks
  }

  if(!pfps.size()){
    std::cerr << "\nSkip event: No PFParticles found.";
    return; // Skip event if there are no reconstructed particles
  }

  std::cout << "PFParticles: " << pfps.size();
  fNPFParticles = pfps.size();

  // Get the vector or vectors of tracks for each PFParticle
  // The vector size of associated tracks to a single PFParticle should be 0 or 1: why?
  art::FindManyP<recob::Track> trackAssn(pfps, evt, fTrackLabel);

  // Find the neutrino ID
  for(const art::Ptr<recob::PFParticle> &pfp : pfps){

    fParticleID->push_back(pfp->Self());
    fParticlePDG->push_back(pfp->PdgCode());
    fNDaughthers->push_back(pfp->NumDaughters());

    if(pfp->IsPrimary()){
    	fIsPrimary->push_back(true);
    }
    else{
    	fIsPrimary->push_back(false);
    }

    // Check if there is an associated track to this particle
    std::vector<art::Ptr<recob::Track>> tracks = trackAssn.at(pfp.key());
    if(!tracks.empty()){
      assert(tracks.size()==1);
      std::cout << "Found a track!";
      fLengths->push_back(tracks.at(0)->Length());
      fVPoints->push_back(tracks.at(0)->CountValidPoints());        
      fStartDirPhi->push_back(tracks.at(0)->StartDirection().Phi());
      fStartDirMag->push_back(tracks.at(0)->StartDirection().Mag2());
      fStartDirZ->push_back(tracks.at(0)->StartDirection().Z());
      fStartX->push_back(tracks.at(0)->Start().X());
      fStartY->push_back(tracks.at(0)->Start().Y());
      fStartZ->push_back(tracks.at(0)->Start().Z());
      fEndX->push_back(tracks.at(0)->End().X());
      fEndY->push_back(tracks.at(0)->End().Y());
      fEndZ->push_back(tracks.at(0)->End().Z());
    } // if there is a track
    else{
      assert(tracks.size()==0);
      fLengths->emplace_back();
      fVPoints->emplace_back();        
      fStartDirPhi->emplace_back();
      fStartDirMag->emplace_back();
      fStartDirZ->emplace_back();
      fStartX->emplace_back();
      fStartY->emplace_back();
      fStartZ->emplace_back();
      fEndX->emplace_back();
      fEndY->emplace_back();
      fEndZ->emplace_back();      
    } // no track

  } // pfps vector

/*
    // Check if there is an associated track to this particle
    std::vector< art::Ptr<recob::Track> > this_tracks = trackAssn.at(pfp.key());
    if(!this_tracks.empty()){
      assert(this_tracks.size()==1);
      std::cout << "PFParticle has a track!" << std::endl;

      for(const art::Ptr<recob::Track> &track : this_tracks){

        fLengths->push_back(track->Length());
        fVPoints->push_back(track->CountValidPoints());        

        fStartDirPhi->push_back(track->StartDirection().Phi());
        fStartDirMag->push_back(track->StartDirection().Mag2());
        fStartDirZ->push_back(track->StartDirection().Z());

        fStartX->push_back(track->Start().X());
        fStartY->push_back(track->Start().Y());
        fStartZ->push_back(track->Start().Z());

        fEndX->push_back(track->End().X());
        fEndY->push_back(track->End().Y());
        fEndZ->push_back(track->End().Z());

      }
    } // nTracks > 0
  } // pfps vector
*/
  // Fill the output tree with all the relevant variables
  fTree->Fill();
}

void sbnd::ValidateTracks::beginJob()
{
  // Implementation of optional member function here.
  // The TFileService is used to define the TTree and writing it to the output file
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analysis Output Tree");

  //Add branches to out tree
  fTree->Branch("eventID",      &fEventID, "eventID/i");
  // Truth
  fTree->Branch("truePDG",    &fTruePDG);
  fTree->Branch("positionT",  &fPositionT);
  fTree->Branch("momentumE",  &fMomentumE);
  // Reco
  fTree->Branch("nPFParticles", &fNPFParticles, "nPFParticles/i");
  fTree->Branch("isPrimary",    &fIsPrimary);
  fTree->Branch("nDaughters",   &fNDaughthers);
  fTree->Branch("particleID",   &fParticleID);
  fTree->Branch("particlePDG",  &fParticlePDG);
  fTree->Branch("Lengths",      &fLengths);
  fTree->Branch("ValidPoints",  &fVPoints);
  fTree->Branch("StartDirPhi", 	&fStartDirPhi);
  fTree->Branch("StartDirZ",    &fStartDirZ);
  fTree->Branch("StartDirMag",  &fStartDirMag);
  fTree->Branch("StartX",       &fStartX);
  fTree->Branch("StartY",       &fStartY);
  fTree->Branch("StartZ",       &fStartZ);
  fTree->Branch("EndX",         &fEndX);
  fTree->Branch("EndY",         &fEndY);
  fTree->Branch("EndZ",         &fEndZ);

}

void sbnd::ValidateTracks::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::ValidateTracks)
