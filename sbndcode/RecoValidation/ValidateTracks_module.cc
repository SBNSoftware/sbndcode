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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
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
  // --- tracks
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
  // --- pfparticles
  unsigned int fNPFParticles;
  std::vector<int> *fNDaughthers;
  std::vector<int> *fParticleID;
  std::vector<int> *fParticlePDG;
  std::vector<bool> *fIsPrimary;
  // --- tracks
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
  // --- hits
  std::vector<unsigned int> *fNHits;
  // std::vector<raw::TDCtick_t> *fStartTick;
  // std::vector<raw::TDCtick_t> *fEndTick;
  // std::vector<raw::ChannelID_t> *fHitChannel;
  // std::vector<geo::View_t> *fHitView;
  // std::vector<geo::WireID> *fHitWireID;
  // std::vector<int> *fHitPlaneID;
  // std::vector<float> *fHitADC;
  // std::vector<float> *fHitRMS;
  std::vector< std::vector<raw::TDCtick_t>> *fStartTick;
  std::vector< std::vector<raw::TDCtick_t>> *fEndTick;
  std::vector<raw::ChannelID_t> *fHitChannel;
  std::vector<geo::View_t> *fHitView;
  std::vector<geo::WireID> *fHitWireID;
  std::vector< std::vector<int>> *fHitPlaneID;
  std::vector< std::vector<float>> *fHitADC;
  std::vector<float> *fHitRMS;
  // --- cals should be vector of vector
  std::vector<std::vector<float>> *fTrackdQdx;
  std::vector<std::vector<float>> *fTrackdEdx;
  std::vector<std::vector<float>> *fTrackResidualRange;
  // --- chi2 should be vector of vector
  std::vector<double> *fChi2Kaon;
  std::vector<double> *fChi2Muon;
  std::vector<double> *fChi2Pion;
  std::vector<double> *fChi2Proton;
  std::vector<double> *fChi2NDoF;
  std::vector<double> *fPIDA;


  // Module labels
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fHitLabel;
  std::string fG4Label;
  std::string fGenLabel;
  std::string fCalLabel;
  std::string fPIDLabel;

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
  // --- tracks
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
  fEndZ(nullptr),
  // --- hits
  fNHits(nullptr),
  fStartTick(nullptr),
  fEndTick(nullptr),
  fHitChannel(nullptr),
  fHitView(nullptr),
  fHitWireID(nullptr),
  fHitPlaneID(nullptr),
  fHitADC(nullptr),
  fHitRMS(nullptr),
  // --- cals
  fTrackdQdx(nullptr),
  fTrackdEdx(nullptr),
  fTrackResidualRange(nullptr),
  // --- chi2
  fChi2Kaon(nullptr),
  fChi2Muon(nullptr),
  fChi2Pion(nullptr),
  fChi2Proton(nullptr),
  fChi2NDoF(nullptr),
  fPIDA(nullptr)

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fHitLabel = p.get<std::string>("HitLabel");
  fG4Label = p.get<std::string>("G4Label");
  fGenLabel = p.get<std::string>("GenLabel");
  fCalLabel = p.get<std::string>("CalLabel");
  fPIDLabel = p.get<std::string>("PIDLabel");
}

void sbnd::ValidateTracks::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  // Define out event ID variable
  fEventID = evt.id().event();
  std::cout << "\n\n====================\n";
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
  fNHits        ->clear();
  fStartTick    ->clear();
  fEndTick      ->clear();
  fHitPlaneID   ->clear();
  fHitADC       ->clear();
  fTrackdQdx          ->clear();
  fTrackdEdx          ->clear();
  fTrackResidualRange ->clear();
  fChi2Kaon   ->clear();
  fChi2Muon   ->clear();
  fChi2Pion   ->clear();
  fChi2Proton ->clear();
  fChi2NDoF   ->clear();
  fPIDA       ->clear();
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
  // Handles
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  art::Handle< std::vector<recob::Track> > trackHandle;
  art::Handle< std::vector<recob::Hit> > hitHandle;
  // Object vectors
  std::vector< art::Ptr<recob::PFParticle> > pfps;
  std::vector< art::Ptr<recob::Track> > tracks;
  std::vector< art::Ptr<recob::Hit> > hits;

  if(evt.getByLabel(fPFParticleLabel, pfpHandle)){ // Make sure the handle is valid
    art::fill_ptr_vector(pfps, pfpHandle); // Fill the vector with the art::Ptr PFParticles
  }
  if(evt.getByLabel(fTrackLabel, trackHandle)){
    art::fill_ptr_vector(tracks, trackHandle);
  }
  if(evt.getByLabel(fHitLabel, hitHandle)){
    art::fill_ptr_vector(hits, hitHandle);
  }

  if(!pfps.size()){
    std::cerr << "\nSkip event: No PFParticles found.";
    return; // Skip event if there are no reconstructed particles
  }

  std::cout << "\nPFParticles: " << pfps.size() << std::endl;
  fNPFParticles = pfps.size();

  // Get the vector of tracks for each PFParticle via association
  // The vector size of associated tracks to a single PFParticle
  // should be 0 or 1 as a PFP will be either a shower or a track
  art::FindManyP<recob::Track> trackAssn(pfps, evt, fTrackLabel);
  art::FindManyP<recob::Hit> hitAssn(tracks, evt, fTrackLabel);
  art::FindManyP<anab::Calorimetry> calAssn(tracks, evt, fCalLabel);
  art::FindManyP<anab::ParticleID> pidAssn(tracks, evt, fPIDLabel);
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
    if(tracks.empty()){
      continue;
      /*
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
      */
    }
    else{
      assert(tracks.size()==1);
      std::cout << "Found a track!" << std::endl;
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

      // Get the hits associated to the track
      std::vector<art::Ptr<recob::Hit>> trackhits = hitAssn.at(tracks[0].key());
      fNHits->push_back(trackhits.size());
      if(trackhits.empty()) continue;
      std::vector<raw::TDCtick_t> temp_start;
      std::vector<raw::TDCtick_t> temp_end;
      std::vector<int> temp_planeid;
      std::vector<float> temp_adc;
      for (auto const& hit: trackhits){
        temp_start.push_back(hit->StartTick());
        temp_end.push_back(hit->EndTick());
        temp_planeid.push_back(hit->WireID().Plane);
        temp_adc.push_back(hit->SummedADC());
      /*
        fStartTick->push_back(hit->StartTick());
        fEndTick->push_back(hit->EndTick());
        fHitPlaneID->push_back(hit->WireID().Plane);
        fHitADC->push_back(hit->SummedADC());
        // fHitChannel->push_back(hit->Channel());
        // fHitView->push_back(hit->View());
        // fHitWireId->push_back(hit->WireID());
        // fHitRMS->push_back(hit->RMS());
        */
      } // hits
      fStartTick  ->push_back(temp_start);
      fEndTick    ->push_back(temp_end);
      fHitPlaneID ->push_back(temp_planeid);
      fHitADC     ->push_back(temp_adc);  

      // Get the calorimetry information associated to this track
      std::vector< art::Ptr<anab::Calorimetry> > trackcals = calAssn.at(tracks[0].key());
      if(trackcals.empty()) continue;
      for(auto &cal : trackcals){
        if(!cal->PlaneID().isValid) continue;
        // Only look at the collection plane, since this is where the dEdx is acquired
        int calplaneid = cal->PlaneID().Plane;
        if(calplaneid!=2) continue;
        fTrackdQdx->push_back(cal->dQdx());
        fTrackdEdx->push_back(cal->dEdx());
        fTrackResidualRange->push_back(cal->ResidualRange());
      } // calorimetry

      // Get the chi2 information associated to this track
      std::vector< art::Ptr<anab::ParticleID> > trackpids = pidAssn.at(tracks[0].key());
      if(trackpids.empty()) continue;
      for(auto &pid : trackpids){
        if(pid->PlaneID()){
          assert(pid->PlaneID().deepestIndex()<3);
          fChi2Kaon->push_back(pid->Chi2Kaon());
          fChi2Muon->push_back(pid->Chi2Muon());
          fChi2Pion->push_back(pid->Chi2Pion());
          fChi2Proton->push_back(pid->Chi2Proton());
          fChi2NDoF->push_back(pid->Ndf());
          fPIDA->push_back(pid->PIDA());
        }

      }

    } // tracks
  } // pfps

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
  fTree->Branch("nHits",        &fNHits);
  fTree->Branch("StartTick",    &fStartTick);
  fTree->Branch("EndTick",      &fEndTick);
  fTree->Branch("HitPlaneID",   &fHitPlaneID);
  fTree->Branch("fHitADC",      &fHitADC);
  fTree->Branch("fTrackdQdx",           &fTrackdQdx);  
  fTree->Branch("fTrackdEdx",           &fTrackdEdx);
  fTree->Branch("fTrackResidualRange",  &fTrackResidualRange);
  fTree->Branch("fChi2Kaon", &fChi2Kaon);
  fTree->Branch("fChi2Muon", &fChi2Muon);
  fTree->Branch("fChi2Pion", &fChi2Pion);
  fTree->Branch("fChi2Proton", &fChi2Proton);
  fTree->Branch("fChi2NDoF", &fChi2NDoF);
  fTree->Branch("fPIDA", &fPIDA);

}

void sbnd::ValidateTracks::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::ValidateTracks)
