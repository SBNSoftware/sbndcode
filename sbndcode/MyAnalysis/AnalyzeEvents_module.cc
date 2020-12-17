////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (art v3_05_01)
// File:        AnalyzeEvents_module.cc
//
// Generated at Thu Dec  3 10:52:26 2020 by Vu Chi Lan Nguyen using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

// Default framework includes 
 
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
#include "canvas/Persistency/Common/FindManyP.h" //Find associations as pointers
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ROOT includes
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

// C++ includes 
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>


namespace test {
  class AnalyzeEvents;
}

//define class

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

  //Output histogram
  TH1D *fTrackLengthHist;
  TH2D *fdEdxvsResidualRange;  
  TH2D *fdQdxvsResidualRange; 

  // Declare tree
  TTree *fTree;

  // Variables in tree
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  unsigned int fNPrimaryDaughters;

  std::vector<float> 			*fDaughterTrackLengths;  
  std::vector< std::vector<float> > 	*fDaughterTrackdEdx;
  std::vector< std::vector<float> > 	*fDaughterTrackResidualRange;
  std::vector<bool> 			*fDaughterLongestTrack;
  std::vector<int> 			*fDaughterTrackTruePDG; 
  std::vector< std::vector<float> >	*fDaughterTrackdQdx;
  std::vector<float>			*fDaughterTrackHits;
  std::vector<float> 			*fX;
  std::vector<float>			*fY;
  std::vector<float>			*fZ;
//  std::vector< geo::Point_t >		*fPosition;
 
  //labels
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fSpacePointLabel;

  //Additional member functions
  size_t GetNeutrinoID(const std::vector< art::Ptr<recob::PFParticle> > &pfps, unsigned int &nPrimaries) const;
  
  std::vector< art::Ptr<recob::PFParticle> > GetNeutrinoDaughters(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const size_t &neutrinoID) const;
  
  std::vector< art::Ptr<recob::Track> > GetNeutrinoDaughterTracks(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const art::FindManyP<recob::Track> &trackAssoc) const;
  
  int GetLongestTrackID(const std::vector< art::Ptr<recob::Track> > &trks) const;

}; //End of class AnalyzeEvents


// class AnalyzeEvents constructor

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTrackLengthHist(nullptr),
    fdEdxvsResidualRange(nullptr),
    fdQdxvsResidualRange(nullptr),
    fTree(nullptr),
    fEventID(99999),
    fNPFParticles(99999),
    fNPrimaries(99999), 
    fNPrimaryDaughters(99999),
    fDaughterTrackLengths(nullptr),
    fDaughterTrackdEdx(nullptr),
    fDaughterTrackResidualRange(nullptr),
    fDaughterLongestTrack(nullptr),
    fDaughterTrackTruePDG(nullptr),
    fDaughterTrackdQdx(nullptr),
    fDaughterTrackHits(nullptr),
    fX(nullptr),
    fY(nullptr),
    fZ(nullptr)
  //Initialising all variables to unphysical values 
  //Can check when something hasn't been filled
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  // Get the producer name for the PFParticles from the config file  
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fCalorimetryLabel = p.get<std::string>("CalorimetryLabel");
  fSpacePointLabel = p.get<std::string>("SpacePointLabel");
} //End of constructor


// analyze() function to goes through each event
void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
 
  //services
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
   
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  // Define event ID
  fEventID = e.id().event();

  //clear
  fDaughterTrackLengths->clear();
  fDaughterTrackdEdx->clear();
  fDaughterTrackResidualRange->clear(); 
  fDaughterLongestTrack->clear();
  fDaughterTrackTruePDG->clear();
  fDaughterTrackdQdx->clear();
  fDaughterTrackHits->clear();
 
  // Counter
  fNPFParticles = 0;
  fNPrimaries = 0;
  fNPrimaryDaughters = 0;

  // Access PFParticles from Pandora
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;
  if(e.getByLabel(fPFParticleLabel, pfparticleHandle)){  	//Check if handle is valie 
    art::fill_ptr_vector(pfparticleVect, pfparticleHandle); 	//Fill vector with art::Ptr PFParticle
  } 

  // Access Tracks from pandoraTrack
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > trackVect;
  if(e.getByLabel(fTrackLabel, trackHandle)){
    art::fill_ptr_vector(trackVect, trackHandle);
  } 

  //If no reco events, skip
  if(pfparticleVect.empty()){
    std::cerr<<"Error: No reco'ed particles found in event. Skipping...\n";
    return;
  }

  // Get # of pfp
  fNPFParticles = pfparticleVect.size();

  // Get the neutrino ID using GetNeutrinoID()
  size_t neutrinoID = this->GetNeutrinoID(pfparticleVect, fNPrimaries);

  //If no primary neutrinos found, skip the event
  if(neutrinoID == 99999){
    std::cerr<<"Error: No reco'ed neutrino found in event. Skipping...\n";
    return;
  }
 
  //Get a list of neutrino daughters as a vector of PFParticle using GetNeutrinoDaughters()
  std::vector< art::Ptr<recob::PFParticle> > neutrinoDaughters = this->GetNeutrinoDaughters(pfparticleVect, neutrinoID);

  //If no daughter partickes, skip ths event
  if(neutrinoDaughters.empty()){
    std::cerr<<"Error: No reco'ed neutrino daughters found in event. Skipping...\n";
    return;
  }

  //Get # of neutrino daughters
  fNPrimaryDaughters = neutrinoDaughters.size();

  //Get vector of vector of tracks for each PFParticle & calorimetry for each track


  //-----Define association-----//
  //
  //Vector size of associated track to a single PFParticle is 0 or 1
  art::FindManyP<recob::Track> trackAssoc (pfparticleVect, e, fTrackLabel);

  //track-calo vector size is 0-3, 1 for each plane
  art::FindManyP<anab::Calorimetry> calorimetryAssoc(trackVect, e, fCalorimetryLabel);
    
  art::FindManyP<recob::Hit> hitAssoc(trackVect, e, fTrackLabel);

  art::FindManyP<recob::SpacePoint> spacepointAssoc(trackVect, e, fSpacePointLabel);

  //Get neutrino daughter trakcs using GetNeutrinoDaughterTracks
  std::vector< art::Ptr<recob::Track> > neutrinoDaughterTracks = this->GetNeutrinoDaughterTracks(neutrinoDaughters, trackAssoc);

  //If no tracks found, skip
  if(neutrinoDaughterTracks.empty()){
    std::cerr<<"Error: No tracks asso'ed to reco'ed neutrino daughters found in event. Skipping...\n";
    return;
  } 

  //Get ID of longest track. Later add this add a flag in the tree
  int longestID = this->GetLongestTrackID(neutrinoDaughterTracks);

  if(longestID == -1){
    std::cerr<<"Error: Invalid longest track ID. Skipping...\n";
    return;
  }

  //-------Loop over track----------//
  
  for(const art::Ptr<recob::Track> &trk: neutrinoDaughterTracks){
    
    fDaughterTrackLengths->push_back(trk->Length());

    if(trk->ID() == longestID) fDaughterLongestTrack->push_back(true);
    else fDaughterLongestTrack->push_back(false);

    //Fill histogram
    fTrackLengthHist->Fill(trk->Length());

    //Get association with tracks
    
    std::vector< art::Ptr<anab::Calorimetry> > trackCalo = calorimetryAssoc.at(trk.key());
    std::vector< art::Ptr<recob::Hit> > trackHits = hitAssoc.at(trk.key());
    std::vector< art::Ptr<recob::SpacePoint> > trackSpacePoint = spacepointAssoc.at(trk.key());

    if(trackCalo.empty()) continue;

    for(const art::Ptr<anab::Calorimetry> &cal: trackCalo){
      if(!cal->PlaneID().isValid) continue;

      //Get plane number. Collection plane = 2
      int planenum = cal->PlaneID().Plane;

      if(planenum!=2) continue;
   
      fDaughterTrackdEdx->push_back(cal->dEdx());
      fDaughterTrackResidualRange->push_back(cal->ResidualRange()); 
      fDaughterTrackdQdx->push_back(cal->dQdx()); 
 
      //Fill histogram
      for(unsigned int i = 0; i < cal->dEdx().size(); i++){
        fdEdxvsResidualRange->Fill(cal->ResidualRange()[i], cal->dEdx()[i]);
        fdQdxvsResidualRange->Fill(cal->ResidualRange()[i], cal->dQdx()[i]);
      }
          
    } //End of Calorimetry
 

    //Get track hit truth info

    //Get ID code of true particle contributing the most energy to the track
    int trkidtruth = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trackHits, true);
    //Use Particle Inventory Service to get the MC particle contributing to this ID
    const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trkidtruth);
    //Save the PDG code of the track
    fDaughterTrackTruePDG->push_back(particle->PdgCode());

    //Get associated hit of track
    for(const art::Ptr<recob::Hit> &hit: trackHits){
      if(!hit->WireID().isValid) continue;

      int planenum = hit->WireID().Plane;
      
      if(planenum!=2) continue;
     
      fDaughterTrackHits->push_back(hit->Integral());
    
      for(const art::Ptr<recob::SpacePoint> &sp: trackSpacePoint){
      //  int ID = sp->ID();
        fX->push_back(sp->XYZ()[0]);
	fY->push_back(sp->XYZ()[1]);
        fZ->push_back(sp->XYZ()[2]);
       // fPosition->push_back(sp->position());
      }

    }//End of Hits 


  } //End of Tracks

  //Fill tree
  fTree->Fill();
  
} //End of analyze()

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  
  // TFileService define and write tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","Analyser Output Tree");
  
  // Add branch
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/i");
  fTree->Branch("daughterTrackLengths",&fDaughterTrackLengths);
  fTree->Branch("daughterTrackdEdx",&fDaughterTrackdEdx);
  fTree->Branch("daughterTrackResidualRange",&fDaughterTrackResidualRange);
  fTree->Branch("daughterLongestTrack",&fDaughterLongestTrack);
  fTree->Branch("daughterTrackTruePDG",&fDaughterTrackTruePDG);
  fTree->Branch("daughterTrackdQdx", &fDaughterTrackdQdx);
  fTree->Branch("daughterTrackHits", &fDaughterTrackHits);
  fTree->Branch("X",&fX);
  fTree->Branch("Y",&fY);
  fTree->Branch("Z",&fZ);
 // fTree->Branch("Position",&fPosition); 
 
  // Histograms
  fTrackLengthHist = tfs->make<TH1D>("trackLengthHist", "Reconstructed track lengths",70,0,500 );
  fdEdxvsResidualRange = tfs->make<TH2D>("dEdxvsResidualRange", "dEdx vs Residual Range", 200, 0, 50, 200, 0, 30);
  fdQdxvsResidualRange = tfs->make<TH2D>("dQdxvsResidualRange", "dQdx vs Residual Range", 200, 0, 50, 200, 0, 6500);
} //End of beginJob()

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
} //End of endJob()


size_t test::AnalyzeEvents::GetNeutrinoID(const std::vector< art::Ptr<recob::PFParticle> > &pfps, unsigned int &nPrimaries) const{
  
  //NOTE: This example only search for a single neutrino in each event. For real data-like analysis, should look for a list of neutrino in each event

  //Initiate non-physical ID for checking later 
  size_t neutrinoID = 99999;

  //Find the NeutrinoID
  for(const art::Ptr<recob::PFParticle> &pfp: pfps){
    //check if primary particle if neutrino, else skip
    if(!(pfp->IsPrimary() &&(std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode())==12))) continue;
    nPrimaries++;
    neutrinoID = pfp->Self();
  }
  return neutrinoID;
}//End of GetNeutrinoID()

std::vector< art::Ptr<recob::PFParticle> > test::AnalyzeEvents::GetNeutrinoDaughters(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const size_t &neutrinoID) const{
  //vector to hold daughters and return them at the end
  std::vector< art::Ptr<recob::PFParticle> > daughters;
  daughters.clear();

  for(const art::Ptr<recob::PFParticle> &pfp: pfps){
    //Check if pfparticle is daughter of neutrino, else skip
    if(pfp->Parent() != neutrinoID) continue;
    daughters.push_back(pfp);
  }
  return daughters;
}//End of GetNeutrinoDaughters()

std::vector< art::Ptr<recob::Track> > test::AnalyzeEvents::GetNeutrinoDaughterTracks(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const art::FindManyP<recob::Track> &trackAssoc) const{
  //vector to hold the daughter tracks
  std::vector< art::Ptr<recob::Track> > daughterTracks;
  daughterTracks.clear();

  for(const art::Ptr<recob::PFParticle> &pfp: pfps){
    //get track
    std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
    if(pfpTracks.empty()) continue;
    
    for(const art::Ptr<recob::Track> &trk: pfpTracks){
      daughterTracks.push_back(trk);
    }  
  }
  return daughterTracks;
} //End of GetNeutrinoDaughterTracks

int test::AnalyzeEvents::GetLongestTrackID(const std::vector< art::Ptr<recob::Track> > &trks) const{

  float longest = -99999;
  int longestID = -1;

  //Iterate over track length and get ID of longest track
  for(const art::Ptr<recob::Track> &trk: trks){
    if(trk->Length() > longest){
      longest = trk->Length();
      longestID = trk->ID();
    }
  }
  return longestID;
}//End of GetLongestTrackID

DEFINE_ART_MODULE(test::AnalyzeEvents)
