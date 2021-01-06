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

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
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

  int NumberOfSharedHits(int TrackID, detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits); //Returns the number of hits in the vector that are matched to a MC track.

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;
  // Services
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  art::ServiceHandle<geo::Geometry> geom;
  // Variables to fill the output tree with
  unsigned int fEventID;
  std::vector<int> *mctruth_pdg;
  // Truth
  std::vector<bool> *mctruth_isnc;
  std::vector<bool> *mctruth_iscc;
  // MCParticles
  std::vector<double> *mcpart_trackid;
  std::vector<int>    *mcpart_pdg;
  std::vector<int>    *mcpart_ntpoints;
  std::vector<double> *mcpart_position_t;
  std::vector<double> *mcpart_momentum_e;
  std::vector<double> *mcpart_momentum_p;
  std::vector<double> *mcpart_momentum_pt;
  std::vector<double> *mcpart_momentum_mass;
  std::vector<double> *mcpart_vx;
  std::vector<double> *mcpart_vy;
  std::vector<double> *mcpart_vz;
  std::vector<double> *mcpart_endx;
  std::vector<double> *mcpart_endy;
  std::vector<double> *mcpart_endz;

  // Reco
  // --- pfparticles
  unsigned int pfpart_number;
  std::vector<int>  *pfpart_ndaughthers;
  std::vector<int>  *pfpart_id;
  std::vector<int>  *pfpart_pdg;
  std::vector<bool> *pfpart_isprimary;
  // // -- clusters
  std::vector<int>   *cluster_id;
  std::vector<int>   *cluster_nhits;
  std::vector<float> *cluster_width;
  std::vector<int>   *cluster_plane;
  // --- tracks
  std::vector<int>   *track_nspoints;
  std::vector<int>   *track_nclusters;
  std::vector<float> *track_lenght;
  std::vector<float> *track_vpoints;
  std::vector<float> *track_startdirphi;
  std::vector<float> *track_startdirz;
  std::vector<float> *track_startdirmag;
  std::vector<float> *track_startx;
  std::vector<float> *track_starty;
  std::vector<float> *track_startz;
  std::vector<float> *track_endx;
  std::vector<float> *track_endy;
  std::vector<float> *track_endz;
  std::vector<int>   *track_nhits;
  // --- hits
  std::vector< std::vector<raw::TDCtick_t>>   *hits_starttick;
  std::vector< std::vector<raw::TDCtick_t>>   *hits_endtick;
  std::vector< std::vector<raw::ChannelID_t>> *hits_channel;
  std::vector< std::vector<geo::View_t>>      *hits_view;
  std::vector< std::vector<geo::WireID>>      *hits_wireid;
  std::vector< std::vector<int>>              *hits_planeid;
  std::vector< std::vector<float>>            *hits_adc;
  std::vector<float> *hits_rms;
  // --- cals should be vector of vector
  std::vector<std::vector<float>> *track_dqdx;
  std::vector<std::vector<float>> *track_dedx;
  std::vector<std::vector<float>> *track_resrange;
  // --- chi2 should be vector of vector
  std::vector<double> *track_chi2kaon;
  std::vector<double> *track_chi2muon;
  std::vector<double> *track_chi2pion;
  std::vector<double> *track_chi2proton;
  std::vector<double> *track_chi2ndof;
  std::vector<double> *track_pida;

  // Module labels
  std::string fPFParticleLabel;
  std::string fSpacePointLabel;
  std::string fClusterLabel;
  std::string fTrackLabel;
  std::string fHitLabel;
  std::string fG4Label;
  std::string fGenLabel;
  std::string fCalLabel;
  std::string fPIDLabel;
  bool fVerbose;
  // Additional member functions
};


sbnd::ValidateTracks::ValidateTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // All vectors must be initialized in the class constructer
  // True
  mctruth_pdg(nullptr),
  mctruth_isnc(nullptr),
  mctruth_iscc(nullptr),
  // MC Particle
  mcpart_trackid(nullptr),
  mcpart_pdg(nullptr),
  mcpart_ntpoints(nullptr),
  mcpart_position_t(nullptr),
  mcpart_momentum_e(nullptr),
  mcpart_momentum_p(nullptr),
  mcpart_momentum_pt(nullptr),
  mcpart_momentum_mass(nullptr),
  mcpart_vx(nullptr),
  mcpart_vy(nullptr),
  mcpart_vz(nullptr),
  mcpart_endx(nullptr),
  mcpart_endy(nullptr),
  mcpart_endz(nullptr),
  // Reco
  pfpart_ndaughthers(nullptr),
  pfpart_id(nullptr),
  pfpart_pdg(nullptr),
  pfpart_isprimary(nullptr),
  // // --- clusters
  cluster_id(nullptr),
  cluster_nhits(nullptr),
  cluster_width(nullptr),
  cluster_plane(nullptr),
  // --- tracks
  track_nspoints(nullptr),
  track_nclusters(nullptr),
  track_lenght(nullptr),
  track_vpoints(nullptr),
  track_startdirphi(nullptr),
  track_startdirz(nullptr),
  track_startdirmag(nullptr),
  track_startx(nullptr),
  track_starty(nullptr),
  track_startz(nullptr),
  track_endx(nullptr),
  track_endy(nullptr),
  track_endz(nullptr),
  track_nhits(nullptr),
  // --- hits
  hits_starttick(nullptr),
  hits_endtick(nullptr),
  hits_channel(nullptr),
  hits_view(nullptr),
  hits_wireid(nullptr),
  hits_planeid(nullptr),
  hits_adc(nullptr),
  hits_rms(nullptr),
  // --- cals
  track_dqdx(nullptr),
  track_dedx(nullptr),
  track_resrange(nullptr),
  // --- chi2
  track_chi2kaon(nullptr),
  track_chi2muon(nullptr),
  track_chi2pion(nullptr),
  track_chi2proton(nullptr),
  track_chi2ndof(nullptr),
  track_pida(nullptr)

  {
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fSpacePointLabel = p.get<std::string>("SpacePointLabel");
  fClusterLabel = p.get<std::string>("ClusterLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fHitLabel = p.get<std::string>("HitLabel");
  fG4Label = p.get<std::string>("G4Label");
  fGenLabel = p.get<std::string>("GenLabel");
  fCalLabel = p.get<std::string>("CalLabel");
  fPIDLabel = p.get<std::string>("PIDLabel");

  fVerbose = p.get<bool>("Verbose");
}

void sbnd::ValidateTracks::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  // Define out event ID variable
  fEventID = evt.id().event();
  std::cout << "\n\n====================\n";
  std::cout << "Event ID: " << fEventID << std::endl;

  // Make sure the vectors are empty and counters set to zero.
  // Truth
  mctruth_pdg->clear();
  mctruth_isnc->clear();
  mctruth_iscc->clear();
  // Matched MC Particles
  mcpart_trackid          ->clear();
  mcpart_pdg              ->clear();
  mcpart_ntpoints         ->clear();
  mcpart_position_t       ->clear();
  mcpart_momentum_e       ->clear();
  mcpart_momentum_p       ->clear();
  mcpart_momentum_pt      ->clear();
  mcpart_momentum_mass    ->clear();
  mcpart_vx               ->clear();
  mcpart_vy               ->clear();
  mcpart_vz               ->clear();
  mcpart_endx             ->clear();
  mcpart_endy             ->clear();
  mcpart_endz             ->clear();
  // Reco
  pfpart_number = 0;
  pfpart_ndaughthers ->clear();
  pfpart_id          ->clear();
  pfpart_pdg         ->clear();
  pfpart_isprimary   ->clear();
  cluster_id    ->clear();
  cluster_nhits ->clear();
  cluster_width ->clear();
  cluster_plane ->clear();

  track_lenght       ->clear();
  track_vpoints      ->clear();
  track_startdirphi  ->clear();
  track_startdirz    ->clear();
  track_startdirmag  ->clear();
  track_startx       ->clear();
  track_starty       ->clear();
  track_startz       ->clear();
  track_endx         ->clear();
  track_endy         ->clear();
  track_endz         ->clear();
  track_nspoints     ->clear();
  track_nclusters    ->clear();
  track_nhits        ->clear();

  hits_starttick ->clear();
  hits_endtick   ->clear();
  hits_planeid   ->clear();
  hits_adc       ->clear();

  track_dqdx       ->clear();
  track_dedx       ->clear();
  track_resrange   ->clear();
  track_chi2kaon   ->clear();
  track_chi2muon   ->clear();
  track_chi2pion   ->clear();
  track_chi2proton ->clear();
  track_chi2ndof   ->clear();
  track_pida       ->clear();

  // =========================================================
  // Services
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // =========================================================
  // Truth stuff
  // Handles
  art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
  art::Handle<std::vector<simb::MCParticle>> mcpartHandle;
  // Object vectors
  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  std::vector<art::Ptr<simb::MCParticle>> mcparts;

  if(fVerbose){std::cout << "declared mctruth and mcpart hangles and vectors" << std::endl;}

  if(evt.getByLabel(fGenLabel, mctruthHandle)){
    art::fill_ptr_vector(mctruths, mctruthHandle);
  }
  if(evt.getByLabel(fG4Label, mcpartHandle)){
    art::fill_ptr_vector(mcparts, mcpartHandle);
  }

  if(fVerbose){std::cout << "mctruths size: " << mctruths.size() << ", mcparts size: " << mcparts.size() << std::endl;}
  //const simb::MCNeutrino& nu = mctruths[0]->GetNeutrino();
  for(const auto& mctruth: mctruths){
  	const simb::MCNeutrino nu = mctruth->GetNeutrino();
  	mctruth_pdg->push_back(nu.Nu().PdgCode());
  	mctruth_isnc->push_back(nu.CCNC());
  	mctruth_iscc->push_back(!nu.CCNC());
  }
  // //List the particles in the event
  // const sim::ParticleList& particles = particleInventory->ParticleList();
  // for(const auto& particle: particles) {
  //   const simb::MCParticle* this_particle = particle.second;
  // }

  // =========================================================
  // Reco
  // Handles
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  art::Handle< std::vector<recob::SpacePoint> > spointHandle;
  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  art::Handle< std::vector<recob::Track> > trackHandle;
  art::Handle< std::vector<recob::Hit> > hitHandle;
  // Object vectors
  std::vector< art::Ptr<recob::PFParticle> > pfps;
  std::vector< art::Ptr<recob::SpacePoint> > spoints;
  std::vector< art::Ptr<recob::Cluster> > clusters;
  std::vector< art::Ptr<recob::Track> > tracks;
  std::vector< art::Ptr<recob::Hit> > hits;

  if(evt.getByLabel(fPFParticleLabel, pfpHandle)){ // Make sure the handle is valid
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "pfpHandle not valid" << std::endl;}
        return;
    }
    art::fill_ptr_vector(pfps, pfpHandle); // Fill the vector with the art::Ptr PFParticles
  }
  if(evt.getByLabel(fSpacePointLabel, spointHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "spointHandle not valid" << std::endl;}
        return;
    }
    art::fill_ptr_vector(spoints, spointHandle);
  }
  if(evt.getByLabel(fClusterLabel, clusterHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "clusterHandle not valid" << std::endl;}
        return;
    }
  	art::fill_ptr_vector(clusters, clusterHandle);
  }
  if(evt.getByLabel(fTrackLabel, trackHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "trackHandle not valid" << std::endl;}
        return;
    }
    art::fill_ptr_vector(tracks, trackHandle);
  }
  if(evt.getByLabel(fHitLabel, hitHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "hitHandle not valid" << std::endl;}
        return;
    }
    art::fill_ptr_vector(hits, hitHandle);
  }

  if(!pfps.size()){
    if(fVerbose){std::cout << "\nSkip event: No PFParticles found.";}
    return; // Skip event if there are no reconstructed particles
  }

  if(fVerbose){std::cout << "\nPFParticles: " << pfps.size() << std::endl;}
  pfpart_number = pfps.size();

  // Get the vector of tracks for each PFParticle via association
  // The vector size of associated tracks to a single PFParticle
  // should be 0 or 1 as a PFP will be either a shower or a track
  art::FindManyP<recob::SpacePoint> spointAssn(pfps, evt, fSpacePointLabel);
  art::FindManyP<recob::Cluster>    clusterAssn(pfps, evt, fPFParticleLabel);
  art::FindManyP<recob::Hit>        clusterhitAssn(clusters, evt, fHitLabel);
  art::FindManyP<recob::Track>      trackAssn(pfps, evt, fTrackLabel);
  art::FindManyP<recob::Hit>        hitAssn(tracks, evt, fTrackLabel);
  art::FindManyP<anab::Calorimetry> calAssn(tracks, evt, fCalLabel);
  art::FindManyP<anab::ParticleID>  pidAssn(tracks, evt, fPIDLabel);
  // Find the neutrino ID

  for(const art::Ptr<recob::PFParticle> &pfp : pfps){

    pfpart_id->push_back(pfp->Self());
    pfpart_pdg->push_back(pfp->PdgCode());
    pfpart_ndaughthers->push_back(pfp->NumDaughters());

    if(pfp->IsPrimary()){
    	pfpart_isprimary->push_back(true);
    }
    else{
    	pfpart_isprimary->push_back(false);
    }

    // Check associated space points
    const std::vector<art::Ptr<recob::SpacePoint>> spoints = spointAssn.at(pfp.key());
    if(spoints.empty()){
    	track_nspoints->push_back(0);
    	if(fVerbose){std::cout << "Empty space point" << std::endl;}
    }
    else{
    	track_nspoints->push_back(spoints.size());
    	if(fVerbose){std::cout << "Found space points. spoints.size(): " << spoints.size() << std::endl;}
    }

    // Get the hits from the track by accessing the clusters
    std::vector<art::Ptr<recob::Hit> > cluster_hits;
    const std::vector< art::Ptr< recob::Cluster> > thisclusters = clusterAssn.at(pfp.key());
    if(fVerbose){std::cout << "thisclusters.size(): " << thisclusters.size() << std::endl;}
    for (const auto& cluster: thisclusters){
    	cluster_id->push_back(cluster->ID());
    	cluster_nhits->push_back(cluster->NHits());
    	cluster_width->push_back(cluster->Width());
    	if(cluster->hasPlane()){
    		cluster_plane->push_back(cluster->Plane().Plane);
    	}
    }
    track_nclusters->push_back(thisclusters.size());


    // Check if there is an associated track to this particle    
    const std::vector<art::Ptr<recob::Track>> tracks = trackAssn.at(pfp.key());
    if(tracks.empty()){
    	if(fVerbose){std::cout << "No tracks." << std::endl;}
    	continue;
      /*
      assert(tracks.size()==0);
      track_lenght->emplace_back();
      */
    }
    else{
      assert(tracks.size()==1);
      if(fVerbose){std::cout << "Found a track!" << std::endl;}
      track_lenght->push_back(tracks.at(0)->Length());
      track_vpoints->push_back(tracks.at(0)->CountValidPoints());        
      track_startdirphi->push_back(tracks.at(0)->StartDirection().Phi());
      track_startdirmag->push_back(tracks.at(0)->StartDirection().Mag2());
      track_startdirz->push_back(tracks.at(0)->StartDirection().Z());
      track_startx->push_back(tracks.at(0)->Start().X());
      track_starty->push_back(tracks.at(0)->Start().Y());
      track_startz->push_back(tracks.at(0)->Start().Z());
      track_endx->push_back(tracks.at(0)->End().X());
      track_endy->push_back(tracks.at(0)->End().Y());
      track_endz->push_back(tracks.at(0)->End().Z());

      // Get the hits associated to the track
      const std::vector<art::Ptr<recob::Hit>> trackhits = hitAssn.at(tracks[0].key());
      track_nhits->push_back(trackhits.size());
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
      } // hits
      hits_starttick ->push_back(temp_start);
      hits_endtick   ->push_back(temp_end);
      hits_planeid   ->push_back(temp_planeid);
      hits_adc       ->push_back(temp_adc);  

      // Use the hits to get the G4ID
      auto mc_matched_id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clock_data, trackhits, true); // rollupUnsavedIDs
      bool valid_match = TruthMatchUtils::Valid(mc_matched_id);
      if (!valid_match){
      	std::cout << "Unable to find MCParticle matched to this track." << std::endl;
      	continue;
      }
      // Use particle inventory service to get the mcparticle contributing to this track
      // For the truth matching, get the MCParticle which trackid is the same as the g4id
      const simb::MCParticle *particle = particleInventory->TrackIdToParticle_P(mc_matched_id);
      mcpart_trackid->push_back(particle->TrackId());
      mcpart_pdg->push_back(particle->PdgCode());
      mcpart_ntpoints->push_back(int(particle->NumberTrajectoryPoints()));
      mcpart_position_t->push_back(particle->Position().T());
      mcpart_momentum_e->push_back(particle->Momentum().E());
      mcpart_momentum_p->push_back(particle->P());
      mcpart_momentum_pt->push_back(particle->Pt());
      mcpart_momentum_mass->push_back(particle->Mass());
      mcpart_vx->push_back(particle->Vx());
      mcpart_vy->push_back(particle->Vy());
      mcpart_vz->push_back(particle->Vz());
      mcpart_endx->push_back(particle->EndY());
      mcpart_endy->push_back(particle->EndY());
      mcpart_endz->push_back(particle->EndZ());

      int nsharedhits = NumberOfSharedHits(mc_matched_id, clock_data, trackhits);
      // Returns the number of hits in the vector that are associated to the MC track.
      int ntruehits = nsharedhits;
      int nrecohits = trackhits.size();
      float hitcompleteness = -9999.;
      float hitpurity = -9999.;
      if(nrecohits!=0)
      	hitpurity = float(nsharedhits)/float(nrecohits);
      if(ntruehits!=0)
      	hitcompleteness = float(nsharedhits)/float(ntruehits);
      std::cout << "shared hits: " << nsharedhits << ", recohits: " << nrecohits << ", truehits: " << ntruehits << std::endl;
      std::cout << "hitpurity: " << hitpurity << ", hitcompleteness: " << hitcompleteness << std::endl;

      // Get the calorimetry information associated to this track
      std::vector< art::Ptr<anab::Calorimetry> > trackcals = calAssn.at(tracks[0].key());
      if(trackcals.empty()) continue;
      for(auto &cal : trackcals){
        if(!cal->PlaneID().isValid) continue;
        // Only look at the collection plane, since this is where the dEdx is acquired
        int calplaneid = cal->PlaneID().Plane;
        if(calplaneid!=2) continue;
        track_dqdx->push_back(cal->dQdx());
        track_dedx->push_back(cal->dEdx());
        track_resrange->push_back(cal->ResidualRange());
        // assert(calplaneid==2);
        // track_dqdx->push_back(cal->dQdx().at(2));
        // track_dedx->push_back(cal->dEdx().at(2));
        // track_resrange->push_back(cal->ResidualRange().at(2));
      } // calorimetry

      // Get the chi2 information associated to this track
      std::vector< art::Ptr<anab::ParticleID> > trackpids = pidAssn.at(tracks[0].key());
      if(trackpids.empty()) continue;
      for(auto &pid : trackpids){
        if(pid->PlaneID()){
          assert(pid->PlaneID().deepestIndex()<3);
          track_chi2kaon->push_back(pid->Chi2Kaon());
          track_chi2muon->push_back(pid->Chi2Muon());
          track_chi2pion->push_back(pid->Chi2Pion());
          track_chi2proton->push_back(pid->Chi2Proton());
          track_chi2ndof->push_back(pid->Ndf());
          track_pida->push_back(pid->PIDA());
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
  fTree->Branch("mctruth_pdg", 	        &mctruth_pdg);
  fTree->Branch("mctruth_isnc",         &mctruth_isnc);
  fTree->Branch("mctruth_iscc",         &mctruth_iscc);
  // Matched MCParticle
  fTree->Branch("mcpart_trackid",       &mcpart_trackid);
  fTree->Branch("mcpart_pdg",           &mcpart_pdg);
  fTree->Branch("mcpart_ntpoints",      &mcpart_ntpoints);
  fTree->Branch("mcpart_position_t",    &mcpart_position_t);
  fTree->Branch("mcpart_momentum_e",    &mcpart_momentum_e);
  fTree->Branch("mcpart_momentum_p",    &mcpart_momentum_p);
  fTree->Branch("mcpart_momentum_pt",   &mcpart_momentum_pt);  
  fTree->Branch("mcpart_momentum_mass", &mcpart_momentum_mass);  
  fTree->Branch("mcpart_vx",        &mcpart_vx);
  fTree->Branch("mcpart_vy",        &mcpart_vy);
  fTree->Branch("mcpart_vz",        &mcpart_vz);
  fTree->Branch("mcpart_endx",      &mcpart_endx);
  fTree->Branch("mcpart_endy",      &mcpart_endy);
  fTree->Branch("mcpart_endz",      &mcpart_endz);
  // Reco
  fTree->Branch("pfpart_number",      &pfpart_number, "nPFParticles/i");
  fTree->Branch("pfpart_isprimary",   &pfpart_isprimary);
  fTree->Branch("pfpart_ndaughthers", &pfpart_ndaughthers);
  fTree->Branch("pfpart_id",          &pfpart_id);
  fTree->Branch("pfpart_pdg",         &pfpart_pdg);
  fTree->Branch("cluster_id",        &cluster_id);
  fTree->Branch("cluster_nhits",     &cluster_nhits);
  fTree->Branch("cluster_width",     &cluster_width);
  fTree->Branch("cluster_plane",     &cluster_plane);
  fTree->Branch("track_nspoints",    &track_nspoints);
  fTree->Branch("track_nclusters",   &track_nclusters);
  fTree->Branch("track_lenght",      &track_lenght);
  fTree->Branch("track_vpoints",     &track_vpoints);
  fTree->Branch("track_startdirphi", &track_startdirphi);
  fTree->Branch("track_startdirz",   &track_startdirz);
  fTree->Branch("track_startdirmag", &track_startdirmag);
  fTree->Branch("track_startx",      &track_startx);
  fTree->Branch("track_starty",      &track_starty);
  fTree->Branch("track_startz",      &track_startz);
  fTree->Branch("track_endx",        &track_endx);
  fTree->Branch("track_endy",        &track_endy);
  fTree->Branch("track_endz",        &track_endz);
  fTree->Branch("track_nhits",       &track_nhits);
  fTree->Branch("hits_starttick",    &hits_starttick);
  fTree->Branch("hits_endtick",      &hits_endtick);
  fTree->Branch("hits_planeid",      &hits_planeid);
  fTree->Branch("hits_adc",          &hits_adc);
  fTree->Branch("track_dqdx",        &track_dqdx);  
  fTree->Branch("track_dedx",        &track_dedx);
  fTree->Branch("track_resrange",    &track_resrange);
  fTree->Branch("track_chi2kaon",    &track_chi2kaon);
  fTree->Branch("track_chi2muon",    &track_chi2muon);
  fTree->Branch("track_chi2pion",    &track_chi2pion);
  fTree->Branch("track_chi2proton",  &track_chi2proton);
  fTree->Branch("track_chi2ndof",    &track_chi2ndof);
  fTree->Branch("track_pida",        &track_pida);

}

int sbnd::ValidateTracks::NumberOfSharedHits(int trackID, detinfo::DetectorClocksData const& clock_data, const std::vector<art::Ptr<recob::Hit>> &hits){

  int nsharedhits = 0;
  for(auto hit : hits){
    auto hit_matched_id = TruthMatchUtils::TrueParticleID(clock_data, hit, true);
    if(hit_matched_id == trackID)
      nsharedhits++;
  }

  return nsharedhits;
}


void sbnd::ValidateTracks::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::ValidateTracks)
