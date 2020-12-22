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

private:

  // Declare member data here.
  // Output tree declaration
  TTree *fTree;
  // Services
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  // Variables to fill the output tree with
  unsigned int fEventID;
  int mctruth_pdg;
  bool mctruth_isnc;
  bool mctruth_iscc;

  // Truth
  // --- tracks
  // Truth matching by hand. Does not save all the MCParticles
  std::vector<int>    *mcpart_pdg_hm;
  std::vector<double> *mcpart_positiont_hm;
  std::vector<double> *mcpart_momentume_hm;
  std::vector<double> *mcpart_trackid_hm;
  std::vector<double> *mcpart_vx_hm;
  std::vector<double> *mcpart_vy_hm;
  std::vector<double> *mcpart_vz_hm;
  std::vector<double> *mcpart_endx_hm;
  std::vector<double> *mcpart_endy_hm;
  std::vector<double> *mcpart_endz_hm;
  // Truth matching with particle service
  std::vector<int>    *mcpart_pdg;
  std::vector<double> *mcpart_position;
  std::vector<double> *mcpart_positiont;
  std::vector<double> *mcpart_endposition;
  std::vector<double> *mcpart_endpositiont;
  std::vector<double> *mcpart_momentum;
  std::vector<double> *mcpart_momentum_e;
  std::vector<double> *mcpart_momentum_p;
  std::vector<double> *mcpart_momentum_pt;
  std::vector<double> *mcpart_momentum_mass;
  std::vector<double> *mcpart_endmomentum;
  std::vector<double> *mcpart_endmomentum_e;
  std::vector<double> *mcpart_endmomentum_p;
  std::vector<double> *mcpart_endmomentum_pt;
  std::vector<double> *mcpart_endmomentum_mass;
  std::vector<int>    *mcpart_trackid;
  std::vector<int>    *mcpart_g4id;

  // Reco
  // --- pfparticles
  unsigned int pfpart_number;
  std::vector<int>  *pfpart_ndaughthers;
  std::vector<int>  *pfpart_id;
  std::vector<int>  *pfpart_pdg;
  std::vector<bool> *pfpart_isprimary;
  // --- tracks
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
  std::vector<int>  *track_nspoints;
  std::vector<int>  *track_nclusters;
  std::vector<int>  *track_nhits;
  std::vector<int>  *cluster_nhits;
  // --- hits
  // std::vector<raw::TDCtick_t> *hits_starttick;
  // std::vector<raw::TDCtick_t> *hits_endtick;
  // std::vector<raw::ChannelID_t> *hits_channel;
  // std::vector<geo::View_t> *hits_view;
  // std::vector<geo::WireID> *hits_wireid;
  // std::vector<int> *hits_planeid;
  // std::vector<float> *hits_adc;
  // std::vector<float> *hits_rms;
  std::vector< std::vector<raw::TDCtick_t>> *hits_starttick;
  std::vector< std::vector<raw::TDCtick_t>> *hits_endtick;
  std::vector<raw::ChannelID_t> *hits_channel;
  std::vector<geo::View_t> *hits_view;
  std::vector<geo::WireID> *hits_wireid;
  std::vector< std::vector<int>> *hits_planeid;
  std::vector< std::vector<float>> *hits_adc;
  std::vector<float> *hits_rms;
  // --- cals should be vector of vector
  std::vector<std::vector<float>> *track_dqdx;
  std::vector<std::vector<float>> *track_dedx;
  std::vector<std::vector<float>> *track_resrange;
  // std::vector<float> *track_dqdx;
  // std::vector<float> *track_dedx;
  // std::vector<float> *track_resrange;
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
  // mctruth_pdg(nullptr),
  // mctruth_isnc(nullptr),
  // mctruth_iscc(nullptr),
  mcpart_pdg_hm(nullptr),
  mcpart_positiont_hm(nullptr),
  mcpart_momentume_hm(nullptr),
  mcpart_trackid_hm(nullptr),
  mcpart_vx_hm(nullptr),
  mcpart_vy_hm(nullptr),
  mcpart_vz_hm(nullptr),
  mcpart_endx_hm(nullptr),
  mcpart_endy_hm(nullptr),
  mcpart_endz_hm(nullptr),
  mcpart_pdg(nullptr),
  mcpart_position(nullptr),
  mcpart_positiont(nullptr),
  mcpart_endposition(nullptr),
  mcpart_endpositiont(nullptr),
  mcpart_momentum(nullptr),
  mcpart_momentum_e(nullptr),
  mcpart_momentum_p(nullptr),
  mcpart_momentum_pt(nullptr),
  mcpart_momentum_mass(nullptr),
  mcpart_endmomentum(nullptr),
  mcpart_endmomentum_e(nullptr),
  mcpart_endmomentum_p(nullptr),
  mcpart_endmomentum_pt(nullptr),
  mcpart_endmomentum_mass(nullptr),
  mcpart_trackid(nullptr),
  mcpart_g4id(nullptr),
  // Reco
  pfpart_ndaughthers(nullptr),
  pfpart_id(nullptr),
  pfpart_pdg(nullptr),
  pfpart_isprimary(nullptr),
  // --- tracks
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
  track_nspoints(nullptr),
  track_nclusters(nullptr),
  track_nhits(nullptr),
  cluster_nhits(nullptr),
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
  // --- clusters

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
  /*
  // Make sure the vectors are empty and counters set to zero.
  // Truth
  fPDG->clear();
  mcpart_position->clear();
  mcpart_positiont->clear();
  mcpart_endposition->clear();
  mcpart_endpositiont->clear();
  mcpart_momentum->clear();
  mcpart_momentume->clear();
  mcpart_momentump->clear();
  mcpart_momentumpt->clear();
  mcpart_momentummass->clear();
  mcpart_momentum->clear();
  mcpart_momentume->clear();
  mcpart_momentump->clear();
  mcpart_momentumpt->clear();
  mcpart_momentummass->clear();
 */
  mctruth_pdg = 0;
  mctruth_isnc = false;
  mctruth_iscc = false;
  std::cout << "set mctruth bools";
  
  mcpart_pdg_hm->clear();
  mcpart_positiont_hm->clear();
  mcpart_momentume_hm->clear();
  mcpart_trackid_hm  ->clear();
  mcpart_vx_hm  ->clear();
  mcpart_vy_hm  ->clear();
  mcpart_vz_hm  ->clear();
  mcpart_endx_hm  ->clear();
  mcpart_endy_hm  ->clear();
  mcpart_endz_hm  ->clear();
  mcpart_pdg          ->clear();
  mcpart_positiont    ->clear();
  mcpart_momentum_e    ->clear();
  mcpart_momentum_p    ->clear();
  mcpart_momentum_pt   ->clear();
  mcpart_momentum_mass ->clear();
  mcpart_trackid      ->clear();
  mcpart_g4id         ->clear();
  std::cout << "cleared mcpart";

  // Reco
  pfpart_number = 0;
  pfpart_ndaughthers ->clear();
  pfpart_id   ->clear();
  pfpart_pdg  ->clear();
  pfpart_isprimary   ->clear();

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

  cluster_nhits  ->clear();
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
  const simb::MCNeutrino& nu = mctruths[0]->GetNeutrino();
  mctruth_isnc = nu.CCNC();
  mctruth_iscc = !nu.CCNC();
  mctruth_pdg = nu.Nu().PdgCode();
/*
  for(const art::Ptr<simb::MCParticle> &part : mcparts){
    mcpart_pdg->push_back(part->PdgCode());
    mcpart_positiont->push_back(part->Position().T());
    mcpart_momentume->push_back(part->Momentum().E());
    mcpart_trackid->push_back(part->TrackId());
  }
*/
  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();
  for(const auto& particle: particles) {
    const simb::MCParticle* this_particle = particle.second;
    mcpart_trackid->push_back(this_particle->TrackId());
    mcpart_pdg->push_back(this_particle->PdgCode());
    mcpart_positiont->push_back(this_particle->Position().T());
    mcpart_momentum_e->push_back(this_particle->E());
    mcpart_momentum_p->push_back(this_particle->P());
    mcpart_momentum_pt->push_back(this_particle->Pt());
    mcpart_momentum_mass->push_back(this_particle->Mass());
  }

  // =========================================================
  // Reco
  // Handles
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  art::Handle< std::vector<recob::Track> > trackHandle;
  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  art::Handle< std::vector<recob::Hit> > hitHandle;
  art::Handle< std::vector<recob::SpacePoint> > spointHandle;
  // Object vectors
  std::vector< art::Ptr<recob::PFParticle> > pfps;
  std::vector< art::Ptr<recob::Track> > tracks;
  std::vector< art::Ptr<recob::Cluster> > clusters;
  std::vector< art::Ptr<recob::Hit> > hits;
  std::vector< art::Ptr<recob::SpacePoint> > spoints;

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
  if(evt.getByLabel(fTrackLabel, trackHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "trackHandle not valid" << std::endl;}
        return;
    }
    art::fill_ptr_vector(tracks, trackHandle);
  }
  if(evt.getByLabel(fClusterLabel, clusterHandle)){
  	if (!pfpHandle.isValid()) {
        if(fVerbose){std::cout << "clusterHandle not valid" << std::endl;}
        return;
    }
  	art::fill_ptr_vector(clusters, clusterHandle);
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
  art::FindManyP<recob::Track> trackAssn(pfps, evt, fTrackLabel);
  art::FindManyP<recob::SpacePoint> spointAssn(pfps, evt, fSpacePointLabel);
  art::FindManyP<recob::Cluster> clusterAssn(pfps, evt, fPFParticleLabel);
  art::FindManyP<recob::Hit> hitAssn(tracks, evt, fTrackLabel);
  art::FindManyP<recob::Hit> hitClusterAssn(clusters, evt, fHitLabel);
  art::FindManyP<anab::Calorimetry> calAssn(tracks, evt, fCalLabel);
  art::FindManyP<anab::ParticleID> pidAssn(tracks, evt, fPIDLabel);
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
    	const std::vector< art::Ptr< recob::Hit> > thishits = hitClusterAssn.at(cluster.key());
    	if(fVerbose){std::cout << "thishits.size(): " << thishits.size() << std::endl;}
    	cluster_hits.insert(cluster_hits.end(), thishits.begin(), thishits.end());
    	cluster_nhits->push_back(thishits.size());
    }
    track_nhits->push_back(cluster_hits.size());
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

      // const std::vector<art::Ptr<recob::SpacePoint>> spoints = spointAssn.at(tracks[0].key());
      // if(spoints.empty()){
      // 	track_nspoints->push_back(0);
      // 	std::cout << "empty space point";
      // }
      // else{
      // 	track_nspoints->push_back(spoints.size());
      // 	std::cout << "found space point";
      // }

      // // Get the hits from the track by accessing the clusters
      // std::vector<art::Ptr<recob::Hit> > cluster_hits;
      // const std::vector< art::Ptr< recob::Cluster> > thisclusters = clusterAssn.at(tracks[0].key());
      // for (const auto& cluster: thisclusters){
      //   const std::vector< art::Ptr< recob::Hit> > thishits = hitAssn.at(cluster.key());
      //   cluster_hits.insert(cluster_hits.end(), thishits.begin(), thishits.end());
      //   cluster_nhits->push_back(thishits.size());
      // }
      // track_nhits->push_back(cluster_hits.size());
      // track_nclusters->push_back(thisclusters.size());

      // Get the hits associated to the track directly
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
      /*
        hits_starttick->push_back(hit->StartTick());
        hits_endtick->push_back(hit->EndTick());
        hits_planeid->push_back(hit->WireID().Plane);
        hits_adc->push_back(hit->SummedADC());
        // hits_channel->push_back(hit->Channel());
        // hits_view->push_back(hit->View());
        // hits_wireid->push_back(hit->WireID());
        // hits_rms->push_back(hit->RMS());
        */
      } // hits
      hits_starttick  ->push_back(temp_start);
      hits_endtick    ->push_back(temp_end);
      hits_planeid ->push_back(temp_planeid);
      hits_adc     ->push_back(temp_adc);  

      // Use the hits to get the G4ID
      auto this_g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clock_data, trackhits, true); // rollupUnsavedIDs
      bool valid_g4id = TruthMatchUtils::Valid(this_g4id);
      if (!valid_g4id) this_g4id = -5;
      mcpart_g4id->push_back(this_g4id);
      // For the truth matching, get the MCParticle which trackId is the same as the G4ID
      for(const art::Ptr<simb::MCParticle> &part : mcparts){
        auto temp_trackid = part->TrackId();
        if(temp_trackid!=this_g4id) continue;
        mcpart_trackid_hm->push_back(part->TrackId());
        mcpart_pdg_hm->push_back(part->PdgCode());
        mcpart_positiont_hm->push_back(part->Position().T());
        mcpart_momentume_hm->push_back(part->Momentum().E());
        mcpart_vx_hm->push_back(part->Vx());
        mcpart_vy_hm->push_back(part->Vy());
        mcpart_vz_hm->push_back(part->Vz());
        mcpart_endx_hm->push_back(part->EndY());
        mcpart_endy_hm->push_back(part->EndY());
        mcpart_endz_hm->push_back(part->EndZ());
      }

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
  fTree->Branch("mcpart_pdg",           &mcpart_pdg);
  fTree->Branch("mcpart_positiont",     &mcpart_positiont);
  fTree->Branch("mcpart_momentum_e",    &mcpart_momentum_e);
  fTree->Branch("mcpart_momentum_p",    &mcpart_momentum_p);
  fTree->Branch("mcpart_momentum_pt",   &mcpart_momentum_pt);  
  fTree->Branch("mcpart_momentum_mass", &mcpart_momentum_mass);  
  fTree->Branch("mcpart_trackid",       &mcpart_trackid);
  fTree->Branch("mcpart_g4id",          &mcpart_g4id);
  fTree->Branch("mcpart_pdg_hm",       &mcpart_pdg_hm);
  fTree->Branch("mcpart_positiont_hm", &mcpart_positiont_hm);
  fTree->Branch("mcpart_momentume_hm", &mcpart_momentume_hm);
  fTree->Branch("mcpart_trackid_hm",   &mcpart_trackid_hm);
  fTree->Branch("mcpart_vx_hm",        &mcpart_vx_hm);
  fTree->Branch("mcpart_vy_hm",        &mcpart_vy_hm);
  fTree->Branch("mcpart_vz_hm",        &mcpart_vz_hm);
  fTree->Branch("mcpart_endx_hm",      &mcpart_endx_hm);
  fTree->Branch("mcpart_endy_hm",      &mcpart_endy_hm);
  fTree->Branch("mcpart_endz_hm",      &mcpart_endz_hm);
  // Reco
  fTree->Branch("pfpart_number",      &pfpart_number, "nPFParticles/i");
  fTree->Branch("pfpart_isprimary",   &pfpart_isprimary);
  fTree->Branch("pfpart_ndaughthers", &pfpart_ndaughthers);
  fTree->Branch("pfpart_id",          &pfpart_id);
  fTree->Branch("pfpart_pdg",         &pfpart_pdg);
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
  fTree->Branch("cluster_hits",      &cluster_nhits);
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

void sbnd::ValidateTracks::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::ValidateTracks)
