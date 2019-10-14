////////////////////////////////////////////////////////////////////////
// Class:       NuMuRecoTree
// Module Type: analyzer
// File:        NuMuRecoTree_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <algorithm>

namespace sbnd {

  class NuMuRecoTree : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> TrackModuleLabel {
        Name("TrackModuleLabel"),
        Comment("tag of track producer data product")
      };

      fhicl::Atom<art::InputTag> HitModuleLabel {
        Name("HitModuleLabel"),
        Comment("tag of hit producer data product")
      };

      fhicl::Atom<art::InputTag> ShowerModuleLabel {
        Name("ShowerModuleLabel"),
        Comment("tag of shower producer data product")
      };

      fhicl::Atom<art::InputTag> PidModuleLabel {
        Name("PidModuleLabel"),
        Comment("tag of PID producer data product")
      };

      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of calorimetry producer data product")
      };

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
        Name("fitter"),
      };

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit NuMuRecoTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetMuonVars();

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Hit> > HitVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fHitModuleLabel; ///< name of TPC track producer
    art::InputTag fShowerModuleLabel; ///< name of TPC track producer
    art::InputTag fPidModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of TPC track producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    TPCGeoAlg fTpcGeo;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    
    // Tree (One entry per primary muon)
    TTree *fMuonTree;

    //Muon tree parameters
    bool is_cosmic;         // True origin of PFP is cosmic
    bool is_nu;             // True origin of PFP is nu in AV
    double time;
    double vtx_x;
    double vtx_y;
    double vtx_z;
    double end_x;
    double end_y;
    double end_z;
    double length;
    double contained_length;
    double momentum;
    double theta;
    double phi;
    double e_dep;
    int n_hits;
    double hit_area;
    int n_cr_tracks;
    int n_nu_showers;
    int n_nu_tracks;
    double track_vtx_x;
    double track_vtx_y;
    double track_vtx_z;
    double track_end_x;
    double track_end_y;
    double track_end_z;
    double track_length;
    double track_range_momentum;
    double track_mcs_momentum;
    double track_theta;
    double track_phi;
    double track_e_dep;
    int track_n_hits;
    int track_n_hits_true;
    double track_hit_area;
    double track_hit_area_true;
    double track_mu_chi2;
    double track_pi_chi2;
    double track_p_chi2;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);
    
    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

  }; // class NuMuRecoTree


  // Constructor
  NuMuRecoTree::NuMuRecoTree(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTrackModuleLabel     (config().TrackModuleLabel())
    , fHitModuleLabel       (config().HitModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPidModuleLabel       (config().PidModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fMcsFitter            (config().fitter)
  {

  } // NuMuRecoTree()


  void NuMuRecoTree::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fMuonTree = tfs->make<TTree>("muons", "muons");

    fMuonTree->Branch("is_cosmic",       &is_cosmic,       "is_cosmic/O");
    fMuonTree->Branch("is_nu",           &is_nu,           "is_nu/O");
    fMuonTree->Branch("time",       &time,       "time/D");
    fMuonTree->Branch("vtx_x",       &vtx_x,       "vtx_x/D");
    fMuonTree->Branch("vtx_y",       &vtx_y,       "vtx_y/D");
    fMuonTree->Branch("vtx_z",       &vtx_z,       "vtx_z/D");
    fMuonTree->Branch("end_x",       &end_x,       "end_x/D");
    fMuonTree->Branch("end_y",       &end_y,       "end_y/D");
    fMuonTree->Branch("end_z",       &end_z,       "end_z/D");
    fMuonTree->Branch("length",       &length,       "length/D");
    fMuonTree->Branch("contained_length",       &contained_length,       "contained_length/D");
    fMuonTree->Branch("momentum",       &momentum,       "momentum/D");
    fMuonTree->Branch("theta",       &theta,       "theta/D");
    fMuonTree->Branch("phi",       &phi,       "phi/D");
    fMuonTree->Branch("e_dep",       &e_dep,       "e_dep/D");
    fMuonTree->Branch("n_hits",       &n_hits,       "n_hits/I");
    fMuonTree->Branch("hit_area",       &hit_area,       "hit_area/D");
    fMuonTree->Branch("n_cr_tracks",       &n_cr_tracks,       "n_cr_tracks/I");
    fMuonTree->Branch("n_nu_showers",       &n_nu_showers,       "n_nu_showers/I");
    fMuonTree->Branch("n_nu_tracks",       &n_nu_tracks,       "n_nu_tracks/I");
    fMuonTree->Branch("track_vtx_x",       &track_vtx_x,       "track_vtx_x/D");
    fMuonTree->Branch("track_vtx_y",       &track_vtx_y,       "track_vtx_y/D");
    fMuonTree->Branch("track_vtx_z",       &track_vtx_z,       "track_vtx_z/D");
    fMuonTree->Branch("track_end_x",       &track_end_x,       "track_end_x/D");
    fMuonTree->Branch("track_end_y",       &track_end_y,       "track_end_y/D");
    fMuonTree->Branch("track_end_z",       &track_end_z,       "track_end_z/D");
    fMuonTree->Branch("track_length",         &track_length,       "track_length/D");
    fMuonTree->Branch("track_range_momentum", &track_range_momentum,       "track_range_momentum/D");
    fMuonTree->Branch("track_mcs_momentum",   &track_mcs_momentum,       "track_mcs_momentum/D");
    fMuonTree->Branch("track_theta",       &track_theta,       "track_theta/D");
    fMuonTree->Branch("track_phi",       &track_phi,       "track_phi/D");
    fMuonTree->Branch("track_e_dep",       &track_e_dep,       "track_e_dep/D");
    fMuonTree->Branch("track_n_hits",       &track_n_hits,       "track_n_hits/I");
    fMuonTree->Branch("track_n_hits_true",       &track_n_hits_true,       "track_n_hits_true/I");
    fMuonTree->Branch("track_hit_area",       &track_hit_area,       "track_hit_area/D");
    fMuonTree->Branch("track_hit_area_true",       &track_hit_area_true,       "track_hit_area_true/D");
    fMuonTree->Branch("track_mu_chi2",       &track_mu_chi2,       "track_mu_chi2/D");
    fMuonTree->Branch("track_pi_chi2",       &track_pi_chi2,       "track_pi_chi2/D");
    fMuonTree->Branch("track_p_chi2",       &track_p_chi2,       "track_p_chi2/D");
   

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Ana Module -------------------"<<std::endl;

  }// NuMuRecoTree::beginJob()


  void NuMuRecoTree::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get truth info and matching
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get PFParticles from pandora
    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);
    if( !pfParticleHandle.isValid() ){
      if(fVerbose) std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

    // Get cosmic and neutrino particles
    PFParticleVector cr_particles;
    PFParticleVector nu_particles;
    this->GetFinalStatePFParticleVectors(pfParticleMap, cr_particles, nu_particles);

    // Get neutrino tracks and showers
    TrackVector nu_tracks;
    ShowerVector nu_showers;
    this->CollectTracksAndShowers(nu_particles, pfParticleHandle, event, nu_tracks, nu_showers);

    // Get cosmic ray tracks and showers
    TrackVector cr_tracks;
    ShowerVector cr_showers;
    this->CollectTracksAndShowers(cr_particles, pfParticleHandle, event, cr_tracks, cr_showers);

    // Get PFParticle to track associations
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTrackModuleLabel);
    
    // Get track handle
    auto trackHandle = event.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(trackHandle, event, fTrackModuleLabel);
    // Get track to PID associations
    art::FindMany<anab::ParticleID> findManyPid(trackHandle, event, fPidModuleLabel);
    // Get track to calorimetry associations
    art::FindMany<anab::Calorimetry> findManyCalo(trackHandle, event, fCaloModuleLabel);

    // Get shower handle
    auto showerHandle = event.getValidHandle<std::vector<recob::Shower>>(fShowerModuleLabel);
    // Get shower to hit associations
    art::FindManyP<recob::Hit> findManyHitsShower(showerHandle, event, fShowerModuleLabel);

    // Get hit handle
    auto hitHandle = event.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                        MUON RECO ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    // Truth matching
    std::map<int, TrackVector> match_nu_tracks;
    for(auto const& track : nu_tracks){
      HitVector hits = findManyHits.at(track->ID());
      int id = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      match_nu_tracks[id].push_back(track);
    }

    std::map<int, TrackVector> match_cr_tracks;
    for(auto const& track : cr_tracks){
      HitVector hits = findManyHits.at(track->ID());
      int id = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      match_cr_tracks[id].push_back(track);
    }

    std::map<int, ShowerVector> match_nu_showers;
    for(auto const& shower : nu_showers){
      HitVector hits = findManyHitsShower.at(shower.key());
      int id = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      match_nu_showers[id].push_back(shower);
    }

    std::map<int, HitVector> match_hits;
    for(size_t i = 0; i < hitHandle->size(); i++){
      const art::Ptr<recob::Hit> hit(hitHandle, i);
      int id = RecoUtils::TrueParticleID(hit, false);
      match_hits[id].push_back(hit);
    }

    // Put them in a map for easier access
    for (auto const& particle: (*particleHandle)){
      // Store muon if primary and stable and inside TPC
      if(std::abs(particle.PdgCode()) != 13) continue;
      if(particle.Mother() != 0) continue;
      if(particle.StatusCode() != 1) continue;
      if(!fTpcGeo.InVolume(particle)) continue;

      ResetMuonVars();

      // True variables
      int id = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(id);
      if(truth->Origin() == simb::kBeamNeutrino) is_nu = true;
      if(truth->Origin() == simb::kCosmicRay) is_cosmic = true;

      time = particle.T();

      vtx_x = particle.Vx();
      vtx_y = particle.Vy();
      vtx_z = particle.Vz();

      end_x = particle.EndX();
      end_y = particle.EndY();
      end_z = particle.EndZ();

      length = particle.Trajectory().TotalLength();
      contained_length = fTpcGeo.TpcLength(particle);
      momentum = particle.P();
      std::pair<TVector3, TVector3> se = fTpcGeo.CrossingPoints(particle);
      theta = (se.second-se.first).Theta();
      phi = (se.second-se.first).Phi();

      e_dep = 0;
      for(size_t i = 0; i < particle.NumberTrajectoryPoints() - 1; i++){
        geo::Point_t pos {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
        if(fTpcGeo.InFiducial(pos, 0)) e_dep += particle.E(i) - particle.E(i+1);
      }

      if(match_hits.find(id) != match_hits.end()){
        n_hits = match_hits[id].size();
        hit_area = 0;
        for(int i = 0; i < n_hits; i++){
          hit_area += match_hits[id][i]->Integral();
        }
      }

      // Number of cr tracks
      if(match_cr_tracks.find(id) != match_cr_tracks.end()){
        n_cr_tracks = match_cr_tracks[id].size();
      }
      // Number of nu showers
      if(match_nu_showers.find(id) != match_nu_showers.end()){
        n_nu_showers = match_nu_showers[id].size();
      }
      // Number of nu tracks
      recob::Track track;
      if(match_nu_tracks.find(id) != match_nu_tracks.end()){
        n_nu_tracks = match_nu_tracks[id].size();
        double longest_length = -99999;
        for(int i = 0; i < n_nu_tracks; i++){
          if(match_nu_tracks[id][i]->Length() > longest_length){
            track = *match_nu_tracks[id][i];
            longest_length = match_nu_tracks[id][i]->Length();
          }
        }
      }

      // Track variables
      if(n_nu_tracks < 1){ 
        fMuonTree->Fill();
        continue;
      }

      track_vtx_x = track.Vertex().X();
      track_vtx_y = track.Vertex().Y();
      track_vtx_z = track.Vertex().Z();

      track_end_x = track.End().X();
      track_end_y = track.End().Y();
      track_end_z = track.End().Z();

      track_length = track.Length();
      track_range_momentum = fRangeFitter.GetTrackMomentum(track_length, 13);
      track_mcs_momentum = fMcsFitter.fitMcs(track).bestMomentum();
      track_theta = track.Theta();
      track_phi = track.Phi();
       
      std::vector<const anab::Calorimetry*> calos = findManyCalo.at(track.ID());
      for(size_t i = 0; i < calos.size(); i++){
        if(calos[i]->PlaneID().Plane != 2) continue;
        track_e_dep = calos[i]->KineticEnergy();
      }

      HitVector hits = findManyHits.at(track.ID());
      track_n_hits = hits.size();
      track_hit_area = 0;
      track_hit_area_true = 0;
      for(int i = 0; i < track_n_hits; i++){
        track_hit_area += hits[i]->Integral();
        int id_true = RecoUtils::TrueParticleID(hits[i], false);
        if(id_true == id){
          track_n_hits_true++;
          track_hit_area_true += hits[i]->Integral();
        }
      }

      std::vector<const anab::ParticleID*> pids = findManyPid.at(track.ID());
      for(size_t i = 0; i < pids.size(); i++){
        if(pids[i]->PlaneID().Plane != 2) continue;
        track_mu_chi2 = pids[i]->Chi2Muon();
        track_pi_chi2 = pids[i]->Chi2Pion();
        track_p_chi2 = pids[i]->Chi2Proton();
      }

      fMuonTree->Fill();
    }
    
  } // NuMuRecoTree::analyze()


  void NuMuRecoTree::endJob(){

  } // NuMuRecoTree::endJob()

  void NuMuRecoTree::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void NuMuRecoTree::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles){

      for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
          const art::Ptr<recob::PFParticle> pParticle(it->second);

          // Only look for primary particles
          if (!pParticle->IsPrimary()) continue;

          // Check if this particle is identified as the neutrino
          const int pdg(pParticle->PdgCode());
          const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

          // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
          if (!isNeutrino){
              crParticles.push_back(pParticle);
              continue;
          }

          // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
          for (const size_t daughterId : pParticle->Daughters()){
              if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                  std::cout << "  Invalid PFParticle collection!" <<"\n";

              nuParticles.push_back(pfParticleMap.at(daughterId));
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void NuMuRecoTree::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
  {
      // Get the associations between PFParticles and tracks/showers from the event
      art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, fTrackModuleLabel);
      art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, evt, fShowerModuleLabel);
     
      for (const art::Ptr<recob::PFParticle> &pParticle : particles){
          const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
          const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
          const unsigned int nTracks(associatedTracks.size());
          const unsigned int nShowers(associatedShowers.size());

          // Check if the PFParticle has no associated tracks or showers
          if (nTracks == 0 && nShowers == 0){
              std::cout << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
              continue;
          }

          // Check if there is an associated track
          if (nTracks == 1 && nShowers == 0){
              tracks.push_back(associatedTracks.front());
              continue;
          }

          // Check if there is an associated shower
          if (nTracks == 0 && nShowers == 1){
              showers.push_back(associatedShowers.front());
              continue;
          }

          std::cout << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self() << "\n";
      }
  }

  // Reset the tree variables
  void NuMuRecoTree::ResetMuonVars(){
    is_cosmic = false;
    is_nu = false;
    time = -99999;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    end_x = -99999;
    end_y = -99999;
    end_z = -99999;
    length = -99999;
    contained_length = -99999;
    momentum = -99999;
    theta = -99999;
    phi = -99999;
    e_dep = -99999;
    n_hits = 0;
    hit_area = -99999;
    n_cr_tracks = 0;
    n_nu_showers = 0;
    n_nu_tracks = 0;
    track_vtx_x = -99999;
    track_vtx_y = -99999;
    track_vtx_z = -99999;
    track_end_x = -99999;
    track_end_y = -99999;
    track_end_z = -99999;
    track_length = -99999;
    track_range_momentum = -99999;
    track_mcs_momentum = -99999;
    track_theta = -99999;
    track_phi = -99999;
    track_e_dep = -99999;
    track_n_hits = 0;
    track_n_hits_true = 0;
    track_hit_area = -99999;
    track_hit_area_true = -99999;
    track_mu_chi2 = -99999;
    track_pi_chi2 = -99999;
    track_p_chi2 = -99999;
    
  }
    
  
  //------------------------------------------------------------------------------------------------------------------------------------------

  DEFINE_ART_MODULE(NuMuRecoTree)
} // namespace sbnd

