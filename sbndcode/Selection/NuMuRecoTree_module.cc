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
#include "sbndcode/CosmicId/Algs/StoppingParticleCosmicIdAlg.h"

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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

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

      fhicl::Table<StoppingParticleCosmicIdAlg::Config> SPTagAlg {
        Name("SPTagAlg"),
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
    void ResetParticleVars();
    void ResetNuVars();

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
    StoppingParticleCosmicIdAlg  fStopTagger;

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorClocks const* fDetectorClocks;
    detinfo::DetectorProperties const* fDetectorProperties;
    double readoutWindow;
    double driftTime;
    
    // Tree (One entry per primary muon)
    TTree *fParticleTree;

    //Muon tree parameters
    bool is_cosmic;         // True origin of PFP is cosmic
    bool is_nu;             // True origin of PFP is nu
    bool is_cc;             // True origin of PFP is CC nu
    int nu_pdg;
    int pdg;
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
    std::string end_process;
    std::string fin_process;
    double vtx_x_tpc;
    double vtx_y_tpc;
    double vtx_z_tpc;
    double end_x_tpc;
    double end_y_tpc;
    double end_z_tpc;
    double stopping_chi2;
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
    double mu_range_momentum;
    double mu_mcs_momentum;
    double track_dca;
    double track_ave_angle;
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
    double track_flat_chi2;

    // Tree (One entry per reconstructed numuCC)
    TTree *fNuMuTree;
    double nu_vtx_x;
    double nu_vtx_y;
    double nu_vtx_z;
    double true_nu_energy;
    double true_had_energy;
    double true_mu_momentum;
    int nu_n_tracks;
    int nu_n_showers;
    bool mu_contained;
    bool tracks_contained;
    double reco_mu_momentum;
    double reco_track_energy;
    double reco_track_hit_area;
    double reco_shower_hit_area;
    double reco_track_hit_energy;
    double reco_shower_hit_energy;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);
    
    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

    std::pair<double, double> XLimitsTPC(const simb::MCParticle& particle);

    std::string FinalProcess(int id, const std::map<int, simb::MCParticle>& particleMap) const;

    double AverageDCA(const recob::Track& track);

    double HitEnergy(art::Ptr<recob::Hit> hit);

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
    , fStopTagger           (config().SPTagAlg())
  {

  } // NuMuRecoTree()


  void NuMuRecoTree::beginJob()
  {

    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    readoutWindow  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    driftTime = (2.*fGeometryService->DetHalfWidth())/fDetectorProperties->DriftVelocity(); // [us]

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fParticleTree = tfs->make<TTree>("particles", "particles");

    fParticleTree->Branch("is_cosmic",        &is_cosmic);
    fParticleTree->Branch("is_nu",            &is_nu);
    fParticleTree->Branch("is_cc",            &is_cc);
    fParticleTree->Branch("nu_pdg",           &nu_pdg);
    fParticleTree->Branch("pdg",              &pdg);
    fParticleTree->Branch("time",             &time);
    fParticleTree->Branch("vtx_x",            &vtx_x);
    fParticleTree->Branch("vtx_y",            &vtx_y);
    fParticleTree->Branch("vtx_z",            &vtx_z);
    fParticleTree->Branch("end_x",            &end_x);
    fParticleTree->Branch("end_y",            &end_y);
    fParticleTree->Branch("end_z",            &end_z);
    fParticleTree->Branch("length",           &length);
    fParticleTree->Branch("contained_length", &contained_length);
    fParticleTree->Branch("momentum",         &momentum);
    fParticleTree->Branch("theta",            &theta);
    fParticleTree->Branch("phi",              &phi);
    fParticleTree->Branch("end_process",      &end_process);
    fParticleTree->Branch("fin_process",      &fin_process);
    fParticleTree->Branch("vtx_x_tpc",        &vtx_x_tpc);
    fParticleTree->Branch("vtx_y_tpc",        &vtx_y_tpc);
    fParticleTree->Branch("vtx_z_tpc",        &vtx_z_tpc);
    fParticleTree->Branch("end_x_tpc",        &end_x_tpc);
    fParticleTree->Branch("end_y_tpc",        &end_y_tpc);
    fParticleTree->Branch("end_z_tpc",        &end_z_tpc);
    fParticleTree->Branch("stopping_chi2",    &stopping_chi2);
    fParticleTree->Branch("e_dep",            &e_dep);
    fParticleTree->Branch("n_hits",           &n_hits);
    fParticleTree->Branch("hit_area",         &hit_area);
    fParticleTree->Branch("n_cr_tracks",      &n_cr_tracks);
    fParticleTree->Branch("n_nu_showers",     &n_nu_showers);
    fParticleTree->Branch("n_nu_tracks",      &n_nu_tracks);
    fParticleTree->Branch("track_vtx_x",      &track_vtx_x);
    fParticleTree->Branch("track_vtx_y",      &track_vtx_y);
    fParticleTree->Branch("track_vtx_z",      &track_vtx_z);
    fParticleTree->Branch("track_end_x",      &track_end_x);
    fParticleTree->Branch("track_end_y",      &track_end_y);
    fParticleTree->Branch("track_end_z",      &track_end_z);
    fParticleTree->Branch("track_length",         &track_length);
    fParticleTree->Branch("track_range_momentum", &track_range_momentum);
    fParticleTree->Branch("track_mcs_momentum",   &track_mcs_momentum);
    fParticleTree->Branch("mu_range_momentum", &mu_range_momentum);
    fParticleTree->Branch("mu_mcs_momentum",   &mu_mcs_momentum);
    fParticleTree->Branch("track_dca",        &track_dca);
    fParticleTree->Branch("track_ave_angle",  &track_ave_angle);
    fParticleTree->Branch("track_theta",      &track_theta);
    fParticleTree->Branch("track_phi",        &track_phi);
    fParticleTree->Branch("track_e_dep",      &track_e_dep);
    fParticleTree->Branch("track_n_hits",     &track_n_hits);
    fParticleTree->Branch("track_n_hits_true",&track_n_hits_true);
    fParticleTree->Branch("track_hit_area",   &track_hit_area);
    fParticleTree->Branch("track_hit_area_true",  &track_hit_area_true);
    fParticleTree->Branch("track_mu_chi2",    &track_mu_chi2);
    fParticleTree->Branch("track_pi_chi2",    &track_pi_chi2);
    fParticleTree->Branch("track_p_chi2",     &track_p_chi2);
    fParticleTree->Branch("track_flat_chi2",     &track_flat_chi2);
   
    fNuMuTree = tfs->make<TTree>("numu", "numu");

    fNuMuTree->Branch("nu_vtx_x",          &nu_vtx_x);
    fNuMuTree->Branch("nu_vtx_y",          &nu_vtx_y);
    fNuMuTree->Branch("nu_vtx_z",          &nu_vtx_z);
    fNuMuTree->Branch("true_nu_energy",    &true_nu_energy);
    fNuMuTree->Branch("true_mu_momentum",  &true_mu_momentum);
    fNuMuTree->Branch("true_had_energy",   &true_had_energy);
    fNuMuTree->Branch("nu_n_tracks",       &nu_n_tracks);
    fNuMuTree->Branch("nu_n_showers",      &nu_n_showers);
    fNuMuTree->Branch("mu_contained",      &mu_contained);
    fNuMuTree->Branch("tracks_contained",  &tracks_contained);
    fNuMuTree->Branch("reco_mu_momentum",  &reco_mu_momentum);
    fNuMuTree->Branch("reco_track_energy", &reco_track_energy);
    fNuMuTree->Branch("reco_track_hit_area", &reco_track_hit_area);
    fNuMuTree->Branch("reco_shower_hit_area", &reco_shower_hit_area);
    fNuMuTree->Branch("reco_track_hit_energy", &reco_track_hit_energy);
    fNuMuTree->Branch("reco_shower_hit_energy", &reco_shower_hit_energy);

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
    art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, event, fShowerModuleLabel);
    
    // Get track handle
    auto trackHandle = event.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(trackHandle, event, fTrackModuleLabel);
    // Get track to PID associations
    art::FindMany<anab::ParticleID> findManyPid(trackHandle, event, fPidModuleLabel);
    // Get track to calorimetry associations
    art::FindManyP<anab::Calorimetry> findManyCalo(trackHandle, event, fCaloModuleLabel);

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


    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      // Only interested in track-like muons, pions and protons
      pdg = particle.PdgCode();
      if(!(std::abs(pdg) == 13||std::abs(pdg) == 211||std::abs(pdg) == 2212)) continue;
      // Only want primary particles
      if(particle.Mother() != 0) continue;
      // Only want stable particles (post fsi)
      if(particle.StatusCode() != 1) continue;
      // Only want particles that are inside the TPC
      if(!fTpcGeo.InVolume(particle)) continue;
      // Check that it crosses in the reconstructable window
      double time = particle.T() * 1e-3; //[us]
      if (time <= -driftTime || time >= readoutWindow) continue;
      // Get the minimum and maximum |x| position in the TPC
      std::pair<double, double> xLimits = XLimitsTPC(particle);
      // Calculate the expected time of arrival of those points
      double minTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.second)/fDetectorProperties->DriftVelocity(); 
      double maxTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.first)/fDetectorProperties->DriftVelocity(); 
      // If both times are below or above the readout window time then skip
      if((minTime < 0 && maxTime < 0) || (minTime > readoutWindow && maxTime < readoutWindow)) continue;
      int id = particle.TrackId();
      particles[id] = particle;
    }

    // Put them in a map for easier access
    for (auto const& part: particles){
      simb::MCParticle particle = part.second;

      ResetParticleVars();
      pdg = particle.PdgCode();

      // True variables
      int id = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(id);
      if(truth->Origin() == simb::kBeamNeutrino){ 
        is_nu = true;
        nu_pdg = truth->GetNeutrino().Nu().PdgCode();
        if(truth->GetNeutrino().CCNC() == simb::kCC) is_cc = true;
      }
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

      end_process = particle.EndProcess();
      fin_process = FinalProcess(id, particles);

      vtx_x_tpc = se.first.X();
      vtx_y_tpc = se.first.Y();
      vtx_z_tpc = se.first.Z();

      end_x_tpc = se.second.X();
      end_y_tpc = se.second.Y();
      end_z_tpc = se.second.Z();

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
        fParticleTree->Fill();
        continue;
      }

      track_vtx_x = track.Vertex().X();
      track_vtx_y = track.Vertex().Y();
      track_vtx_z = track.Vertex().Z();

      track_end_x = track.End().X();
      track_end_y = track.End().Y();
      track_end_z = track.End().Z();

      track_length = track.Length();
      track_range_momentum = fRangeFitter.GetTrackMomentum(track_length, pdg);
      track_mcs_momentum = fMcsFitter.fitMcs(track, pdg).bestMomentum();
      mu_range_momentum = fRangeFitter.GetTrackMomentum(track_length, 13);
      mu_mcs_momentum = fMcsFitter.fitMcs(track, 13).bestMomentum();
      track_dca = AverageDCA(track); 
      std::vector<float> angles = fMcsFitter.fitMcs(track, 13).scatterAngles();
      track_ave_angle = std::accumulate(angles.begin(), angles.end(), 0)/angles.size();
      track_theta = track.Theta();
      track_phi = track.Phi();
       
      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
      // Loop over planes (Y->V->U) and choose the next plane's calorimetry if there are 1.5x more points (collection plane more reliable)
      if(calos.size()==0) continue;
      size_t nhits = 0;
      art::Ptr<anab::Calorimetry> calo = calos[0];
      size_t best_plane = 0;
      for( size_t i = calos.size(); i > 0; i--){
        if(calos[i-1]->dEdx().size() > nhits*1.5){
          nhits = calos[i-1]->dEdx().size();
          calo = calos[i-1];
          best_plane = i-1;
        }
      }

      stopping_chi2 = fStopTagger.StoppingChiSq(track.End(), calos);

      track_e_dep = calo->KineticEnergy();

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
      std::vector<double> flat_chi2 = fStopTagger.FlatChi2(calos);
      for(size_t i = 0; i < pids.size(); i++){
        if(pids[i]->PlaneID().Plane != best_plane) continue;
        track_mu_chi2 = pids[i]->Chi2Muon();
        track_pi_chi2 = pids[i]->Chi2Pion();
        track_p_chi2 = pids[i]->Chi2Proton();
        if(i < flat_chi2.size()) track_flat_chi2 = flat_chi2[i];
      }

      fParticleTree->Fill();
    }

    //----------------------------------------------------------------------------------------------------------
    //                                        NEUTRINO ENERGY RECO
    //----------------------------------------------------------------------------------------------------------

    // PFParticle truth matching
    std::map<int, PFParticleVector> match_pfps;
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Check if this particle is identified as the neutrino
      if (!pParticle->IsPrimary()) continue;
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || 
                            std::abs(pdg) == pandora::NU_MU || 
                            std::abs(pdg) == pandora::NU_TAU);
      if(!isNeutrino) continue;

      for (const size_t daughterId : pParticle->Daughters()){

        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue;

        // Get the first associated track
        recob::Track tpcTrack = *associatedTracks.front();
        HitVector hits = findManyHits.at(tpcTrack.ID());
        int id = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        match_pfps[id].push_back(pDaughter);
      }
    }

    // Loop over all the particles
    for (auto const& part: particles){
      simb::MCParticle particle = part.second;
      // Find any primary muons
      pdg = particle.PdgCode();
      if(std::abs(pdg) != 13) continue;
      // Get the MCTruth associations
      int id = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(id);
      // Check the truth is a beam neutrino
      if(truth->Origin() != simb::kBeamNeutrino) continue;
      // Check that it's a numuCC interaction
      if(std::abs(truth->GetNeutrino().Nu().PdgCode()) != 14) continue;
      if(truth->GetNeutrino().CCNC() != simb::kCC) continue;
      // Check vertex in AV
      geo::Point_t vtx {particle.Vx(), particle.Vy(), particle.Vz()};
      if(!fTpcGeo.InFiducial(vtx, 0)) continue;

      ResetNuVars();
      // Save true information
      nu_vtx_x = vtx.X();
      nu_vtx_y = vtx.Y();
      nu_vtx_z = vtx.Z();
      true_nu_energy = truth->GetNeutrino().Nu().E();
      true_mu_momentum = particle.P();
      true_had_energy = true_nu_energy - std::sqrt(std::pow(true_mu_momentum, 2) + std::pow(0.10566, 2));

      // Loop over associated PFParticles
      if(match_pfps.find(id) == match_pfps.end()) continue;
      double longest = -99999;
      recob::Track mu_track;
      size_t parent_id = -99999;
      size_t mu_id = -99999;
      for(size_t i = 0; i < match_pfps[id].size(); i++){
        // Get associated track
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(match_pfps[id][i].key()));
        if(associatedTracks.size() != 1) continue;
        recob::Track track = *associatedTracks.front();
        // If track is longest then get the parent PFParticle ID
        if(track.Length() > longest){
          longest = track.Length();
          mu_track = track;
          parent_id = match_pfps[id][i]->Parent();
          mu_id = match_pfps[id][i]->Self();
        }
      }
      if(longest == -99999) continue;

      // For the daughter associated with the primary muon save reco kinematics
      mu_contained = fTpcGeo.InFiducial(mu_track.End(), 1.5);
      if(mu_contained) reco_mu_momentum = fRangeFitter.GetTrackMomentum(mu_track.Length(), 13);
      else reco_mu_momentum = fMcsFitter.fitMcs(mu_track, 13).bestMomentum();

      reco_track_energy = 0;
      reco_track_hit_area = 0;
      reco_shower_hit_area = 0;
      reco_track_hit_energy = 0;
      reco_shower_hit_energy = 0;

      // Get all of the daughters in the PFParticle slice
      if(pfParticleMap.find(parent_id) == pfParticleMap.end()) continue;
      art::Ptr<recob::PFParticle> pParent = pfParticleMap.at(parent_id);
      for (const size_t daughter_id : pParent->Daughters()){
        if (daughter_id == mu_id) continue;

        // Get all associated tracks and showers for each daughter
        if(pfParticleMap.find(daughter_id) == pfParticleMap.end()) continue;
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughter_id);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pDaughter.key()));

        // Add up all of the deposited energy
        for(size_t i = 0; i < associatedTracks.size(); i++){
          nu_n_tracks++;
          // Check if tracks exit
          if(!fTpcGeo.InFiducial(associatedTracks[i]->End(), 1.5)) tracks_contained = false;

          // Loop over planes (Y->V->U) and choose the next plane's calorimetry if there are 1.5x more points (collection plane more reliable)
          std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(associatedTracks[i]->ID());
          if(calos.size()==0) continue;
          size_t nhits = 0;
          size_t best_plane = 0;
          for( size_t j = calos.size(); j > 0; j--){
            if(calos[j-1]->dEdx().size() > nhits*1.5){
              nhits = calos[j-1]->dEdx().size();
              best_plane = j-1;
            }
          }
          // Get the track kinetic energy returned by the calorimetry module
          if(calos[best_plane]->KineticEnergy() > 0) reco_track_energy += calos[best_plane]->KineticEnergy()/1e3; //[GeV]

          // Get the area and energy of hits associated with the tracks
          HitVector hits = findManyHits.at(associatedTracks[i]->ID());
          for(size_t j = 0; j < hits.size(); j++){
            if(hits[j]->WireID().Plane == best_plane){ 
              reco_track_hit_area += hits[j]->Integral();
              reco_track_hit_energy += HitEnergy(hits[j]);
            }
          }
        }

        for(size_t i = 0; i < associatedShowers.size(); i++){
          nu_n_showers++;
          HitVector hits = findManyHitsShower.at(associatedShowers[i]->ID());
          for(size_t j = 0; j < hits.size(); j++){
            // Assume collection plane is best for showers
            if(hits[j]->WireID().Plane == 2){ 
              reco_shower_hit_area += hits[j]->Integral();
              reco_shower_hit_energy += HitEnergy(hits[j]);
            }
          }
        }
      }
      // Save if secondary particles contained
      fNuMuTree->Fill();
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
  void NuMuRecoTree::ResetParticleVars(){
    is_cosmic = false;
    is_nu = false;
    is_cc = false;
    nu_pdg = -99999;
    pdg = -99999;
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
    end_process = "unknown";
    fin_process = "unknown";
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
    mu_range_momentum = -99999;
    mu_mcs_momentum = -99999;
    track_dca = -99999;
    track_ave_angle = -99999;
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
    track_flat_chi2 = -99999;
    
  }

  // Reset the tree variables
  void NuMuRecoTree::ResetNuVars(){
    nu_vtx_x = -99999;
    nu_vtx_y = -99999;
    nu_vtx_z = -99999;
    true_nu_energy = -99999;
    true_mu_momentum = -99999;
    true_had_energy = -99999;
    nu_n_tracks = 0;
    nu_n_showers = 0;
    mu_contained = true;
    tracks_contained = true;
    reco_mu_momentum = -99999;
    reco_track_energy = -99999;
    reco_track_hit_area = -99999;
    reco_shower_hit_area = -99999;
    reco_track_hit_energy = -99999;
    reco_shower_hit_energy = -99999;
  }
    
  std::pair<double, double> NuMuRecoTree::XLimitsTPC(const simb::MCParticle& particle){
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();

    double minimum = 99999;
    double maximum = -99999;

    int nTrajPoints = particle.NumberTrajectoryPoints();
    for (int traj_i = 0; traj_i < nTrajPoints; traj_i++){
      TVector3 trajPoint(particle.Vx(traj_i), particle.Vy(traj_i), particle.Vz(traj_i));
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
        if(std::abs(trajPoint[0]) < minimum) minimum = std::abs(trajPoint[0]);
        if(std::abs(trajPoint[0]) > maximum) maximum = std::abs(trajPoint[0]);
      }
    }
    return std::make_pair(minimum, maximum);                                                                                                                  
  }
  //------------------------------------------------------------------------------------------------------------------------------------------

  // Follow processes down to final particle of same type
  std::string NuMuRecoTree::FinalProcess(int id, const std::map<int, simb::MCParticle>& particles) const{
    int parent_pdg = particles.at(id).PdgCode();
    std::string process = particles.at(id).EndProcess();
    int n_daughters = particles.at(id).NumberDaughters();
    double max_p = -99999;
    int max_id = -99999;

    // Find daughter particle with same PDG and highest momentum
    for(int i = 0; i < n_daughters; i++){
      int daughter_id = particles.at(id).Daughter(i);
      if(particles.find(daughter_id) == particles.end()) continue;
      simb::MCParticle daughter = particles.at(daughter_id);
      int daughter_pdg = daughter.PdgCode();
      if(daughter_pdg != parent_pdg) continue;
      if(daughter.P() < max_p) continue;
      max_p = daughter.P();
      max_id = daughter_id;
      process = daughter.EndProcess();
    }

    if(max_id != -99999) return FinalProcess(max_id, particles);
    return process;
  }
    
  double NuMuRecoTree::AverageDCA(const recob::Track& track){

    TVector3 start = track.Vertex<TVector3>();
    TVector3 end = track.End<TVector3>();
    double denominator = (end - start).Mag();
    size_t npts = track.NumberTrajectoryPoints();
    double aveDCA = 0;
    int usedPts = 0;
    for(size_t i = 0; i < npts; i++){
      TVector3 point = track.LocationAtPoint<TVector3>(i);
      if(!track.HasValidPoint(i)) continue;
      aveDCA += (point - start).Cross(point - end).Mag()/denominator;
      usedPts++;
    }
    return aveDCA/usedPts;
  }

  double NuMuRecoTree::HitEnergy(art::Ptr<recob::Hit> hit){

    double ADCtoEl = 0.02354; //FIXME from calorimetry_sbnd.fcl
    double time = hit->PeakTime();
    double timetick = fDetectorProperties->SamplingRate()*1e-3;
    double presamplings = fDetectorProperties->TriggerOffset();
    time -= presamplings;
    time = time*timetick;
    double tau = fDetectorProperties->ElectronLifetime();
    double correction = std::exp(time/tau);
    return fDetectorProperties->ModBoxCorrection((hit->Integral()/ADCtoEl)*correction)/1e3; //[GeV]
    
  }

  DEFINE_ART_MODULE(NuMuRecoTree)
} // namespace sbnd

