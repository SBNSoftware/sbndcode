////////////////////////////////////////////////////////////////////////
// Class:       CosmicIdTree
// Module Type: analyzer
// File:        CosmicIdTree_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/CosmicId/Algs/CosmicIdAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace sbnd {

  class CosmicIdTree : public art::EDAnalyzer {
  public:

    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("")
      };

    };

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel"),
        Comment("tag of CRT hit producer data product")
      };

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of tpc track producer data product")
      };

      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of tpc calorimetry data product")
      };

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

      fhicl::Table<CosmicIdAlg::Config> CosIdAlg {
        Name("CosIdAlg"),
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CosmicIdTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetTrackVars();
    void ResetPfpVars();

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;

  private:

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);
    
    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTHitLabel;   ///< name of CRT producer
    art::InputTag fCRTTrackLabel;   ///< name of CRT producer
    art::InputTag fTPCTrackLabel; ///< name of CRT producer
    art::InputTag fCaloModuleLabel; ///< name of CRT producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on
    double fBeamTimeMin;
    double fBeamTimeMax;

    TPCGeoAlg fTpcGeo;

    CRTBackTracker fCrtBackTrack;

    CosmicIdAlg fCosId;

    // Trees
    TTree *fTrackTree;
    TTree *fPfpTree;

    // Track tree parameters
    std::string track_type;
    bool track_pfp_nu;
    int track_nu_tpc;
    int track_pdg;
    double track_time;
    double track_length;
    double track_momentum;
    double track_theta;
    double track_phi;
    bool track_crt_hit_true_match;
    double track_crt_hit_dca;
    bool track_crt_track_true_match;
    double track_crt_track_dca;
    double track_crt_track_angle;
    bool track_stops;
    double track_stop_ratio_start;
    double track_stop_ratio_end;
    double track_fiducial_dist_start;
    double track_fiducial_dist_end;
    int track_tpc;
    bool track_apa_cross;
    double track_apa_dist;
    double track_apa_min_dist;
    double track_pandora_nu_score;

    // PFParticle tree parameters
    std::string pfp_type;
    bool pfp_nu;
    int pfp_nu_tpc;
    int pfp_n_tracks;
    int pfp_pdg;
    double pfp_time;
    double pfp_length;
    double pfp_momentum;
    double pfp_theta;
    double pfp_phi;
    double pfp_tracks_angle;
    double pfp_second_length;
    bool pfp_crt_hit_true_match;
    double pfp_crt_hit_dca;
    double pfp_sec_crt_hit_dca;
    bool pfp_crt_track_true_match;
    double pfp_crt_track_dca;
    double pfp_crt_track_angle;
    bool pfp_stops;
    double pfp_stop_ratio_start;
    double pfp_stop_ratio_end;
    double pfp_sec_stop_ratio_start;
    double pfp_sec_stop_ratio_end;
    double pfp_fiducial_dist_start;
    double pfp_fiducial_dist_end;
    double pfp_sec_fiducial_dist_start;
    double pfp_sec_fiducial_dist_end;
    int pfp_tpc;
    bool pfp_apa_cross;
    double pfp_apa_dist;
    double pfp_apa_min_dist;
    double pfp_sec_apa_min_dist;
    double pfp_pandora_nu_score;

  }; // class CosmicIdTree


  // Constructor
  CosmicIdTree::CosmicIdTree(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTHitLabel         (config().CRTHitLabel())
    , fCRTTrackLabel       (config().CRTTrackLabel())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fCaloModuleLabel     (config().CaloModuleLabel())
    , fPandoraLabel        (config().PandoraLabel())
    , fVerbose             (config().Verbose())
    , fBeamTimeMin         (config().BeamTimeLimits().BeamTimeMin())
    , fBeamTimeMax         (config().BeamTimeLimits().BeamTimeMax())
    , fCrtBackTrack        (config().CrtBackTrack())
    , fCosId               (config().CosIdAlg())
  {

  } //CosmicIdTree()


  void CosmicIdTree::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    // Track tree
    fTrackTree = tfs->make<TTree>("tracks", "tracks");

    fTrackTree->Branch("track_type",                 &track_type);
    fTrackTree->Branch("track_pfp_nu",               &track_pfp_nu, "track_pfp_nu/O");
    fTrackTree->Branch("track_nu_tpc",               &track_nu_tpc, "track_nu_tpc/I");
    fTrackTree->Branch("track_pdg",                  &track_pdg, "track_pdg/I");
    fTrackTree->Branch("track_time",                 &track_time, "track_time/D");
    fTrackTree->Branch("track_length",               &track_length, "track_length/D");
    fTrackTree->Branch("track_momentum",             &track_momentum, "track_momentum/D");
    fTrackTree->Branch("track_theta",                &track_theta, "track_theta/D");
    fTrackTree->Branch("track_phi",                  &track_phi, "track_phi/D");
    fTrackTree->Branch("track_crt_hit_true_match",   &track_crt_hit_true_match, "track_crt_hit_true_match/O");
    fTrackTree->Branch("track_crt_hit_dca",          &track_crt_hit_dca, "track_crt_hit_dca/D");
    fTrackTree->Branch("track_crt_track_true_match", &track_crt_track_true_match, "track_crt_track_true_match/O");
    fTrackTree->Branch("track_crt_track_dca",        &track_crt_track_dca, "track_crt_track_dca/D");
    fTrackTree->Branch("track_crt_track_angle",      &track_crt_track_angle, "track_crt_track_angle/D");
    fTrackTree->Branch("track_stops",                &track_stops, "track_stops/O");
    fTrackTree->Branch("track_stop_ratio_start",     &track_stop_ratio_start, "track_stop_ratio_start/D");
    fTrackTree->Branch("track_stop_ratio_end",       &track_stop_ratio_end, "track_stop_ratio_end/D");
    fTrackTree->Branch("track_fiducial_dist_start",  &track_fiducial_dist_start, "track_fiducial_dist_start/D");
    fTrackTree->Branch("track_fiducial_dist_end",    &track_fiducial_dist_end, "track_fiducial_dist_end/D");
    fTrackTree->Branch("track_tpc",                  &track_tpc, "track_tpc/I");
    fTrackTree->Branch("track_apa_cross",            &track_apa_cross, "track_apa_cross/O");
    fTrackTree->Branch("track_apa_dist",             &track_apa_dist, "track_apa_dist/D");
    fTrackTree->Branch("track_apa_min_dist",         &track_apa_min_dist, "track_apa_min_dist/D");
    fTrackTree->Branch("track_pandora_nu_score",     &track_pandora_nu_score, "track_pandora_nu_score/D");

    // PFParticle tree
    fPfpTree = tfs->make<TTree>("pfps", "pfps");

    fPfpTree->Branch("pfp_type",                 &pfp_type);
    fPfpTree->Branch("pfp_nu",                   &pfp_nu, "pfp_nu/O");
    fPfpTree->Branch("pfp_nu_tpc",               &pfp_nu_tpc, "pfp_nu_tpc/I");
    fPfpTree->Branch("pfp_n_tracks",             &pfp_n_tracks, "pfp_n_tracks/I");
    fPfpTree->Branch("pfp_pdg",                  &pfp_pdg, "pfp_pdg/I");
    fPfpTree->Branch("pfp_time",                 &pfp_time, "pfp_time/D");
    fPfpTree->Branch("pfp_length",               &pfp_length, "pfp_length/D");
    fPfpTree->Branch("pfp_momentum",             &pfp_momentum, "pfp_momentum/D");
    fPfpTree->Branch("pfp_theta",                &pfp_theta, "pfp_theta/D");
    fPfpTree->Branch("pfp_phi",                  &pfp_phi, "pfp_phi/D");
    fPfpTree->Branch("pfp_tracks_angle",         &pfp_tracks_angle, "pfp_tracks_angle/D");
    fPfpTree->Branch("pfp_second_length",        &pfp_second_length, "pfp_second_length/D");
    fPfpTree->Branch("pfp_crt_hit_true_match",   &pfp_crt_hit_true_match, "pfp_crt_hit_true_match/O");
    fPfpTree->Branch("pfp_crt_hit_dca",          &pfp_crt_hit_dca, "pfp_crt_hit_dca/D");
    fPfpTree->Branch("pfp_sec_crt_hit_dca",      &pfp_sec_crt_hit_dca, "pfp_sec_crt_hit_dca/D");
    fPfpTree->Branch("pfp_crt_track_true_match", &pfp_crt_track_true_match, "pfp_crt_track_true_match/O");
    fPfpTree->Branch("pfp_crt_track_dca",        &pfp_crt_track_dca, "pfp_crt_track_dca/D");
    fPfpTree->Branch("pfp_crt_track_angle",      &pfp_crt_track_angle, "pfp_crt_track_angle/D");
    fPfpTree->Branch("pfp_stops",                &pfp_stops, "pfp_stops/O");
    fPfpTree->Branch("pfp_stop_ratio_start",     &pfp_stop_ratio_start, "pfp_stop_ratio_start/D");
    fPfpTree->Branch("pfp_stop_ratio_end",       &pfp_stop_ratio_end, "pfp_stop_ratio_end/D");
    fPfpTree->Branch("pfp_sec_stop_ratio_start", &pfp_sec_stop_ratio_start, "pfp_sec_stop_ratio_start/D");
    fPfpTree->Branch("pfp_sec_stop_ratio_end",   &pfp_sec_stop_ratio_end, "pfp_sec_stop_ratio_end/D");
    fPfpTree->Branch("pfp_fiducial_dist_start",  &pfp_fiducial_dist_start, "pfp_fiducial_dist_start/D");
    fPfpTree->Branch("pfp_fiducial_dist_end",    &pfp_fiducial_dist_end, "pfp_fiducial_dist_end/D");
    fPfpTree->Branch("pfp_sec_fiducial_dist_start", &pfp_sec_fiducial_dist_start, "pfp_sec_fiducial_dist_start/D");
    fPfpTree->Branch("pfp_sec_fiducial_dist_end",  &pfp_sec_fiducial_dist_end, "pfp_sec_fiducial_dist_end/D");
    fPfpTree->Branch("pfp_tpc",                  &pfp_tpc, "pfp_tpc/I");
    fPfpTree->Branch("pfp_apa_cross",            &pfp_apa_cross, "pfp_apa_cross/O");
    fPfpTree->Branch("pfp_apa_dist",             &pfp_apa_dist, "pfp_apa_dist/D");
    fPfpTree->Branch("pfp_apa_min_dist",         &pfp_apa_min_dist, "pfp_apa_min_dist/D");
    fPfpTree->Branch("pfp_sec_apa_min_dist",     &pfp_sec_apa_min_dist, "pfp_sec_apa_min_dist/D");
    fPfpTree->Branch("pfp_pandora_nu_score",     &pfp_pandora_nu_score, "pfp_pandora_nu_score/D");

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Tree Module -------------------"<<std::endl;

  } // CosmicIdTree::beginJob()


  void CosmicIdTree::analyze(const art::Event& event)
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
    // Get g4 particles
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get CRT hits from the event
    art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
    if (event.getByLabel(fCRTHitLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitList, crtHitHandle);

    // Initialize the CRT backtracker for speed
    fCrtBackTrack.Initialize(event);

    // Loop over the CRT hits and match them to true particles
    std::vector<sbn::crt::CRTHit> crtHits;
    std::map<int, int> numHitMap;
    int hit_i = 0;
    for(auto const& hit : (crtHitList)){
      // Don't try to match CRT hits in time with the beam
      double hitTime = hit->ts1_ns * 1e-3;
      if(hitTime > fBeamTimeMin && hitTime < fBeamTimeMax) continue;
      crtHits.push_back(*hit);
      int hitTrueID = fCrtBackTrack.TrueIdFromHitId(event, hit_i);
      hit_i++;
      numHitMap[hitTrueID]++;
    }

    // Get CRT tracks from the event
    art::Handle< std::vector<sbn::crt::CRTTrack>> crtTrackHandle;
    std::vector<art::Ptr<sbn::crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackHandle);

    // Loop over the CRT tracks and match them to true particles
    std::vector<sbn::crt::CRTTrack> crtTracks;
    std::map<int, int> numTrackMap;
    int track_i = 0;
    for(auto const& track : (crtTrackList)){
      // Don't try to match CRT tracks in time with the beam
      double trackTime = track->ts1_ns * 1e-3;
      if(trackTime > fBeamTimeMin && trackTime < fBeamTimeMax) continue;
      crtTracks.push_back(*track);
      int trackTrueID = fCrtBackTrack.TrueIdFromTrackId(event, track_i);
      track_i++;
      numTrackMap[trackTrueID]++;
    }

    // Get reconstructed tracks from the event and hit/calorimetry associations
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    // Get PFParticles from pandora
    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);
    if( !pfParticleHandle.isValid() ){
      if(fVerbose) std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
    // Get PFParticle to track associations
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTPCTrackLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> findManyPFPMetadata(pfParticleHandle,
        event, fPandoraLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    
    // Sort the true particles by type
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;
    std::vector<int> dirtParticleIds;
    std::vector<int> crParticleIds;
    // Record where the beam activity occurs
    int nuTpc = -2;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
      parts.push_back(particle);

      // Get MCTruth
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partID);
      int pdg = std::abs(particle.PdgCode());
      double time = particle.T() * 1e-3; //[us]

      // If origin is a neutrino
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        // If neutrino vertex is not inside the TPC then call it a dirt particle
        if(!fTpcGeo.InFiducial(vtx, 0.)){ 
          dirtParticleIds.push_back(partID);
        }
        // If it's a primary muon
        else if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partID);
        }
        // Other nu particles
        else{
          nuParticleIds.push_back(partID);
        }
      }

      // If origin is a cosmic ray
      else if(truth->Origin() == simb::kCosmicRay){
        crParticleIds.push_back(partID);
        
      }

      // If particle is primary, charged and in time with the beam get the TPC
      if(particle.Mother()==0 && 
         (pdg==13||pdg==111||pdg==211||pdg==2212||pdg==11) && 
         time > fBeamTimeMin && time < fBeamTimeMax){
        std::pair<TVector3, TVector3> cross = fTpcGeo.CrossingPoints(particle);
        if(cross.first.X() != cross.second.X()){
          if(cross.first.X() < 0 && cross.second.X() < 0 && (nuTpc == 0 || nuTpc == -2)) nuTpc = 0;
          else if(cross.first.X() > 0 && cross.second.X() > 0 && (nuTpc == 1 || nuTpc == -2)) nuTpc = 1;
          else nuTpc = -1;
        }
      }
      
    }

    //----------------------------------------------------------------------------------------------------------
    //                                      FAKE PDS RECONSTRUCTION
    //----------------------------------------------------------------------------------------------------------

    // Create fake flashes in each tpc
    std::pair<std::vector<double>, std::vector<double>> fakeFlashes = CosmicIdUtils::FakeTpcFlashes(parts);
    std::vector<double> fakeTpc0Flashes = fakeFlashes.first;
    std::vector<double> fakeTpc1Flashes = fakeFlashes.second;
    bool tpc0BeamFlash = CosmicIdUtils::BeamFlash(fakeTpc0Flashes, fBeamTimeMin, fBeamTimeMax);
    bool tpc1BeamFlash = CosmicIdUtils::BeamFlash(fakeTpc1Flashes, fBeamTimeMin, fBeamTimeMax);

    // If there are no flashes in time with the beam then ignore the event
    if(!tpc0BeamFlash && !tpc1BeamFlash) return;

    //----------------------------------------------------------------------------------------------------------
    //                                FILLING THE PFPARTICLE TREE
    //----------------------------------------------------------------------------------------------------------

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

    //Loop over the pfparticle map
    std::map<int, bool> isPfpNu;
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;

      ResetPfpVars();
      pfp_nu_tpc = nuTpc;

      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

      // FIXME Won't ever look at cosmic pfps as they don't have daughters
      pfp_nu = isNeutrino;

      std::vector<recob::Track> nuTracks;
      // Loop over daughters of pfparticle
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue;

        // Get the first associated track
        recob::Track tpcTrack = *associatedTracks.front();
        nuTracks.push_back(tpcTrack);

        isPfpNu[tpcTrack.ID()] = isNeutrino;

        // Truth match muon tracks and pfps
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ 
          pfp_type = "NuMu";
        }
        else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){ 
          if(pfp_type != "NuMu") pfp_type = "Nu";
        }
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()){ 
          if(pfp_type != "NuMu" && pfp_type != "Nu") pfp_type = "Dirt";
        }
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()){
          if(pfp_type != "NuMu" && pfp_type != "Nu" && pfp_type != "Dirt") pfp_type = "Cr";
        }

        // Get the TPC the pfp was detected in
        std::vector<art::Ptr<recob::Hit>> tpcHits = findManyHits.at(tpcTrack.ID());
        int tpc = fTpcGeo.DetectedInTPC(tpcHits);
        if(tpc == pfp_tpc || pfp_tpc == -99999) pfp_tpc = tpc;
        else if(pfp_tpc != tpc) pfp_tpc = -1;
      }

      // Don't consider PFParticles with no associated tracks
      if(nuTracks.size() == 0) continue;

      pfp_n_tracks = nuTracks.size();

      // Sort tracks by length
      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
                return left.Length() > right.Length();});

      // Choose longest track as cosmic muon candidate
      recob::Track tpcTrack = nuTracks[0];
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());

      // Get truth variables
      if(particles.find(trueId) != particles.end()){
        pfp_pdg = particles[trueId].PdgCode();
        pfp_momentum = particles[trueId].P();
        pfp_time = particles[trueId].T();
        // Does the true particle match any CRT hits or tracks?
        if(numHitMap.find(trueId) != numHitMap.end()){
          if(numHitMap[trueId] > 0) pfp_crt_hit_true_match = true;
        }
        if(numTrackMap.find(trueId) != numTrackMap.end()){
          if(numTrackMap[trueId] > 0) pfp_crt_track_true_match = true;
        }
        // Does particle stop in the TPC?
        geo::Point_t end {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
        if(fTpcGeo.InFiducial(end, 0.)) pfp_stops = true;
        // Does the true particle cross the APA?
        pfp_apa_cross = fTpcGeo.CrossesApa(particles[trueId]);
        // Distance from the APA of the reco track at the true time
        pfp_apa_dist = fCosId.ApaAlg().ApaDistance(detProp, tpcTrack, pfp_time/1e3, hits);
      }

      pfp_length = tpcTrack.Length();
      pfp_theta = tpcTrack.Theta();
      pfp_phi = tpcTrack.Phi();

      // Get the second longest track originating from the same vertex
      recob::Track secTrack = tpcTrack;
      bool useSecTrack = false;
      if(nuTracks.size() > 1){
        TVector3 start = tpcTrack.Vertex<TVector3>();
        TVector3 end = tpcTrack.End<TVector3>();

        for(size_t i = 1; i < nuTracks.size(); i++){
          recob::Track track2 = nuTracks[i];
          TVector3 start2 = track2.Vertex<TVector3>();
          TVector3 end2 = track2.End<TVector3>();
          // Do they share the same vertex? (no delta rays)
          if((start-start2).Mag() > 5.) continue;
          // Get the angle between the two longest tracks
          pfp_tracks_angle = (end - start).Angle(end2 - start2);
          pfp_second_length = track2.Length();
          useSecTrack = true;
          secTrack = track2;
          break;
        }

      }

      // CRT hit cut - get the distance of closest approach for the nearest CRT hit
      std::pair<sbn::crt::CRTHit, double> closestHit = fCosId.CrtHitAlg().T0Alg().ClosestCRTHit(detProp, tpcTrack, crtHits, event);
      pfp_crt_hit_dca = closestHit.second;
      if(useSecTrack){
        std::pair<sbn::crt::CRTHit, double> closestSecHit = fCosId.CrtHitAlg().T0Alg().ClosestCRTHit(detProp, secTrack, crtHits, event);
        pfp_sec_crt_hit_dca = closestHit.second;
      }

      // CRT track cut - get the average distance of closest approach and angle between tracks for the nearest CRT track
      std::pair<sbn::crt::CRTTrack, double> closestTrackDca = fCosId.CrtTrackAlg().TrackAlg().ClosestCRTTrackByDCA(detProp, tpcTrack, crtTracks, event);
      pfp_crt_track_dca = closestTrackDca.second;
      std::pair<sbn::crt::CRTTrack, double> closestTrackAngle = fCosId.CrtTrackAlg().TrackAlg().ClosestCRTTrackByAngle(detProp, tpcTrack, crtTracks, event);
      pfp_crt_track_angle = closestTrackAngle.second;

      // Stopping cut - get the chi2 ratio of the start and end of the track
      pfp_stop_ratio_start = fCosId.StoppingAlg().StoppingChiSq(tpcTrack.Vertex(), calos);
      pfp_stop_ratio_end = fCosId.StoppingAlg().StoppingChiSq(tpcTrack.End(), calos);
      if(useSecTrack){
        std::vector<art::Ptr<anab::Calorimetry>> secCalos = findManyCalo.at(secTrack.ID());
        pfp_sec_stop_ratio_start = fCosId.StoppingAlg().StoppingChiSq(secTrack.Vertex(), secCalos);
        pfp_sec_stop_ratio_end = fCosId.StoppingAlg().StoppingChiSq(secTrack.End(), secCalos);
      }

      // Fiducial cut - Get the fiducial volume the start and end points are contained in
      pfp_fiducial_dist_start = fTpcGeo.MinDistToWall(tpcTrack.Vertex());
      pfp_fiducial_dist_end = fTpcGeo.MinDistToWall(tpcTrack.End());
      if(useSecTrack){
        pfp_sec_fiducial_dist_start = fTpcGeo.MinDistToWall(secTrack.Vertex());
        pfp_sec_fiducial_dist_end = fTpcGeo.MinDistToWall(secTrack.End());
      }

      // APA cut - get the minimum distance to the APA at all PDS times
      std::pair<double, double> ApaMin = fCosId.ApaAlg().MinApaDistance(detProp, tpcTrack, hits, fakeTpc0Flashes, fakeTpc1Flashes);
      pfp_apa_min_dist = ApaMin.first;
      if(useSecTrack){
        std::vector<art::Ptr<recob::Hit>> secHits = findManyHits.at(secTrack.ID());
        std::pair<double, double> ApaMin = fCosId.ApaAlg().MinApaDistance(detProp, secTrack, hits, fakeTpc0Flashes, fakeTpc1Flashes);
        pfp_sec_apa_min_dist = ApaMin.first;
      }

      // Get the PFParticle Nu Score for the PFP Neutrino
      pfp_pandora_nu_score = fCosId.PandoraNuScoreAlg().GetPandoraNuScore(*pParticle, findManyPFPMetadata);

      // Fill the PFParticle tree
      fPfpTree->Fill();

    }

    //----------------------------------------------------------------------------------------------------------
    //                                      FILLING THE TRACK TREE
    //----------------------------------------------------------------------------------------------------------

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){

      ResetTrackVars();
      track_nu_tpc = nuTpc;

      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());

      // Determine the type of the track
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) track_type = "NuMu";
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) track_type = "Nu";
      if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()) track_type = "Cr";
      if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()) track_type = "Dirt";

      // Check if it belongs to a neutrino tagged PFP
      if(isPfpNu.find(tpcTrack.ID()) != isPfpNu.end()){
        if(isPfpNu[tpcTrack.ID()]) track_pfp_nu = true;
      }

      // Get truth variables
      if(particles.find(trueId) != particles.end()){
        track_pdg = particles[trueId].PdgCode();
        track_momentum = particles[trueId].P();
        track_time = particles[trueId].T();
        // Does the true particle match any CRT hits or tracks?
        if(numHitMap.find(trueId) != numHitMap.end()){
          if(numHitMap[trueId] > 0) track_crt_hit_true_match = true;
        }
        if(numTrackMap.find(trueId) != numTrackMap.end()){
          if(numTrackMap[trueId] > 0) track_crt_track_true_match = true;
        }
        // Does particle stop in the TPC?
        geo::Point_t end {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
        if(fTpcGeo.InFiducial(end, 0.)) track_stops = true;
        // Does the true particle cross the APA?
        track_apa_cross = fTpcGeo.CrossesApa(particles[trueId]);
        // Distance from the APA of the reco track at the true time
        track_apa_dist = fCosId.ApaAlg().ApaDistance(detProp, tpcTrack, track_time/1e3, hits);
      }

      track_length = tpcTrack.Length();
      track_theta = tpcTrack.Theta();
      track_phi = tpcTrack.Phi();

      // CRT hit cut - get the distance of closest approach for the nearest CRT hit
      std::pair<sbn::crt::CRTHit, double> closestHit = fCosId.CrtHitAlg().T0Alg().ClosestCRTHit(detProp, tpcTrack, crtHits, event);
      track_crt_hit_dca = closestHit.second;

      // CRT track cut - get the average distance of closest approach and angle between tracks for the nearest CRT track
      std::pair<sbn::crt::CRTTrack, double> closestTrackDca = fCosId.CrtTrackAlg().TrackAlg().ClosestCRTTrackByDCA(detProp, tpcTrack, crtTracks, event);
      track_crt_track_dca = closestTrackDca.second;
      std::pair<sbn::crt::CRTTrack, double> closestTrackAngle = fCosId.CrtTrackAlg().TrackAlg().ClosestCRTTrackByAngle(detProp, tpcTrack, crtTracks, event);
      track_crt_track_angle = closestTrackAngle.second;

      // Stopping cut - get the chi2 ratio of the start and end of the track
      track_stop_ratio_start = fCosId.StoppingAlg().StoppingChiSq(tpcTrack.Vertex(), calos);
      track_stop_ratio_end = fCosId.StoppingAlg().StoppingChiSq(tpcTrack.End(), calos);

      // Fiducial cut - Get the fiducial volume the start and end points are contained in
      track_fiducial_dist_start = fTpcGeo.MinDistToWall(tpcTrack.Vertex());
      track_fiducial_dist_end = fTpcGeo.MinDistToWall(tpcTrack.End());

      // Geometry cut - get the TPC the track was detected in
      track_tpc = fTpcGeo.DetectedInTPC(hits);

      // APA cut - get the minimum distance to the APA at all PDS times
      std::pair<double, double> ApaMin = fCosId.ApaAlg().MinApaDistance(detProp, tpcTrack, hits, fakeTpc0Flashes, fakeTpc1Flashes);
      track_apa_min_dist = ApaMin.first;

      // The PFP Nu Score only exists for PFP Neutrinos
      if (track_pfp_nu){
        for(auto const &pfp : (*pfParticleHandle)){
          // Get the associated track if there is one
          const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pfp.Self()));
          if(associatedTracks.size() != 1) continue;
          recob::Track trk = *associatedTracks.front();
          if(trk.ID() != tpcTrack.ID()) continue;

          recob::PFParticle PFPNeutrino = fCosId.PandoraNuScoreAlg().GetPFPNeutrino(pfp, (*pfParticleHandle));
          track_pandora_nu_score = fCosId.PandoraNuScoreAlg().GetPandoraNuScore(PFPNeutrino, findManyPFPMetadata);
          break;
        }
      }
      // Fill the Track tree
      fTrackTree->Fill();
    }
    
  } // CosmicIdTree::analyze()


  void CosmicIdTree::endJob(){
  
    
  } // CosmicIdTree::endJob()

  void CosmicIdTree::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  void CosmicIdTree::ResetTrackVars(){
    track_type = "none";
    track_pfp_nu = false;
    track_nu_tpc = -99999;
    track_pdg = -99999;
    track_length = -99999;
    track_momentum = -99999;
    track_theta = -99999;
    track_phi = -99999;
    track_crt_hit_true_match = false;
    track_crt_hit_dca = -99999;
    track_crt_track_true_match = false;
    track_crt_track_dca = -99999;
    track_crt_track_angle = -99999;
    track_stops = false;
    track_stop_ratio_start = -99999;
    track_stop_ratio_end = -99999;
    track_fiducial_dist_start = -99999;
    track_fiducial_dist_end = -99999;
    track_tpc = -99999;
    track_apa_cross = false;
    track_apa_dist = -99999;
    track_apa_min_dist = -99999;
    track_pandora_nu_score = -99999;
  }

  void CosmicIdTree::ResetPfpVars(){
    pfp_type = "none";
    pfp_nu = false;
    pfp_nu_tpc = -99999;
    pfp_n_tracks = 0;
    pfp_pdg = -99999;
    pfp_length = -99999;
    pfp_momentum = -99999;
    pfp_theta = -99999;
    pfp_phi = -99999;
    pfp_tracks_angle = -99999;
    pfp_second_length = -99999;
    pfp_crt_hit_true_match = false;
    pfp_crt_hit_dca = -99999;
    pfp_sec_crt_hit_dca = -99999;
    pfp_crt_track_true_match = false;
    pfp_crt_track_dca = -99999;
    pfp_crt_track_angle = -99999;
    pfp_stops = false;
    pfp_stop_ratio_start = -99999;
    pfp_stop_ratio_end = -99999;
    pfp_sec_stop_ratio_start = -99999;
    pfp_sec_stop_ratio_end = -99999;
    pfp_fiducial_dist_start = -99999;
    pfp_fiducial_dist_end = -99999;
    pfp_sec_fiducial_dist_start = -99999;
    pfp_sec_fiducial_dist_end = -99999;
    pfp_tpc = -99999;
    pfp_apa_cross = false;
    pfp_apa_dist = -99999;
    pfp_apa_min_dist = -99999;
    pfp_sec_apa_min_dist = -99999;
    pfp_pandora_nu_score = -99999;
  }
  
  DEFINE_ART_MODULE(CosmicIdTree)
} // namespace sbnd
