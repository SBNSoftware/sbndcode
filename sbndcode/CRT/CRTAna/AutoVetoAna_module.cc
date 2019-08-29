////////////////////////////////////////////////////////////////////////
// Class:       AutoVetoAna
// Module Type: analyzer
// File:        AutoVetoAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
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

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace sbnd {

  class AutoVetoAna : public art::EDAnalyzer {
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
      fhicl::Atom<art::InputTag> GenModuleLabel {
        Name("GenModuleLabel"),
        Comment("tag of generator level data product")
      };

      fhicl::Atom<art::InputTag> CosModuleLabel {
        Name("CosModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> G4ModuleLabel {
        Name("G4ModuleLabel"),
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

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

      fhicl::Table<CRTTrackMatchAlg::Config> CrtTrackMatch {
        Name("CrtTrackMatch"),
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit AutoVetoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    void ResetVars();
    void ResetEventVars();
    std::vector<const simb::MCParticle*> InterestingParticles(std::vector<const simb::MCParticle*> particles);

  private:

    // fcl file parameters
    art::InputTag fGenModuleLabel;      ///< name of gen producer
    art::InputTag fCosModuleLabel;      ///< name of detsim producer
    art::InputTag fG4ModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTHitLabel;   ///< name of CRT producer
    art::InputTag fCRTTrackLabel;   ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on
    double fBeamTimeMin;
    double fBeamTimeMax;

    TPCGeoAlg fTpcGeo;

    CRTBackTracker fCrtBackTrack;
    CRTTrackMatchAlg fCrtTrackMatch;

    // Tree
    TTree *fMCTruthTree;

    std::string type;
    double time;
    double energy;
    double vtx_x;
    double vtx_y;
    double vtx_z;
    int num_crt_hits;
    int n_top_high;
    int n_top_low;
    int n_bot;
    int n_side_left;
    int n_side_right;
    int n_face_back;
    int n_face_front;
    int num_crt_tracks;
    int num_crt_through;
    bool hit_cut;
    bool track_cut;
    bool through_cut;
    bool in_tpc;
    double track_x1;
    double track_y1;
    double track_z1;
    double track_x2;
    double track_y2;
    double track_z2;
    int run;
    int subrun;
    int evt;

    TTree *fEventTree;

    bool nu_in_av;
    bool nu_in_fv;
    bool dirt_in_av;
    bool cosmic_in_av;
    bool triggered;
    bool hit_veto;
    int hit_veto_nu;
    int hit_veto_cosmic;
    int hit_veto_dirt;
    int hit_veto_other;
    bool track_veto;
    int track_veto_nu;
    int track_veto_cosmic;
    int track_veto_dirt;
    int track_veto_other;
    bool through_veto;
    int through_veto_nu;
    int through_veto_cosmic;
    int through_veto_dirt;
    int through_veto_other;
    int evt_run;
    int evt_subrun;
    int evt_evt;


  }; // class AutoVetoAna


  // Constructor
  AutoVetoAna::AutoVetoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fGenModuleLabel      (config().GenModuleLabel())
    , fCosModuleLabel      (config().CosModuleLabel())
    , fG4ModuleLabel       (config().G4ModuleLabel())
    , fCRTHitLabel         (config().CRTHitLabel())
    , fCRTTrackLabel       (config().CRTTrackLabel())
    , fVerbose             (config().Verbose())
    , fBeamTimeMin         (config().BeamTimeLimits().BeamTimeMin())
    , fBeamTimeMax         (config().BeamTimeLimits().BeamTimeMax())
    , fCrtBackTrack        (config().CrtBackTrack())
    , fCrtTrackMatch       (config().CrtTrackMatch())
  {

  } //AutoVetoAna()


  void AutoVetoAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    // Tree
    fMCTruthTree = tfs->make<TTree>("mctruth", "mctruth tree");

    fMCTruthTree->Branch("type", &type);
    fMCTruthTree->Branch("time", &time, "time/D");
    fMCTruthTree->Branch("energy", &energy, "energy/D");
    fMCTruthTree->Branch("vtx_x", &vtx_x, "vtx_x/D");
    fMCTruthTree->Branch("vtx_y", &vtx_y, "vtx_y/D");
    fMCTruthTree->Branch("vtx_z", &vtx_z, "vtx_z/D");
    fMCTruthTree->Branch("num_crt_hits", &num_crt_hits, "num_crt_hits/I");
    fMCTruthTree->Branch("n_top_high", &n_top_high, "n_top_high/I");
    fMCTruthTree->Branch("n_top_low", &n_top_low, "n_top_low/I");
    fMCTruthTree->Branch("n_bot", &n_bot, "n_bot/I");
    fMCTruthTree->Branch("n_side_left", &n_side_left, "n_side_left/I");
    fMCTruthTree->Branch("n_side_right", &n_side_right, "n_side_right/I");
    fMCTruthTree->Branch("n_face_front", &n_face_front, "n_face_front/I");
    fMCTruthTree->Branch("n_face_back", &n_face_back, "n_face_back/I");
    fMCTruthTree->Branch("num_crt_tracks", &num_crt_tracks, "num_crt_tracks/I");
    fMCTruthTree->Branch("num_crt_through", &num_crt_through, "num_crt_through/I");
    fMCTruthTree->Branch("hit_cut", &hit_cut, "hit_cut/O");
    fMCTruthTree->Branch("track_cut", &track_cut, "track_cut/O");
    fMCTruthTree->Branch("through_cut", &through_cut, "through_cut/O");
    fMCTruthTree->Branch("in_tpc", &in_tpc, "in_tpc/O");
    fMCTruthTree->Branch("track_x1", &track_x1, "track_x1/D");
    fMCTruthTree->Branch("track_y1", &track_y1, "track_y1/D");
    fMCTruthTree->Branch("track_z1", &track_z1, "track_z1/D");
    fMCTruthTree->Branch("track_x2", &track_x2, "track_x2/D");
    fMCTruthTree->Branch("track_y2", &track_y2, "track_y2/D");
    fMCTruthTree->Branch("track_z2", &track_z2, "track_z2/D");
    fMCTruthTree->Branch("run", &run, "run/I");
    fMCTruthTree->Branch("subrun", &subrun, "subrun/I");
    fMCTruthTree->Branch("evt", &evt, "evt/I");

    fEventTree = tfs->make<TTree>("event", "event tree");

    fEventTree->Branch("nu_in_av", &nu_in_av, "nu_in_av/O");
    fEventTree->Branch("nu_in_fv", &nu_in_fv, "nu_in_fv/O");
    fEventTree->Branch("cosmic_in_av", &cosmic_in_av, "cosmic_in_av/O");
    fEventTree->Branch("dirt_in_av", &dirt_in_av, "dirt_in_av/O");
    fEventTree->Branch("triggered", &triggered, "triggered/O");
    fEventTree->Branch("hit_veto", &hit_veto, "hit_veto/O");
    fEventTree->Branch("hit_veto_nu", &hit_veto_nu, "hit_veto_nu/I");
    fEventTree->Branch("hit_veto_cosmic", &hit_veto_cosmic, "hit_veto_cosmic/I");
    fEventTree->Branch("hit_veto_dirt", &hit_veto_dirt, "hit_veto_dirt/I");
    fEventTree->Branch("hit_veto_other", &hit_veto_other, "hit_veto_other/I");
    fEventTree->Branch("track_veto", &track_veto, "track_veto/O");
    fEventTree->Branch("track_veto_nu", &track_veto_nu, "track_veto_nu/I");
    fEventTree->Branch("track_veto_cosmic", &track_veto_cosmic, "track_veto_cosmic/I");
    fEventTree->Branch("track_veto_dirt", &track_veto_dirt, "track_veto_dirt/I");
    fEventTree->Branch("track_veto_other", &track_veto_other, "track_veto_other/I");
    fEventTree->Branch("through_veto", &through_veto, "through_veto/O");
    fEventTree->Branch("through_veto_nu", &through_veto_nu, "through_veto_nu/I");
    fEventTree->Branch("through_veto_cosmic", &through_veto_cosmic, "through_veto_cosmic/I");
    fEventTree->Branch("through_veto_dirt", &through_veto_dirt, "through_veto_dirt/I");
    fEventTree->Branch("through_veto_other", &through_veto_other, "through_veto_other/I");
    fEventTree->Branch("evt_run", &evt_run, "evt_run/I");
    fEventTree->Branch("evt_subrun", &evt_subrun, "evt_subrun/I");
    fEventTree->Branch("evt_evt", &evt_evt, "evt_evt/I");

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Tree Module -------------------"<<std::endl;

  } // AutoVetoAna::beginJob()


  void AutoVetoAna::analyze(const art::Event& event)
  {

    ResetEventVars();

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }
    evt_run = event.run();
    evt_subrun = event.subRun();
    evt_evt = event.id().event();

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Get genie particles
    art::Handle<std::vector<simb::MCTruth>> genHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthList;
    if(event.getByLabel(fGenModuleLabel, genHandle)) art::fill_ptr_vector(mctruthList, genHandle);
    // Get corsika particles
    art::Handle<std::vector<simb::MCTruth>> cosHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthCosList;
    if(event.getByLabel(fCosModuleLabel, cosHandle)) art::fill_ptr_vector(mctruthCosList, cosHandle);

    // Get all the MCParticles from g4
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fG4ModuleLabel);
    // Make a map between ID and type
    std::map<int, std::string> particleMap;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        if(!fTpcGeo.InFiducial(vtx, 0.)) particleMap[partId] = "Dirt";
        else particleMap[partId] = "Nu";
      }
      else if(truth->Origin() == simb::kCosmicRay) particleMap[partId] = "Cosmic";
    }

    // Get CRT hits from the event
    art::Handle< std::vector<crt::CRTHit>> crtHitHandle;
    std::vector<art::Ptr<crt::CRTHit> > crtHitList;
    if (event.getByLabel(fCRTHitLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitList, crtHitHandle);

    // Initialize the CRT backtracker for speed
    fCrtBackTrack.Initialize(event);

    // Loop over the CRT hits and match them to true particles
    std::vector<std::pair<int, crt::CRTHit>> crtHits;
    std::vector<std::pair<int, crt::CRTHit>> crtHitsBeam;
    int hit_i = 0;
    // FIXME hack for when CRT hits only generated in beam window
    double minHitTime = 99999;
    double maxHitTime = -99999;
    for(auto const& hit : (crtHitList)){
      int hitTrueID = fCrtBackTrack.TrueIdFromHitId(event, hit_i);
      hit_i++;
      crtHits.push_back(std::make_pair(hitTrueID, *hit));
      // If hit's in time with beam record hit and true particle
      double hitTime = hit->ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;
      if(hitTime > fBeamTimeMin && hitTime < fBeamTimeMax && hit->tagger.find("FaceBack")){
        crtHitsBeam.push_back(std::make_pair(hitTrueID, *hit));
      }
    }

    // Get CRT tracks from the event
    art::Handle< std::vector<crt::CRTTrack>> crtTrackHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackHandle);

    // Loop over the CRT tracks and match them to true particles
    std::vector<std::pair<int, crt::CRTTrack>> crtTracks;
    std::vector<std::pair<int, crt::CRTTrack>> crtTracksBeam;
    std::vector<std::pair<int, crt::CRTTrack>> crtThroughTracksBeam;
    int track_i = 0;
    for(auto const& track : (crtTrackList)){
      int trackTrueID = fCrtBackTrack.TrueIdFromTrackId(event, track_i);
      track_i++;
      crtTracks.push_back(std::make_pair(trackTrueID, *track));
      // Don't try to match CRT tracks in time with the beam
      double trackTime = track->ts1_ns * 1e-3;
      if(trackTime > fBeamTimeMin && trackTime < fBeamTimeMax){
        if(!fCrtTrackMatch.CrossesTPC(*track)) continue;
        crtTracksBeam.push_back(std::make_pair(trackTrueID, *track));
        if(track->complete) crtThroughTracksBeam.push_back(std::make_pair(trackTrueID, *track));
      }
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          AUTO VETO CUTS
    //----------------------------------------------------------------------------------------------------------
    if(crtHitsBeam.size() > 0) hit_veto = true;
    if(crtTracksBeam.size() > 0) track_veto = true;
    if(crtThroughTracksBeam.size() > 0) through_veto = true;

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------

    // Loop over MCTruth objects
    for (size_t i = 0; i < mctruthList.size(); i++){
      // Reset all the tree variables
      ResetVars();
      run = evt_run;
      subrun = evt_subrun;
      evt = evt_evt;
      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> mctruth = mctruthList[i];
      // Get type Cosmic/NuMuCC/NuMuNC/NuECC/NuENC/Other
      type = "Other";
      if(mctruth->Origin() == simb::kBeamNeutrino){
        // Get vertex
        geo::Point_t vertex {mctruth->GetNeutrino().Nu().Vx(), 
                             mctruth->GetNeutrino().Nu().Vy(), 
                             mctruth->GetNeutrino().Nu().Vz()};
        vtx_x = vertex.X();
        vtx_y = vertex.Y();
        vtx_z = vertex.Z();
        // Get energy
        energy = mctruth->GetNeutrino().Nu().E();
        // Get time
        time = mctruth->GetNeutrino().Nu().T()*1e-3;
        if(!fTpcGeo.InFiducial(vertex, 0.)){
          type = "Dirt";
        }
        else{
          if(std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 14){
            if(mctruth->GetNeutrino().CCNC() == simb::kCC) type = "NuMuCC";
            else type = "NuMuNC";
          }
          else{
            if(mctruth->GetNeutrino().CCNC() == simb::kCC) type = "NuECC";
            else type = "NuENC";
          }
        }
      }
      if(type == "Other"){
        // Check if it has any particles
        if(mctruth->NParticles() > 0){
          energy = mctruth->GetParticle(0).E();
          time = mctruth->GetParticle(0).T()*1e-3;
        }
      }
      if(time < minHitTime || time > maxHitTime) continue;
      // Get cuts
      hit_cut = hit_veto;
      track_cut = track_veto;
      through_cut = through_veto;
      // Get associated MCParticles and loop over
      std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(mctruth);
      for(size_t j = 0; j < parts.size(); j++){
        int id = parts[j]->TrackId();
        // Get number of CRT hits/tracks
        for(auto const& hit : crtHits){
          if(hit.first == id){ 
            num_crt_hits++;
            if(hit.second.tagger.find("TopHigh") != std::string::npos) n_top_high++;
            if(hit.second.tagger.find("TopLow") != std::string::npos) n_top_low++;
            if(hit.second.tagger.find("Bot") != std::string::npos) n_bot++;
            if(hit.second.tagger.find("SideLeft") != std::string::npos) n_side_left++;
            if(hit.second.tagger.find("SideRight") != std::string::npos) n_side_right++;
            if(hit.second.tagger.find("FaceFront") != std::string::npos) n_face_front++;
            if(hit.second.tagger.find("FaceBack") != std::string::npos) n_face_back++;
          }
        }
        for(auto const& track : crtTracks){
          if(track.first == id){ 
            num_crt_tracks++;
            if(track.second.complete){ 
              num_crt_through++;
              track_x1 = track.second.x1_pos;
              track_y1 = track.second.y1_pos;
              track_z1 = track.second.z1_pos;
              track_x2 = track.second.x2_pos;
              track_y2 = track.second.y2_pos;
              track_z2 = track.second.z2_pos;
            }
          }
        }
      }
      // Get the interesting particles
      std::vector<const simb::MCParticle*> particles = InterestingParticles(parts);
      for(size_t j = 0; j < particles.size(); j++){
        if(fTpcGeo.EntersVolume(*particles[j])) in_tpc = true;
      }
      //FIXME need to investigate nu through going tracks
      fMCTruthTree->Fill();

      if(type.find("Nu") != std::string::npos){ 
        nu_in_av = true;
        if(std::abs(vtx_x) < 183.5 && std::abs(vtx_y) < 185 && vtx_z > 30 && vtx_z < 450) nu_in_fv = true;
      }
      if(type == "Dirt" && in_tpc) dirt_in_av = true;
    }

    // Loop over MCTruth objects
    for (size_t i = 0; i < mctruthCosList.size(); i++){
      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> mctruth = mctruthCosList[i];
      std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(mctruth);
      std::vector<const simb::MCParticle*> particles = InterestingParticles(parts);
      for(size_t j = 0; j < particles.size(); j++){
        // Reset all the tree variables
        ResetVars();
        run = evt_run;
        subrun = evt_subrun;
        evt = evt_evt;
        // Get type Cosmic/NuMuCC/NuMuNC/NuECC/NuENC/Other
        type = "Cosmic";
        // Check if it has any particles
        energy = particles[j]->E();
        time = particles[j]->T()*1e-3;
        if(time < minHitTime || time > maxHitTime) continue;
        // Get cuts
        hit_cut = hit_veto;
        track_cut = track_veto;
        through_cut = through_veto;
        int id = particles[j]->TrackId();
        // Get number of CRT hits/tracks
        for(auto const& hit : crtHits){
          if(hit.first == id){ 
            num_crt_hits++;
            if(hit.second.tagger.find("TopHigh") != std::string::npos) n_top_high++;
            if(hit.second.tagger.find("TopLow") != std::string::npos) n_top_low++;
            if(hit.second.tagger.find("Bot") != std::string::npos) n_bot++;
            if(hit.second.tagger.find("SideLeft") != std::string::npos) n_side_left++;
            if(hit.second.tagger.find("SideRight") != std::string::npos) n_side_right++;
            if(hit.second.tagger.find("FaceFront") != std::string::npos) n_face_front++;
            if(hit.second.tagger.find("FaceBack") != std::string::npos) n_face_back++;
          }
        }
        for(auto const& track : crtTracks){
          if(track.first == id){ 
            num_crt_tracks++;
            if(track.second.complete) num_crt_through++;
          }
        }
        if(fTpcGeo.EntersVolume(*particles[j])) in_tpc = true;
        fMCTruthTree->Fill();

        if(in_tpc) cosmic_in_av = true;
      }
    }

    if(nu_in_av || dirt_in_av || cosmic_in_av) triggered = true;
    for(auto const& hit : crtHitsBeam){
      int id = hit.first;
      if(particleMap.find(id) != particleMap.end()){
        if(particleMap[id] == "Nu") hit_veto_nu++;
        if(particleMap[id] == "Cosmic") hit_veto_cosmic++;
        if(particleMap[id] == "Dirt") hit_veto_dirt++;
      }
      else hit_veto_other++;
    }
    for(auto const& track : crtTracksBeam){
      int id = track.first;
      if(particleMap.find(id) != particleMap.end()){
        if(particleMap[id] == "Nu") track_veto_nu++;
        if(particleMap[id] == "Cosmic") track_veto_cosmic++;
        if(particleMap[id] == "Dirt") track_veto_dirt++;
      }
      else track_veto_other++;
    }
    for(auto const& track : crtThroughTracksBeam){
      int id = track.first;
      if(particleMap.find(id) != particleMap.end()){
        if(particleMap[id] == "Nu") through_veto_nu++;
        if(particleMap[id] == "Cosmic") through_veto_cosmic++;
        if(particleMap[id] == "Dirt") through_veto_dirt++;
      }
      else through_veto_other++;
    }

    fEventTree->Fill();

  } // AutoVetoAna::analyze()


  void AutoVetoAna::endJob(){
  
    
  } // AutoVetoAna::endJob()

  // Reset variables and counters
  void AutoVetoAna::ResetVars(){
    type = "Other";
    time = -99999;
    energy = -99999;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    num_crt_hits = 0;
    n_top_high = 0;
    n_top_low = 0;
    n_bot = 0;
    n_side_left = 0;
    n_side_right = 0;
    n_face_front = 0;
    n_face_back = 0;
    num_crt_tracks = 0;
    num_crt_through = 0;
    hit_cut = false;
    track_cut = false;
    through_cut = false;
    in_tpc = false;
    track_x1 = -99999;
    track_y1 = -99999;
    track_z1 = -99999;
    track_x2 = -99999;
    track_y2 = -99999;
    track_z2 = -99999;
    run = -99999;
    subrun = -99999;
    evt = -99999;
  }

  void AutoVetoAna::ResetEventVars(){
    nu_in_av = false;
    nu_in_fv = false;
    dirt_in_av = false;
    cosmic_in_av = false;
    triggered = false;
    hit_veto = false;
    hit_veto_nu = 0;
    hit_veto_cosmic = 0;
    hit_veto_dirt = 0;
    hit_veto_other = 0;
    track_veto = false;
    track_veto_nu = 0;
    track_veto_cosmic = 0;
    track_veto_dirt = 0;
    track_veto_other = 0;
    through_veto = false;
    through_veto_nu = 0;
    through_veto_cosmic = 0;
    through_veto_dirt = 0;
    through_veto_other = 0;
    evt_run = -99999;
    evt_subrun = -99999;
    evt_evt = -99999;
  }

  // Give us a list of the stable, primary particles that we're interested in
  std::vector<const simb::MCParticle*> AutoVetoAna::InterestingParticles(std::vector<const simb::MCParticle*> particles){

    std::vector<const simb::MCParticle*> interesting;

    // Loop over all of the particles
    for(size_t j = 0; j < particles.size(); j++){
      // Only consider stable final states particles
      if(particles[j]->StatusCode() != 1) continue;
      // Only want primary particles
      if(particles[j]->Mother() != 0) continue;

      // Only consider muons, pi0, charged pi and protons TODO for now...
      int pdg = std::abs(particles[j]->PdgCode());
      if(!(pdg == 13 || pdg == 111 || pdg == 211 || pdg == 2212 || pdg == 11)) continue;
      interesting.push_back(particles[j]);
    }

    return interesting;

  } // AutoVetoAna::InterestingParticles()

  DEFINE_ART_MODULE(AutoVetoAna)
} // namespace sbnd


