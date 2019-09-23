////////////////////////////////////////////////////////////////////////
// Class:       SelectionComparison
// Module Type: analyzer
// File:        SelectionComparison_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/CosmicId/Algs/CosmicIdAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
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

  class SelectionComparison : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<art::InputTag> PidModuleLabel {
        Name("PidModuleLabel"),
        Comment("tag of PID producer data product")
      };

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CosmicIdAlg::Config> CosIdAlg {
        Name("CosIdAlg"),
      };

      fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
        Name("fitter"),
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit SelectionComparison(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetPfpVars();

    // Apply the proposal selection
    std::pair<bool, recob::Track> ProposalSelection(std::vector<recob::Track> tracks);
    // Apply Rhiannon's selection
    std::pair<bool, recob::Track> RhiSelection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid);

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fPidModuleLabel; ///< name of TPC track producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on
    double        fBeamTimeMin;
    double        fBeamTimeMax;

    CosmicIdAlg cosIdAlg;
    TPCGeoAlg fTpcGeo;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    
    // Tree (One entry per reconstructed pfp)
    TTree *fPfpTree;

    //Pfp tree parameters
    bool is_cosmic;         // True origin of PFP is cosmic
    bool is_dirt;           // True origin of PFP is dirt interaction
    bool is_nu;             // True origin of PFP is nu in AV
    int nu_pdg;             // Pdg of neutrino if not cosmic
    bool is_cc;             // Is interaction CC if not cosmic
    int nu_int;             // Interaction type of neutrino if not cosmic
    bool cosmic_id;         // ID'd as a cosmic
    int n_tracks;           // Number of reconstructed tracks
    double nu_energy;       // Energy of true neutrino
    bool mu_cont;           // Is true muon contained
    double mu_length;       // True contained length of muon if true numuCC
    double mu_mom;          // True momentum of muon if true numuCC
    double mu_theta;        // True theta of muon if true numuCC
    double mu_phi;          // True phi of muon if true numuCC
    bool prop_in_fv;             // True vertex in fiducial volume
    bool prop_selected;          // Selected as numuCC
    int prop_selected_pdg;       // PDG of particle prop_selected as muon
    bool prop_selected_cont;     // Is true prop_selected particle contained 
    double prop_selected_length; // True contained length of prop_selected particle
    double prop_selected_mom;    // True momentum of prop_selected particle
    double prop_selected_theta;  // True theta of prop_selected particle
    double prop_selected_phi;    // True phi of prop_selected particle
    bool prop_track_cont;        // Is reconstructed prop_track contained
    double prop_track_length;    // Selected muon prop_track length
    double prop_track_mom;       // Selected muon prop_track momentum
    double prop_track_theta;     // Selected muon prop_track theta
    double prop_track_phi;       // Selected muon prop_track phi
    bool rhi_in_fv;             // True vertex in fiducial volume
    bool rhi_selected;          // Selected as numuCC
    int rhi_selected_pdg;       // PDG of particle rhi_selected as muon
    bool rhi_selected_cont;     // Is true rhi_selected particle contained 
    double rhi_selected_length; // True contained length of rhi_selected particle
    double rhi_selected_mom;    // True momentum of rhi_selected particle
    double rhi_selected_theta;  // True theta of rhi_selected particle
    double rhi_selected_phi;    // True phi of rhi_selected particle
    bool rhi_track_cont;        // Is reconstructed rhi_track contained
    double rhi_track_length;    // Selected muon rhi_track length
    double rhi_track_mom;       // Selected muon rhi_track momentum
    double rhi_track_theta;     // Selected muon rhi_track theta
    double rhi_track_phi;       // Selected muon rhi_track phi

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

  }; // class SelectionComparison


  // Constructor
  SelectionComparison::SelectionComparison(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fPidModuleLabel       (config().PidModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fBeamTimeMin          (config().BeamTimeLimits().BeamTimeMin())
    , fBeamTimeMax          (config().BeamTimeLimits().BeamTimeMax())
    , cosIdAlg              (config().CosIdAlg())
    , fMcsFitter            (config().fitter)
  {

  } // SelectionComparison()


  void SelectionComparison::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fPfpTree = tfs->make<TTree>("pfps", "pfps");

    fPfpTree->Branch("is_cosmic",       &is_cosmic,       "is_cosmic/O");
    fPfpTree->Branch("is_dirt",         &is_dirt,         "is_dirt/O");
    fPfpTree->Branch("is_nu",           &is_nu,           "is_nu/O");
    fPfpTree->Branch("nu_pdg",          &nu_pdg,          "nu_pdg/I");
    fPfpTree->Branch("is_cc",           &is_cc,           "is_cc/O");
    fPfpTree->Branch("nu_int",          &nu_int,          "nu_int/I");
    fPfpTree->Branch("cosmic_id",       &cosmic_id,       "cosmic_id/O");
    fPfpTree->Branch("n_tracks",        &n_tracks,        "n_tracks/I");
    fPfpTree->Branch("nu_energy",       &nu_energy,       "nu_energy/D");
    fPfpTree->Branch("mu_cont",         &mu_cont,         "mu_cont/O");
    fPfpTree->Branch("mu_length",       &mu_length,       "mu_length/D");
    fPfpTree->Branch("mu_mom",          &mu_mom,          "mu_mom/D");
    fPfpTree->Branch("mu_theta",        &mu_theta,        "mu_theta/D");
    fPfpTree->Branch("mu_phi",          &mu_phi,          "mu_phi/D");
    fPfpTree->Branch("prop_in_fv",           &prop_in_fv,           "prop_in_fv/O");
    fPfpTree->Branch("prop_selected",        &prop_selected,        "prop_selected/O");
    fPfpTree->Branch("prop_selected_pdg",    &prop_selected_pdg,    "prop_selected_pdg/I");
    fPfpTree->Branch("prop_selected_cont",   &prop_selected_cont,   "prop_selected_cont/O");
    fPfpTree->Branch("prop_selected_length", &prop_selected_length, "prop_selected_length/D");
    fPfpTree->Branch("prop_selected_mom",    &prop_selected_mom,    "prop_selected_mom/D");
    fPfpTree->Branch("prop_selected_theta",  &prop_selected_theta,  "prop_selected_theta/D");
    fPfpTree->Branch("prop_selected_phi",    &prop_selected_phi,    "prop_selected_phi/D");
    fPfpTree->Branch("prop_track_cont",      &prop_track_cont,      "prop_track_cont/O");
    fPfpTree->Branch("prop_track_length",    &prop_track_length,    "prop_track_length/D");
    fPfpTree->Branch("prop_track_mom",       &prop_track_mom,       "prop_track_mom/D");
    fPfpTree->Branch("prop_track_theta",     &prop_track_theta,     "prop_track_theta/D");
    fPfpTree->Branch("prop_track_phi",       &prop_track_phi,       "prop_track_phi/D");
    fPfpTree->Branch("rhi_in_fv",           &rhi_in_fv,           "rhi_in_fv/O");
    fPfpTree->Branch("rhi_selected",        &rhi_selected,        "rhi_selected/O");
    fPfpTree->Branch("rhi_selected_pdg",    &rhi_selected_pdg,    "rhi_selected_pdg/I");
    fPfpTree->Branch("rhi_selected_cont",   &rhi_selected_cont,   "rhi_selected_cont/O");
    fPfpTree->Branch("rhi_selected_length", &rhi_selected_length, "rhi_selected_length/D");
    fPfpTree->Branch("rhi_selected_mom",    &rhi_selected_mom,    "rhi_selected_mom/D");
    fPfpTree->Branch("rhi_selected_theta",  &rhi_selected_theta,  "rhi_selected_theta/D");
    fPfpTree->Branch("rhi_selected_phi",    &rhi_selected_phi,    "rhi_selected_phi/D");
    fPfpTree->Branch("rhi_track_cont",      &rhi_track_cont,      "rhi_track_cont/O");
    fPfpTree->Branch("rhi_track_length",    &rhi_track_length,    "rhi_track_length/D");
    fPfpTree->Branch("rhi_track_mom",       &rhi_track_mom,       "rhi_track_mom/D");
    fPfpTree->Branch("rhi_track_theta",     &rhi_track_theta,     "rhi_track_theta/D");
    fPfpTree->Branch("rhi_track_phi",       &rhi_track_phi,       "rhi_track_phi/D");

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Ana Module -------------------"<<std::endl;

  }// SelectionComparison::beginJob()


  void SelectionComparison::analyze(const art::Event& event)
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
    // Put them in a map for easier access
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    for (auto const& particle: (*particleHandle)){
      // Store particle
      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
    }

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
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);
    
    // Get track to hit and colorimetry associations
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindMany<anab::ParticleID> findManyPid(tpcTrackHandle, event, fPidModuleLabel);

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
    //                                        RUNNING SELECTIONS
    //----------------------------------------------------------------------------------------------------------

    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || 
                            std::abs(pdg) == pandora::NU_MU || 
                            std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;

      ResetPfpVars();

      std::vector<recob::Track> nuTracks;

      // Loop over daughters of pfparticle and do some truth matching
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue; //TODO check how often this occurs

        // Get the first associated track
        recob::Track tpcTrack = *associatedTracks.front();
        nuTracks.push_back(tpcTrack);

        // Truth match muon tracks and pfps
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);

        // Skip if no corresponding true particle
        if(particles.find(trueId) == particles.end()) continue;
        int pdg = std::abs(particles[trueId].PdgCode());

        // Get the origin of the particle
        art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(trueId);
        if(truth->Origin() == simb::kBeamNeutrino){
          // Save neutrino interaction info
          nu_pdg = truth->GetNeutrino().Nu().PdgCode();
          if(truth->GetNeutrino().CCNC() == simb::kCC) is_cc = true;
          nu_int = truth->GetNeutrino().InteractionType();
          nu_energy = truth->GetNeutrino().Nu().E();

          // If neutrino vertex is not inside the TPC then call it a dirt particle
          geo::Point_t vtx {truth->GetNeutrino().Nu().Vx(), 
                            truth->GetNeutrino().Nu().Vy(), 
                            truth->GetNeutrino().Nu().Vz()};
          if(!fTpcGeo.InFiducial(vtx, 0, 0)) is_dirt = true;
          else is_nu = true;

          // If it's a numuCC then save the muon kinematics
          if(std::abs(nu_pdg) == 14 && is_cc && pdg == 13 && particles[trueId].Mother() == 0){
            TVector3 start(particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz());
            TVector3 end(particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
            mu_length = fTpcGeo.TpcLength(particles[trueId]);
            mu_mom = particles[trueId].P();
            mu_theta = (end-start).Theta();
            mu_phi = (end-start).Phi();
            mu_cont = fTpcGeo.IsContained(particles[trueId]);
          }

          // Selection specific FV definitions
          if(fTpcGeo.InFiducial(vtx, 16.5, 15., 15., 16.5, 15., 80.)) prop_in_fv = true;
          if(fTpcGeo.InFiducial(vtx, 8.25, 15., 15., 8.25, 15., 85., 8.25)) rhi_in_fv = true;
        }
        else if(truth->Origin() == simb::kCosmicRay) is_cosmic = true;
      }

      n_tracks = nuTracks.size();

      // Skip any PFPs without any tracks in them // TODO check how many numuCC this misses
      if(n_tracks == 0) continue;

      // Does pfp look like a cosmic?
      cosmic_id = cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes);

      // -------------------------------------- PROPOSAL SELECTION ---------------------------------------
      // Apply the proposal selection
      std::pair<bool, recob::Track> proposal = ProposalSelection(nuTracks);
      prop_selected = proposal.first;

      // Calculate kinematic variables for prop_selected prop_track
      prop_track_cont = fTpcGeo.InFiducial(proposal.second.End(), 5.);
      prop_track_length = proposal.second.Length();
      prop_track_mom = 0.;
      if(prop_track_cont){
        prop_track_mom = fRangeFitter.GetTrackMomentum(prop_track_length, 13);
      }
      else{
        recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(proposal.second);
        prop_track_mom = mcsResult.bestMomentum();
      }
      prop_track_theta = proposal.second.Theta();
      prop_track_phi = proposal.second.Phi();

      // Get the true prop_selected particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(proposal.second.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      // Get the true kinematic variables
      if(particles.find(trueId) != particles.end()){
        TVector3 start(particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz());
        TVector3 end(particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
        prop_selected_length = fTpcGeo.TpcLength(particles[trueId]);
        prop_selected_mom = particles[trueId].P();
        prop_selected_theta = (end-start).Theta();
        prop_selected_phi = (end-start).Phi();
        prop_selected_cont = fTpcGeo.IsContained(particles[trueId]);
        prop_selected_pdg = particles[trueId].PdgCode();
      }

      // -------------------------------------- RHIANNON'S SELECTION ---------------------------------------
      // Apply the proposal selection
      std::pair<bool, recob::Track> rhiannon = RhiSelection(nuTracks, findManyPid);
      rhi_selected = rhiannon.first;

      // Calculate kinematic variables for rhi_selected rhi_track
      rhi_track_cont = fTpcGeo.InFiducial(rhiannon.second.End(), 5.);
      rhi_track_length = rhiannon.second.Length();
      rhi_track_mom = 0.;
      if(rhi_track_cont){
        rhi_track_mom = fRangeFitter.GetTrackMomentum(rhi_track_length, 13);
      }
      else{
        recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(rhiannon.second);
        rhi_track_mom = mcsResult.bestMomentum();
      }
      rhi_track_theta = rhiannon.second.Theta();
      rhi_track_phi = rhiannon.second.Phi();

      // Get the true rhi_selected particle
      std::vector<art::Ptr<recob::Hit>> rhi_hits = findManyHits.at(rhiannon.second.ID());
      int rhi_trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(rhi_hits, false);
      // Get the true kinematic variables
      if(particles.find(rhi_trueId) != particles.end()){
        TVector3 start(particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz());
        TVector3 end(particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
        rhi_selected_length = fTpcGeo.TpcLength(particles[trueId]);
        rhi_selected_mom = particles[trueId].P();
        rhi_selected_theta = (end-start).Theta();
        rhi_selected_phi = (end-start).Phi();
        rhi_selected_cont = fTpcGeo.IsContained(particles[trueId]);
        rhi_selected_pdg = particles[trueId].PdgCode();
      }

      fPfpTree->Fill();

    }  
    
  } // SelectionComparison::analyze()


  void SelectionComparison::endJob(){

  } // SelectionComparison::endJob()

  void SelectionComparison::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  // Reset the tree variables
  void SelectionComparison::ResetPfpVars(){
    is_cosmic = false;
    is_dirt = false;
    is_nu = false;
    nu_pdg = -99999;
    is_cc = false;
    nu_int = -99999;
    cosmic_id = false;
    n_tracks = 0;
    nu_energy = -99999;
    mu_cont = false;
    mu_length = -99999;
    mu_mom = -99999;
    mu_theta = -99999;
    mu_phi = -99999;
    prop_in_fv = false;
    prop_selected = false;
    prop_selected_pdg = -99999;
    prop_selected_cont = false;
    prop_selected_length = -99999;
    prop_selected_mom = -99999;
    prop_selected_theta = -99999;
    prop_selected_phi = -99999;
    prop_track_cont = false;
    prop_track_length = -99999;
    prop_track_mom = -99999;
    prop_track_theta = -99999;
    prop_track_phi = -99999;
    rhi_in_fv = false;
    rhi_selected = false;
    rhi_selected_pdg = -99999;
    rhi_selected_cont = false;
    rhi_selected_length = -99999;
    rhi_selected_mom = -99999;
    rhi_selected_theta = -99999;
    rhi_selected_phi = -99999;
    rhi_track_cont = false;
    rhi_track_length = -99999;
    rhi_track_mom = -99999;
    rhi_track_theta = -99999;
    rhi_track_phi = -99999;
  }
    
  // Apply the proposal selection
  std::pair<bool, recob::Track> SelectionComparison::ProposalSelection(std::vector<recob::Track> tracks){

    bool is_selected = false;

    // Get the longest track and ID as muon
    std::sort(tracks.begin(), tracks.end(), [](auto& left, auto& right){
              return left.Length() > right.Length();});
    recob::Track track = tracks[0];
    double length = track.Length();

    // Check if the track is contained
    bool track_contained = fTpcGeo.InFiducial(track.End(), 5.);

    // Apply a fiducial volume cut to the vertex (start of track) TODO CPA cut
    bool vertex_contained = fTpcGeo.InFiducial(track.Start(), 16.5, 15., 15., 16.5, 15., 80.);

    // Apply selection criteria
    if(track_contained){
      if(vertex_contained && length > 50.) is_selected = true;
    }
    else{
      if(vertex_contained && length > 100.) is_selected = true;
    }

    return std::make_pair(is_selected, track);
  }

  // Apply rhiannon's selection
  std::pair<bool, recob::Track> SelectionComparison::RhiSelection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid){

    bool is_selected = false;
    bool has_candidate = false;
    recob::Track candidate = tracks[0];

    // Loop over tracks and count how many escape
    int n_escape = 0;
    double longest_escape = 0;
    int longest_i = -99999;
    double longest_first = 0;
    double longest_second = 0;
    for(size_t i = 0; i < tracks.size(); i++){
      if(!fTpcGeo.InFiducial(tracks[i].End(), 0.)){ //TODO containment def 
        n_escape++;
        double length = fTpcGeo.LengthInFiducial(tracks[i], 10, 20, 10, 10, 20, 10);
        if(length > longest_escape){ 
          longest_escape = length;
          longest_i = i;
        }
      }

      if(tracks[i].Length() > longest_first){
        longest_second = longest_first;
        longest_first = tracks[i].Length();
      }
      else if(tracks[i].Length() > longest_second){
        longest_second = tracks[i].Length();
      }
    }

    // Case 1: 1 escaping track
    if(n_escape == 1){
      // If track longer than 100 cm then ID as muon
      // If more than 1 escaping track > 100 cm then choose longest
      if(longest_escape > 100){
        has_candidate = true;
        candidate = tracks[longest_i];
      }
      // Else don't select ? TODO what if long contained track but short exiting one
    }

    // Case 2: everything else
    else{
      std::vector<std::pair<recob::Track, double>> candidates;
      for(auto const& track : tracks){
        std::vector<const anab::ParticleID*> pids = fmpid.at(track.ID());
        bool is_proton = false;
        double muon_chi = 99999;
        for(size_t i = 0; i < pids.size(); i++){
          // Only use the collection plane
          if(pids[i]->PlaneID().Plane != 2) continue;
          // Run proton chi^2 to cut out stopping protons
          if(pids[i]->Chi2Proton() < 80){ 
            is_proton = true;
            continue;
          }
          // Muon candidates = longest by 1.5x and any with low muon chi^2
          muon_chi = pids[i]->Chi2Muon();
        }
        if(!is_proton && (muon_chi<16 || track.Length() >= 1.5*longest_second)) 
          candidates.push_back(std::make_pair(track, muon_chi));
      }
      // 0 candidates = don't select
      // 1 candidate = select as muon
      if(candidates.size() == 1){
        has_candidate = true;
        candidate = candidates[0].first;
      }
      // >1 candidate = Select longest if it was a candidate or track with lowest muon chi^2
      else if(candidates.size() > 1){
        // Sort by length
        std::sort(candidates.begin(), candidates.end(), [](auto& left, auto& right){
              return left.first.Length() > right.first.Length();});
        if(candidates[0].first.Length() >= 1.5*longest_second){
          has_candidate = true;
          candidate = candidates[0].first;
        }
        else{
          // Sort by muon chi2
          std::sort(candidates.begin(), candidates.end(), [](auto& left, auto& right){
                return left.second < right.second;});
          has_candidate = true;
          candidate = candidates[0].first;
        }
      }
    }

    // Check vertex (start of muon candidate) in FV
    // Fiducial definitions 8.25 cm from X (inc CPA), 15 cm from Y and front, 85 cm from back
    bool vertex_contained = fTpcGeo.InFiducial(candidate.Start(), 8.25, 15., 15., 8.25, 15., 85., 8.25);
    if(vertex_contained && has_candidate) is_selected = true;

    return std::make_pair(is_selected, candidate);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  DEFINE_ART_MODULE(SelectionComparison)
} // namespace sbnd

