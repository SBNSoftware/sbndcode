////////////////////////////////////////////////////////////////////////
// Class:       SelectionTree
// Module Type: analyzer
// File:        SelectionTree_module.cc
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

  class SelectionTree : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> GenModuleLabel {
        Name("GenModuleLabel"),
        Comment("tag of generator data product")
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
    explicit SelectionTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetPfpVars();
    void ResetNuMuVars();

    // Apply the proposal selection
    std::pair<bool, recob::Track> ProposalSelection(std::vector<recob::Track> tracks);
    // Apply Rhiannon's selection
    std::pair<bool, recob::Track> RhiSelection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid);

    void FillSelectionTree(std::string selection, std::pair<bool, recob::Track> selected, int trueId, std::map<int, simb::MCParticle> particles);

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fGenModuleLabel;      ///< name of detsim producer
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

    std::vector<std::string> selections {"prop", "rhi"};
    
    // Tree (One entry per reconstructed pfp)
    TTree *fPfpTree;

    //Pfp tree parameters
    bool is_cosmic;         // True origin of PFP is cosmic
    bool is_dirt;           // True origin of PFP is dirt interaction
    bool is_nu;             // True origin of PFP is nu in AV
    int nu_pdg;             // Pdg of neutrino if not cosmic
    bool is_cc;             // Is interaction CC if not cosmic
    int nu_int;             // Interaction type of neutrino if not cosmic
    double vtx_x;
    double vtx_y;
    double vtx_z;
    bool cosmic_id;         // ID'd as a cosmic
    int n_tracks;           // Number of reconstructed tracks
    double nu_energy;       // Energy of true neutrino
    bool mu_cont;           // Is true muon contained
    double mu_length;       // True contained length of muon if true numuCC
    double mu_mom;          // True momentum of muon if true numuCC
    double mu_theta;        // True theta of muon if true numuCC
    double mu_phi;          // True phi of muon if true numuCC
    std::map<std::string, bool> selected;      // Selected as numuCC?
    std::map<std::string, int> true_pdg;       // PDG of particle prop_true as muon
    std::map<std::string, bool> true_cont;     // Is selected true particle contained 
    std::map<std::string, double> true_length; // True contained length of selected particle
    std::map<std::string, double> true_mom;    // True momentum of selected particle
    std::map<std::string, double> true_theta;  // True theta of selected particle
    std::map<std::string, double> true_phi;    // True phi of selected particle
    std::map<std::string, bool> reco_cont;     // Is reconstructed track contained
    std::map<std::string, double> reco_length; // Selected muon reco length
    std::map<std::string, double> reco_mom;    // Selected muon reco momentum
    std::map<std::string, double> reco_theta;  // Selected muon reco theta
    std::map<std::string, double> reco_phi;    // Selected muon reco phi
    std::map<std::string, double> reco_nu_e;   // Reconstructed neutrino energy assuming nu_mu CC

    // Tree (one entry per numu CC)
    TTree *fNuMuTree;

    // NuMu tree parameters;
    double nu_vtx_x;
    double nu_vtx_y;
    double nu_vtx_z;
    double nu_nu_energy;
    double nu_mu_length;
    double nu_mu_mom;
    double nu_mu_theta;
    double nu_mu_phi;
    bool nu_mu_cont;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

  }; // class SelectionTree


  // Constructor
  SelectionTree::SelectionTree(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fGenModuleLabel       (config().GenModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fPidModuleLabel       (config().PidModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fBeamTimeMin          (config().BeamTimeLimits().BeamTimeMin())
    , fBeamTimeMax          (config().BeamTimeLimits().BeamTimeMax())
    , cosIdAlg              (config().CosIdAlg())
    , fMcsFitter            (config().fitter)
  {

  } // SelectionTree()


  void SelectionTree::beginJob()
  {

    // Maps need initializing
    for(auto const& sel : selections){
      selected[sel] = false;
      true_pdg[sel] = -99999;
      true_cont[sel] = false;
      true_length[sel] = -99999;
      true_mom[sel] = -99999;
      true_theta[sel] = -99999;
      true_phi[sel] = -99999;
      reco_cont[sel] = false;
      reco_length[sel] = -99999;
      reco_mom[sel] = -99999;
      reco_theta[sel] = -99999;
      reco_phi[sel] = -99999;
    }

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fPfpTree = tfs->make<TTree>("pfps", "pfps");

    fPfpTree->Branch("is_cosmic",       &is_cosmic,       "is_cosmic/O");
    fPfpTree->Branch("is_dirt",         &is_dirt,         "is_dirt/O");
    fPfpTree->Branch("is_nu",           &is_nu,           "is_nu/O");
    fPfpTree->Branch("nu_pdg",          &nu_pdg,          "nu_pdg/I");
    fPfpTree->Branch("is_cc",           &is_cc,           "is_cc/O");
    fPfpTree->Branch("nu_int",          &nu_int,          "nu_int/I");
    fPfpTree->Branch("vtx_x",           &vtx_x,           "vtx_x/D");
    fPfpTree->Branch("vtx_y",           &vtx_y,           "vtx_y/D");
    fPfpTree->Branch("vtx_z",           &vtx_z,           "vtx_z/D");
    fPfpTree->Branch("cosmic_id",       &cosmic_id,       "cosmic_id/O");
    fPfpTree->Branch("n_tracks",        &n_tracks,        "n_tracks/I");
    fPfpTree->Branch("nu_energy",       &nu_energy,       "nu_energy/D");
    fPfpTree->Branch("mu_cont",         &mu_cont,         "mu_cont/O");
    fPfpTree->Branch("mu_length",       &mu_length,       "mu_length/D");
    fPfpTree->Branch("mu_mom",          &mu_mom,          "mu_mom/D");
    fPfpTree->Branch("mu_theta",        &mu_theta,        "mu_theta/D");
    fPfpTree->Branch("mu_phi",          &mu_phi,          "mu_phi/D");
    for(auto const& sel : selections){
      fPfpTree->Branch((sel+"_selected").c_str(),   &selected[sel],   (sel+"_selected/O").c_str());
      fPfpTree->Branch((sel+"_true_pdg").c_str(),   &true_pdg[sel],   (sel+"_true_pdg/I").c_str());
      fPfpTree->Branch((sel+"_true_cont").c_str(),  &true_cont[sel],  (sel+"_true_cont/O").c_str());
      fPfpTree->Branch((sel+"_true_length").c_str(),&true_length[sel],(sel+"_true_length/D").c_str());
      fPfpTree->Branch((sel+"_true_mom").c_str(),   &true_mom[sel],   (sel+"_true_mom/D").c_str());
      fPfpTree->Branch((sel+"_true_theta").c_str(), &true_theta[sel], (sel+"_true_theta/D").c_str());
      fPfpTree->Branch((sel+"_true_phi").c_str(),   &true_phi[sel],   (sel+"_true_phi/D").c_str());
      fPfpTree->Branch((sel+"_reco_cont").c_str(),  &reco_cont[sel],  (sel+"_reco_cont/O").c_str());
      fPfpTree->Branch((sel+"_reco_length").c_str(),&reco_length[sel],(sel+"_reco_length/D").c_str());
      fPfpTree->Branch((sel+"_reco_mom").c_str(),   &reco_mom[sel],   (sel+"_reco_mom/D").c_str());
      fPfpTree->Branch((sel+"_reco_theta").c_str(), &reco_theta[sel], (sel+"_reco_theta/D").c_str());
      fPfpTree->Branch((sel+"_reco_phi").c_str(),   &reco_phi[sel],   (sel+"_reco_phi/D").c_str());
      fPfpTree->Branch((sel+"_reco_nu_e").c_str(),  &reco_nu_e[sel],  (sel+"_reco_nu_e/D").c_str());
    }

    fNuMuTree = tfs->make<TTree>("numu", "numu");

    fNuMuTree->Branch("nu_vtx_x", &nu_vtx_x, "nu_vtx_x/D");
    fNuMuTree->Branch("nu_vtx_y", &nu_vtx_y, "nu_vtx_y/D");
    fNuMuTree->Branch("nu_vtx_z", &nu_vtx_z, "nu_vtx_z/D");
    fNuMuTree->Branch("nu_nu_energy", &nu_nu_energy, "nu_nu_energy/D");
    fNuMuTree->Branch("nu_mu_length", &nu_mu_length, "nu_mu_length/D");
    fNuMuTree->Branch("nu_mu_mom", &nu_mu_mom, "nu_mu_mom/D");
    fNuMuTree->Branch("nu_mu_theta", &nu_mu_theta, "nu_mu_theta/D");
    fNuMuTree->Branch("nu_mu_phi", &nu_mu_phi, "nu_mu_phi/D");
    fNuMuTree->Branch("nu_mu_cont", &nu_mu_cont, "nu_mu_cont/O");

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Ana Module -------------------"<<std::endl;

  }// SelectionTree::beginJob()


  void SelectionTree::analyze(const art::Event& event)
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

    art::Handle<std::vector<simb::MCTruth>> genHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthList;
    if(event.getByLabel(fGenModuleLabel, genHandle)) art::fill_ptr_vector(mctruthList, genHandle);


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
    //                                      FILLING THE TRUTH TREE
    //----------------------------------------------------------------------------------------------------------

    for (size_t i = 0; i < mctruthList.size(); i++){

      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> truth = mctruthList.at(i);

      if(truth->Origin() != simb::kBeamNeutrino) continue;

      // Push back all unique neutrino energies
      bool numucc = false;
      if(truth->GetNeutrino().CCNC() == simb::kCC && std::abs(truth->GetNeutrino().Nu().PdgCode()) == 14) numucc = true;

      // Get truth info if numuCC in AV
      if(!numucc) continue;
      geo::Point_t vtx {truth->GetNeutrino().Nu().Vx(), 
                        truth->GetNeutrino().Nu().Vy(), 
                        truth->GetNeutrino().Nu().Vz()};
      if(!fTpcGeo.InFiducial(vtx, 0.)) continue;

      ResetNuMuVars();

      nu_vtx_x = vtx.X();
      nu_vtx_y = vtx.Y();
      nu_vtx_z = vtx.Z();
      nu_nu_energy = truth->GetNeutrino().Nu().E();

      // Get the primary muon
      std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(truth);
      for(auto const& part : parts){
        if(std::abs(part->PdgCode()) != 13) continue;
        if(part->Mother() != 0) continue;
        if(part->StatusCode() != 1) continue;

        nu_mu_length = fTpcGeo.TpcLength(*part);
        nu_mu_mom = part->P();
        TVector3 start(part->Vx(), part->Vy(), part->Vz());
        TVector3 end(part->EndX(), part->EndY(), part->EndZ());
        nu_mu_theta = (end-start).Theta();
        nu_mu_phi = (end-start).Phi();
        nu_mu_cont = fTpcGeo.IsContained(*part);
      }

      fNuMuTree->Fill();
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
    //                                     FILLING THE SELECTION TREE
    //----------------------------------------------------------------------------------------------------------

    //Loop over the pfparticle map
    std::vector<double> used_nus;
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
      // Assign labels based on the particle constributing the most hits
      std::vector<art::Ptr<recob::Hit>> all_hits;
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue; //TODO check how often this occurs

        // Get the first associated track
        recob::Track tpcTrack = *associatedTracks.front();
        nuTracks.push_back(tpcTrack);

        // Truth match the pfps using all hits associated to all tracks associated to neutrino pfp
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        all_hits.insert(all_hits.end(), hits.begin(), hits.end());

      }

      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(all_hits, false);

      // Skip if no corresponding true particle
      if(particles.find(trueId) == particles.end()) continue;

      // Get the origin of the particle
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(trueId);
      if(truth->Origin() == simb::kBeamNeutrino){
        // Save neutrino interaction info
        nu_pdg = truth->GetNeutrino().Nu().PdgCode();
        if(truth->GetNeutrino().CCNC() == simb::kCC) is_cc = true;
        nu_int = truth->GetNeutrino().InteractionType();
        nu_energy = truth->GetNeutrino().Nu().E();

        // Avoid double counting neutrinos
        // FIXME if this ever happens need better way of deciding which pfp to keep
        if(std::find(used_nus.begin(), used_nus.end(), nu_energy) != used_nus.end()) continue;
        used_nus.push_back(nu_energy);

        // If neutrino vertex is not inside the TPC then call it a dirt particle
        geo::Point_t vtx {truth->GetNeutrino().Nu().Vx(), 
                          truth->GetNeutrino().Nu().Vy(), 
                          truth->GetNeutrino().Nu().Vz()};
        vtx_x = vtx.X();
        vtx_y = vtx.Y();
        vtx_z = vtx.Z();
        if(!fTpcGeo.InFiducial(vtx, 0, 0)) is_dirt = true;
        else is_nu = true;

        // If it's a numuCC then save the muon kinematics
        // Get the primary muon
        std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(truth);
        for(auto const& part : parts){
          if(std::abs(part->PdgCode()) != 13) continue;
          if(part->Mother() != 0) continue;
          if(part->StatusCode() != 1) continue;
          TVector3 start(part->Vx(), part->Vy(), part->Vz());
          TVector3 end(part->EndX(), part->EndY(), part->EndZ());
          mu_length = fTpcGeo.TpcLength(*part);
          mu_mom = part->P();
          mu_theta = (end-start).Theta();
          mu_phi = (end-start).Phi();
          mu_cont = fTpcGeo.IsContained(*part);
        }

      }
      else if(truth->Origin() == simb::kCosmicRay) is_cosmic = true;

      n_tracks = nuTracks.size();

      // Skip any PFPs without any tracks in them // TODO check how many numuCC this misses
      if(n_tracks == 0) continue;

      // Does pfp look like a cosmic?
      cosmic_id = cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes);

      // -------------------------------------- APPLY SELECTIONS ---------------------------------------
      for(auto const& sel : selections){
        if(sel == "prop"){
          std::pair<bool, recob::Track> sel_track = ProposalSelection(nuTracks);
          std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(sel_track.second.ID());
          int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
          FillSelectionTree(sel, sel_track, trueId, particles);
        }
        if(sel == "rhi"){
          std::pair<bool, recob::Track> sel_track = RhiSelection(nuTracks, findManyPid);
          std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(sel_track.second.ID());
          int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
          FillSelectionTree(sel, sel_track, trueId, particles);
        } 
      }

      fPfpTree->Fill();

    }  
    
  } // SelectionTree::analyze()


  void SelectionTree::endJob(){

  } // SelectionTree::endJob()

  void SelectionTree::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  // Reset the tree variables
  void SelectionTree::ResetPfpVars(){
    is_cosmic = false;
    is_dirt = false;
    is_nu = false;
    nu_pdg = -99999;
    is_cc = false;
    nu_int = -99999;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    cosmic_id = false;
    n_tracks = 0;
    nu_energy = -99999;
    mu_cont = false;
    mu_length = -99999;
    mu_mom = -99999;
    mu_theta = -99999;
    mu_phi = -99999;
    for(auto const& sel : selections){
      selected[sel] = false;
      true_pdg[sel] = -99999;
      true_cont[sel] = false;
      true_length[sel] = -99999;
      true_mom[sel] = -99999;
      true_theta[sel] = -99999;
      true_phi[sel] = -99999;
      reco_cont[sel] = false;
      reco_length[sel] = -99999;
      reco_mom[sel] = -99999;
      reco_theta[sel] = -99999;
      reco_phi[sel] = -99999;
      reco_nu_e[sel] = -99999;
    }
  }

  void SelectionTree::ResetNuMuVars(){
    nu_vtx_x = -99999;
    nu_vtx_y = -99999;
    nu_vtx_z = -99999;
    nu_nu_energy = -99999;
    nu_mu_mom = -99999;
    nu_mu_theta = -99999;
    nu_mu_phi = -99999;
    nu_mu_cont = false;
  }

  // Apply the proposal selection
  std::pair<bool, recob::Track> SelectionTree::ProposalSelection(std::vector<recob::Track> tracks){

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
  std::pair<bool, recob::Track> SelectionTree::RhiSelection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid){

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
      if(!fTpcGeo.InFiducial(tracks[i].End(), 2.)){ //TODO containment def 
        n_escape++;
        //double length = fTpcGeo.LengthInFiducial(tracks[i], 10, 20, 10, 10, 20, 10);
        double length = tracks[i].Length();
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
    else if(n_escape == 0){
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
    
  void SelectionTree::FillSelectionTree(std::string selection, std::pair<bool, recob::Track> sel_track, int trueId, std::map<int, simb::MCParticle> particles){

    selected[selection] = sel_track.first;

    // Calculate kinematic variables for prop_sel_track prop_track
    reco_cont[selection] = fTpcGeo.InFiducial(sel_track.second.End(), 5.);
    reco_length[selection] = sel_track.second.Length();
    reco_mom[selection] = 0.;
    if(reco_cont[selection]){
      reco_mom[selection] = fRangeFitter.GetTrackMomentum(reco_length[selection], 13);
    }
    else{
      recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(sel_track.second);
      reco_mom[selection] = mcsResult.bestMomentum();
    }
    reco_theta[selection] = sel_track.second.Theta();
    reco_phi[selection] = sel_track.second.Phi();

    // Get the true kinematic variables
    if(particles.find(trueId) != particles.end()){
      TVector3 start(particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz());
      TVector3 end(particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
      true_length[selection] = fTpcGeo.TpcLength(particles[trueId]);
      true_mom[selection] = particles[trueId].P();
      true_theta[selection] = (end-start).Theta();
      true_phi[selection] = (end-start).Phi();
      true_cont[selection] = fTpcGeo.IsContained(particles[trueId]);
      true_pdg[selection] = particles[trueId].PdgCode();
    }

  }

  DEFINE_ART_MODULE(SelectionTree)
} // namespace sbnd

