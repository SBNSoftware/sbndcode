////////////////////////////////////////////////////////////////////////
// Class:       NumuccTree
// Module Type: analyzer
// File:        NumuccTree_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/CosmicId/Algs/CosmicIdAlg.h"
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"

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

  class NumuccTree : public art::EDAnalyzer {
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

      fhicl::Table<StoppingParticleCosmicIdAlg::Config> SPTagAlg {
        Name("SPTagAlg"),
      };

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit NumuccTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called for every sub run
    virtual void beginSubRun(const art::SubRun& subrun) override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetNuMuVars();

    std::pair<std::pair<bool, bool>, recob::Track> Selection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid, art::FindManyP<anab::Calorimetry> fmcalo);

    double AverageDCA(const recob::Track& track);

    double HitEnergy(art::Ptr<recob::Hit> hit);
    double NeutrinoEnergy(std::map<int, double> track_energies, double shower_energy, int id);

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fGenModuleLabel;      ///< name of detsim producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fShowerModuleLabel; ///< name of TPC track producer
    art::InputTag fPidModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of TPC track producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on
    double        fBeamTimeMin;
    double        fBeamTimeMax;

    CosmicIdAlg cosIdAlg;
    TPCGeoAlg fTpcGeo;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    StoppingParticleCosmicIdAlg  fStopTagger;
    detinfo::DetectorProperties const* fDetectorProperties;

    // Tree (one entry per numu CC)
    TTree *fNuMuTree;

    // NuMu tree parameters;
    int nu_pdg;             // Pdg of neutrino if not cosmic
    bool is_cc;             // Is interaction CC if not cosmic
    int nu_int;             // Interaction type of neutrino if not cosmic
    double true_vtx_x;
    double true_vtx_y;
    double true_vtx_z;
    double true_nu_energy;
    double true_length;
    double true_mom;
    double true_theta;
    double true_phi;
    int true_n_tracks;
    bool true_cont;
    bool reconstructed;
    bool cosmic_id;
    bool in_fv;
    bool selected;
    double reco_vtx_x;      // Reconstructed neutrino vertex X
    double reco_vtx_y;      // Reconstructed neutrino vertex Y
    double reco_vtx_z;      // Reconstructed neutrino vertex Z
    double reco_nu_energy;  // Reconstructed neutrino energy assuming nu_mu CC
    double reco_length;     // Selected muon reco length
    double reco_mom;        // Selected muon reco momentum
    double reco_theta;      // Selected muon reco theta
    double reco_phi;        // Selected muon reco phi
    int reco_n_tracks;
    bool reco_cont;         // Is reconstructed track contained

    TTree *fPotTree;
    double pot;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

  }; // class NumuccTree


  // Constructor
  NumuccTree::NumuccTree(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fGenModuleLabel       (config().GenModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fShowerModuleLabel    (config().ShowerModuleLabel())
    , fPidModuleLabel       (config().PidModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , fBeamTimeMin          (config().BeamTimeLimits().BeamTimeMin())
    , fBeamTimeMax          (config().BeamTimeLimits().BeamTimeMax())
    , cosIdAlg              (config().CosIdAlg())
    , fMcsFitter            (config().fitter)
    , fStopTagger           (config().SPTagAlg())
  {

  } // NumuccTree()


  void NumuccTree::beginJob()
  {

    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fNuMuTree = tfs->make<TTree>("numu", "numu");

    fNuMuTree->Branch("nu_pdg",    &nu_pdg);
    fNuMuTree->Branch("is_cc",     &is_cc);
    fNuMuTree->Branch("nu_int",    &nu_int);
    fNuMuTree->Branch("true_vtx_x", &true_vtx_x);
    fNuMuTree->Branch("true_vtx_y", &true_vtx_y);
    fNuMuTree->Branch("true_vtx_z", &true_vtx_z);
    fNuMuTree->Branch("true_nu_energy", &true_nu_energy);
    fNuMuTree->Branch("true_length", &true_length);
    fNuMuTree->Branch("true_mom", &true_mom);
    fNuMuTree->Branch("true_theta", &true_theta);
    fNuMuTree->Branch("true_phi", &true_phi);
    fNuMuTree->Branch("true_cont", &true_cont);
    fNuMuTree->Branch("true_n_tracks", &true_n_tracks);
    fNuMuTree->Branch("reconstructed", &reconstructed);
    fNuMuTree->Branch("cosmic_id", &cosmic_id);
    fNuMuTree->Branch("in_fv", &in_fv);
    fNuMuTree->Branch("selected", &selected);
    fNuMuTree->Branch("reco_vtx_x", &reco_vtx_x);
    fNuMuTree->Branch("reco_vtx_y", &reco_vtx_y);
    fNuMuTree->Branch("reco_vtx_z", &reco_vtx_z);
    fNuMuTree->Branch("reco_nu_energy", &reco_nu_energy);
    fNuMuTree->Branch("reco_length", &reco_length);
    fNuMuTree->Branch("reco_mom", &reco_mom);
    fNuMuTree->Branch("reco_theta", &reco_theta);
    fNuMuTree->Branch("reco_phi", &reco_phi);
    fNuMuTree->Branch("reco_cont", &reco_cont);
    fNuMuTree->Branch("reco_n_tracks", &reco_n_tracks);

    fPotTree = tfs->make<TTree>("pots", "pots");
    fPotTree->Branch("pot",  &pot);

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Ana Module -------------------"<<std::endl;

  }// NumuccTree::beginJob()

  // Called for every sub run
  void NumuccTree::beginSubRun(const art::SubRun& subrun){

    art::Handle< sumdata::POTSummary > potHandle;
    subrun.getByLabel(fGenModuleLabel, potHandle);
    const sumdata::POTSummary& potSum = (*potHandle);
    pot = potSum.totpot;

    fPotTree->Fill();

    return;
  }

  void NumuccTree::analyze(const art::Event& event)
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
    for (auto const& particle: (*particleHandle)){
      // Store particle
      int partId = particle.TrackId();
      particles[partId] = particle;
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
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    // Get shower handle
    auto showerHandle = event.getValidHandle<std::vector<recob::Shower>>(fShowerModuleLabel);
    // Get PFParticle to shower associations
    art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, event, fShowerModuleLabel);
    // Get shower to hit associations
    art::FindManyP<recob::Hit> findManyHitsShower(showerHandle, event, fShowerModuleLabel);

    
    //----------------------------------------------------------------------------------------------------------
    //                                     NUMUCC SELECTION
    //----------------------------------------------------------------------------------------------------------

    // Create list of primary muons and primary pfps which contain the most of them
    std::map<int, std::pair<size_t, int>> muPfpMap;
    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg = std::abs(pParticle->PdgCode());
      if(!(pdg == pandora::NU_E || pdg == pandora::NU_MU || pdg == pandora::NU_TAU)) continue;

      size_t pfp_key = pParticle.key();
      // Loop over daughters of pfparticle and do some truth matching
      // Assign labels based on the particle constributing the most hits
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        std::vector<art::Ptr<recob::Hit>> all_hits;
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));

        // Add up track and shower energy
        for(size_t i = 0; i < associatedTracks.size(); i++){
          std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(associatedTracks[i]->ID());
          all_hits.insert(all_hits.end(), hits.begin(), hits.end());
        }

        // Get the true particle associated with the daughter
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(all_hits, false);
        if(particles.find(trueId) == particles.end()) continue;
        // Only care about primary muons
        if(std::abs(particles[trueId].PdgCode()) != 13) continue;
        if(particles[trueId].Mother() != 0) continue;
        if(particles[trueId].StatusCode() != 1) continue;

        // Match the mother PFParticle to the primary muons if muon not matched already
        if(muPfpMap.find(trueId) == muPfpMap.end()) muPfpMap[trueId] = std::make_pair(pfp_key, all_hits.size());
        // If muon is matched, match this pfp if it shares more hits with the primary muon
        else if(muPfpMap[trueId].second < (int)all_hits.size()) std::make_pair(pfp_key, all_hits.size());
      }
    }

    for (size_t i = 0; i < mctruthList.size(); i++){

      ResetNuMuVars();

      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> truth = mctruthList.at(i);

      if(truth->Origin() != simb::kBeamNeutrino) continue;

      // Push back all unique neutrino energies
      nu_pdg = truth->GetNeutrino().Nu().PdgCode();
      if(!(truth->GetNeutrino().CCNC() == simb::kCC && std::abs(nu_pdg) == 14)) continue;
      is_cc = true;
      nu_int = truth->GetNeutrino().Mode();

      // Get truth info if numuCC in AV
      geo::Point_t vtx {truth->GetNeutrino().Nu().Vx(), 
                        truth->GetNeutrino().Nu().Vy(), 
                        truth->GetNeutrino().Nu().Vz()};
      if(!fTpcGeo.InFiducial(vtx, 0.)) continue;

      true_vtx_x = vtx.X();
      true_vtx_y = vtx.Y();
      true_vtx_z = vtx.Z();
      true_nu_energy = truth->GetNeutrino().Nu().E();

      // Get the primary muon
      std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(truth);
      int mu_id = -99999;
      for(auto const& part : parts){
        int part_pdg = std::abs(part->PdgCode());
        if((part_pdg == 211 || part_pdg == 321 || part_pdg == 2212) && part->E()>0.01) true_n_tracks++; 
        // Find the primary muon
        if(part_pdg != 13) continue;
        if(part->Mother() != 0) continue;
        if(part->StatusCode() != 1) continue;

        true_length = fTpcGeo.TpcLength(*part);
        true_mom = part->P();
        TVector3 start(part->Vx(), part->Vy(), part->Vz());
        TVector3 end(part->EndX(), part->EndY(), part->EndZ());
        true_theta = (end-start).Theta();
        true_phi = (end-start).Phi();
        true_cont = fTpcGeo.IsContained(*part);
        mu_id = part->TrackId();
      }

      if(muPfpMap.find(mu_id) == muPfpMap.end()){
        reconstructed = false;
        fNuMuTree->Fill();
        continue;
      }

      size_t pfp_key = muPfpMap[mu_id].first;

      for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

        const art::Ptr<recob::PFParticle> pParticle(it->second);
        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;
        // Check if this particle is identified as the neutrino
        const int pdg = std::abs(pParticle->PdgCode());
        if(!(pdg == pandora::NU_E || pdg == pandora::NU_MU || pdg == pandora::NU_TAU)) continue;

        if(pParticle.key() != pfp_key) continue;

        std::vector<recob::Track> nuTracks;
        std::map<int, double> track_energies;
        double shower_energy = 0;
       
        // Loop over daughters of pfparticle and do some truth matching
        // Assign labels based on the particle constributing the most hits
        for (const size_t daughterId : pParticle->Daughters()){
       
          // Get tracks associated with daughter
          art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
          const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
          const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pDaughter.key()));
       
          // Add up track and shower energy
          for(size_t i = 0; i < associatedTracks.size(); i++){
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
            track_energies[associatedTracks[i]->ID()] = calos[best_plane]->KineticEnergy()/1e3;
          }
          for(size_t i = 0; i < associatedShowers.size(); i++){
            std::vector<art::Ptr<recob::Hit>> hits = findManyHitsShower.at(associatedShowers[i]->ID());
            for(size_t j = 0; j < hits.size(); j++){
              // Assume collection plane is best for showers
              if(hits[j]->WireID().Plane == 2){ 
                shower_energy += HitEnergy(hits[j]);
              }
            }
          }
       
          if(associatedTracks.size() != 1) continue;
       
          // Get the first associated track
          recob::Track tpcTrack = *associatedTracks.front();
          nuTracks.push_back(tpcTrack);
       
        }
       
        reco_n_tracks = nuTracks.size() - 1;
       
        // Skip any PFPs without any tracks in them
        if(nuTracks.size() > 0) reconstructed = true;
        else{
          reconstructed = false;
          continue;
        }
       
        // Does pfp look like a cosmic?
        cosmic_id = cosIdAlg.CosmicId(*pParticle, pfParticleMap, event);
       
        // -------------------------------------- APPLY SELECTION ---------------------------------------
        std::pair<std::pair<bool, bool>, recob::Track> sel_track = Selection(nuTracks, findManyPid, findManyCalo);

        selected = sel_track.first.first;
        in_fv = sel_track.first.second;
       
        // Calculate kinematic variables for prop_sel_track prop_track
        reco_cont = fTpcGeo.InFiducial(sel_track.second.End(), 1.5);
        reco_length = sel_track.second.Length();
        reco_mom = 0.;
        if(reco_cont){
          reco_mom = fRangeFitter.GetTrackMomentum(reco_length, 13);
        }
        else{
          recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(sel_track.second);
          reco_mom = mcsResult.bestMomentum();
        }
        reco_theta = sel_track.second.Theta();
        reco_phi = sel_track.second.Phi();
        reco_vtx_x = sel_track.second.Start().X();
        reco_vtx_y = sel_track.second.Start().Y();
        reco_vtx_z = sel_track.second.Start().Z();
        reco_nu_energy = NeutrinoEnergy(track_energies, shower_energy, sel_track.second.ID());

      }

      fNuMuTree->Fill();
    }

  } // NumuccTree::analyze()


  void NumuccTree::endJob(){

  } // NumuccTree::endJob()

  void NumuccTree::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  void NumuccTree::ResetNuMuVars(){
    nu_pdg = -99999;
    is_cc = false;
    nu_int = -99999;
    true_vtx_x = -99999;
    true_vtx_y = -99999;
    true_vtx_z = -99999;
    true_nu_energy = -99999;
    true_length = -99999;
    true_mom = -99999;
    true_theta = -99999;
    true_phi = -99999;
    true_cont = false;
    true_n_tracks = 0;
    reconstructed = false;
    cosmic_id = false;
    in_fv = false;
    selected = false;
    reco_vtx_x = -99999;
    reco_vtx_y = -99999;
    reco_vtx_z = -99999;
    reco_nu_energy = -99999;
    reco_length = -99999;
    reco_mom = -99999;
    reco_theta = -99999;
    reco_phi = -99999;
    reco_cont = false;
    reco_n_tracks = 0;

  }

  // Apply my selection
  std::pair<std::pair<bool, bool>, recob::Track> NumuccTree::Selection(std::vector<recob::Track> tracks, art::FindMany<anab::ParticleID> fmpid, art::FindManyP<anab::Calorimetry> fmcalo){

    bool is_selected = false;
    bool has_candidate = false;
    recob::Track candidate = tracks[0];

    // Loop over tracks and count how many escape
    // For contained tracks apply PID cuts to only retain muon-like tracks
    int n_escape = 0;
    double longest_escape = 0;
    std::vector<recob::Track> long_tracks;
    for(size_t i = 0; i < tracks.size(); i++){
      // Find longest escaping track and don't apply track cuts if escaping
      if(!fTpcGeo.InFiducial(tracks[i].End(), 1.5)){
        n_escape++;
        double length = tracks[i].Length();
        if(length > longest_escape){ 
          longest_escape = length;
        }
        long_tracks.push_back(tracks[i]);
        continue;
      }

      // Select if longer than 150 cm
      if(tracks[i].Length() > 100.){
        long_tracks.push_back(tracks[i]);
        continue;
      }

      // Loop over planes (Y->V->U) and choose the next plane's calorimetry if there are 1.5x more points (collection plane more reliable)
      std::vector<art::Ptr<anab::Calorimetry>> calos = fmcalo.at(tracks[i].ID());
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

      // Get rid of any protons using chi2
      std::vector<const anab::ParticleID*> pids = fmpid.at(tracks[i].ID());
      bool is_proton = false;
      for(size_t i = 0; i < pids.size(); i++){
        // Only use the collection plane
        if(pids[i]->PlaneID().Plane != best_plane) continue;
        // If minimum chi2 is proton then ignore
        if(pids[i]->Chi2Proton() < pids[i]->Chi2Muon() && pids[i]->Chi2Proton() < pids[i]->Chi2Pion()){ 
          is_proton = true;
          continue;
        }
      }
      if(is_proton) continue;

      // Get rid of tracks which don't scatter like muons
      std::vector<float> angles = fMcsFitter.fitMcs(tracks[i], 13).scatterAngles();
      double ave_angle = std::accumulate(angles.begin(), angles.end(), 0)/angles.size();
      if(AverageDCA(tracks[i]) < 0.2 || ave_angle < 30) continue;

      // Get rid of any contained particles which don't stop (most muons do)
      double stop_chi2 = fStopTagger.StoppingChiSq(tracks[i].End(), calos);
      if(stop_chi2 < 1.2) continue;

      // Reject any tracks shorter than 25 cm
      if(tracks[i].Length() < 25) continue;

      // Check momentum reconstruction quality
      double range_mom = fRangeFitter.GetTrackMomentum(tracks[i].Length(), 13);
      double mcs_mom = fMcsFitter.fitMcs(tracks[i], 13).bestMomentum();
      double mom_diff = (mcs_mom - range_mom)/range_mom;
      double mom_diff_limit = 0.5 + std::exp(-(tracks[i].Length()-15.)/10.);
      if(mom_diff > mom_diff_limit) continue;

      long_tracks.push_back(tracks[i]);
    }

    std::sort(long_tracks.begin(), long_tracks.end(), [](auto& left, auto& right){
              return left.Length() > right.Length();});

    // Case 1: 1 escaping track
    if(n_escape == 1 && long_tracks.size() > 0){
      // If escaping track is longest and length > 50 cm then ID as muon
      if(longest_escape == long_tracks[0].Length() && longest_escape > 50){
        has_candidate = true;
        candidate = long_tracks[0];
      }
    }

    // Case 2: All tracks contained
    else if(n_escape == 0 && long_tracks.size() > 0){
      // Select the longest track
      has_candidate = true;
      candidate = long_tracks[0];
    }

    // Check vertex (start of muon candidate) in FV
    // Fiducial definition 50 cm from back, 10 cm from left, right and bottom, 15 cm from front, 20 cm from top 
    // 5 cm either side of CPA, 2.5 cm either side of APA gap
    bool in_fiducial = fTpcGeo.InFiducial(candidate.Start(), 10., 10., 15., 10., 20., 50., 5., 2.5);
    if(has_candidate) is_selected = true;

    return std::make_pair(std::make_pair(is_selected, in_fiducial), candidate);

  }

 double NumuccTree::AverageDCA(const recob::Track& track){

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

  double NumuccTree::HitEnergy(art::Ptr<recob::Hit> hit){

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

  double NumuccTree::NeutrinoEnergy(std::map<int, double> track_energies, double shower_energy, int id){
    double nu_e = 0;
    if(shower_energy>0 && shower_energy<5) nu_e = shower_energy;
    for(auto const& kv : track_energies){
      if(kv.first != id && kv.second>0 && kv.second<5) nu_e += kv.second;
    }
    return nu_e;
  }

  DEFINE_ART_MODULE(NumuccTree)
} // namespace sbnd

