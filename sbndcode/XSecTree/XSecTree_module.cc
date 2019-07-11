////////////////////////////////////////////////////////////////////////
// Class:       XSecTree
// Module Type: analyzer
// File:        XSecTree_module.cc
//
// Analysis module for selecting cross section distributions from truth
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom2.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace sbnd {

  class XSecTree : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> GenModuleLabel {
        Name("GenModuleLabel"),
        Comment("tag of generator level data product")
      };

      fhicl::Atom<art::InputTag> G4ModuleLabel {
        Name("G4ModuleLabel"),
        Comment("tag of geant4 level data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<double> WallCut {
        Name("WallCut"),
        Comment("Fiducial cut from all walls but back [cm]")
      };

      fhicl::Atom<double> BackCut {
        Name("BackCut"),
        Comment("Fiducial cut from all back wall [cm]")
      };

      fhicl::Atom<double> MinContainedLength {
        Name("MinContainedLength"),
        Comment("Minimum length of longest particle if contained [cm]")
      };

      fhicl::Atom<double> MinExitingLength {
        Name("MinExitingLength"),
        Comment("Minimum length of longest particle if exiting [cm]")
      };

      fhicl::Atom<double> MuonThreshold {
        Name("MuonThreshold"),
        Comment("Momentum threshold for reconstructing muons [GeV]")
      };

      fhicl::Atom<double> Pi0Threshold {
        Name("Pi0Threshold"),
        Comment("Momentum threshold for reconstructing pi0s [GeV]")
      };

      fhicl::Atom<double> PhotonThreshold {
        Name("PhotonThreshold"),
        Comment("Momentum threshold for reconstructing photons [GeV]")
      };

      fhicl::Atom<double> PionThreshold {
        Name("PionThreshold"),
        Comment("Momentum threshold for reconstructing charged pions [GeV]")
      };

      fhicl::Atom<double> ProtonThreshold {
        Name("ProtonThreshold"),
        Comment("Momentum threshold for reconstructing protons [GeV]")
      };

      fhicl::Atom<double> MuonEff {
        Name("MuonEff"),
        Comment("Efficiency for reconstructing muon [%]")
      };

      fhicl::Atom<double> Pi0Eff {
        Name("Pi0Eff"),
        Comment("Efficiency for reconstructing pi0s [%]")
      };

      fhicl::Atom<double> PhotonEff {
        Name("PhotonEff"),
        Comment("Efficiency for reconstructing photons [%]")
      };

      fhicl::Atom<double> PionEff {
        Name("PionEff"),
        Comment("Efficiency for reconstructing pions [%]")
      };

      fhicl::Atom<double> ProtonEff {
        Name("ProtonEff"),
        Comment("Efficiency for reconstructing protons [%]")
      };
      
      fhicl::Atom<double> PionPidEff {
        Name("PionPidEff"),
        Comment("PID efficiency for pions [%]")
      };

      fhicl::Atom<double> ProtonPidEff {
        Name("ProtonPidEff"),
        Comment("PID efficiency for protons [%]")
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit XSecTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called for every sub run
    virtual void beginSubRun(const art::SubRun& subrun) override;

    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset all counters/variables to their null values
    void ResetVars();

    // Give us a list of the stable, primary particles that we're interested in
    std::vector<const simb::MCParticle*> InterestingParticles(std::vector<const simb::MCParticle*> particles);

    // Apply reconstruction efficiencies to true particles
    std::vector<const simb::MCParticle*> RecoParticles(std::vector<const simb::MCParticle*> particles);

    // Smear momentum for exiting particles using MCS based method
    double SmearMcsMomentum(double momentum);

    // Smear momentum for contained particles using range based method
    double SmearRangeMomentum(double momentum);

    // Calculate the visible energy as neutrino energy estimator
    double VisibleEnergy(std::vector<const simb::MCParticle*> particles);

    // Calculate transverse variables (https://link.aps.org/accepted/10.1103/PhysRevC.94.015503)
    TVector3 DeltaPT(TVector3 mu_mom, TVector3 pr_mom);

    // Inclusive charged current selection
    bool IsCCInc(geo::Point_t vertex, std::vector<const simb::MCParticle*> reco_particles);

  private:

    // fcl file parameters
    art::InputTag fGenModuleLabel;      ///< name of gen producer
    art::InputTag fG4ModuleLabel;       ///< name of g4 producer
    bool          fVerbose;             ///< print information about what's going on

    double        fWallCut;
    double        fBackCut;

    double fMinContainedLength;
    double fMinExitingLength;

    double fMuonThreshold;
    double fPi0Threshold;
    double fPhotonThreshold;
    double fPionThreshold;
    double fProtonThreshold;

    double fMuonEff;
    double fPi0Eff;
    double fPhotonEff;
    double fPionEff;
    double fProtonEff;

    double fPionPidEff;
    double fProtonPidEff;

    TPCGeoAlg fTPCGeo;
    trkf::TrackMomentumCalculator fRangeFitter;

    TRandom2 *fRandom;

    // Global variables
    int lep_j;
    int longest_j;

    // Tree
    TTree *fXSecTree;
    TTree *fMetaDataTree;

    // XSec true tree parameters
    bool true_particles_contained;
    bool true_lep_contained;
    int true_cc;
    int true_nu_pdg;
    int true_int_type;
    unsigned int true_n_pipm;
    unsigned int true_n_pi0;
    unsigned int true_n_pr;
    double true_nu_energy;
    double true_lep_mom;
    double true_lep_theta;
    double true_pr1_mom;
    double true_pr1_theta;
    double true_lep_pr1_angle;
    double true_pipm1_mom;
    double true_pipm1_theta;
    double true_lep_pipm1_angle;
    double true_delta_pt;
    double true_delta_alphat;
    double true_delta_phit;

    // XSec reco tree parameters
    bool reco_particles_contained;
    bool reco_lep_contained;
    int reco_cc;
    int reco_nu_pdg;
    int reco_int_type;
    unsigned int reco_n_pipm;
    unsigned int reco_n_pi0;
    unsigned int reco_n_pr;
    double reco_nu_energy;
    double reco_lep_mom;
    double reco_lep_theta;
    double reco_pr1_mom;
    double reco_pr1_theta;
    double reco_lep_pr1_angle;
    double reco_pipm1_mom;
    double reco_pipm1_theta;
    double reco_lep_pipm1_angle;
    double reco_delta_pt;
    double reco_delta_alphat;
    double reco_delta_phit;

    // MetaData tree parameters
    double pot;

  }; // class XSecTree


  // Constructor
  XSecTree::XSecTree(Parameters const& config)
    : EDAnalyzer(config)
    , fGenModuleLabel       (config().GenModuleLabel())
    , fG4ModuleLabel        (config().G4ModuleLabel())
    , fVerbose              (config().Verbose())
    , fWallCut              (config().WallCut())
    , fBackCut              (config().BackCut())
    , fMinContainedLength   (config().MinContainedLength())
    , fMinExitingLength     (config().MinExitingLength())
    , fMuonThreshold        (config().MuonThreshold())
    , fPi0Threshold         (config().Pi0Threshold())
    , fPhotonThreshold      (config().PhotonThreshold())
    , fPionThreshold        (config().PionThreshold())
    , fProtonThreshold      (config().ProtonThreshold())
    , fMuonEff              (config().MuonEff())
    , fPi0Eff               (config().Pi0Eff())
    , fPhotonEff            (config().PhotonEff())
    , fPionEff              (config().PionEff())
    , fProtonEff            (config().ProtonEff())
    , fPionPidEff           (config().PionPidEff())
    , fProtonPidEff         (config().ProtonPidEff())
  {

  }


  void XSecTree::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    // Define histograms
    fXSecTree = tfs->make<TTree>("interaction", "xsec tree");

    fXSecTree->Branch("true_particles_contained", &true_particles_contained, "true_particles_contained/O");
    fXSecTree->Branch("true_lep_contained", &true_lep_contained, "true_lep_contained/O");
    fXSecTree->Branch("true_cc", &true_cc, "true_cc/I");
    fXSecTree->Branch("true_nu_pdg", &true_nu_pdg, "true_nu_pdg/I");
    fXSecTree->Branch("true_int_type", &true_int_type, "true_int_type/I");
    fXSecTree->Branch("true_n_pipm", &true_n_pipm, "true_n_pipm/i");
    fXSecTree->Branch("true_n_pi0", &true_n_pi0, "true_n_pi0/i");
    fXSecTree->Branch("true_n_pr", &true_n_pr, "true_n_pr/i");
    fXSecTree->Branch("true_nu_energy", &true_nu_energy, "true_nu_energy/D");
    fXSecTree->Branch("true_lep_mom", &true_lep_mom, "true_lep_mom/D");
    fXSecTree->Branch("true_lep_theta", &true_lep_theta, "true_lep_theta/D");
    fXSecTree->Branch("true_pr1_mom", &true_pr1_mom, "true_pr1_mom/D");
    fXSecTree->Branch("true_pr1_theta", &true_pr1_theta, "true_pr1_theta/D");
    fXSecTree->Branch("true_lep_pr1_angle", &true_lep_pr1_angle, "true_lep_pr1_angle/D");
    fXSecTree->Branch("true_pipm1_mom", &true_pipm1_mom, "true_pipm1_mom/D");
    fXSecTree->Branch("true_pipm1_theta", &true_pipm1_theta, "true_pipm1_theta/D");
    fXSecTree->Branch("true_lep_pipm1_angle", &true_lep_pipm1_angle, "true_lep_pipm1_angle/D");
    fXSecTree->Branch("true_delta_pt", &true_delta_pt, "true_delta_pt/D");
    fXSecTree->Branch("true_delta_alphat", &true_delta_alphat, "true_delta_alphat/D");
    fXSecTree->Branch("true_delta_phit", &true_delta_phit, "true_delta_phit/D");

    fXSecTree->Branch("reco_particles_contained", &reco_particles_contained, "reco_particles_contained/O");
    fXSecTree->Branch("reco_lep_contained", &reco_lep_contained, "reco_lep_contained/O");
    fXSecTree->Branch("reco_cc", &reco_cc, "reco_cc/I");
    fXSecTree->Branch("reco_nu_pdg", &reco_nu_pdg, "reco_nu_pdg/I");
    fXSecTree->Branch("reco_int_type", &reco_int_type, "reco_int_type/I");
    fXSecTree->Branch("reco_n_pipm", &reco_n_pipm, "reco_n_pipm/i");
    fXSecTree->Branch("reco_n_pi0", &reco_n_pi0, "reco_n_pi0/i");
    fXSecTree->Branch("reco_n_pr", &reco_n_pr, "reco_n_pr/i");
    fXSecTree->Branch("reco_nu_energy", &reco_nu_energy, "reco_nu_energy/D");
    fXSecTree->Branch("reco_lep_mom", &reco_lep_mom, "reco_lep_mom/D");
    fXSecTree->Branch("reco_lep_theta", &reco_lep_theta, "reco_lep_theta/D");
    fXSecTree->Branch("reco_pr1_mom", &reco_pr1_mom, "reco_pr1_mom/D");
    fXSecTree->Branch("reco_pr1_theta", &reco_pr1_theta, "reco_pr1_theta/D");
    fXSecTree->Branch("reco_lep_pr1_angle", &reco_lep_pr1_angle, "reco_lep_pr1_angle/D");
    fXSecTree->Branch("reco_pipm1_mom", &reco_pipm1_mom, "reco_pipm1_mom/D");
    fXSecTree->Branch("reco_pipm1_theta", &reco_pipm1_theta, "reco_pipm1_theta/D");
    fXSecTree->Branch("reco_lep_pipm1_angle", &reco_lep_pipm1_angle, "reco_lep_pipm1_angle/D");
    fXSecTree->Branch("reco_delta_pt", &reco_delta_pt, "reco_delta_pt/D");
    fXSecTree->Branch("reco_delta_alphat", &reco_delta_alphat, "reco_delta_alphat/D");
    fXSecTree->Branch("reco_delta_phit", &reco_delta_phit, "reco_delta_phit/D");

    fMetaDataTree = tfs->make<TTree>("metadata", "xsec tree");
    fMetaDataTree->Branch("pot", &pot, "pot/D");

    // Initial output
    std::cout<<"----------------- XSec Tree Module -------------------"<<std::endl;

    fRandom = new TRandom2();

  } // XSecTree::beginJob()


  // Called for every sub run
  void XSecTree::beginSubRun(const art::SubRun& subrun){

    art::Handle< sumdata::POTSummary > potHandle;
    subrun.getByLabel(fGenModuleLabel, potHandle);
    const sumdata::POTSummary& potSum = (*potHandle);
    pot = potSum.totpot;

    fMetaDataTree->Fill();
    return;
  } // XSecTree::beginSubRun()


  void XSecTree::analyze(const art::Event& event)
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
    // Retrieve all the truth info in the events
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::Handle<std::vector<simb::MCTruth>> genHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthList;
    if(event.getByLabel(fGenModuleLabel, genHandle)) art::fill_ptr_vector(mctruthList, genHandle);

    //----------------------------------------------------------------------------------------------------------
    //                                           FILLING THE TREE
    //----------------------------------------------------------------------------------------------------------
    // Loop over all the neutrino interactions
    for (size_t i = 0; i < mctruthList.size(); i++){
      if(fVerbose) std::cout<<"\n\nNeutrino: "<<i<<"\n";
      // Reset all the tree variables
      ResetVars();

      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> mctruth = mctruthList[i];

      // Check the interaction is within the TPC
      geo::Point_t vertex {mctruth->GetNeutrino().Nu().Vx(), 
                           mctruth->GetNeutrino().Nu().Vy(), 
                           mctruth->GetNeutrino().Nu().Vz()};
      if(fVerbose) std::cout<<"->Vertex: ("<<vertex.X()<<", "<<vertex.Y()<<", "<<vertex.Z()<<")\n";
      if(!fTPCGeo.InFiducial(vertex, 0.)) continue;

      //---------------------- FILLING ALL THE TRUE INTERACTION PARAMETERS ------------------------------------

      // Fill all the neutrino parameters
      if(mctruth->GetNeutrino().CCNC() == simb::kCC) true_cc = 1;
      else true_cc = 0;
      true_nu_pdg = mctruth->GetNeutrino().Nu().PdgCode();
      true_int_type = mctruth->GetNeutrino().Mode();
      true_nu_energy = mctruth->GetNeutrino().Nu().E();

      // Get all the particles propagated by geant 4
      std::vector<const simb::MCParticle*> parts = pi_serv->MCTruthToParticles_Ps(mctruth);
      // Get only the interesting ones
      std::vector<const simb::MCParticle*> particles = InterestingParticles(parts);

      // Lepton stuff
      for(size_t j = 0; j < particles.size(); j++){
        int pdg = std::abs(particles[j]->PdgCode());
        // Only consider the primary muon
        if(!(pdg == 13)) continue;
        lep_j = j;
        true_lep_mom = particles[j]->P();
        TVector3 start = particles[j]->Position().Vect();
        TVector3 end = particles[j]->EndPosition().Vect();
        true_lep_theta = (end - start).Theta();
        true_lep_contained = fTPCGeo.IsContained(*particles[j]);
      }

      // Secondary particle stuff
      for(size_t j = 0; j < particles.size(); j++){

        if((int)j == lep_j) continue;

        // Only consider pi0, charged pi and protons
        int pdg = std::abs(particles[j]->PdgCode());

        TVector3 start = particles[j]->Position().Vect();
        TVector3 end = particles[j]->EndPosition().Vect();

        // Count particle numbers
        if(pdg == 111) true_n_pi0++;
        if(pdg == 211){ 
          true_n_pipm++;
          if(particles[j]->P() < true_pipm1_mom) continue;
          true_pipm1_mom = particles[j]->P();
          true_pipm1_theta = (end - start).Theta();
          if(!fTPCGeo.IsContained(*particles[j])) true_particles_contained = false;

          // Calculate things relative to the lepton if there is one
          if(lep_j == -1) continue;
          TVector3 lep_end = particles[lep_j]->EndPosition().Vect();
          true_lep_pipm1_angle = lep_end.Angle(end);
        }
        if(pdg == 2212){ 
          true_n_pr++;

          if(particles[j]->P() < true_pr1_mom) continue;
          true_pr1_mom = particles[j]->P();
          true_pr1_theta = (end - start).Theta();
          if(!fTPCGeo.IsContained(*particles[j])) true_particles_contained = false;

          // Calculate things relative to the lepton if there is one
          if(lep_j == -1) continue;
          TVector3 lep_end = particles[lep_j]->EndPosition().Vect();
          true_lep_pr1_angle = lep_end.Angle(end);

          // Transverse variables
          TVector3 delta_pt = DeltaPT(particles[lep_j]->Momentum().Vect(), particles[j]->Momentum().Vect());
          true_delta_pt = delta_pt.Mag();
          true_delta_alphat = delta_pt.Theta()*TMath::RadToDeg();
          true_delta_phit = delta_pt.Phi()*TMath::RadToDeg();
        }

      }

      //---------------------------- CC INCLUSIVE SELECTION --------------------------------------

      // TODO Make a better neutrino energy calculator
      reco_nu_energy = VisibleEnergy(particles);

      // Apply reconstruction efficiencies to the true particles
      std::vector<const simb::MCParticle*> reco_particles = RecoParticles(parts);
      if(fVerbose) std::cout<<"Number of reconstructed particles = "<<reco_particles.size()<<"\n";

      bool cc_selected = IsCCInc(vertex, reco_particles); 

      if(cc_selected){ 
        reco_cc = 1;
        reco_nu_pdg = 14;
      }
      else{ 
        reco_cc = -1;
      }

      //------------------------------------ RECO FSI ------------------------------------------

      for(size_t j = 0; j < reco_particles.size(); j++){
        // Don't look at the particle selected as the muon
        if(longest_j == (int)j) continue;

        int pdg = std::abs(reco_particles[j]->PdgCode());

        double contained_length = fTPCGeo.TpcLength(*reco_particles[j]);
        TVector3 start = reco_particles[j]->Position().Vect();
        TVector3 end = reco_particles[j]->EndPosition().Vect();

        // Apply PID estimation and count reco particles
        if(pdg == 111) reco_n_pi0++;

        // Treat muons and pions as the same
        if(pdg == 13 || pdg == 211){
          // Check that particle is contained
          if(!fTPCGeo.IsContained(*reco_particles[j])) reco_particles_contained = false;

          double rand_pid = fRandom->Rndm();
          // ID as pion, calculate leading pion variables
          if(rand_pid < fPionPidEff){ 
            reco_n_pipm++;

            // Leading pion momentum and angle
            double pipm1_mom = fRangeFitter.GetTrackMomentum(contained_length, 13);
            if(pipm1_mom < reco_pipm1_mom) continue;
            // FIXME no idea how to estimate pion momentum
            reco_pipm1_mom = pipm1_mom;
            reco_pipm1_theta = (end - start).Theta();

            // Angle between leading pion and lepton candidate
            if(longest_j == -1) continue;
            TVector3 lep_end = reco_particles[longest_j]->EndPosition().Vect();
            reco_lep_pipm1_angle = lep_end.Angle(end);
          }

          // ID as proton, calculate leading proton variables
          else{
            reco_n_pr++;

            // Leading proton momentum and angle
            double pr1_mom = fRangeFitter.GetTrackMomentum(contained_length, 2212);
            if(pr1_mom < reco_pr1_mom) continue;
            reco_pr1_mom = pr1_mom;
            reco_pr1_theta = (end - start).Theta();

            // Angle between leading proton and lepton candidate
            if(longest_j == -1) continue;
            TVector3 lep_end = reco_particles[longest_j]->EndPosition().Vect();
            reco_lep_pr1_angle = lep_end.Angle(end);

            // Transverse variables
            TVector3 delta_pt = DeltaPT(reco_particles[longest_j]->Momentum().Vect(), reco_particles[j]->Momentum().Vect());
            reco_delta_pt = delta_pt.Mag();
            reco_delta_alphat = delta_pt.Theta()*TMath::RadToDeg();
            reco_delta_phit = delta_pt.Phi()*TMath::RadToDeg();
          }
        }

        if(pdg == 2212){
          // Check that particle is contained
          if(!fTPCGeo.IsContained(*reco_particles[j])) reco_particles_contained = false;

          double rand_pid = fRandom->Rndm();
          // ID as proton, calculate leading proton variables
          if(rand_pid < fProtonPidEff){ 
            reco_n_pr++;

            // Leading proton momentum and angle
            double pr1_mom = fRangeFitter.GetTrackMomentum(contained_length, 2212);
            if(pr1_mom < reco_pr1_mom) continue;
            reco_pr1_mom = pr1_mom;
            reco_pr1_theta = (end - start).Theta();

            // Angle between leading proton and lepton candidate
            if(longest_j == -1) continue;
            TVector3 lep_end = reco_particles[longest_j]->EndPosition().Vect();
            reco_lep_pr1_angle = lep_end.Angle(end);

            // Transverse variables
            TVector3 delta_pt = DeltaPT(reco_particles[longest_j]->Momentum().Vect(), reco_particles[j]->Momentum().Vect());
            reco_delta_pt = delta_pt.Mag();
            reco_delta_alphat = delta_pt.Theta()*TMath::RadToDeg();
            reco_delta_phit = delta_pt.Phi()*TMath::RadToDeg();
          }

          // ID as pion, calculate leading pion variables
          else{ 
            reco_n_pipm++;

            // Leading pion momentum and angle
            double pipm1_mom = fRangeFitter.GetTrackMomentum(contained_length, 13);
            if(pipm1_mom < reco_pipm1_mom) continue;
            reco_pipm1_mom = pipm1_mom;
            reco_pipm1_theta = (end - start).Theta();

            // Angle between leading pion and lepton candidate
            if(longest_j == -1) continue;
            TVector3 lep_end = reco_particles[longest_j]->EndPosition().Vect();
            reco_lep_pipm1_angle = lep_end.Angle(end);
          }
        }
          
      }

      
      if(fVerbose) std::cout<<"->true_cc:                  "<<true_cc<<"\n"
                            <<"->true_nu_pdg:              "<<true_nu_pdg<<"\n"
                            <<"->true_nu_energy:           "<<true_nu_energy<<"\n"
                            <<"->true_int_type:            "<<true_int_type<<"\n"
                            <<"->true_lep_contained:       "<<true_lep_contained<<"\n"
                            <<"-->true_lep_mom:            "<<true_lep_mom<<"\n"
                            <<"-->true_lep_theta:          "<<true_lep_theta<<"\n"
                            <<"->true_particles_contained: "<<true_particles_contained<<"\n"
                            <<"->true_n_pi0:               "<<true_n_pi0<<"\n"
                            <<"->true_n_pipm:              "<<true_n_pipm<<"\n"
                            <<"-->true_pipm1_mom:          "<<true_pipm1_mom<<"\n"
                            <<"-->true_pipm1_theta:        "<<true_pipm1_theta<<"\n"
                            <<"-->true_lep_pipm1_angle:    "<<true_lep_pipm1_angle<<"\n"
                            <<"->true_n_pr:                "<<true_n_pr<<"\n"
                            <<"-->true_pr1_mom:            "<<true_pr1_mom<<"\n"
                            <<"-->true_pr1_theta:          "<<true_pr1_theta<<"\n"
                            <<"-->true_lep_pr1_angle:      "<<true_lep_pr1_angle<<"\n"
                            <<"-->true_delta_pt:           "<<true_delta_pt<<"\n"
                            <<"-->true_delta_alphat:       "<<true_delta_alphat<<"\n"
                            <<"-->true_delta_phit:         "<<true_delta_phit<<"\n"
                            <<"\n"
                            <<"->reco_cc:                  "<<reco_cc<<"\n"
                            <<"->reco_nu_pdg:              "<<reco_nu_pdg<<"\n"
                            <<"->reco_nu_energy:           "<<reco_nu_energy<<"\n"
                            <<"->reco_int_type:            "<<reco_int_type<<"\n"
                            <<"->reco_lep_contained:       "<<reco_lep_contained<<"\n"
                            <<"-->reco_lep_mom:            "<<reco_lep_mom<<"\n"
                            <<"-->reco_lep_theta:          "<<reco_lep_theta<<"\n"
                            <<"->reco_particles_contained: "<<reco_particles_contained<<"\n"
                            <<"->reco_n_pi0:               "<<reco_n_pi0<<"\n"
                            <<"->reco_n_pipm:              "<<reco_n_pipm<<"\n"
                            <<"-->reco_pipm1_mom:          "<<reco_pipm1_mom<<"\n"
                            <<"-->reco_pipm1_theta:        "<<reco_pipm1_theta<<"\n"
                            <<"-->reco_lep_pipm1_angle:    "<<reco_lep_pipm1_angle<<"\n"
                            <<"->reco_n_pr:                "<<reco_n_pr<<"\n"
                            <<"-->reco_pr1_mom:            "<<reco_pr1_mom<<"\n"
                            <<"-->reco_pr1_theta:          "<<reco_pr1_theta<<"\n"
                            <<"-->reco_lep_pr1_angle:      "<<reco_lep_pr1_angle<<"\n"
                            <<"-->reco_delta_pt:           "<<reco_delta_pt<<"\n"
                            <<"-->reco_delta_alphat:       "<<reco_delta_alphat<<"\n"
                            <<"-->reco_delta_phit:         "<<reco_delta_phit<<"\n";


      fXSecTree->Fill();
    }

  } // XSecTree::analyze()


  void XSecTree::endJob(){

  } // XSecTree::endJob()


  // Reset variables and counters
  void XSecTree::ResetVars(){

    lep_j = -1;
    longest_j = -1;
    
    true_particles_contained = true;
    true_lep_contained       = false;
    true_cc                  = -1;
    true_nu_pdg              = -99999;
    true_int_type            = -99999;
    true_n_pipm              = 0;
    true_n_pi0               = 0;
    true_n_pr                = 0;
    true_nu_energy           = -99999;
    true_lep_mom             = -99999;
    true_lep_theta           = -99999;
    true_pr1_mom             = -99999;
    true_pr1_theta           = -99999;
    true_lep_pr1_angle       = -99999;
    true_pipm1_mom           = -99999;
    true_pipm1_theta         = -99999;
    true_lep_pipm1_angle     = -99999;
    true_delta_pt            = -99999;
    true_delta_alphat        = -99999;
    true_delta_phit          = -99999;

    reco_particles_contained = true;
    reco_lep_contained       = false;
    reco_cc                  = -1;
    reco_nu_pdg              = -99999;
    reco_int_type            = -99999;
    reco_n_pipm              = 0;
    reco_n_pi0               = 0;
    reco_n_pr                = 0;
    reco_nu_energy           = -99999;
    reco_lep_mom             = -99999;
    reco_lep_theta           = -99999;
    reco_pr1_mom             = -99999;
    reco_pr1_theta           = -99999;
    reco_lep_pr1_angle       = -99999;
    reco_pipm1_mom           = -99999;
    reco_pipm1_theta         = -99999;
    reco_lep_pipm1_angle     = -99999;
    reco_delta_pt            = -99999;
    reco_delta_alphat        = -99999;
    reco_delta_phit          = -99999;

  } // XSecTree::ResetVars

  
  // Give us a list of the stable, primary particles that we're interested in
  std::vector<const simb::MCParticle*> XSecTree::InterestingParticles(std::vector<const simb::MCParticle*> particles){

    std::vector<const simb::MCParticle*> interesting;

    // Loop over all of the particles
    for(size_t j = 0; j < particles.size(); j++){
      // Only consider stable final states particles
      if(particles[j]->StatusCode() != 1) continue;
      // Only want primary particles
      if(particles[j]->Mother() != 0) continue;

      // Only consider muons, pi0, charged pi and protons TODO for now...
      int pdg = std::abs(particles[j]->PdgCode());
      if(!(pdg == 13 || pdg == 111 || pdg == 211 || pdg == 2212)) continue;
      interesting.push_back(particles[j]);
    }

    return interesting;

  } // XSecTree::InterestingParticles()


  // Apply reconstruction efficiencies to true particles
  std::vector<const simb::MCParticle*> XSecTree::RecoParticles(std::vector<const simb::MCParticle*> particles){

    std::vector<const simb::MCParticle*> reco_particles;

    for(size_t j = 0; j < particles.size(); j++){
      // Only consider stable final states particles
      if(particles[j]->StatusCode() != 1) continue;
      // Only want primary particles
      if(particles[j]->Mother() != 0) continue;

      int pdg = std::abs(particles[j]->PdgCode());
      if(!(pdg == 13 || pdg == 111 || pdg == 211 || pdg == 2212)) continue;

      // Apply efficiency cut (energy thresholds + flat efficiency)
      double rand_eff = fRandom->Rndm();
      if(pdg == 13 && (particles[j]->P() < fMuonThreshold || rand_eff > fMuonEff)) continue;
      // Consider the secondary photons when reconstructing pi0
      if(pdg == 111){// && (particles[j]->P() < fPi0Threshold || rand_eff > fPi0Eff)) continue;
        int n_photons = 0;
        // Get the pi0 daughters
        std::vector<int> d_trk_ids;
        for(int k = 0; k < particles[j]->NumberDaughters(); k++){
          d_trk_ids.push_back(particles[j]->Daughter(k));
        }
        if(d_trk_ids.size() < 2) continue;
        // Match by track ID to g4 particles and apply threshold and effciencies
        for(size_t i = 0; i < particles.size(); i++){
          for(size_t k = 0; k < d_trk_ids.size(); k++){
            if(particles[i]->TrackId() != d_trk_ids[k]) continue;
            double rand_ph_eff = fRandom->Rndm();
            if(particles[i]->PdgCode() == 22 && 
             particles[i]->P() > fPhotonThreshold && 
             rand_ph_eff < fPhotonEff &&
             fTPCGeo.IsContained(*particles[i])) n_photons++;
          }
        }
        // Need to reconstruct 2 photons to reconstruct pi0
        if(n_photons < 2) continue;
      }
      if(pdg == 211 && (particles[j]->P() < fPionThreshold || rand_eff > fPionEff)) continue;
      if(pdg == 2212 && (particles[j]->P() < fProtonThreshold || rand_eff > fProtonEff)) continue;

      reco_particles.push_back(particles[j]);
    }

    return reco_particles;

  } // XSecTree::RecoParticles()


  // Smear momentum for exiting particles using MCS based method
  double XSecTree::SmearMcsMomentum(double momentum){
    //For exiting muons use multiple coulomb scattering bias and resolution
    //Values from Fig 12 of https://arxiv.org/pdf/1703.06187.pdf
    double bias[] = {0.0273,0.0409,0.0352,0.0250,0.0227,0.0068,0.0364,0.0273,0.0227};
    double resolution[] = {0.127,0.145,0.143,0.141,0.164,0.177,0.250,0.266,0.341};
    int pos = 8; 
    for(int i=0; i<9; i++){
      if(momentum<(0.34+0.41*(i+1))) {pos = i; break;}
    }
    double momentum_smear = fRandom->Gaus(momentum, resolution[pos]*momentum) + bias[pos]*momentum;
    if(momentum_smear<0) {momentum_smear = 0;}

    return momentum_smear;

  } // XSecTree::SmearMcsMomentum()


  // Smear momentum for contained particles using range based method
  double XSecTree::SmearRangeMomentum(double momentum){
    //For contained muons use range based bias and resolution
    //Values from Fig 5 of https://arxiv.org/pdf/1703.06187.pdf
    double bias[] = {-0.0035,-0.0059,-0.0047,-0.0059,-0.0035,-0.0029,-0.0076,-0.0059,0.0006};
    double resolution[] = {0.017,0.021,0.023,0.026,0.025,0.030,0.030,0.040,0.032};
    int pos = 8; 
    for(int i=0; i<9; i++){
      if(momentum<(0.33+0.186*(i+1))) {pos = i; break;}
    }
    double momentum_smear = fRandom->Gaus(momentum, resolution[pos]*momentum)+bias[pos]*momentum;
    if(momentum_smear<0){momentum_smear = 0;}

    return momentum_smear;

  } // XSecTree::SmearRangeMomentum()


  // Calculate the visible energy as neutrino energy estimator
  // TODO poorly copied from Gray's stuff, will likely change anyway
  double XSecTree::VisibleEnergy(std::vector<const simb::MCParticle*> particles){
    //
    double visible_E = 0;

    for(size_t j = 0; j < particles.size(); j++){
      // Only consider stable final states particles
      if(particles[j]->StatusCode() != 1) continue;
      // Only want primary particles
      if(particles[j]->Mother() != 0) continue;

      int pdg = std::abs(particles[j]->PdgCode());
      if(pdg == 111) continue;

      if(pdg == 13){
        double smearing_percentage;
        if(fTPCGeo.IsContained(*particles[j])) smearing_percentage = 0.02;
        else smearing_percentage = -0.102 * TMath::Log(0.000612*fTPCGeo.TpcLength(*particles[j]));
        double lep_E = particles[j]->E();
        visible_E += std::max(fRandom->Gaus(lep_E, smearing_percentage*lep_E), 0.);
        continue;
      }

      double mass = particles[j]->Mass();
      double this_visible_energy = (particles[j]->E() - mass); // GeV
      this_visible_energy = fRandom->Gaus(this_visible_energy, 0.05*this_visible_energy);
      this_visible_energy = std::max(this_visible_energy, 0.);
      if(this_visible_energy > 0.021) visible_E += this_visible_energy;
    }

    return visible_E;

  } // XSecTree::VisibleEnergy()


  // Calculate transverse variables (https://link.aps.org/accepted/10.1103/PhysRevC.94.015503)
  // TODO should this be just proton momentum or all other tracks?
  TVector3 XSecTree::DeltaPT(TVector3 mu_mom, TVector3 pr_mom){
    // Assume neutrino is directly along Z, magnitude unimportant
    TVector3 nu_mom(0, 0, 1);

    // Get the transverse momentum for the muon and proton
    TVector3 mu_rot(mu_mom);
    mu_rot.Rotate(TMath::Pi(), nu_mom);
    TVector3 mu_pt = (mu_mom - mu_rot)*0.5;

    TVector3 pr_rot(mu_mom);
    pr_rot.Rotate(TMath::Pi(), nu_mom);
    TVector3 pr_pt = (pr_mom - pr_rot)*0.5;

    // Calculate the variables
    TVector3 tot_pt = mu_pt + pr_pt;
    double phi = TMath::ACos(mu_pt.Dot(pr_pt)*(-1)/mu_pt.Mag()/pr_pt.Mag());
    double theta = TMath::ACos(tot_pt.Dot(mu_pt)*(-1)/tot_pt.Mag()/mu_pt.Mag());

    TVector3 delta_pt;
    delta_pt.SetMagThetaPhi(tot_pt.Mag(), theta, phi);

    return delta_pt;

  } // XSecTree::DeltaPT()


  // Inclusive charged current selection
  bool XSecTree::IsCCInc(geo::Point_t vertex, std::vector<const simb::MCParticle*> reco_particles){
    bool cc_selected = true;
    // Check vertex is inside the fiducial volume
    if(!fTPCGeo.InFiducial(vertex, fWallCut, fWallCut, fWallCut, fWallCut, fWallCut, fBackCut)){ 
      if(fVerbose) std::cout<<"Not in fiducial\n";
      cc_selected = false; 
    }

    // Loop over the mu/pi/pr secondary particles and find the longest
    double max_contained_length = 0;
    int max_pdg = -99999;
    for(size_t j = 0; j < reco_particles.size(); j++){
      // Only consider track-like particles
      int pdg = std::abs(reco_particles[j]->PdgCode());
      if(!(pdg == 13 || pdg == 211 || pdg == 2212)) continue;

      double contained_length = fTPCGeo.TpcLength(*reco_particles[j]);
      if(contained_length < max_contained_length) continue;
      max_contained_length = contained_length;
      max_pdg = pdg;
      longest_j = j;

      // Check if particle is contained
      reco_lep_contained = fTPCGeo.IsContained(*reco_particles[j]);

      TVector3 start = reco_particles[j]->Position().Vect();
      TVector3 end = reco_particles[j]->EndPosition().Vect();
      reco_lep_theta = (end - start).Theta();

      // Smear momentum based on whether particle is contained or not
      if(reco_lep_contained) reco_lep_mom = fRangeFitter.GetTrackMomentum(contained_length, 13);
      else reco_lep_mom = SmearMcsMomentum(reco_particles[j]->P());
    }

    if(fVerbose) std::cout<<"max contained length = "<<max_contained_length<<" pdg = "<<max_pdg<<"\n";

    // Check length of particle
    if(reco_lep_contained && max_contained_length < fMinContainedLength){ 
      cc_selected = false;
      if(fVerbose) std::cout<<"contained length too short\n";
    }
    if(!reco_lep_contained && max_contained_length < fMinExitingLength){ 
      cc_selected = false;
      if(fVerbose) std::cout<<"exiting length too short\n";
    }

    return cc_selected;

  }

  DEFINE_ART_MODULE(XSecTree)
} // namespace sbnd
