////////////////////////////////////////////////////////////////////////
// Class:       CosmicIdAna
// Module Type: analyzer
// File:        CosmicIdAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CosmicIdAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

// C++ includes
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace {
  // Local namespace for local functions
  // Declare here, define later

}

namespace sbnd {

  class CosmicIdAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
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

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CosmicIdAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    CosmicIdAlg cosIdAlg;
    CRTTruthRecoAlg truthAlg;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;
    
    // histograms
    std::vector<std::string> trueCategories{"NuMu","Dirt","Cr","Other"};
    size_t nTC = trueCategories.size();
    std::vector<std::string> recoCategories{"NuMu","Dirt","Cr","Other","OtherNu"};
    size_t nRC = recoCategories.size();
    std::vector<std::string> cuts{"None","FV","SP","Geo","CC","AC","CT","CH","PT","Tot","Remain"};
    size_t nCuts = cuts.size();

    TH1D *hTrueMom[4][11];
    TH1D *hTrueLength[4][11];
    TH1D *hTrueTheta[4][11];
    TH1D *hTruePhi[4][11];

    TH1D *hRecoMom[5][11];
    TH1D *hRecoLength[5][11];
    TH1D *hRecoTheta[5][11];
    TH1D *hRecoPhi[5][11];

    // performance counters

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

  }; // class CosmicIdAna


  // Constructor
  CosmicIdAna::CosmicIdAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fPandoraLabel         (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , cosIdAlg              (config().CosIdAlg())
    , fMcsFitter            (config().fitter)
  {

  } // CosmicIdAna()


  void CosmicIdAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    for(size_t i = 0; i < nTC; i++){
      for(size_t j = 0; j < nCuts; j++){
        TString hMom_name     = Form("hTrueMom%s_%s", trueCategories[i].c_str(), cuts[j].c_str());
        hTrueMom[i][j]        = tfs->make<TH1D>(hMom_name,    "", 20, 0,    2);
        TString hLength_name  = Form("hTrueLength%s_%s", trueCategories[i].c_str(), cuts[j].c_str());
        hTrueLength[i][j]     = tfs->make<TH1D>(hLength_name, "", 20, 0,    500);
        TString hTheta_name   = Form("hTrueTheta%s_%s", trueCategories[i].c_str(), cuts[j].c_str());
        hTrueTheta[i][j]      = tfs->make<TH1D>(hTheta_name,  "", 20, 0,    3.2);
        TString hPhi_name     = Form("hTruePhi%s_%s", trueCategories[i].c_str(), cuts[j].c_str());
        hTruePhi[i][j]        = tfs->make<TH1D>(hPhi_name,    "", 20, -3.2, 3.2);
      }
    }

    for(size_t i = 0; i < nRC; i++){
      for(size_t j = 0; j < nCuts; j++){
        TString hMom_name     = Form("hRecoMom%s_%s", recoCategories[i].c_str(), cuts[j].c_str());
        hRecoMom[i][j]        = tfs->make<TH1D>(hMom_name,    "", 20, 0,    2);
        TString hLength_name  = Form("hRecoLength%s_%s", recoCategories[i].c_str(), cuts[j].c_str());
        hRecoLength[i][j]     = tfs->make<TH1D>(hLength_name, "", 20, 0,    500);
        TString hTheta_name   = Form("hRecoTheta%s_%s", recoCategories[i].c_str(), cuts[j].c_str());
        hRecoTheta[i][j]      = tfs->make<TH1D>(hTheta_name,  "", 20, 0,    3.2);
        TString hPhi_name     = Form("hRecoPhi%s_%s", recoCategories[i].c_str(), cuts[j].c_str());
        hRecoPhi[i][j]        = tfs->make<TH1D>(hPhi_name,    "", 20, -3.2, 3.2);
      }
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- Cosmic ID Ana Module -------------------"<<std::endl;

  }// CosmicIdAna::beginJob()


  void CosmicIdAna::analyze(const art::Event& event)
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
    // Get PFParticle to track associations
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);
    
    // Get track to hit and colorimetry associations
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------

    // Record all true particles and sort by type
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;
    std::vector<int> dirtParticleIds;
    std::vector<int> crParticleIds;

    // Loop over all true particles
    for (auto const& particle: (*particleHandle)){
      // Store particle
      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
      // Get MCTruth
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      int pdg = std::abs(particle.PdgCode());

      // If origin is a neutrino
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        // If neutrino vertex is not inside the TPC then call it a dirt particle
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)){ 
          dirtParticleIds.push_back(partId);
        }
        // If it's a primary muon
        else if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
        }
        // Other nu particles
        else{
          nuParticleIds.push_back(partId);
        }
      }

      // If origin is a cosmic ray
      else if(truth->Origin() == simb::kCosmicRay){
        crParticleIds.push_back(partId);
      }
    }

    //----------------------------------------------------------------------------------------------------------
    //                                      FAKE PDS RECONSTRUCTION
    //----------------------------------------------------------------------------------------------------------

    // Create fake flashes in each tpc
    std::pair<std::vector<double>, std::vector<double>> fakeFlashes = CosmicRemovalUtils::FakeTpcFlashes(parts);
    std::vector<double> fakeTpc0Flashes = fakeFlashes.first;
    std::vector<double> fakeTpc1Flashes = fakeFlashes.second;
    bool tpc0BeamFlash = CosmicRemovalUtils::BeamFlash(fakeTpc0Flashes, 4.);
    bool tpc1BeamFlash = CosmicRemovalUtils::BeamFlash(fakeTpc1Flashes, 4.);

    // If there are no flashes in time with the beam then ignore the event
    if(!tpc0BeamFlash && !tpc1BeamFlash) return;

    //----------------------------------------------------------------------------------------------------------
    //                                          COSMIC ID - CALCULATING CUTS
    //----------------------------------------------------------------------------------------------------------

    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;

      int pfpType = 3;
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

        // Truth match muon tracks and pfps
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        int trackType = 3;
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ 
          trackType = 0;
          pfpType = 0;
        }
        else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){ 
          if(pfpType != 0) pfpType = 4;
        }
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()){ 
          trackType = 1;
          if(pfpType != 0 && pfpType != 4) pfpType = 1;
        }
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()){
          trackType = 2;
          if(pfpType != 0 && pfpType != 4 && pfpType != 1) pfpType = 2;
        }
        
        // Fill some histograms
        if(particles.find(trueId) != particles.end()){
          if(std::abs(particles[trueId].PdgCode()) == 13){
            std::pair<TVector3, TVector3> se = truthAlg.TpcCrossPoints(particles[trueId]);
            double momentum = particles[trueId].P();
            double length = truthAlg.TpcLength(particles[trueId]);
            double theta = (se.second-se.first).Theta();
            double phi = (se.second-se.first).Phi();
            for(size_t j = 0; j < nCuts; j++){
              bool plot = false;
              if(j == 0) plot = true;
              if(j == 1){
                cosIdAlg.SetCuts(true, false, false, false, false, false, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 2){
                cosIdAlg.SetCuts(false, true, false, false, false, false, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 3){
                cosIdAlg.SetCuts(false, false, true, false, false, false, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 4){
                cosIdAlg.SetCuts(false, false, false, true, false, false, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 5){
                cosIdAlg.SetCuts(false, false, false, false, true, false, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 6){
                cosIdAlg.SetCuts(false, false, false, false, false, true, false, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 7){
                cosIdAlg.SetCuts(false, false, false, false, false, false, true, false);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 8){
                cosIdAlg.SetCuts(false, false, false, false, false, false, false, true);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 9){
                cosIdAlg.SetCuts(true, true, true, true, true, true, true, true);
                if(cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              }
              if(j == 10 && !cosIdAlg.CosmicId(tpcTrack, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
              if(!plot) continue;
              hTrueMom[trackType][j]->Fill(momentum);
              hTrueLength[trackType][j]->Fill(length);
              hTrueTheta[trackType][j]->Fill(theta);
              hTruePhi[trackType][j]->Fill(phi);
            }
          }
        }

      }

      if(nuTracks.size() == 0) continue;

      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
                return left.Length() > right.Length();});

      recob::Track nuTrack = nuTracks[0];
      double recoMuMomentum = 0.;
      bool exits = CosmicRemovalUtils::InFiducial(nuTrack.End(), 5., 5.);
      double length = nuTrack.Length();
      double theta = nuTrack.Theta();
      double phi = nuTrack.Phi();
      if(exits){
        recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(nuTrack);
        recoMuMomentum = mcsResult.bestMomentum();
      }
      else{
        recoMuMomentum = fRangeFitter.GetTrackMomentum(length, 13);
      }
      for(size_t j = 0; j < nCuts; j++){
        bool plot = false;
        if(j == 0) plot = true;
        if(j == 1){
          cosIdAlg.SetCuts(true, false, false, false, false, false, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){
            plot = true;
          }
        }
        if(j == 2){
          cosIdAlg.SetCuts(false, true, false, false, false, false, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 3){
          cosIdAlg.SetCuts(false, false, true, false, false, false, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 4){
          cosIdAlg.SetCuts(false, false, false, true, false, false, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 5){
          cosIdAlg.SetCuts(false, false, false, false, true, false, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 6){
          cosIdAlg.SetCuts(false, false, false, false, false, true, false, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 7){
          cosIdAlg.SetCuts(false, false, false, false, false, false, true, false);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 8){
          cosIdAlg.SetCuts(false, false, false, false, false, false, false, true);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
            plot = true;
          }
        }
        if(j == 9){
          cosIdAlg.SetCuts(true, true, true, true, true, true, true, true);
          if(cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)) plot = true;
        }
        if(j == 10 && !cosIdAlg.CosmicId(*pParticle, pfParticleMap, event, fakeTpc0Flashes, fakeTpc1Flashes)){ 
          plot = true;
        }
        if(!plot) continue;
        hRecoMom[pfpType][j]->Fill(recoMuMomentum);
        hRecoLength[pfpType][j]->Fill(length);
        hRecoTheta[pfpType][j]->Fill(theta);
        hRecoPhi[pfpType][j]->Fill(phi);
      }

    }  
    
  } // CosmicIdAna::analyze()


  void CosmicIdAna::endJob(){

  } // CosmicIdAna::endJob()

  void CosmicIdAna::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  DEFINE_ART_MODULE(CosmicIdAna)
} // namespace sbnd

