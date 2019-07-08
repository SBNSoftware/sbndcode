////////////////////////////////////////////////////////////////////////
// Class:       PandoraNuMuRecoAna
// Module Type: analyzer
// File:        PandoraNuMuRecoAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
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

  class PandoraNuMuRecoAna : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Inputs {
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

      fhicl::Atom<art::InputTag> HitModuleLabel {
        Name("HitModuleLabel"),
        Comment("tag of hit producer data product")
      };

      fhicl::Atom<art::InputTag> ShowerModuleLabel {
        Name("ShowerModuleLabel"),
        Comment("tag of shower producer data product")
      };
      
      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Inputs

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Table<PandoraNuMuRecoAna::Inputs> inputs {
        Name("inputs"),
      };
      
      fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
        Name("fitter"),
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit PandoraNuMuRecoAna(Parameters const& config);
 
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
    art::InputTag fHitModuleLabel; ///< name of calorimetry producer
    art::InputTag fShowerModuleLabel;
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    // histograms
    std::vector<std::string> trueCategories{"","NuTrack","NuShower","CrTrack","Mix","NoReco"};
    size_t nTC = trueCategories.size();

    TH1D *hTrueMuLength[6];
    TH1D *hTrueMuMom[6];
    TH1D *hTrueMuTheta[6];
    TH1D *hTrueMuPhi[6];

    TH1D *hMuCompletenessEdep;
    TH1D *hMuCompletenessHitN;
    TH1D *hMuEdep;
    TH1D *hMuHitN;

    TH1D *hMuLengthDiff;
    TH1D *hMuThetaDiff;
    TH1D *hMuPhiDiff;
    TH1D *hMuAngle;

    TH1D *hMuVertexDiff;
    TH1D *hMuVertexXDiff;
    TH1D *hMuVertexYDiff;
    TH1D *hMuVertexZDiff;

    TH1D *hMcsDiff;
    TH1D *hMcsDiffCut;
    TH1D *hMcsDiffAll;
    TH1D *hMcsDiffAllCut;
    TH1D *hRangeDiff;
    TH1D *hRangeDiffCut;
    TH1D *hBothDiff;
    TH1D *hBothDiffCut;

    TH2D *hLengthTrue;
    TH2D *hLengthTrueDiff;
    TH2D *hThetaTrue;
    TH2D *hThetaTrueDiff;
    TH2D *hPhiTrue;
    TH2D *hPhiTrueDiff;

    TH2D *hMcsTrue;
    TH2D *hMcsTrueDiff;
    TH2D *hMcsTrueCut;
    TH2D *hMcsTrueDiffCut;
    TH2D *hMcsTrueAll;
    TH2D *hMcsTrueDiffAll;
    TH2D *hMcsTrueAllCut;
    TH2D *hMcsTrueDiffAllCut;

    TH2D *hRangeTrue;
    TH2D *hRangeTrueDiff;
    TH2D *hRangeTrueCut;
    TH2D *hRangeTrueDiffCut;

    TH2D *hBothTrue;
    TH2D *hBothTrueDiff;
    TH2D *hBothTrueCut;
    TH2D *hBothTrueDiffCut;

    TPCGeoAlg fTpcGeo;
    // Momentum fitters
    trkf::TrajectoryMCSFitter     fMcsFitter; 
    trkf::TrackMomentumCalculator fRangeFitter;

    int nNuMuCC = 0;
    int nMu = 0;
    int nMuIsNuTrack = 0;
    int nMuIsNuShower = 0;
    int nMuIsCrTrack = 0;
    int nMuNoReco = 0;
    int nMuIsMix = 0;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);

  }; // class PandoraNuMuRecoAna


  // Constructor
  PandoraNuMuRecoAna::PandoraNuMuRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().inputs().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().inputs().TpcTrackModuleLabel())
    , fHitModuleLabel       (config().inputs().HitModuleLabel())
    , fShowerModuleLabel    (config().inputs().ShowerModuleLabel())
    , fPandoraLabel         (config().inputs().PandoraLabel())
    , fVerbose              (config().inputs().Verbose())
    , fMcsFitter            (config().fitter)
  {

  } // PandoraNuMuRecoAna()


  void PandoraNuMuRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    for(size_t i = 0; i < nTC; i++){
      TString hLength_name  = Form("hTrueMuLength%s", trueCategories[i].c_str());
      hTrueMuLength[i]      = tfs->make<TH1D>(hLength_name,     "", 20, 0,    500);
      TString hMom_name     = Form("hTrueMuMom%s", trueCategories[i].c_str());
      hTrueMuMom[i]         = tfs->make<TH1D>(hMom_name,        "", 20, 0,    2);
      TString hTheta_name   = Form("hTrueMuTheta%s", trueCategories[i].c_str());
      hTrueMuTheta[i]       = tfs->make<TH1D>(hTheta_name,      "", 20, 0,    3.2);
      TString hPhi_name     = Form("hTrueMuPhi%s", trueCategories[i].c_str());
      hTrueMuPhi[i]         = tfs->make<TH1D>(hPhi_name,        "", 20, -3.2, 3.2);
    }

    hMuCompletenessEdep   = tfs->make<TH1D>("hMuCompletenessEdep",    "", 24, 0, 1.2);
    hMuCompletenessHitN   = tfs->make<TH1D>("hMuCompletenessHitN",    "", 24, 0, 1.2);
    hMuEdep               = tfs->make<TH1D>("hMuEdep",                "", 40, 0, 2);
    hMuHitN               = tfs->make<TH1D>("hMuHitN",                "", 40, 0, 2);

    hMuLengthDiff         = tfs->make<TH1D>("hMuLengthDiff",          "", 50, -50,  50);
    hMuThetaDiff          = tfs->make<TH1D>("hMuThetaDiff",           "", 50, -3.2, 3.2);
    hMuPhiDiff            = tfs->make<TH1D>("hMuPhiDiff",             "", 50, -3.2, 3.2);
    hMuAngle              = tfs->make<TH1D>("hMuAngle",               "", 50, 0,    3.2);

    hMuVertexDiff         = tfs->make<TH1D>("hMuVertexDiff",          "", 30, 0,   15);
    hMuVertexXDiff        = tfs->make<TH1D>("hMuVertexXDiff",         "", 60, -15, 15);
    hMuVertexYDiff        = tfs->make<TH1D>("hMuVertexYDiff",         "", 60, -15, 15);
    hMuVertexZDiff        = tfs->make<TH1D>("hMuVertexZDiff",         "", 60, -15, 15);

    hMcsDiff              = tfs->make<TH1D>("hMcsDiff",               "", 50, -1,  1);
    hMcsDiffCut           = tfs->make<TH1D>("hMcsDiffCut",            "", 50, -1,  1);
    hMcsDiffAll           = tfs->make<TH1D>("hMcsDiffAll",            "", 50, -1,  1);
    hMcsDiffAllCut        = tfs->make<TH1D>("hMcsDiffAllCut",         "", 50, -1,  1);
    hRangeDiff            = tfs->make<TH1D>("hRangeDiff",             "", 50, -1,  1);
    hRangeDiffCut         = tfs->make<TH1D>("hRangeDiffCut",          "", 50, -1,  1);
    hBothDiff             = tfs->make<TH1D>("hBothDiff",              "", 50, -1,  1);
    hBothDiffCut          = tfs->make<TH1D>("hBothDiffCut",           "", 50, -1,  1);

    hLengthTrue           = tfs->make<TH2D>("hLengthTrue",            "", 50, 0,    500, 50, 0,    500);
    hLengthTrueDiff       = tfs->make<TH2D>("hLengthTrueDiff",        "", 50, 0,    500, 50, -1,   1);
    hThetaTrue            = tfs->make<TH2D>("hThetaTrue",             "", 50, 0,    3.2, 50, 0,    3.2);
    hThetaTrueDiff        = tfs->make<TH2D>("hThetaTrueDiff",         "", 50, 0,    3.2, 50, -1,   1);
    hPhiTrue              = tfs->make<TH2D>("hPhiTrue",               "", 50, -3.2, 3.2, 50, -3.2, 3.2);
    hPhiTrueDiff          = tfs->make<TH2D>("hPhiTrueDiff",           "", 50, -3.2, 3.2, 50, -1,   1);

    hMcsTrue              = tfs->make<TH2D>("hMcsTrue",               "", 50, 0, 2, 50, 0,  2);
    hMcsTrueDiff          = tfs->make<TH2D>("hMcsTrueDiff",           "", 50, 0, 2, 50, -1, 1);
    hMcsTrueCut           = tfs->make<TH2D>("hMcsTrueCut",            "", 50, 0, 2, 50, 0,  2);
    hMcsTrueDiffCut       = tfs->make<TH2D>("hMcsTrueDiffCut",        "", 50, 0, 2, 50, -1, 1);
    hMcsTrueAll           = tfs->make<TH2D>("hMcsTrueAll",            "", 50, 0, 2, 50, 0,  2);
    hMcsTrueDiffAll       = tfs->make<TH2D>("hMcsTrueDiffAll",        "", 50, 0, 2, 50, -1, 1);
    hMcsTrueAllCut        = tfs->make<TH2D>("hMcsTrueAllCut",         "", 50, 0, 2, 50, 0,  2);
    hMcsTrueDiffAllCut    = tfs->make<TH2D>("hMcsTrueDiffAllCut",     "", 50, 0, 2, 50, -1, 1);

    hRangeTrue            = tfs->make<TH2D>("hRangeTrue",             "", 50, 0, 2, 50, 0,  2);
    hRangeTrueDiff        = tfs->make<TH2D>("hRangeTrueDiff",         "", 50, 0, 2, 50, -1, 1);
    hRangeTrueCut         = tfs->make<TH2D>("hRangeTrueCut",          "", 50, 0, 2, 50, 0,  2);
    hRangeTrueDiffCut     = tfs->make<TH2D>("hRangeTrueDiffCut",      "", 50, 0, 2, 50, -1, 1);

    hBothTrue             = tfs->make<TH2D>("hBothTrue",              "", 50, 0, 2, 50, 0,  2);
    hBothTrueDiff         = tfs->make<TH2D>("hBothTrueDiff",          "", 50, 0, 2, 50, -1, 1);
    hBothTrueCut          = tfs->make<TH2D>("hBothTrueCut",           "", 50, 0, 2, 50, 0,  2);
    hBothTrueDiffCut      = tfs->make<TH2D>("hBothTrueDiffCut",       "", 50, 0, 2, 50, -1, 1);

    // Initial output
    if(fVerbose) std::cout<<"----------------- Pandora NuMu Reco Ana Module -------------------"<<std::endl;

  }// PandoraNuMuRecoAna::beginJob()


  void PandoraNuMuRecoAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);

    if( !pfParticleHandle.isValid() ){
      std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }

    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

    std::vector<art::Ptr<recob::PFParticle>> crParticles;
    std::vector<art::Ptr<recob::PFParticle>> nuParticles;

    this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, nuParticles);

    std::vector<art::Ptr<recob::Track>> tracks;
    std::vector<art::Ptr<recob::Shower>> showers;
    this->CollectTracksAndShowers(nuParticles, pfParticleHandle, event, tracks, showers);

    std::vector<art::Ptr<recob::Track>> crTracks;
    std::vector<art::Ptr<recob::Shower>> crShowers;
    this->CollectTracksAndShowers(crParticles, pfParticleHandle, event, crTracks, crShowers);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<int> lepParticleIds;

    std::vector<double> usedNuVtx;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      simb::MCNeutrino mcNu = truth->GetNeutrino();
      bool isNuMuCC = false;
      if(std::abs(mcNu.Lepton().PdgCode())==13 && mcNu.CCNC() == simb::kCC) isNuMuCC = true;
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(mcNu.Nu().Vx()); vtx.SetY(mcNu.Nu().Vy()); vtx.SetZ(mcNu.Nu().Vz());

        if(!CosmicIdUtils::InFiducial(vtx, 0, 0)) continue;

        if(isNuMuCC && std::find(usedNuVtx.begin(), usedNuVtx.end(), vtx.X())==usedNuVtx.end()){ 
          nNuMuCC++;
          usedNuVtx.push_back(vtx.X());
        }

        if(std::abs(particle.PdgCode())==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partId);
          nMu++;
          hTrueMuLength[0]->Fill(fTpcGeo.TpcLength(particle));
          hTrueMuMom[0]->Fill(particle.P());
          std::pair<TVector3, TVector3> se = fTpcGeo.CrossingPoints(particle);
          double theta = (se.second-se.first).Theta();
          double phi = (se.second-se.first).Phi();
          hTrueMuTheta[0]->Fill(theta);
          hTrueMuPhi[0]->Fill(phi);
        }

      }
    }

    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    auto tpcShowerHandle = event.getValidHandle<std::vector<recob::Shower>>(fShowerModuleLabel);
    art::FindManyP<recob::Hit> findManyHitsShower(tpcShowerHandle, event, fShowerModuleLabel);

    auto hitHandle = event.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);

    std::vector<int> nuRecoTrackIds;
    std::vector<int> nuRecoShowerIds;
    std::vector<int> crRecoTrackIds;
    std::vector<int> crRecoShowerIds;

    for(auto const& track : tracks){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      nuRecoTrackIds.push_back(trueId);
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        simb::MCParticle part = particles[trueId];

        // Vertex differences
        std::pair<TVector3, TVector3> se = fTpcGeo.CrossingPoints(part);
        TVector3 trueVtx = se.first;
        TVector3 trueEnd = se.second;
        TVector3 recoVtx = track->Vertex<TVector3>();
        TVector3 recoEnd = track->End<TVector3>();
        hMuVertexDiff->Fill((recoVtx - trueVtx).Mag());
        hMuVertexXDiff->Fill(recoVtx.X() - trueVtx.X());
        hMuVertexYDiff->Fill(recoVtx.Y() - trueVtx.Y());
        hMuVertexZDiff->Fill(recoVtx.Z() - trueVtx.Z());

        // Length differences
        double trueLength = fTpcGeo.TpcLength(part);
        double recoLength = track->Length();
        hMuLengthDiff->Fill(recoLength - trueLength);
        hLengthTrue->Fill(trueLength, recoLength);
        hLengthTrueDiff->Fill(trueLength, (recoLength-trueLength)/trueLength);

        // Angle differences
        double trueTheta = (trueEnd - trueVtx).Theta();
        double truePhi = (trueEnd - trueVtx).Phi();
        double recoTheta = track->Theta();
        double recoPhi = track->Phi();
        hMuThetaDiff->Fill(TMath::ATan2(TMath::Sin(recoTheta-trueTheta), TMath::Cos(recoTheta-trueTheta)));
        hThetaTrue->Fill(trueTheta, recoTheta);
        hThetaTrueDiff->Fill(trueTheta, (recoTheta-trueTheta)/trueTheta);
        hMuPhiDiff->Fill(TMath::ATan2(TMath::Sin(recoPhi-truePhi), TMath::Cos(recoPhi-truePhi)));
        hPhiTrue->Fill(truePhi, recoPhi);
        hPhiTrueDiff->Fill(truePhi, (recoPhi-truePhi)/truePhi);
        hMuAngle->Fill((recoEnd - recoVtx).Angle((trueEnd - trueVtx)));

        // Completeness and contamination
        double trueHitN = 0;
        double trueEdep = 0;
        // Get all hits associated with true particle
        for (unsigned int i = 0; i < hitHandle->size(); ++i){
          const art::Ptr<recob::Hit> hit(hitHandle, i);
          if(RecoUtils::TrueParticleID(hit, false) == trueId){
            // Add up number and energy (area)
            trueHitN++;
            trueEdep += hit->Integral();
          }
        }
        double recoHitN = 0;
        double recoEdep = 0;
        double recoHitNTot = 0;
        double recoEdepTot = 0;
        // Get all hits associated with track
        for(auto const& hit : hits){
          if(RecoUtils::TrueParticleID(hit, false) == trueId){
            // Add up number and energy (area)
            recoHitN++;
            recoEdep += hit->Integral();
          }
          recoHitNTot++;
          recoEdepTot += hit->Integral();
        }
        hMuCompletenessHitN->Fill(recoHitN/trueHitN);
        hMuCompletenessEdep->Fill(recoEdep/trueEdep);
        hMuHitN->Fill(recoHitNTot/trueHitN);
        hMuEdep->Fill(recoEdepTot/trueEdep);

        // Momentum reconstruction
        double trueMom = part.P();
        bool exits = false;
        TVector3 outEnd (part.EndX(), part.EndY(), part.EndZ());
        if(trueEnd != outEnd) exits = true;
        double recoMom = 0.;
        if(!exits){
          recoMom = fRangeFitter.GetTrackMomentum(recoLength, 13);
          double diff = recoMom - trueMom;

          hRangeDiff->Fill(diff);
          hRangeTrue->Fill(trueMom, recoMom);
          hRangeTrueDiff->Fill(trueMom, diff/trueMom);

          hBothDiff->Fill(diff);
          hBothTrue->Fill(trueMom, recoMom);
          hBothTrueDiff->Fill(trueMom, diff/trueMom);

          if(recoLength > 50){
            hRangeTrueCut->Fill(trueMom, recoMom);
            hRangeTrueDiffCut->Fill(trueMom, diff/trueMom);
            hRangeDiffCut->Fill(diff);

            hBothTrueCut->Fill(trueMom, recoMom);
            hBothTrueDiffCut->Fill(trueMom, diff/trueMom);
            hBothDiffCut->Fill(diff);
          }

          recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(*track);
          double mcsMom = mcsResult.bestMomentum();
          double mcsDiff = mcsMom - trueMom;

          hMcsTrueAll->Fill(trueMom, mcsMom);
          hMcsTrueDiffAll->Fill(trueMom, mcsDiff/trueMom);
          hMcsDiffAll->Fill(mcsDiff);

          if(recoLength > 100){
            hMcsTrueAllCut->Fill(trueMom, mcsMom);
            hMcsTrueDiffAllCut->Fill(trueMom, mcsDiff/trueMom);
            hMcsDiffAllCut->Fill(mcsDiff);
          }
        }

        else{
          recob::MCSFitResult mcsResult = fMcsFitter.fitMcs(*track);
          recoMom = mcsResult.bestMomentum();
          double diff = recoMom - trueMom;

          hMcsTrue->Fill(trueMom, recoMom);
          hMcsTrueDiff->Fill(trueMom, diff/trueMom);
          hMcsDiff->Fill(diff);

          hMcsTrueAll->Fill(trueMom, recoMom);
          hMcsTrueDiffAll->Fill(trueMom, diff/trueMom);
          hMcsDiffAll->Fill(diff);

          hBothTrue->Fill(trueMom, recoMom);
          hBothTrueDiff->Fill(trueMom, diff/trueMom);
          hBothDiff->Fill(diff);

          if(recoLength > 100){
            hMcsTrueCut->Fill(trueMom, recoMom);
            hMcsTrueDiffCut->Fill(trueMom, diff/trueMom);
            hMcsDiffCut->Fill(diff);

            hMcsTrueAllCut->Fill(trueMom, recoMom);
            hMcsTrueDiffAllCut->Fill(trueMom, diff/trueMom);
            hMcsDiffAllCut->Fill(diff);

            hBothTrueCut->Fill(trueMom, recoMom);
            hBothTrueDiffCut->Fill(trueMom, diff/trueMom);
            hBothDiffCut->Fill(diff);
          }
        }
      }
    }

    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);
      std::vector<art::Ptr<recob::Hit>> hits = findManyHitsShower.at(shower.key());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      nuRecoShowerIds.push_back(trueId);
    }

    for(auto const& crTrack : crTracks){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(crTrack->ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      crRecoTrackIds.push_back(trueId);
    }

    for(size_t i = 0; i < lepParticleIds.size(); i++){
      int lepId = lepParticleIds[i];
      std::cout<<"True particle "<<lepId<<":\n";
      simb::MCParticle particle = particles[lepId];
      double length = fTpcGeo.TpcLength(particle);
      double mom = particle.P();
      std::pair<TVector3, TVector3> se = fTpcGeo.CrossingPoints(particle);
      double theta = (se.second-se.first).Theta();
      double phi = (se.second-se.first).Phi();

      int trueType = 5;
      if(std::find(nuRecoTrackIds.begin(), nuRecoTrackIds.end(), lepId) != nuRecoTrackIds.end()){
        if(fVerbose) std::cout<<"ID AS NU TRACK\n";
        if(trueType == 5) trueType = 1;
        else trueType = 4;
      }
      if(std::find(nuRecoShowerIds.begin(), nuRecoShowerIds.end(), lepId) != nuRecoShowerIds.end()){
        if(fVerbose) std::cout<<"ID AS NU SHOWER\n";
        if(trueType == 5) trueType = 2;
        else trueType = 4;
      }
      if(std::find(crRecoTrackIds.begin(), crRecoTrackIds.end(), lepId) != crRecoTrackIds.end()){
        if(fVerbose) std::cout<<"ID AS CR TRACK\n";
        if(trueType == 5) trueType = 3;
        else trueType = 4;
      }

      hTrueMuLength[trueType]->Fill(length);
      hTrueMuMom[trueType]->Fill(mom);
      hTrueMuTheta[trueType]->Fill(theta);
      hTrueMuPhi[trueType]->Fill(phi);

      if(trueType == 5) nMuNoReco++;
      else if(trueType == 1) nMuIsNuTrack++;
      else if(trueType == 2) nMuIsNuShower++;
      else if(trueType == 3) nMuIsCrTrack++;
      else if(trueType == 4) nMuIsMix++;
    }

  } // PandoraNuMuRecoAna::analyze()


  void PandoraNuMuRecoAna::endJob(){

    std::cout<<"Number of NuMuCC interactions in TPC = "<<nNuMuCC<<"\n"
             <<"Total true mu particles     = "<<nMu<<"\n"
             <<"Reco as nu track            = "<<nMuIsNuTrack<<"\n"
             <<"Reco as nu shower           = "<<nMuIsNuShower<<"\n"
             <<"Reco as cr track            = "<<nMuIsCrTrack<<"\n"
             <<"Reco as mix                 = "<<nMuIsMix<<"\n"
             <<"Not reconstructed           = "<<nMuNoReco<<"\n";

  } // PandoraNuMuRecoAna::endJob()

  void PandoraNuMuRecoAna::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
      
  void PandoraNuMuRecoAna::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles){

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
      
  void PandoraNuMuRecoAna::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
  {
      // Get the associations between PFParticles and tracks/showers from the event
      art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, fTpcTrackModuleLabel);
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


  DEFINE_ART_MODULE(PandoraNuMuRecoAna)
} // namespace sbnd

