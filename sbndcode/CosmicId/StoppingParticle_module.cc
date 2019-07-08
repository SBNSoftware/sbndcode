////////////////////////////////////////////////////////////////////////
// Class:       StoppingParticle
// Module Type: analyzer
// File:        StoppingParticleAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

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

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TGraphAsymmErrors.h"

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

  class StoppingParticle : public art::EDAnalyzer {
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
      
      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of calorimetry producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<int> TrackID {
        Name("TrackID"),
        Comment("Track ID")
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit StoppingParticle(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
    bool          fVerbose;             ///< print information about what's going on
    int           fTrackID;

    // histograms
    // True time of neutrino particles to get beam window
    TH1D* hStoppingNuMu;
    TH1D* hStoppingCos;
    TH1D* hSlopeNuMu;
    TH1D* hSlopeCos;
    TH1D* hChi2NuMu;
    TH1D* hChi2Cos;
    TH1D* hChi2NuMu5;
    TH1D* hChi2Cos5;
    TH1D* hChi2NuMu1;
    TH1D* hChi2Cos1;
    TH1D* hChi2NuMu_15;
    TH1D* hChi2Cos_15;
    TH1D* hChi2NuMu_10;
    TH1D* hChi2Cos_10;
    TH1D* hPolChi2NuMu;
    TH1D* hPolChi2Cos;
    TH1D* hAveNuMu;
    TH1D* hAveCos;
    //2d hist of slope and chi2

    // Other variables shared between different methods.
    CRTTrackMatchAlg trackAlg;

    std::map<int, std::vector<double>> chi2numu;
    std::map<int, std::vector<double>> chi2cos;

    // Performance Counters
    int nCosmicTracks = 0;
    int nCosmicTracksVal = 0;
    int nLepTracks = 0;
    int nLepTracksVal = 0;
    int nSlopeNuMu = 0;
    int nSlopeCos = 0;
    int nChi2NuMu = 0;
    int nChi2Cos = 0;
    int nAveNuMu = 0;
    int nAveCos = 0;
    int nNuMuRemoved = 0;
    int nCosRemoved = 0;

  }; // class StoppingParticle


  // Constructor
  StoppingParticle::StoppingParticle(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fVerbose              (config().Verbose())
    , fTrackID              (config().TrackID())
  {


  } // StoppingParticle()


  void StoppingParticle::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    hStoppingNuMu = tfs->make<TH1D>("hStoppingNuMu", "", 10, 0, 10);
    hStoppingCos = tfs->make<TH1D>("hStoppingCos", "", 10, 0, 10);
    hSlopeNuMu = tfs->make<TH1D>("hSlopeNuMu", "", 40, -0.5, 0.5);
    hSlopeCos = tfs->make<TH1D>("hSlopeCos", "", 40, -0.5, 0.5);
    hChi2NuMu = tfs->make<TH1D>("hChi2NuMu", "", 50, 0, 2.5);
    hChi2Cos = tfs->make<TH1D>("hChi2Cos", "", 50, 0, 2.5);
    hChi2NuMu5 = tfs->make<TH1D>("hChi2NuMu5", "", 50, 0, 2.5);
    hChi2Cos5 = tfs->make<TH1D>("hChi2Cos5", "", 50, 0, 2.5);
    hChi2NuMu1 = tfs->make<TH1D>("hChi2NuMu1", "", 50, 0, 2.5);
    hChi2Cos1 = tfs->make<TH1D>("hChi2Cos1", "", 50, 0, 2.5);
    hChi2NuMu_15 = tfs->make<TH1D>("hChi2NuMu_15", "", 50, 0, 2.5);
    hChi2Cos_15 = tfs->make<TH1D>("hChi2Cos_15", "", 50, 0, 2.5);
    hChi2NuMu_10 = tfs->make<TH1D>("hChi2NuMu_10", "", 50, 0, 2.5);
    hChi2Cos_10 = tfs->make<TH1D>("hChi2Cos_10", "", 50, 0, 2.5);
    hPolChi2NuMu = tfs->make<TH1D>("hPolChi2NuMu", "", 50, 0, 2.5);
    hPolChi2Cos = tfs->make<TH1D>("hPolChi2Cos", "", 50, 0, 2.5);
    hAveNuMu = tfs->make<TH1D>("hAveNuMu", "", 40, 0, 5);
    hAveCos = tfs->make<TH1D>("hAveCos", "", 40, 0, 5);

    // Initial output
    if(fVerbose) std::cout<<"----------------- TPC Cosmic Removal Ana Module -------------------"<<std::endl;

  }// StoppingParticle::beginJob()


  void StoppingParticle::analyze(const art::Event& event)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

    // Loop over true particles in readout window
    std::map<int, simb::MCParticle> particles;
    std::vector<int> cosParticleIds;
    std::vector<int> lepParticleIds;

    for (auto const& particle: (*particleHandle)){

      int partId = particle.TrackId();
      particles[partId] = particle;
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partId);
      geo::Point_t start {particle.Vx(), particle.Vy(), particle.Vz()};
      bool sif = CosmicIdUtils::InFiducial(start, 0, 0);
      geo::Point_t end {particle.EndX(), particle.EndY(), particle.EndZ()};
      bool eif = CosmicIdUtils::InFiducial(end, 0, 0);
      if((sif && !eif) || (eif && !sif)){
        if(truth->Origin() == simb::kBeamNeutrino){
          geo::Point_t vtx;
          vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
          double time = particle.T() * 1e-3; // [us]
          if(fVerbose && particle.Mother()==0) 
            std::cout<<"Nu VTX = "<<vtx<<" ID = "<<partId<<" pdg = "<<particle.PdgCode()<<" time = "<<time<<" length = "<<particle.Trajectory().TotalLength()
                     <<" start = ("<<particle.Vx()<<", "<<particle.Vy()<<", "<<particle.Vz()<<") end = ("<<particle.EndX()<<", "<<particle.EndY()<<", "<<particle.EndZ()<<")\n";
          if(!CosmicIdUtils::InFiducial(vtx, 0, 0)) continue;
          if(std::abs(particle.PdgCode())==13 && particle.Mother()==0) lepParticleIds.push_back(partId);
       
        }
        else if(particle.Trajectory().TotalLength() > 500){
          cosParticleIds.push_back(partId);
        }
      }

    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<"\n\n";
    int nStopLep = 0;
    int nStopCos = 0;

    for(auto const& tpcTrack : (*tpcTrackHandle)){

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if (particles.find(trueId) == particles.end()){ 
        continue; 
      }

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());
      if(calos.size()==0) continue;

      geo::Point_t trueStart {particles[trueId].Vx(), particles[trueId].Vy(), particles[trueId].Vz()};
      geo::Point_t trueEnd {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
      geo::Point_t end = tpcTrack.End();

      if(CosmicIdUtils::InFiducial(trueStart, 0., 0.)){
        if(std::sqrt(std::pow(trueStart.Y()-tpcTrack.Vertex().Y(),2)+std::pow(trueStart.Z()-tpcTrack.Vertex().Z(),2)) > 20 && std::sqrt(std::pow(trueStart.Y()-tpcTrack.End().Y(),2)+std::pow(trueStart.Z()-tpcTrack.End().Z(),2)) > 20) continue;
        if(std::sqrt(std::pow(trueStart.Y()-tpcTrack.Vertex().Y(),2)+std::pow(trueStart.Z()-tpcTrack.Vertex().Z(),2)) < std::sqrt(std::pow(trueStart.Y()-tpcTrack.End().Y(),2)+std::pow(trueStart.Z()-tpcTrack.End().Z(),2))) end = tpcTrack.Vertex();
      }
      if(CosmicIdUtils::InFiducial(trueEnd, 0., 0.)){
        if(std::sqrt(std::pow(trueEnd.Y()-tpcTrack.Vertex().Y(),2)+std::pow(trueEnd.Z()-tpcTrack.Vertex().Z(),2)) > 20 && std::sqrt(std::pow(trueEnd.Y()-tpcTrack.End().Y(),2)+std::pow(trueEnd.Z()-tpcTrack.End().Z(),2)) > 20) continue;
        if(std::sqrt(std::pow(trueEnd.Y()-tpcTrack.Vertex().Y(),2)+std::pow(trueEnd.Z()-tpcTrack.Vertex().Z(),2)) < std::sqrt(std::pow(trueEnd.Y()-tpcTrack.End().Y(),2)+std::pow(trueEnd.Z()-tpcTrack.End().Z(),2))) end = tpcTrack.Vertex();
      }

      // Is track from a neutrino interaction
      bool isLep = false;
      bool isCos = false;
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){
        // Set neutrino flag
        isLep = true;
        nLepTracks++;
        nStopLep++;
        if(fVerbose) std::cout<<"\nPrimary lepton!\n";
      }
      // Else is it from a primary cosmic (true length > 500)
      if(std::find(cosParticleIds.begin(), cosParticleIds.end(), trueId) != cosParticleIds.end()){
        // Set primary cosmic flag
        isCos = true;
        nCosmicTracks++;
        nStopCos++;
        if(fVerbose) std::cout<<"\nPrimary cosmic!\n";
      }
      if(!isCos && !isLep) continue;
      if(fVerbose) std::cout<<"Track ID:"<<tpcTrack.ID()<<":\n";

      if(fVerbose) std::cout<<"calo size = "<<calos[0]->dEdx().size()<<" "<<calos[1]->dEdx().size()<<" "<<calos[2]->dEdx().size()<<"\n";
      size_t nhits = 0;
      art::Ptr<anab::Calorimetry> calo = calos[0];
      for( size_t i = calos.size(); i > 0; i--){
        if(calos[i-1]->dEdx().size() > nhits*1.5){
          nhits = calos[i-1]->dEdx().size();
          calo = calos[i-1];
        }
      }

      double distStart = (calo->XYZ()[0] - end).Mag2();
      double distEnd = (calo->XYZ()[nhits-1] - end).Mag2();

      std::vector<double> v_resrg;
      std::vector<double> v_dedx;
      std::vector<double> v_resrg5;
      std::vector<double> v_dedx5;
      std::vector<double> v_resrg1;
      std::vector<double> v_dedx1;

      std::vector<double> v_resrgFull;
      std::vector<double> v_dedxFull;

      double maxResRg = std::max(calo->ResidualRange()[0], calo->ResidualRange()[nhits-1]);
      double maxBegin = maxResRg * 0.2;
      if(maxResRg>20) maxBegin = 10.;
      double dedxBegin = 0;
      int nptsBegin = 0;
      double minEnd = maxResRg - maxBegin;
      double dedxEnd = 0;
      int nptsEnd = 0;

      double maxDedx = 0;
      double resrgStart = 0;
      double maxDedx5 = 0;
      double resrgStart5 = 0;
      double maxDedx1 = 0;
      double resrgStart1 = 0;

      for(size_t i = 0; i < nhits; i++){
        double dedx = calo->dEdx()[i];
        double resrg = calo->ResidualRange()[i];
        
        if(dedx < 100){
          v_resrgFull.push_back(resrg);
          v_dedxFull.push_back(dedx);
        }

        if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
        if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

        if(resrg < 10 && dedx > maxDedx && dedx < 30){
          maxDedx = dedx;
          resrgStart = resrg;
        }
        if(resrg < 5 && dedx > maxDedx5 && dedx < 30){
          maxDedx5 = dedx;
          resrgStart5 = resrg;
        }
        if(resrg < 1 && dedx > maxDedx1 && dedx < 30){
          maxDedx1 = dedx;
          resrgStart1 = resrg;
        }

        if(resrg > 0.4 && resrg < maxBegin){
          dedxBegin += dedx;
          nptsBegin++;
        }

        if(resrg > minEnd){
          dedxEnd += dedx;
          nptsEnd++;
        }

      }

      for(size_t i = 0; i < nhits; i++){
        double dedx = calo->dEdx()[i];
        double resrg = calo->ResidualRange()[i];

        if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
        if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

        if(resrg > resrgStart && resrg < resrgStart+30 && dedx < 30){
          v_resrg.push_back(resrg);
          v_dedx.push_back(dedx);
        }
        if(resrg > resrgStart5 && resrg < resrgStart5+30 && dedx < 30){
          v_resrg5.push_back(resrg);
          v_dedx5.push_back(dedx);
        }
        if(resrg >  resrgStart1 && resrg < resrgStart1+30 && dedx < 30){
          v_resrg1.push_back(resrg);
          v_dedx1.push_back(dedx);
        }
      }

      dedxBegin = dedxBegin/nptsBegin;
      dedxEnd = dedxEnd/nptsEnd;

      if(v_dedx.size() < 10){if(fVerbose) std::cout<<"size = "<<v_dedx.size()<<" full size = "<<v_dedxFull.size()<<" Too small\n"; continue;}

      TGraph *gdedx = new TGraph(v_dedx.size(), &v_resrg[0], &v_dedx[0]);
      TGraph *gdedxFull = new TGraph(v_dedxFull.size(), &v_resrgFull[0], &v_dedxFull[0]);

      try{ gdedx->Fit("pol1", "Q"); } catch(...){ std::cout<<"pol1 fit failed\n"; continue; }
      TF1* polfit = gdedx->GetFunction("pol1");
      double slope = polfit->GetParameter(1);
      double chi2 = polfit->GetChisquare();

      try{ gdedx->Fit("pol0", "Q", "", 0, resrgStart+20); } catch(...){ if(fVerbose) std::cout<<"pol0 fit failed\n"; continue; }
      TF1* pol0fit = gdedx->GetFunction("pol0");
      double pol0chi2 = pol0fit->GetChisquare();

      try{ gdedx->Fit("expo", "Q", "", 0, resrgStart+20); } catch(...){ if(fVerbose) std::cout<<"Exp fit failed\n"; continue; }
      TF1* expfit = gdedx->GetFunction("expo");
      double expchi2 = expfit->GetChisquare();

      TGraph *gdedx5 = new TGraph(v_dedx5.size(), &v_resrg5[0], &v_dedx5[0]);
      try{ gdedx5->Fit("pol0", "Q", "", 0, resrgStart5+20); } catch(...){ if(fVerbose) std::cout<<"pol0 fit failed\n"; continue; }
      TF1* pol0fit5 = gdedx5->GetFunction("pol0");
      double pol0chi25 = pol0fit5->GetChisquare();
      try{ gdedx5->Fit("expo", "Q", "", 0, resrgStart5+20); } catch(...){ if(fVerbose) std::cout<<"Exp fit failed\n"; continue; }
      TF1* expfit5 = gdedx5->GetFunction("expo");
      double expchi25 = expfit5->GetChisquare();

      try{ gdedx5->Fit("pol0", "Q", "", 0, resrgStart5+15); } catch(...){ if(fVerbose) std::cout<<"pol0 fit failed\n"; continue; }
      TF1* pol0fit_15 = gdedx5->GetFunction("pol0");
      double pol0chi2_15 = pol0fit_15->GetChisquare();
      try{ gdedx5->Fit("expo", "Q", "", 0, resrgStart5+15); } catch(...){ if(fVerbose) std::cout<<"Exp fit failed\n"; continue; }
      TF1* expfit_15 = gdedx5->GetFunction("expo");
      double expchi2_15 = expfit_15->GetChisquare();

      try{ gdedx5->Fit("pol0", "Q", "", 0, resrgStart5+10); } catch(...){ if(fVerbose) std::cout<<"pol0 fit failed\n"; continue; }
      TF1* pol0fit_10 = gdedx5->GetFunction("pol0");
      double pol0chi2_10 = pol0fit_10->GetChisquare();
      try{ gdedx5->Fit("expo", "Q", "", 0, resrgStart5+10); } catch(...){ if(fVerbose) std::cout<<"Exp fit failed\n"; continue; }
      TF1* expfit_10 = gdedx5->GetFunction("expo");
      double expchi2_10 = expfit_10->GetChisquare();

      TGraph *gdedx1 = new TGraph(v_dedx1.size(), &v_resrg1[0], &v_dedx1[0]);
      try{ gdedx1->Fit("pol0", "Q", "", 0, resrgStart1+20); } catch(...){ if(fVerbose) std::cout<<"pol0 fit failed\n"; continue; }
      TF1* pol0fit1 = gdedx1->GetFunction("pol0");
      double pol0chi21 = pol0fit1->GetChisquare();
      try{ gdedx1->Fit("expo", "Q", "", 0, resrgStart1+20); } catch(...){ if(fVerbose) std::cout<<"Exp fit failed\n"; continue; }
      TF1* expfit1 = gdedx1->GetFunction("expo");
      double expchi21 = expfit1->GetChisquare();

      if(isCos) nCosmicTracksVal++;
      if(isLep) nLepTracksVal++;

      if(tpcTrack.ID()==fTrackID){
        TCanvas *c3 = new TCanvas("c3","",700,700);
        gdedx->SetMarkerStyle(3);
        gdedx->Draw("ap");
        c3->SaveAs("gdedx.root");
        TCanvas *c4 = new TCanvas("c4","",700,700);
        gdedxFull->SetMarkerStyle(3);
        gdedxFull->Draw("ap");
        c4->SaveAs("gdedxFull.root");
      }

      if(fVerbose) std::cout<<"resrg size = "<<v_resrg.size()<<" dedx size = "<<v_dedx.size()<<" slope = "<<slope<<" chi2 = "<<chi2<<" ratio = "<<pol0chi2/chi2<<"\n";
      if(fVerbose) std::cout<<"dedx begin = "<<dedxBegin<<" dedx end = "<<dedxEnd<<" exp chi2 = "<<expchi2<<" pol0chi2 = "<<pol0chi2<<" ratio = "<<pol0chi2/expchi2<<"\n";
/*
      double lim = 6;
      for(int i = 0; i < 10; i++){
        //try{ gdedx->Fit("pol0", "Q", "", 0, resrgStart+lim); } catch(...){ lim+=2; continue; }
        try{ gdedx->Fit("pol0", "Q", "", 0, resrgStart+7); } catch(...){ std::cout<<"pol0 fit failed\n"; continue; }
        TF1* pol0fitlim = gdedx->GetFunction("pol0");
        double pol0chi2lim = pol0fitlim->GetChisquare();
        //try{ gdedx->Fit("expo", "Q", "", 0, resrgStart+lim); } catch(...){ lim += 2; continue; }
      try{ gdedx->Fit("expo", "Q", "", 0, resrgStart+7); } catch(...){ std::cout<<"Exp fit failed\n"; continue; }
        TF1* expfitlim = gdedx->GetFunction("expo");
        double expchi2lim = expfitlim->GetChisquare();
        if(isLep) chi2numu[lim].push_back(pol0chi2lim/expchi2lim);
        if(isCos) chi2cos[lim].push_back(pol0chi2lim/expchi2lim);
        lim += 2;
      }
*/
      if(isLep){
        hSlopeNuMu->Fill(slope);
        hChi2NuMu->Fill(pol0chi2/expchi2);
        hChi2NuMu5->Fill(pol0chi25/expchi25);
        hChi2NuMu1->Fill(pol0chi21/expchi21);
        hChi2NuMu_15->Fill(pol0chi2_15/expchi2_15);
        hChi2NuMu_10->Fill(pol0chi2_10/expchi2_10);
        hPolChi2NuMu->Fill(pol0chi2/chi2);
        hAveNuMu->Fill(dedxBegin/dedxEnd);
        if(slope < -0.1){ if(fVerbose) std::cout<<"REMOVED: slope cut\n"; nSlopeNuMu++; }
        if((pol0chi2)/expchi2>1.2){ if(fVerbose) std::cout<<"REMOVED: chi2 cut\n"; nChi2NuMu++; }
        if((dedxBegin/dedxEnd) > 1.4){ if(fVerbose) std::cout<<"REMOVED: ave cut\n"; nAveNuMu++; }
        if((dedxBegin/dedxEnd) > 1.4||(pol0chi2/expchi2)>1.3) nNuMuRemoved++;
      }

      if(isCos){
        hSlopeCos->Fill(slope);
        hChi2Cos->Fill(pol0chi2/expchi2);
        hChi2Cos5->Fill(pol0chi25/expchi25);
        hChi2Cos1->Fill(pol0chi21/expchi21);
        hChi2Cos_10->Fill(pol0chi2_10/expchi2_10);
        hChi2Cos_15->Fill(pol0chi2_15/expchi2_15);
        hPolChi2Cos->Fill(pol0chi2/chi2);
        hAveCos->Fill(dedxBegin/dedxEnd);
        if(slope < -0.1){ if(fVerbose) std::cout<<"REMOVED: slope cut\n"; nSlopeCos++; }
        if((pol0chi2)/expchi2>1.2){ if(fVerbose) std::cout<<"REMOVED: chi2 cut\n"; nChi2Cos++; }
        if((dedxBegin/dedxEnd) > 1.4){ if(fVerbose) std::cout<<"REMOVED: ave cut\n"; nAveCos++; }
        if((dedxBegin/dedxEnd) > 1.4||(pol0chi2/expchi2)>1.3) nCosRemoved++;
      }

    }

    hStoppingNuMu->Fill(nStopLep);
    hStoppingCos->Fill(nStopCos);

  } // StoppingParticle::analyze()


  void StoppingParticle::endJob(){
/*
    double lim = 6;
    for(int i = 0; i < 10; i++){
      std::vector<double> v_chi2numu = chi2numu[lim];
      std::sort(v_chi2numu.begin(), v_chi2numu.end());
      size_t numusize = v_chi2numu.size();
      int index = std::floor((double)numusize*0.9);
      double chi2lim = v_chi2numu[index];
      std::vector<double> v_chi2cos = chi2cos[lim];
      std::sort(v_chi2cos.begin(), v_chi2cos.end());
      size_t cossize = v_chi2cos.size();
      size_t nremain = 0;
      bool stop = false;
      for(size_t j = 0; j < cossize; j++){
        if(v_chi2cos[j] > chi2lim && !stop){
          nremain = j;
          stop = true;
        }
      }
      std::cout<<"Limit = "<<lim<<" chi2 limit = "<<chi2lim<<" n numu = "<<numusize<<" remaining = "<<index<<" n cos = "<<cossize<<" remaining = "<<nremain<<"\n";
      lim += 2;
    }
*/
    std::cout<<"Total lepton tracks       = "<<nLepTracks<<"\n"
             <<"Valid tracks              = "<<nLepTracksVal<<"\n"
             <<"Removed by slope cut      = "<<nSlopeNuMu<<"\n"
             <<"Removed by chi2 cut       = "<<nChi2NuMu<<"\n"
             <<"Removed by ave cut        = "<<nAveNuMu<<"\n"
             <<"Total removed             = "<<nNuMuRemoved<<"\n"
             <<"-----------------------------------------------------------\n"
             <<"Total cosmic tracks       = "<<nCosmicTracks<<"\n"
             <<"Valid tracks              = "<<nCosmicTracksVal<<"\n"
             <<"Removed by slope cut      = "<<nSlopeCos<<"\n"
             <<"Removed by chi2 cut       = "<<nChi2Cos<<"\n"
             <<"Removed by ave cut        = "<<nAveCos<<"\n"
             <<"Total removed             = "<<nCosRemoved<<"\n";

  } // StoppingParticle::endJob()

  DEFINE_ART_MODULE(StoppingParticle)
} // namespace sbnd

