////////////////////////////////////////////////////////////////////////
// Class:       CRTFullRecoAna
// Module Type: analyzer
// File:        CRTFullRecoAna_module.cc
//
// Analysis module for evaluating CRT reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
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

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TF1.h"
#include "TString.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace sbnd {

  struct CRTTruthMatch{
    int partID;
    double trueT0;
    std::map<std::string, geo::Point_t> trueCrosses;
    std::map<std::string, bool> validCrosses;
    std::map<std::string, std::vector<crt::CRTHit>> crtHits;
    std::vector<crt::CRTTrack> crtTracks;
    bool hasTpcTrack;
    std::vector<double> hitT0s;
    std::vector<double> trackT0s;
  };

  class CRTFullRecoAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CRTSimLabel {
        Name("CRTSimLabel"),
        Comment("tag of CRT simulation data product")
      };

      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel"),
        Comment("tag of CRT simulation data product")
      };

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel"),
        Comment("tag of CRT simulation data product")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of TPC track data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTT0MatchAlg::Config> CRTT0Alg {
        Name("CRTT0Alg"),
      };

      fhicl::Table<CRTTrackMatchAlg::Config> CRTTrackAlg {
        Name("CRTTrackAlg"),
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTFullRecoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTSimLabel;         ///< name of CRT producer
    art::InputTag fCRTHitLabel;         ///< name of CRT producer
    art::InputTag fCRTTrackLabel;       ///< name of CRT producer
    art::InputTag fTPCTrackLabel;       ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on
    
    // Histograms
    TH1D* hCRTHitTimes;
    std::map<std::string, TH3D*> hTaggerXYZResolution;

    std::vector<std::string> stage{"CRTHit", "CRTTrack", "HitT0", "TrackT0"};
    std::vector<std::string> category{"PossMatch", "ShouldMatch", "CRTVol", "TPCVol", "CRTCross", "TPCCross"};
    std::vector<std::string> level{"Total", "Matched"};
    TH1D* hEffMom[4][6][2];
    TH1D* hEffTheta[4][6][2];
    TH1D* hEffPhi[4][6][2];

    std::vector<std::string> trackType{"Complete", "Incomplete", "Both"};
    TH1D* hTrackThetaDiff[3];
    TH1D* hTrackPhiDiff[3];
    TH2D* hTrackTheta[3];
    TH2D* hTrackPhi[3];

    std::vector<std::string> t0Type{"HitT0", "TrackT0"};
    std::vector<std::string> purity{"Matched", "Correct"};
    TH1D* hPurityMom[2][2];
    TH1D* hPurityTheta[2][2];
    TH1D* hPurityPhi[2][2];

    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

    CRTT0MatchAlg crtT0Alg;
    CRTTrackMatchAlg crtTrackAlg;

    CRTBackTracker fCrtBackTrack;

  }; // class CRTFullRecoAna

  // Constructor
  CRTFullRecoAna::CRTFullRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTSimLabel          (config().CRTSimLabel())
    , fCRTHitLabel          (config().CRTHitLabel())
    , fCRTTrackLabel        (config().CRTTrackLabel())
    , fTPCTrackLabel        (config().TPCTrackLabel())
    , fVerbose              (config().Verbose())
    , crtT0Alg              (config().CRTT0Alg())
    , crtTrackAlg           (config().CRTTrackAlg())
    , fCrtBackTrack         (config().CrtBackTrack())
  {
  }

  void CRTFullRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    hCRTHitTimes = tfs->make<TH1D>("CRTHitTimes","",100, -5000, 5000);

    for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
      std::string taggerName = fCrtGeo.GetTagger(i).name;
      TString histName = Form("%s_resolution", taggerName.c_str());
      hTaggerXYZResolution[taggerName] = tfs->make<TH3D>(histName, "", 50, -25, 25, 50, -25, 25, 50, -25, 25);
    }

    for(size_t si = 0; si < stage.size(); si++){
      for(size_t ci = 0; ci < category.size(); ci++){
        for(size_t li = 0; li < level.size(); li++){
          TString momName = Form("MuMom%s_%s_%s", stage[si].c_str(), category[ci].c_str(), level[li].c_str());
          hEffMom[si][ci][li] = tfs->make<TH1D>(momName, "", 20, 0, 10);
          TString thetaName = Form("MuTheta%s_%s_%s", stage[si].c_str(), category[ci].c_str(), level[li].c_str());
          hEffTheta[si][ci][li] = tfs->make<TH1D>(thetaName, "", 20, 0, 180);
          TString phiName = Form("MuPhi%s_%s_%s", stage[si].c_str(), category[ci].c_str(), level[li].c_str());
          hEffPhi[si][ci][li] = tfs->make<TH1D>(phiName, "", 20, -180, 180);
        }
      }
    }

    for(size_t ti = 0; ti < trackType.size(); ti++){
      TString thetaDiffName = Form("TrackThetaDiff%s", trackType[ti].c_str());
      hTrackThetaDiff[ti] = tfs->make<TH1D>(thetaDiffName, "", 50, -0.2, 0.2);
      TString phiDiffName = Form("TrackPhiDiff%s", trackType[ti].c_str());
      hTrackPhiDiff[ti] = tfs->make<TH1D>(phiDiffName, "", 50, -0.2, 0.2);
      TString thetaName = Form("TrackTheta%s", trackType[ti].c_str());
      hTrackTheta[ti] = tfs->make<TH2D>(thetaName, "", 40, 0, 180, 40, 0, 180);
      TString phiName = Form("TrackPhi%s", trackType[ti].c_str());
      hTrackPhi[ti] = tfs->make<TH2D>(phiName, "", 40, -180, 0, 40, -180, 0);
    }

    for(size_t ti = 0; ti < t0Type.size(); ti++){
      for(size_t pi = 0; pi < purity.size(); pi++){
        TString momName = Form("PurityMom%s_%s", t0Type[ti].c_str(), purity[pi].c_str());
        hPurityMom[ti][pi] = tfs->make<TH1D>(momName, "", 20, 0, 10);
        TString thetaName = Form("PurityTheta%s_%s", t0Type[ti].c_str(), purity[pi].c_str());
        hPurityTheta[ti][pi] = tfs->make<TH1D>(thetaName, "", 20, 0, 180);
        TString phiName = Form("PurityPhi%s_%s", t0Type[ti].c_str(), purity[pi].c_str());
        hPurityPhi[ti][pi] = tfs->make<TH1D>(phiName, "", 20, -180, 180);
      }
    }

    // Initial output
    std::cout<<"----------------- Full CRT Reconstruction Analysis Module -------------------"<<std::endl;

  }// CRTFullRecoAna::beginJob()

  void CRTFullRecoAna::analyze(const art::Event& event)
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
    art::Handle< std::vector<crt::CRTHit>> crtHitHandle;
    std::vector<art::Ptr<crt::CRTHit> > crtHitList;
    if (event.getByLabel(fCRTHitLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitList, crtHitHandle);

    // Get CRT tracks from the event
    art::Handle< std::vector<crt::CRTTrack>> crtTrackHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackHandle);
    art::FindManyP<crt::CRTHit> findManyCrtHits(crtTrackHandle, event, fCRTTrackLabel);

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    fCrtBackTrack.Initialize(event);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
     
    // Only evaluate performance over the time window that hits have been reconstructed
    double minHitTime = 99999;
    double maxHitTime = -99999;
    for(auto const& hit : (*crtHitHandle)){
      double hitTime = (double)(int)hit.ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;
    }

    std::map<int, simb::MCParticle> particles;
    std::map<int, CRTTruthMatch> truthMatching;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      CRTTruthMatch truthMatch;
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
      
      // Only consider primary muons
      if(!(std::abs(particle.PdgCode()) == 13 && particle.Mother()==0)) continue;

      double time = particle.T() * 1e-3;
      if(time < minHitTime || time > maxHitTime) continue;

      truthMatch.partID = partID;
      truthMatch.trueT0 = time;
      truthMatch.hasTpcTrack = false;

      // Loop over the taggers
      for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
        std::string taggerName = fCrtGeo.GetTagger(i).name;
        // Find the intersections between the true particle and each tagger (if exists)
        geo::Point_t trueCP = fCrtGeo.TaggerCrossingPoint(taggerName, particle);
        // Make map for each tagger with true ID as key and value an XYZ position
        if(trueCP.X() != -99999){ 
          truthMatch.trueCrosses[taggerName] = trueCP;
        }
        bool validCP = fCrtGeo.ValidCrossingPoint(taggerName, particle);
        truthMatch.validCrosses[taggerName] = validCP;
      } 
      
      truthMatching[partID] = truthMatch;

    }

    //--------------------------------------------- CRT HIT MATCHING -------------------------------------------

    // Loop over CRT hits
    std::vector<crt::CRTHit> crtHits;
    int hit_i = 0;
    for (auto const& crtHit: (*crtHitHandle)){
      crtHits.push_back(crtHit);
      // Get tagger of CRT hit
      std::string taggerName = crtHit.tagger;

      hCRTHitTimes->Fill(crtHit.ts1_ns*1e-3);

      // Get associated true particle
      int partID = fCrtBackTrack.TrueIdFromHitId(event, hit_i);
      hit_i++;
      if(truthMatching.find(partID) == truthMatching.end()) continue;

      truthMatching[partID].crtHits[taggerName].push_back(crtHit);
    }

    //------------------------------------------- CRT TRACK MATCHING -------------------------------------------

    // Loop over CRT tracks
    std::vector<crt::CRTTrack> crtTracks;
    int track_i = 0;
    for (auto const& crtTrack : (*crtTrackHandle)){
      crtTracks.push_back(crtTrack);

      // Get associated true particle
      int partID = fCrtBackTrack.TrueIdFromTrackId(event, track_i);
      track_i++;
      if(truthMatching.find(partID) == truthMatching.end()) continue;

      truthMatching[partID].crtTracks.push_back(crtTrack);

    }

    //------------------------------------------- CRT T0 MATCHING ----------------------------------------------

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      // Get the associated true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int partID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if(truthMatching.find(partID) == truthMatching.end()) continue;
      truthMatching[partID].hasTpcTrack = true;

      // Calculate t0 from CRT Hit matching
      //int tpc = fTpcGeo.DetectedInTPC(hits);
      double hitT0 = crtT0Alg.T0FromCRTHits(tpcTrack, crtHits, event);
      if(hitT0 != -99999) truthMatching[partID].hitT0s.push_back(hitT0);

      // Calculate t0 from CRT Track matching
      double trackT0 = crtTrackAlg.T0FromCRTTracks(tpcTrack, crtTracks, event);
      if(trackT0 != -99999) truthMatching[partID].trackT0s.push_back(trackT0);

    }

    //----------------------------------------------------------------------------------------------------------
    //                                        PERFORMANCE METRICS
    //----------------------------------------------------------------------------------------------------------
    // Loop over all the true particles in the event inside the reco window
    for(auto const& truthMatch : truthMatching){
      CRTTruthMatch match = truthMatch.second;
      int partID = match.partID;

      double momentum = particles[partID].P();
      TVector3 start(particles[partID].Vx(), particles[partID].Vy(), particles[partID].Vz());
      TVector3 end(particles[partID].EndX(), particles[partID].EndY(), particles[partID].EndZ());
      double theta = (end-start).Theta();
      double phi = (end-start).Phi();

      //--------------------------------- CRT HIT RECONSTRUCTION PERFORMANCE -------------------------------------
      // Loop over the true crossing points
      for(auto const& trueCross : match.trueCrosses){
        // Find the closest matched reconstructed hit on that tagger
        std::string taggerName = trueCross.first;
        geo::Point_t trueCP = trueCross.second;

        double minDist = 99999;
        crt::CRTHit closestHit;
        // Loop over the associated CRTHits and find the closest
        if(match.crtHits.find(taggerName) == match.crtHits.end()) continue;
        for(auto const& crtHit : match.crtHits[taggerName]){
          double dist = pow(crtHit.x_pos-trueCP.X(), 2) + pow(crtHit.y_pos-trueCP.Y(), 2) + pow(crtHit.z_pos-trueCP.Z(), 2);
          if(dist < minDist){
            closestHit = crtHit;
            minDist = dist;
          }
        }
        // RESOLUTIONS
        // Fill the histogram of resolutions for that tagger
        if(minDist != 99999){
          hTaggerXYZResolution[taggerName]->Fill(closestHit.x_pos-trueCP.X(), closestHit.y_pos-trueCP.Y(), closestHit.z_pos-trueCP.Z());
        }
      }
      
      //------------------------------- CRT TRACK RECONSTRUCTION PERFORMANCE -------------------------------------
     
      // RESOLUTIONS
      double minThetaDiff = 99999;
      double minPhiDiff = 99999;
      double minTheta = 99999;
      double minPhi = 99999;
      bool complete = true;
      // Get true theta and phi in CRT volume
      simb::MCParticle particle = particles[match.partID];
      TVector3 trueStart(particle.Vx(), particle.Vy(), particle.Vz());
      TVector3 trueEnd(particle.EndX(), particle.EndY(), particle.EndZ());
      if(trueEnd.Y() > trueStart.Y()) std::swap(trueStart, trueEnd);
      double trueTheta = (trueEnd - trueStart).Theta();
      double truePhi = (trueEnd - trueStart).Phi();
      // Loop over the associated CRT tracks and find the one with the closest angle
      for(auto const& crtTrack : match.crtTracks){
        // Take largest y point as the start
        TVector3 recoStart(crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
        TVector3 recoEnd(crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
        if(recoEnd.Y() > recoStart.Y()) std::swap(recoStart, recoEnd);
        double recoTheta = (recoEnd - recoStart).Theta();
        double recoPhi = (recoEnd - recoStart).Phi();
        double thetaDiff = TMath::ATan2(TMath::Sin(recoTheta-trueTheta), TMath::Cos(recoTheta-trueTheta));
        double phiDiff = TMath::ATan2(TMath::Sin(recoPhi-truePhi), TMath::Cos(recoPhi-truePhi));
        if(thetaDiff < minThetaDiff){
          minThetaDiff = thetaDiff;
          minPhiDiff = phiDiff;
          minTheta = recoTheta;
          minPhi = recoPhi;
          complete = crtTrack.complete;
        }
      }
      // Fill theta, phi difference plots for complete/incomplete tracks (maybe consider entry/exit points)
      if(minThetaDiff != 99999){
        if(complete){
          hTrackThetaDiff[0]->Fill(minThetaDiff);
          hTrackPhiDiff[0]->Fill(minPhiDiff);
          hTrackTheta[0]->Fill(trueTheta*180./TMath::Pi(), minTheta*180./TMath::Pi());
          hTrackPhi[0]->Fill(truePhi*180./TMath::Pi(), minPhi*180./TMath::Pi());
        }
        else{
          hTrackThetaDiff[1]->Fill(minThetaDiff);
          hTrackPhiDiff[1]->Fill(minPhiDiff);
          hTrackTheta[1]->Fill(trueTheta*180./TMath::Pi(), minTheta*180./TMath::Pi());
          hTrackPhi[1]->Fill(truePhi*180./TMath::Pi(), minPhi*180./TMath::Pi());
        }
        hTrackThetaDiff[2]->Fill(minThetaDiff);
        hTrackPhiDiff[2]->Fill(minPhiDiff);
        hTrackTheta[2]->Fill(trueTheta*180./TMath::Pi(), minTheta*180./TMath::Pi());
        hTrackPhi[2]->Fill(truePhi*180./TMath::Pi(), minPhi*180./TMath::Pi());
      }
      
      //--------------------------------------- CRT HIT T0 PERFORMANCE -------------------------------------------

      // PURITY
      // Loop over the associated T0s and find the one from the longest reco track
      // If it is within a certain time then fill plots with angles, momentum, reco length
      double bestHitT0Diff = 99999;
      for(auto const& hitT0 : match.hitT0s){
        double hitT0Diff = std::abs(match.trueT0 - hitT0);
        if(hitT0Diff < bestHitT0Diff) bestHitT0Diff = hitT0Diff;
      }
      if(bestHitT0Diff != 99999){
        hPurityMom[0][0]->Fill(momentum);
        hPurityTheta[0][0]->Fill(theta*180./TMath::Pi());
        hPurityPhi[0][0]->Fill(phi*180./TMath::Pi());
      }
      if(bestHitT0Diff < 2){
        hPurityMom[0][1]->Fill(momentum);
        hPurityTheta[0][1]->Fill(theta*180./TMath::Pi());
        hPurityPhi[0][1]->Fill(phi*180./TMath::Pi());
      }
     
      //------------------------------------- CRT TRACK T0 PERFORMANCE -------------------------------------------
      // PURITY
      // Loop over the associated T0s and find the one from the longest reco track
      // If it is within a certain time then fill plots with angles, momentum, reco length
      double bestTrackT0Diff = 99999;
      for(auto const& trackT0 : match.trackT0s){
        double trackT0Diff = std::abs(match.trueT0 - trackT0);
        if(trackT0Diff < bestTrackT0Diff) bestTrackT0Diff = trackT0Diff;
      }
      if(bestTrackT0Diff != 99999){
        hPurityMom[1][0]->Fill(momentum);
        hPurityTheta[1][0]->Fill(theta*180./TMath::Pi());
        hPurityPhi[1][0]->Fill(phi*180./TMath::Pi());
      }
      if(bestTrackT0Diff < 2){
        hPurityMom[1][1]->Fill(momentum);
        hPurityTheta[1][1]->Fill(theta*180./TMath::Pi());
        hPurityPhi[1][1]->Fill(phi*180./TMath::Pi());
      }

      // --------------------------------- EFFICIENCIES FOR EVERYTHING -------------------------------------------
      // If there are any matched T0s then fill plots with angles, momentum, reco length for different categories
      bool matchesCRTHit = (match.crtHits.size() > 0);
      bool matchesCRTTrack = (match.crtTracks.size() > 0);
      bool matchesHitT0 = (match.hitT0s.size() > 0);
      bool matchesTrackT0 = (match.trackT0s.size() > 0);
      std::vector<bool> matches{matchesCRTHit, matchesCRTTrack, matchesHitT0, matchesTrackT0};

      std::vector<bool> crtHitCategories;
      std::vector<bool> crtTrackCategories;
      std::vector<bool> hitT0Categories;
      std::vector<bool> trackT0Categories;

      // - Possible match (hits & hitT0: cross >= 1 tagger, tracks & trackT0: cross >= 2 taggers)
      crtHitCategories.push_back((match.trueCrosses.size() > 0));
      crtTrackCategories.push_back((match.trueCrosses.size() > 1));
      hitT0Categories.push_back((match.trueCrosses.size() > 0) && match.hasTpcTrack);
      trackT0Categories.push_back((match.trueCrosses.size() > 1) && match.hasTpcTrack);

      // - Should match: (hits: cross 2 perp strips in tagger, tracks: cross 2 perp strips in 2 taggers, 
      //                  hitT0: hit and tpc track associated, trackT0: crt and tpc track associated)
      bool validHit = false;
      for(auto const& valid : match.validCrosses){
        if(valid.second) validHit = true;
      }
      crtHitCategories.push_back(validHit);
      crtTrackCategories.push_back(match.crtHits.size() > 1);
      hitT0Categories.push_back((match.crtHits.size() > 0 && match.hasTpcTrack));
      trackT0Categories.push_back((match.crtTracks.size() > 0 && match.hasTpcTrack));

      // - Enters CRT volume
      bool entersCRT = fCrtGeo.EntersVolume(particles[partID]);
      crtHitCategories.push_back(entersCRT);
      crtTrackCategories.push_back(entersCRT);
      hitT0Categories.push_back((entersCRT && match.hasTpcTrack));
      trackT0Categories.push_back((entersCRT && match.hasTpcTrack));

      // - Enters TPC volume
      bool entersTPC = fTpcGeo.EntersVolume(particles[partID]) && entersCRT;
      crtHitCategories.push_back(entersTPC);
      crtTrackCategories.push_back(entersTPC);
      hitT0Categories.push_back((entersTPC && match.hasTpcTrack));
      trackT0Categories.push_back((entersTPC && match.hasTpcTrack));

      // - Crosses CRT volume
      bool crossesCRT = fCrtGeo.CrossesVolume(particles[partID]);
      crtHitCategories.push_back(crossesCRT);
      crtTrackCategories.push_back(crossesCRT);
      hitT0Categories.push_back((crossesCRT && match.hasTpcTrack));
      trackT0Categories.push_back((crossesCRT && match.hasTpcTrack));

      // - Crosses TPC volume
      bool crossesTPC = fTpcGeo.CrossesVolume(particles[partID]) && entersCRT;
      crtHitCategories.push_back(crossesTPC);
      crtTrackCategories.push_back(crossesTPC);
      hitT0Categories.push_back((crossesTPC && match.hasTpcTrack));
      trackT0Categories.push_back((crossesTPC && match.hasTpcTrack));

      std::vector<std::vector<bool>> matchCategories{crtHitCategories, crtTrackCategories, hitT0Categories, trackT0Categories};

      for(size_t si = 0; si < stage.size(); si++){
        for(size_t ci = 0; ci < category.size(); ci++){
          if(!matchCategories[si][ci]) continue;
          hEffMom[si][ci][0]->Fill(momentum);
          hEffTheta[si][ci][0]->Fill(theta*180./TMath::Pi());
          hEffPhi[si][ci][0]->Fill(phi*180./TMath::Pi());

          if(!matches[si]) continue;
          hEffMom[si][ci][1]->Fill(momentum);
          hEffTheta[si][ci][1]->Fill(theta*180./TMath::Pi());
          hEffPhi[si][ci][1]->Fill(phi*180./TMath::Pi());
        }
      }
    }

  } // CRTFullRecoAna::analyze()

  void CRTFullRecoAna::endJob(){

  } // CRTFullRecoAna::endJob()
  
  
  DEFINE_ART_MODULE(CRTFullRecoAna)
} // namespace sbnd


