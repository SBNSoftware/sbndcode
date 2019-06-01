////////////////////////////////////////////////////////////////////////
// Class:       CRTHitRecoAna
// Module Type: analyzer
// File:        CRTHitRecoAna_module.cc
//
// Analysis module for evaluating CRT hit reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTEventDisplay.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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
#include "TGraphAsymmErrors.h"
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace {
  // Local namespace for local functions
  // Declare here, define later

}

namespace sbnd {

  class CRTHitRecoAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel"),
        Comment("tag of CRT hit product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<bool> VeryVerbose {
        Name("VeryVerbose"),
        Comment("Print detailed information about every track")
      };

      fhicl::Atom<bool> Plot {
        Name("Plot"),
        Comment("Plot tracks")
      };

      fhicl::Atom<int> PlotTrackID {
        Name("PlotTrackID"),
        Comment("ID of track to plot")
      };

      fhicl::Atom<double> MinAngleNpePlot {
        Name("MinAngleNpePlot"),
        Comment("Minimum angle for plotting Npe per hit")
      };

      fhicl::Table<CRTEventDisplay::Config> Evd {
        Name("Evd"),
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTHitRecoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTHitLabel;         ///< name of CRT hit producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot
    double        fMinAngleNpePlot;     ///< Maximum angle for plotting Npe per hit

    // histograms
    TH1D* hNpeAngleCut;

    std::map<std::string, TH1D*> hHitTime;
    std::map<std::string, TH1D*> hNpe;
    std::map<std::string, TH1D*> hAngle;
    std::map<std::string, TH1D*> hRecoSipmDist;
    std::map<std::string, TH1D*> hTrueSipmDist;

    std::map<std::string, TH1D*> hEffWidthTotal;
    std::map<std::string, TH1D*> hEffWidthReco;
    std::map<std::string, TH1D*> hEffLengthTotal;
    std::map<std::string, TH1D*> hEffLengthReco;

    std::map<std::string,TH2D*> hTrueRecoSipmDist;
    std::map<std::string,TH2D*> hNpeAngle;
    std::map<std::string,TH2D*> hNpeSipmDist;
    std::map<std::string,TH2D*> hNpeStripDist;

    // CRT helpers
    CRTGeoAlg fCrtGeo;
    CRTEventDisplay evd;
    CRTBackTracker fCrtBackTrack;

  }; // class CRTHitRecoAna

  // Constructor
  CRTHitRecoAna::CRTHitRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTHitLabel          (config().CRTHitLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , fMinAngleNpePlot      (config().MinAngleNpePlot())
    , evd                   (config().Evd())
    , fCrtBackTrack         (config().CrtBackTrack())
  {

  }

  void CRTHitRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    // Define histograms
    hNpeAngleCut         = tfs->make<TH1D>("NpeAngleCut", "", 40, 0,  600);
    for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
      std::string tagger = fCrtGeo.GetTagger(i).name;
      hRecoSipmDist[tagger] = tfs->make<TH1D>(Form("RecoSipmDist_%s", tagger.c_str()), "", 40, -2,    13);
      hTrueSipmDist[tagger] = tfs->make<TH1D>(Form("TrueSipmDist_%s", tagger.c_str()), "", 40, -2,    13);
      hHitTime[tagger]      = tfs->make<TH1D>(Form("HitTime_%s", tagger.c_str()),      "", 40, -2000, 3000);
      hNpe[tagger]          = tfs->make<TH1D>(Form("Npe_%s", tagger.c_str()),          "", 40, 0,     600);
      hAngle[tagger]        = tfs->make<TH1D>(Form("Angle_%s", tagger.c_str()),        "", 40, 0,     1.6);

      hEffWidthTotal[tagger]   = tfs->make<TH1D>(Form("EffWidthTotal_%s", tagger.c_str()),  "",  20, 0, 11.2);
      hEffWidthReco[tagger]    = tfs->make<TH1D>(Form("EffWidthReco_%s", tagger.c_str()),   "",  20, 0, 11.2);
      hEffLengthTotal[tagger]  = tfs->make<TH1D>(Form("EffLengthTotal_%s", tagger.c_str()), "",  20, 0,  450);
      hEffLengthReco[tagger]   = tfs->make<TH1D>(Form("EffLengthReco_%s", tagger.c_str()),  "",  20, 0,  450);

      hTrueRecoSipmDist[tagger] = tfs->make<TH2D>(Form("TrueRecoSipmDist_%s", tagger.c_str()), "", 30, -2, 13,   30, -2, 13);
      hNpeAngle[tagger]         = tfs->make<TH2D>(Form("NpeAngle_%s", tagger.c_str()),         "", 30, 0,  1.6,  30, 0,  600);
      hNpeSipmDist[tagger]      = tfs->make<TH2D>(Form("NpeSipmDist_%s", tagger.c_str()),      "", 30, 0,  11.2, 30, 0,  600);
      hNpeStripDist[tagger]     = tfs->make<TH2D>(Form("NpeStripDist_%s", tagger.c_str()),     "", 30, 0,  450,  30, 0,  600);
    }

    // Initial output
    std::cout<<"----------------- CRT Hit Reco Ana Module -------------------"<<std::endl;

  }// CRTHitRecoAna::beginJob()

  void CRTHitRecoAna::analyze(const art::Event& event)
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
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get all the CRT hits
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);

    // Get hit to data associations
    art::FindManyP<crt::CRTData> findManyData(crtHitHandle, event, fCRTHitLabel);

    fCrtBackTrack.Initialize(event);
    
    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    std::map<int, std::vector<crt::CRTHit>> crtHits;
    double minHitTime = 99999;
    double maxHitTime = -99999;
    int ht_i = 0;
    for(auto const& hit : (*crtHitHandle)){
      int trueId = fCrtBackTrack.TrueIdFromHitId(event, ht_i);
      ht_i++;
      if(trueId == -99999) continue;
      crtHits[trueId].push_back(hit);
      double hitTime = (double)(int)hit.ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          EFFICIENCIES
    //----------------------------------------------------------------------------------------------------------
    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;

      // Only consider particles within the time limit of the generated CRT hits
      double time = particle.T() * 1e-3;
      if(time < minHitTime || time > maxHitTime) continue;

      if(!(std::abs(particle.PdgCode()) == 13 && particle.Mother()==0)) continue;

      std::vector<std::string> stripNames = fCrtGeo.CrossesStrips(particle);

      for(auto const& stripName : stripNames){
        std::string tagger = fCrtGeo.GetTaggerName(stripName);
        if(!fCrtGeo.ValidCrossingPoint(tagger, particle)) continue;
        geo::Point_t cross = fCrtGeo.StripCrossingPoint(stripName, particle);
        double sipmDist = fCrtGeo.DistanceBetweenSipms(cross, stripName);
        double stripDist = fCrtGeo.DistanceDownStrip(cross, stripName);
        hEffWidthTotal[tagger]->Fill(sipmDist);
        hEffLengthTotal[tagger]->Fill(stripDist);

        //Loop over crt hits
        if(crtHits.find(partId) == crtHits.end()) continue;
        bool match = false;
        for(size_t i = 0; i < crtHits[partId].size(); i++){
          //If you find a hit on the same tagger then fill
          if(crtHits[partId][i].tagger == tagger) match = true;
        }

        if(!match) continue;
        hEffWidthReco[tagger]->Fill(sipmDist);
        hEffLengthReco[tagger]->Fill(stripDist);
      }
    }

    //----------------------------------------------------------------------------------------------------------
    //                                        CRT HIT ANALYSIS
    //----------------------------------------------------------------------------------------------------------
    int hit_i = 0;
    for(auto const& hit : (*crtHitHandle)){
      std::string tagger = hit.tagger;

      // Calculate the time of the CRT hit in us
      double time = (double)(int)hit.ts1_ns * 1e-3;
      hHitTime[tagger]->Fill(time);

      geo::Point_t hitPos {hit.x_pos, hit.y_pos, hit.z_pos};

      // Get the data associated with the hit
      std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(hit_i);

      // Get all the strips from the crt hit
      std::vector<std::string> stripNames;
      for(auto const& dat : data){
        std::string name = fCrtGeo.ChannelToStripName(dat->Channel());
        if(std::find(stripNames.begin(), stripNames.end(), name) == stripNames.end()) 
          stripNames.push_back(name);
      }

      // Match hit to true track
      int trueId = fCrtBackTrack.TrueIdFromHitId(event, hit_i);
      hit_i++;
      if(particles.find(trueId) == particles.end()) continue;
      if(!(std::abs(particles[trueId].PdgCode()) == 13 && particles[trueId].Mother()==0)) continue;

      // Calculate the angle to the tagger
      double angle = TMath::Pi()/2. - fCrtGeo.AngleToTagger(tagger, particles[trueId]);
      if(angle > fMinAngleNpePlot) hNpeAngleCut->Fill(hit.peshit);
      hNpe[tagger]->Fill(hit.peshit);
      hAngle[tagger]->Fill(angle);
      hNpeAngle[tagger]->Fill(angle, hit.peshit);

      // Calculate the true distance from the channel
      for(auto const& stripName : stripNames){
        // Calculate the reconstructed position between sipms
        double hitDist = fCrtGeo.DistanceBetweenSipms(hitPos, stripName);
        hRecoSipmDist[tagger]->Fill(hitDist);

        // Calculate the true position between the sipms
        geo::Point_t truePos = fCrtGeo.StripCrossingPoint(stripName, particles[trueId]);
        double trueDist = fCrtGeo.DistanceBetweenSipms(truePos, stripName);
        hTrueSipmDist[tagger]->Fill(trueDist);
        hTrueRecoSipmDist[tagger]->Fill(trueDist, hitDist);

        // Fill the number of pe as a function of distances
        hNpeSipmDist[tagger]->Fill(trueDist, hit.peshit);
        double stripDist = fCrtGeo.DistanceDownStrip(truePos, stripName);
        hNpeStripDist[tagger]->Fill(stripDist, hit.peshit);
      }

    }

    if(fPlot){
      evd.SetDrawCrtData(true);
      evd.SetDrawCrtHits(true);
      evd.SetDrawCrtTracks(true);
      if(fVeryVerbose) evd.SetPrint(true);
      if(fPlotTrackID != -99999) evd.SetTrueId(fPlotTrackID);
      evd.Draw(event);
    }

  } // CRTHitRecoAna::analyze()

  void CRTHitRecoAna::endJob(){

  } // CRTHitRecoAna::endJob()

  DEFINE_ART_MODULE(CRTHitRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


