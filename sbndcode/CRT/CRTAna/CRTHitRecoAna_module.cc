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

      fhicl::Atom<art::InputTag> CRTModuleLabel {
        Name("CRTModuleLabel"),
        Comment("tag of CRT simulation data product")
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
    art::InputTag fCRTModuleLabel;      ///< name of CRT producer
    art::InputTag fCRTHitLabel;         ///< name of CRT hit producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot

    // histograms
    TH1D* hHitTime;
    TH1D* hNpe;
    TH1D* hRecoSipmDist;
    TH1D* hTrueSipmDist;

    TH1D* hEffWidthTotal;
    TH1D* hEffWidthReco;
    TH1D* hEffLengthTotal;
    TH1D* hEffLengthReco;

    TH2D* hTrueRecoSipmDist;
    TH2D* hNpeAngle;
    TH2D* hNpeSipmDist;
    TH2D* hNpeStripDist;

    TGraphAsymmErrors* gEffWidth;
    TGraphAsymmErrors* gEffLength;

    // CRT helpers
    CRTGeoAlg fCrtGeo;
    CRTEventDisplay evd;
    CRTBackTracker fCrtBackTrack;

    int nValidCrosses = 0;
    int nRecoHits = 0;

  }; // class CRTHitRecoAna

  // Constructor
  CRTHitRecoAna::CRTHitRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTModuleLabel       (config().CRTModuleLabel())
    , fCRTHitLabel          (config().CRTHitLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , evd                   (config().Evd())
    , fCrtBackTrack         (config().CrtBackTrack())
  {

  }

  void CRTHitRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    // Define histograms

    hRecoSipmDist   = tfs->make<TH1D>("RecoSipmDist", ";Distance from SiPM (cm);N Tracks",  60, -2,    13);
    hTrueSipmDist   = tfs->make<TH1D>("TrueSipmDist", ";Distance from SiPM (cm);N Tracks",  60, -2,    13);
    hHitTime    = tfs->make<TH1D>("HitTime",  ";Hit time (s);N Tracks",             50, -2000, 3000);
    hNpe        = tfs->make<TH1D>("Npe",      ";Num PE;N Tracks",                   60, 0,     300);

    hEffWidthTotal   = tfs->make<TH1D>("EffWidthTotal", ";Distance from SiPM (cm)",  20, 0, 11.2);
    hEffWidthReco    = tfs->make<TH1D>("EffWidthReco", ";Distance from SiPM (cm)",   20, 0, 11.2);
    hEffLengthTotal  = tfs->make<TH1D>("EffLengthTotal", ";Distance along strip (cm)", 20, 0,  450);
    hEffLengthReco   = tfs->make<TH1D>("EffLengthReco", ";Distance along strip (cm)",  20, 0,  450);

    hTrueRecoSipmDist = tfs->make<TH2D>("TrueRecoSipmDist", ";True SiPM dist (cm);Reco SiPM dist (cm)",  30, -2, 13,  30, -2, 13);
    hNpeAngle         = tfs->make<TH2D>("NpeAngle",         ";Angle to tagger (rad);Num PE", 30, 0, 1.6, 30, 0,  300);
    hNpeSipmDist      = tfs->make<TH2D>("NpeSipmDist",      ";SiPM dist (cm);Num PE", 30, 0, 11.2, 30, 0,  300);
    hNpeStripDist      = tfs->make<TH2D>("NpeStripDist",    ";Strip dist (cm);Num PE", 30, 0, 450, 30, 0,  300);

    gEffWidth = tfs->makeAndRegister<TGraphAsymmErrors>("EffWidth", ";SiPM dist (cm)");
    gEffLength = tfs->makeAndRegister<TGraphAsymmErrors>("EffLength", ";Strip dist (cm)");
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

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Get all the CRT hits
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
    std::map<int, std::vector<crt::CRTHit>> crtHits;
    double minHitTime = 99999;
    double maxHitTime = -99999;
    for(auto const& hit : (*crtHitHandle)){
      nRecoHits++;
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, hit);
      if(trueId == -99999) continue;
      crtHits[trueId].push_back(hit);
      double hitTime = (double)(int)hit.ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;
    }

    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;

      // Only consider particles within the time limit of the generated CRT hits
      double time = particle.T() * 1e-3;
      if(time < minHitTime || time > maxHitTime) continue;

      if(!(std::abs(particle.PdgCode()) == 13 && particle.Mother()==0)) continue;

      for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
        if(fCrtGeo.ValidCrossingPoint(fCrtGeo.GetTagger(i).name, particle)) nValidCrosses++;
      }

      //if(!fCrtGeo.EntersVolume(particle)) continue;

      //FIXME should probably do with auxdet IDEs
      std::vector<std::string> stripNames = fCrtGeo.CrossesStrips(particle);

      for(auto const& stripName : stripNames){
        if(!fCrtGeo.ValidCrossingPoint(fCrtGeo.GetTaggerName(stripName), particle)) continue;
        geo::Point_t cross = fCrtGeo.StripCrossingPoint(stripName, particle);
        double sipmDist = fCrtGeo.DistanceBetweenSipms(cross, stripName);
        double stripDist = fCrtGeo.DistanceDownStrip(cross, stripName);
        hEffWidthTotal->Fill(sipmDist);
        hEffLengthTotal->Fill(stripDist);

        //Loop over crt hits
        if(crtHits.find(partId) == crtHits.end()) continue;
        bool match = false;
        for(size_t i = 0; i < crtHits[partId].size(); i++){
          //If you find a hit on the same tagger then fill
          if(crtHits[partId][i].tagger == fCrtGeo.GetTaggerName(stripName)) match = true;
        }

        if(!match) continue;
        hEffWidthReco->Fill(sipmDist);
        hEffLengthReco->Fill(stripDist);
      }
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    if(crtHitHandle.isValid()){
      art::FindManyP<crt::CRTData> findManyData(crtHitHandle, event, fCRTHitLabel);
      int hit_i = 0;
      for(auto const& hit : (*crtHitHandle)){
        double time = (double)(int)hit.ts1_ns * 1e-3;
        hHitTime->Fill(time);
        //hNpe->Fill(hit.peshit);
        geo::Point_t hitPos {hit.x_pos, hit.y_pos, hit.z_pos};

        // Get the data associated with the hit
        std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(hit_i);
        hit_i++;
        std::vector<std::string> stripNames;
        for(auto const& dat : data){
          std::string name = fCrtGeo.ChannelToStripName(dat->Channel());
          if(std::find(stripNames.begin(), stripNames.end(), name) == stripNames.end()) 
            stripNames.push_back(name);
        }

        // Match hit to true track
        int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, hit);
        if(particles.find(trueId) == particles.end()) continue;

        // Calculate the angle to the tagger
        double angle = TMath::Pi()/2. - fCrtGeo.AngleToTagger(hit.tagger, particles[trueId]);
        hNpeAngle->Fill(angle, hit.peshit);

        if(angle>1.3) hNpe->Fill(hit.peshit);

        // Calculate the true distance from the channel
        for(auto const& stripName : stripNames){
          double hitDist = fCrtGeo.DistanceBetweenSipms(hitPos, stripName);
          hRecoSipmDist->Fill(hitDist);
          geo::Point_t truePos = fCrtGeo.StripCrossingPoint(stripName, particles[trueId]);
          double trueDist = fCrtGeo.DistanceBetweenSipms(truePos, stripName);
          hTrueSipmDist->Fill(trueDist);
          hTrueRecoSipmDist->Fill(trueDist, hitDist);

          hNpeSipmDist->Fill(trueDist, hit.peshit);
          double stripDist = fCrtGeo.DistanceDownStrip(truePos, stripName);
          hNpeStripDist->Fill(stripDist, hit.peshit);
        }

      }
    }

    if(fPlot){
      evd.SetDrawCrtData(true);
      evd.SetDrawCrtHits(true);
      if(fVeryVerbose) evd.SetPrint(true);
      if(fPlotTrackID != -99999) evd.SetTrueId(fPlotTrackID);
      evd.Draw(event);
    }

  } // CRTHitRecoAna::analyze()

  void CRTHitRecoAna::endJob(){

    std::cout<<"True crosses = "<<nValidCrosses<<", num hits = "<<nRecoHits<<" ("<<(double)nRecoHits/nValidCrosses<<")\n";

    gEffLength->BayesDivide(hEffLengthReco, hEffLengthTotal);
    gEffLength->Draw("ap");
    gEffWidth->BayesDivide(hEffWidthReco, hEffWidthTotal);
    gEffWidth->Draw("ap");

  } // CRTHitRecoAna::endJob()

    
  DEFINE_ART_MODULE(CRTHitRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


