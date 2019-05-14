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
    TH1D* hSipmDist;

    TH2D* hTrueRecoSipmDist;
    TH2D* hNpeAngle;

    // CRT helpers
    CRTGeoAlg fCrtGeo;
    CRTEventDisplay evd;
    CRTBackTracker fCrtBackTrack;

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

    hSipmDist   = tfs->make<TH1D>("SipmDist", ";Distance from SiPM (cm);N Tracks",  50, 0,    13);
    hHitTime    = tfs->make<TH1D>("HitTime",  ";Hit time (s);N Tracks",             50, -5000, 5000);
    hNpe        = tfs->make<TH1D>("Npe",      ";Num PE;N Tracks",                   50, 0,     300);

    hTrueRecoSipmDist = tfs->make<TH2D>("TrueRecoSipmDist", ";True SiPM dist (cm);Reco SiPM dist (cm)",  20, -10,   13,  20, 0, 13);
    hNpeAngle         = tfs->make<TH2D>("NpeAngle",         ";Angle to tagger (rad);Num PE",             40, 0, 3.2, 10, 0,  300);
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
    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
    if(crtHitHandle.isValid()){
      art::FindManyP<crt::CRTData> findManyData(crtHitHandle, event, fCRTHitLabel);
      int hit_i = 0;
      for(auto const& hit : (*crtHitHandle)){
        double time = (double)(int)hit.ts1_ns * 1e-3;
        hHitTime->Fill(time);
        hNpe->Fill(hit.peshit);
        geo::Point_t hitPos {hit.x_pos, hit.y_pos, hit.z_pos};

        // Get the data associated with the hit
        std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(hit_i);
        hit_i++;
        std::vector<int> chs;
        for(auto const& dat : data){
          if(!(dat->Channel() % 2)) chs.push_back(dat->Channel());
        }
        if(chs.size() == 0) continue;
        int ch1 = chs[0];
        int ch2 = chs[0];
        if(chs.size() == 2) ch2 = chs[1];

        // Work out how far along the strips the hit position is
        double hitDist1 = std::abs(fCrtGeo.DistanceBetweenSipms(hitPos, ch1));
        hSipmDist->Fill(hitDist1);
        double hitDist2 = std::abs(fCrtGeo.DistanceBetweenSipms(hitPos, ch2));
        if(chs.size() == 2){
          hSipmDist->Fill(hitDist2);
        }

        // Match hit to true track
        int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, hit);
        if(particles.find(trueId) == particles.end()) continue;

        // Calculate the angle to the tagger
        double angle = fCrtGeo.AngleToTagger(hit.tagger, particles[trueId]);
        hNpeAngle->Fill(angle, hit.peshit);

        // Calculate the true distance from the channel
        geo::Point_t truePos = fCrtGeo.TaggerCrossingPoint(hit.tagger, particles[trueId]);
        double trueDist1 = std::abs(fCrtGeo.DistanceBetweenSipms(truePos, ch1));
        hTrueRecoSipmDist->Fill(trueDist1, hitDist1);
        if(chs.size() == 2){
          double trueDist2 = std::abs(fCrtGeo.DistanceBetweenSipms(truePos, ch2));
          hTrueRecoSipmDist->Fill(trueDist2, hitDist2);
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

  } // CRTHitRecoAna::endJob()

    
  DEFINE_ART_MODULE(CRTHitRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


