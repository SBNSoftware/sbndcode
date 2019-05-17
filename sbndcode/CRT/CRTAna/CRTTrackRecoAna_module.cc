////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackRecoAna
// Module Type: analyzer
// File:        CRTTrackRecoAna_module.cc
//
// Analysis module for evaluating CRT track reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
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
#include <cmath>
#include <algorithm>

namespace sbnd {

  class CRTTrackRecoAna : public art::EDAnalyzer {
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
        Comment("tag of CRT hit data product")
      };

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel"),
        Comment("tag of CRT track data product")
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
    explicit CRTTrackRecoAna(Parameters const& config);
 
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
    art::InputTag fCRTTrackLabel;       ///< name of CRT track producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot
    
    // n-tuples
    std::map<std::string, TH1D*> hTrackDist;

    std::map<std::string, TH1D*> hCrossDistance;

    std::map<std::string, TH1D*> hEffMomTotal;
    std::map<std::string, TH1D*> hEffMomReco;
    std::map<std::string, TH1D*> hEffThetaTotal;
    std::map<std::string, TH1D*> hEffThetaReco;
    std::map<std::string, TH1D*> hEffPhiTotal;
    std::map<std::string, TH1D*> hEffPhiReco;

    // Other variables shared between different methods.
    CRTGeoAlg fCrtGeo;

    CRTEventDisplay evd;
    CRTBackTracker fCrtBackTrack;

    // Performance Counters

    // Truth Counters
    int nParts = 0;
    int nSinglePlane = 0;
    int nMatchTrack = 0;
    int nTwoStripPlanes = 0;

  }; // class CRTTrackRecoAna

  // Constructor
  CRTTrackRecoAna::CRTTrackRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTHitLabel          (config().CRTHitLabel())
    , fCRTTrackLabel        (config().CRTTrackLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , evd                   (config().Evd())
    , fCrtBackTrack         (config().CrtBackTrack())
  {

  }

  void CRTTrackRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    // Define histograms
    for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
      std::string tagger = fCrtGeo.GetTagger(i).name;
      hCrossDistance[tagger] = tfs->make<TH1D>(Form("CrossDistance_%s", tagger.c_str()), "", 40, 0, 100);
      hTrackDist[tagger]     = tfs->make<TH1D>(Form("TrackDist_%s", tagger.c_str()),     "", 40, 0, 200);

      hEffMomTotal[tagger]   = tfs->make<TH1D>(Form("EffMomTotal_%s", tagger.c_str()),   "", 20, 0,    10);
      hEffMomReco[tagger]    = tfs->make<TH1D>(Form("EffMomReco_%s", tagger.c_str()),    "", 20, 0,    10);
      hEffThetaTotal[tagger] = tfs->make<TH1D>(Form("EffThetaTotal_%s", tagger.c_str()), "", 20, 0,    3.2);
      hEffThetaReco[tagger]  = tfs->make<TH1D>(Form("EffThetaReco_%s", tagger.c_str()),  "", 20, 0,    3.2);
      hEffPhiTotal[tagger]   = tfs->make<TH1D>(Form("EffPhiTotal_%s", tagger.c_str()),   "", 20, -3.2, 3.2);
      hEffPhiReco[tagger]    = tfs->make<TH1D>(Form("EffPhiReco_%s", tagger.c_str()),    "", 20, -3.2, 3.2);
    }

    // Initial output
    std::cout<<"----------------- CRT Track Reco Ana Module -------------------"<<std::endl;

    // Take a position that is known to be the center of a strip


  }// CRTTrackRecoAna::beginJob()

  void CRTTrackRecoAna::analyze(const art::Event& event)
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

    // Get all the CRT tracks
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCRTTrackLabel);

    // Get hit to data associations
    art::FindManyP<crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    std::map<int, std::vector<crt::CRTTrack>> crtTracks;
    double minTrackTime = 99999;
    double maxTrackTime = -99999;
    for(auto const& track : (*crtTrackHandle)){
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, track);
      if(trueId == -99999) continue;
      crtTracks[trueId].push_back(track);
      double trackTime = (double)(int)track.ts1_ns * 1e-3;
      if(trackTime < minTrackTime) minTrackTime = trackTime;
      if(trackTime > maxTrackTime) maxTrackTime = trackTime;
    }

    std::map<int, std::vector<crt::CRTHit>> crtHits;
    for(auto const& hit : (*crtHitHandle)){
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, hit);
      if(trueId == -99999) continue;
      crtHits[trueId].push_back(hit);
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
      if(time < minTrackTime || time > maxTrackTime) continue;

      if(!(std::abs(particle.PdgCode()) == 13 && particle.Mother()==0)) continue;

      // Only consider particles which generated more than 1 crt hit
      if(crtHits.find(partId) == crtHits.end()) continue;
      int nPlanesHit = 0;
      std::vector<std::string> hitTaggers;
      for(auto const& hit : crtHits[partId]){
        if(std::find(hitTaggers.begin(), hitTaggers.end(), hit.tagger) != hitTaggers.end()) continue;
        hitTaggers.push_back(hit.tagger);
        nPlanesHit++;
      }
      if(nPlanesHit < 2) continue;

      double momentum = particle.P();
      TVector3 start (particle.Vx(), particle.Vy(), particle.Vz());
      TVector3 end (particle.EndX(), particle.EndY(), particle.EndZ());
      double theta = (end-start).Theta();
      double phi = (end-start).Phi();
      for(auto const& tagger : hitTaggers){
        hEffMomTotal[tagger]->Fill(momentum);
        hEffThetaTotal[tagger]->Fill(theta);
        hEffPhiTotal[tagger]->Fill(phi);
      }

      if(crtTracks.find(partId) == crtTracks.end()) continue;
      for(auto const& tagger : hitTaggers){
        hEffMomReco[tagger]->Fill(momentum);
        hEffThetaReco[tagger]->Fill(theta);
        hEffPhiReco[tagger]->Fill(phi);
      }
    }

    //----------------------------------------------------------------------------------------------------------
    //                                        CRT TRACK ANALYSIS
    //----------------------------------------------------------------------------------------------------------
    int track_i = 0;
    for(auto const& track : (*crtTrackHandle)){

      std::vector<art::Ptr<crt::CRTHit>> hits = findManyHits.at(track_i);
      track_i++;

      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, track);
      if(particles.find(trueId) == particles.end()) continue;

      // Calculate the average distance over the length of the track 
      TVector3 start {track.x1_pos, track.y1_pos, track.z1_pos};
      TVector3 end {track.x2_pos, track.y2_pos, track.z2_pos};
      double denominator = (end - start).Mag();

      std::vector<double> crtLims = fCrtGeo.CRTLimits();
      simb::MCParticle particle = particles[trueId];

      int nTraj = particle.NumberTrajectoryPoints();
      double aveDCA = 0;
      double npts = 0;
      for(int j = 0; j < nTraj; j++){
        TVector3 pt (particle.Vx(j), particle.Vy(j), particle.Vz(j));
        if (pt.X() >= crtLims[0] && pt.X() <= crtLims[3] && pt.Y() >= crtLims[1] && pt.Y() <= crtLims[4] && pt.Z() >= crtLims[2] && pt.Z() <= crtLims[5]){
          double numerator = ((pt - start).Cross(pt - end)).Mag();
          aveDCA += numerator/denominator;
          npts++;
        }
      }

      aveDCA = aveDCA/npts;

      // Find the taggers and positions of the true crossing points
      crt::CRTHit startHit, endHit;
      for(auto const& hit : hits){
        geo::Point_t trueCross = fCrtGeo.TaggerCrossingPoint(hit->tagger, particles[trueId]);
        if(trueCross.X() == -99999) continue;
        // For each tagger calculate the distance between the CRT track and true track
        double dist = std::sqrt(std::pow(hit->x_pos - trueCross.X(), 2)
                                + std::pow(hit->y_pos - trueCross.Y(), 2)
                                + std::pow(hit->z_pos - trueCross.Z(), 2));
        hCrossDistance[hit->tagger]->Fill(dist);
        hTrackDist[hit->tagger]->Fill(aveDCA);
      }

    }
    
    if(fPlot){
      evd.SetDrawCrtHits(true);
      evd.SetDrawCrtTracks(true);
      if(fVeryVerbose) evd.SetPrint(true);
      if(fPlotTrackID != -99999) evd.SetTrueId(fPlotTrackID);
      evd.Draw(event);
    }


  } // CRTTrackRecoAna::analyze()

  void CRTTrackRecoAna::endJob(){  

  } // CRTTrackRecoAna::endJob()


  DEFINE_ART_MODULE(CRTTrackRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


