////////////////////////////////////////////////////////////////////////
// Class:       RunningT0TaggingAna
// Module Type: analyzer
// File:        RunningT0TaggingAna_module.cc
//
// Example code for running t0 tagging within an analysis script
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// C++ includes
#include <map>
#include <vector>

namespace sbnd {

  class RunningT0Tagging : public art::EDAnalyzer {
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
        Comment("tag of CRT hit producer data product")
      };

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTT0MatchAlg::Config> CRTHitMatch {
        Name("CRTHitMatch"),
      };

      fhicl::Table<CRTTrackMatchAlg::Config> CRTTrackMatch {
        Name("CRTTrackMatch"),
      };

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit RunningT0Tagging(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;  ///< name of detsim producer
    art::InputTag fCRTHitLabel;     ///< name of CRT hit producer
    art::InputTag fCRTTrackLabel;   ///< name of CRT track producer
    art::InputTag fTPCTrackLabel;   ///< name of TPC track producer
    bool          fVerbose;         ///< print information about what's going on

    CRTT0MatchAlg    fCRTHitMatch;   ///< CRT hit - TPC track matching algorithms
    CRTTrackMatchAlg fCRTTrackMatch; ///< CRT track - TPC track matching algorithms

    // Performance counters
    int nTotal = 0;
    int nCrtHitMatched = 0;
    int nCrtHitCorrect = 0;
    int nCrtTrackMatched = 0;
    int nCrtTrackCorrect = 0;

  }; // class RunningT0Tagging


  // Constructor
  RunningT0Tagging::RunningT0Tagging(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTHitLabel         (config().CRTHitLabel())
    , fCRTTrackLabel       (config().CRTTrackLabel())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fVerbose             (config().Verbose())
    , fCRTHitMatch         (config().CRTHitMatch())
    , fCRTTrackMatch       (config().CRTTrackMatch())
  {

  } //RunningT0Tagging()


  void RunningT0Tagging::beginJob()
  {

    // Initial output
    if(fVerbose) std::cout<<"----------------- Running CRT t0 tagging example -------------------"<<std::endl;

  } // RunningT0Tagging::beginJob()


  void RunningT0Tagging::analyze(const art::Event& event)
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
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Put into a map for easier access
    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partID = particle.TrackId();
      particles[partID] = particle;
    }

    // Get CRT hits from the event
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
    std::vector<crt::CRTHit> crtHits;
    for(auto const& hit : (*crtHitHandle)){
      crtHits.push_back(hit);
    }

    // Get CRT tracks from the event
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCRTTrackLabel);
    std::vector<crt::CRTTrack> crtTracks;
    for(auto const& track : (*crtTrackHandle)){
      crtTracks.push_back(track);
    }

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);

    // Get the associated hits
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                        RUNNING THE T0 TAGGING
    //----------------------------------------------------------------------------------------------------------

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

      // Get the true particle
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      if(particles.find(trackTrueID) == particles.end()) continue;

      // Only consider primary muons
      if(!(std::abs(particles[trackTrueID].PdgCode()) == 13 && particles[trackTrueID].Mother() == 0)) continue;

      // Only consider particles inside the reco hit time window
      double trueTime = particles[trackTrueID].T() * 1e-3; // [us]
      if(fVerbose) std::cout<<"True track time = "<<trueTime<<" us\n";

      nTotal++;

      // Calculate t0 from CRT hit matching
      double hitT0 = fCRTHitMatch.T0FromCRTHits(tpcTrack, crtHits, event);
      if(hitT0 == -99999){
        if(fVerbose) std::cout<<"Couldn't match to CRT hit.\n";
      }
      else{
        nCrtHitMatched++;
        if(fVerbose) std::cout<<"Hit matched time = "<<hitT0<<" us\n";
        if(std::abs(trueTime - hitT0) < 2) nCrtHitCorrect++;
      }
      // Calculate t0 from CRT track matching
      double trackT0 = fCRTTrackMatch.T0FromCRTTracks(tpcTrack, crtTracks, event);
      if(trackT0 == -99999){
        if(fVerbose) std::cout<<"Couldn't match to CRT track.\n";
      }
      else{
        nCrtTrackMatched++;
        if(fVerbose) std::cout<<"Track matched time = "<<trackT0<<" us\n";
        if(std::abs(trueTime - trackT0) < 2) nCrtTrackCorrect++;
      }

      if(fVerbose) std::cout<<"\n";

    }

    
  } // RunningT0Tagging::analyze()


  void RunningT0Tagging::endJob(){

    std::cout<<"Running CRT t0 tagging:\n"
             <<"-- CRT hit matching:\n"
             <<"---- Efficiency = "<<((double)nCrtHitMatched/nTotal)*100<<"%\n"
             <<"---- Purity     = "<<((double)nCrtHitCorrect/nCrtHitMatched)*100<<"%\n"
             <<"-- CRT track matching:\n"
             <<"---- Efficiency = "<<((double)nCrtTrackMatched/nTotal)*100<<"%\n"
             <<"---- Purity     = "<<((double)nCrtTrackCorrect/nCrtTrackMatched)*100<<"%\n";
    
  } // RunningT0Tagging::endJob()

  
  DEFINE_ART_MODULE(RunningT0Tagging)
} // namespace sbnd


