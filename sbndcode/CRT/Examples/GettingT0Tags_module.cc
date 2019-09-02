////////////////////////////////////////////////////////////////////////
// Class:       GettingT0Tags
// Module Type: analyzer
// File:        GettingT0Tags_module.cc
//
// Example code for getting CRT t0 tags from reco file
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h" //FIXME Why is this needed to compile?!

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
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

  class GettingT0Tags : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CRTHitT0Label {
        Name("CRTHitT0Label"),
        Comment("tag of CRT hit matching data product")
      };

      fhicl::Atom<art::InputTag> CRTTrackT0Label {
        Name("CRTTrackT0Label"),
        Comment("tag of CRT track matching data product")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit GettingT0Tags(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;  ///< name of detsim producer
    art::InputTag fCRTHitT0Label;   ///< name of CRT hit matching module
    art::InputTag fCRTTrackT0Label; ///< name of CRT track matching module
    art::InputTag fTPCTrackLabel;   ///< name of TPC track producer
    bool          fVerbose;         ///< print information about what's going on

    // Performance counters
    int nTotal = 0;
    int nCrtHitMatched = 0;
    int nCrtHitCorrect = 0;
    int nCrtTrackMatched = 0;
    int nCrtTrackCorrect = 0;

  }; // class GettingT0Tags


  // Constructor
  GettingT0Tags::GettingT0Tags(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTHitT0Label       (config().CRTHitT0Label())
    , fCRTTrackT0Label     (config().CRTTrackT0Label())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fVerbose             (config().Verbose())
  {

  } //GettingT0Tags()


  void GettingT0Tags::beginJob()
  {

    // Initial output
    if(fVerbose) std::cout<<"----------------- Getting CRT t0 tags example -------------------"<<std::endl;

  } // GettingT0Tags::beginJob()


  void GettingT0Tags::analyze(const art::Event& event)
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

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);

    // Get the associated hits
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    // Get the associated T0 from CRT hit matching
    art::FindManyP<anab::T0> findManyHitT0(tpcTrackHandle, event, fCRTHitT0Label);

    // Get the associated T0 from CRT track matching
    art::FindManyP<anab::T0> findManyTrackT0(tpcTrackHandle, event, fCRTTrackT0Label);

    //----------------------------------------------------------------------------------------------------------
    //                                        GETTING THE T0 TAGS
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

      // Retrieve the t0 from CRT hit matching
      std::vector<art::Ptr<anab::T0>> hitT0s = findManyHitT0.at(tpcTrack.ID());
      if(hitT0s.size() < 1){
        if(fVerbose) std::cout<<"No matching CRT hit.\n";
      }
      else{
        nCrtHitMatched++;
        double hitT0 = hitT0s[0]->Time() * 1e-3; // [us]
        if(fVerbose) std::cout<<"Hit matched time = "<<hitT0<<" us\n";
        if(std::abs(trueTime - hitT0) < 2) nCrtHitCorrect++;
      }

      // Retrieve the t0 from CRT track matching
      std::vector<art::Ptr<anab::T0>> trackT0s = findManyTrackT0.at(tpcTrack.ID());
      if(trackT0s.size() < 1){
        if(fVerbose) std::cout<<"No matching CRT track.\n";
      }
      else{
        nCrtTrackMatched++;
        double trackT0 = trackT0s[0]->Time() * 1e-3; // [us]
        if(fVerbose) std::cout<<"Track matched time = "<<trackT0<<" us\n";
        if(std::abs(trueTime - trackT0) < 2) nCrtTrackCorrect++;
      }

      if(fVerbose) std::cout<<"\n";

    }

    
  } // GettingT0Tags::analyze()


  void GettingT0Tags::endJob(){

    std::cout<<"Getting CRT t0 tags from reco file:\n"
             <<"-- CRT hit matching:\n"
             <<"---- Efficiency = "<<((double)nCrtHitMatched/nTotal)*100<<"%\n"
             <<"---- Purity     = "<<((double)nCrtHitCorrect/nCrtHitMatched)*100<<"%\n"
             <<"-- CRT track matching:\n"
             <<"---- Efficiency = "<<((double)nCrtTrackMatched/nTotal)*100<<"%\n"
             <<"---- Purity     = "<<((double)nCrtTrackCorrect/nCrtTrackMatched)*100<<"%\n";
   
  } // GettingT0Tags::endJob()

  
  DEFINE_ART_MODULE(GettingT0Tags)
} // namespace sbnd


