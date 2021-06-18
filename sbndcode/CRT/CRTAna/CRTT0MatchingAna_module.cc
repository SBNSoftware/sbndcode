
////////////////////////////////////////////////////////////////////////
// Class:       CRTT0MatchingAnaAna
// Module Type: analyzer
// File:        CRTT0MatchingAnaAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTHit_Legacy.hh"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

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
#include "cetlib/pow.h" // cet::sum_of_squares()

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

namespace sbnd {

  class CRTT0MatchingAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of tpc track producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTT0MatchAlg::Config> CRTT0Alg {
        Name("CRTT0Alg"),
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTT0MatchingAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Calculate the distance from the track crossing point to CRT overlap coordinates
    double DistToCrtHit(TVector3 trackPos, sbn::crt::CRTHit crtHit);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTHitLabel;   ///< name of CRT producer
    art::InputTag fTPCTrackLabel; ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on

    CRTT0MatchAlg t0Alg;

    CRTGeoAlg fCrtGeo;
    TPCGeoAlg fTpcGeo;

    CRTBackTracker fCrtBackTrack;

    // Histograms
    std::map<std::string, TH1D*> hDCA;
    std::map<std::string, TH1D*> hMatchDCA;
    std::map<std::string, TH1D*> hNoMatchDCA;

    std::map<std::string, TH1D*> hDoL;
    std::map<std::string, TH1D*> hMatchDoL;
    std::map<std::string, TH1D*> hNoMatchDoL;

    std::map<std::string, TH1D*> hT0;
    std::map<std::string, TH1D*> hMatchT0;
    std::map<std::string, TH1D*> hNoMatchT0;

    std::map<std::string, TH1D*> hEffDCATotal;
    std::map<std::string, TH1D*> hEffDCAReco;
    std::map<std::string, TH1D*> hEffDoLTotal;
    std::map<std::string, TH1D*> hEffDoLReco;
    std::map<std::string, TH1D*> hEffLengthTotal;
    std::map<std::string, TH1D*> hEffLengthReco;

    std::map<std::string, TH1D*> hPurityDCATotal;
    std::map<std::string, TH1D*> hPurityDCAReco;
    std::map<std::string, TH1D*> hPurityDoLTotal;
    std::map<std::string, TH1D*> hPurityDoLReco;
    std::map<std::string, TH1D*> hPurityLengthTotal;
    std::map<std::string, TH1D*> hPurityLengthReco;

  }; // class CRTT0MatchingAna


  // Constructor
  CRTT0MatchingAna::CRTT0MatchingAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTHitLabel         (config().CRTHitLabel())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fVerbose             (config().Verbose())
    , t0Alg                (config().CRTT0Alg())
    , fCrtBackTrack        (config().CrtBackTrack())
  {

  } //CRTT0MatchingAna()


  void CRTT0MatchingAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i = 0; i < fCrtGeo.NumTaggers() + 1; i++){
      std::string tagger = "All";
      if(i < fCrtGeo.NumTaggers()){
        tagger = fCrtGeo.GetTagger(i).name;
      }
      hDCA[tagger]        = tfs->make<TH1D>(Form("DCA_%s", tagger.c_str()),        "", 50, 0, 100);
      hMatchDCA[tagger]   = tfs->make<TH1D>(Form("MatchDCA_%s", tagger.c_str()),   "", 50, 0, 100);
      hNoMatchDCA[tagger] = tfs->make<TH1D>(Form("NoMatchDCA_%s", tagger.c_str()), "", 50, 0, 100);

      hDoL[tagger]        = tfs->make<TH1D>(Form("DoL_%s", tagger.c_str()),        "", 100, 0, 0.25);
      hMatchDoL[tagger]   = tfs->make<TH1D>(Form("MatchDoL_%s", tagger.c_str()),   "", 100, 0, 0.25);
      hNoMatchDoL[tagger] = tfs->make<TH1D>(Form("NoMatchDoL_%s", tagger.c_str()), "", 100, 0, 0.25);
      
      hT0[tagger]        = tfs->make<TH1D>(Form("T0_%s", tagger.c_str()),        "", 600, -3000, 3000);
      hMatchT0[tagger]   = tfs->make<TH1D>(Form("MatchT0_%s", tagger.c_str()),   "", 600, -3000, 3000);
      hNoMatchT0[tagger] = tfs->make<TH1D>(Form("NoMatchT0_%s", tagger.c_str()), "", 600, -3000, 3000);
      
      hEffDCATotal[tagger] = tfs->make<TH1D>(Form("EffDCATotal_%s", tagger.c_str()), "", 50, 0, 100);
      hEffDCAReco[tagger]  = tfs->make<TH1D>(Form("EffDCAReco_%s", tagger.c_str()),  "", 50, 0, 100);
      hEffDoLTotal[tagger] = tfs->make<TH1D>(Form("EffDoLTotal_%s", tagger.c_str()), "", 100, 0, 0.25);
      hEffDoLReco[tagger]  = tfs->make<TH1D>(Form("EffDoLReco_%s", tagger.c_str()),  "", 100, 0, 0.25);
      hEffLengthTotal[tagger] = tfs->make<TH1D>(Form("EffLengthTotal_%s", tagger.c_str()), "", 20, 0, 600);
      hEffLengthReco[tagger]  = tfs->make<TH1D>(Form("EffLengthReco_%s", tagger.c_str()),  "", 20, 0, 600);
    
      hPurityDCATotal[tagger] = tfs->make<TH1D>(Form("PurityDCATotal_%s", tagger.c_str()), "", 50, 0, 100);
      hPurityDCAReco[tagger]  = tfs->make<TH1D>(Form("PurityDCAReco_%s", tagger.c_str()),  "", 50, 0, 100);
      hPurityDoLTotal[tagger] = tfs->make<TH1D>(Form("PurityDoLTotal_%s", tagger.c_str()), "", 100, 0, 0.25);
      hPurityDoLReco[tagger]  = tfs->make<TH1D>(Form("PurityDoLReco_%s", tagger.c_str()),  "", 100, 0, 0.25);
      hPurityLengthTotal[tagger] = tfs->make<TH1D>(Form("PurityLengthTotal_%s", tagger.c_str()), "", 20, 0, 600);
      hPurityLengthReco[tagger]  = tfs->make<TH1D>(Form("PurityLengthReco_%s", tagger.c_str()),  "", 20, 0, 600);
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT T0 Matching Ana Module -------------------"<<std::endl;

  } // CRTT0MatchingAna::beginJob()


  void CRTT0MatchingAna::analyze(const art::Event& event)
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
    art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
    if (event.getByLabel(fCRTHitLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitList, crtHitHandle);

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    fCrtBackTrack.Initialize(event);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    
    std::map<int, simb::MCParticle> particles;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
      
    }

    std::map<int, std::vector<std::string>> crtTaggerMap;
    std::vector<sbn::crt::CRTHit> crtHits;
    int hit_i = 0;
    double minHitTime = 99999;
    double maxHitTime = -99999;

    for(auto const& hit : (*crtHitHandle)){
      double hitTime = (double)(int)hit.ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;

      crtHits.push_back(hit);
      int trueId = fCrtBackTrack.TrueIdFromHitId(event, hit_i);
      hit_i++;
      if(trueId == -99999) continue;
      if(crtTaggerMap.find(trueId) == crtTaggerMap.end()){
        crtTaggerMap[trueId].push_back(hit.tagger);
      }
      else if(std::find(crtTaggerMap[trueId].begin(), crtTaggerMap[trueId].end(), hit.tagger) == crtTaggerMap[trueId].end()){
        crtTaggerMap[trueId].push_back(hit.tagger);
      }
    }

    //    std::cout << " New event" << std::endl;
    //----------------------------------------------------------------------------------------------------------
    //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
      if(particles.find(trackTrueID) == particles.end()) continue;
      // Only consider primary muons
      if(!(std::abs(particles[trackTrueID].PdgCode()) == 13 && particles[trackTrueID].Mother() == 0)) continue;

      // Only consider particles inside the reco hit time window
      double trueTime = particles[trackTrueID].T() * 1e-3;
      if(trueTime < minHitTime || trueTime > maxHitTime) continue;

      //      std::cout << "new track " << trueTime << std::endl;
      // Calculate t0 from CRT Hit matching
      matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
      //      std::cout << "closest match " << closest.t0 << std::endl;
      double sin_angle = -99999;
      if(closest.dca != -99999){ 
        hDCA[closest.thishit.tagger]->Fill(closest.dca);
        hDCA["All"]->Fill(closest.dca);
	sin_angle = closest.dca/closest.extrapLen;
        hDoL[closest.thishit.tagger]->Fill(sin_angle);
        hDoL["All"]->Fill(sin_angle);	
	hT0[closest.thishit.tagger]->Fill(closest.t0);
	hT0["All"]->Fill(closest.t0);


      // Is hit matched to that track
        int hitTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, closest.thishit);
        if(hitTrueID == trackTrueID && hitTrueID != -99999){
          hMatchDCA[closest.thishit.tagger]->Fill(closest.dca);
          hMatchDCA["All"]->Fill(closest.dca);
	  hMatchDoL[closest.thishit.tagger]->Fill(sin_angle);
	  hMatchDoL["All"]->Fill(sin_angle);	
	  hMatchT0[closest.thishit.tagger]->Fill(closest.t0);
	  hMatchT0["All"]->Fill(closest.t0);
        }
        else{
          hNoMatchDCA[closest.thishit.tagger]->Fill(closest.dca);
          hNoMatchDCA["All"]->Fill(closest.dca);
	  hNoMatchDoL[closest.thishit.tagger]->Fill(sin_angle);
	  hNoMatchDoL["All"]->Fill(sin_angle);	
	  hNoMatchT0[closest.thishit.tagger]->Fill(closest.t0);
	  hNoMatchT0["All"]->Fill(closest.t0);
        }
      }

      int nbins = hEffDCATotal.begin()->second->GetNbinsX();
      for(int i = 0; i < nbins; i++){
        double DCAcut = hEffDCATotal.begin()->second->GetBinCenter(i);

        // Fill total efficiency histogram with each cut if track matches any hits
        if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
          for(auto const& tagger : crtTaggerMap[trackTrueID]){
            hEffDCATotal[tagger]->Fill(DCAcut);

            // If closest hit is below limit and track matches any hits then fill efficiency
            if(closest.dca < DCAcut && closest.dca != -99999){
              hEffDCAReco[tagger]->Fill(DCAcut);
            }
          }
          // Fill total efficiency histograms
          hEffDCATotal["All"]->Fill(DCAcut);
          if(closest.dca < DCAcut && closest.dca != -99999){
            hEffDCAReco["All"]->Fill(DCAcut);
          }
        }

        // Fill total purity histogram with each cut if closest hit is below limit
        if(closest.dca < DCAcut && closest.dca != -99999){
          hPurityDCATotal[closest.thishit.tagger]->Fill(DCAcut);
          hPurityDCATotal["All"]->Fill(DCAcut);

          // If closest hit is below limit and matched time is correct then fill purity
          double hitTime = closest.thishit.ts1_ns * 1e-3;
          if(particles.find(trackTrueID) != particles.end()){
            if(std::abs(hitTime - trueTime) < 2.){
              hPurityDCAReco[closest.thishit.tagger]->Fill(DCAcut);
              hPurityDCAReco["All"]->Fill(DCAcut);
            }
          }

        }
      }

      nbins = hEffDoLTotal.begin()->second->GetNbinsX();
      
      for(int i = 0; i < nbins; i++){
        double DCAcut = hEffDoLTotal.begin()->second->GetBinCenter(i);

        // Fill total efficiency histogram with each cut if track matches any hits
        if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
          for(auto const& tagger : crtTaggerMap[trackTrueID]){
            hEffDoLTotal[tagger]->Fill(DCAcut);

            // If closest hit is below limit and track matches any hits then fill efficiency
            if(sin_angle < DCAcut && closest.dca != -99999){
              hEffDoLReco[tagger]->Fill(DCAcut);
            }
          }
          // Fill total efficiency histograms
          hEffDoLTotal["All"]->Fill(DCAcut);
          if(sin_angle < DCAcut && closest.dca != -99999){
            hEffDoLReco["All"]->Fill(DCAcut);
          }
        }

        // Fill total purity histogram with each cut if closest hit is below limit
        if(sin_angle < DCAcut && closest.dca != -99999){
          hPurityDoLTotal[closest.thishit.tagger]->Fill(DCAcut);
          hPurityDoLTotal["All"]->Fill(DCAcut);

          // If closest hit is below limit and matched time is correct then fill purity
          double hitTime = closest.thishit.ts1_ns * 1e-3;
          if(particles.find(trackTrueID) != particles.end()){
            if(std::abs(hitTime - trueTime) < 2.){
              hPurityDoLReco[closest.thishit.tagger]->Fill(DCAcut);
              hPurityDoLReco["All"]->Fill(DCAcut);
            }
          }

        }
      }

      double fixedCut = 30.;

      // Fill total efficiency histogram with each cut if track matches any hits
      if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
        for(auto const& tagger : crtTaggerMap[trackTrueID]){
          hEffLengthTotal[tagger]->Fill(tpcTrack.Length());

          // If closest hit is below limit and track matches any hits then fill efficiency
          if(closest.dca < fixedCut && closest.dca >=0 ){
            hEffLengthReco[tagger]->Fill(tpcTrack.Length());
          }
        }
        // Fill total efficiency histograms
        hEffLengthTotal["All"]->Fill(tpcTrack.Length());
        if(closest.dca < fixedCut && closest.dca >=0){
          hEffLengthReco["All"]->Fill(tpcTrack.Length());
        }
      }

      // Fill total purity histogram with each cut if closest hit is below limit
      if(closest.dca < fixedCut && closest.dca >= 0){
        hPurityLengthTotal[closest.thishit.tagger]->Fill(tpcTrack.Length());
        hPurityLengthTotal["All"]->Fill(tpcTrack.Length());

        // If closest hit is below limit and matched time is correct then fill purity
        double hitTime = closest.thishit.ts1_ns * 1e-3;
        if(particles.find(trackTrueID) != particles.end()){
          double trueTime = particles[trackTrueID].T() * 1e-3;
          if(std::abs(hitTime - trueTime) < 2.){
            hPurityLengthReco[closest.thishit.tagger]->Fill(tpcTrack.Length());
            hPurityLengthReco["All"]->Fill(tpcTrack.Length());
          }
        }

      }

    }

    
  } // CRTT0MatchingAna::analyze()


  void CRTT0MatchingAna::endJob(){
  
    
  } // CRTT0MatchingAna::endJob()

  
  DEFINE_ART_MODULE(CRTT0MatchingAna)
} // namespace sbnd
