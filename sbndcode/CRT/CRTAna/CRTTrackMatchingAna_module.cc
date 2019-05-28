////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackMatchingAnaAna
// Module Type: analyzer
// File:        CRTTrackMatchingAnaAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/CRT/CRTUtils/CRTEventDisplay.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
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

  class CRTTrackMatchingAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel"),
        Comment("tag of tpc track producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<CRTTrackMatchAlg::Config> CRTTrackAlg {
        Name("CRTTrackAlg"),
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

      fhicl::Table<CRTEventDisplay::Config> Evd {
        Name("Evd"),
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTTrackMatchingAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Calculate the distance from the track crossing point to CRT overlap coordinates
    double DistToCrtHit(TVector3 trackPos, crt::CRTHit crtHit);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTTrackLabel;   ///< name of CRT producer
    art::InputTag fTPCTrackLabel; ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on

    CRTTrackMatchAlg trackAlg;

    CRTGeoAlg fCrtGeo;
    TPCGeoAlg fTpcGeo;

    CRTBackTracker fCrtBackTrack;
    CRTEventDisplay evd;

    // Histograms
    TH1D* hAngle;
    TH1D* hMatchAngle;
    TH1D* hNoMatchAngle;

    TH1D* hDCA;
    TH1D* hMatchDCA;
    TH1D* hNoMatchDCA;

    TH2D* hMatchAngleDCA;
    TH2D* hNoMatchAngleDCA;

    TH1D* hEffAngleTotal;
    TH1D* hEffAngleReco;
    TH1D* hPurityAngleTotal;
    TH1D* hPurityAngleReco;

    TH1D* hEffDCATotal;
    TH1D* hEffDCAReco;
    TH1D* hPurityDCATotal;
    TH1D* hPurityDCAReco;

  }; // class CRTTrackMatchingAna


  // Constructor
  CRTTrackMatchingAna::CRTTrackMatchingAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTTrackLabel       (config().CRTTrackLabel())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fVerbose             (config().Verbose())
    , trackAlg             (config().CRTTrackAlg())
    , fCrtBackTrack        (config().CrtBackTrack())
    , evd                   (config().Evd())
  {

  } //CRTTrackMatchingAna()


  void CRTTrackMatchingAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    
    hAngle         = tfs->make<TH1D>("Angle",         "", 50, 0, 2);
    hMatchAngle    = tfs->make<TH1D>("MatchAngle",    "", 50, 0, 2);
    hNoMatchAngle  = tfs->make<TH1D>("NoMatchAngle",  "", 50, 0, 2);

    hDCA        = tfs->make<TH1D>("DCA",        "", 50, 0, 150);
    hMatchDCA   = tfs->make<TH1D>("MatchDCA",   "", 50, 0, 150);
    hNoMatchDCA = tfs->make<TH1D>("NoMatchDCA", "", 50, 0, 150);

    hMatchAngleDCA   = tfs->make<TH2D>("MatchAngleDCA",   "", 20, 0, 2, 20, 0, 150);
    hNoMatchAngleDCA = tfs->make<TH2D>("NoMatchAngleDCA", "", 20, 0, 2, 20, 0, 150);

    hEffAngleTotal    = tfs->make<TH1D>("EffAngleTotal",    "", 20, 0, 1);
    hEffAngleReco     = tfs->make<TH1D>("EffAngleReco",     "", 20, 0, 1);
    hPurityAngleTotal = tfs->make<TH1D>("PurityAngleTotal", "", 20, 0, 1);
    hPurityAngleReco  = tfs->make<TH1D>("PurityAngleReco",  "", 20, 0, 1);

    hEffDCATotal    = tfs->make<TH1D>("EffDCATotal",    "", 20, 0, 100);
    hEffDCAReco     = tfs->make<TH1D>("EffDCAReco",     "", 20, 0, 100);
    hPurityDCATotal = tfs->make<TH1D>("PurityDCATotal", "", 20, 0, 100);
    hPurityDCAReco  = tfs->make<TH1D>("PurityDCAReco",  "", 20, 0, 100);

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT T0 Matching Ana Module -------------------"<<std::endl;

  } // CRTTrackMatchingAna::beginJob()


  void CRTTrackMatchingAna::analyze(const art::Event& event)
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

    // Get CRT hits from the event
    art::Handle< std::vector<crt::CRTTrack>> crtTrackHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackHandle);

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

    std::vector<crt::CRTTrack> crtTracks;
    std::map<int, std::vector<crt::CRTTrack>> crtTrackMap;
    int track_i = 0;
    double minTrackTime = 99999;
    double maxTrackTime = -99999;

    for(auto const& track : (*crtTrackHandle)){
      crtTracks.push_back(track);
      int trueID = fCrtBackTrack.TrueIdFromTrackId(event, track_i);
      track_i++;
      if(trueID != -99999) crtTrackMap[trueID].push_back(track);

      double trackTime = (double)(int)track.ts1_ns * 1e-3;
      if(trackTime < minTrackTime) minTrackTime = trackTime;
      if(trackTime > maxTrackTime) maxTrackTime = trackTime;
    }

    //----------------------------------------------------------------------------------------------------------
    //                            ANGLE AND DISTANCE OF CLOSEST APPROACH ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);

      if(particles.find(trackTrueID) == particles.end()) continue;
      // Only consider primary muons
      if(!(std::abs(particles[trackTrueID].PdgCode()) == 13 && particles[trackTrueID].Mother() == 0)) continue;

      // Only consider particles withing time window of reco track
      double trueTime = particles[trackTrueID].T() * 1e-3;
      if(trueTime < minTrackTime || trueTime > maxTrackTime) continue;

      //----------------------------------------------------------------------------------------------------------
      //                                        SINGLE ANGLE CUT ANALYSIS
      //----------------------------------------------------------------------------------------------------------
      // Find the closest track by angle
      std::pair<crt::CRTTrack, double> closestAngle = trackAlg.ClosestCRTTrackByAngle(tpcTrack, crtTracks, event);
      if(closestAngle.second != -99999){ 
        hAngle->Fill(closestAngle.second);
      }
      // Is crt track matched to that tpc track
      if(closestAngle.second != -99999){
        int crtTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, closestAngle.first);
        if(crtTrueID == trackTrueID && crtTrueID != -99999){
          hMatchAngle->Fill(closestAngle.second);
        }
        else{
          hNoMatchAngle->Fill(closestAngle.second);
        }
      }

      int nbinsAngle = hEffAngleTotal->GetNbinsX();
      for(int i = 0; i < nbinsAngle; i++){
        double angleCut = hEffAngleTotal->GetBinCenter(i);

        // Fill total efficiency histogram with each cut if track matches any hits
        if(crtTrackMap.find(trackTrueID) != crtTrackMap.end()){
          hEffAngleTotal->Fill(angleCut);

          // If closest hit is below limit and track matches any hits then fill efficiency
          if(closestAngle.second < angleCut && closestAngle.second != -99999){
            hEffAngleReco->Fill(angleCut);
          }
        }

        // Fill total purity histogram with each cut if closest hit is below limit
        if(closestAngle.second < angleCut && closestAngle.second != -99999){
          hPurityAngleTotal->Fill(angleCut);

          // If closest hit is below limit and matched time is correct then fill purity
          double trackTime = closestAngle.first.ts1_ns * 1e-3;
          if(particles.find(trackTrueID) != particles.end()){
            if(std::abs(trackTime - trueTime) < 2.){
              hPurityAngleReco->Fill(angleCut);
            }
          }
        }
      }

      //----------------------------------------------------------------------------------------------------------
      //                                        SINGLE DCA CUT ANALYSIS
      //----------------------------------------------------------------------------------------------------------
      // Find the closest track by average distance of closest approach
      std::pair<crt::CRTTrack, double> closestDCA = trackAlg.ClosestCRTTrackByDCA(tpcTrack, crtTracks, event);
      if(closestDCA.second != -99999){
        hDCA->Fill(closestDCA.second);
      }
      if(closestDCA.second != -99999){
        int crtTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, closestDCA.first);
        if(crtTrueID == trackTrueID && crtTrueID != -99999){
          hMatchDCA->Fill(closestDCA.second);
        }
        else{
          hNoMatchDCA->Fill(closestDCA.second);
        }
      }

      int nbinsDCA = hEffDCATotal->GetNbinsX();
      for(int i = 0; i < nbinsDCA; i++){
        double DCAcut = hEffDCATotal->GetBinCenter(i);

        // Fill total efficiency histogram with each cut if track matches any hits
        if(crtTrackMap.find(trackTrueID) != crtTrackMap.end()){
          hEffDCATotal->Fill(DCAcut);

          // If closest hit is below limit and track matches any hits then fill efficiency
          if(closestDCA.second < DCAcut && closestDCA.second != -99999){
            hEffDCAReco->Fill(DCAcut);
          }
        }

        // Fill total purity histogram with each cut if closest hit is below limit
        if(closestDCA.second < DCAcut && closestDCA.second != -99999){
          hPurityDCATotal->Fill(DCAcut);

          // If closest hit is below limit and matched time is correct then fill purity
          double trackTime = closestDCA.first.ts1_ns * 1e-3;
          if(particles.find(trackTrueID) != particles.end()){
            if(std::abs(trackTime - trueTime) < 2.){
              hPurityDCAReco->Fill(DCAcut);
            }
          }
        }
      }

      
      //----------------------------------------------------------------------------------------------------------
      //                                    JOINT DCA AND ANGLE CUT ANALYSIS
      //----------------------------------------------------------------------------------------------------------
      std::vector<crt::CRTTrack> possTracks = trackAlg.AllPossibleCRTTracks(tpcTrack, crtTracks, event);
      for(auto const& possTrack : possTracks){
        int crtTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, possTrack);
        double angle = trackAlg.AngleBetweenTracks(tpcTrack, possTrack);
        double DCA = trackAlg.AveDCABetweenTracks(tpcTrack, possTrack, event);
        if(crtTrueID == trackTrueID && crtTrueID != -99999){
          hMatchAngleDCA->Fill(angle, DCA);
        }
        else{
          hNoMatchAngleDCA->Fill(angle, DCA);
        }
      }

    }

   /* 
    evd.SetDrawCrtTracks(true);
    evd.SetDrawCrtData(true);
    evd.SetDrawTpcTracks(true);
    evd.SetDrawTrueTracks(true);
    evd.SetDrawTpc(true);
    evd.SetTrueId(114);
    evd.Draw(event);
*/
    
  } // CRTTrackMatchingAna::analyze()


  void CRTTrackMatchingAna::endJob(){
  
    
  } // CRTTrackMatchingAna::endJob()

  
  DEFINE_ART_MODULE(CRTTrackMatchingAna)
} // namespace sbnd


