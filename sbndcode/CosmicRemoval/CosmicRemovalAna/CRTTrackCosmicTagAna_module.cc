////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackCosmicTagAna
// Module Type: analyzer
// File:        CRTTrackCosmicTagAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CrtTrackCosmicTagAlg.h"
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

  class CRTTrackCosmicTagAna : public art::EDAnalyzer {
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

      fhicl::Table<CRTTrackMatchAlg::Config> TrackMatchAlg {
        Name("TrackMatchAlg"),
      };

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };

      fhicl::Table<CrtTrackCosmicTagAlg::Config> CTTagAlg {
        Name("CTTagAlg"),
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTTrackCosmicTagAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

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
    CrtTrackCosmicTagAlg  ctTag;

    // Histograms
    std::vector<std::string> types {"NuMu", "Cr", "Dirt", "Nu"};
    std::map<std::string, TH1D*> hNumTrueMatches;
    std::map<std::string, TH1D*> hMatchDCA;
    std::map<std::string, TH1D*> hNoMatchDCA;
    std::map<std::string, TH1D*> hMatchAngle;
    std::map<std::string, TH1D*> hNoMatchAngle;
    std::map<std::string, TH1D*> hDCATotal;
    std::map<std::string, TH1D*> hDCATag;
    std::map<std::string, TH1D*> hLengthTotal;
    std::map<std::string, TH1D*> hLengthTag;

  }; // class CRTTrackCosmicTagAna


  // Constructor
  CRTTrackCosmicTagAna::CRTTrackCosmicTagAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCRTTrackLabel       (config().CRTTrackLabel())
    , fTPCTrackLabel       (config().TPCTrackLabel())
    , fVerbose             (config().Verbose())
    , trackAlg             (config().TrackMatchAlg())
    , fCrtBackTrack        (config().CrtBackTrack())
    , ctTag                (config().CTTagAlg())
  {

  } //CRTTrackCosmicTagAna()


  void CRTTrackCosmicTagAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    for(auto type : types){
      hNumTrueMatches[type] = tfs->make<TH1D>(Form("%sNumTrueMatches", type.c_str()), "", 10, 0, 10);
      hMatchDCA[type]       = tfs->make<TH1D>(Form("%sMatchDCA", type.c_str()),       "", 60, 0, 100);
      hNoMatchDCA[type]     = tfs->make<TH1D>(Form("%sNoMatchDCA", type.c_str()),     "", 60, 0, 100);
      hMatchAngle[type]     = tfs->make<TH1D>(Form("%sMatchAngle", type.c_str()),     "", 60, 0, 3.2);
      hNoMatchAngle[type]   = tfs->make<TH1D>(Form("%sNoMatchAngle", type.c_str()),   "", 60, 0, 3.2);
      hDCATotal[type]       = tfs->make<TH1D>(Form("%sDCATotal", type.c_str()),       "", 20, 0, 80);
      hDCATag[type]         = tfs->make<TH1D>(Form("%sDCATag", type.c_str()),         "", 20, 0, 80);
      hLengthTotal[type]    = tfs->make<TH1D>(Form("%sLengthTotal", type.c_str()),    "", 20, 0, 400);
      hLengthTag[type]      = tfs->make<TH1D>(Form("%sLengthTag", type.c_str()),      "", 20, 0, 400);
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT T0 Matching Ana Module -------------------"<<std::endl;

  } // CRTTrackCosmicTagAna::beginJob()


  void CRTTrackCosmicTagAna::analyze(const art::Event& event)
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
    art::Handle< std::vector<crt::CRTTrack>> crtTrackHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackHandle);

    fCrtBackTrack.Initialize(event);
    std::vector<crt::CRTTrack> crtTracks;
    std::map<int, int> numCrtTrackMap;
    int track_i = 0;
    for(auto const& track : (crtTrackList)){
      crtTracks.push_back(*track);
      int trackTrueID = fCrtBackTrack.TrueIdFromTrackId(event, track_i);
      track_i++;
      numCrtTrackMap[trackTrueID]++;
    }

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    std::map<int, int> numTpcTrackMap;
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      numTpcTrackMap[trackTrueID]++;
    }


    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    
    std::map<int, simb::MCParticle> particles;
    std::vector<int> nuParticleIds;
    std::vector<int> lepParticleIds;
    std::vector<int> dirtParticleIds;
    std::vector<int> crParticleIds;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;

      // Get MCTruth
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(partID);
      int pdg = std::abs(particle.PdgCode());

      // If origin is a neutrino
      std::string type = "none";
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        // If neutrino vertex is not inside the TPC then call it a dirt particle
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)){ 
          dirtParticleIds.push_back(partID);
          type = "Dirt";
        }
        // If it's a primary muon
        else if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partID);
          type = "NuMu";
        }
        // Other nu particles
        else{
          nuParticleIds.push_back(partID);
          type = "Nu";
        }
      }

      // If origin is a cosmic ray
      else if(truth->Origin() == simb::kCosmicRay){
        crParticleIds.push_back(partID);
        type = "Cr";
      }

      if(type == "none") continue;
      if(numTpcTrackMap.find(partID) == numTpcTrackMap.end()) continue;
      if(numCrtTrackMap.find(partID) != numCrtTrackMap.end()){
        hNumTrueMatches[type]->Fill(numCrtTrackMap[partID]);
      }
      else{
        hNumTrueMatches[type]->Fill(0);
      }
      
    }

    //----------------------------------------------------------------------------------------------------------
    //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    // Loop over reconstructed tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
      std::string type = "none";
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trackTrueID) != lepParticleIds.end()) type = "NuMu";
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trackTrueID) != nuParticleIds.end()) type = "Nu";
      if(std::find(crParticleIds.begin(), crParticleIds.end(), trackTrueID) != crParticleIds.end()) type = "Cr";
      if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trackTrueID) != dirtParticleIds.end()) type = "Dirt";
      if(type == "none") continue;

      // Calculate t0 from CRT Hit matching
      std::pair<crt::CRTTrack, double> closestAngle = trackAlg.ClosestCRTTrackByAngle(tpcTrack, crtTracks, event);
      std::pair<crt::CRTTrack, double> closestDCA = trackAlg.ClosestCRTTrackByDCA(tpcTrack, crtTracks, event);

      if(closestAngle.second != -99999){
        int crtTrackTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, closestAngle.first);
        if(crtTrackTrueID == trackTrueID && crtTrackTrueID != -99999){
          hMatchAngle[type]->Fill(closestAngle.second);
        }
        else{
          hNoMatchAngle[type]->Fill(closestAngle.second);
        }
      }

      if(closestDCA.second != -99999){
        int crtTrackTrueID = fCrtBackTrack.TrueIdFromTotalEnergy(event, closestDCA.first);
        if(crtTrackTrueID == trackTrueID && crtTrackTrueID != -99999){
          hMatchDCA[type]->Fill(closestDCA.second);
        }
        else{
          hNoMatchDCA[type]->Fill(closestDCA.second);
        }
      }
        

      int nbins = hDCATotal.begin()->second->GetNbinsX();
      for(int i = 0; i < nbins; i++){
        double DCAcut = hDCATotal.begin()->second->GetBinCenter(i);

        hDCATotal[type]->Fill(DCAcut);

        // If closest hit is below limit and track matches any hits then fill efficiency
        if(closestDCA.second < DCAcut && closestDCA.second != -99999){
          double trackTime = closestDCA.first.ts1_ns * 1e-3;
          if(trackTime > 0 && trackTime < 4) continue;
          hDCATag[type]->Fill(DCAcut);
        }

      }

      hLengthTotal[type]->Fill(tpcTrack.Length());
      if(ctTag.CrtTrackCosmicTag(tpcTrack, crtTracks, event)){
        hLengthTag[type]->Fill(tpcTrack.Length());
      }
    }

    
  } // CRTTrackCosmicTagAna::analyze()


  void CRTTrackCosmicTagAna::endJob(){
  
    
  } // CRTTrackCosmicTagAna::endJob()

  
  DEFINE_ART_MODULE(CRTTrackCosmicTagAna)
} // namespace sbnd


