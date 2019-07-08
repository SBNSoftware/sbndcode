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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
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

#include "Pandora/PdgTable.h"

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

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
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

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    
  private:

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);
    
    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTTrackLabel;   ///< name of CRT producer
    art::InputTag fTPCTrackLabel; ///< name of CRT producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    CRTTrackMatchAlg trackAlg;

    CRTGeoAlg fCrtGeo;
    TPCGeoAlg fTpcGeo;

    CRTBackTracker fCrtBackTrack;
    CrtTrackCosmicTagAlg  ctTag;

    // Histograms
    std::vector<std::string> types {"NuMuTrack", "CrTrack", "DirtTrack", "NuTrack", "NuMuPfp", "CrPfp", "DirtPfp", "NuPfp"};
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
    , fPandoraLabel        (config().PandoraLabel())
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
      double trackTime = track->ts1_ns * 1e-3;
      if(trackTime > 0 && trackTime < 4) continue;
      numCrtTrackMap[trackTrueID]++;
    }

    // Get reconstructed tracks from the event
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    // Get PFParticles from pandora
    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);
    if( !pfParticleHandle.isValid() ){
      if(fVerbose) std::cout<<"Failed to find the PFParticles."<<std::endl;
      return;
    }
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
    // Get PFParticle to track associations
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTPCTrackLabel);

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
      if(truth->Origin() == simb::kBeamNeutrino){
        geo::Point_t vtx;
        vtx.SetX(truth->GetNeutrino().Nu().Vx()); vtx.SetY(truth->GetNeutrino().Nu().Vy()); vtx.SetZ(truth->GetNeutrino().Nu().Vz());
        // If neutrino vertex is not inside the TPC then call it a dirt particle
        if(!CosmicRemovalUtils::InFiducial(vtx, 0, 0)){ 
          dirtParticleIds.push_back(partID);
        }
        // If it's a primary muon
        else if(pdg==13 && particle.Mother()==0){ 
          lepParticleIds.push_back(partID);
        }
        // Other nu particles
        else{
          nuParticleIds.push_back(partID);
        }
      }

      // If origin is a cosmic ray
      else if(truth->Origin() == simb::kCosmicRay){
        crParticleIds.push_back(partID);
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
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trackTrueID) != lepParticleIds.end()) type = "NuMuTrack";
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trackTrueID) != nuParticleIds.end()) type = "NuTrack";
      if(std::find(crParticleIds.begin(), crParticleIds.end(), trackTrueID) != crParticleIds.end()) type = "CrTrack";
      if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trackTrueID) != dirtParticleIds.end()) type = "DirtTrack";
      if(type == "none") continue;

      if(numCrtTrackMap.find(trackTrueID) != numCrtTrackMap.end()){
        hNumTrueMatches[type]->Fill(numCrtTrackMap[trackTrueID]);
      }
      else{
        hNumTrueMatches[type]->Fill(0);
      }

      // Calculate t0 from CRT track matching
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


    //Loop over the pfparticle map
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){

      const art::Ptr<recob::PFParticle> pParticle(it->second);
      // Only look for primary particles
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;

      std::string type = "none";
      std::vector<recob::Track> nuTracks;

      // Loop over daughters of pfparticle
      for (const size_t daughterId : pParticle->Daughters()){

        // Get tracks associated with daughter
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pDaughter.key()));
        if(associatedTracks.size() != 1) continue;

        // Get the first associated track
        recob::Track tpcTrack = *associatedTracks.front();
        nuTracks.push_back(tpcTrack);

        // Truth match muon tracks and pfps
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);
        if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()){ 
          type = "NuMuPfp";
        }
        else if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()){ 
          if(type != "NuMuPfp") type = "NuPfp";
        }
        else if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()){ 
          if(type != "NuMuPfp" && type != "NuPfp") type = "DirtPfp";
        }
        else if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()){
          if(type != "NuMuPfp" && type != "NuPfp" && type != "DirtPfp") type = "CrPfp";
        }
      }

      if(nuTracks.size() == 0) continue;

      std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
                return left.Length() > right.Length();});

      recob::Track tpcTrack = nuTracks[0];
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(hits, false);

      if(numCrtTrackMap.find(trackTrueID) != numCrtTrackMap.end()){
        hNumTrueMatches[type]->Fill(numCrtTrackMap[trackTrueID]);
      }
      else{
        hNumTrueMatches[type]->Fill(0);
      }

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


  void CRTTrackCosmicTagAna::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  
  DEFINE_ART_MODULE(CRTTrackCosmicTagAna)
} // namespace sbnd


