////////////////////////////////////////////////////////////////////////
// Class:       StoppingCosmicIdAna
// Module Type: analyzer
// File:        StoppingCosmicIdAnaAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CosmicId/Algs/StoppingParticleCosmicIdAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <algorithm>

namespace sbnd {

  class StoppingCosmicIdAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };
      
      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of calorimetry producer data product")
      };

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Table<StoppingParticleCosmicIdAlg::Config> SPTagAlg {
        Name("SPTagAlg"),
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit StoppingCosmicIdAna(Parameters const& config);
 
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
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    art::InputTag fCaloModuleLabel; ///< name of calorimetry producer
    art::InputTag fPandoraLabel;
    bool          fVerbose;             ///< print information about what's going on

    StoppingParticleCosmicIdAlg spTag;
    TPCGeoAlg fTPCGeo;

    // histograms
    // True time of neutrino particles to get beam window
    std::vector<std::string> types {"NuMuTrack", "CrTrack", "DirtTrack", "NuTrack", "NuMuPfp", "CrPfp", "DirtPfp", "NuPfp"};
    std::map<std::string, TH1D*> hRatioTotal;
    std::map<std::string, TH1D*> hRatioTag;
    std::map<std::string, TH1D*> hStopLength;
    std::map<std::string, TH1D*> hNoStopLength;
    std::map<std::string, TH1D*> hStopChiSq;
    std::map<std::string, TH1D*> hNoStopChiSq;

  }; // class StoppingCosmicIdAna


  // Constructor
  StoppingCosmicIdAna::StoppingCosmicIdAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fCaloModuleLabel      (config().CaloModuleLabel())
    , fPandoraLabel        (config().PandoraLabel())
    , fVerbose              (config().Verbose())
    , spTag                 (config().SPTagAlg())
  {


  } // StoppingCosmicIdAna()


  void StoppingCosmicIdAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    for(auto type : types){
      hRatioTotal[type] = tfs->make<TH1D>(Form("%sRatioTotal", type.c_str()), "", 20, 1, 3);
      hRatioTag[type] = tfs->make<TH1D>(Form("%sRatioTag", type.c_str()), "", 20, 1, 3);
      hStopLength[type] = tfs->make<TH1D>(Form("%sStopLength", type.c_str()), "", 50, 0, 500);
      hNoStopLength[type] = tfs->make<TH1D>(Form("%sNoStopLength", type.c_str()), "", 50, 0, 500);
      hStopChiSq[type] = tfs->make<TH1D>(Form("%sStopChiSq", type.c_str()), "", 50, 0, 5);
      hNoStopChiSq[type] = tfs->make<TH1D>(Form("%sNoStopChiSq", type.c_str()), "", 50, 0, 5);
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- Stopping Particle Cosmic ID Ana Module -------------------"<<std::endl;

  }// StoppingCosmicIdAna::beginJob()


  void StoppingCosmicIdAna::analyze(const art::Event& event)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

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

    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

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
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);

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
        if(!fTPCGeo.InFiducial(vtx, 0, 0)){ 
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
    //                                    STOPPING CHI2 ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

    for(auto const& tpcTrack : (*tpcTrackHandle)){

      // Match to the true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
      std::string type = "none";
      if(std::find(lepParticleIds.begin(), lepParticleIds.end(), trueId) != lepParticleIds.end()) type = "NuMuTrack";
      if(std::find(nuParticleIds.begin(), nuParticleIds.end(), trueId) != nuParticleIds.end()) type = "NuTrack";
      if(std::find(crParticleIds.begin(), crParticleIds.end(), trueId) != crParticleIds.end()) type = "CrTrack";
      if(std::find(dirtParticleIds.begin(), dirtParticleIds.end(), trueId) != dirtParticleIds.end()) type = "DirtTrack";
      if(type == "none") continue;

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());
      if(calos.size()==0) continue;

      // Only focus on muon tracks or it'll be too hard
      if(std::abs(particles[trueId].PdgCode()) != 13) continue;

      geo::Point_t trueEnd {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
      bool stops = false;

      if(fTPCGeo.InFiducial(trueEnd, 0.)) stops = true;

      if(stops) hStopLength[type]->Fill(tpcTrack.Length());
      else hNoStopLength[type]->Fill(tpcTrack.Length());

      TVector3 trueEndVec (particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
      TVector3 start = tpcTrack.Vertex<TVector3>();
      TVector3 end = tpcTrack.End<TVector3>();
      geo::Point_t recoEnd = tpcTrack.End();
      if((start-trueEndVec).Mag() < (end-trueEndVec).Mag()) recoEnd = tpcTrack.Vertex();

      double chi2 = spTag.StoppingChiSq(recoEnd, calos);
      if(stops) hStopChiSq[type]->Fill(chi2);
      else hNoStopChiSq[type]->Fill(chi2);

      double startChi2 = spTag.StoppingChiSq(tpcTrack.Vertex(), calos);
      double endChi2 = spTag.StoppingChiSq(tpcTrack.End(), calos);

      int nbins = hRatioTotal.begin()->second->GetNbinsX();
      for(int i = 0; i < nbins; i++){
        double ratioCut = hRatioTotal.begin()->second->GetBinCenter(i);

        hRatioTotal[type]->Fill(ratioCut);

        // If closest hit is below limit and track matches any hits then fill efficiency
        if((startChi2 > ratioCut && !fTPCGeo.InFiducial(tpcTrack.End(), 5.) && tpcTrack.End().Y() > 0)
           || (endChi2 > ratioCut && !fTPCGeo.InFiducial(tpcTrack.Vertex(), 5.) && tpcTrack.Vertex().Y() > 0)){
          hRatioTag[type]->Fill(ratioCut);
        }

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
        int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
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
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);

      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(tpcTrack.ID());
      if(calos.size()==0) continue;

      // Only focus on muon tracks or it'll be too hard
      if(std::abs(particles[trueId].PdgCode()) != 13) continue;

      geo::Point_t trueEnd {particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ()};
      bool stops = false;

      if(fTPCGeo.InFiducial(trueEnd, 0.)) stops = true;

      if(stops) hStopLength[type]->Fill(tpcTrack.Length());
      else hNoStopLength[type]->Fill(tpcTrack.Length());

      TVector3 trueEndVec (particles[trueId].EndX(), particles[trueId].EndY(), particles[trueId].EndZ());
      TVector3 start = tpcTrack.Vertex<TVector3>();
      TVector3 end = tpcTrack.End<TVector3>();
      geo::Point_t recoEnd = tpcTrack.End();
      if((start-trueEndVec).Mag() < (end-trueEndVec).Mag()) recoEnd = tpcTrack.Vertex();

      double chi2 = spTag.StoppingChiSq(recoEnd, calos);
      if(stops) hStopChiSq[type]->Fill(chi2);
      else hNoStopChiSq[type]->Fill(chi2);

      double startChi2 = spTag.StoppingChiSq(tpcTrack.Vertex(), calos);
      double endChi2 = spTag.StoppingChiSq(tpcTrack.End(), calos);

      int nbins = hRatioTotal.begin()->second->GetNbinsX();
      for(int i = 0; i < nbins; i++){
        double ratioCut = hRatioTotal.begin()->second->GetBinCenter(i);

        hRatioTotal[type]->Fill(ratioCut);

        // If closest hit is below limit and track matches any hits then fill efficiency
        if((startChi2 > ratioCut && !fTPCGeo.InFiducial(tpcTrack.End(), 5.) && tpcTrack.End().Y() > 0)
           || (endChi2 > ratioCut && !fTPCGeo.InFiducial(tpcTrack.Vertex(), 5.) && tpcTrack.Vertex().Y() > 0)){
          hRatioTag[type]->Fill(ratioCut);
        }

      }

    }


  } // StoppingCosmicIdAna::analyze()


  void StoppingCosmicIdAna::endJob(){


  } // StoppingCosmicIdAna::endJob()

  void StoppingCosmicIdAna::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }

  DEFINE_ART_MODULE(StoppingCosmicIdAna)
} // namespace sbnd
