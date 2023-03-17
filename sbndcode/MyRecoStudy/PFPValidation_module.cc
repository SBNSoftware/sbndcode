////////////////////////////////////////////////////////////////////////
// Class:       PFPValidation
// Plugin Type: analyzer (art v3_02_06)
// File:        PFPValidation_module.cc
//
// 8th March 2023: Lan Nguyen
// Generated at Wed Oct  2 03:27:09 2019 by Edward Tyley using cetskelgen
// from cetlib version v3_07_02.
// 09/10/2021 Adopted for ICARUS (Sergey Martynenko smartynen@bnl.gov)  
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbnci/Common/Modules/MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <numeric> 

namespace sbnd {
class PFPValidation;
}

class sbnd::PFPValidation : public art::EDAnalyzer {
  public:
  explicit PFPValidation(fhicl::ParameterSet const& pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPValidation(PFPValidation const&) = delete;
  PFPValidation(PFPValidation&&) = delete;
  PFPValidation& operator=(PFPValidation const&) = delete;
  PFPValidation& operator=(PFPValidation&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  void beginJob() override;

  void clearPFPTree();
  void clearPosTree();
  void clearTrueTree();

  std::map<int, int> GetTruePrimaryHits(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<recob::Hit>>& allHits) const;

  std::map<int, float> GetTruePrimaryEnergies(
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<sim::SimChannel>>& simchannels) const;

  std::map<int, float> GetTruePrimaryHitEnergies(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<recob::Hit>>& allHits) const;
  
  float GetTotalEnergyInHits(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit>>& hits) const;

  std::map<int, std::vector<std::vector<double>>> GetTruePrimaryHitXYZ(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<sim::SimChannel>>& simchannels,
      const std::vector<art::Ptr<recob::Hit>>& allHits) const;

  private:
  std::string fPFPLabel; 
  std::string fHitLabel;
  std::string fSimChannelLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fCaloLabel;
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
 
  bool fVerbose;

  TTree* pfpTree;
  TTree* posTree;
  TTree* trueTree;

  // Fill the tree once per particle
  unsigned int fEventID;

  //pfp
  int fpfpID, pfpTrueID;
  int pfpPdg, pfpTruePdg, pfpNumHits, eventNumPFP;
  std::vector<double> fVertex;

  //reco-to-truth 
  float pfpTrueEnergy, pfpTrueMomentum, pfpTrueDepositedEnergy, pfpTrackScore;
  float pfpHitPurity, pfpHitComp, pfpEnergyPurity, pfpEnergyComp, pfpHitSPRatio;
  float G4TrueDepositedEnergy;
  std::vector<double> pfpTrueVertex;
  float trueAngleBShower;

  //hit ans space point
  std::vector<double> spX, spY, spZ;
  int pfptrueID1, pfptrueID2;
  std::vector<double> pfptrueX1, pfptrueY1, pfptrueZ1; //matched to true shower 1 -> this is the matched shower for 1-to-1 matching
  std::vector<double> pfptrueX2, pfptrueY2, pfptrueZ2; //matched to true shower 2
  
  //true stuff
  int trueID;
  float trueE;
  std::vector<double> trueVertex;
  std::vector<double> trueX, trueY, trueZ;
 
  //track info
  float trackE, trackLength;
  std::vector<double> trackStart;
  std::vector<double> trackDir;

  //shower info
  float showerOpenAngle, showerE, showerLength;
  std::vector<double> showerStart;
  std::vector<double> showerDir;

};

sbnd::PFPValidation::PFPValidation(fhicl::ParameterSet const& pset)
    : EDAnalyzer { pset }
    , fPFPLabel(pset.get<std::string>("PFPLabel"))
    , fHitLabel(pset.get<std::string>("HitLabel"))
    , fSimChannelLabel(pset.get<std::string>("SimChannelLabel"))
    , fTrackLabel(pset.get<std::string>("TrackLabel"))
    , fShowerLabel(pset.get<std::string>("ShowerLabel"))
    , fCaloLabel(pset.get<std::string>("CaloLabel"))
    , fVerbose(pset.get<bool>("Verbose", false))
{
}

void sbnd::PFPValidation::analyze(art::Event const& evt)
{

  fEventID = evt.id().event();
  std::cout << std::setprecision(2) << std::fixed;

  std::cout << std::endl << ">>>>>>>>>>>>>>>>>>>Processing Event " << fEventID << std::endl;

  // Get the true g4 particles and make a map form trackId
  std::map<int, const simb::MCParticle*> trueParticles;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  for (auto const& particleIt : particles) {
    const simb::MCParticle* particle = particleIt.second;
    trueParticles[particle->TrackId()] = particle;
  }

  std::map<int, std::vector<int>> showerMothers = ShowerUtils::GetShowerMothersCandidates(trueParticles);

  std::map<int, std::vector<int>> truePrimaries;
  for (auto const& [trackId, particle] : trueParticles) {
    if (abs(particle->PdgCode()) == 11 || abs(particle->PdgCode()) == 22) {
      auto const& showerMotherIter(showerMothers.find(trackId));
      if (showerMotherIter != showerMothers.end()) {
        truePrimaries[trackId] = showerMotherIter->second;
      }
    } 
    else {
      truePrimaries[trackId] = { trackId };
    }
  }
  
  if (fVerbose) {
    std::cout << "Event True Particle Map" << std::endl;
    for (auto const& part: trueParticles) {
      std::cout << "particle g4id: " << part.first << ", pdg code: "<< part.second->PdgCode() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Event True Primaries Map" << std::endl;
    for (auto const& prim: truePrimaries) {
      std::cout << "primary g4id: " << prim.first << " had daugh g4 id : ";
      for (auto const& daug: prim.second) {
	std::cout << daug << " "  ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  // SIM CHANNELS
  auto const simChannelHandle(evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel));
  std::vector<art::Ptr<sim::SimChannel>> simChannels;
  art::fill_ptr_vector(simChannels, simChannelHandle);

  // HITS
  auto const hitHandle(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel));
  std::vector<art::Ptr<recob::Hit>> allHits;

  art::fill_ptr_vector(allHits, hitHandle);
  
  // Associations to Hits
  art::FindManyP<recob::SpacePoint> fHitSP(allHits, evt, fPFPLabel);

  // Get map of true primary particle to number of reco hits / energy in reco hits
  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt));
  const std::map<int, int> truePrimaryHits(this->GetTruePrimaryHits(clockData, truePrimaries, allHits));
  const std::map<int, float> truePrimaryEnergies(this->GetTruePrimaryEnergies( truePrimaries, simChannels));
  const std::map<int, float> truePrimaryHitEnergies(this->GetTruePrimaryHitEnergies(clockData, truePrimaries, allHits));

  const std::map<int, std::vector<std::vector<double>>> truePrimaryHitXYZ(this->GetTruePrimaryHitXYZ(clockData, truePrimaries, simChannels, allHits));

  // PFPs
  auto const pfpHandle(evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  // Associations to PFPs
  art::FindManyP<recob::SpacePoint> fPFPSP(pfps, evt, fPFPLabel);
  art::FindManyP<recob::Cluster> fPFPClus(pfps, evt, fPFPLabel);
  art::FindManyP<recob::Track> fPFPTrack(pfps, evt, fTrackLabel);
  art::FindManyP<recob::Shower> fPFPShower(pfps, evt, fShowerLabel);
  art::FindManyP<recob::Vertex> fPFPVertex(pfps, evt, fPFPLabel);

  // Tracks 
  auto const trackHandle(evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel));
  std::vector<art::Ptr<recob::Track>> allTracks;
  art::fill_ptr_vector(allTracks, trackHandle);

  // Association to Tracks
  art::FindManyP<anab::Calorimetry> fTrackCalo(allTracks, evt, fCaloLabel);
  art::FindManyP<recob::Hit> fTrackHit(allTracks, evt, fTrackLabel);

  // Showers 
  auto const showerHandle(evt.getValidHandle<std::vector<recob::Shower>>(fShowerLabel));
  std::vector<art::Ptr<recob::Shower>> allShowers;
  art::fill_ptr_vector(allShowers, showerHandle);

  // Association to Showers
  art::FindManyP<recob::SpacePoint> fShowerSP(allShowers, evt, fShowerLabel);

  // CLUSTERS
  auto const clusterHandle(evt.getValidHandle<std::vector<recob::Cluster>>(fPFPLabel));
  std::vector<art::Ptr<recob::Cluster>> clusters;
  art::fill_ptr_vector(clusters, clusterHandle);

  //Associations to Clusters
  art::FindManyP<recob::Hit> fClusHit(clusters, evt, fPFPLabel);

  // Create a map between PFParticles and their IDs
  art::FindManyP<larpandoraobj::PFParticleMetadata> fpfpmd(pfps, evt, fPFPLabel);
  if (!fpfpmd.isValid() || fpfpmd.size() == 0) return;
  
  std::map<long unsigned int, art::Ptr<recob::PFParticle>> pfpMap;
  std::map<long unsigned int, float> pfpTrackScoreMap;
  
  for (auto const& part: truePrimaryEnergies) {
  
    if ((part.first == 3) || (part.first == 4)) {
      clearTrueTree();
       
      trueID = part.first;
      trueE = part.second;
  
      const simb::MCParticle* trueP(trueParticles.at(trueID));
      trueVertex.push_back(trueP->Vx());
      trueVertex.push_back(trueP->Vy());
      trueVertex.push_back(trueP->Vz());
  
      //ALL true hits positions
      for (auto const it: truePrimaryHitXYZ) {
        if (it.first == trueID) {
          for (auto const pos: it.second) {
            trueX.push_back(pos[0]);
            trueY.push_back(pos[1]);
            trueZ.push_back(pos[2]);
          }
        }
      }

      trueTree->Fill();

      if (fVerbose) std::cout << "Filling particle g4id: " << part.first <<  ", has # hits = " <<  trueX.size() << std::endl;
    }
  }

  for (auto const& pfp : pfps) {
    long unsigned int pfpId(pfp->Self());
    pfpMap[pfpId] = pfp;
    const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec(fpfpmd.at(pfpId));
    for (auto const& pfpMeta : pfpMetaVec) {
      const larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMeta->GetPropertiesMap());
      auto const pfpTrackScoreIter(propertiesMap.find("TrackScore"));
      pfpTrackScoreMap[pfpId] = pfpTrackScoreIter == propertiesMap.end() ? -999 : pfpTrackScoreIter->second;
    }
  }

  std::map<int, double> pfpHitMap;

  for (auto const& pfp : pfps) {
    // Get all the associations with PFParticle
    const std::vector<art::Ptr<recob::Cluster>>& clusters(fPFPClus.at(pfp.key()));

    // Get the hits from the PFParticle
    std::vector<art::Ptr<recob::Hit>> pfpHits;

    for (const auto& cluster : clusters) {
      const std::vector<art::Ptr<recob::Hit>>& hits(fClusHit.at(cluster.key()));
      pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
    }
    pfpHitMap[pfp->Self()] = pfpHits.size();

  }

//  for (auto const it: pfpHitMap ) {
//    std::cout << "pfp id = " << it.first << " has #hits = " << it.second << std::endl;
//  }

  eventNumPFP = pfps.size() - 1;

  for (auto const& pfp : pfps) {

    clearPFPTree();
    clearPosTree();

    // Get all the associations with PFParticle
    const std::vector<art::Ptr<recob::Cluster>>& clusters(fPFPClus.at(pfp.key()));
    const std::vector<art::Ptr<recob::SpacePoint>>& sps(fPFPSP.at(pfp.key()));

    // Get the hits from the PFParticle
    std::vector<art::Ptr<recob::Hit>> pfpHits;

    for (const auto& cluster : clusters) {
      const std::vector<art::Ptr<recob::Hit>>& hits(fClusHit.at(cluster.key()));
      pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
    }

    if (pfpHits.empty()) continue;

    //saving top 2 pfparticle for now
    fpfpID = pfp->Self();
    if (fpfpID == 2) continue;

    pfpPdg = pfp->PdgCode();
    pfpTrackScore = pfpTrackScoreMap.at(pfp->Self());
    pfpNumHits = pfpHits.size();
    pfpHitSPRatio = (float)sps.size() / pfpHits.size();

    if (fVerbose) {
      std::cout << std::endl << "------------PFP ID: " << pfp->Self() << " , pdg = " << pfp->PdgCode() << std::endl;
      std::cout << "PFP # of sp = " << sps.size() << std::endl;
    }

    //2 = collection plane
    const std::pair<int, double> trueId(ShowerUtils::TrueParticleIDFromTrueChain(clockData, truePrimaries, pfpHits, 2));
    if (trueId.first == -99999) {
      pfpTree->Fill();
      continue;
    }

    //true ID : # true hits
    const std::map<int, int> pfpTrueHitsMap(GetTruePrimaryHits(clockData, truePrimaries, pfpHits));
    //true ID : # true energy
    const std::map<int, float> pfpTrueEnergyMap(GetTruePrimaryHitEnergies(clockData, truePrimaries, pfpHits));
    //true ID: vector<XYZ>
    const std::map<int, std::vector<std::vector<double>>> pfpTruePrimaryHitXYZ(this->GetTruePrimaryHitXYZ(clockData, truePrimaries, simChannels, pfpHits));

    const simb::MCParticle* trueParticle(trueParticles.at(trueId.first));

    // Some hard-coded shit here, I can't think
    // Matched to either G4ID 3 or 4
    pfptrueID1 = trueId.first;
    if (pfptrueID1 == 3) pfptrueID2 = 4;
    if (pfptrueID1 == 4) pfptrueID2 = 3;
  

    pfpTrueID = trueId.first;
    pfpTruePdg = trueParticle->PdgCode();
    pfpTrueEnergy = trueParticle->E();
    pfpTrueMomentum = trueParticle->P();

    //# of reco hits that are true
    const int pfpHitsTrueHits(pfpTrueHitsMap.at(trueId.first));
    const float pfpHitsTrueEnergy(pfpTrueEnergyMap.at(trueId.first));
    const float pfpHitsTotalEnergy(GetTotalEnergyInHits(clockData, pfpHits));

    // # of reco hits that matched to true / # of reco hits
    pfpHitPurity = (float)pfpHitsTrueHits / pfpHits.size();

    // # of reco hits that are true / # of true hits
    pfpHitComp = (float)pfpHitsTrueHits / truePrimaryHits.at(trueId.first);

    pfpEnergyPurity = pfpHitsTrueEnergy / pfpHitsTotalEnergy;
    pfpEnergyComp = pfpHitsTrueEnergy / truePrimaryHitEnergies.at(trueId.first);

    // Energy is in all planes, hence divide deposited energy by 3
    pfpTrueDepositedEnergy = truePrimaryHitEnergies.at(trueId.first) / 3;
    G4TrueDepositedEnergy = truePrimaryEnergies.at(trueId.first) / 3;
    
    //Save SpacePoint location
    for(const auto& sp: sps) {
      spX.push_back(sp->XYZ()[0]);
      spY.push_back(sp->XYZ()[1]);
      spZ.push_back(sp->XYZ()[2]);
    } 

    //Hit-Matched-To-PFP positions: true shower 1
    for (auto const it: pfpTruePrimaryHitXYZ) {
      if (it.first == pfptrueID1) {
        for (auto const pos: it.second) {
          pfptrueX1.push_back(pos[0]);
          pfptrueY1.push_back(pos[1]);
          pfptrueZ1.push_back(pos[2]);
        }
      }
    }
    
    //Hit-Matched-To-PFP positions: true shower 2
    for (auto const it: pfpTruePrimaryHitXYZ) {
      if (it.first == pfptrueID2) {
        for (auto const pos: it.second) {
          pfptrueX2.push_back(pos[0]);
          pfptrueY2.push_back(pos[1]);
          pfptrueZ2.push_back(pos[2]);
        }
      }
    }

    // True Vertex: pi0 g4id is 2
    const simb::MCParticle* truePi0(trueParticles.at(2));
    pfpTrueVertex.push_back(truePi0->Vx());
    pfpTrueVertex.push_back(truePi0->Vy());
    pfpTrueVertex.push_back(truePi0->Vz());

    const simb::MCParticle* truePh1(trueParticles.at(3));
    const simb::MCParticle* truePh2(trueParticles.at(4));
    auto const vec1 =  truePh1->Momentum().Vect();
    auto const vec2 =  truePh2->Momentum().Vect();

    float trueCosAngleBShower = (vec1.Dot(vec2)) /(vec1.Mag() * vec2.Mag());  
    trueAngleBShower = acos(trueCosAngleBShower)* 180.0 / M_PI;
    std::cout << "true angle between shower " << trueAngleBShower << std::endl;

    // Reco Vertex: 1 vtx to 1 PFP 
    const std::vector<art::Ptr<recob::Vertex>>& pfpVertex(fPFPVertex.at(pfp.key()));

    if (pfpVertex.size() == 1) {
      const art::Ptr<recob::Vertex>& vtx(pfpVertex.front());
      for (auto const& pos: vtx->position()) fVertex.push_back(pos);
    }
    
    // Track: 1 track to 1 PFP 
    const std::vector<art::Ptr<recob::Track>>& pfpTrack(fPFPTrack.at(pfp.key()));

    if (pfpTrack.size() == 1) {
      const art::Ptr<recob::Track>& track(pfpTrack.front());

      trackLength = track->Length();      
      trackStart.push_back(track->Start().X());
      trackStart.push_back(track->Start().Y());
      trackStart.push_back(track->Start().Z());
      trackDir.push_back(track->StartDirection().X());
      trackDir.push_back(track->StartDirection().Y());
      trackDir.push_back(track->StartDirection().Z());

      //Track Calo
      const std::vector<art::Ptr<anab::Calorimetry>>& calos(fTrackCalo.at(track.key()));
//      unsigned int bestPlane = 999;
      unsigned int numHits = 0;    

      for (auto const& calo: calos) {
//        unsigned int thisPlane = calo->PlaneID().Plane;
        unsigned int thisPlaneHits = calo->dEdx().size();
        if(thisPlaneHits > numHits) { 
//          bestPlane = thisPlane;
          numHits = thisPlaneHits;
          trackE = calo->KineticEnergy();
        }
      }

//      //Track Hit
//      const std::vector<art::Ptr<recob::Hit>>& trackHits(fTrackHit.at(track.key()));
//      unsigned int trackSPNum = 0;
//      for (auto const& tHit: trackHits){
//        const std::vector<art::Ptr<recob::SpacePoint>>& trackSPs(fHitSP.at(tHit.key()));
//        trackSPNum += trackSPs.size();
//        for(const auto& sp: trackSPs) {
//          trackX.push_back(sp->XYZ()[0]);
//          trackY.push_back(sp->XYZ()[1]);
//          trackZ.push_back(sp->XYZ()[2]);
//        } 
//      }
//
//      if (fVerbose) { 
//        std::cout << std::endl;
//        std::cout << "This pfp has 1 association with track" << std::endl;  
//        std::cout << "# of SP points = " << trackSPNum << std::endl; 
//        std::cout << "Track E = " << trackE << " on plane " << bestPlane << std::endl;
//      }

    } else {
//      if (fVerbose) std::cout << "This pfp has " << pfpTrack.size() << " association with track" << std::endl;  
    }

    // Shower: 1 track to 1 PFP
    const std::vector<art::Ptr<recob::Shower>>& pfpShower(fPFPShower.at(pfp.key()));
 
    if (pfpShower.size() == 1) { 
      const art::Ptr<recob::Shower>& shower(pfpShower.front());
 
      showerOpenAngle = shower->OpenAngle();
      showerLength = shower->Length();
      int planeID = shower->best_plane();
      showerE = shower->Energy()[planeID];
      showerStart.push_back(shower->ShowerStart().X());
      showerStart.push_back(shower->ShowerStart().Y());
      showerStart.push_back(shower->ShowerStart().Z());
      showerDir.push_back(shower->Direction().X()); 
      showerDir.push_back(shower->Direction().Y()); 
      showerDir.push_back(shower->Direction().Z()); 
      
//      const std::vector<art::Ptr<recob::SpacePoint>>& showerSPs(fShowerSP.at(shower.key()));
//      for(const auto& sp: showerSPs) {
//        showerX.push_back(sp->XYZ()[0]);
//        showerY.push_back(sp->XYZ()[1]);
//        showerZ.push_back(sp->XYZ()[2]);
//      } 
//
//      if (fVerbose) {
//        std::cout << std::endl;
//        std::cout << "This pfp has 1 association with shower" << std::endl;  
//        std::cout << "# of shower sp = " << showerSPs.size() << std::endl;
//        std::cout << "Shower E = " << showerE << " on plane " << planeID << std::endl;
//      }
    } else {
//      if (fVerbose) std::cout << "This pfp has " << pfpShower.size() << " association with shower" << std::endl;  
    } 

    pfpTree->Fill();
    posTree->Fill();

  } //End PFParticle

} //End of analyze

std::map<int, int> sbnd::PFPValidation::GetTruePrimaryHits(
    const detinfo::DetectorClocksData& clockData,
    const std::map<int, std::vector<int>>& truePrimaries,
    const std::vector<art::Ptr<recob::Hit>>& allHits) const
{

  std::map<int, int> trueParticleHits;
  for (const auto& hit : allHits) {
    const int trackID(TruthMatchUtils::TrueParticleID(clockData, hit, true));
    ++trueParticleHits[trackID];
  }

  std::map<int, int> truePrimaryHits;
  for (const auto& truePrimary : truePrimaries) {
    for (const auto& trueDaughter : truePrimary.second) {
      truePrimaryHits[truePrimary.first] += trueParticleHits[trueDaughter];
    }
  }
  return truePrimaryHits;
}

std::map<int, float> sbnd::PFPValidation::GetTruePrimaryEnergies(
    const std::map<int, std::vector<int>>& truePrimaries,
    const std::vector<art::Ptr<sim::SimChannel>>& simChannels) const
{

  std::map<int, float> trueParticleEnergies;
  for (const auto& simChannel : simChannels) {
    const auto tdcideMap(simChannel->TDCIDEMap());
    for (const auto& [tdc, ideVec] : tdcideMap) {
      for (const auto& ide : ideVec) {
        trueParticleEnergies[std::abs(ide.trackID)] += ide.energy;
      }
    }
  }

  std::map<int, float> truePrimaryEnergies;
  for (const auto& truePrimary : truePrimaries) {
    for (const auto& trueDaughter : truePrimary.second) {
      truePrimaryEnergies[truePrimary.first] += trueParticleEnergies[trueDaughter];
    }
  }
  return truePrimaryEnergies;
}

std::map<int, float> sbnd::PFPValidation::GetTruePrimaryHitEnergies(
    const detinfo::DetectorClocksData& clockData,
    const std::map<int, std::vector<int>>& truePrimaries,
    const std::vector<art::Ptr<recob::Hit>>& allHits) const
{

  TruthMatchUtils::IDToEDepositMap idToEDepMap;
  for (const art::Ptr<recob::Hit>& pHit : allHits)
    TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepMap, clockData, pHit, true);

  std::map<int, float> truePrimaryHitEnergies;
  for (const auto& truePrimary : truePrimaries) {
    for (const auto& trueDaughter : truePrimary.second) {
      truePrimaryHitEnergies[truePrimary.first] += idToEDepMap[trueDaughter];
    }
  }
  return truePrimaryHitEnergies;
}

float sbnd::PFPValidation::GetTotalEnergyInHits(
    const detinfo::DetectorClocksData& clockData,
    const std::vector<art::Ptr<recob::Hit>>& hits) const
{

  TruthMatchUtils::IDToEDepositMap idToEDepMap;
  for (const art::Ptr<recob::Hit>& pHit : hits)
    TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepMap, clockData, pHit, true);

  return std::accumulate(idToEDepMap.cbegin(), idToEDepMap.cend(), 0.f,
      [](float sum, auto const& iter) { return sum + iter.second; });
}

std::map<int, std::vector<std::vector<double>>> sbnd::PFPValidation::GetTruePrimaryHitXYZ(
    const detinfo::DetectorClocksData& clockData,
    const std::map<int, std::vector<int>>& truePrimaries,
    const std::vector<art::Ptr<sim::SimChannel>>& simChannels,
    const std::vector<art::Ptr<recob::Hit>>& allHits
    ) const
{  

  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  std::map<int, std::vector<std::vector<double>>> idToXYZMap;
//  std::map<int, float> idToEMap;

  unsigned int i = 0; 
  for (const art::Ptr<recob::Hit>& pHit : allHits){
   
    //Check if hit has valid Channel
    bool foundChan = false;
    for (auto const& sc: simChannels){
      if(sc->Channel() == pHit->Channel()) {
        foundChan = true;
        continue;
      }
    }

    if (foundChan) {
      const std::vector<const sim::IDE*> IDEs(btServ->HitToSimIDEs_Ps(clockData, pHit));
      if (IDEs.empty()) continue;
     
      for (const auto& ide: IDEs) {
        std::vector<double> XYZ;
	XYZ.push_back(ide->x);
	XYZ.push_back(ide->y);
	XYZ.push_back(ide->z);
        idToXYZMap[std::abs(ide->trackID)].push_back(XYZ);
//        idToEMap[std::abs(ide->trackID)] += ide->energy;
      }
    }
    i++;
  }

  std::map<int, std::vector<std::vector<double>>> truePrimaryHitXYZ;
  std::map<int, float> truePrimaryHitE;
  for (const auto& truePrimary : truePrimaries) {
    for (const auto& trueDaughter : truePrimary.second) {

//      truePrimaryHitE[truePrimary.first] += idToEMap[trueDaughter];

      for (const auto& xyz: idToXYZMap[trueDaughter]) {
        truePrimaryHitXYZ[truePrimary.first].push_back(xyz);
      }
    }
  }

//  std::cout << "My True E Map" << std::endl;
//  for (auto const& prim: idToEMap) {
//    if(prim.second > 0) std::cout << "primary g4id: " << prim.first << " had E : " << prim.second << std::endl;
//  }
//  std::cout << std::endl;
//
//  std::cout << "My True Primaries E" << std::endl;
//  for (auto const& prim: truePrimaryHitE) {
//    if(prim.second > 0) std::cout << "primary g4id: " << prim.first << " had E : " << prim.second << std::endl;
//  }
//  std::cout << std::endl;
//
//  std::cout << "My True XYZ Map" << std::endl;
//  for (auto const& prim: idToXYZMap) {
//    if(prim.second.size() > 0 )std::cout << "primary g4id: " << prim.first << " had #hits : " << prim.second.size() << std::endl;
//  }
//  std::cout << std::endl;
//
//  std::cout << "My True Primary XYZ" << std::endl;
//  for (auto const& prim: truePrimaryHitXYZ) {
//    if(prim.second.size()) std::cout << "primary g4id: " << prim.first << " had #hits : " << prim.second.size() << std::endl;
//  }
//  std::cout << std::endl;

  return truePrimaryHitXYZ;

}


void sbnd::PFPValidation::beginJob()
{
  pfpTree = tfs->make<TTree>("pfpTree", "Tree with metrics for each pfp");

  pfpTree->Branch("eventID", &fEventID);
  pfpTree->Branch("eventNumPFP", &eventNumPFP); 
  pfpTree->Branch("pfpTrueID", &pfpTrueID);
  pfpTree->Branch("pfpID", &fpfpID);
  pfpTree->Branch("pfpPdg", &pfpPdg);
  pfpTree->Branch("pfpTruePdg", &pfpTruePdg);
  pfpTree->Branch("pfpNumHits", &pfpNumHits);
  pfpTree->Branch("pfpTrueEnergy", &pfpTrueEnergy);
  pfpTree->Branch("pfpTrueMomentum", &pfpTrueMomentum);
  pfpTree->Branch("pfpTrueDepositedEnergy", &pfpTrueDepositedEnergy);
  pfpTree->Branch("G4TrueDepositedEnergy", &G4TrueDepositedEnergy);
  pfpTree->Branch("pfpTrackScore", &pfpTrackScore);
  pfpTree->Branch("pfpHitPurity", &pfpHitPurity);
  pfpTree->Branch("pfpHitComp", &pfpHitComp);
  pfpTree->Branch("pfpEnergyPurity", &pfpEnergyPurity);
  pfpTree->Branch("pfpEnergyComp", &pfpEnergyComp);
  pfpTree->Branch("pfpHitSPRatio", &pfpHitSPRatio);
  pfpTree->Branch("trueAngleBShower", &trueAngleBShower);
  pfpTree->Branch("showerOpenAngle", &showerOpenAngle);
  pfpTree->Branch("showerE", &showerE);
  pfpTree->Branch("showerLength", &showerLength);
  pfpTree->Branch("showerStart", &showerStart);
  pfpTree->Branch("showerDir", &showerDir);
  pfpTree->Branch("trackE", &trackE);
  pfpTree->Branch("trackLength", &trackLength);
  pfpTree->Branch("trackStart", &trackStart);
  pfpTree->Branch("trackDir", &trackDir);

  posTree = tfs->make<TTree>("posTree", "Tree with position info for each pfp");
  posTree->Branch("eventID", &fEventID);
  posTree->Branch("eventNumPFP", &eventNumPFP); 
  posTree->Branch("pfpID", &fpfpID);
  posTree->Branch("pfptrueID1", &pfptrueID1);
  posTree->Branch("pfptrueID2", &pfptrueID2);
  posTree->Branch("pfptrueVertex", &pfpTrueVertex);
  posTree->Branch("pfpVertex", &fVertex);
  posTree->Branch("pfpspX", &spX);
  posTree->Branch("pfpspY", &spY);
  posTree->Branch("pfpspZ", &spZ);
  posTree->Branch("pfptrueX1", &pfptrueX1);
  posTree->Branch("pfptrueY1", &pfptrueY1);
  posTree->Branch("pfptrueZ1", &pfptrueZ1);
  posTree->Branch("pfptrueX2", &pfptrueX2);
  posTree->Branch("pfptrueY2", &pfptrueY2);
  posTree->Branch("pfptrueZ2", &pfptrueZ2);

  trueTree = tfs->make<TTree>("trueTree", "Tree with info for each true particle");
  trueTree->Branch("eventID", &fEventID);
  trueTree->Branch("trueID", &trueID);
  trueTree->Branch("trueE", &trueE);
  trueTree->Branch("trueVertex", &trueVertex);
  trueTree->Branch("trueX", &trueX);
  trueTree->Branch("trueY", &trueY);
  trueTree->Branch("trueZ", &trueZ);
}

void sbnd::PFPValidation::clearPFPTree()
{
//  fEventID = -999;
//  eventNumPFP = -999; 
  pfpTrueID = -999;
  fpfpID = -999;
  pfpPdg = -999;
  pfpTruePdg = -999;
  pfpTrueEnergy = -999;
  pfpTrueMomentum = -999;
  pfpTrueDepositedEnergy = -999;
  pfpTrackScore = -999;
  pfpNumHits = -999;
  pfpHitPurity = -999;
  pfpHitComp = -999;
  pfpEnergyPurity = -999;
  pfpEnergyComp = -999;
  pfpHitSPRatio = -999;
  trueAngleBShower = -999;
  showerOpenAngle = -999;
  showerE = -999;
  showerLength = -999;
  showerStart.clear();
  showerDir.clear();
  trackE = -999;
  trackLength = -999;
  trackStart.clear();
  trackDir.clear();
}

void sbnd::PFPValidation::clearPosTree()
{
//  fEventID = -999;
//  eventNumPFP = -999; 
  fpfpID = -999;
  pfptrueID1 =  -999;
  pfptrueID2 =  -999;
  fVertex.clear();
  pfpTrueVertex.clear();
  spX.clear();
  spY.clear();
  spZ.clear();
  pfptrueX1.clear();
  pfptrueY1.clear();
  pfptrueZ1.clear();
  pfptrueX2.clear();
  pfptrueY2.clear();
  pfptrueZ2.clear();
}

void sbnd::PFPValidation::clearTrueTree()
{
//  fEventID = -999;
  trueID = -999;
  trueE = -999;
  trueVertex.clear();
  trueX.clear();
  trueY.clear();
  trueZ.clear();
}
DEFINE_ART_MODULE(sbnd::PFPValidation)
