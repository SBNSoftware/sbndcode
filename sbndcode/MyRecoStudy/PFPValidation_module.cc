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

  std::map<int, int> GetTruePrimaryHits(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, const simb::MCParticle*>& trueParticles,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<recob::Hit>>& allHits) const;

  std::map<int, float> GetTruePrimaryEnergies(
      const std::map<int, const simb::MCParticle*>& trueParticles,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<sim::SimChannel>>& simchannels) const;

  std::map<int, float> GetTruePrimaryHitEnergies(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, const simb::MCParticle*>& trueParticles,
      const std::map<int, std::vector<int>>& truePrimaries,
      const std::vector<art::Ptr<recob::Hit>>& allHits) const;
  
  float GetTotalEnergyInHits(const detinfo::DetectorClocksData& clockData, const std::vector<art::Ptr<recob::Hit>>& hits) const;

  std::map<int, std::vector<std::vector<double>>> GetTruePrimaryHitXYZ(
      const detinfo::DetectorClocksData& clockData,
      const std::map<int, const simb::MCParticle*>& trueParticles,
      const std::map<int, std::vector<int>>& truePrimaries,
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

  // Fill the tree once per particle
  unsigned int fEventID;
  unsigned int fpfpID, pfpTrueID;

  //reco-to-truth 
  int pfpPdg, pfpTruePdg, pfpNumHits, eventNumPFP;
  float pfpTrueEnergy, pfpTrueMomentum, pfpTrueDepositedEnergy, pfpTrackScore;
  float pfpHitPurity, pfpHitComp, pfpEnergyPurity, pfpEnergyComp, pfpHitSPRatio;
  float G4TrueDepositedEnergy;
  
  //hit ans space point
  std::vector<double> spX, spY, spZ;
  std::vector<double> trackX, trackY, trackZ;
  std::vector<double> showerX, showerY, showerZ;
 
  //track info
  float trackE, trackLength;
  std::vector<double> trackDir;

  //shower info
  float showerOpenAngle, showerE, showerLength;
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

  std::cout << std::endl << "Processing Event " << fEventID << std::endl;

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
  
  // SIM CHANNELS
  auto const simChannelHandle(evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel));
  std::vector<art::Ptr<sim::SimChannel>> simChannels;
  art::fill_ptr_vector(simChannels, simChannelHandle);

  // HITS
  auto const hitHandle(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel));
  std::vector<art::Ptr<recob::Hit>> allHits;
  art::fill_ptr_vector(allHits, hitHandle);
  
  // Associations to PFPs
  art::FindManyP<recob::SpacePoint> fHitSP(allHits, evt, fPFPLabel);

  // Get map of true primary particle to number of reco hits / energy in reco hits
  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt));
  const std::map<int, int> truePrimaryHits(this->GetTruePrimaryHits(clockData, trueParticles, truePrimaries, allHits));
  const std::map<int, float> truePrimaryEnergies(this->GetTruePrimaryEnergies(trueParticles, truePrimaries, simChannels));
  const std::map<int, float> truePrimaryHitEnergies(this->GetTruePrimaryHitEnergies(clockData, trueParticles, truePrimaries, allHits));
//  const std::map<int, std::vector<std::vector<double>>> truePrimaryHitXYZ(this->GetTruePrimaryHitXYZ(clockData, trueParticles, truePrimaries, allHits));

  // PFPs
  auto const pfpHandle(evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  // Associations to PFPs
  art::FindManyP<recob::SpacePoint> fPFPSP(pfps, evt, fPFPLabel);
  art::FindManyP<recob::Cluster> fPFPClus(pfps, evt, fPFPLabel);
  art::FindManyP<recob::Track> fPFPTrack(pfps, evt, fTrackLabel);
  art::FindManyP<recob::Shower> fPFPShower(pfps, evt, fShowerLabel);

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

  eventNumPFP = pfps.size();

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

    fpfpID = pfp->Self();
    pfpPdg = pfp->PdgCode();
    pfpTrackScore = pfpTrackScoreMap.at(pfp->Self());
    pfpNumHits = pfpHits.size();
    pfpHitSPRatio = (float)sps.size() / pfpHits.size();

    if (fVerbose) {
      std::cout << std::endl << "PFP ID: " << pfp->Self() << " , pdg = " << pfp->PdgCode() << std::endl;
      std::cout << "PFP # of sp = " << sps.size() << std::endl;
    }

    //2 = collection plane
    const std::pair<int, double> trueId(ShowerUtils::TrueParticleIDFromTrueChain(clockData, truePrimaries, pfpHits, 2));
    if (trueId.first == -99999) {
      pfpTree->Fill();
      continue;
    }

    //true ID : # true hits
    const std::map<int, int> pfpTrueHitsMap(GetTruePrimaryHits(clockData, trueParticles, truePrimaries, pfpHits));
    //true ID : # true energy
    const std::map<int, float> pfpTrueEnergyMap(GetTruePrimaryHitEnergies(clockData, trueParticles, truePrimaries, pfpHits));
    //true ID: vector<XYZ>
    const std::map<int, std::vector<std::vector<double>>> pfpTruePrimaryHitXYZ(this->GetTruePrimaryHitXYZ(clockData, trueParticles, truePrimaries, pfpHits));

    const simb::MCParticle* trueParticle(trueParticles.at(trueId.first));

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

    // Track: 1 track to 1 PFP 
    const std::vector<art::Ptr<recob::Track>>& pfpTrack(fPFPTrack.at(pfp.key()));

    if (pfpTrack.size() == 1) {
      const art::Ptr<recob::Track>& track(pfpTrack.front());

      trackLength = track->Length();      
      trackDir.push_back(track->StartDirection().X());
      trackDir.push_back(track->StartDirection().Y());
      trackDir.push_back(track->StartDirection().Z());

      //Track Calo
      const std::vector<art::Ptr<anab::Calorimetry>>& calos(fTrackCalo.at(track.key()));
      unsigned int bestPlane = 999;
      unsigned int numHits = 0;    

      for (auto const& calo: calos) {
        unsigned int thisPlane = calo->PlaneID().Plane;
        unsigned int thisPlaneHits = calo->dEdx().size();
        if(thisPlaneHits > numHits) { 
          bestPlane = thisPlane;
          numHits = thisPlaneHits;
          trackE = calo->KineticEnergy();
        }
      }

      //Track Hit
      const std::vector<art::Ptr<recob::Hit>>& trackHits(fTrackHit.at(track.key()));
      unsigned int trackSPNum = 0;
      for (auto const& tHit: trackHits){
        const std::vector<art::Ptr<recob::SpacePoint>>& trackSPs(fHitSP.at(tHit.key()));
        trackSPNum += trackSPs.size();
        for(const auto& sp: trackSPs) {
          trackX.push_back(sp->XYZ()[0]);
          trackY.push_back(sp->XYZ()[1]);
          trackZ.push_back(sp->XYZ()[2]);
        } 

      }

      if (fVerbose) { 
        std::cout << std::endl;
        std::cout << "This pfp has 1 association with track" << std::endl;  
        std::cout << "# of SP points = " << trackSPNum << std::endl; 
        std::cout << "Track E = " << trackE << " on plane " << bestPlane << std::endl;
      }

    } else {
      if (fVerbose) std::cout << "This pfp has " << pfpTrack.size() << " association with track" << std::endl;  
    }

    // Shower: 1 track to 1 PFP
    const std::vector<art::Ptr<recob::Shower>>& pfpShower(fPFPShower.at(pfp.key()));
 
    if (pfpShower.size() == 1) { 
      const art::Ptr<recob::Shower>& shower(pfpShower.front());
 
      showerOpenAngle = shower->OpenAngle();
      showerLength = shower->Length();
      int planeID = shower->best_plane();
      showerE = shower->Energy()[planeID];
      showerDir.push_back(shower->Direction().X()); 
      showerDir.push_back(shower->Direction().Y()); 
      showerDir.push_back(shower->Direction().Z()); 
      
      const std::vector<art::Ptr<recob::SpacePoint>>& showerSPs(fShowerSP.at(shower.key()));
      for(const auto& sp: showerSPs) {
        showerX.push_back(sp->XYZ()[0]);
        showerY.push_back(sp->XYZ()[1]);
        showerZ.push_back(sp->XYZ()[2]);
      } 

      if (fVerbose) {
        std::cout << std::endl;
        std::cout << "This pfp has 1 association with shower" << std::endl;  
        std::cout << "# of shower sp = " << showerSPs.size() << std::endl;
        std::cout << "Shower E = " << showerE << " on plane " << planeID << std::endl;
      }
    } else {
      if (fVerbose) std::cout << "This pfp has " << pfpShower.size() << " association with shower" << std::endl;  
    } 

    pfpTree->Fill();

  } //End PFParticle

} //End of analyze

std::map<int, int> sbnd::PFPValidation::GetTruePrimaryHits(
    const detinfo::DetectorClocksData& clockData,
    const std::map<int, const simb::MCParticle*>& trueParticles,
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
    const std::map<int, const simb::MCParticle*>& trueParticles,
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
    const std::map<int, const simb::MCParticle*>& trueParticles,
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
    const std::map<int, const simb::MCParticle*>& trueParticles,
    const std::map<int, std::vector<int>>& truePrimaries,
    const std::vector<art::Ptr<recob::Hit>>& allHits) const
{  

  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  std::map<int, std::vector<std::vector<double>>> idToXYZMap;

  for (const art::Ptr<recob::Hit>& pHit : allHits){

    const std::vector<const sim::IDE*> IDEs(btServ->HitToSimIDEs_Ps(clockData, pHit));
    if (IDEs.empty()) continue;
     
    for (const auto& ide: IDEs) {
      std::vector<double> XYZ{ide->x, ide->y, ide->z};
      idToXYZEMap[std::abs(ide->trackID)].push_back(XYZ);
    }
  }

  std::map<int, std::vector<std::vector<double>>> truePrimaryHitXYZ;
  for (const auto& truePrimary : truePrimaries) {
    for (const auto& trueDaughter : truePrimary.second) {
      truePrimaryHitXYZ[truePrimary.first] = idToXYZMap[trueDaughter];
    }
  }
 
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
  pfpTree->Branch("showerOpenAngle", &showerOpenAngle);
  pfpTree->Branch("showerE", &showerE);
  pfpTree->Branch("showerLength", &showerLength);
  pfpTree->Branch("showerDir", &showerDir);
  pfpTree->Branch("trackE", &trackE);
  pfpTree->Branch("trackLength", &trackLength);
  pfpTree->Branch("trackDir", &trackDir);

  posTree = tfs->make<TTree>("posTree", "Tree with position info for each pfp");
  posTree->Branch("pfpspX", &spX);
  posTree->Branch("pfpspY", &spY);
  posTree->Branch("pfpspZ", &spZ);
  posTree->Branch("pfpshowerX", &spX);
  posTree->Branch("pfpshowerY", &spY);
  posTree->Branch("pfpshowerZ", &spZ);
//  posTree->Branch("pfptrueX", &pfptrueX);
//  posTree->Branch("pfptrueY", &pfptrueY);
//  posTree->Branch("pfptrueZ", &pfptrueZ);
//  posTree->Branch("trueX", &trueX);
//  posTree->Branch("trueY", &trueY);
//  posTree->Branch("trueZ", &trueZ);
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
  showerOpenAngle = -999;
  showerE = -999;
  showerLength = -999;
  showerDir.clear();
  trackE = -999;
  trackLength = -999;
  trackDir.clear();
}

void sbnd::PFPValidation::clearPosTree()
{
  spX.clear();
  spY.clear();
  spZ.clear();
  showerX.clear();
  showerY.clear();
  showerZ.clear();
  trackX.clear();
  trackY.clear();
  trackZ.clear();
}

DEFINE_ART_MODULE(sbnd::PFPValidation)
