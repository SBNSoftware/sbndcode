////////////////////////////////////////////////////////////////////////
// Class:       PFPValidation
// Plugin Type: analyzer (art v3_02_06)
// File:        PFPValidation_module.cc
//
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
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  // Declare member data here.
  unsigned int maxLabelLength;

  TTree* trueTree;
  TTree* pfpTree;
  TTree* eventTree;

  // Fill the tree once per particle
  unsigned int fEventID;
  unsigned int fpfpID;

  std::string pfpModuleLabel;
  int pfpPdg, pfpTruePdg, pfpNumHits;
  float pfpTrueEnergy, pfpTrueMomentum, pfpTrueDepositedEnergy, pfpTrackScore;
  float pfpHitPurity, pfpHitComp, pfpEnergyPurity, pfpEnergyComp, pfpHitSPRatio;
  float G4TrueDepositedEnergy;
  
  int eventNumPFP;
  std::vector<double> spX, spY, spZ;
  std::vector<double> pfptrueX, pfptrueY, pfptrueZ;
  std::vector<double> trueX, trueY, trueZ;
};

sbnd::PFPValidation::PFPValidation(fhicl::ParameterSet const& pset)
    : EDAnalyzer { pset }
    , fPFPLabel(pset.get<std::string>("PFPLabel"))
    , fHitLabel(pset.get<std::string>("HitLabel"))
    , fSimChannelLabel(pset.get<std::string>("SimChannelLabel"))
    , maxLabelLength(0)
{
}

void sbnd::PFPValidation::analyze(art::Event const& evt)
{

  fEventID = evt.id().event();
  std::cout << std::setprecision(2) << std::fixed;

  std::cout << std::endl; 
  std::cout << "Processing Event " << fEventID << std::endl;

  // Get the true g4 particles and make a map form trackId
  std::map<int, const simb::MCParticle*> trueParticles;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  for (auto const& particleIt : particles) {
    const simb::MCParticle* particle = particleIt.second;
    trueParticles[particle->TrackId()] = particle;
  }

  std::cout << "trueParticles map size = " << trueParticles.size() << std::endl;
  std::cout << "--------------------------" << std::endl;
 
   
  for (const auto& tp : trueParticles) {
    std::cout << "trueParticle trk ID: " << tp.first <<std::endl;
    std::cout << "trueParticle pdg: " << tp.second->PdgCode() << std::endl;
  }
  std::cout << "--------------------------" << std::endl;
  

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
  
  std::cout << "truePrimaries map size = " << trueParticles.size() << std::endl;

  for (const auto& truePrimary : truePrimaries) {
    std::cout << "truePrimary trk ID: " << truePrimary.first <<std::endl;
    for (const auto& trueDaughter : truePrimary.second) {
      std::cout << "trueDaughter trk ID: " << trueDaughter <<std::endl;
    }
    std::cout << "--------------------------" << std::endl;
  }
  // SIM CHANNELS
  auto const simChannelHandle(evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel));
  std::vector<art::Ptr<sim::SimChannel>> simChannels;
  art::fill_ptr_vector(simChannels, simChannelHandle);

  // HITS
  auto const hitHandle(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel));
  std::vector<art::Ptr<recob::Hit>> allHits;
  art::fill_ptr_vector(allHits, hitHandle);

  // Get map of true primary particle to number of reco hits / energy in reco hits
  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt));
  const std::map<int, int> truePrimaryHits(this->GetTruePrimaryHits(clockData, trueParticles, truePrimaries, allHits));
  const std::map<int, float> truePrimaryEnergies(this->GetTruePrimaryEnergies(trueParticles, truePrimaries, simChannels));
  const std::map<int, float> truePrimaryHitEnergies(this->GetTruePrimaryHitEnergies(clockData, trueParticles, truePrimaries, allHits));
  const std::map<int, std::vector<std::vector<double>>> truePrimaryHitXYZ(this->GetTruePrimaryHitXYZ(clockData, trueParticles, truePrimaries, allHits));

  // PFPs
  auto const pfpHandle(evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);
  
  // Associations to PFPs
  art::FindManyP<recob::SpacePoint> fmPFPSP(pfps, evt, fPFPLabel);
  art::FindManyP<recob::Cluster> fmPFPClus(pfps, evt, fPFPLabel);
  
  // CLUSTERS
  auto const clusterHandle(evt.getValidHandle<std::vector<recob::Cluster>>(fPFPLabel));
  std::vector<art::Ptr<recob::Cluster>> clusters;
  art::fill_ptr_vector(clusters, clusterHandle);

  //Associations to Clusters
  art::FindManyP<recob::Hit> fmClusHit(clusters, evt, fPFPLabel);

  // Create a map between PFParticles and their IDs
  art::FindManyP<larpandoraobj::PFParticleMetadata> fmpfpmd(pfps, evt, fPFPLabel);
  if (!fmpfpmd.isValid() || fmpfpmd.size() == 0) return;
  
  std::map<long unsigned int, art::Ptr<recob::PFParticle>> pfpMap;
  std::map<long unsigned int, float> pfpTrackScoreMap;

  for (auto const& pfp : pfps) {
    long unsigned int pfpId(pfp->Self());
    pfpMap[pfpId] = pfp;
    const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec(fmpfpmd.at(pfpId));
    for (auto const& pfpMeta : pfpMetaVec) {
      const larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMeta->GetPropertiesMap());
      auto const pfpTrackScoreIter(propertiesMap.find("TrackScore"));
      pfpTrackScoreMap[pfpId] = pfpTrackScoreIter == propertiesMap.end() ? -999 : pfpTrackScoreIter->second;
    }
  }

  eventNumPFP = pfps.size();

  for (auto const& pfp : pfps) {

    clearPFPTree();

    // Get the hits from the PFParticle
    std::vector<art::Ptr<recob::Hit>> pfpHits;
    const std::vector<art::Ptr<recob::Cluster>>& clusters(fmPFPClus.at(pfp.key()));
    for (const auto& cluster : clusters) {
      const std::vector<art::Ptr<recob::Hit>>& hits(fmClusHit.at(cluster.key()));
      pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
    }

    const std::vector<art::Ptr<recob::SpacePoint>>& sps(fmPFPSP.at(pfp.key()));

    if (pfpHits.empty()) continue;

    fpfpID = pfp->Self();
    pfpPdg = pfp->PdgCode();
    pfpTrackScore = pfpTrackScoreMap.at(pfp->Self());
    pfpNumHits = pfpHits.size();
    pfpHitSPRatio = (float)sps.size() / pfpHits.size();

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

    //Save all true hit xyz 
    
    //Save pfp-true hit xyz 
    for (const auto& pfpTrue: pfpTruePrimaryHitXYZ) {
      const auto g4ID = pfpTrue.first;
      if(g4ID != trueId.first) continue;      
      for (const auto& XYZ : pfpTrue.second) {
        pfptrueX.push_back(XYZ[0]);
	pfptrueY.push_back(XYZ[1]);
        pfptrueZ.push_back(XYZ[2]);
      }
    }

    //Save pfp-true hit xyz 
    for (const auto& pfpTrue: truePrimaryHitXYZ) {
      const auto g4ID = pfpTrue.first;
      if(g4ID != trueId.first) continue;      
      for (const auto& XYZ : pfpTrue.second) {
        trueX.push_back(XYZ[0]);
	trueY.push_back(XYZ[1]);
        trueZ.push_back(XYZ[2]);
      }
    }
    
    //Save SpacePoint location
    for(const auto& sp: sps) {
      spX.push_back(sp->XYZ()[0]);
      spY.push_back(sp->XYZ()[1]);
      spZ.push_back(sp->XYZ()[2]);
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

    const std::vector<sim::TrackIDE> trackIDEs(btServ->HitToTrackIDEs(clockData, pHit));

    if (trackIDEs.empty()) continue;

    for (const sim::TrackIDE& trackIDE : trackIDEs) {
      
      //rollup is on therefore use abs
      const int g4ID(static_cast<int>(std::abs(trackIDE.trackID)));
      std::vector<const sim::IDE*> IDEs(btServ->TrackIdToSimIDEs_Ps(g4ID));
      std::vector<double> XYZ(btServ->SimIDEsToXYZ(IDEs));

      idToXYZMap[g4ID].push_back(XYZ); 
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
  pfpTree->Branch("pfpspX", &spX);
  pfpTree->Branch("pfpspY", &spY);
  pfpTree->Branch("pfpspZ", &spZ);
  pfpTree->Branch("pfptrueX", &pfptrueX);
  pfpTree->Branch("pfptrueY", &pfptrueY);
  pfpTree->Branch("pfptrueZ", &pfptrueZ);
  pfpTree->Branch("trueX", &trueX);
  pfpTree->Branch("trueY", &trueY);
  pfpTree->Branch("trueZ", &trueZ);
}

void sbnd::PFPValidation::clearPFPTree()
{
//  fEventID = -999;
  eventNumPFP = -999; 
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
  spX.clear();
  spY.clear();
  spZ.clear();
  pfptrueX.clear();
  pfptrueY.clear();
  pfptrueZ.clear();
  trueX.clear();
  trueY.clear();
  trueZ.clear();
}

DEFINE_ART_MODULE(sbnd::PFPValidation)
