////////////////////////////////////////////////////////////////////////
// Class:       PFPValidationCI
// Plugin Type: analyzer (art v3_02_06)
// File:        PFPValidationCI_module.cc
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
#include "sbndcode/Common/Modules/MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TTree.h"
#include <iostream>
#include <vector>
#include <numeric> 

namespace sbndcode {
class PFPValidationCI;
}

class sbndcode::PFPValidationCI : public art::EDAnalyzer {
  public:
  explicit PFPValidationCI(fhicl::ParameterSet const& pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPValidationCI(PFPValidationCI const&) = delete;
  PFPValidationCI(PFPValidationCI&&) = delete;
  PFPValidationCI& operator=(PFPValidationCI const&) = delete;
  PFPValidationCI& operator=(PFPValidationCI&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  void beginJob() override;

  void clearPFPTree(const std::string& PFParticleLabel);
  void clearTrueTree();
  void clearEventTree();

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

  template <class T>
  void initTree(TTree* Tree,
      const std::string& branchName,
      std::map<std::string, T>& Metric,
      const std::vector<std::string>& fPFPLabels);

  struct TruthMatch {

    TruthMatch()
        : mRecoId(0)
        , mTrueId(-999)
        , mRecoPdg(-999)
        , mNumHits(-999)
        , mHitComp(-999)
        , mHitPurity(-999)
        , mEnergyComp(-999)
        , mEnergyPurity(-999)
        , mHitSPRatio(-999)
        , mTrackScore(-999)
    {
    }

    TruthMatch(long unsigned int recoId, int trueId, int recoPdg, int numHits,
        float hitComp, float hitPurity, float energyComp, float energyPurity,
        float hitSPRatio, float trackScore)
        : mRecoId(recoId)
        , mTrueId(trueId)
        , mRecoPdg(recoPdg)
        , mNumHits(numHits)
        , mHitComp(hitComp)
        , mHitPurity(hitPurity)
        , mEnergyComp(energyComp)
        , mEnergyPurity(energyPurity)
        , mHitSPRatio(hitSPRatio)
        , mTrackScore(trackScore)
    {
    }

    long unsigned int mRecoId;
    int mTrueId, mRecoPdg, mNumHits;
    float mHitComp, mHitPurity, mEnergyComp, mEnergyPurity;
    float mHitSPRatio, mTrackScore;

    bool operator>(const TruthMatch& match) const { return mHitComp > match.mHitComp; }
  };

  private:
  std::vector<std::string> fPFPLabels; 
  std::vector<std::string> fHitLabels;
  std::string fSimChannelLabel;
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  // Declare member data here.
  unsigned int maxLabelLength;

  TTree* trueTree;
  TTree* pfpTree;
  TTree* eventTree;

  // Fill the tree once per particle
  int truePdg, motherPdg, numTrueHits;
  float trueEnergy, trueMomentum, trueDepositedEnergy, trueDepositedEnergyInHits;
  std::string trueProcess;
  bool isPrimary;
  std::map<std::string, int> recoPdg, recoPFPShowers, recoHits, recoPFPs, recoPFPTracks;
  std::map<std::string, float> hitPurity, energyPurity, hitComp, energyComp, hitSPRatio, recoTrackScore;

  std::string pfpModuleLabel;
  std::map<std::string, int> pfpPdg, pfpTruePdg, pfpNumHits;
  std::map<std::string, float> pfpTrueEnergy, pfpTrueMomentum, pfpTrueDepositedEnergy, pfpTrackScore;
  std::map<std::string, float> pfpHitPurity, pfpHitComp, pfpEnergyPurity, pfpEnergyComp, pfpHitSPRatio;

  std::map<std::string, int> eventNumPFP, eventNumPFPNu, eventNumPFPShower, eventNumPFPTrack;
};

sbndcode::PFPValidationCI::PFPValidationCI(fhicl::ParameterSet const& pset)
    : EDAnalyzer { pset }
    , fPFPLabels(pset.get<std::vector<std::string>>("PFPLabels"))
    , fHitLabels(pset.get<std::vector<std::string>>("HitLabels"))
    , fSimChannelLabel(pset.get<std::string>("SimChannelLabel"))
    , maxLabelLength(0)
{
  for (std::string const& fPFPLabel : fPFPLabels) {
    if (fPFPLabel.length() > maxLabelLength) {
      maxLabelLength = fPFPLabel.length();
    }
  }
}

void sbndcode::PFPValidationCI::beginJob()
{

  pfpTree = tfs->make<TTree>("pfpTree", "Tree with metrics for each pfp");
  trueTree = tfs->make<TTree>("trueTree", "Tree with metrics for each true particle");
  eventTree = tfs->make<TTree>("eventTree", "Tree with metrics for each event");

  trueTree->Branch("truePdg", &truePdg);
  trueTree->Branch("motherPdg", &motherPdg);
  trueTree->Branch("numTrueHits", &numTrueHits);
  trueTree->Branch("trueEnergy", &trueEnergy);
  trueTree->Branch("trueMomentum", &trueMomentum);
  trueTree->Branch("trueDepositedEnergy", &trueDepositedEnergy);
  trueTree->Branch("trueDepositedEnergyInHits", &trueDepositedEnergyInHits);
  trueTree->Branch("trueProcess", &trueProcess);
  trueTree->Branch("isPrimary", &isPrimary);

  initTree(trueTree, "recoPdg", recoPdg, fPFPLabels);
  initTree(trueTree, "recoTrackScore", recoTrackScore, fPFPLabels);
  initTree(trueTree, "recoHits", recoHits, fPFPLabels);
  initTree(trueTree, "recoPFPs", recoPFPs, fPFPLabels);
  initTree(trueTree, "recoPFPTracks", recoPFPTracks, fPFPLabels);
  initTree(trueTree, "recoPFPShowers", recoPFPShowers, fPFPLabels);
  initTree(trueTree, "hitPurity", hitPurity, fPFPLabels);
  initTree(trueTree, "hitComp", hitComp, fPFPLabels);
  initTree(trueTree, "energyPurity", energyPurity, fPFPLabels);
  initTree(trueTree, "energyComp", energyComp, fPFPLabels);
  initTree(trueTree, "hitSPRatio", hitSPRatio, fPFPLabels);

  initTree(pfpTree, "pfpPdg", pfpPdg, fPFPLabels);
  initTree(pfpTree, "pfpTruePdg", pfpTruePdg, fPFPLabels);
  initTree(pfpTree, "pfpNumHits", pfpNumHits, fPFPLabels);
  initTree(pfpTree, "pfpTrueEnergy", pfpTrueEnergy, fPFPLabels);
  initTree(pfpTree, "pfpTrueMomentum", pfpTrueMomentum, fPFPLabels);
  initTree(pfpTree, "pfpTrueDepositedEnergy", pfpTrueDepositedEnergy, fPFPLabels);
  initTree(pfpTree, "pfpTrackScore", pfpTrackScore, fPFPLabels);
  initTree(pfpTree, "pfpHitPurity", pfpHitPurity, fPFPLabels);
  initTree(pfpTree, "pfpHitComp", pfpHitComp, fPFPLabels);
  initTree(pfpTree, "pfpEnergyPurity", pfpEnergyPurity, fPFPLabels);
  initTree(pfpTree, "pfpEnergyComp", pfpEnergyComp, fPFPLabels);
  initTree(pfpTree, "pfpHitSPRatio", pfpHitSPRatio, fPFPLabels);

  initTree(eventTree, "numPFP", eventNumPFP, fPFPLabels);
  initTree(eventTree, "numPFPNu", eventNumPFPNu, fPFPLabels);
  initTree(eventTree, "numPFPShower", eventNumPFPShower, fPFPLabels);
  initTree(eventTree, "numPFPTrack", eventNumPFPTrack, fPFPLabels);
}

void sbndcode::PFPValidationCI::analyze(art::Event const& evt)
{

  std::cout << std::setprecision(2) << std::fixed;

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

  // Initialse some stuff
  auto const simChannelHandle(evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel));

  // Get all the hits
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > allHits;
  for(auto const& fHitModuleLabel: fHitLabels){

    if(evt.getByLabel(fHitModuleLabel,hitListHandle)) {
      art::fill_ptr_vector(allHits, hitListHandle);
    }
  }

  std::vector<art::Ptr<sim::SimChannel>> simChannels;
  art::fill_ptr_vector(simChannels, simChannelHandle);

  // Get map of true primary particle to number of reco hits / energy in reco hits
  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt));
  const std::map<int, int> truePrimaryHits(this->GetTruePrimaryHits(clockData, trueParticles, truePrimaries, allHits));
  const std::map<int, float> truePrimaryEnergies(this->GetTruePrimaryEnergies(trueParticles, truePrimaries, simChannels));
  const std::map<int, float> truePrimaryHitEnergies(this->GetTruePrimaryHitEnergies(clockData, trueParticles, truePrimaries, allHits));

  std::map<std::string, std::map<long unsigned int, art::Ptr<recob::PFParticle>>> pfpMap;
  std::map<std::string, std::map<long unsigned int, float>> pfpTrackScoreMap;

  std::map<std::string, std::vector<TruthMatch>> pfpTruthMatchMap;

  for (auto const fPFPLabel : fPFPLabels) {

    // Get all the PFPs
    auto const pfpHandle(evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
    std::vector<art::Ptr<recob::PFParticle>> pfps;
    art::fill_ptr_vector(pfps, pfpHandle);

    auto const spHandle(evt.getValidHandle<std::vector<recob::SpacePoint>>(fPFPLabel));
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    art::fill_ptr_vector(allSpacePoints, spHandle);

    auto const clusterHandle(evt.getValidHandle<std::vector<recob::Cluster>>(fPFPLabel));
    std::vector<art::Ptr<recob::Cluster>> clusters;
    art::fill_ptr_vector(clusters, clusterHandle);

    art::FindManyP<recob::SpacePoint> fmPFPSP(pfpHandle, evt, fPFPLabel);
    if (!fmPFPSP.isValid() || fmPFPSP.size() == 0) {
      std::cout << "PFP-SpacePoint assns not valid" << std::endl;
      continue;
    }

    art::FindManyP<recob::Cluster> fmPFPClus(pfpHandle, evt, fPFPLabel);
    if (!fmPFPClus.isValid() || fmPFPClus.size() == 0) {
      std::cout << "PFP-Cluster assns not valid" << std::endl;
      continue;
    }

    art::FindManyP<recob::Hit> fmClusHit(clusterHandle, evt, fPFPLabel);
    if (!fmClusHit.isValid() || fmClusHit.size() == 0) {
      std::cout << "Cluster-Hit assns not valid" << std::endl;
      continue;
    }

    // Create a map between PFParticles and their IDs
    art::FindManyP<larpandoraobj::PFParticleMetadata> fmpfpmd(pfps, evt, fPFPLabel);
    if (!fmpfpmd.isValid() || fmpfpmd.size() == 0) {
      std::cout << "PFP-MetaData assns not valid" << std::endl;
      return;
    }

    // std::map<long unsigned int, float > pfpScoreMap;
    for (auto const& pfp : pfps) {
      long unsigned int pfpId(pfp->Self());
      pfpMap[fPFPLabel][pfpId] = pfp;
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec(fmpfpmd.at(pfpId));
      for (auto const& pfpMeta : pfpMetaVec) {
        const larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMeta->GetPropertiesMap());
        auto const pfpTrackScoreIter(propertiesMap.find("TrackScore"));
        pfpTrackScoreMap[fPFPLabel][pfpId] = pfpTrackScoreIter == propertiesMap.end() ? -999 : pfpTrackScoreIter->second;
      }
    }

    clearEventTree();
    eventNumPFP[fPFPLabel] = pfps.size();

    for (auto const& pfp : pfps) {

      clearPFPTree(fPFPLabel);

      // Update event wide metrics
      if (pfp->PdgCode() == 11) {
        ++eventNumPFPShower[fPFPLabel];
      } else if (pfp->PdgCode() == 13) {
        ++eventNumPFPTrack[fPFPLabel];
      } else {
        ++eventNumPFPNu[fPFPLabel];
      }

      // Get the hits from the PFParticle
      std::vector<art::Ptr<recob::Hit>> pfpHits;
      const std::vector<art::Ptr<recob::Cluster>>& clusters(fmPFPClus.at(pfp.key()));
      for (const auto& cluster : clusters) {
        const std::vector<art::Ptr<recob::Hit>>& hits(fmClusHit.at(cluster.key()));
        pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
      }

      const std::vector<art::Ptr<recob::SpacePoint>>& sps(fmPFPSP.at(pfp.key()));

      if (pfpHits.empty())
        continue;

      pfpPdg[fPFPLabel] = pfp->PdgCode();
      pfpTrackScore[fPFPLabel] = pfpTrackScoreMap[fPFPLabel].at(pfp->Self());
      pfpNumHits[fPFPLabel] = pfpHits.size();
      pfpHitSPRatio[fPFPLabel] = (float)sps.size() / pfpHits.size();

      const std::pair<int, double> trueId(ShowerUtils::TrueParticleIDFromTrueChain(clockData, truePrimaries, pfpHits, 2));
      if (trueId.first == -99999) {
        pfpTree->Fill();
        continue;
      }

      const std::map<int, int> pfpTrueHitsMap(GetTruePrimaryHits(clockData, trueParticles, truePrimaries, pfpHits));
      const std::map<int, float> pfpTrueEnergyMap(GetTruePrimaryHitEnergies(clockData, trueParticles, truePrimaries, pfpHits));

      const simb::MCParticle* trueParticle(trueParticles.at(trueId.first));

      pfpTruePdg[fPFPLabel] = trueParticle->PdgCode();
      pfpTrueEnergy[fPFPLabel] = trueParticle->E();
      pfpTrueMomentum[fPFPLabel] = trueParticle->P();

      const int pfpHitsTrueHits(pfpTrueHitsMap.at(trueId.first));
      const float pfpHitsTrueEnergy(pfpTrueEnergyMap.at(trueId.first));
      const float pfpHitsTotalEnergy(GetTotalEnergyInHits(clockData, pfpHits));

      pfpHitPurity[fPFPLabel] = (float)pfpHitsTrueHits / pfpHits.size();
      pfpHitComp[fPFPLabel] = (float)pfpHitsTrueHits / truePrimaryHits.at(trueId.first);
      pfpEnergyPurity[fPFPLabel] = pfpHitsTrueEnergy / pfpHitsTotalEnergy;
      pfpEnergyComp[fPFPLabel] = pfpHitsTrueEnergy / truePrimaryHitEnergies.at(trueId.first);
      // Energy is in all planes, hence divide deposited energy by 3
      pfpTrueDepositedEnergy[fPFPLabel] = truePrimaryHitEnergies.at(trueId.first) / 3;

      pfpTree->Fill();

      TruthMatch match(pfp->Self(), trueId.first, pfpPdg[fPFPLabel], pfpNumHits[fPFPLabel],
          pfpHitComp[fPFPLabel], pfpHitPurity[fPFPLabel],
          pfpEnergyComp[fPFPLabel], pfpEnergyPurity[fPFPLabel],
          pfpHitSPRatio[fPFPLabel], pfpTrackScore[fPFPLabel]);

      pfpTruthMatchMap[fPFPLabel].push_back(match);
    }
    eventTree->Fill();
  }

  for (auto const& [trueId, trueDaughters] : truePrimaries) {

    clearTrueTree();

    numTrueHits = truePrimaryHits.at(trueId);
    if (!numTrueHits)
      continue;

    const simb::MCParticle* trueParticle(trueParticles.at(trueId));
    truePdg = trueParticle->PdgCode();
    trueEnergy = (trueParticle->P() * trueParticle->P()) / (2 * trueParticle->Mass());
    trueMomentum = trueParticle->P();
    trueDepositedEnergy = truePrimaryEnergies.at(trueId) / 3;
    trueDepositedEnergyInHits = truePrimaryHitEnergies.at(trueId) / 3;
    trueProcess = trueParticle->Process();
    isPrimary = trueProcess == "primary";

    auto const trueMother(trueParticles.find(trueParticle->Mother()));
    if (trueMother != trueParticles.end()) {
      motherPdg = trueMother->second->PdgCode();
    }

    if (trueProcess == "primary") {
      std::cout << "True Particle: " << truePdg << " with true deposited energy: " << trueDepositedEnergy << std::endl;
    }
    for (auto const fPFPLabel : fPFPLabels) {
      std::vector<TruthMatch> pfpMatches;
      unsigned int pfpTracks(0), pfpShowers(0);
      for (TruthMatch const& truthMatch : pfpTruthMatchMap[fPFPLabel]) {
        // Check the pfp is matched to the true particle we care about
        if (trueId == truthMatch.mTrueId) {
          pfpMatches.push_back(truthMatch);
          if (truthMatch.mRecoPdg == 11)
            ++pfpShowers;
          if (truthMatch.mRecoPdg == 13)
            ++pfpTracks;
        }
      }
      recoPFPs[fPFPLabel] = pfpTracks + pfpShowers;
      recoPFPTracks[fPFPLabel] = pfpTracks;
      recoPFPShowers[fPFPLabel] = pfpShowers;

      if (pfpMatches.size() > 0) {
        TruthMatch bestMatch;
        // Find the "Best Match" as the pfp with the highest completeness
        for (TruthMatch const& pfpMatch : pfpMatches) {
          if (pfpMatch > bestMatch) {
            bestMatch = pfpMatch;
          }
        }

        recoPdg[fPFPLabel] = bestMatch.mRecoPdg;
        recoHits[fPFPLabel] = bestMatch.mNumHits;
        hitComp[fPFPLabel] = bestMatch.mHitComp;
        hitPurity[fPFPLabel] = bestMatch.mHitPurity;
        energyComp[fPFPLabel] = bestMatch.mEnergyComp;
        energyPurity[fPFPLabel] = bestMatch.mEnergyPurity;
        hitSPRatio[fPFPLabel] = bestMatch.mHitSPRatio;
        recoTrackScore[fPFPLabel] = bestMatch.mTrackScore;
        if (trueProcess == "primary") {
          std::cout << "  " << std::setw(maxLabelLength) << fPFPLabel
                    << ": particles: " << pfpTracks + pfpShowers << "(" << pfpTracks << "+" << pfpShowers << ")"
                    << " and reco pdg: " << bestMatch.mRecoPdg
                    << " and Comp: " << bestMatch.mHitComp
                    << " and Purity: " << bestMatch.mHitPurity
                    << " and Track Score: " << bestMatch.mTrackScore
                    << std::endl;
        }
      }
    }
    trueTree->Fill();
  }
  std::cout << "\n"
            << std::endl;
}

std::map<int, int> sbndcode::PFPValidationCI::GetTruePrimaryHits(
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

std::map<int, float> sbndcode::PFPValidationCI::GetTruePrimaryEnergies(
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

std::map<int, float> sbndcode::PFPValidationCI::GetTruePrimaryHitEnergies(
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

float sbndcode::PFPValidationCI::GetTotalEnergyInHits(
    const detinfo::DetectorClocksData& clockData,
    const std::vector<art::Ptr<recob::Hit>>& hits) const
{

  TruthMatchUtils::IDToEDepositMap idToEDepMap;
  for (const art::Ptr<recob::Hit>& pHit : hits)
    TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepMap, clockData, pHit, true);

  return std::accumulate(idToEDepMap.cbegin(), idToEDepMap.cend(), 0.f,
      [](float sum, auto const& iter) { return sum + iter.second; });
}

template <class T>
void sbndcode::PFPValidationCI::initTree(TTree* Tree,
    const std::string& branchName,
    std::map<std::string, T>& Metric,
    const std::vector<std::string>& fPFPLabels)
{

  for (auto const& fPFPLabel : fPFPLabels) {
    const std::string branchString(branchName + "_" + fPFPLabel);
    Tree->Branch(branchString.c_str(), &Metric[fPFPLabel], 32000, 0);
  }
}

void sbndcode::PFPValidationCI::clearPFPTree(const std::string& PFParticleLabel)
{
  pfpPdg[PFParticleLabel] = -999;
  pfpTruePdg[PFParticleLabel] = -999;
  pfpTrueEnergy[PFParticleLabel] = -999;
  pfpTrueMomentum[PFParticleLabel] = -999;
  pfpTrueDepositedEnergy[PFParticleLabel] = -999;
  pfpTrackScore[PFParticleLabel] = -999;
  pfpNumHits[PFParticleLabel] = -999;
  pfpHitPurity[PFParticleLabel] = -999;
  pfpHitComp[PFParticleLabel] = -999;
  pfpEnergyPurity[PFParticleLabel] = -999;
  pfpEnergyComp[PFParticleLabel] = -999;
  pfpHitSPRatio[PFParticleLabel] = -999;
}

void sbndcode::PFPValidationCI::clearTrueTree()
{
  truePdg = -999;
  motherPdg = -999;
  numTrueHits = -999;
  trueEnergy = -999;
  trueMomentum = -999;
  trueDepositedEnergy = -999;
  trueDepositedEnergyInHits = -999;
  trueProcess = "";
  isPrimary = false;
  for (auto const& fPFPLabel : fPFPLabels) {
    recoPFPs[fPFPLabel] = -999;
    recoPFPTracks[fPFPLabel] = -999;
    recoPFPShowers[fPFPLabel] = -999;
    recoPdg[fPFPLabel] = -999;
    recoTrackScore[fPFPLabel] = -999;
    recoHits[fPFPLabel] = -999;
    hitPurity[fPFPLabel] = -999;
    energyPurity[fPFPLabel] = -999;
    hitComp[fPFPLabel] = -999;
    energyComp[fPFPLabel] = -999;
    hitSPRatio[fPFPLabel] = -999;
  }
}

void sbndcode::PFPValidationCI::clearEventTree()
{
  for (auto const& fPFPLabel : fPFPLabels) {
    eventNumPFP[fPFPLabel] = 0;
    eventNumPFPNu[fPFPLabel] = 0;
    eventNumPFPTrack[fPFPLabel] = 0;
    eventNumPFPShower[fPFPLabel] = 0;
  }
}

DEFINE_ART_MODULE(sbndcode::PFPValidationCI)
