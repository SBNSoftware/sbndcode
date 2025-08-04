////////////////////////////////////////////////////////////////////////
// Class:       POTprintout
// Plugin Type: analyzer (Unknown Unknown)
// File:        POTprintout_module.cc
//
// Generated on 1st July 2025 by Anna Beever
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileService.h"

// Additional LArSoft includes
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// ROOT includes
#include <TTree.h>

constexpr double double_default = -999.0;

namespace nuclearFragments {
  class POTprintout;
}


class nuclearFragments::POTprintout : public art::EDAnalyzer {
public:
  explicit POTprintout(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  POTprintout(POTprintout const&) = delete;
  POTprintout(POTprintout&&) = delete;
  POTprintout& operator=(POTprintout const&) = delete;
  POTprintout& operator=(POTprintout&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun (const art::SubRun& subrun);

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Create output TTree for subrun
  TTree *fSubRunTree;

  // Tree variables
  unsigned int fRun;
  unsigned int fSubRun;
  /*unsigned int fReco_nPFParticles;
  unsigned int fReco_nPrimaryChildren;
  bool fIsCCQE = 0;
  bool fHasCluster = 0;*/
  double POT;

  /*std::vector<double> fReco_childTrackLengths;
  std::vector<double> fReco_childTrackCompleteness;
  std::vector<double> fReco_childTrackPurity;
  std::vector<std::vector<double>> fReco_childTrackdEdx;
  std::vector<std::vector<double>> fReco_childTrackResRange;
  std::vector<int> fReco_truthMatchedTrackID;
  std::vector<int> fReco_truthMatchedPDG;
  std::vector<double> fReco_truthMatchedKE;
  std::vector<double> fReco_truthMatchedE;
  std::vector<double> fReco_truthMatchedP;
  std::vector<int> fMC_particlePDG;
  std::vector<double> fMC_particleMass;
  std::vector<double> fMC_particleE;
  std::vector<double> fMC_particleP;
  std::vector<double> fMC_particleKE;
  std::vector<bool> fMC_isReconstructed;
  std::vector<int> fMC_trackID;*/

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;
  std::string fNuGenLabel;
  std::string fLArGeantLabel;

  /*// Maps
  std::map<int,int> fTrackHitsMap;
  std::map<int,int> fTrackIDtoTruthPDGMap;
  std::map<int,double> fTrackIDtoTruthKEMap;
  std::map<int,double> fTrackIDtoTruthEMap;
  std::map<int,double> fTrackIDtoTruthPMap;*/

  // Functions 
  void ResetVariables();
  //float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  //float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  //detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


nuclearFragments::POTprintout::POTprintout(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fHitLabel(p.get<std::string>("HitLabel")),
  fNuGenLabel(p.get<std::string>("NuGenLabel")),
  fLArGeantLabel(p.get<std::string>("LArGeantLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void nuclearFragments::POTprintout::analyze(art::Event const& e)
{
  
}

void nuclearFragments::POTprintout::beginJob()
{
  // Get TFileService to create output TTree
  art::ServiceHandle<art::TFileService> tfs;
  fSubRunTree = tfs->make<TTree>("pottree", "Output TTree");

  // Add branches to TTRee
  fSubRunTree->Branch("run", &fRun);
  fSubRunTree->Branch("subRun", &fSubRun);
  /*fTree->Branch("isCCQE", &fIsCCQE);
  fTree->Branch("hasCluster", &fHasCluster);
  fTree->Branch("reco_nPFParticles", &fReco_nPFParticles);
  fTree->Branch("reco_nPrimaryChildren", &fReco_nPrimaryChildren);
  fTree->Branch("reco_childTrackLengths", &fReco_childTrackLengths);
  fTree->Branch("reco_childTrackCompleteness", &fReco_childTrackCompleteness);
  fTree->Branch("reco_childTrackPurity", &fReco_childTrackPurity);
  fTree->Branch("reco_childTrackdEdx", &fReco_childTrackdEdx);
  fTree->Branch("reco_childTrackResRange", &fReco_childTrackResRange);
  fTree->Branch("reco_truthMatchedTrackID", &fReco_truthMatchedTrackID);
  fTree->Branch("reco_truthMatchedPDG", &fReco_truthMatchedPDG);
  fTree->Branch("reco_truthMatchedE", &fReco_truthMatchedE);
  fTree->Branch("reco_truthMatchedP", &fReco_truthMatchedP);
  fTree->Branch("reco_truthMatchedKE", &fReco_truthMatchedKE);
  fTree->Branch("MC_particlePDG", &fMC_particlePDG);
  fTree->Branch("MC_particleMass", &fMC_particleMass);
  fTree->Branch("MC_particleE", &fMC_particleE);
  fTree->Branch("MC_particleP", &fMC_particleP);
  fTree->Branch("MC_particleKE", &fMC_particleKE);
  fTree->Branch("MC_isReconstructed", &fMC_isReconstructed);
  fTree->Branch("MC_trackID", &fMC_trackID);*/
  fSubRunTree->Branch("POT", &POT);
}

void nuclearFragments::POTprintout::beginSubRun(const art::SubRun& subrun)
{
  // Go through POTSummary objects
    fRun = subrun.run();
    fSubRun = subrun.subRun();
    art::Handle<sumdata::POTSummary> potHandle;
    subrun.getByLabel(fNuGenLabel, potHandle);
    const sumdata::POTSummary& potSum = (*potHandle);
    POT = potSum.totpot;
    std::cout << "POT: " << POT << std::endl;
    fSubRunTree->Fill();
    return;
}

void nuclearFragments::POTprintout::endJob()
{
  // Implementation of optional member function here.
}

void nuclearFragments::POTprintout::ResetVariables()
{
  fRun = -99999;
  fSubRun = -99999;
  // Set all trackCounters to zero for the current event
  //POT = -999;
  /*fIsCCQE = 0;
  fHasCluster = 0;
  fTrackHitsMap.clear();
  fTrackIDtoTruthPDGMap.clear();
  fTrackIDtoTruthKEMap.clear();
  fReco_nPFParticles = 0;
  fReco_nPrimaryChildren = 0;
  fReco_childTrackLengths.clear();
  fReco_childTrackCompleteness.clear();
  fReco_childTrackPurity.clear();
  fReco_childTrackdEdx.clear();
  fReco_childTrackResRange.clear();
  fReco_truthMatchedTrackID.clear();
  fReco_truthMatchedPDG.clear();
  fReco_truthMatchedE.clear();
  fReco_truthMatchedP.clear();
  fReco_truthMatchedKE.clear();
  fMC_particlePDG.clear();
  fMC_particleMass.clear();
  fMC_particleE.clear();
  fMC_particleP.clear();
  fMC_particleKE.clear();
  fMC_isReconstructed.clear();
  fMC_trackID.clear();*/
}

/*float nuclearFragments::POTprintout::Completeness(std::vector< art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (fTrackHitsMap[trackID] == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(fTrackHitsMap[trackID]);
}

float nuclearFragments::POTprintout::Purity(std::vector<art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (objectHits.size() == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}*/

DEFINE_ART_MODULE(nuclearFragments::POTprintout)
