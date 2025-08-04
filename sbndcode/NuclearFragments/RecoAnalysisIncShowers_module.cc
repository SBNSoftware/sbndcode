////////////////////////////////////////////////////////////////////////
// Class:       RecoAnalysisIncShowers
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoAnalysisIncShowers_module.cc
//
// Generated on Friday 11th July 2025 by Anna Beever
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FileBlock.h"
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
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include <TTree.h>

constexpr double double_default = -999.0;

namespace nuclearFragments {
  class RecoAnalysisIncShowers;
}


class nuclearFragments::RecoAnalysisIncShowers : public art::EDAnalyzer {
public:
  explicit RecoAnalysisIncShowers(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoAnalysisIncShowers(RecoAnalysisIncShowers const&) = delete;
  RecoAnalysisIncShowers(RecoAnalysisIncShowers&&) = delete;
  RecoAnalysisIncShowers& operator=(RecoAnalysisIncShowers const&) = delete;
  RecoAnalysisIncShowers& operator=(RecoAnalysisIncShowers&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void respondToOpenInputFile(const art::FileBlock& fb);

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Input file
  std::string fInputFile;

  // Create output TTree
  TTree *fEventTree;
  //TTree *fHitTree;

  // Tree variables
  unsigned int fEventID;
  unsigned int fReco_nPFParticles;
  unsigned int fReco_nPrimaryChildren;
  bool fIsCCQE = 0;
  bool fHasCluster = 0;
  double fReco_nuVertexX = double_default;
  double fReco_nuVertexY = double_default;
  double fReco_nuVertexZ = double_default;

  std::vector<double> fReco_childTrackLengths;
  std::vector<double> fReco_childTrackCompleteness;
  std::vector<double> fReco_childTrackPurity;
  std::vector<double> fReco_childTrackStartX;
  std::vector<double> fReco_childTrackStartY;
  std::vector<double> fReco_childTrackStartZ;
  std::vector<double> fReco_childTrackEndX;
  std::vector<double> fReco_childTrackEndY;
  std::vector<double> fReco_childTrackEndZ;
  std::vector<std::vector<double>> fReco_childTrackdEdx;
  std::vector<std::vector<double>> fReco_childTrackResRange;
  std::vector<double> fReco_childShowerLengths;
  std::vector<double> fReco_childShowerCompleteness;
  std::vector<double> fReco_childShowerPurity;
  std::vector<double> fReco_childShowerStartX;
  std::vector<double> fReco_childShowerStartY;
  std::vector<double> fReco_childShowerStartZ;
  std::vector<std::vector<double>> fReco_childShowerdEdx;
  //std::vector<std::vector<double>> fReco_childShowerResRange;
  std::vector<int> fReco_truthMatchedShowerID;
  std::vector<int> fReco_truthMatchedTrackID;
  std::vector<int> fReco_truthMatchedPFPg4ID;
  std::vector<int> fReco_truthMatchedTrackPDG;
  std::vector<double> fReco_truthMatchedTrackKE;
  std::vector<double> fReco_truthMatchedTrackE;
  std::vector<double> fReco_truthMatchedTrackP;
  std::vector<int> fMC_particlePDG;
  std::vector<double> fMC_particleMass;
  std::vector<double> fMC_particleE;
  std::vector<double> fMC_particleP;
  std::vector<double> fMC_particleKE;
  std::vector<double> fMC_particleVx;
  std::vector<double> fMC_particleVy;
  std::vector<double> fMC_particleVz;
  std::vector<double> fMC_particleEndX;
  std::vector<double> fMC_particleEndY;
  std::vector<double> fMC_particleEndZ;
  std::vector<bool> fMC_isReconstructedAsTrack;
  std::vector<bool> fMC_isReconstructedAsShower;
  std::vector<bool> fMC_isReconstructedAsPFP;
  //std::vector<bool> fMC_isShower;
  //std::vector<bool> fMC_isTrack;
  std::vector<int> fMC_trackID;

  std::vector<int> fReco_hitNumber;
  std::vector<int> fReco_hitG4ID;
  std::vector<bool> fReco_hitHasTrackObject;
  std::vector<bool> fReco_hitHasShowerObject;
  std::vector<bool> fReco_hitHasPFPObject;
  std::vector<int> fReco_hitTrackObjectG4ID;
  std::vector<int> fReco_hitShowerObjectG4ID;
  std::vector<int> fReco_hitPFPObjectG4ID;
  std::vector<std::vector<int>> fReco_hitG4IDs;
  std::vector<std::vector<double>> fReco_hitEdeps;
  

  // New
  std::vector<int> fReco_hitPlane;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  //std::string fShowerLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;
  std::string fVertexLabel;
  std::string fClusterLabel;
  std::string fRecoHitLabel;
  std::string fLArGeantLabel;
  std::string fNuGenLabel;

  // Maps
  std::map<int,int> fTrackHitsMap;
  std::map<int,int> fTrackIDtoTruthPDGMap;
  std::map<int,double> fTrackIDtoTruthKEMap;
  std::map<int,double> fTrackIDtoTruthEMap;
  std::map<int,double> fTrackIDtoTruthPMap;

  // Functions 
  void ResetVariables();
  float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  void AnalyseTrack(const art::Ptr<recob::Track> &track, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, std::vector<art::Ptr<recob::Hit>> hitVector);
  void AnalyseShower(const art::Ptr<recob::Shower> &shower, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc, std::vector<art::Ptr<recob::Hit>> hitVector);
  void AnalyseTruNu(const art::Ptr<simb::MCTruth> &truNu, art::FindManyP<simb::MCParticle> truNuParticleAssoc);
  void AnalyseMCTruth(art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle, art::FindManyP<simb::MCParticle> truNuParticleAssoc);
  void TruthMatch();
  void AnalyseReco(art::ValidHandle<std::vector<recob::Hit>> hitHandle, art::ValidHandle<std::vector<recob::Slice>> sliceHandle, art::ValidHandle<std::vector<recob::Track>> trackHandle, art::FindManyP<recob::PFParticle> slicePFPAssoc, art::FindManyP<recob::Track> pfpTrackAssoc, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, art::FindOneP<recob::Vertex> pfpVertexAssoc, art::ValidHandle<std::vector<recob::Shower>> showerHandle, art::FindManyP<recob::Shower> pfpShowerAssoc, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc, art::FindManyP<recob::Cluster> pfpClusterAssoc, art::FindManyP<recob::Hit> clusterHitAssoc);
  void shout(std::string message);

  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


nuclearFragments::RecoAnalysisIncShowers::RecoAnalysisIncShowers(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  //fShowerLabel(p.get<std::string>("ShowerLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fHitLabel(p.get<std::string>("HitLabel")),
  fVertexLabel(p.get<std::string>("VertexLabel")),
  fClusterLabel(p.get<std::string>("ClusterLabel")),
  fRecoHitLabel(p.get<std::string>("RecoHitLabel")),
  fLArGeantLabel(p.get<std::string>("LArGeantLabel")),
  fNuGenLabel(p.get<std::string>("NuGenLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void nuclearFragments::RecoAnalysisIncShowers::analyze(art::Event const& e)
{

  ResetVariables();

  // Set the event ID
  fEventID = e.id().event();

  //shout("====================== new event ===============================");

  // ========================================================================
  // Handles
  // ========================================================================

  // Get hit labelled objects for event
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);

  // Get nugen labelled objects
  art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fNuGenLabel);

  // Get track labelled objects
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);

  // Get shower labelled objects
  art::ValidHandle<std::vector<recob::Shower>> showerHandle = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);

  // Get PFP labelled objects
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

  // Get slice labelled objects
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);

  // Get cluster labelled objects
  art::ValidHandle<std::vector<recob::Cluster>> clusterHandle = e.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel);

  // ========================================================================
  // Associations
  // ========================================================================

  // Get associations between slices and PFParticles
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

  // Get tracks associated with the PFPs
  art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);

  // Associate tracks with hits
  art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackLabel);

  // Associate calorimetry info with tracks
  art::FindManyP<anab::Calorimetry> trackCaloAssoc(trackHandle, e, fCalorimetryLabel);

  // Get showers associated with the PFPs
  art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, e, fShowerLabel);

  // Associate showers with hits
  art::FindManyP<recob::Hit> showerHitAssoc(showerHandle,e,fShowerLabel);

  // Associate calorimetry info with showers
  art::FindManyP<anab::Calorimetry> showerCaloAssoc(showerHandle, e, fCalorimetryLabel);

  // Get associations between MC Particles and true neutrinos
  art::FindManyP<simb::MCParticle> truNuParticleAssoc(truNuHandle, e, fLArGeantLabel);

  // Get associations between PFPs and vertices
  art::FindOneP<recob::Vertex>pfpVertexAssoc(pfpHandle, e, fVertexLabel);

  // Get associations between PFPs and their clusters
  art::FindManyP<recob::Cluster>pfpClusterAssoc(pfpHandle, e, fClusterLabel);

  // Get associations between clusters and hits
  art::FindManyP<recob::Hit>clusterHitAssoc(clusterHandle, e, fRecoHitLabel);

  // ========================================================================
  // Main analysis
  // ========================================================================


  AnalyseReco(hitHandle, sliceHandle, trackHandle, slicePFPAssoc, pfpTrackAssoc, trackHitAssoc, trackCaloAssoc, pfpVertexAssoc, showerHandle, pfpShowerAssoc, showerHitAssoc, showerCaloAssoc, pfpClusterAssoc, clusterHitAssoc);

  AnalyseMCTruth(truNuHandle, truNuParticleAssoc);

  TruthMatch();

  // Fill event tree
  fEventTree->Fill();

}

void nuclearFragments::RecoAnalysisIncShowers::respondToOpenInputFile(const art::FileBlock& fb)
{
  fInputFile = fb.fileName();
  shout("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
  shout(fInputFile);
  shout("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
}

void nuclearFragments::RecoAnalysisIncShowers::beginJob()
{
  // Get TFileService to create output TTree
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("eventTree", "Output TTree");

  // Add branches to TTree
  fEventTree->Branch("eventID", &fEventID);
  fEventTree->Branch("inputFile", &fInputFile);
  fEventTree->Branch("isCCQE", &fIsCCQE);
  fEventTree->Branch("hasCluster", &fHasCluster);
  fEventTree->Branch("reco_nPFParticles", &fReco_nPFParticles);
  fEventTree->Branch("reco_nPrimaryChildren", &fReco_nPrimaryChildren);
  fEventTree->Branch("reco_childTrackLengths", &fReco_childTrackLengths);
  fEventTree->Branch("reco_childTrackCompleteness", &fReco_childTrackCompleteness);
  fEventTree->Branch("reco_childTrackPurity", &fReco_childTrackPurity);
  fEventTree->Branch("reco_childTrackdEdx", &fReco_childTrackdEdx);
  fEventTree->Branch("reco_childTrackResRange", &fReco_childTrackResRange);
  fEventTree->Branch("reco_childTrackStartX", &fReco_childTrackStartX);
  fEventTree->Branch("reco_childTrackStartY", &fReco_childTrackStartY);
  fEventTree->Branch("reco_childTrackStartZ", &fReco_childTrackStartZ);
  fEventTree->Branch("reco_childTrackEndX", &fReco_childTrackEndX);
  fEventTree->Branch("reco_childTrackEndY", &fReco_childTrackEndY);
  fEventTree->Branch("reco_childTrackEndZ", &fReco_childTrackEndZ);
  fEventTree->Branch("reco_childShowerLengths", &fReco_childShowerLengths);
  fEventTree->Branch("reco_childShowerCompleteness", &fReco_childShowerCompleteness);
  fEventTree->Branch("reco_childShowerPurity", &fReco_childShowerPurity);
  fEventTree->Branch("reco_childShowerdEdx", &fReco_childShowerdEdx);
  //fEventTree->Branch("reco_childShowerResRange", &fReco_childShowerResRange);
  fEventTree->Branch("reco_childShowerStartX", &fReco_childShowerStartX);
  fEventTree->Branch("reco_childShowerStartY", &fReco_childShowerStartY);
  fEventTree->Branch("reco_childShowerStartZ", &fReco_childShowerStartZ);
  fEventTree->Branch("reco_truthMatchedTrackID", &fReco_truthMatchedTrackID);
  fEventTree->Branch("reco_truthMatchedShowerID", &fReco_truthMatchedShowerID);
  fEventTree->Branch("reco_truthMatchedPFPg4ID", &fReco_truthMatchedPFPg4ID);
  fEventTree->Branch("reco_truthMatchedTrackPDG", &fReco_truthMatchedTrackPDG);
  fEventTree->Branch("reco_truthMatchedTrackE", &fReco_truthMatchedTrackE);
  fEventTree->Branch("reco_truthMatchedTrackP", &fReco_truthMatchedTrackP);
  fEventTree->Branch("reco_truthMatchedTrackKE", &fReco_truthMatchedTrackKE);
  fEventTree->Branch("reco_nuVertexX", &fReco_nuVertexX);
  fEventTree->Branch("reco_nuVertexY", &fReco_nuVertexY);
  fEventTree->Branch("reco_nuVertexZ", &fReco_nuVertexZ);
  fEventTree->Branch("reco_hitNumber", &fReco_hitNumber);
  fEventTree->Branch("reco_hitG4ID", &fReco_hitG4ID);
  fEventTree->Branch("reco_hitHasTrackObject", &fReco_hitHasTrackObject);
  fEventTree->Branch("reco_hitHasShowerObject", &fReco_hitHasShowerObject);
  fEventTree->Branch("reco_hitHasPFPObject", &fReco_hitHasPFPObject);
  fEventTree->Branch("reco_hitTrackObjectG4ID", &fReco_hitTrackObjectG4ID);
  fEventTree->Branch("reco_hitShowerObjectG4ID", &fReco_hitShowerObjectG4ID);
  fEventTree->Branch("reco_hitPFPObjectG4ID", &fReco_hitPFPObjectG4ID);
  fEventTree->Branch("reco_hitG4IDs", &fReco_hitG4IDs);
  fEventTree->Branch("reco_hitEdeps", &fReco_hitEdeps);
  fEventTree->Branch("MC_particlePDG", &fMC_particlePDG);
  fEventTree->Branch("MC_particleMass", &fMC_particleMass);
  fEventTree->Branch("MC_particleE", &fMC_particleE);
  fEventTree->Branch("MC_particleP", &fMC_particleP);
  fEventTree->Branch("MC_particleKE", &fMC_particleKE);
  fEventTree->Branch("MC_particleVx", &fMC_particleVx);
  fEventTree->Branch("MC_particleVy", &fMC_particleVy);
  fEventTree->Branch("MC_particleVz", &fMC_particleVz);
  fEventTree->Branch("MC_particleEndX", &fMC_particleEndX);
  fEventTree->Branch("MC_particleEndY", &fMC_particleEndY);
  fEventTree->Branch("MC_particleEndZ", &fMC_particleEndZ);
  fEventTree->Branch("MC_isReconstructedAsTrack", &fMC_isReconstructedAsTrack);
  fEventTree->Branch("MC_isReconstructedAsShower", &fMC_isReconstructedAsShower);
  fEventTree->Branch("MC_isReconstructedAsPFP", &fMC_isReconstructedAsPFP);
  //fEventTree->Branch("MC_isShower", &fMC_isShower);
  //fEventTree->Branch("MC_isTrack", &fMC_isTrack);
  fEventTree->Branch("MC_trackID", &fMC_trackID);
}

void nuclearFragments::RecoAnalysisIncShowers::endJob()
{
  // Implementation of optional member function here.
}

void nuclearFragments::RecoAnalysisIncShowers::ResetVariables()
{
  // Set all trackCounters to zero for the current event
  //fInputFile = "empty";
  fEventID = 0;
  fIsCCQE = 0;
  fHasCluster = 0;
  fReco_nuVertexX = double_default;
  fReco_nuVertexY = double_default;
  fReco_nuVertexZ = double_default;
  fTrackHitsMap.clear();
  fTrackIDtoTruthPDGMap.clear();
  fTrackIDtoTruthKEMap.clear();
  fTrackIDtoTruthEMap.clear();
  fTrackIDtoTruthPMap.clear();
  fReco_nPFParticles = 0;
  fReco_nPrimaryChildren = 0;
  fReco_childTrackLengths.clear();
  fReco_childTrackCompleteness.clear();
  fReco_childTrackPurity.clear();
  fReco_childTrackdEdx.clear();
  fReco_childTrackResRange.clear();
  fReco_childTrackStartX.clear();
  fReco_childTrackStartY.clear();
  fReco_childTrackStartZ.clear();
  fReco_childTrackEndX.clear();
  fReco_childTrackEndY.clear();
  fReco_childTrackEndZ.clear();
  fReco_childShowerLengths.clear();
  fReco_childShowerCompleteness.clear();
  fReco_childShowerPurity.clear();
  fReco_childShowerdEdx.clear();
  //fReco_childShowerResRange.clear();
  fReco_childShowerStartX.clear();
  fReco_childShowerStartY.clear();
  fReco_childShowerStartZ.clear();
  fReco_truthMatchedTrackID.clear();
  fReco_truthMatchedShowerID.clear();
  fReco_truthMatchedPFPg4ID.clear();
  fReco_truthMatchedTrackPDG.clear();
  fReco_truthMatchedTrackE.clear();
  fReco_truthMatchedTrackP.clear();
  fReco_truthMatchedTrackKE.clear();
  fReco_hitNumber.clear();
  fReco_hitG4ID.clear();
  fReco_hitHasTrackObject.clear();
  fReco_hitHasShowerObject.clear();
  fReco_hitHasPFPObject.clear();
  fReco_hitTrackObjectG4ID.clear();
  fReco_hitShowerObjectG4ID.clear();
  fReco_hitPFPObjectG4ID.clear();
  fReco_hitG4IDs.clear();
  fReco_hitEdeps.clear();
  fMC_particlePDG.clear();
  fMC_particleMass.clear();
  fMC_particleE.clear();
  fMC_particleP.clear();
  fMC_particleKE.clear();
  fMC_particleVx.clear();
  fMC_particleVy.clear();
  fMC_particleVz.clear();
  fMC_particleEndX.clear();
  fMC_particleEndY.clear();
  fMC_particleEndZ.clear();
  fMC_isReconstructedAsTrack.clear();
  fMC_isReconstructedAsShower.clear();
  fMC_isReconstructedAsPFP.clear();
  //fMC_isShower.clear();
  //fMC_isTrack.clear();
  fMC_trackID.clear();
}

float nuclearFragments::RecoAnalysisIncShowers::Completeness(std::vector< art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (fTrackHitsMap[trackID] == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(fTrackHitsMap[trackID]);
}

float nuclearFragments::RecoAnalysisIncShowers::Purity(std::vector<art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (objectHits.size() == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}

void nuclearFragments::RecoAnalysisIncShowers::AnalyseTrack(const art::Ptr<recob::Track> &track, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, std::vector<art::Ptr<recob::Hit>> hitVector)
{

  // Fill in track variables
  fReco_childTrackStartX.push_back((track->Start()).X());
  fReco_childTrackStartY.push_back((track->Start()).Y());
  fReco_childTrackStartZ.push_back((track->Start()).Z());

  fReco_childTrackEndX.push_back((track->End()).X());
  fReco_childTrackEndY.push_back((track->End()).Y());
  fReco_childTrackEndZ.push_back((track->End()).Z());

  // Get hits associated with this track
  std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());

  // Get just plane 2 hits
  std::vector<art::Ptr<recob::Hit>> trackHitsPlane2;

  for(unsigned int j = 0; j < trackHits.size(); j++)
  {
    if(trackHits[j]->WireID().Plane == 2)
    {
      trackHitsPlane2.push_back(trackHits[j]);
    }
  }

  // Fill track variables
  //shout(("trackhits has size " + std::to_string(trackHits.size())).c_str());
  //shout(("trackhits plane 2 has size " + std::to_string(trackHitsPlane2.size())).c_str());
  int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
  //shout(("this track has g4ID of " + std::to_string(trackID)).c_str());
  
  //int trackIDPlane2 = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHitsPlane2,true);
  //shout(("this track plane 2 has g4ID of " + std::to_string(trackIDPlane2)).c_str());
  fReco_truthMatchedTrackID.push_back(trackID);
  fReco_childTrackCompleteness.push_back(Completeness(trackHits,trackID));
  fReco_childTrackPurity.push_back(Purity(trackHits,trackID));
  fReco_childTrackLengths.push_back(track->Length());

  // Compare hits
  for(std::size_t i = 0; i < hitVector.size(); i++){
    int isInTrack = std::count(trackHits.begin(),trackHits.end(),hitVector[i]);
    if(isInTrack > 1){
      std::cout << "!!!!!!!!!!!!!!!weird!!!!!!!!!!!!!!!!!!!1" << std::endl;
    }
    if(isInTrack==1){
      fReco_hitHasTrackObject[i] = true;
      fReco_hitTrackObjectG4ID[i] = trackID;
    }
  }

  // Get the calorimetry objects associated with this track
  std::vector<art::Ptr<anab::Calorimetry>> calos = trackCaloAssoc.at(track.key());

  // Fill the track's calorimetry info
  for(const art::Ptr<anab::Calorimetry> &calo : calos)
  {
    const int plane = calo->PlaneID().Plane;
    if(plane != 2)
    {
      continue;
    }
    std::vector<double> dEdx (calo->dEdx().begin(), calo->dEdx().end());
    std::vector<double> resRange (calo->ResidualRange().begin(), calo->ResidualRange().end());
    fReco_childTrackdEdx.push_back(dEdx);
    fReco_childTrackResRange.push_back(resRange);
  }

}


void nuclearFragments::RecoAnalysisIncShowers::AnalyseShower(const art::Ptr<recob::Shower> &shower, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc, std::vector<art::Ptr<recob::Hit>> hitVector)
{

  // Fill in shower variables
  fReco_childShowerStartX.push_back((shower->ShowerStart()).X());
  fReco_childShowerStartY.push_back((shower->ShowerStart()).Y());
  fReco_childShowerStartZ.push_back((shower->ShowerStart()).Z());

  // Get hits associated with this shower
  std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssoc.at(shower.key());

  // Fill shower variables
  int showerID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,showerHits,true);
  fReco_truthMatchedShowerID.push_back(showerID);
  fReco_childShowerCompleteness.push_back(Completeness(showerHits,showerID));
  fReco_childShowerPurity.push_back(Purity(showerHits,showerID));
  fReco_childShowerLengths.push_back(shower->Length());

  // Compare hits
  for(std::size_t i = 0; i < hitVector.size(); i++){
    int isInShower = std::count(showerHits.begin(),showerHits.end(),hitVector[i]);
    if(isInShower > 1){
      std::cout << "!!!!!!!!!!!!!!!weird!!!!!!!!!!!!!!!!!!!1" << std::endl;
    }
    if(isInShower==1){
      fReco_hitHasShowerObject[i] = true;
      fReco_hitShowerObjectG4ID[i] = showerID;
    }
  }

  
  fReco_childShowerdEdx.push_back(shower->dEdx());
  //fReco_childShowerResRange.push_back(resRange);
  

}

void nuclearFragments::RecoAnalysisIncShowers::AnalyseTruNu(const art::Ptr<simb::MCTruth> &truNu, art::FindManyP<simb::MCParticle> truNuParticleAssoc)
{

  // Get neutrino
  const simb::MCNeutrino nu = truNu->GetNeutrino();

  // Check interaction type (only works for GENIE events for now)
  if(nu.InteractionType() == simb::kCCQE){
    fIsCCQE = 1;
    std::cout << "CCQE" << std::endl;
  }

  // Get particles associated with true neutrino
  std::vector<art::Ptr<simb::MCParticle>> particles = truNuParticleAssoc.at(truNu.key());

  // Loop over the neutrino's particles 
  for(const art::Ptr<simb::MCParticle> &particle : particles){

    // Only get primary particles - the ones leaving the nucleus (?)
    if(particle->Process() != "primary" || particle->StatusCode() != 1){
      continue;
    }

    //shout(("particle of G4ID " + std::to_string(particle->TrackId())).c_str());
    //shout(("particle of PDG code " + std::to_string(particle->PdgCode())).c_str());

    // Fill MC particle truth info

    fMC_trackID.push_back(particle->TrackId());
    fMC_particlePDG.push_back(particle->PdgCode());
    fMC_particleMass.push_back(particle->Mass());
    fMC_particleE.push_back(particle->E());
    fMC_particleP.push_back(particle->P());
    double KE = particle->E() - particle->Mass();
    fMC_particleKE.push_back(KE);
    fMC_particleVx.push_back(particle->Vx());
    fMC_particleVy.push_back(particle->Vy());
    fMC_particleVz.push_back(particle->Vz());
    fMC_particleEndX.push_back(particle->EndX());
    fMC_particleEndY.push_back(particle->EndY());
    fMC_particleEndZ.push_back(particle->EndZ());
    fTrackIDtoTruthPDGMap.insert({particle->TrackId(), particle->PdgCode()});
    fTrackIDtoTruthEMap.insert({particle->TrackId(), particle->E()});
    fTrackIDtoTruthPMap.insert({particle->TrackId(), particle->P()});
    fTrackIDtoTruthKEMap.insert({particle->TrackId(), KE});


    int checkIfReconstructedAsTrack = std::count(fReco_truthMatchedTrackID.begin(), fReco_truthMatchedTrackID.end(), particle->TrackId());
    fMC_isReconstructedAsTrack.push_back(checkIfReconstructedAsTrack!=0);

    int checkIfReconstructedAsShower = std::count(fReco_truthMatchedShowerID.begin(), fReco_truthMatchedShowerID.end(), particle->TrackId());
    fMC_isReconstructedAsShower.push_back(checkIfReconstructedAsShower!=0);

    int checkIfReconstructedAsPFP = std::count(fReco_truthMatchedPFPg4ID.begin(), fReco_truthMatchedPFPg4ID.end(), particle->TrackId());
    fMC_isReconstructedAsPFP.push_back(checkIfReconstructedAsPFP!=0);

    // check for clusters
    if(particle->PdgCode() == 1000010020 || particle->PdgCode() == 1000010030 || particle->PdgCode() == 1000020030 || particle->PdgCode() == 1000020040){
      fHasCluster = 1;
    }

  }
  
}

void nuclearFragments::RecoAnalysisIncShowers::AnalyseMCTruth(art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle, art::FindManyP<simb::MCParticle> truNuParticleAssoc)
{
  // Get true neutrinos
  std::vector<art::Ptr<simb::MCTruth>> truNuVector;

  if(truNuHandle.isValid())
  {
    art::fill_ptr_vector(truNuVector,truNuHandle);
  }

  // Loop over the true neutrinos 

  for(const art::Ptr<simb::MCTruth> &truNu : truNuVector){
    if(truNu.isNull()){
      continue;
    }

    AnalyseTruNu(truNu, truNuParticleAssoc);
  }

}

void nuclearFragments::RecoAnalysisIncShowers::TruthMatch()
{
  for(std::size_t i=0; i<fReco_truthMatchedTrackID.size(); i++){
    fReco_truthMatchedTrackPDG.push_back(fTrackIDtoTruthPDGMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedTrackE.push_back(fTrackIDtoTruthEMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedTrackP.push_back(fTrackIDtoTruthPMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedTrackKE.push_back(fTrackIDtoTruthKEMap[fReco_truthMatchedTrackID[i]]);
  }
}

void nuclearFragments::RecoAnalysisIncShowers::AnalyseReco(art::ValidHandle<std::vector<recob::Hit>> hitHandle, art::ValidHandle<std::vector<recob::Slice>> sliceHandle, art::ValidHandle<std::vector<recob::Track>> trackHandle, art::FindManyP<recob::PFParticle> slicePFPAssoc, art::FindManyP<recob::Track> pfpTrackAssoc, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, art::FindOneP<recob::Vertex> pfpVertexAssoc, art::ValidHandle<std::vector<recob::Shower>> showerHandle, art::FindManyP<recob::Shower> pfpShowerAssoc, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc, art::FindManyP<recob::Cluster> pfpClusterAssoc, art::FindManyP<recob::Hit> clusterHitAssoc)
{
  // Fill hits
  std::vector<art::Ptr<recob::Hit>> hitVector;

  if(hitHandle.isValid())
  {
    art::fill_ptr_vector(hitVector,hitHandle);
  }

  // Track-Hit map setup

  int hitCount = 0;
  for(const art::Ptr<recob::Hit> &hit : hitVector){
    fTrackHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
    fReco_hitNumber.push_back(hitCount);
    fReco_hitHasTrackObject.push_back(false);
    fReco_hitHasShowerObject.push_back(false);
    fReco_hitHasPFPObject.push_back(false);
    fReco_hitTrackObjectG4ID.push_back(-999);
    fReco_hitShowerObjectG4ID.push_back(-999);
    fReco_hitPFPObjectG4ID.push_back(-999);
    fReco_hitG4ID.push_back(TruthMatchUtils::TrueParticleID(clockData,hit,true));
    hitCount++;
    TruthMatchUtils::IDToEDepositMap idToEDepMap;
    std::vector<int> G4IDsForThisHit;
    std::vector<double> eDepsForThisHit;
    TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepMap,clockData,hit,true);
    for(const auto & [G4ID,eDep] : idToEDepMap)
    {
      G4IDsForThisHit.push_back(G4ID);
      eDepsForThisHit.push_back(eDep);
    }
    fReco_hitG4IDs.push_back(G4IDsForThisHit);
    fReco_hitEdeps.push_back(eDepsForThisHit);
  }


  // Fill slices vector
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if(sliceHandle.isValid())
  {
    art::fill_ptr_vector(sliceVector,sliceHandle);
  }

  //Filling neutrino hierachy variables
  int nuID = -1;
  int nuSliceKey = -1;

  // Loop over slices

  for(const art::Ptr<recob::Slice> &slice : sliceVector)
  {
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

    // Loop over PFPs

    for(const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs)
    {

      // Find primary neutrino
      const bool isPrimary(slicePFP->IsPrimary());
      const bool isNeutrino((std::abs(slicePFP->PdgCode() == 12) || (std::abs(slicePFP->PdgCode() == 14))));

      if(!(isPrimary && isNeutrino))
      {
        continue;
      }

      //shout("new PFP");

      const art::Ptr<recob::Vertex> vertex = pfpVertexAssoc.at(slicePFP.key());
      geo::Point_t vertexPos = vertex.isNonnull() ? vertex->position() : geo::Point_t(double_default, double_default, double_default);
      fReco_nuVertexX = vertexPos.X();
      fReco_nuVertexY = vertexPos.Y();
      fReco_nuVertexZ = vertexPos.Z();

      //locate neutrino 

      nuSliceKey = slice.key();
      nuID = slicePFP->Self();

      // fill general variables 
      fReco_nPFParticles = slicePFPs.size();
      fReco_nPrimaryChildren = slicePFP->NumDaughters();

      // Get PFPs from the neutrino

      std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));
      
      // Actual analysis loop
      for(const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
      {
        if(nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
        {
          continue;
        }

        fReco_truthMatchedPFPg4ID.push_back(-999);

        // Check if PFP is reconstructed (as anything)
        // ============= clusterAna ==============

        std::vector<art::Ptr<recob::Cluster>> clusters(pfpClusterAssoc.at(nuSlicePFP.key()));

        std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit>>> planeHits;
        std::vector<art::Ptr<recob::Hit>> allHits;

        //Loop over the clusters in the plane and get the hits
        for (const art::Ptr<recob::Cluster> &cluster : clusters) {

          //Get the hits
          std::vector<art::Ptr<recob::Hit>> hits = clusterHitAssoc.at(cluster.key());

          //Get the plane.
          const geo::PlaneID::PlaneID_t plane(cluster->Plane().Plane);

          planeHits[plane].insert(planeHits[plane].end(), hits.begin(), hits.end());
          allHits.insert(allHits.end(), hits.begin(), hits.end());

        }

        // Now work out G4ID of the hits 

        if(!allHits.empty())
        {
          //shout(("allHits has size " + std::to_string(allHits.size())).c_str());
          int PFPg4ID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,allHits,true);
          fReco_truthMatchedPFPg4ID.back() = PFPg4ID;
          //shout(("this PFP has g4ID of " + std::to_string(PFPg4ID)).c_str());
          for(std::size_t i = 0; i < hitVector.size(); i++){
            int isInPFP = std::count(allHits.begin(),allHits.end(),hitVector[i]);
            if(isInPFP > 1){
              std::cout << "!!!!!!!!!!!!!!!weird!!!!!!!!!!!!!!!!!!!1" << std::endl;
            }
            if(isInPFP==1){
              fReco_hitHasPFPObject[i] = true;
              fReco_hitPFPObjectG4ID[i] = PFPg4ID;
            }
          }
        }

        // Do same G4ID but with a separate thing for which plane?

        


        // ============ end clusterAna ==============
        
        

        // Get tracks associated with this PFParticle
        std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

        // Check is one track per PFP
        if(tracks.size() != 1)
        {
          continue;
        }

        art::Ptr<recob::Track> track = tracks.at(0);
        
        AnalyseTrack(track, trackHitAssoc, trackCaloAssoc, hitVector);

        // Get shower associated with this PFParticle
        std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

        // Check is one shower per PFP
        if(showers.size() != 1)
        {
          continue;
        }

        art::Ptr<recob::Shower> shower = showers.at(0);
        
        AnalyseShower(shower, showerHitAssoc, showerCaloAssoc, hitVector);

      }

      break; //==========why?
    }

    // Assume only one neutrino per slice here - meed to fix this!
    if(nuID >=0)
    {
      break;
    }
  }
}

void nuclearFragments::RecoAnalysisIncShowers::shout(std::string message)
{
  std::cout << message << std::endl;
}


DEFINE_ART_MODULE(nuclearFragments::RecoAnalysisIncShowers)
