////////////////////////////////////////////////////////////////////////
// Class:       PidAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        PidAnalysis_module.cc
//
// Written on Monday 19th January 2026 by Anna Beever
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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include <TTree.h>

constexpr double double_default = -999.0;

namespace nuclearFragments {
  class PidAnalysis;
}


class nuclearFragments::PidAnalysis : public art::EDAnalyzer {
public:
  explicit PidAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PidAnalysis(PidAnalysis const&) = delete;
  PidAnalysis(PidAnalysis&&) = delete;
  PidAnalysis& operator=(PidAnalysis const&) = delete;
  PidAnalysis& operator=(PidAnalysis&&) = delete;

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
  bool fIsClearCosmic = 0;
  double fReco_nuVertexX = double_default;
  double fReco_nuVertexY = double_default;
  double fReco_nuVertexZ = double_default;
  double fMC_nuVertexX = double_default;
  double fMC_nuVertexY = double_default;
  double fMC_nuVertexZ = double_default;

  std::vector<int> fRecoPFP_selfID;
  std::vector<int> fRecoPFP_parentID;
  std::vector<double> fReco_trackLength;
  std::vector<double> fReco_trackCompleteness;
  std::vector<double> fReco_trackPurity;
  std::vector<double> fReco_trackStartX;
  std::vector<double> fReco_trackStartY;
  std::vector<double> fReco_trackStartZ;
  std::vector<double> fReco_trackEndX;
  std::vector<double> fReco_trackEndY;
  std::vector<double> fReco_trackEndZ;
  std::vector<double> fReco_trackChisqMuon;
  std::vector<double> fReco_trackChisqProton;
  std::vector<double> fReco_trackChisqPion;
  std::vector<double> fReco_trackChisqKaon;
  std::vector<double> fReco_showerLengths;
  std::vector<double> fReco_showerCompleteness;
  std::vector<double> fReco_showerPurity;
  std::vector<int> fReco_pandoraPDG;
  std::vector<int> fReco_truthMatchedShowerG4ID;
  std::vector<int> fReco_truthMatchedTrackG4ID;
  std::vector<int> fRecoPFP_truthMatchedG4ID;
  std::vector<int> fMC_allParticleTrackG4ID;
  std::vector<int> fMC_allParticlePDG;
  std::vector<int> fMC_allParticleMotherID;
  std::vector<double> fMC_allParticleP;
  std::vector<double> fMC_allParticleEk;
  std::vector<int> fMC_primaryParticleTrackG4ID;
  std::vector<int> fMC_primaryParticleParentG4ID;
  std::vector<double> fMC_primaryParticleVx;
  std::vector<double> fMC_primaryParticleVy;
  std::vector<double> fMC_primaryParticleVz;
  std::vector<double> fMC_primaryParticleEndX;
  std::vector<double> fMC_primaryParticleEndY;
  std::vector<double> fMC_primaryParticleEndZ;
  std::vector<bool> fMC_primaryIsReconstructedAsTrack;
  std::vector<bool> fMC_primaryIsReconstructedAsShower;
  std::vector<bool> fMC_primaryIsReconstructedAsPFP;  
  std::vector<bool> fMC_allParticleIsReconstructedAsPFP;  

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fChisqLabel;
  std::string fShowerLabel;
  //std::string fShowerLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;
  std::string fVertexLabel;
  std::string fClusterLabel;
  std::string fSpacePointLabel;
  std::string fRecoHitLabel;
  std::string fLArGeantLabel;
  std::string fNuGenLabel;

  // Maps
  std::map<int,int> fTrackHitsMap;

  // Functions 
  void ResetVariables();
  float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  void AnalyseTrack(const art::Ptr<recob::Track> &track, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, art::FindManyP<anab::ParticleID> trackChisqAssoc);
  void AnalyseShower(const art::Ptr<recob::Shower> &shower, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc);
  void AnalyseTruNu(const art::Ptr<simb::MCTruth> &truNu, art::FindManyP<simb::MCParticle> truNuParticleAssoc);
  void AnalyseMCTruth(art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle, art::FindManyP<simb::MCParticle> truNuParticleAssoc);
  void AnalyseReco(art::ValidHandle<std::vector<recob::Hit>> hitHandle, art::ValidHandle<std::vector<recob::Slice>> sliceHandle, art::ValidHandle<std::vector<recob::Track>> trackHandle, art::FindManyP<recob::PFParticle> slicePFPAssoc, art::FindManyP<recob::Track> pfpTrackAssoc, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, art::FindOneP<recob::Vertex> pfpVertexAssoc, art::ValidHandle<std::vector<recob::Shower>> showerHandle, art::FindManyP<recob::Shower> pfpShowerAssoc, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc, art::FindManyP<recob::Cluster> pfpClusterAssoc, art::FindManyP<recob::Hit> clusterHitAssoc, art::FindManyP<recob::SpacePoint> pfpSpacePointAssoc);
  void shout(std::string message);

  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


nuclearFragments::PidAnalysis::PidAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fChisqLabel(p.get<std::string>("ChisqLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  //fShowerLabel(p.get<std::string>("ShowerLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fHitLabel(p.get<std::string>("HitLabel")),
  fVertexLabel(p.get<std::string>("VertexLabel")),
  fClusterLabel(p.get<std::string>("ClusterLabel")),
  fSpacePointLabel(p.get<std::string>("SpacePointLabel")),
  fRecoHitLabel(p.get<std::string>("RecoHitLabel")),
  fLArGeantLabel(p.get<std::string>("LArGeantLabel")),
  fNuGenLabel(p.get<std::string>("NuGenLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void nuclearFragments::PidAnalysis::analyze(art::Event const& e)
{

  ResetVariables();

  // Set the event ID
  fEventID = e.id().event();

  std::cout << "=============== new event ===============" << std::endl;

  // ========================================================================
  // Handles
  // ========================================================================

  // Get hit labelled objects for event
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);

  // Get spacepoint labelled objects for event
  art::ValidHandle<std::vector<recob::SpacePoint>> spacepointHandle = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);

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

  // Get associations between PFPs and spacepoints
  art::FindManyP<recob::SpacePoint>pfpSpacePointAssoc(pfpHandle, e, fSpacePointLabel);

  // Get associations between spacepoints and hits
  art::FindManyP<recob::Hit>spacepointHitAssoc(spacepointHandle, e, fRecoHitLabel);

  // Get associations between tracks and chisquares
  art::FindManyP<anab::ParticleID>trackChisqAssoc(trackHandle, e, fChisqLabel);

  // ========================================================================
  // Main analysis
  // ========================================================================

  // Fill hits
  std::vector<art::Ptr<recob::Hit>> hitVector;

  if(hitHandle.isValid())
  {
    art::fill_ptr_vector(hitVector,hitHandle);
  }

  // Track-Hit map setup

  for(const art::Ptr<recob::Hit> &hit : hitVector){
    fTrackHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
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

      const art::Ptr<recob::Vertex> vertex = pfpVertexAssoc.at(slicePFP.key());
      geo::Point_t vertexPos = vertex.isNonnull() ? vertex->position() : geo::Point_t(double_default, double_default, double_default);
      fReco_nuVertexX = vertexPos.X();
      fReco_nuVertexY = vertexPos.Y();
      fReco_nuVertexZ = vertexPos.Z();

      nuSliceKey = slice.key();
      nuID = slicePFP->Self();

      fReco_nPFParticles = slicePFPs.size();
      fReco_nPrimaryChildren = slicePFP->NumDaughters();

      // Get PFPs from the neutrino

      std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));
      
      for(const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
      {
        fRecoPFP_selfID.push_back(nuSlicePFP->Self());
        fRecoPFP_parentID.push_back(nuSlicePFP->Parent());
        fRecoPFP_truthMatchedG4ID.push_back(-999);

        //Get the clusters
        std::vector<art::Ptr<recob::Cluster>> clusters(pfpClusterAssoc.at(nuSlicePFP.key()));

        std::vector<art::Ptr<recob::Hit>> allHits;

        //Loop over the clusters get the hits
        for (const art::Ptr<recob::Cluster> &cluster : clusters) {

          //Get the hits
          std::vector<art::Ptr<recob::Hit>> hits = clusterHitAssoc.at(cluster.key());

          allHits.insert(allHits.end(), hits.begin(), hits.end());

        }

        // Now work out G4ID of the hits 
        
        if(!allHits.empty())
        {

          int PFPg4ID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,allHits,true);
          fRecoPFP_truthMatchedG4ID.back() = PFPg4ID;
  
        }

        // Get tracks associated with this PFParticle
        std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

        // Check is one track per PFP
        if(tracks.size() != 1)
        {
          continue;
        }

        fReco_pandoraPDG.push_back(nuSlicePFP->PdgCode());

        art::Ptr<recob::Track> track = tracks.at(0);
        
        AnalyseTrack(track, trackHitAssoc, trackCaloAssoc, trackChisqAssoc);

        // Get shower associated with this PFParticle
        std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

        // Check is one shower per PFP
        if(showers.size() != 1)
        {
          continue;
        }

        art::Ptr<recob::Shower> shower = showers.at(0);
        
        AnalyseShower(shower, showerHitAssoc, showerCaloAssoc);

      }

      break;
    }

    // Assume only one neutrino per slice here - meed to fix this!
    if(nuID >=0)
    {
      break;
    }
  }

  AnalyseMCTruth(truNuHandle, truNuParticleAssoc);

  // Fill event tree
  fEventTree->Fill();

}

void nuclearFragments::PidAnalysis::respondToOpenInputFile(const art::FileBlock& fb)
{
  fInputFile = fb.fileName();
}

void nuclearFragments::PidAnalysis::beginJob()
{
  // Get TFileService to create output TTree
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("eventTree", "Output TTree");

  // Add branches to TTree
  fEventTree->Branch("eventID", &fEventID);
  fEventTree->Branch("inputFile", &fInputFile);
  fEventTree->Branch("reco_nPFParticles", &fReco_nPFParticles);
  fEventTree->Branch("reco_nPrimaryChildren", &fReco_nPrimaryChildren);
  fEventTree->Branch("reco_trackLength", &fReco_trackLength);
  fEventTree->Branch("reco_trackCompleteness", &fReco_trackCompleteness);
  fEventTree->Branch("reco_trackPurity", &fReco_trackPurity);
  fEventTree->Branch("reco_trackStartX", &fReco_trackStartX);
  fEventTree->Branch("reco_trackStartY", &fReco_trackStartY);
  fEventTree->Branch("reco_trackStartZ", &fReco_trackStartZ);
  fEventTree->Branch("reco_trackEndX", &fReco_trackEndX);
  fEventTree->Branch("reco_trackEndY", &fReco_trackEndY);
  fEventTree->Branch("reco_trackEndZ", &fReco_trackEndZ);
  fEventTree->Branch("reco_trackChisqMuon", &fReco_trackChisqMuon);
  fEventTree->Branch("reco_trackChisqProton", &fReco_trackChisqProton);
  fEventTree->Branch("reco_trackChisqPion", &fReco_trackChisqPion);
  fEventTree->Branch("reco_trackChisqKaon", &fReco_trackChisqKaon);
  fEventTree->Branch("reco_showerLengths", &fReco_showerLengths);
  fEventTree->Branch("reco_showerCompleteness", &fReco_showerCompleteness);
  fEventTree->Branch("reco_showerPurity", &fReco_showerPurity);
  fEventTree->Branch("reco_pandoraPDG", &fReco_pandoraPDG);
  fEventTree->Branch("reco_truthMatchedTrackG4ID", &fReco_truthMatchedTrackG4ID);
  fEventTree->Branch("reco_truthMatchedShowerG4ID", &fReco_truthMatchedShowerG4ID);
  fEventTree->Branch("recoPFP_truthMatchedG4ID", &fRecoPFP_truthMatchedG4ID);
  fEventTree->Branch("recoPFP_selfID", &fRecoPFP_selfID);
  fEventTree->Branch("recoPFP_parentID", &fRecoPFP_parentID);
  fEventTree->Branch("reco_nuVertexX", &fReco_nuVertexX);
  fEventTree->Branch("reco_nuVertexY", &fReco_nuVertexY);
  fEventTree->Branch("reco_nuVertexZ", &fReco_nuVertexZ);
  fEventTree->Branch("MC_nuVertexX", &fMC_nuVertexX);
  fEventTree->Branch("MC_nuVertexY", &fMC_nuVertexY);
  fEventTree->Branch("MC_nuVertexZ", &fMC_nuVertexZ);
  fEventTree->Branch("MC_primaryParticleTrackG4ID", &fMC_primaryParticleTrackG4ID);
  fEventTree->Branch("MC_allParticleTrackG4ID", &fMC_allParticleTrackG4ID);
  fEventTree->Branch("MC_allParticlePDG", &fMC_allParticlePDG);
  fEventTree->Branch("MC_allParticleMotherID", &fMC_allParticleMotherID);
  fEventTree->Branch("MC_allParticleP", &fMC_allParticleP);
  fEventTree->Branch("MC_allParticleEk", &fMC_allParticleEk);
  fEventTree->Branch("MC_primaryParticleParentG4ID", &fMC_primaryParticleParentG4ID);
  fEventTree->Branch("MC_primaryParticleVx", &fMC_primaryParticleVx);
  fEventTree->Branch("MC_primaryParticleVy", &fMC_primaryParticleVy);
  fEventTree->Branch("MC_primaryParticleVz", &fMC_primaryParticleVz);
  fEventTree->Branch("MC_primaryParticleEndX", &fMC_primaryParticleEndX);
  fEventTree->Branch("MC_primaryParticleEndY", &fMC_primaryParticleEndY);
  fEventTree->Branch("MC_primaryParticleEndZ", &fMC_primaryParticleEndZ);
  fEventTree->Branch("MC_primaryIsReconstructedAsTrack", &fMC_primaryIsReconstructedAsTrack);
  fEventTree->Branch("MC_primaryIsReconstructedAsShower", &fMC_primaryIsReconstructedAsShower);
  fEventTree->Branch("MC_primaryIsReconstructedAsPFP", &fMC_primaryIsReconstructedAsPFP);
  fEventTree->Branch("MC_allParticleIsReconstructedAsPFP", &fMC_allParticleIsReconstructedAsPFP);
}

void nuclearFragments::PidAnalysis::endJob()
{
  // Implementation of optional member function here.
}

void nuclearFragments::PidAnalysis::ResetVariables()
{
  fTrackHitsMap.clear();
  fEventID = 0;
  fReco_nuVertexX = double_default;
  fReco_nuVertexY = double_default;
  fReco_nuVertexZ = double_default;
  fMC_nuVertexX = double_default;
  fMC_nuVertexY = double_default;
  fMC_nuVertexZ = double_default;
  fReco_nPFParticles = 0;
  fReco_nPrimaryChildren = 0;
  fReco_trackLength.clear();
  fReco_trackCompleteness.clear();
  fReco_trackPurity.clear();
  fReco_trackStartX.clear();
  fReco_trackStartY.clear();
  fReco_trackStartZ.clear();
  fReco_trackEndX.clear();
  fReco_trackEndY.clear();
  fReco_trackEndZ.clear();
  fReco_trackChisqMuon.clear();
  fReco_trackChisqProton.clear();
  fReco_trackChisqPion.clear();
  fReco_trackChisqKaon.clear();
  fReco_showerLengths.clear();
  fReco_showerCompleteness.clear();
  fReco_showerPurity.clear();
  fReco_pandoraPDG.clear();
  fReco_truthMatchedTrackG4ID.clear();
  fReco_truthMatchedShowerG4ID.clear();
  fRecoPFP_truthMatchedG4ID.clear();
  fRecoPFP_selfID.clear();
  fRecoPFP_parentID.clear();
  fMC_allParticleTrackG4ID.clear();
  fMC_allParticlePDG.clear();
  fMC_allParticleMotherID.clear();
  fMC_allParticleP.clear();
  fMC_allParticleEk.clear();
  fMC_primaryParticleTrackG4ID.clear();
  fMC_primaryParticleParentG4ID.clear();
  fMC_primaryParticleVx.clear();
  fMC_primaryParticleVy.clear();
  fMC_primaryParticleVz.clear();
  fMC_primaryParticleEndX.clear();
  fMC_primaryParticleEndY.clear();
  fMC_primaryParticleEndZ.clear();
  fMC_primaryIsReconstructedAsTrack.clear();
  fMC_primaryIsReconstructedAsShower.clear();
  fMC_primaryIsReconstructedAsPFP.clear();
  fMC_allParticleIsReconstructedAsPFP.clear();
}

float nuclearFragments::PidAnalysis::Completeness(std::vector< art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (fTrackHitsMap[trackID] == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(fTrackHitsMap[trackID]);
}

float nuclearFragments::PidAnalysis::Purity(std::vector<art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (objectHits.size() == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}

void nuclearFragments::PidAnalysis::AnalyseTrack(const art::Ptr<recob::Track> &track, art::FindManyP<recob::Hit> trackHitAssoc, art::FindManyP<anab::Calorimetry> trackCaloAssoc, art::FindManyP<anab::ParticleID> trackChisqAssoc)
{
  // Fill in track variables
  fReco_trackStartX.push_back((track->Start()).X());
  fReco_trackStartY.push_back((track->Start()).Y());
  fReco_trackStartZ.push_back((track->Start()).Z());

  fReco_trackEndX.push_back((track->End()).X());
  fReco_trackEndY.push_back((track->End()).Y());
  fReco_trackEndZ.push_back((track->End()).Z());

  std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());

  int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);

  fReco_truthMatchedTrackG4ID.push_back(trackID);
  std::cout << "track ID: " << trackID << std::endl;
  fReco_trackCompleteness.push_back(Completeness(trackHits,trackID));
  fReco_trackPurity.push_back(Purity(trackHits,trackID));
  fReco_trackLength.push_back(track->Length());
  std::cout << "length is: " << track->Length() << std::endl;
  fReco_trackChisqMuon.push_back(-1.);
  fReco_trackChisqProton.push_back(-1.);
  fReco_trackChisqPion.push_back(-1.);
  fReco_trackChisqKaon.push_back(-1.);

  // Get the calorimetry objects associated with this track
  std::vector<art::Ptr<anab::Calorimetry>> calos = trackCaloAssoc.at(track.key());

  if(calos.size() != 3){
    //std::cout << "not enough calo entries, it has size: " << calos.size() << std::endl;
  }
  const size_t maxHits = calos.size() != 3 ? -1 : std::max({calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size()});
  const int bestPlane  = calos.size() != 3 ? -1 : (calos[2]->dEdx().size() == maxHits) ? 2 : (calos[0]->dEdx().size() == maxHits) ? 0 : (calos[1]->dEdx().size() == maxHits) ? 1 : -1;

  const std::vector<art::Ptr<anab::ParticleID>> chi2s = trackChisqAssoc.at(track.key());
  const std::vector<anab::sParticleIDAlgScores> AlgScoresVec = chi2s[bestPlane]->ParticleIDAlgScores();

  //std::cout << "maxHits: " << maxHits << std::endl;
  //std::cout << "alg scores vec size: " << AlgScoresVec.size() << std::endl;

  std::vector<std::pair<int, float>> chi2vec;

  for(size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
  {
    const anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

    //std::cout << AlgScore.fAlgName << std::endl;

    if(AlgScore.fAlgName == "Chi2")
      {
        chi2vec.push_back({AlgScore.fAssumedPdg, AlgScore.fValue});
        //std::cout << "Assumed PDG: " << AlgScore.fAssumedPdg << std::endl;
        //std::cout << "Value: " << AlgScore.fValue << std::endl;
 
        switch(AlgScore.fAssumedPdg)
          {
          case 13:
            fReco_trackChisqMuon.back() = AlgScore.fValue;
            break;
          case 211:
            fReco_trackChisqPion.back() = AlgScore.fValue;
            break;
          case 321:
            fReco_trackChisqKaon.back() = AlgScore.fValue;
            break;
          case 2212:
            fReco_trackChisqProton.back() = AlgScore.fValue;
            break;
          }
      }
  }

  //std::cout << "muon PDG score: " << fReco_trackChisqMuon.back() << std::endl;
  //std::cout << "proton PDG score: " << fReco_trackChisqProton.back() << std::endl;
  //std::cout << "pion PDG score: " << fReco_trackChisqPion.back() << std::endl;
  //std::cout << "kaon PDG score: " << fReco_trackChisqKaon.back() << std::endl;

  /*if(fReco_trackChisqMuon.back() == -1){
    for(const art::Ptr<anab::Calorimetry> &calo : calos)
    {
      //const int plane = calo->PlaneID().Plane;
      //std::cout << "plane: " << plane << std::endl;
      //std::vector<double> dEdx (calo->dEdx().begin(), calo->dEdx().end());
      //std::vector<double> resRange (calo->ResidualRange().begin(), calo->ResidualRange().end());
      //for(size_t i=0; i < dEdx.size(); i++){
        //std::cout << dEdx[i] << std::endl;
        //std::cout << resRange[i] << std::endl;
      //}
    }
  }*/



}


void nuclearFragments::PidAnalysis::AnalyseShower(const art::Ptr<recob::Shower> &shower, art::FindManyP<recob::Hit> showerHitAssoc, art::FindManyP<anab::Calorimetry> showerCaloAssoc)
{

  // Get hits associated with this shower
  std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssoc.at(shower.key());

  // Fill shower variables
  int showerID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,showerHits,true);
  fReco_truthMatchedShowerG4ID.push_back(showerID);
  fReco_showerCompleteness.push_back(Completeness(showerHits,showerID));
  fReco_showerPurity.push_back(Purity(showerHits,showerID));
  fReco_showerLengths.push_back(shower->Length());
  

}

void nuclearFragments::PidAnalysis::AnalyseTruNu(const art::Ptr<simb::MCTruth> &truNu, art::FindManyP<simb::MCParticle> truNuParticleAssoc)
{

  // Get neutrino
  const simb::MCNeutrino mcn = truNu->GetNeutrino();
  const simb::MCParticle nu = mcn.Nu();

  fMC_nuVertexX = nu.Vx();
  fMC_nuVertexY = nu.Vy();
  fMC_nuVertexZ = nu.Vz();

  // Get particles associated with true neutrino
  std::vector<art::Ptr<simb::MCParticle>> particles = truNuParticleAssoc.at(truNu.key());


  // Loop over the neutrino's particles 
  for(const art::Ptr<simb::MCParticle> &particle : particles){


    if(particle->StatusCode() != 1){
      continue;
    }

    //truth info on all particles propagated by Geant4
    fMC_allParticleTrackG4ID.push_back(particle->TrackId());
    fMC_allParticlePDG.push_back(particle->PdgCode());
    fMC_allParticleMotherID.push_back(particle->Mother());
    fMC_allParticleP.push_back(particle->P());
    double Ek = particle->E() - particle->Mass();
    fMC_allParticleEk.push_back(Ek);

    int checkIfReconstructedAsPFP_allParticles = std::count(fRecoPFP_truthMatchedG4ID.begin(), fRecoPFP_truthMatchedG4ID.end(), particle->TrackId());
    fMC_allParticleIsReconstructedAsPFP.push_back(checkIfReconstructedAsPFP_allParticles!=0);

    // Only get primary particles - the ones leaving the nucleus (?)
    if(particle->Process() != "primary" || particle->StatusCode() != 1){
      continue;
    }

    //shout(("particle of G4ID " + std::to_string(particle->TrackId())).c_str());
    //shout(("particle of PDG code " + std::to_string(particle->PdgCode())).c_str());

    // Fill MC particle truth info

    fMC_primaryParticleTrackG4ID.push_back(particle->TrackId());
    fMC_primaryParticleParentG4ID.push_back(particle->Mother());
    fMC_primaryParticleVx.push_back(particle->Vx());
    fMC_primaryParticleVy.push_back(particle->Vy());
    fMC_primaryParticleVz.push_back(particle->Vz());
    fMC_primaryParticleEndX.push_back(particle->EndX());
    fMC_primaryParticleEndY.push_back(particle->EndY());
    fMC_primaryParticleEndZ.push_back(particle->EndZ());
    
    int checkIfReconstructedAsTrack = std::count(fReco_truthMatchedTrackG4ID.begin(), fReco_truthMatchedTrackG4ID.end(), particle->TrackId());
    fMC_primaryIsReconstructedAsTrack.push_back(checkIfReconstructedAsTrack!=0);

    int checkIfReconstructedAsShower = std::count(fReco_truthMatchedShowerG4ID.begin(), fReco_truthMatchedShowerG4ID.end(), particle->TrackId());
    fMC_primaryIsReconstructedAsShower.push_back(checkIfReconstructedAsShower!=0);

    int checkIfReconstructedAsPFP = std::count(fRecoPFP_truthMatchedG4ID.begin(), fRecoPFP_truthMatchedG4ID.end(), particle->TrackId());
    fMC_primaryIsReconstructedAsPFP.push_back(checkIfReconstructedAsPFP!=0);

  }
  
}

void nuclearFragments::PidAnalysis::AnalyseMCTruth(art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle, art::FindManyP<simb::MCParticle> truNuParticleAssoc)
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

void nuclearFragments::PidAnalysis::shout(std::string message)
{
  std::cout << message << std::endl;
}


DEFINE_ART_MODULE(nuclearFragments::PidAnalysis)
