////////////////////////////////////////////////////////////////////////
// Class:       RecoAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoAnalysis_module.cc
//
// Generated on Wednesday 12th Feb 10:55:37 2025 by Anna Beever
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

// ROOT includes
#include <TTree.h>

constexpr double double_default = -999.0;

namespace nuclearFragments {
  class RecoAnalysis;
}


class nuclearFragments::RecoAnalysis : public art::EDAnalyzer {
public:
  explicit RecoAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoAnalysis(RecoAnalysis const&) = delete;
  RecoAnalysis(RecoAnalysis&&) = delete;
  RecoAnalysis& operator=(RecoAnalysis const&) = delete;
  RecoAnalysis& operator=(RecoAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Create output TTree
  TTree *fTree;

  // Tree variables
  unsigned int fEventID;
  unsigned int fReco_nPFParticles;
  unsigned int fReco_nPrimaryChildren;
  bool fIsCCQE = 0;
  bool fHasCluster = 0;

  std::vector<double> fReco_childTrackLengths;
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
  std::vector<int> fMC_trackID;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;
  std::string fNuGenLabel;
  std::string fLArGeantLabel;

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
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


nuclearFragments::RecoAnalysis::RecoAnalysis(fhicl::ParameterSet const& p)
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

void nuclearFragments::RecoAnalysis::analyze(art::Event const& e)
{
  // Set the event ID
  fEventID = e.id().event();

  ResetVariables();

  // Track-Hit map setup  
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
  std::vector<art::Ptr<recob::Hit>> hitVector;

  if(hitHandle.isValid())
  {
    art::fill_ptr_vector(hitVector,hitHandle);
  }

  for(const art::Ptr<recob::Hit> &hit : hitVector){
    fTrackHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
  }

  // Get event slices
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if(sliceHandle.isValid())
  {
    art::fill_ptr_vector(sliceVector,sliceHandle);
  }

  // Get associations between slices and PFParticles
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

  //Filling neutrino hierachy variables
  int nuID = -1;
  int nuSliceKey = -1;

  for(const art::Ptr<recob::Slice> &slice : sliceVector)
  {
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

    for(const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs)
    {
      const bool isPrimary(slicePFP->IsPrimary());
      const bool isNeutrino((std::abs(slicePFP->PdgCode() == 12) || (std::abs(slicePFP->PdgCode() == 14))));

      if(!(isPrimary && isNeutrino))
      {
        continue;
      }

      nuSliceKey = slice.key();
      nuID = slicePFP->Self();
      fReco_nPFParticles = slicePFPs.size();
      fReco_nPrimaryChildren = slicePFP->NumDaughters();

      // Finding nu slices to loop through
      art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
      art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);

      std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

      int trackCounter = 0;
      
      // Actual analysis loop
      for(const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
      {
        if(nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
        {
          continue;
        }

        // Get tracks associated with this PFParticle
        std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

        if(tracks.size() != 1)
        {
          continue;
        }

        art::Ptr<recob::Track> track = tracks.at(0);

        trackCounter++;

        const art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
        art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackLabel);

        std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());
        int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
        fReco_truthMatchedTrackID.push_back(trackID);
        fReco_childTrackCompleteness.push_back(Completeness(trackHits,trackID));
        fReco_childTrackPurity.push_back(Purity(trackHits,trackID));
        fReco_childTrackLengths.push_back(track->Length());

        art::FindManyP<anab::Calorimetry> trackCaloAssoc(trackHandle, e, fCalorimetryLabel);

        // Get the calorimetry object
        std::vector<art::Ptr<anab::Calorimetry>> calos = trackCaloAssoc.at(track.key());

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

      break;
    }

    // Assume only one neutrino per slice here
    if(nuID >=0)
    {
      break;
    }
  }

  //Truth level
  art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fNuGenLabel);
  std::vector<art::Ptr<simb::MCTruth>> truNuVector;

  if(truNuHandle.isValid())
  {
    art::fill_ptr_vector(truNuVector,truNuHandle);
  }

  art::FindManyP<simb::MCParticle> truNuParticleAssoc(truNuHandle, e, fLArGeantLabel);

  for(const art::Ptr<simb::MCTruth> &truNu : truNuVector){
    if(truNu.isNull()){
      continue;
    }

    const simb::MCNeutrino nu = truNu->GetNeutrino();
    if(nu.InteractionType() == simb::kCCQE){
      fIsCCQE = 1;
      std::cout << "CCQE" << std::endl;
    }

    std::vector<art::Ptr<simb::MCParticle>> particles = truNuParticleAssoc.at(truNu.key());

    for(const art::Ptr<simb::MCParticle> &particle : particles){
      if(particle->Process() != "primary" || particle->StatusCode() != 1){
        continue;
      }

      fMC_trackID.push_back(particle->TrackId());
      fMC_particlePDG.push_back(particle->PdgCode());
      fMC_particleMass.push_back(particle->Mass());
      fMC_particleE.push_back(particle->E());
      fMC_particleP.push_back(particle->P());
      double KE = particle->E() - particle->Mass();
      fMC_particleKE.push_back(KE);
      fTrackIDtoTruthPDGMap.insert({particle->TrackId(), particle->PdgCode()});
      fTrackIDtoTruthEMap.insert({particle->TrackId(), particle->E()});
      fTrackIDtoTruthPMap.insert({particle->TrackId(), particle->P()});
      fTrackIDtoTruthKEMap.insert({particle->TrackId(), KE});


      int checkIfReconstructed = std::count(fReco_truthMatchedTrackID.begin(), fReco_truthMatchedTrackID.end(), particle->TrackId());
      fMC_isReconstructed.push_back(checkIfReconstructed!=0);

      //filter for clusters
      if(particle->PdgCode() == 1000010020 || particle->PdgCode() == 1000010030 || particle->PdgCode() == 1000020030 || particle->PdgCode() == 1000020040){
        fHasCluster = 1;
      }

    }
  }

  for(std::size_t i=0; i<fReco_truthMatchedTrackID.size(); i++){
    fReco_truthMatchedPDG.push_back(fTrackIDtoTruthPDGMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedE.push_back(fTrackIDtoTruthEMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedP.push_back(fTrackIDtoTruthPMap[fReco_truthMatchedTrackID[i]]);
    fReco_truthMatchedKE.push_back(fTrackIDtoTruthKEMap[fReco_truthMatchedTrackID[i]]);
  }

  //True event type
  /*
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          ++_n_nu;
        }
    }

    ResizeVectors(nuVars, _n_nu);

    int nuCounter = 0;

    for(auto const& MCTruthHandle : MCTruthHandles)
      {
        std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

        for(auto const& mct : MCTruthVec)
          {
            if(mct->Origin() != 1)
              continue;

            AnalyseMCTruth(e, nuVars, mct, nuCounter, "nu");
            AnalyseNuReco(nuCounter);

            ++nuCounter;
          }
      }
  }*/


  // Fill tree
  fTree->Fill();

}

void nuclearFragments::RecoAnalysis::beginJob()
{
  // Get TFileService to create output TTree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  // Add branches to TTRee
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("isCCQE", &fIsCCQE);
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
  fTree->Branch("MC_trackID", &fMC_trackID);
}

void nuclearFragments::RecoAnalysis::endJob()
{
  // Implementation of optional member function here.
}

void nuclearFragments::RecoAnalysis::ResetVariables()
{
  // Set all trackCounters to zero for the current event
  fIsCCQE = 0;
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
  fMC_trackID.clear();
}

float nuclearFragments::RecoAnalysis::Completeness(std::vector< art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (fTrackHitsMap[trackID] == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(fTrackHitsMap[trackID]);
}

float nuclearFragments::RecoAnalysis::Purity(std::vector<art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }

  return (objectHits.size() == 0) ? double_default : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}

DEFINE_ART_MODULE(nuclearFragments::RecoAnalysis)
