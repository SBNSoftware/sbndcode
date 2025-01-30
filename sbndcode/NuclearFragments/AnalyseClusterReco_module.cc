////////////////////////////////////////////////////////////////////////
// Class:       AnalyseClusterReco
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyseClusterReco_module.cc
//
// Generated at Thu Jan 16 10:55:37 2025 by Anna Beever using cetskelgen
// from cetlib version 3.18.02.
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

constexpr float float_default = -999.0f;

namespace testReco {
  class AnalyseClusterReco;
}


class testReco::AnalyseClusterReco : public art::EDAnalyzer {
public:
  explicit AnalyseClusterReco(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseClusterReco(AnalyseClusterReco const&) = delete;
  AnalyseClusterReco(AnalyseClusterReco&&) = delete;
  AnalyseClusterReco& operator=(AnalyseClusterReco const&) = delete;
  AnalyseClusterReco& operator=(AnalyseClusterReco&&) = delete;

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
  unsigned int fNPFParticles;
  unsigned int fNPrimaryChildren;
  bool fThreeTracks;

  std::vector<double> fChildTrackLengths;
  std::vector<float> fChildTrackCompleteness;
  std::vector<float> fChildTrackPurity;
  std::vector<std::vector<float>> fChildTrackdEdx;
  std::vector<std::vector<float>> fChildTrackResRange;
  std::vector<bool> fChildTrackIsLongest;
  std::vector<bool> fChildTrackIsMiddle;
  std::vector<bool> fChildTrackIsShortest;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;

  // Maps
  std::map<int,int> fTrackHitsMap;

  // Functions 
  void ResetVariables();
  float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


testReco::AnalyseClusterReco::AnalyseClusterReco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fHitLabel(p.get<std::string>("HitLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void testReco::AnalyseClusterReco::analyze(art::Event const& e)
{
  // Set the event ID
  fEventID = e.id().event();

  ResetVariables();

  // Track-Hit map setup  
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
  std::vector<art::Ptr<recob::Hit>> hitVector;

  if(hitHandle.isValid())
  {
    std::cout << "valid hit handle" << std::endl;
    art::fill_ptr_vector(hitVector,hitHandle);
  }

  int hitCount = 0;
  for(const art::Ptr<recob::Hit> &hit : hitVector){
    hitCount++;
    //if(TruthMatchUtils::TrueParticleID(clockData,hit,true) >= 0){
    fTrackHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
    std::cout << TruthMatchUtils::TrueParticleID(clockData,hit,true) << ",";
    //}
  }
  std::cout << "end" << std::endl;
  std::cout << fTrackHitsMap[0] << std::endl;
  std::cout << fTrackHitsMap[1] << std::endl;
  std::cout << fTrackHitsMap[2] << std::endl;
  std::cout << fTrackHitsMap[3] << std::endl;
  std::cout << hitCount << std::endl;

  // Get event slices
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if(sliceHandle.isValid())
  {
    std::cout << "valid slice handle" << std::endl;
    art::fill_ptr_vector(sliceVector,sliceHandle);
  }

  // Get associations between slices and PFParticles
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

  //Filling neutrino hierachy variables
  int nuID = -1;
  int nuSliceKey = -1;

  for(const art::Ptr<recob::Slice> &slice : sliceVector)
  {
    std::cout << "new slice" << std::endl;
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

    for(const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs)
    {
      const bool isPrimary(slicePFP->IsPrimary());
      const bool isNeutrino((std::abs(slicePFP->PdgCode() == 12) || (std::abs(slicePFP->PdgCode() == 14))));

      if(!(isPrimary && isNeutrino))
      {
        continue;
      }
      std::cout << "found neutrino" << std::endl;

      nuSliceKey = slice.key();
      nuID = slicePFP->Self();
      fNPFParticles = slicePFPs.size();
      fNPrimaryChildren = slicePFP->NumDaughters();

      // Finding nu slices to loop through
      art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
      art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);

      std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

      // Identifying track length hierachy
      int longestID = std::numeric_limits<int>::lowest();
      int middleID = std::numeric_limits<int>::lowest();
      int shortestID = std::numeric_limits<int>::lowest();
      double longestLength = std::numeric_limits<double>::lowest();
      double shortestLength = std::numeric_limits<double>::max();
      int trackCounter = 0;

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

        // Find the middle length (if any)
        if(trackCounter == 3){
          if(track->Length() > shortestLength && track->Length() < longestLength){
            middleID = track->ID();
            //middleLength = track->Length();
          }
          else if(track->Length() < shortestLength){
            middleID = shortestID;
            //middleLength = shortestLength;
          }
          else if(track->Length() > longestLength){
            middleID = longestID;
            //middleLength = longestLength;
          }
        }

        // Check if this track is longer than the current longest
        if(track->Length() > longestLength)
        {
          longestID = track->ID();
          longestLength = track->Length();
        }
        if(track->Length() < shortestLength)
        {
          shortestID = track->ID();
          shortestLength = track->Length();
        }

      }

      if(trackCounter == 3){
        fThreeTracks = 1;
      }
      
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

        art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
        art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackLabel);

        std::vector<art::Ptr<recob::Hit> > trackHits = trackHitAssoc.at(track.key());
        int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
        fChildTrackCompleteness.push_back(Completeness(trackHits,trackID));
        float comp = Completeness(trackHits, trackID);
        std::cout << "completeness: " << comp << std::endl;
        fChildTrackPurity.push_back(Purity(trackHits,trackID));
        float pur = Purity(trackHits, trackID);
        std::cout << "purity: " << pur << std::endl;

        fChildTrackLengths.push_back(track->Length());
        //std::cout << "pushed back length: " << track->Length() << std::endl;
        fChildTrackIsLongest.push_back(track->ID() == longestID);
        fChildTrackIsMiddle.push_back(track->ID() == middleID);
        fChildTrackIsShortest.push_back(track->ID() == shortestID);
        //std::cout << "isLongest: " << (track->ID() == longestID) << std::endl;
        //std::cout << "isMiddle: " << (track->ID() == middleID) << std::endl;
        //std::cout << "isShortest: " << (track->ID() == shortestID) << std::endl;
        //std::cout << "size of longest vector: " << fChildTrackIsLongest.size() << std::endl;

        //art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
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
          fChildTrackdEdx.push_back(calo->dEdx());
          fChildTrackResRange.push_back(calo->ResidualRange());
          std::vector<float> resRange = calo->ResidualRange();
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

  // Fill tree
  fTree->Fill();

}

void testReco::AnalyseClusterReco::beginJob()
{
  // Get TFileService to create output TTree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  // Add branches to TTRee
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaryChildren", &fNPrimaryChildren);
  fTree->Branch("childTrackLengths", &fChildTrackLengths);
  fTree->Branch("childTrackCompleteness", &fChildTrackCompleteness);
  fTree->Branch("childTrackPurity", &fChildTrackPurity);
  fTree->Branch("childTrackIsLongest", &fChildTrackIsLongest);
  fTree->Branch("childTrackIsMiddle", &fChildTrackIsMiddle);
  fTree->Branch("childTrackIsShortest", &fChildTrackIsShortest);
  fTree->Branch("threeTracks", &fThreeTracks);
  fTree->Branch("childTrackdEdx", &fChildTrackdEdx);
  fTree->Branch("childTrackResRange", &fChildTrackResRange);
}

void testReco::AnalyseClusterReco::endJob()
{
  // Implementation of optional member function here.
}

void testReco::AnalyseClusterReco::ResetVariables()
{
  // Set all trackCounters to zero for the current event
  fNPFParticles = 0;
  fNPrimaryChildren = 0;
  fThreeTracks = 0;
  fChildTrackLengths.clear();
  fChildTrackCompleteness.clear();
  fChildTrackPurity.clear();
  fChildTrackIsLongest.clear();
  fChildTrackIsMiddle.clear();
  fChildTrackIsShortest.clear();
  fChildTrackdEdx.clear();
  fChildTrackResRange.clear();
  fTrackHitsMap.clear();
}

float testReco::AnalyseClusterReco::Completeness(std::vector< art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  //std::cout << "COMPLETENESS" << std::endl;
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  //std::cout << "number of hits: " << objectHits.size() << std::endl;

  std::cout << "number of hits measured: " << objectHitsMap[trackID] << std::endl;
  std::cout << "object hits map size: " << objectHitsMap.size() << std::endl;

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    //if(TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true) >= 0){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
    std::cout << TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true) << ",";
    //}
  }
  std::cout << "end" << std::endl;

  std::cout << "object hits map size: " << objectHitsMap.size() << std::endl;
  //std::cout << "number of object hits: " << objectHitsMap[trackID] << std::endl;
  //std::cout << "track ID: " << trackID << std::endl;
  //std::cout << fTrackHitsMap[trackID] << std::endl;

  std::cout << "number of hits in object: " << objectHits.size() << std::endl;
  std::cout << "number of hits measured: " << objectHitsMap[trackID] << std::endl;
  std::cout << "number of hits: " << fTrackHitsMap[trackID] << std::endl;

  return (fTrackHitsMap[trackID] == 0) ? float_default : objectHitsMap[trackID]/static_cast<float>(fTrackHitsMap[trackID]);
}

float testReco::AnalyseClusterReco::Purity(std::vector<art::Ptr<recob::Hit>> const &objectHits, int const &trackID)
{
  //std::cout << "PURITY" << std::endl;
  std::map<int,int> objectHitsMap;
  objectHitsMap.clear();

  //std::cout << "number of hits: " << objectHits.size() << std::endl;

  for(unsigned int i = 0; i < objectHits.size(); ++i){
    //if(TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true) >= 0){
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
    //}
  }

  //std::cout << "track ID: " << trackID << std::endl;
  //std::cout << "number of hits from the true track: " << objectHitsMap[trackID] << std::endl;

  if(objectHitsMap[trackID] > static_cast<float>(objectHits.size())){
    std::cout << "number of hits: " << objectHits.size() << std::endl;
    std::cout << "number of hits from the true track: " << objectHitsMap[trackID] << std::endl;
  }

  return (objectHits.size() == 0) ? float_default : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

DEFINE_ART_MODULE(testReco::AnalyseClusterReco)
