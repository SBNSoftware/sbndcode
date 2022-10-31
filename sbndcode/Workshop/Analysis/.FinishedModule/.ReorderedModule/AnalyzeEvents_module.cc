////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Wed Oct 13 08:28:13 2021 by Edward Tyley using cetskelgen
// from  version .
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

// Additional Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Additional LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Slice.h"

// ROOT includes
#include <TH1F.h>
#include <TTree.h>

// STL includes
#include <string>
#include <vector>

namespace test {
class AnalyzeEvents;
}

class test::AnalyzeEvents : public art::EDAnalyzer {
  public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  private:
  // Declare member data here.

  // Create out output tree
  TTree* fTree;
  // Create out output histogram
  TH1F* fTrackLengthHist;

  // Tree Variables
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;
  float fT0;

  std::vector<float> fDaughterTrackLengths;
  std::vector<bool> fDaughterLongestTrack;
  std::vector<int> fDaughterTrackTruePDG;

  std::vector<std::vector<float>> fDaughterTrackdEdx;
  std::vector<std::vector<float>> fDaughterTrackResidualRange;

  // Define input labels
  const std::string fPFParticleLabel;
  const std::string fTrackLabel;
  const std::string fCaloLabel;
  const std::string fSliceLabel;
  const std::string fOptLabel;

  // Write some member functions to make our code more readable
  void ResetTree();
  const int GetNeutrinoID(const std::vector<art::Ptr<recob::PFParticle>>& pfpVec);
  const std::vector<art::Ptr<recob::Track>> GetDaughterTracks(
      const size_t neutrinoID,
      const std::vector<art::Ptr<recob::PFParticle>>& pfpVec,
      const art::FindManyP<recob::Track>& pfpTrackAssns) const;
  const art::Ptr<recob::Slice> GetNeutrinoSlice(
      const size_t neutrinoID,
      const std::vector<art::Ptr<recob::PFParticle>>& pfpVec,
      const art::FindManyP<recob::Slice>& pfpSliceAssns) const;
  const int GetLongestTrackID(const std::vector<art::Ptr<recob::Track>>& daughterTracks) const;
};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
    : EDAnalyzer { p }
    // Initialise out input labels by reading the fhicl parameters
    , fPFParticleLabel(p.get<std::string>("PFParticleLabel"))
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fCaloLabel(p.get<std::string>("CalorimetryLabel"))
    , fSliceLabel(p.get<std::string>("SliceLabel"))
    , fOptLabel(p.get<std::string>("OptLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{

  // Initialise the servies we need for getting truth information
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventoryService;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // Implementation of required member function here.
  fEventID = e.id().event();

  // Reset all of our variables to 0 or empty vectors
  this->ResetTree();

  // Load the PFParticles from pandora
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  // Load all of the tracks from pandoraTrack
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> trackVec;
  if (e.getByLabel(fTrackLabel, trackHandle))
    art::fill_ptr_vector(trackVec, trackHandle);

  // Load all of the slices from pandora
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  if (e.getByLabel(fSliceLabel, sliceHandle))
    art::fill_ptr_vector(sliceVec, sliceHandle);

  // Load the associations between PFPs, Tracks and Calorimetries
  // Load the associations between PFPs, Slices and T0
  art::FindManyP<recob::Track> pfpTrackAssns(pfpVec, e, fTrackLabel);
  art::FindManyP<anab::Calorimetry> trackCaloAssns(trackVec, e, fCaloLabel);
  art::FindManyP<recob::Hit> trackHitAssns(trackVec, e, fTrackLabel);
  art::FindManyP<recob::Slice> pfpSliceAssns(pfpVec, e, fSliceLabel);
  art::FindManyP<anab::T0> sliceT0Assns(sliceVec, e, fOptLabel);

  // If there are no PFParticles then give up and skip the event
  if (pfpVec.empty())
    return;

  // Loop over the PFParticles and find the neutrino
  const size_t neutrinoID(this->GetNeutrinoID(pfpVec));

  // Check that we found a reconstructed neutrino, if not skip the event
  if (neutrinoID == std::numeric_limits<size_t>::max())
    return;

  // Collect the daughter tracks
  const std::vector<art::Ptr<recob::Track>> daughterTracks(this->GetDaughterTracks(neutrinoID, pfpVec, pfpTrackAssns));

  // Check that we found a reconstructed neutrino, if not skip the event
  if (daughterTracks.empty())
    return;

  // Search for the longest daughter track ID
  const int longestID(this->GetLongestTrackID(daughterTracks));
  // Check that we found a longest track
  if (longestID == -1)
    return;

  for (const art::Ptr<recob::Track>& track : daughterTracks) {
    // Add the parameters of the track to the vectors stored in the tree
    fDaughterTrackLengths.push_back(track->Length());
    fDaughterLongestTrack.push_back(track->ID() == longestID);

    // Fill the histogram with the length of this track
    fTrackLengthHist->Fill(track->Length());

    // Get the calorimetry objects associated with the track
    // In this case we expect 3 Calorimetries per track (one for each plane)
    const std::vector<art::Ptr<anab::Calorimetry>> trackCalos(trackCaloAssns.at(track.key()));

    // Iterate through the calorimetry objects and select the one on the collection plane (plane 2)
    for (const art::Ptr<anab::Calorimetry>& calo : trackCalos) {

      // Get the plane number in a simple format
      const int planeNum(calo->PlaneID().Plane);
      // If it is not on the collection plane, skip it
      if (planeNum != 2)
        continue;

      // Store the calorimetry data in the vectors stored in the tree
      // Note here we are pushing back a vector into a vector
      fDaughterTrackdEdx.push_back(calo->dEdx());
      fDaughterTrackResidualRange.push_back(calo->ResidualRange());
    }

    // Get the hits from the track
    const std::vector<art::Ptr<recob::Hit>> trackHits(trackHitAssns.at(track.key()));
    // Use some utility functions to access the backtracker and find the matching particle ID
    const int trackTrueID(TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trackHits, true));

    // Check we found a valid match
    if (TruthMatchUtils::Valid(trackTrueID)) {

      // Get the True particle
      const simb::MCParticle* particle(particleInventoryService->TrackIdToParticle_P(trackTrueID));
      // Save the PDG code of the track
      fDaughterTrackTruePDG.push_back(particle->PdgCode());
    }
  }
  
  // Access the neutrino slice to recover the T0
  const art::Ptr<recob::Slice> neutrinoSlice(this->GetNeutrinoSlice(neutrinoID, pfpVec, pfpSliceAssns));
  
  // Check that we found a reconstructed neutrino slice, if not skip the event
  if (!neutrinoSlice)
    return;

  // Get the T0 object associated with the slice
  const std::vector<art::Ptr<anab::T0>> sliceT0s(sliceT0Assns.at(neutrinoSlice.key()));

  // There should only be 1 T0 per slice
  if (sliceT0s.size() == 1) {
    const art::Ptr<anab::T0>& t0(sliceT0s.front());
    fT0 = t0->Time();
  } // T0s

  // Store the outputs in the TTree
  fTree->Fill();
}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Get the TFileService to create out output histogram for us
  // We know the maximum track length is about 350cm, so lets set 70 bins (5cm each)
  fTrackLengthHist = tfs->make<TH1F>("trackLengthHist", "Reconstructed Track Lengths;Track Length [cm]", 70, 0, 350);

  // Add branches to the TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaries", &fNPrimaries);
  fTree->Branch("nPrimaryDaughters", &fNPrimaryDaughters);

  fTree->Branch("daughterTrackLengths", &fDaughterTrackLengths);
  fTree->Branch("daughterLongestTrack", &fDaughterLongestTrack);
  fTree->Branch("daughterTrackTruePDG", &fDaughterTrackTruePDG);

  fTree->Branch("daughterTrackdEdx", &fDaughterTrackdEdx);
  fTree->Branch("daughterTrackResidualRange", &fDaughterTrackResidualRange);
  
  fTree->Branch("t0", &fT0);
}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

void test::AnalyzeEvents::ResetTree()
{
  // Reset all of our variables to 0 or empty vectors
  // This ensures things are not kept from the previous event
  fNPFParticles      = 0;
  fNPrimaries        = 0;
  fNPrimaryDaughters = 0;
  fT0                = 0.;
  fDaughterTrackLengths.clear();
  fDaughterLongestTrack.clear();
  fDaughterTrackTruePDG.clear();
  fDaughterTrackdEdx.clear();
  fDaughterTrackResidualRange.clear();
}

const int test::AnalyzeEvents::GetNeutrinoID(const std::vector<art::Ptr<recob::PFParticle>>& pfpVec)
{
  // Initialise the neutrino ID
  size_t neutrinoID(std::numeric_limits<size_t>::max());

  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    fNPFParticles++;

    // Check that we are looking at a primary that has a neutrino PDG code, if not move on to the next PFP
    if (!(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12)))
      continue;

    neutrinoID = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;
  }

  return neutrinoID;
}

const std::vector<art::Ptr<recob::Track>> test::AnalyzeEvents::GetDaughterTracks(
    const size_t neutrinoID,
    const std::vector<art::Ptr<recob::PFParticle>>& pfpVec,
    const art::FindManyP<recob::Track>& pfpTrackAssns) const
{
  std::vector<art::Ptr<recob::Track>> daughterTracks;

  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    // We are only interested in the daughter particles
    if (pfp->Parent() != neutrinoID)
      continue;

    // Get the tracks associated to the PFParticle
    const std::vector<art::Ptr<recob::Track>> pfpTracks(pfpTrackAssns.at(pfp.key()));

    // There should only ever be 0 or 1 tracks associated to the PFP
    if (pfpTracks.size() == 1)
      daughterTracks.push_back(pfpTracks.front()); // As we have one entry in the vector, take the first (only) element
  }
  return daughterTracks;
}

const art::Ptr<recob::Slice> test::AnalyzeEvents::GetNeutrinoSlice(
    const size_t neutrinoID,
    const std::vector<art::Ptr<recob::PFParticle>>& pfpVec,
    const art::FindManyP<recob::Slice>& pfpSliceAssns) const
{
  art::Ptr<recob::Slice> pfpSlice;
  // Now access the slices and corresponding timing information
  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    // Start by assessing the neutrino PFParticle itself
    if(pfp->Self() != neutrinoID) continue;

    // Get the slices associated with the current PFParticle
    const std::vector<art::Ptr<recob::Slice>> pfpSlices(pfpSliceAssns.at(pfp.key()));

    // There should only ever be 0 or 1 slices associated to the neutrino PFP
    if (pfpSlices.size() == 1) {
      // Get the first (only) element of the vector
      pfpSlice = pfpSlices.front();
    }
  } // PFParticle vector

  // Will be null if != 1 slice found
  return pfpSlice;
}

const int test::AnalyzeEvents::GetLongestTrackID(const std::vector<art::Ptr<recob::Track>>& trackVec) const
{
  int longestID(-1);
  float longestLength(std::numeric_limits<float>::lowest());

  for (const art::Ptr<recob::Track>& track : trackVec) {
    // If this is the longest track, save the track ID and the new max length
    if (track->Length() > longestLength) {
      longestID = track->ID();
      longestLength = track->Length();
    }
  }
  return longestID;
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
