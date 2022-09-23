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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

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

  std::vector<float> fDaughterTrackLengths;
  std::vector<bool> fDaughterLongestTrack;

  std::vector<std::vector<float>> fDaughterTrackdEdx;
  std::vector<std::vector<float>> fDaughterTrackResidualRange;

  // Define input labels
  const std::string fPFParticleLabel;
  const std::string fTrackLabel;
  const std::string fCaloLabel;
};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
    : EDAnalyzer { p }
    // Initialise out input labels by reading the fhicl parameters
    , fPFParticleLabel(p.get<std::string>("PFParticleLabel"))
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fCaloLabel(p.get<std::string>("CalorimetryLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();

  // Reset all of our variables to 0 or empty vectors
  // This ensures things are not kept from the previous event
  fNPFParticles = 0;
  fNPrimaries = 0;
  fNPrimaryDaughters = 0;
  fDaughterTrackLengths.clear();
  fDaughterLongestTrack.clear();
  fDaughterTrackdEdx.clear();
  fDaughterTrackResidualRange.clear();

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

  // If there are no PFParticles then give up and skip the event
  if (pfpVec.empty())
    return;

  // Initialise neutrino ID
  size_t neutrinoID(std::numeric_limits<size_t>::max());

  // Loop over the PFParticles and find the neutrino
  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    fNPFParticles++;

    // Check that we are looking at a primary that has a neutrino PDG code, if not move on to the next PFP
    if (!(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12)))
      continue;

    neutrinoID = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;
  }

  // Check that we found a reconstructed neutrino, if not skip the event
  if (neutrinoID == std::numeric_limits<size_t>::max())
    return;

  // Load the associations between PFPs, Tracks and Calorimetries
  art::FindManyP<recob::Track> pfpTrackAssns(pfpVec, e, fTrackLabel);
  art::FindManyP<anab::Calorimetry> trackCaloAssns(trackVec, e, fCaloLabel);

  // Search for the longest daughter track ID
  int longestID(-1);
  float longestLength(std::numeric_limits<float>::lowest());

  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    // We are only interested in the daughter particles
    if (pfp->Parent() != neutrinoID)
      continue;

    // Get the tracks associated with the PFP
    const std::vector<art::Ptr<recob::Track>> pfpTracks(pfpTrackAssns.at(pfp.key()));

    // There should only ever be 0 or 1 tracks associated to the PFP
    if (pfpTracks.size() == 1) {
      // Get the first (only) element of the vector
      const art::Ptr<recob::Track>& pfpTrack(pfpTracks.front());

      // If this is the longest track, save the track ID and the new max length
      if (pfpTrack->Length() > longestLength) {
        longestID = pfpTrack->ID();
        longestLength = pfpTrack->Length();
      }
    }
  }

  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    // We are only interested in the daughter particles
    if (pfp->Parent() != neutrinoID)
      continue;

    // Get the tracks associated with the PFP
    const std::vector<art::Ptr<recob::Track>> pfpTracks(pfpTrackAssns.at(pfp.key()));

    // There should only ever be 0 or 1 tracks associated to the PFP
    if (pfpTracks.size() == 1) {
      // Get the first (only) element of the vector
      const art::Ptr<recob::Track>& pfpTrack(pfpTracks.front());

      // Add the parameters of the track to the vectors stored in the tree
      fDaughterTrackLengths.push_back(pfpTrack->Length());
      fDaughterLongestTrack.push_back(pfpTrack->ID() == longestID);

      // Fill the histogram with the length of this track
      fTrackLengthHist->Fill(pfpTrack->Length());

      // Get the calorimetry objects associated with the track
      // In this case we expect 3 Calorimetries per track (one for each plane)
      const std::vector<art::Ptr<anab::Calorimetry>> trackCalos(trackCaloAssns.at(pfpTrack.key()));

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
    }
  }
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

  fTree->Branch("daughterTrackdEdx", &fDaughterTrackdEdx);
  fTree->Branch("daughterTrackResidualRange", &fDaughterTrackResidualRange);
}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
