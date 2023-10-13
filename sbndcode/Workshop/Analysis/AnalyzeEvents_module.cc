////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Thu Sep 28 08:45:26 2023 by Isobel Mawby using cetskelgen
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

// Additional LArSoft includes
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Additional framework includes
#include "art_root_io/TFileService.h"

// ROOT includes
#include <TTree.h>
#include <TH1F.h>

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
  // Create output TTree
  TTree *fTree;

  // Create output histogram
  TH1F *fTrackLengthHist;

  // Tree variables
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaryChildren;

  std::vector<float> fChildTrackLengths;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    fSliceLabel(p.get<std::string>("SliceLabel")),
    fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
    fTrackLabel(p.get<std::string>("TrackLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Set the event ID
  fEventID = e.id().event();

  // Set all counters to 0 for the current event 
  fNPFParticles = 0;
  fNPrimaryChildren = 0;
  fChildTrackLengths.clear();

  // Get event slices
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if (sliceHandle.isValid())
      art::fill_ptr_vector(sliceVector, sliceHandle);

  // Get associations between slices and pfparticles
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

  // Filling our neutrino hierarchy variables
  int nuID = -1;
  int nuSliceKey = -1;

  std::cout << "NEW EVENT" << std::endl;

  for (const art::Ptr<recob::Slice> &slice : sliceVector)
  {
      std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

      for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs)
      {
          const bool isPrimary(slicePFP->IsPrimary());
          const bool isNeutrino((std::abs(slicePFP->PdgCode()) == 12) || (std::abs(slicePFP->PdgCode()) == 14));

          if (!(isPrimary && isNeutrino))
              continue;

          // We have found our neutrino!
          nuSliceKey = slice.key();
          nuID = slicePFP->Self();
          fNPFParticles = slicePFPs.size();
          fNPrimaryChildren = slicePFP->NumDaughters();

          break;
      }

      if (nuID >= 0)
          break;
  }

  // Now let's look at our tracks
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);

  std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

  for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
  {
      // We are only interested in neutrino children particles
      if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
          continue;

      // Get tracks associated with this PFParticle
      std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

      // There should only be 0 or 1 tracks associated with a PFP
      if (tracks.size() != 1)
          continue;

      // Get the track
      art::Ptr<recob::Track> track = tracks.at(0);

      // Add parameters from the track to the branch vector
      fChildTrackLengths.push_back(track->Length());

      // Fill the track length histogram with this entry
      fTrackLengthHist->Fill(track->Length());
  }

  // Fill tree
  fTree->Fill();
}

void test::AnalyzeEvents::beginJob()
{
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");
  fTrackLengthHist = tfs->make<TH1F>("trackLengthHist", "Reconstructed Track Lengths; Track Length (cm);N Tracks", 20, 0, 350);

  // Add branches to TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaryChildren", &fNPrimaryChildren);
  fTree->Branch("childTrackLengths", &fChildTrackLengths);
}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
