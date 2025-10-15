///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//   ██████╗░░█████╗░██████╗░██╗░░░██╗██╗░██████╗  ███████╗██╗██████╗░░██████╗████████╗  //
//   ██╔══██╗██╔══██╗██╔══██╗╚██╗░██╔╝╚█║██╔════╝  ██╔════╝██║██╔══██╗██╔════╝╚══██╔══╝  //
//   ██████╦╝███████║██████╦╝░╚████╔╝░░╚╝╚█████╗░  █████╗░░██║██████╔╝╚█████╗░░░░██║░░░  //
//   ██╔══██╗██╔══██║██╔══██╗░░╚██╔╝░░░░░░╚═══██╗  ██╔══╝░░██║██╔══██╗░╚═══██╗░░░██║░░░  //
//   ██████╦╝██║░░██║██████╦╝░░░██║░░░░░░██████╔╝  ██║░░░░░██║██║░░██║██████╔╝░░░██║░░░  //
//   ╚═════╝░╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░░░░╚═════╝░  ╚═╝░░░░░╚═╝╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░  //
//                                                                                       //
//   ░█████╗░███╗░░██╗░█████╗░██╗░░░░░██╗░░░██╗███████╗███████╗██████╗░                  //
//   ██╔══██╗████╗░██║██╔══██╗██║░░░░░╚██╗░██╔╝╚════██║██╔════╝██╔══██╗                  //
//   ███████║██╔██╗██║███████║██║░░░░░░╚████╔╝░░░███╔═╝█████╗░░██████╔╝                  //
//   ██╔══██║██║╚████║██╔══██║██║░░░░░░░╚██╔╝░░██╔══╝░░██╔══╝░░██╔══██╗                  //
//   ██║░░██║██║░╚███║██║░░██║███████╗░░░██║░░░███████╗███████╗██║░░██║                  //
//   ╚═╝░░╚═╝╚═╝░░╚══╝╚═╝░░╚═╝╚══════╝░░░╚═╝░░░╚══════╝╚══════╝╚═╝░░╚═╝                  //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

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
#include "lardataobj/AnalysisBase/Calorimetry.h"

// SBN(D) includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

// Additional framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

// ROOT includes
#include <TTree.h>

namespace test {
  class AnalyseEvents;
}


class test::AnalyseEvents : public art::EDAnalyzer {
public:
  explicit AnalyseEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseEvents(AnalyseEvents const&) = delete;
  AnalyseEvents(AnalyseEvents&&) = delete;
  AnalyseEvents& operator=(AnalyseEvents const&) = delete;
  AnalyseEvents& operator=(AnalyseEvents&&) = delete;

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
  double fOpT0;

  std::vector<float> fChildTrackLengths;
  std::vector<bool> fChildTrackIsLongest;
  std::vector<std::vector<float>> fChildTrackdEdx;
  std::vector<std::vector<float>> fChildTrackResRange;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fOpT0FinderLabel;
};

test::AnalyseEvents::AnalyseEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fOpT0FinderLabel(p.get<std::string>("OpT0FinderLabel"))
  {
    // Call appropriate consumes<>() for any products to be retrieved by this module.
  }

void test::AnalyseEvents::analyze(art::Event const& e)
{
  // Set the event ID
  fEventID = e.id().event();

  // Prepare variables for new event (reset counters to 0 / set default values / empty vectors)
  fNPFParticles = 0;
  fNPrimaryChildren = 0;
  fOpT0 = std::numeric_limits<double>::lowest();
  fChildTrackLengths.clear();
  fChildTrackIsLongest.clear();
  fChildTrackdEdx.clear();
  fChildTrackResRange.clear();

  // Get event slices
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if (sliceHandle.isValid())
    art::fill_ptr_vector(sliceVector, sliceHandle);

  // Get associations between slices and pfparticles & opt0 results
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);
  art::FindManyP<sbn::OpT0Finder> sliceOpT0Assoc(sliceHandle, e, fOpT0FinderLabel);

  // Filling our neutrino hierarchy variables
  int nuID = -1;
  int nuSliceKey = -1;

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

          // Get any OpT0Finder results associated with our slice
          std::vector<art::Ptr<sbn::OpT0Finder>> opT0s = sliceOpT0Assoc.at(nuSliceKey);

          // Occasionally there may be multiple results, let's use the one with the best score
          std::sort(opT0s.begin(), opT0s.end(),
                    [](auto const& a, auto const& b)
                    { return a->score > b->score; });

          // The best score will now be at the front of the vector (if there were any)
          if (opT0s.size() != 0)
            fOpT0 = opT0s[0]->time;

          break;
        }

      if (nuID >= 0)
        break;
    }

  if(nuID < 0)
    return;

  // Now let's look at our tracks
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle =
    e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  art::ValidHandle<std::vector<recob::Track>> trackHandle =
    e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);

  art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);
  art::FindManyP<anab::Calorimetry> trackCaloAssoc(trackHandle, e, fCalorimetryLabel);

  std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

  // Let's find the longest track before we progress with filling the track variables
  int longestID = std::numeric_limits<int>::lowest();
  float longestLength = std::numeric_limits<float>::lowest();

  for(const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
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

      // Check if this track is longer than the current longest
      if(track->Length() > longestLength)
        {
          // If yes, then overwrite the variables to reflect the new longest track
          longestID = track->ID();
          longestLength = track->Length();
        }
    }

  // Now loop through the PFPs again to fill the track variables for the tree
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
      // Was this track the one we found to be the longest earlier?
      fChildTrackIsLongest.push_back(track->ID() == longestID);

      // Get the calorimetry object
      std::vector<art::Ptr<anab::Calorimetry>> calos = trackCaloAssoc.at(track.key());

      for(auto const& calo : calos)
        {
          const int plane = calo->PlaneID().Plane;

          // Only interested in the collection plane (2)
          if(plane != 2)
            continue;

          fChildTrackdEdx.push_back(calo->dEdx());
          fChildTrackResRange.push_back(calo->ResidualRange());
        }
    }

  // Fill tree
  fTree->Fill();
}

void test::AnalyseEvents::beginJob()
{
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  // Add branches to TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaryChildren", &fNPrimaryChildren);
  fTree->Branch("opT0", &fOpT0);
  fTree->Branch("childTrackLengths", &fChildTrackLengths);
  fTree->Branch("childTrackIsLongest", &fChildTrackIsLongest);
  fTree->Branch("childTrackdEdx", &fChildTrackdEdx);
  fTree->Branch("childTrackResRange", &fChildTrackResRange);
}

void test::AnalyseEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyseEvents)
