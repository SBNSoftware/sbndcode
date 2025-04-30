////////////////////////////////////////////////////////////////////////
// Class:       RecoFilter
// Plugin Type: filter (Unknown Unknown)
// File:        RecoFilter_module.cc
//
// Generated at Tue Apr 29 19:39:15 2025 by Anna Beever using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

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

namespace nuclearFragments {
  class RecoFilter;
}


class nuclearFragments::RecoFilter : public art::EDFilter {
public:
  explicit RecoFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoFilter(RecoFilter const&) = delete;
  RecoFilter(RecoFilter&&) = delete;
  RecoFilter& operator=(RecoFilter const&) = delete;
  RecoFilter& operator=(RecoFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  bool fHasRecoDeuteronWithDesiredP;
  int fCountRecoDeuteronWithDesiredP;
  bool fHasNonRecoDeuteronWithDesiredP;
  int fCountNonRecoDeuteronWithDesiredP;
  bool fHasAlphaWithDesiredP;
  std::vector<bool> fMC_isReconstructed;
  std::vector<int> fReco_truthMatchedTrackID;

  // Define input labels
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fNuGenLabel;
  std::string fLArGeantLabel;

  //Functions
  void ResetVariables();
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

};


nuclearFragments::RecoFilter::RecoFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fNuGenLabel(p.get<std::string>("NuGenLabel")),
  fLArGeantLabel(p.get<std::string>("LArGeantLabel"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool nuclearFragments::RecoFilter::filter(art::Event& e)
{ 

  ResetVariables();

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

      // Finding nu slices to loop through
      art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
      art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, e, fTrackLabel);

      std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));
      
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

        const art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
        art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackLabel);

        std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());
        int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
        fReco_truthMatchedTrackID.push_back(trackID);
        }

      }

      break;
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

    std::vector<art::Ptr<simb::MCParticle>> particles = truNuParticleAssoc.at(truNu.key());

    for(const art::Ptr<simb::MCParticle> &particle : particles){
      if(particle->Process() != "primary" || particle->StatusCode() != 1){
        continue;
      }

      int checkIfReconstructed = std::count(fReco_truthMatchedTrackID.begin(), fReco_truthMatchedTrackID.end(), particle->TrackId());
      fMC_isReconstructed.push_back(checkIfReconstructed!=0);

      //filter for clusters

      if(particle->PdgCode() == 1000010020 && particle->P() > 0.8 && particle->P() < 0.9 && checkIfReconstructed!=0){
        fHasRecoDeuteronWithDesiredP = 1;
        fCountRecoDeuteronWithDesiredP++;
      }

      if(particle->PdgCode() == 1000010020 && particle->P() > 0.8 && particle->P() < 0.9 && checkIfReconstructed == 0){
        fHasNonRecoDeuteronWithDesiredP = 1;
        fCountNonRecoDeuteronWithDesiredP++;
      }

      if(particle->PdgCode() == 1000020040 && particle->P() > 1.0){
        fHasAlphaWithDesiredP = 1;
      }

    }
  }

  
  if((fHasNonRecoDeuteronWithDesiredP && fCountNonRecoDeuteronWithDesiredP < 10) || (fHasRecoDeuteronWithDesiredP && fCountRecoDeuteronWithDesiredP < 10) || (fHasAlphaWithDesiredP) )
  {
    std::cout << "event passed filters" << std::endl;
    return true;
  }

  return false;

}

void nuclearFragments::RecoFilter::ResetVariables()
{
  // Set all trackCounters to zero for the current event
  fHasRecoDeuteronWithDesiredP = 0;
  fHasNonRecoDeuteronWithDesiredP = 0;
  fHasAlphaWithDesiredP = 0;
}

DEFINE_ART_MODULE(nuclearFragments::RecoFilter)
