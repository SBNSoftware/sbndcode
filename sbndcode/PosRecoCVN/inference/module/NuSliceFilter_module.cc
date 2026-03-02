// ============================================================================
// NuSliceFilter_module.cc
//
// art EDFilter that selects events containing at least one high-purity
// neutrino slice from Pandora reconstruction.
//
// Cuts applied:
//   - is_clear_cosmic == 0  (no PFP in the slice has IsClearCosmic flag)
//   - nu_score > NuScoreCut (Pandora NuScore from PFParticleMetadata)
//   - Pandora vertex inside Fiducial Volume
//
// Produces:
//   - int "bestSliceID": ID of the best neutrino slice (-1 if event is rejected)
//     Downstream modules (e.g. NuSliceAnalyzer) read this to avoid repeating
//     the selection logic.
//
// Events passing the cuts continue to the ResNet (PosRecoCVNProducer).
// Events not passing are dropped, so the ResNet only runs on neutrino-rich events.
// ============================================================================

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// C++
#include <memory>
#include <vector>
#include <cmath>
#include <string>

namespace opdet {

class NuSliceFilter : public art::EDFilter {
public:
  explicit NuSliceFilter(fhicl::ParameterSet const& p);
  bool filter(art::Event& e) override;

private:
  // FCL parameters
  art::InputTag fPandoraLabel;
  float         fNuScoreCut;
  double        fFVXmin, fFVXmax;
  double        fFVYmin, fFVYmax;
  double        fFVZmin, fFVZmax;
  int           fVerbosity;

  bool InFV(double x, double y, double z) const;
};

// ----------------------------------------------------------------------------
NuSliceFilter::NuSliceFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
  , fPandoraLabel( p.get<art::InputTag>("PandoraLabel", "pandora") )
  , fNuScoreCut(   p.get<float>("NuScoreCut", 0.5) )
  , fVerbosity(    p.get<int>("Verbosity", 0) )
{
  auto fvx = p.get<std::vector<double>>("FiducialVolumeX", {-180.0, 180.0});
  auto fvy = p.get<std::vector<double>>("FiducialVolumeY", {-185.0, 185.0});
  auto fvz = p.get<std::vector<double>>("FiducialVolumeZ", {10.0, 490.0});
  fFVXmin = fvx.at(0); fFVXmax = fvx.at(1);
  fFVYmin = fvy.at(0); fFVYmax = fvy.at(1);
  fFVZmin = fvz.at(0); fFVZmax = fvz.at(1);

  produces<int>("bestSliceID");
}

// ----------------------------------------------------------------------------
bool NuSliceFilter::InFV(double x, double y, double z) const {
  return (x > fFVXmin && x < fFVXmax &&
          y > fFVYmin && y < fFVYmax &&
          z > fFVZmin && z < fFVZmax);
}

// ----------------------------------------------------------------------------
bool NuSliceFilter::filter(art::Event& e) {

  // Default: no valid slice found
  int  bestSliceID = -1;
  bool foundSlice  = false;

  // --- Load Pandora Slices
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fPandoraLabel, sliceHandle);
  if (!sliceHandle.isValid()) {
    mf::LogWarning("NuSliceFilter") << "No slices found with label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  // --- Load Pandora PFParticles
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPandoraLabel, pfpHandle);
  if (!pfpHandle.isValid()) {
    mf::LogWarning("NuSliceFilter") << "No PFParticles found with label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  // --- Associations
  art::FindManyP<recob::PFParticle>                 slice_pfp_assns(sliceHandle, e, fPandoraLabel);
  art::FindManyP<recob::Vertex>                     pfp_vtx_assns(pfpHandle,   e, fPandoraLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_meta_assns(pfpHandle,  e, fPandoraLabel);

  if (!slice_pfp_assns.isValid() || !pfp_vtx_assns.isValid() || !pfp_meta_assns.isValid()) {
    mf::LogWarning("NuSliceFilter") << "Invalid associations for pandora label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  std::vector<art::Ptr<recob::Slice>> sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);

  float bestNuScore = -999.f;

  for (auto const& slice : sliceVect) {

    std::vector<art::Ptr<recob::PFParticle>> pfpVec = slice_pfp_assns.at(slice.key());

    bool   isClearCosmic = false;
    float  nuScore       = -999.f;
    bool   hasNuPFP      = false;
    bool   hasVertex     = false;
    double vtxX = 0, vtxY = 0, vtxZ = 0;

    for (auto const& pfp : pfpVec) {

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metaVec =
        pfp_meta_assns.at(pfp.key());

      for (auto const& meta : metaVec) {
        const larpandoraobj::PFParticleMetadata::PropertiesMap& props =
          meta->GetPropertiesMap();

        if (props.count("IsClearCosmic")) { isClearCosmic = true; }

        if (pfp->IsPrimary() &&
            (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14)) {
          if (props.count("NuScore")) {
            nuScore  = props.at("NuScore");
            hasNuPFP = true;
          }
          const std::vector<art::Ptr<recob::Vertex>> vtxVec = pfp_vtx_assns.at(pfp.key());
          if (!vtxVec.empty()) {
            const geo::Point_t pos = vtxVec.front()->position();
            vtxX = pos.X(); vtxY = pos.Y(); vtxZ = pos.Z();
            hasVertex = true;
          }
        }
      }

      if (isClearCosmic) break;
    }

    // --- Apply cuts
    if (isClearCosmic) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: IsClearCosmic";
      continue;
    }
    if (!hasNuPFP)                  continue;
    if (nuScore < fNuScoreCut) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: NuScore=" << nuScore;
      continue;
    }
    if (!hasVertex)                 continue;
    if (!InFV(vtxX, vtxY, vtxZ)) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: vertex outside FV";
      continue;
    }

    if (fVerbosity > 0)
      mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " PASSED: NuScore=" << nuScore
                                   << " vtx=(" << vtxX << "," << vtxY << "," << vtxZ << ")";

    if (nuScore > bestNuScore) {
      bestNuScore  = nuScore;
      bestSliceID  = slice->ID();
      foundSlice   = true;
    }
  }

  if (fVerbosity > 0)
    mf::LogInfo("NuSliceFilter") << "Event " << e.event() << ": "
                                 << (foundSlice ? "PASSED sliceID=" + std::to_string(bestSliceID) : "REJECTED");

  e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
  return foundSlice;
}

DEFINE_ART_MODULE(NuSliceFilter)

} // namespace opdet
