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

// --- art framework ---
// EDFilter: base class for modules that accept or reject events.
//   filter() must return true (keep event) or false (drop event).
// ModuleMacros: provides DEFINE_ART_MODULE() macro to register this class with art.
// Event: the object representing one physics event; used to read and write data products.
// Handle: a smart pointer with validity check used to access data collections in the event.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// FindManyP: reads many-to-many association tables stored in the ROOT file by Pandora.
//   Usage: FindManyP<B> assns(handleA, event, label)
//          assns.at(a.key()) -> vector of B objects associated to object a
// FindOneP: same but for one-to-one associations (not used here, kept for reference).
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// ParameterSet: reads FCL configuration parameters at job startup.
// MessageLogger: logging system (mf::LogInfo, mf::LogWarning, etc.).
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// --- LArSoft data products produced by Pandora ---
// Slice: a group of TPC hits that Pandora believes come from the same physical origin
//   (one neutrino interaction, one cosmic muon, etc.). Each event has several slices.
// PFParticle: a reconstructed particle within a slice, organized in a parent-child tree.
//   The root PFParticle (IsPrimary()==true) with PDG 12 or 14 represents the neutrino hypothesis.
// Vertex: the 3D interaction point associated to a PFParticle.
// PFParticleMetadata: a key-value map attached to each PFParticle by Pandora with scores:
//   "NuScore"       -> float [0,1], probability that this PFParticle is a neutrino
//   "IsClearCosmic" -> presence of this key (any value) flags the PFP as a clear cosmic
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
  // --- FCL parameters (read once at job startup, fixed for all events) ---
  art::InputTag fPandoraLabel;  // label of the Pandora process that produced slices/PFPs
  float         fNuScoreCut;   // minimum NuScore to accept a slice as neutrino candidate
  double        fFVXmin, fFVXmax; // fiducial volume limits in X (cm)
  double        fFVYmin, fFVYmax; // fiducial volume limits in Y (cm)
  double        fFVZmin, fFVZmax; // fiducial volume limits in Z (cm)
  int           fVerbosity;       // 0=silent, 1=per-event summary, 2=per-slice details

  // Returns true if the point (x,y,z) is inside the configured fiducial volume.
  bool InFV(double x, double y, double z) const;
};

// ----------------------------------------------------------------------------
// Constructor: runs once at job startup.
// p.get<T>("Key", default) reads a parameter from the FCL file; if the key is
// absent the default value is used, so the FCL only needs to specify overrides.
// ----------------------------------------------------------------------------
NuSliceFilter::NuSliceFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
  , fPandoraLabel( p.get<art::InputTag>("PandoraLabel", "pandora") )
  , fNuScoreCut(   p.get<float>("NuScoreCut", 0.5) )
  , fVerbosity(    p.get<int>("Verbosity", 0) )
{
  // The FV is given as [min, max] pairs in the FCL because FHiCL has no pair type.
  // We unpack them into 6 separate doubles so InFV() does direct comparisons.
  auto fvx = p.get<std::vector<double>>("FiducialVolumeX", {-180.0, 180.0});
  auto fvy = p.get<std::vector<double>>("FiducialVolumeY", {-185.0, 185.0});
  auto fvz = p.get<std::vector<double>>("FiducialVolumeZ", {10.0, 490.0});
  fFVXmin = fvx.at(0); fFVXmax = fvx.at(1);
  fFVYmin = fvy.at(0); fFVYmax = fvy.at(1);
  fFVZmin = fvz.at(0); fFVZmax = fvz.at(1);

  // Declare that this module will write one int per event named "bestSliceID".
  // art requires this declaration before any call to e.put(); omitting it is a crash.
  // Every code path in filter() must call e.put() exactly once, even when returning false.
  produces<int>("bestSliceID");
}

// ----------------------------------------------------------------------------
// InFV: returns true if (x,y,z) lies strictly inside the fiducial volume.
// The fiducial volume excludes detector edges where reconstruction quality
// degrades and where cosmic muons entering from the side are more likely.
// ----------------------------------------------------------------------------
bool NuSliceFilter::InFV(double x, double y, double z) const {
  return (x > fFVXmin && x < fFVXmax &&
          y > fFVYmin && y < fFVYmax &&
          z > fFVZmin && z < fFVZmax);
}

// ----------------------------------------------------------------------------
// filter(): called once per event. Returns true to keep the event, false to drop it.
// ----------------------------------------------------------------------------
bool NuSliceFilter::filter(art::Event& e) {

  // Start with "nothing found". If any early return fires, -1 is written to the event.
  int  bestSliceID = -1;
  bool foundSlice  = false;

  // --- Load the list of Pandora slices from the event ---
  // A slice is a group of hits that Pandora reconstructed as coming from the same
  // physical source. Typically one event has several slices: one neutrino + several cosmics.
  // Handle is a smart pointer; isValid() returns false if Pandora did not run or failed.
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fPandoraLabel, sliceHandle);
  if (!sliceHandle.isValid()) {
    mf::LogWarning("NuSliceFilter") << "No slices found with label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  // --- Load the list of Pandora PFParticles from the event ---
  // Each slice contains several PFParticles: one root (the neutrino or cosmic hypothesis)
  // and its children (muon, proton, shower, etc.). All have the same Pandora label.
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPandoraLabel, pfpHandle);
  if (!pfpHandle.isValid()) {
    mf::LogWarning("NuSliceFilter") << "No PFParticles found with label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  // --- Build association lookup tables ---
  // Pandora stores connections between objects as separate association tables in the ROOT file.
  // FindManyP loads these tables so we can query them by object index (key()).
  //   slice_pfp_assns : given a Slice, get all PFParticles that belong to it
  //   pfp_vtx_assns   : given a PFParticle, get its reconstructed 3D vertex
  //   pfp_meta_assns  : given a PFParticle, get its metadata (NuScore, IsClearCosmic, ...)
  art::FindManyP<recob::PFParticle>                 slice_pfp_assns(sliceHandle, e, fPandoraLabel);
  art::FindManyP<recob::Vertex>                     pfp_vtx_assns(pfpHandle,   e, fPandoraLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_meta_assns(pfpHandle,  e, fPandoraLabel);

  if (!slice_pfp_assns.isValid() || !pfp_vtx_assns.isValid() || !pfp_meta_assns.isValid()) {
    mf::LogWarning("NuSliceFilter") << "Invalid associations for pandora label " << fPandoraLabel;
    e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
    return false;
  }

  // Convert the Handle into a vector of art::Ptr so we can use .key() to index the associations.
  std::vector<art::Ptr<recob::Slice>> sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);

  // Track the best neutrino candidate across all slices in this event.
  float bestNuScore = -999.f;

  // --- Loop over all slices in the event ---
  for (auto const& slice : sliceVect) {

    // Get all PFParticles that belong to this slice.
    std::vector<art::Ptr<recob::PFParticle>> pfpVec = slice_pfp_assns.at(slice.key());

    // Per-slice flags and values, reset for each slice.
    bool   isClearCosmic = false; // true if ANY PFP in this slice is flagged as cosmic
    float  nuScore       = -999.f;
    bool   hasNuPFP      = false; // true if we found a primary nu PFP with NuScore
    bool   hasVertex     = false; // true if the nu PFP has a reconstructed vertex
    double vtxX = 0, vtxY = 0, vtxZ = 0;

    // --- Loop over all PFParticles in this slice ---
    for (auto const& pfp : pfpVec) {

      // Get the metadata map for this PFParticle.
      // PropertiesMap is a std::map<std::string, float> produced by Pandora.
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metaVec =
        pfp_meta_assns.at(pfp.key());

      for (auto const& meta : metaVec) {
        const larpandoraobj::PFParticleMetadata::PropertiesMap& props =
          meta->GetPropertiesMap();

        // IsClearCosmic: Pandora sets this key on ANY PFP in a slice it is confident
        // is a cosmic. If any PFP in the slice is cosmic, we reject the whole slice.
        if (props.count("IsClearCosmic")) { isClearCosmic = true; }

        // NuScore and vertex are stored only on the primary PFP with PDG 12 (nu_e)
        // or 14 (nu_mu). This is the root of the neutrino hierarchy in Pandora.
        if (pfp->IsPrimary() &&
            (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14)) {
          if (props.count("NuScore")) {
            nuScore  = props.at("NuScore");
            hasNuPFP = true;
          }
          // The vertex is stored as a separate data product associated to the primary PFP.
          // vtxVec should have exactly one element for a well-formed Pandora output.
          const std::vector<art::Ptr<recob::Vertex>> vtxVec = pfp_vtx_assns.at(pfp.key());
          if (!vtxVec.empty()) {
            const geo::Point_t pos = vtxVec.front()->position();
            vtxX = pos.X(); vtxY = pos.Y(); vtxZ = pos.Z();
            hasVertex = true;
          }
        }
      }

      // Once a cosmic PFP is found there is no need to inspect the rest of the slice.
      if (isClearCosmic) break;
    }

    // --- Apply selection cuts to this slice ---
    // Each cut is a separate guard; failing any one skips to the next slice.

    // Cut 1: reject if Pandora flagged any PFP in this slice as a clear cosmic.
    if (isClearCosmic) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: IsClearCosmic";
      continue;
    }
    // Cut 2: reject if there is no primary nu PFP (slice has no neutrino hypothesis).
    if (!hasNuPFP)                  continue;
    // Cut 3: reject if NuScore is below the configured threshold.
    if (nuScore < fNuScoreCut) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: NuScore=" << nuScore;
      continue;
    }
    // Cut 4: reject if there is no reconstructed vertex (degenerate Pandora output).
    if (!hasVertex)                 continue;
    // Cut 5: reject if the vertex is outside the fiducial volume.
    //   Vertices near the edges have worse reconstruction and higher cosmic contamination.
    if (!InFV(vtxX, vtxY, vtxZ)) {
      if (fVerbosity > 1) mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " REJECTED: vertex outside FV";
      continue;
    }

    if (fVerbosity > 0)
      mf::LogInfo("NuSliceFilter") << "  Slice " << slice->ID() << " PASSED: NuScore=" << nuScore
                                   << " vtx=(" << vtxX << "," << vtxY << "," << vtxZ << ")";

    // Keep only the slice with the highest NuScore among all passing candidates.
    // In practice most events have at most one neutrino slice, but we handle the
    // general case to avoid ambiguity.
    if (nuScore > bestNuScore) {
      bestNuScore  = nuScore;
      bestSliceID  = slice->ID();
      foundSlice   = true;
    }
  }

  if (fVerbosity > 0)
    mf::LogInfo("NuSliceFilter") << "Event " << e.event() << ": "
                                 << (foundSlice ? "PASSED sliceID=" + std::to_string(bestSliceID) : "REJECTED");

  // Write the best slice ID into the event so downstream modules (NuSliceAnalyzer)
  // can read it without repeating the selection logic. Value is -1 if event is rejected.
  // This e.put() must always be called because we declared produces<int>("bestSliceID").
  e.put(std::make_unique<int>(bestSliceID), "bestSliceID");
  return foundSlice;
}

DEFINE_ART_MODULE(NuSliceFilter)

} // namespace opdet
