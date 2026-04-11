// ============================================================================
// NuSliceAnalyzer_module.cc
//
// art EDAnalyzer that runs on events passing NuSliceFilter.
// For each event it:
//   1. Reads the bestSliceID produced by NuSliceFilter (no cut re-evaluation)
//   2. Extracts the Pandora vertex and NuScore of the selected slice
//   3. Computes SpacePoint barycenter and PCA (charge-weighted by collection-plane hit integral)
//   4. Reads the ResNet (PosRecoCVNProducer) predictions from PixelMapVars data product
//   5. Fills a TTree with all the above
//
// FCL parameters:
//   PandoraLabel     : input tag for pandora products (default: "pandora")
//   FilterLabel      : input tag for NuSliceFilter (default: "nuselector")
//   ResNetLabel      : input tag for PixelMapVars data product (default: "ResNetInference")
//   Verbosity        : 0=minimal, 1=info, 2=debug
// ============================================================================

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// PosRecoCVN data product
#include "sbndcode/NuIntNNOpReco/3-inference-larsoft-module/module/PixelMapVars.h"

// ROOT
#include "TTree.h"

// Eigen for PCA (same as PosRecoCVNProducer_module)
#include <Eigen/Dense>

// C++
#include <vector>
#include <cmath>
#include <string>

namespace opdet {

class NuSliceAnalyzer : public art::EDAnalyzer {
public:
  explicit NuSliceAnalyzer(fhicl::ParameterSet const& p);
  void beginJob() override;
  void analyze(art::Event const& e) override;

private:
  // FCL parameters
  art::InputTag fPandoraLabel;
  art::InputTag fFilterLabel;
  art::InputTag fResNetLabel;
  int           fVerbosity;

  // TTree
  TTree* fTree;

  // Tree variables
  int    fRun, fSubrun, fEvent;
  int    fSliceID;
  float  fSliceNuScore;
  // Pandora vertex of the selected slice
  double fVtxX, fVtxY, fVtxZ;
  // SpacePoint charge barycenter (collection-plane hit integral weighted)
  double fSpBaryX, fSpBaryY, fSpBaryZ;
  // SpacePoint PCA first principal axis (unit vector, pcaZ >= 0)
  double fSpPCAx, fSpPCAy, fSpPCAz;
  // PCA eigenvalues (λ1≥λ2≥λ3, charge-weighted covariance; λ1 eigenvector = spPCAx/y/z)
  double fSpPCALam1, fSpPCALam2, fSpPCALam3;
  // Secondary (v2) and tertiary (v3) PCA eigenvectors
  double fSpPCAv2x, fSpPCAv2y, fSpPCAv2z;
  double fSpPCAv3x, fSpPCAv3y, fSpPCAv3z;
  int    fNSpacePoints;
  // ResNet CNN predictions (from PixelMapVars data product)
  double fCnnPredX, fCnnPredY, fCnnPredZ;

  void ResetVars();
};

// ----------------------------------------------------------------------------
NuSliceAnalyzer::NuSliceAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fPandoraLabel( p.get<art::InputTag>("PandoraLabel", "pandora") )
  , fFilterLabel(  p.get<art::InputTag>("FilterLabel",  "nuselector") )
  , fResNetLabel(  p.get<art::InputTag>("ResNetLabel",  "ResNetInference") )
  , fVerbosity(    p.get<int>("Verbosity", 0) )
  , fTree(nullptr)
{}

// ----------------------------------------------------------------------------
void NuSliceAnalyzer::ResetVars() {
  fRun = fSubrun = fEvent = 0;
  fSliceID      = -1;
  fSliceNuScore = -999.f;
  fVtxX = fVtxY = fVtxZ         = -999.0;
  fSpBaryX = fSpBaryY = fSpBaryZ = -999.0;
  fSpPCAx  = fSpPCAy  = fSpPCAz  = -999.0;
  fSpPCALam1 = fSpPCALam2 = fSpPCALam3 = -999.0;
  fSpPCAv2x = fSpPCAv2y = fSpPCAv2z = -999.0;
  fSpPCAv3x = fSpPCAv3y = fSpPCAv3z = -999.0;
  fNSpacePoints = 0;
  fCnnPredX = fCnnPredY = fCnnPredZ = -999.0;
}

// ----------------------------------------------------------------------------
void NuSliceAnalyzer::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("nuslicetree", "NuSlice TPC + CNN analysis tree");

  // Event ID
  fTree->Branch("run",    &fRun,    "run/I");
  fTree->Branch("subrun", &fSubrun, "subrun/I");
  fTree->Branch("event",  &fEvent,  "event/I");

  // Selected slice info
  fTree->Branch("sliceID",       &fSliceID,       "sliceID/I");
  fTree->Branch("sliceNuScore",  &fSliceNuScore,  "sliceNuScore/F");

  // Pandora vertex
  fTree->Branch("vtxX", &fVtxX, "vtxX/D");
  fTree->Branch("vtxY", &fVtxY, "vtxY/D");
  fTree->Branch("vtxZ", &fVtxZ, "vtxZ/D");

  // SpacePoint charge barycenter
  fTree->Branch("spBaryX", &fSpBaryX, "spBaryX/D");
  fTree->Branch("spBaryY", &fSpBaryY, "spBaryY/D");
  fTree->Branch("spBaryZ", &fSpBaryZ, "spBaryZ/D");

  // SpacePoint PCA first axis
  fTree->Branch("spPCAx", &fSpPCAx, "spPCAx/D");
  fTree->Branch("spPCAy", &fSpPCAy, "spPCAy/D");
  fTree->Branch("spPCAz", &fSpPCAz, "spPCAz/D");
  // PCA eigenvalues (λ1≥λ2≥λ3, charge-weighted covariance; λ1 eigenvector = spPCAx/y/z)
  fTree->Branch("spPCALam1", &fSpPCALam1, "spPCALam1/D");
  fTree->Branch("spPCALam2", &fSpPCALam2, "spPCALam2/D");
  fTree->Branch("spPCALam3", &fSpPCALam3, "spPCALam3/D");
  // Secondary (v2) and tertiary (v3) PCA eigenvectors
  fTree->Branch("spPCAv2x",  &fSpPCAv2x,  "spPCAv2x/D");
  fTree->Branch("spPCAv2y",  &fSpPCAv2y,  "spPCAv2y/D");
  fTree->Branch("spPCAv2z",  &fSpPCAv2z,  "spPCAv2z/D");
  fTree->Branch("spPCAv3x",  &fSpPCAv3x,  "spPCAv3x/D");
  fTree->Branch("spPCAv3y",  &fSpPCAv3y,  "spPCAv3y/D");
  fTree->Branch("spPCAv3z",  &fSpPCAv3z,  "spPCAv3z/D");

  fTree->Branch("nSpacePoints", &fNSpacePoints, "nSpacePoints/I");

  // ResNet CNN predictions
  fTree->Branch("cnnPredX", &fCnnPredX, "cnnPredX/D");
  fTree->Branch("cnnPredY", &fCnnPredY, "cnnPredY/D");
  fTree->Branch("cnnPredZ", &fCnnPredZ, "cnnPredZ/D");
}

// ----------------------------------------------------------------------------
void NuSliceAnalyzer::analyze(art::Event const& e) {
  ResetVars();

  fRun    = e.run();
  fSubrun = e.subRun();
  fEvent  = e.event();

  // =========================================================================
  // 1. Read ResNet predictions from PixelMapVars data product
  // =========================================================================
  art::Handle<PixelMapVars> pmvHandle;
  e.getByLabel(fResNetLabel, pmvHandle);
  if (pmvHandle.isValid() && pmvHandle->passed_filters) {
    if (!pmvHandle->dEpromx_pred.empty()) fCnnPredX = pmvHandle->dEpromx_pred.at(0);
    if (!pmvHandle->dEpromy_pred.empty()) fCnnPredY = pmvHandle->dEpromy_pred.at(0);
    if (!pmvHandle->dEpromz_pred.empty()) fCnnPredZ = pmvHandle->dEpromz_pred.at(0);
    if (fVerbosity > 0)
      mf::LogInfo("NuSliceAnalyzer") << "  ResNet prediction: (" << fCnnPredX << ", " << fCnnPredY << ", " << fCnnPredZ << ")";
  } else if (fVerbosity > 0) {
    mf::LogInfo("NuSliceAnalyzer") << "  PixelMapVars not available or did not pass ResNet filters";
  }

  // =========================================================================
  // 2. Read bestSliceID produced by NuSliceFilter (no cut re-evaluation)
  // =========================================================================
  art::Handle<int> sidHandle;
  e.getByLabel(art::InputTag(fFilterLabel.label(), "bestSliceID"), sidHandle);
  if (!sidHandle.isValid() || *sidHandle < 0) {
    if (fVerbosity > 0) mf::LogInfo("NuSliceAnalyzer") << "  No valid neutrino slice (bestSliceID not found or -1)";
    fTree->Fill();
    return;
  }
  fSliceID = *sidHandle;

  // =========================================================================
  // 3. Load Pandora products and find the selected slice by ID
  // =========================================================================
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fPandoraLabel, sliceHandle);
  if (!sliceHandle.isValid()) {
    mf::LogWarning("NuSliceAnalyzer") << "No slices found";
    fTree->Fill();
    return;
  }

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPandoraLabel, pfpHandle);
  if (!pfpHandle.isValid()) {
    mf::LogWarning("NuSliceAnalyzer") << "No PFParticles found";
    fTree->Fill();
    return;
  }

  art::FindManyP<recob::PFParticle>                 slice_pfp_assns(sliceHandle, e, fPandoraLabel);
  art::FindManyP<recob::Vertex>                     pfp_vtx_assns(pfpHandle,   e, fPandoraLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_meta_assns(pfpHandle,  e, fPandoraLabel);
  art::FindManyP<recob::SpacePoint>                 pfp_sp_assns(pfpHandle,    e, fPandoraLabel);

  art::Handle<std::vector<recob::SpacePoint>> spHandle;
  e.getByLabel(fPandoraLabel, spHandle);
  art::FindManyP<recob::Hit> sp_hit_assns(spHandle, e, fPandoraLabel);

  std::vector<art::Ptr<recob::Slice>> sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);

  // Find the slice whose ID matches bestSliceID
  art::Ptr<recob::Slice> bestSlice;
  for (auto const& slice : sliceVect) {
    if (slice->ID() == fSliceID) { bestSlice = slice; break; }
  }

  if (bestSlice.isNull()) {
    mf::LogWarning("NuSliceAnalyzer") << "Slice ID=" << fSliceID << " not found in event";
    fTree->Fill();
    return;
  }

  // =========================================================================
  // 4. Extract NuScore and vertex from the primary neutrino PFP
  // =========================================================================
  std::vector<art::Ptr<recob::PFParticle>> pfpVec = slice_pfp_assns.at(bestSlice.key());

  for (auto const& pfp : pfpVec) {
    if (!pfp->IsPrimary() ||
        (std::abs(pfp->PdgCode()) != 12 && std::abs(pfp->PdgCode()) != 14)) continue;

    const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metaVec =
      pfp_meta_assns.at(pfp.key());
    for (auto const& meta : metaVec) {
      const larpandoraobj::PFParticleMetadata::PropertiesMap& props = meta->GetPropertiesMap();
      if (props.count("NuScore")) fSliceNuScore = props.at("NuScore");
    }

    const std::vector<art::Ptr<recob::Vertex>> vtxVec = pfp_vtx_assns.at(pfp.key());
    if (!vtxVec.empty()) {
      const geo::Point_t pos = vtxVec.front()->position();
      fVtxX = pos.X(); fVtxY = pos.Y(); fVtxZ = pos.Z();
    }
    break; // One primary neutrino PFP per slice
  }

  if (fVerbosity > 0)
    mf::LogInfo("NuSliceAnalyzer") << "  Best slice ID=" << fSliceID
                                   << " NuScore=" << fSliceNuScore
                                   << " vtx=(" << fVtxX << "," << fVtxY << "," << fVtxZ << ")";

  // =========================================================================
  // 5. Collect SpacePoints from all PFPs in the best slice
  // =========================================================================
  std::vector<double> spX, spY, spZ, spW;

  for (auto const& pfp : pfpVec) {
    const std::vector<art::Ptr<recob::SpacePoint>> spVec = pfp_sp_assns.at(pfp.key());
    for (auto const& sp : spVec) {
      // Use collection-plane (plane 2) hit integral as charge weight
      double integral = 1.0;
      if (sp_hit_assns.isValid()) {
        const std::vector<art::Ptr<recob::Hit>> hits = sp_hit_assns.at(sp.key());
        for (auto const& h : hits) {
          if (h->WireID().Plane == 2) { integral = h->Integral(); break; }
        }
      }
      spX.push_back(sp->position().X());
      spY.push_back(sp->position().Y());
      spZ.push_back(sp->position().Z());
      spW.push_back(integral);
    }
  }

  fNSpacePoints = (int)spX.size();
  if (fVerbosity > 0) mf::LogInfo("NuSliceAnalyzer") << "  SpacePoints: " << fNSpacePoints;

  // =========================================================================
  // 6. Compute charge barycenter and PCA
  // =========================================================================
  if (!spX.empty()) {
    double totalW = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
    for (size_t i = 0; i < spX.size(); ++i) {
      cx += spW[i] * spX[i]; cy += spW[i] * spY[i]; cz += spW[i] * spZ[i];
      totalW += spW[i];
    }

    if (totalW > 0.0) {
      cx /= totalW; cy /= totalW; cz /= totalW;
      fSpBaryX = cx; fSpBaryY = cy; fSpBaryZ = cz;

      if (fVerbosity > 0)
        mf::LogInfo("NuSliceAnalyzer") << "  Charge barycenter: (" << fSpBaryX << ", " << fSpBaryY << ", " << fSpBaryZ << ")";

      if (spX.size() >= 2) {
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < spX.size(); ++i) {
          Eigen::Vector3d pt(spX[i] - cx, spY[i] - cy, spZ[i] - cz);
          cov += std::abs(spW[i]) * (pt * pt.transpose());
        }
        cov /= totalW;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
        Eigen::Vector3d eigvals = solver.eigenvalues();
        Eigen::Matrix3d eigvecs = solver.eigenvectors();
        int maxIdx;
        eigvals.maxCoeff(&maxIdx);
        Eigen::Vector3d dir = eigvecs.col(maxIdx);
        // Eigen sorts ascending: col(2)=λ1 (dominant), col(1)=λ2, col(0)=λ3
        fSpPCALam1 = eigvals(2); fSpPCALam2 = eigvals(1); fSpPCALam3 = eigvals(0);
        Eigen::Vector3d v2 = eigvecs.col(1), v3 = eigvecs.col(0);
        fSpPCAv2x = v2(0); fSpPCAv2y = v2(1); fSpPCAv2z = v2(2);
        fSpPCAv3x = v3(0); fSpPCAv3y = v3(1); fSpPCAv3z = v3(2);

        // Sign convention: principal axis points in +Z direction (beam direction)
        if (dir(2) < 0) dir = -dir;

        fSpPCAx = dir(0); fSpPCAy = dir(1); fSpPCAz = dir(2);

        if (fVerbosity > 0)
          mf::LogInfo("NuSliceAnalyzer") << "  PCA axis: (" << fSpPCAx << ", " << fSpPCAy << ", " << fSpPCAz << ")";
      }
    }
  }

  fTree->Fill();
}

DEFINE_ART_MODULE(NuSliceAnalyzer)

} // namespace opdet
