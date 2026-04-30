// ============================================================================
// NuSliceAnalyzer_module.cc
//
// art EDAnalyzer that runs on events passing NuSliceFilter.
// For each event it:
//   1. Reads the bestSliceID produced by NuSliceFilter (no cut re-evaluation)
//   2. Extracts the Pandora vertex and NuScore of the selected slice
//   3. Computes SpacePoint barycenter and PCA (charge-weighted by collection-plane hit integral)
//   4. Reads the ResNet (NuIntNNProducer_posdir_module) predictions from PixelMapVars data product
//   5. Fills a TTree with all the above
//
// FCL parameters:
//   PandoraLabel     : input tag for pandora products (default: "pandora")
//   FilterLabel      : input tag for NuSliceFilter (default: "nuselector")
//   ResNetLabel      : input tag for PixelMapVars data product (default: "ResNetInference")
//   Verbosity        : 0=minimal, 1=info, 2=debug
// ============================================================================

// --- art framework ---
// EDAnalyzer: base class for read-only analysis modules. analyze() is called
//   once per event but cannot write data products back (use EDProducer for that).
// ServiceHandle: access to art services like TFileService from within a module.
// TFileService: ROOT I/O service that manages output files and ROOT objects (TTrees, histograms).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// --- LArSoft data products ---
// Slice, PFParticle, Vertex, PFParticleMetadata: produced by Pandora (see NuSliceFilter comments)
// SpacePoint: a 3D position reconstructed by Pandora by combining hits from multiple wire planes.
//   Each SpacePoint has an (X,Y,Z) coordinate in detector space (cm).
// Hit: a reconstructed signal peak on a single wire at a specific drift time.
//   Each SpacePoint is associated to one Hit per wire plane. We use the collection-plane
//   (plane 2) hit integral as the charge weight for barycenter and PCA calculations.
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// PixelMapVars: custom data product written by NuIntNNProducer_posdir_module (ResNetInference module).
//   Contains the CNN position predictions: dEpromx_pred, dEpromy_pred, dEpromz_pred.
#include "sbndcode/NuIntNNOpReco/3-inference-larsoft-module/module/PixelMapVars.h"

// ROOT TTree for output ntuple.
#include "TTree.h"

// Eigen: linear algebra library used for the 3x3 PCA covariance matrix decomposition.
//   SelfAdjointEigenSolver computes eigenvalues and eigenvectors of a symmetric matrix.
#include <Eigen/Dense>

// C++
#include <vector>
#include <cmath>
#include <string>

namespace opdet {

class NuSliceAnalyzer : public art::EDAnalyzer {
public:
  explicit NuSliceAnalyzer(fhicl::ParameterSet const& p);
  void beginJob() override;   // called once before the first event: creates TTree and branches
  void analyze(art::Event const& e) override; // called once per event

private:
  // --- FCL parameters ---
  art::InputTag fPandoraLabel; // label of the Pandora process ("pandora")
  art::InputTag fFilterLabel;  // label of the NuSliceFilter module ("nuselector")
  art::InputTag fResNetLabel;  // label of the ResNet inference module ("ResNetInference")
  int           fVerbosity;    // 0=silent, 1=per-event info, 2=detailed debug

  // --- Output TTree (owned by TFileService, not by us) ---
  TTree* fTree;

  // --- TTree branch variables ---
  // Convention: -999 means "not filled" (event was rejected or data product missing).
  int    fRun, fSubrun, fEvent;  // event identifier
  int    fSliceID;               // ID of the best neutrino slice selected by NuSliceFilter
  float  fSliceNuScore;          // Pandora NuScore of the selected slice [0,1]

  // Pandora 3D vertex of the selected neutrino interaction (cm)
  double fVtxX, fVtxY, fVtxZ;

  // Charge-weighted barycenter of all SpacePoints in the selected slice (cm).
  // Weight = collection-plane hit integral (ADC x ticks), proportional to ionization charge.
  double fSpBaryX, fSpBaryY, fSpBaryZ;

  // First principal axis from PCA of the SpacePoint cloud (unit vector).
  // Sign convention: the Z component is always >= 0 (points in the beam direction).
  double fSpPCAx, fSpPCAy, fSpPCAz;

  // PCA eigenvalues of the charge-weighted covariance matrix (Lam1 >= Lam2 >= Lam3).
  // Lam1 corresponds to the main axis (spPCAx/y/z); large Lam1/Lam2 ratio means track-like.
  double fSpPCALam1, fSpPCALam2, fSpPCALam3;

  // Second and third PCA eigenvectors (unit vectors, orthogonal to each other and to the main axis).
  double fSpPCAv2x, fSpPCAv2y, fSpPCAv2z;
  double fSpPCAv3x, fSpPCAv3y, fSpPCAv3z;

  int    fNSpacePoints; // total number of SpacePoints in the selected slice

  // ResNet CNN predicted position of the neutrino interaction vertex (cm).
  // Read from the PixelMapVars data product produced by NuIntNNProducer_posdir_module.
  double fCnnPredX, fCnnPredY, fCnnPredZ;

  // Sets all branch variables to their default "not filled" values (-999 or 0).
  // Called at the start of each event so stale values from a previous event never leak.
  void ResetVars();
};

// ----------------------------------------------------------------------------
// Constructor: reads FCL parameters. TTree creation is deferred to beginJob()
// because TFileService is not available yet at construction time.
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
// ResetVars: initializes all branch variables to sentinel values.
// -999 is used as "missing data" sentinel so analysis code can distinguish
// events where a quantity was not filled from events where it is genuinely zero.
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
// beginJob: creates the TTree and registers all branches.
// Called once before the first event. TFileService writes the TTree to the
// output ROOT file under the directory named after this module.
// Branch format string "name/T" specifies the ROOT type: I=int, F=float, D=double.
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
  // PCA eigenvalues (Lam1>=Lam2>=Lam3, charge-weighted covariance; Lam1 eigenvector = spPCAx/y/z)
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
// analyze: main per-event function.
// Because this module lives in an end_path (see FCL), art only calls it for
// events that were NOT dropped by the trigger path filter (NuSliceFilter).
// So every event that reaches this function already has a valid neutrino slice.
// ----------------------------------------------------------------------------
void NuSliceAnalyzer::analyze(art::Event const& e) {
  ResetVars();

  fRun    = e.run();
  fSubrun = e.subRun();
  fEvent  = e.event();

  // =========================================================================
  // 1. Read ResNet predictions from PixelMapVars data product
  // =========================================================================
  // PixelMapVars is written by NuIntNNProducer_posdir_module (ResNetInference module) in the same job.
  // It contains the CNN predicted position (dEpromx/y/z_pred) and a flag passed_filters
  // that indicates whether the event passed the optical pre-selection inside the ResNet module.
  // We read this first so that even if the Pandora/TPC part fails, we still record CNN output.
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
  // 2. Read bestSliceID produced by NuSliceFilter
  // =========================================================================
  // NuSliceFilter wrote an int with instance name "bestSliceID" into the event.
  // We read it using an InputTag of the form (moduleLabel, instanceName).
  // If it is -1 or missing, this event has no valid neutrino slice; we still fill the
  // TTree so that the CNN prediction row is not lost (it will have sliceID=-1).
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
  // We need four association tables in addition to the two main collections:
  //   slice_pfp_assns : Slice -> PFParticles in that slice
  //   pfp_vtx_assns   : PFParticle -> its 3D vertex
  //   pfp_meta_assns  : PFParticle -> its metadata (NuScore, etc.)
  //   pfp_sp_assns    : PFParticle -> its SpacePoints
  // Plus a SpacePoint->Hit association to get charge weights.
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

  // SpacePoint->Hit association: used to get per-SpacePoint charge from the collection plane.
  art::Handle<std::vector<recob::SpacePoint>> spHandle;
  e.getByLabel(fPandoraLabel, spHandle);
  art::FindManyP<recob::Hit> sp_hit_assns(spHandle, e, fPandoraLabel);

  std::vector<art::Ptr<recob::Slice>> sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);

  // Find the slice in the collection whose ID matches the one selected by NuSliceFilter.
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
  // The primary PFP with PDG 12 (nu_e) or 14 (nu_mu) is the root of Pandora's
  // neutrino hierarchy. It is the only PFP that carries NuScore in its metadata
  // and has the interaction vertex associated to it.
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
    break; // There is exactly one primary neutrino PFP per slice
  }

  if (fVerbosity > 0)
    mf::LogInfo("NuSliceAnalyzer") << "  Best slice ID=" << fSliceID
                                   << " NuScore=" << fSliceNuScore
                                   << " vtx=(" << fVtxX << "," << fVtxY << "," << fVtxZ << ")";

  // =========================================================================
  // 5. Collect SpacePoints from all PFPs in the best slice
  // =========================================================================
  // SpacePoints are 3D points reconstructed by Pandora by triangulating hits across
  // wire planes. We gather all SpacePoints from all PFPs in the slice (not just the
  // primary), because daughters (tracks, showers) also carry charge information.
  // For each SpacePoint we look up the Hit on the collection plane (plane index 2)
  // and use its integral (area under the Gaussian fit) as the charge weight.
  // If no collection-plane hit is found, we default to weight=1 (unweighted).
  std::vector<double> spX, spY, spZ, spW;

  for (auto const& pfp : pfpVec) {
    const std::vector<art::Ptr<recob::SpacePoint>> spVec = pfp_sp_assns.at(pfp.key());
    for (auto const& sp : spVec) {
      // Default weight = 1 (used if collection-plane hit is not found)
      double integral = 1.0;
      if (sp_hit_assns.isValid()) {
        const std::vector<art::Ptr<recob::Hit>> hits = sp_hit_assns.at(sp.key());
        for (auto const& h : hits) {
          // Plane 2 is the collection plane in SBND (vertical wires, best S/N)
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
  // 6. Compute charge-weighted barycenter and PCA of the SpacePoint cloud
  // =========================================================================
  // Barycenter: sum(w_i * r_i) / sum(w_i), gives the charge centroid of the interaction.
  // PCA: decompose the 3x3 charge-weighted covariance matrix to find the main axis
  //   of the SpacePoint cloud. Useful as a direction estimate and to distinguish
  //   track-like (Lam1 >> Lam2,Lam3) from shower-like (Lam1 ~ Lam2 >> Lam3) events.
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
        // Build the 3x3 charge-weighted covariance matrix: C = sum(w_i * (r_i - mu)(r_i - mu)^T) / sum(w_i)
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < spX.size(); ++i) {
          Eigen::Vector3d pt(spX[i] - cx, spY[i] - cy, spZ[i] - cz);
          cov += std::abs(spW[i]) * (pt * pt.transpose());
        }
        cov /= totalW;

        // SelfAdjointEigenSolver is used because C is symmetric by construction.
        // Eigen sorts eigenvalues in ASCENDING order, so col(2) is the dominant axis.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
        Eigen::Vector3d eigvals = solver.eigenvalues();
        Eigen::Matrix3d eigvecs = solver.eigenvectors();
        int maxIdx;
        eigvals.maxCoeff(&maxIdx);
        Eigen::Vector3d dir = eigvecs.col(maxIdx);
        // Eigen sorts ascending: col(2)=Lam1 (dominant), col(1)=Lam2, col(0)=Lam3
        fSpPCALam1 = eigvals(2); fSpPCALam2 = eigvals(1); fSpPCALam3 = eigvals(0);
        Eigen::Vector3d v2 = eigvecs.col(1), v3 = eigvecs.col(0);
        fSpPCAv2x = v2(0); fSpPCAv2y = v2(1); fSpPCAv2z = v2(2);
        fSpPCAv3x = v3(0); fSpPCAv3y = v3(1); fSpPCAv3z = v3(2);

        // Eigenvectors have an arbitrary sign flip. We fix the sign so that the
        // principal axis always points towards +Z (beam direction) for consistency.
        if (dir(2) < 0) dir = -dir;

        fSpPCAx = dir(0); fSpPCAy = dir(1); fSpPCAz = dir(2);

        if (fVerbosity > 0)
          mf::LogInfo("NuSliceAnalyzer") << "  PCA axis: (" << fSpPCAx << ", " << fSpPCAy << ", " << fSpPCAz << ")";
      }
    }
  }

  // Always fill the TTree, even if some quantities are -999.
  // The analysis can filter on sliceID >= 0 or vtxX != -999 in post-processing.
  fTree->Fill();
}

DEFINE_ART_MODULE(NuSliceAnalyzer)

} // namespace opdet
