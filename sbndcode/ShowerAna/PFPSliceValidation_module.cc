////////////////////////////////////////////////////////////////////////////////
// Class:       PFPSliceValidation                                            //
// Plugin Type: analyzer (art v3_02_06)                                       //
// File:        PFPSliceValidation_module.cc                                  //
//                                                                            //
// Generated at Wed Oct  2 03:27:09 2019 by Edward Tyley using cetskelgen     //
// from cetlib version v3_07_02.                                              //
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Root Includes
#include "TTree.h"
#include <iostream>
#include <vector>

namespace ana {
class PFPSliceValidation;
}

class ana::PFPSliceValidation : public art::EDAnalyzer {
public:
  explicit PFPSliceValidation(fhicl::ParameterSet const &pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPSliceValidation(PFPSliceValidation const &) = delete;
  PFPSliceValidation(PFPSliceValidation &&)      = delete;
  PFPSliceValidation &operator=(PFPSliceValidation const &) = delete;
  PFPSliceValidation &operator=(PFPSliceValidation &&) = delete;

  // Required functions.
  void analyze(art::Event const &evt) override;
  void beginJob() override;

  // Function to get a map of MCTruth to number of hits from that truth
  std::map<art::Ptr<simb::MCTruth>, int>
  GetTruthHitMap(const detinfo::DetectorClocksData &clockData,
                 const sim::ParticleList &trueParticlesMap,
                 const std::map<int, art::Ptr<simb::MCTruth>> &particleTruthMap,
                 const std::vector<art::Ptr<recob::Hit>> &allHits);

  // Function to match a slice, really any selection of hits, back to MCTruth
  // Also calculates purity and completeness of match (by reference)
  art::Ptr<simb::MCTruth>
  GetSliceTruthMatchHits(const detinfo::DetectorClocksData &clockData,
                         const std::vector<art::Ptr<recob::Hit>> &sliceHits,
                         const std::map<int, art::Ptr<simb::MCTruth>> &particleTruthMap,
                         const std::map<art::Ptr<simb::MCTruth>, int> &truthHitMap,
                         float &completeness, float &purity);

  // Functions to reset tree variables
  void ClearTrueTree();
  void ClearEventTree();

  // Function to create branches on tree for a list of labels
  // Branches will have the form branchName_PFParticleLabel
  template <class T>
  void initTree(TTree *Tree, std::string branchName, std::map<std::string, T> &Metric,
                std::vector<std::string> fPFParticleLabels);

  struct SliceMatch {

    // Default constructor initialise all of the things to -999
    SliceMatch()
        : mRecoId(-999)
        , mRecoPdg(-999)
        , mComp(-999)
        , mPurity(-999)
        , mNuScore(-999)
        , mVtxX(-999)
        , mVtxY(-999)
        , mVtxZ(-999){};

    SliceMatch(int recoId, int recoPdg, float comp, float purity, float nuScore, double vtx[3])
        : mRecoId(recoId)
        , mRecoPdg(recoPdg)
        , mComp(comp)
        , mPurity(purity)
        , mNuScore(nuScore)
        , mVtxX(vtx[0])
        , mVtxY(vtx[1])
        , mVtxZ(vtx[2]){};

    int mRecoId, mRecoPdg;
    float mComp, mPurity, mNuScore;
    float mVtxX, mVtxY, mVtxZ;

    // bool operator<(const SliceMatch& match) const {
    //   return mComp < match,mComp;
    // }
  };

private:
  int fVerbose;

  std::string fHitLabel, fGenieGenModuleLabel;
  std::vector<std::string> fPFParticleLabels;

  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  TTree *eventTree;
  TTree *trueTree;

  // Event wide metrics
  int eventTrueNeutrinos;
  std::map<std::string, int> eventPFPSlices, eventPFPNeutrinos;
  std::map<std::string, std::vector<float>> eventCosmicScores, eventNeutrinoScores;

  // Truth-by-truth metrics
  std::map<std::string, bool> nuMatchNeutrino;
  std::map<std::string, int> nuSlices, nuNeutrinos, bestNuPdg;
  std::map<std::string, float> bestNuPurity, bestNuComp, bestNuScore;

  // Neutrino Interaction variables
  int intType, CCNC, neutrinoPDG, numProtons, numNeutrons, numPi, numPi0, numTrueHits;
  float W, X, Y, QSqr, Pt, Theta, neutrinoE, leptonP;
  float trueVertexX, trueVertexY, trueVertexZ;

  // Vertex reco information
  std::map<std::string, float> pfpVertexX, pfpVertexY, pfpVertexZ;
  std::map<std::string, float> pfpVertexDistX, pfpVertexDistY, pfpVertexDistZ, pfpVertexDistMag;
};

ana::PFPSliceValidation::PFPSliceValidation(fhicl::ParameterSet const &pset)
    : EDAnalyzer{pset}
    , fVerbose(pset.get<int>("Verbose", 0))
    , fHitLabel(pset.get<std::string>("HitLabel"))
    , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
    , fPFParticleLabels(pset.get<std::vector<std::string>>("PFParticleLabels")) {}

void ana::PFPSliceValidation::beginJob() {
  trueTree  = tfs->make<TTree>("trueTree", "Tree with true neutrino metrics");
  eventTree = tfs->make<TTree>("eventTree", "Tree with event wide metrics");

  eventTree->Branch("trueNeutrinos", &eventTrueNeutrinos);

  initTree(eventTree, "pfpNeutrinos", eventPFPNeutrinos, fPFParticleLabels);
  initTree(eventTree, "pfpSlices", eventPFPSlices, fPFParticleLabels);
  // Vector of scores of all slices that match to a cosmic
  initTree(eventTree, "cosmicScores", eventCosmicScores, fPFParticleLabels);
  // Vector of socres of all slices that match to a neutrino
  initTree(eventTree, "nuScores", eventNeutrinoScores, fPFParticleLabels);

  // Truth variables from GENIE
  trueTree->Branch("intType", &intType);
  trueTree->Branch("CCNC", &CCNC);
  trueTree->Branch("neutrinoPDG", &neutrinoPDG);
  trueTree->Branch("numProtons", &numProtons); // Note, 21Mev KE cut
  trueTree->Branch("numPi", &numPi);
  trueTree->Branch("numPi0", &numPi0);
  trueTree->Branch("numTrueHits", &numTrueHits);
  trueTree->Branch("W", &W);
  trueTree->Branch("X", &X);
  trueTree->Branch("Y", &Y);
  trueTree->Branch("QSqr", &QSqr);
  trueTree->Branch("Pt", &Pt);
  trueTree->Branch("Theta", &Theta);
  trueTree->Branch("neutrinoE", &neutrinoE);
  trueTree->Branch("leptonP", &leptonP);

  // Total number of all slices, and only neutrino slices matched to true interaction
  initTree(trueTree, "numSlices", nuSlices, fPFParticleLabels);
  initTree(trueTree, "numNeutrinos", nuNeutrinos, fPFParticleLabels);
  // Metrics of only best matched slice (Highest completeness)
  initTree(trueTree, "bestMatchNeutrino", nuMatchNeutrino, fPFParticleLabels);
  initTree(trueTree, "purity", bestNuPurity, fPFParticleLabels);
  initTree(trueTree, "comp", bestNuComp, fPFParticleLabels);
  initTree(trueTree, "score", bestNuScore, fPFParticleLabels);
  initTree(trueTree, "recoPdg", bestNuPdg, fPFParticleLabels);
  // True vertex, needed for FV cuts
  trueTree->Branch("trueVertexX", &trueVertexX);
  trueTree->Branch("trueVertexY", &trueVertexY);
  trueTree->Branch("trueVertexZ", &trueVertexZ);
  // reco vertex of best matched slice, only available for neutrino slices
  initTree(trueTree, "pfpVertexX", pfpVertexX, fPFParticleLabels);
  initTree(trueTree, "pfpVertexY", pfpVertexY, fPFParticleLabels);
  initTree(trueTree, "pfpVertexZ", pfpVertexZ, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistX", pfpVertexDistX, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistY", pfpVertexDistY, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistZ", pfpVertexDistZ, fPFParticleLabels);
  initTree(trueTree, "pfpVertexDistMag", pfpVertexDistMag, fPFParticleLabels);
}

void ana::PFPSliceValidation::analyze(art::Event const &evt) {
  ClearEventTree();

  // Get the truths in the event:
  const std::vector<art::Ptr<simb::MCTruth>> truthVec = particleInventory->MCTruthVector_Ps();

  std::cout << std::setprecision(1) << std::fixed;
  if (fVerbose) {
    for (auto const &truth : truthVec) {
      std::cout << "Truth: " << truth << std::endl;
      if (truth->NeutrinoSet()) {
        const simb::MCNeutrino neutrino = truth->GetNeutrino();
        std::cout << "Neutrino: " << neutrino << std::endl;

        const simb::MCParticle nu = neutrino.Nu();
        std::cout << "X: " << nu.Vx() << " Y: " << nu.Vy() << " Z " << nu.Vz() << std::endl;
      } // truth->NeutrinoSet
    } // fVerbose
  } // truth: truthVec
  std::cout << std::setprecision(2) << std::fixed;

  // Get a map of each true particle to the MC Truth
  std::map<int, art::Ptr<simb::MCTruth>> particleTruthMap;
  const sim::ParticleList &trueParticlesMap = particleInventory->ParticleList();
  for (auto const &[trackId, particle] : trueParticlesMap) {
    particleTruthMap[trackId] = particleInventory->ParticleToMCTruth_P(particle);
  } // [trackId, particle]: trueParticlesMap
  eventTrueNeutrinos = truthVec.size();

  // Get reco stuff initialised
  //std::vector<art::Handle<std::vector<recob::Hit>>> hitHandles;
  //evt.getManyByType(hitHandles);
  auto hitHandles = evt.getMany<std::vector<recob::Hit>>();

  // Set the handles
  art::Handle<std::vector<recob::Hit>> hitHandle;
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  art::Handle<std::vector<recob::Slice>> sliceHandle;

  // Get all the hits
  std::vector<art::Ptr<recob::Hit>> allHits;
  if (evt.getByLabel(fHitLabel, hitHandle))
    art::fill_ptr_vector(allHits, hitHandle);

  // Get map of true primary particle to number of reco hits / energy in reco hits
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  std::map<art::Ptr<simb::MCTruth>, int> truthHitMap =
      GetTruthHitMap(clockData, trueParticlesMap, particleTruthMap, allHits);

  // Create maps to store the best matched slice to each truth
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, unsigned int>> pfpTruthNuCounterMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, unsigned int>> pfpTruthSliceCounterMap;
  std::map<std::string, std::map<art::Ptr<simb::MCTruth>, SliceMatch>> pfpTruthSliceMatchMap;

  // Initialise the counters in the maps
  for (auto const &fPFParticleLabel : fPFParticleLabels) {
    for (auto const &truth : truthVec) {
      pfpTruthNuCounterMap[fPFParticleLabel][truth]    = 0;
      pfpTruthSliceCounterMap[fPFParticleLabel][truth] = 0;
    } // truth: truthVec
  } // fPFParticleLabel: fPFParticleLabels

  for (auto const &fPFParticleLabel : fPFParticleLabels) {

    if (fVerbose) {
      std::cout << "On PFParticleLabel: " << fPFParticleLabel << std::endl;
    } // fVerbose

    // Get all the PFPs
    std::vector<art::Ptr<recob::Slice>> pfpSliceVec;
    if (evt.getByLabel(fPFParticleLabel, sliceHandle))
      art::fill_ptr_vector(pfpSliceVec, sliceHandle);

    std::vector<art::Ptr<recob::PFParticle>> pfps;
    if (evt.getByLabel(fPFParticleLabel, pfpHandle))
      art::fill_ptr_vector(pfps, pfpHandle);

    art::FindOneP<recob::Vertex> fopfv(pfpHandle, evt, fPFParticleLabel);
    if (!fopfv.isValid() || fopfv.size() == 0) {
      std::cout << "FindMany Vertex handle not valid" << std::endl;
        return;
    }

    art::FindManyP<recob::Hit> fmSliceHits(pfpSliceVec, evt, fPFParticleLabel);
    if (!fmSliceHits.isValid() || fmSliceHits.size() == 0) {
      std::cout << "FindMany Slice Hits not valid" << std::endl;
      return;
    }
    art::FindManyP<recob::PFParticle> fmSlicePFPs(pfpSliceVec, evt, fPFParticleLabel);
    if (!fmSlicePFPs.isValid()) {
      std::cout << "FindMany Slice PFPs not valid" << std::endl;
      return;
    }
    // Create a map between PFParticles and their IDs
    art::FindManyP<larpandoraobj::PFParticleMetadata> fmpfpmd(pfps, evt, fPFParticleLabel);
    if (!fmpfpmd.isValid()) {
      std::cout << "PFP MetaData handle not valid" << std::endl;
      return;
    }

    // Get maps of the pfp's Id to the pfp objects
    // And for pfp neutrinos, get the slice id mva score
    std::map<long unsigned int, art::Ptr<recob::PFParticle>> pfpMap;
    std::map<long unsigned int, float> pfpNuScoreMap;
    std::vector<art::Ptr<recob::PFParticle>> pfpNeutrinoVec;

    for (auto const &pfp : pfps) {

      long unsigned int pfpID = pfp->Self();
      pfpMap[pfpID]           = pfp;

      // Select PFP neutrinos
      if (pfp->PdgCode() == 12 || pfp->PdgCode() == 14) {
        pfpNeutrinoVec.push_back(pfp);
        if (fmpfpmd.size() == 0) {
          std::cout << "PFP neutrino has no metadata" << std::endl;
          return;
        } // fmpfpmd.size()

        // Get the pfparticle metadata to get the MVA score for each slice
        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec =
            fmpfpmd.at(pfpID);
        for (auto const pfpMeta : pfpMetaVec) {
          larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap =
              pfpMeta->GetPropertiesMap();
          pfpNuScoreMap[pfpID] = propertiesMap.at("NuScore");
        } // pfpMeta: pfpMetaVec
      } // pfp->PdgCode()==12 || pfp->PdgCode()==14
    } // pfp: pfps

    eventPFPSlices[fPFParticleLabel]    = pfpSliceVec.size();
    eventPFPNeutrinos[fPFParticleLabel] = pfpNeutrinoVec.size();

    for (const auto &pfpSlice : pfpSliceVec) {

      std::vector<art::Ptr<recob::Hit>> sliceHits        = fmSliceHits.at(pfpSlice.key());
      std::vector<art::Ptr<recob::PFParticle>> slicePFPs = fmSlicePFPs.at(pfpSlice.key());

      bool isNeutrinoSlice(false);
      float nuScore(-999);
      long unsigned int pfpNu(-999);
      // Check if it is a neutrino slice, if so get some info from the pfp neutrino
      for (auto const &pfp : slicePFPs) {
        if (pfp->PdgCode() == 12 || pfp->PdgCode() == 14) {
          pfpNu           = pfp->Self();
          nuScore         = pfpNuScoreMap[pfpNu];
          isNeutrinoSlice = true;
          break; // Should only be 1 neutrino per slice
        } // pfp->PdgCode()==12 || pfp->PdgCode()==14
      } // pfp:slicePFPs

      // Find the MCTruth that contains most of the hits from the slice
      float purity(-999), completeness(-999);
      art::Ptr<simb::MCTruth> trueMatch = GetSliceTruthMatchHits(
          clockData, sliceHits, particleTruthMap, truthHitMap, completeness, purity);

      // Check if it matched to anything
      if (trueMatch.isNull())
        continue;

      if (fVerbose && isNeutrinoSlice)
        std::cout << "True Match: " << trueMatch << " with completeness: " << completeness
                  << " and purity: " << purity << " and score: " << nuScore << std::endl;

      // Increment the counters for the true match
      ++pfpTruthSliceCounterMap[fPFParticleLabel][trueMatch];
      if (isNeutrinoSlice) {
        ++pfpTruthNuCounterMap[fPFParticleLabel][trueMatch];
        if (trueMatch->NeutrinoSet()) {
          eventNeutrinoScores[fPFParticleLabel].push_back(nuScore);
        } else { // cosmicTruth
          eventCosmicScores[fPFParticleLabel].push_back(nuScore);
        } // trueMatch->NeutrinoSet()
      } // isNeutrinoSlice

      // Choose the best match slice, defined as the slice with the best completeness
      if (completeness > pfpTruthSliceMatchMap[fPFParticleLabel][trueMatch].mComp) {
        if (isNeutrinoSlice) {

          art::Ptr<recob::PFParticle> pfpNeutrino = pfpMap.at(pfpNu);
          art::Ptr<recob::Vertex> pfpVertex       = fopfv.at(pfpNeutrino.key());
          double pfpVtx[3]                        = {-999, -999, -999};
          pfpVertex->XYZ(pfpVtx);

          SliceMatch match(pfpNu, pfpNeutrino->PdgCode(), completeness, purity, nuScore, pfpVtx);

          pfpTruthSliceMatchMap[fPFParticleLabel][trueMatch] = match;

        } else { // isNeutrinoSlice
          double pfpVtx[3] = {-999, -999, -999};
          SliceMatch match(-999, -999, completeness, purity, -999, pfpVtx);
          pfpTruthSliceMatchMap[fPFParticleLabel][trueMatch] = match;
        } // else isNeutrinoSlice
      } // bestMatch
    } // pfpSlice:pfpSliceVec
  } // auto const fPFParticleLabel: fPFParticleLabels

  eventTree->Fill();

  // Get the true neutrinos
  for (auto const &truth : truthVec) {

    ClearTrueTree();

    // Only fill the tree for truth neutrinos
    if (!truth->NeutrinoSet())
      continue;

    // Get the truth interaction variables
    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    const simb::MCParticle nu       = neutrino.Nu();
    const simb::MCParticle lepton   = neutrino.Lepton();

    intType     = neutrino.Mode();
    CCNC        = neutrino.CCNC();
    neutrinoPDG = nu.PdgCode();
    W           = neutrino.W();
    X           = neutrino.X();
    Y           = neutrino.Y();
    QSqr        = neutrino.QSqr();
    Pt          = neutrino.Pt();
    Theta       = neutrino.Theta();
    neutrinoE   = nu.E();
    leptonP     = lepton.P();

    trueVertexX = nu.Vx();
    trueVertexY = nu.Vy();
    trueVertexZ = nu.Vz();

    // Number of true hits from the slice
    numTrueHits = truthHitMap.at(truth);

    // Calculate the number of direct daughters
    for (auto const &[particleId, truthIter] : particleTruthMap) {
      if (truthIter != truth)
        continue;
      const simb::MCParticle *particle = trueParticlesMap.at(particleId);
      // We only want the primary daughters
      if (particle->Process() != "primary")
        continue;
      // Apply 21MeV KE cut on protons
      if (particle->PdgCode() == 2212 && (particle->E() - particle->Mass()) > 0.021) {
        ++numProtons;
      } else if (std::abs(particle->PdgCode()) == 211) {
        ++numPi;
      } else if (std::abs(particle->PdgCode()) == 111) {
        ++numPi0;
      }
    }

    if (fVerbose) {
      std::cout << "\nTruth: " << truth << std::endl;
    } // fVerbose

    for (auto const &fPFParticleLabel : fPFParticleLabels) {

      // Check we actually match a slice to the truth
      if (pfpTruthSliceCounterMap[fPFParticleLabel].at(truth)) {

        SliceMatch match = pfpTruthSliceMatchMap[fPFParticleLabel][truth];

        nuSlices[fPFParticleLabel]        = pfpTruthSliceCounterMap[fPFParticleLabel][truth];
        nuNeutrinos[fPFParticleLabel]     = pfpTruthNuCounterMap[fPFParticleLabel][truth];
        bestNuComp[fPFParticleLabel]      = match.mComp;
        bestNuPurity[fPFParticleLabel]    = match.mPurity;
        bestNuScore[fPFParticleLabel]     = match.mNuScore;
        bestNuPdg[fPFParticleLabel]       = match.mRecoPdg;
        nuMatchNeutrino[fPFParticleLabel] = (match.mRecoId != -999);

        // If we matched a neutrino slice, get the vertex info
        if (nuMatchNeutrino[fPFParticleLabel]) {

          pfpVertexX[fPFParticleLabel] = match.mVtxX;
          pfpVertexY[fPFParticleLabel] = match.mVtxY;
          pfpVertexZ[fPFParticleLabel] = match.mVtxZ;

          pfpVertexDistX[fPFParticleLabel] = pfpVertexX[fPFParticleLabel] - nu.Vx();
          pfpVertexDistY[fPFParticleLabel] = pfpVertexY[fPFParticleLabel] - nu.Vy();
          pfpVertexDistZ[fPFParticleLabel] = pfpVertexZ[fPFParticleLabel] - nu.Vz();

          pfpVertexDistMag[fPFParticleLabel] =
              std::hypot(pfpVertexDistX[fPFParticleLabel], pfpVertexDistY[fPFParticleLabel],
                         pfpVertexDistZ[fPFParticleLabel]);
        } // nuMatchNeutrino
      } // pfpTruthSliceCounter

      if (fVerbose) {
        std::cout << "PFParticleLabel: " << fPFParticleLabel << std::endl;

        std::cout << "Nu Slices: " << nuSlices[fPFParticleLabel]
                  << " and Nu Neutrinos: " << nuNeutrinos[fPFParticleLabel]
                  << " with best Nu Pdg: " << bestNuPdg[fPFParticleLabel]
                  << "\nCompleteness: " << bestNuComp[fPFParticleLabel]
                  << " and purity: " << bestNuPurity[fPFParticleLabel]
                  << " and score: " << bestNuScore[fPFParticleLabel] << std::endl;
      } // fVerbose
    } // fPFParticleLabel: fPFParticleLabels
    trueTree->Fill();
  } // truth: truthVec
  std::cout << "\n" << std::endl;
} // analyze

std::map<art::Ptr<simb::MCTruth>, int> ana::PFPSliceValidation::GetTruthHitMap(
    const detinfo::DetectorClocksData &clockData, const sim::ParticleList &trueParticlesMap,
    const std::map<int, art::Ptr<simb::MCTruth>> &particleTruthMap,
    const std::vector<art::Ptr<recob::Hit>> &allHits) {

  // Create a map of true particles to number of hits
  std::map<int, int> trueParticleHits;
  for (const auto &hit : allHits) {
    int trackID     = 0;
    float hitEnergy = 0;

    // For each hit, chose the particle that contributed the most energy
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (const auto &ide : trackIDEs) {
      if (ide.energy > hitEnergy) {
        hitEnergy = ide.energy;
        trackID   = std::abs(ide.trackID);
      } // ide.energy > hitEnergy
    } // ide: trackIDEs
    ++trueParticleHits[trackID];
  } // hit: allHits

  // Roll up the particles in their slices
  std::map<art::Ptr<simb::MCTruth>, int> truthHitMap;
  for (const auto &[trueParticle, truth] : particleTruthMap) {
    truthHitMap[truth] += trueParticleHits[trueParticle];
  } // [trueParticle, truth]: particleTruthMap

  return truthHitMap;
} // GetTruthHitMap

art::Ptr<simb::MCTruth> ana::PFPSliceValidation::GetSliceTruthMatchHits(
    const detinfo::DetectorClocksData &clockData,
    const std::vector<art::Ptr<recob::Hit>> &sliceHits,
    const std::map<int, art::Ptr<simb::MCTruth>> &particleTruthMap,
    const std::map<art::Ptr<simb::MCTruth>, int> &truthHitMap, float &completeness, float &purity) {

  // Create a map of true particles to number of hits
  std::map<int, int> trueParticleHits;
  for (const auto &hit : sliceHits) {
    int trackID     = 0;
    float hitEnergy = 0;

    // For each hit, chose the particle that contributed the most energy
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (const auto &ide : trackIDEs) {
      if (ide.energy > hitEnergy) {
        hitEnergy = ide.energy;
        trackID   = std::abs(ide.trackID);
      } // ide.energy > hitEnergy
    } // ide: trackIDEs
    ++trueParticleHits[trackID];
  } // hit: sliceHits

  // Roll up the particles in their slices
  std::map<art::Ptr<simb::MCTruth>, int> sliceTruthHitMap;
  for (const auto &[trueParticle, truth] : particleTruthMap) {
    sliceTruthHitMap[truth] += trueParticleHits[trueParticle];
  } // [trueParticle, truth]: particleTruthMap

  // Choose the truth that contributed the most hits
  int maxHits = 0;
  art::Ptr<simb::MCTruth> bestTruthMatch;
  for (const auto &[truth, truthHits] : sliceTruthHitMap) {
    if (truthHits > maxHits) {
      maxHits        = truthHits;
      bestTruthMatch = truth;
    } // truthHits > maxHuts
  } // [truth, truthHits]: sliceTruthHitMap

  // If we have truth matched the slice, calculate purtity and completeness
  // Note these are passed by referecne
  if (!bestTruthMatch.isNull()) {
    purity       = (float)maxHits / sliceHits.size();
    completeness = (float)maxHits / truthHitMap.at(bestTruthMatch);
  } // !bestTruthMatch.isNull()

  return bestTruthMatch;
} // GetSliceTruthMatchHits

void ana::PFPSliceValidation::ClearTrueTree() {

  intType     = -999;
  CCNC        = -999;
  neutrinoPDG = -999;
  numProtons  = 0;
  numNeutrons = 0;
  numPi       = 0;
  numPi0      = 0;
  numTrueHits = 0;
  W           = -999;
  X           = -999;
  Y           = -999;
  QSqr        = -999;
  Pt          = -999;
  Theta       = -999;
  neutrinoE   = -999;
  leptonP     = -999;

  trueVertexX = -999;
  trueVertexY = -999;
  trueVertexZ = -999;

  for (auto const &fPFParticleLabel : fPFParticleLabels) {

    nuMatchNeutrino[fPFParticleLabel] = false;
    nuSlices[fPFParticleLabel]        = -99999;
    nuNeutrinos[fPFParticleLabel]     = -99999;
    bestNuPurity[fPFParticleLabel]    = -99999;
    bestNuComp[fPFParticleLabel]      = -99999;
    bestNuScore[fPFParticleLabel]     = -99999;
    bestNuPdg[fPFParticleLabel]       = -99999;

    pfpVertexX[fPFParticleLabel] = -99999;
    pfpVertexY[fPFParticleLabel] = -99999;
    pfpVertexZ[fPFParticleLabel] = -99999;

    pfpVertexDistX[fPFParticleLabel]   = -99999;
    pfpVertexDistY[fPFParticleLabel]   = -99999;
    pfpVertexDistZ[fPFParticleLabel]   = -99999;
    pfpVertexDistMag[fPFParticleLabel] = -99999;
  } // fPFParticleLabel: fPFParticleLabels
} // ClearTrueTree

void ana::PFPSliceValidation::ClearEventTree() {
  eventTrueNeutrinos = -999;
  for (auto const &fPFParticleLabel : fPFParticleLabels) {
    eventPFPNeutrinos[fPFParticleLabel] = -999;
    eventPFPSlices[fPFParticleLabel]    = -999;
    eventCosmicScores[fPFParticleLabel].clear();
    eventNeutrinoScores[fPFParticleLabel].clear();
  } // fPFParticleLabel: fPFParticleLabels
} // ClearEventTree

template <class T>
void ana::PFPSliceValidation::initTree(TTree *Tree, std::string branchName,
                                       std::map<std::string, T> &Metric,
                                       std::vector<std::string> fPFParticleLabels) {

  for (auto const &fPFParticleLabel : fPFParticleLabels) {
    std::string branchString = branchName + "_" + fPFParticleLabel;
    const char *branchChar   = branchString.c_str();
    Tree->Branch(branchChar, &Metric[fPFParticleLabel], 32000, 0);
  } // fPFParticleLabel: fPFParticleLabels
} // initTree

DEFINE_ART_MODULE(ana::PFPSliceValidation)
