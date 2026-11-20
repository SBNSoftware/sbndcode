////////////////////////////////////////////////////////////////////////
// Class:       MichelTagging
// Plugin Type: analyzer (Unknown Unknown)
// File:        MichelTagging.cc
//
// Generated at Thu Nov 6 15:03:34 2025 by Sergio Yebes Baena using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/FileBlock.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"


// ROOT and C++ includes
#include <TTree.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <string.h>
#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <optional>

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include <larcorealg/Geometry/Exceptions.h>
#include "larcorealg/Geometry/GeometryCore.h"

// G4 includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/MCBase/MCTrack.h"

// Reco includes
// PDS
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
// TPC
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "sbndcode/BlipRecoSBND/Alg/BlipRecoAlg.h"

// CRT
#include "sbnobj/SBND/CRT/CRTTrack.hh"

// Cosmic rejection includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "lardataobj/AnalysisBase/T0.h"

// Geometry and mapping
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// Michel vars
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelDriver.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelClassifier.hh"


constexpr int defInt = std::numeric_limits<int>::min();

namespace sbnd {
  class MichelTagging;
}


class sbnd::MichelTagging : public art::EDAnalyzer {
public:
  explicit MichelTagging(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelTagging(MichelTagging const&) = delete;
  MichelTagging(MichelTagging&&) = delete;
  MichelTagging& operator=(MichelTagging const&) = delete;
  MichelTagging& operator=(MichelTagging&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // begin-end functions.
  void beginJob() override;
  void endJob() override;

  void beginSubRun(const art::SubRun &sr);
  void endSubRun(const art::SubRun &sr);

  //Files function
  void respondToOpenInputFile(art::FileBlock const&); 

  //TotalEvents function
  int getTotalGenEvents(const art::Event &e);


  // === OTHER FUNCTIONS FOR ANALYSIS === 

  void clearVars();
  void clearSubRunVars();
  void clearBlipVars();

  //void setupHitMaps(const art::Event& e, const art::ValidHandle<std::vector<recob::Hit>>& hitHandle);
  void setSlcHitsMaps(const art::Event &e, const std::vector<art::Ptr<recob::Hit>>& slcHitVec, const art::ValidHandle<std::vector<simb::MCParticle>>& MCPHandle);

  // void setDefaultGenVars();
  // void setDefaultG4Vars();
  // void setDefaultRecoVars();

  // void getNuInfo(const art::Ptr<simb::MCTruth> mct);
  void getEventBlips(art::ValidHandle<std::vector<blip::Blip>> blipHandle);
  void getGenInfo(const art::Ptr<simb::MCTruth> mct);
  void getSliceTrueInfo(const art::Event &e, const std::vector<art::Ptr<recob::Hit>>& slcHitVec, const art::ValidHandle<std::vector<simb::MCTruth>> MCTHandle);
  void getSlcVertex(const art::Event &e, const std::vector<art::Ptr<recob::PFParticle>>& slcPfpVec, const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle);
  std::vector<art::Ptr<recob::PFParticle>> getSlicePrimPfp(const std::vector<art::Ptr<recob::PFParticle>> slcPfpVec);

  // bool slcCheckClearCosmic(std::vector<art::Ptr < recob::PFParticle>> pfpVec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn);
  // double getSliceNuScore(std::vector<art::Ptr<recob::PFParticle>> pfpVec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn);

  void getG4Info(const std::vector<art::Ptr<simb::MCParticle>> mcpVec);
  
  std::array<std::vector<MuHit>, 3> getPfpHits(const art::Ptr<recob::PFParticle>& pfp, const std::vector<art::Ptr<recob::Hit>>& trackHitVec);
  void updateHits(const art::Ptr<recob::PFParticle>& pfp,const std::vector<art::Ptr<recob::Hit>>& trackHitVec);
  std::optional<std::array<std::vector<MuCalo>, 3>> getPfpCalo(
    const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
    const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle
    );
  void updateCalo(
    const std::vector<art::Ptr<anab::Calorimetry>> &caloVec,
    const std::vector<art::Ptr<recob::Hit>>& trackHitVec, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle
    );
  void concatenateMuonPFPs(
    const art::Event &e,
    const art::Ptr<recob::PFParticle>& initialPfp,
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
    const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle,
    const art::ValidHandle<std::vector<recob::Hit>> &hitHandle
    );
  void analyzePFPs(
    const art::Event &e,
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
    const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle,
    const art::ValidHandle<std::vector<recob::Hit>> &hitHandle
    );

  int getPfpG4ID(
    const art::Event &e,
    const std::vector<art::Ptr<recob::Hit>>& hitVec 
    );

  // void analyzeTrack(
  //   const art::Event &e, 
  //   const art::Ptr<recob::Track> &track,
  //   const art::ValidHandle<std::vector<recob::Track>> &trackHandle, 
  //   const std::vector<art::Ptr<recob::Hit>> &trackHitVec,
  //   const detinfo::DetectorClocksData& clockDataEvt
  //   );
  void getUnassociatedHits(
    const art::Event &e, 
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVec,
    const std::vector<art::Ptr<recob::Hit>>& slcHitVec, 
    const art::ValidHandle<std::vector<recob::Hit>>& hitHandle,
    const art::ValidHandle<std::vector<recob::PFParticle>> &pfpHandle, 
    const art::ValidHandle<std::vector<recob::Track>> &trackHandle
  );

  bool isMuonFV();

  int vertex2DriftTick(double vt, double vx);
  std::vector<int> getPointInWires(double x, double y, double z);



private:

  // === HANDLES ===
  art::ServiceHandle<cheat::BackTrackerService> fbtServ;
  art::ServiceHandle<cheat::ParticleInventoryService> fpiServ = art::ServiceHandle<cheat::ParticleInventoryService>();
  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Fhicl parameters set
  fhicl::ParameterSet fMichelPset;
  MichelParams fmichelParams;
  // === OUTPUT TREES ===
  TTree *fTree;

  // === OUTPUT VARIABLES ===

  // Michel variables
  MichelFeatures fmichelFeatures;
  bool ffeaturesComputed;

  // File variables
  std::string ffileName;

  // Event variables
  unsigned int feventID;
  unsigned int frunID;
  unsigned int fsubrunID;

  //MCTruth variables
  //Consult https://internal.dunescience.org/doxygen/classsimb_1_1MCNeutrino.html#a8e213917487f1f5cd52ca5280d33a777 for information on the numbers meaning
  int fnuPDG;
  double fnuE0;
  int fgenIndex;
  int fnuIntMode;
  int fnuIntType;
  int fnuCCNC;

  int fslcPrimID;
  int fslcPrimPDG;

  // True labels
  bool fmuFound;
  bool fmuDecay;
  int fprimMuID;
  double fprimMuLength;

  //Truth-matching variables
  std::map<int,int> fslcHitID2G4IDMap;
  std::map<int,int> fslcHitID2TruePdgMap;

  // --- Reco variables ---

  // Event variables
  unsigned int feventNSlices;

  // Slice variables
  unsigned int fslcID;

  double fslcVtxX, fslcVtxY, fslcVtxZ;
  int fslcVtxWireU, fslcVtxWireV, fslcVtxWireC;
  int fslcVtxChannelU, fslcVtxChannelV, fslcVtxChannelC;
  int fslcVtxTimeTick; 

  double fslcTrueVtxX, fslcTrueVtxY, fslcTrueVtxZ;
  int fslcTrueVtxWireU, fslcTrueVtxWireV, fslcTrueVtxWireC;
  int fslcTrueVtxChannelU, fslcTrueVtxChannelV, fslcTrueVtxChannelC;
  int fslcTrueVtxTimeTick;

  
  std::vector<MuTrack> fmuTracks;
  std::array<std::vector<MuHit>, 3> fmuHits;
  std::array<std::vector<MuCalo>, 3> fmuCalo;
  std::vector<MuBlip> fblips;
  int fnMuTracks;
  float fmichelScore;

  MichelClassifier fmichelBDT;

  // Bragg info
  bool fbraggFound;
  int fbraggIdx;
  int fbraggChannelID;
  double fbraggPeakTime;
  double fbraggTpcID;
  double fbraggResRange;
  double fbraggdEdx;
  double fbraggX;
  double fbraggY;
  double fbraggZ;

  // PCA info
  bool fpcaValid;
  double fpcaDirX;
  double fpcaDirY;

  // Angle vars
  std::vector<double> fangles;

  // Radial vars
  std::vector<double> fradii;
  std::vector<double> fcharge;
  std::vector<double> fdQdr;
  std::vector<double> fsmoothdQdr;
  int fdropIdx;
  double fmindQdr;

  // Hits info
  std::vector<int> fhitPlaneID;
  std::vector<int> fhitTpcID;
  std::vector<int> fhitChannelID;
  std::vector<double> fhitPeakTime; // Already rescaled by 1/4
  std::vector<double> fhitIntegral;
  std::vector<double> fhitRMS;
  std::vector<int> fhitTruePdg;
  std::vector<int> fhitPfpID;
  std::vector<int> fhitTrueID;

  // Calo info
  std::vector<int> fcaloPlaneID;
  std::vector<int> fcaloTpcID;
  std::vector<int> fcaloChannelID;
  std::vector<double> fcaloPeakTime; // Already rescaled by 1/4
  std::vector<int> fcaloTrueID;
  std::vector<int> fcaloTruePdg;
  std::vector<double> fcalodEdx;
  std::vector<double> fcaloResRange;
  std::vector<double> fcaloX, fcaloY, fcaloZ;

  // Detector parameters
  std::vector<double> PosNull = {0., 0., 0.};
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;
  double fOpticalClock;

  // fhicl parameters
  const std::string fMCTruthLabel;
  const std::string fMCParticleLabel;
  const std::string fPOTModuleLabel;
  const std::string fMCT2MCPLabel;
  const std::string fTrack2HitLabel;
  const std::string fShower2HitLabel;
  const std::string fSliceLabel;
  const std::string fSliceModuleLabel;
  const std::string fPFParticleLabel;
  const std::string fPFParticleMetadataLabel;
  const std::string fVertexLabel;
  const std::string fTrackLabel;
  const std::string fShowerLabel;
  const std::string fCalorimetryLabel;
  const std::string fHitLabel;
  const std::string fHit2SpacePointLabel;
  const std::string fBlipLabel;
  const std::string fTPCPMTBarycenterFMModuleLabel;

  // Reco parameters
  float fMaxVtxDist = 5; // cm
  float fDistanceThreshold = 10; // cm
  float fMinMuLength = 20; // cm
  float fFV = 10; // cm

  // Michel BDT configuration
  const std::string fmichelModelConfig;
};