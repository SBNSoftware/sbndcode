////////////////////////////////////////////////////////////////////////
// Class:       FRANSProducer
// Plugin Type: Producer (art v3_05_01)
// File:        FRANSProducer_module.cc
////////////////////////////////////////////////////////////////////////

//#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "lardata/Utilities/AssociationUtil.h"

// services
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// back tracking
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MergeSimSources/MergeSimSources.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// truth includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// reco includes
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

// FRANS includes
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/FRANS/ChargeDensity.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/FRANS/ChargeDensityPset.h"
// TPCLines includes
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesParameters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesAlgo.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/EventHandle/SEventHandle.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleTriangles.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/LambdaAnaTree.h"

// ROOT inclides
#include "TTree.h"
#include "TLatex.h"
#include "TMVA/RBDT.hxx"

// std includes
#include <vector>
#include <limits>
#include <map>

// Michel vars
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelExtractor.hh"

// Lambda products
// #include "sbndcode/QELambdaAnalysis/Objects/LambdaProducts.hh"
#include "sbnobj/SBND/QELambdaReco/LambdaProducts.hh"

#define fXFidCut1 1.5
#define fXFidCut2 190
#define fYFidCut 190
#define fZFidCut1 10
#define fZFidCut2 490

#define fXAVCut 200
#define fYAVCut 200
#define fZAVCut1 0
#define fZAVCut2 500

#define fDefaulNeutrinoID 99999

// inline constexpr int defInt = std::numeric_limits<int>::min();


// struct PreselectionInfo{
//   bool IsClearCosmic;
//   double NuScore;
//   double BcFMScore;
//   bool RecoIsFiducial;
//   int nPrimPfps;
//   int nPrimTrks;
//   int nPrimShwrs;
//   int nTrksExits;
//   double maxGap;
//   bool passPreselection;

//   void clear(){
//     IsClearCosmic = true;
//     NuScore = -1;
//     BcFMScore = -1;
//     RecoIsFiducial = false;
//     nPrimPfps = 0;
//     nPrimTrks = 0;
//     nPrimShwrs = 0;
//     nTrksExits = 0;
//     maxGap = 0;
//     passPreselection = false;
//   }
// };

struct SliceInfo{
  int SliceID;
  std::vector<art::Ptr<recob::Hit>> SlcHits;
  std::vector<art::Ptr<recob::Hit>> SlcRecobHits;
  std::vector<int> SlcRecobHitsClusterId;
  STriangle Triangle;
  LambdaAnaTree LambdaAnaInfo;
  // PreselectionInfo LambdaPreInfo;
  MichelFeatures michelFeatures;
  int Track1PfpID;
  int Track2PfpID;
  int MainTrackPfpID;

  void clear(){
    SliceID = 0;
    SlcHits.clear();
    SlcRecobHits.clear();
    SlcRecobHitsClusterId.clear();
    Triangle = STriangle();
    LambdaAnaInfo.ResetVars();
    // fFRANSslc.LambdaPreInfo.clear();
    michelFeatures.clear();
    Track1PfpID = -1;
    Track2PfpID = -1;
    MainTrackPfpID = -1;
  }
  
};


namespace qelambda {
  class FRANSProducer;
}

class qelambda::FRANSProducer : public art::EDProducer {
public:


  explicit FRANSProducer(fhicl::ParameterSet const& p);


  // Required functions.
  void produce(art::Event &e) override;

  void beginSubRun(const art::SubRun &sr);
  void endSubRun(const art::SubRun& sr);

private:

  // === HANDLES ===
  art::ServiceHandle<cheat::BackTrackerService> fbtServ;
  art::ServiceHandle<cheat::ParticleInventoryService> fpiServ = art::ServiceHandle<cheat::ParticleInventoryService>();
  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // TFileService
  art::ServiceHandle<art::TFileService> tfs;

  // void respondToCloseInputFile(const art::FileBlock& fb);

  // --- Configuration parameters
  std::string fMCTruthLabel;
  std::string fMCParticleLabel;
  std::string fHitLabel;
  std::string fVertexLabel;
  std::string fRecoLabel;
  std::string fCRUMBSLabel;
  std::string fTPCPMTBarycenterFMModuleLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleMetadataLabel;

  int fView;
  float fBcFMScoreCut;
  float fTrkLenCut;
  float fShwrECut;
  float fVtxDistCut;
  float fShwrNHits;
  int fnPrimPfpsCut1;
  int fnPrimPfpsCut2;
  int fnPrimTrksCut1;
  int fnPrimTrksCut2;
  int fnPrimShwrsCut1;
  int fnPrimShwrsCut2;
  int fnPrimTrksExitCut;
  double fminGapCut;

  // --- Configuration parameters FRANS
  fhicl::ParameterSet fFhiclPsetFRANS;
  FRANSPsetType fPsetFRANS;
  
  // --- Configuration parameters TPCLines
  fhicl::ParameterSet fFhiclPsetTrackFinder;
  fhicl::ParameterSet fFhiclPsetVertexFinder;
  fhicl::ParameterSet fFhiclPsetHough;
  fhicl::ParameterSet fFhiclPsetTPCLines;
  fhicl::ParameterSet fFhiclPsetCalo;
  TrackFinderAlgorithmPsetType fPsetTrackFinder;
  HoughAlgorithmPsetType fPsetHough;
  VertexFinderAlgorithmPsetType fPsetVertexFinder;
  TPCLinesAlgoPsetType fPsetTPCLines;
  CaloAlgorithmPsetType fPsetCalo;

  // Michel BDT configuration
  std::string fmichelModelConfig;
  
  // --- TPC Lines algorithm
  TPCLinesAlgo* fTPCLinesAlgo;

  // --- Store the vertex
  VertexView fTrueVertex;
  VertexView fRecoVertex;
  int fTPC;

  // --- Constants
  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;

  // --- Maps
  std::unordered_map<int, art::Ptr<recob::PFParticle>> fPfpIDToPfpMap;
  std::map<art::Ptr<simb::MCTruth>, int> fMCTHitMap;
  std::map<int, int> fTrackIDHitMap;
  std::map<int, art::Ptr<recob::SpacePoint>> fhitToSpacePointMap;

  // --- Geometry
  geo::WireReadoutGeom const& m_readout = art::ServiceHandle<geo::WireReadout>()->Get();

  // --- Main function to get hits for FRANS and TPCLines
  std::vector<SHit> GetFRANSHitsView(int _view, int _tpc, std::vector<art::Ptr<recob::Hit>> hitV,
                                     std::vector<int> hitClusterId,
                                     std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap);

  
  bool CheckClearCosmic(std::vector<art::Ptr<recob::PFParticle>> pfpVec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn);
  bool FindNuPfp(art::Ptr<recob::PFParticle> &nupfp, std::vector<art::Ptr<recob::PFParticle>> pfpVect);

  // --- Fiducial volume functions
  bool PointInFV(double x, double y, double z);
  bool PointInFVNoCathode(double x, double y, double z);
  bool PointInAV(double x, double y, double z);
  
  // --- Reset variables
  void ResetVars();

  void setupHitMaps(
    const art::Event &e, 
    const art::Handle<std::vector<recob::SpacePoint>>& spHandle);
  
  // --- Distance and conversion functions
  bool GetWirePosition(int ch, double xyz_start[3], double xyz_end[3]);
  double DriftTickToVertex(double tt, int tpc);
  int VertexToDriftTick(double vt, double vx);
  VertexView getPointInWires(const double& x, const double& y, const double& z, bool reco);

  std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<int>> GetAssociatedHits(
    const art::Event &e,
    const std::vector<art::Ptr<recob::PFParticle>>& pfpVect,
    const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
    const art::Handle<std::vector<recob::Cluster>>& clusterHandle);

  bool ApplyLambdaPreselection(
    const art::Event &e,
    const art::Ptr<recob::Slice>& slice,
    const std::vector<art::Ptr<recob::PFParticle>>& slcPfpVec,
    const art::Handle<std::vector<recob::Slice>> &sliceHandle,
    const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
    const art::Handle<std::vector<recob::Shower>> &showerHandle,
    art::Ptr<recob::PFParticle>& nupfp);

  bool getFRANSVertices();

  art::Ptr<recob::PFParticle> getMuonPfp(
    const art::Event &e,
    const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  void getVPfps(
    const art::Event &e,
    const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  // --- Lambda info maker
  void updateLambdaInfo();

  // TTree variables
  TTree *fTreeHeader;
  std::string fInputFileName;
  TTree* fTreeAna;
  LambdaAnaTree fAnaTreeHandle;
  
  // --- Output ana tree

  // === OUTPUT VARIABLES ==
  qelambda::QELambdaReco fLambdaInfo = qelambda::QELambdaReco();

  // File variables
  std::string ffileName;

  // Event variables
  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;
  unsigned int fSliceID;

  // FRANS variables  
  bool fvalidFransVertex;
  bool fFransReco;
  double fbestFRANSScore;
  int fbestSlice;
  bool fIsQELambdaSignal;
  std::vector<art::Ptr<recob::Hit>> fLambdaHits;


  // Michel Tagging
  double fmuonTrueEndX;
  double fmuonTrueEndY;
  double fmuonTrueEndZ;
  fhicl::ParameterSet fMichelPset;
  bool ftrueMuonInFV;
  bool fmuDecay;
  float fmichelScore;
  MichelFeatures fmichelFeatures;
  bool ffeaturesComputed = false;

  MichelClassifier fmichelBDT;


  // Reco vertex
  bool fRecoIsFiducial;

  // FRANS starting positions
  double fHyperonVertexX, fHyperonVertexY, fHyperonVertexZ;
  double fbestHyperonVertexX, fbestHyperonVertexY, fbestHyperonVertexZ;
  double fPrimaryLeptonVertexX, fPrimaryLeptonVertexY, fPrimaryLeptonVertexZ;
  double fbestPrimaryLeptonVertexX, fbestPrimaryLeptonVertexY, fbestPrimaryLeptonVertexZ;

  // Lambda preselection variables
  int fnPrimPfps;
  int fnPrimTrks;
  int fnPrimShwrs;
  int fnTrksExits;
  double fmaxGap;

  // Best FRANS slice info
  SliceInfo fFRANSslc;

};

DEFINE_ART_MODULE(qelambda::FRANSProducer)
