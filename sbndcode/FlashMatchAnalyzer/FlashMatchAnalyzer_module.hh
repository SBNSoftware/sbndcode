////////////////////////////////////////////////////////////////////////
// Class:       FlashMatchAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        FlashMatchAnalyzer_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"


#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "TTree.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include "TH1D.h"

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>


#define fXFidCut1 1.5
#define fXFidCut2 190
#define fYFidCut 190
#define fZFidCut1 10
#define fZFidCut2 490


#define fDefaulNeutrinoID 99999

namespace test {
  class FlashMatchAnalyzer
;
}


class test::FlashMatchAnalyzer : public art::EDAnalyzer {
public:
  explicit FlashMatchAnalyzer
(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashMatchAnalyzer(FlashMatchAnalyzer const&) = delete;
  FlashMatchAnalyzer(FlashMatchAnalyzer&&) = delete;
  FlashMatchAnalyzer & operator=(FlashMatchAnalyzer const&) = delete;
  FlashMatchAnalyzer & operator=(FlashMatchAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  void resetVars();
  void resetTrueVars();
  void resetSimVars();
  void resetWireVars();
  void resetRecoVars();
  int VertexToDriftTick(double vt, double vx);
  bool PointInFV(double x, double y, double z);
  void FillReco2(art::Event const& e, std::vector<art::Ptr<recob::PFParticle>> pfpVect, std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap);
  void FillHits(int clusterId, std::vector<art::Ptr<recob::Hit>> hitVect, std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap);
  void AnalyseCRTTracks(const art::Event &, const std::vector<art::Ptr<sbnd::crt::CRTTrack>> &);
  std::string fMCTruthLabel;
  std::string fMCLabel;
  std::string fSimEnergyDepositLabel;
  std::string fSimEnergyDepositInstanceLabel;
  std::string fSimEnergyDepositLabelOut;
  std::string fSimEnergyDepositInstanceLabelOut;
  std::string fSimChannelLabel;
  std::string fRawDigitLabel;
  std::string fRecobWireLabel;
  std::string fHitLabel;
  std::string fReco2Label;
  std::string fTrackLabel;
  std::string fClusterLabel;
  std::string fSpacePointLabel;
  std::string fVertexLabel;
  std::string fT0label;
  std::string fCalorimetryLabel;
  std::string fParticleIDLabel;  
  std::vector<std::string> fOpFlashesModuleLabel;
  std::string fCRTTrackModuleLabel;
  bool fSaveReco1;
  bool fSaveReco2;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fSaveSimEDOut;
  bool fSaveWaveforms;
  bool fSaveWires;
  bool fSaveHitsReco1;
  bool fSaveHitsReco2;
  bool fSaveSpacePoints;
  bool fSaveVertex;
  bool fCreateTPCMap;
  bool fApplyFiducialCut;
  bool fApplyVertexSCE;
  bool fUseSlices;
  bool fUseSimChannels;
  bool fSaveOpHits;
  bool fSaveOpFlashes;
  std::vector<double> fSaveFlashWindow;
  std::vector<double> fSaveCRTWindow;
  bool fSaveOnlyNeutrino;
  bool fSaveCRTTracks;
  bool fComputePMTRatio;
  bool fSaveChargeBarycenter;

  TTree* fTree;
  TTree* fTreeReco1;
  TTree* fOpAnaTree;
  TTree* fCRTTree;
  int fEventID, fRunID, fSubRunID;

  //True variables
  std::vector<int> fTruePrimariesPDG;
  std::vector<double> fTruePrimariesE;
  std::vector<std::vector<double>> fTruePrimariesStartP;
  double fTrueVx;
  double fTrueVy;
  double fTrueVz;
  double fTrueVt;
  int fTrueVU;
  int fTrueVV;
  int fTrueVC;
  int fTrueVTimeTick;
  double fTrueVEnergy;
  int fIntMode;
  int fIntCCNC;
  int fIntNProtons;
  int fIntNNeutrons;
  int fIntNPi0;
  int fIntNPip;
  int fIntNPim;
  int fIntNMuonP;
  int fIntNMuonM;
  int fIntNElectronP;
  int fIntNElectronM;
  int fIntNLambda;
  bool fIntInFV;

  //True SimEnergyDeposits
  std::vector<double> fEnDepE;
  std::vector<double> fEnDepX;
  std::vector<double> fEnDepY;
  std::vector<double> fEnDepZ;
  std::vector<double> fEnDepU;
  std::vector<double> fEnDepV;
  std::vector<double> fEnDepC;
  std::vector<double> fEnDepT;
  std::vector<int>    fEnDepPDG;


  //True SimEnergyDeposits Out
  std::vector<double> fEnDepEOut;
  std::vector<double> fEnDepXOut;
  std::vector<double> fEnDepYOut;
  std::vector<double> fEnDepZOut;
  std::vector<double> fEnDepTOut;

  //Hit variables (reco2 level)
  std::vector<int> fHitsViewReco2;
  std::vector<double> fHitsPeakTimeReco2;
  std::vector<double> fHitsIntegralReco2;
  std::vector<double> fHitsSummedADCReco2;
  std::vector<double> fHitsChannelReco2;
  std::vector<double> fHitsAmplitudeReco2;
  std::vector<double> fHitsRMSReco2;
  std::vector<double> fHitsStartTReco2;
  std::vector<double> fHitsEndTReco2;
  std::vector<double> fHitsWidthReco2;
  std::vector<double> fHitsChi2Reco2;
  std::vector<double> fHitsNDFReco2;
  std::vector<int> fHitsClusterIDReco2;
  std::vector<double> fHitsXReco2;
  std::vector<double> fHitsYReco2;
  std::vector<double> fHitsZReco2;

  //Hit variables (reco1 level)
  std::vector<int> fHitsViewReco1;
  std::vector<double> fHitsPeakTimeReco1;
  std::vector<double> fHitsIntegralReco1;
  std::vector<double> fHitsSummedADCReco1;
  std::vector<double> fHitsChannelReco1;
  std::vector<double> fHitsAmplitudeReco1;
  std::vector<double> fHitsRMSReco1;
  std::vector<double> fHitsStartTReco1;
  std::vector<double> fHitsEndTReco1;
  std::vector<double> fHitsWidthReco1;
  std::vector<double> fHitsChi2Reco1;
  std::vector<double> fHitsNDFReco1;

  // Flash variables
  int _nopflash;
  std::vector<int> _flash_id;
  std::vector<double> _flash_time;
  std::vector<double> _flash_total_pe;
  std::vector<std::vector<double>> _flash_pe_v;
  std::vector<double> _flash_y;
  std::vector<double> _flash_yerr ;
  std::vector<double> _flash_z;
  std::vector<double> _flash_zerr;
  std::vector<double> _flash_x;
  std::vector<double> _flash_xerr;
  std::vector<int> _flash_tpc;
  std::vector<std::vector<double>> _flash_ophit_time;
  std::vector<std::vector<double>> _flash_ophit_risetime;
  std::vector<std::vector<double>> _flash_ophit_starttime;
  std::vector<std::vector<double>>_flash_ophit_amp;
  std::vector<std::vector<double>> _flash_ophit_area;
  std::vector<std::vector<double>> _flash_ophit_width;
  std::vector<std::vector<double>> _flash_ophit_pe;
  std::vector<std::vector<int>> _flash_ophit_ch;
  std::vector<TH1D* > _flash_time_distribution;


  // Maps to store the photon information per box
  std::set<int> fPDSBoxIDs;

  std::map<int, double> fBoxMap_PECoated;
  std::map<int, double> fBoxMap_PEUncoated;
  std::map<int, int> fBoxMap_NCoatedCh;
  std::map<int, int> fBoxMap_NUncoatedCh;

  // PMT Ratio variables
  std::vector<double> fPEUncoated;
  std::vector<double> fPECoated;
  std::vector<double> fPMTRatioPE;
  std::vector<double> fPMTRatioPEPerBox;


  // CRT track variables 
  std::vector<double>                _tr_start_x;
  std::vector<double>                _tr_start_y;
  std::vector<double>                _tr_start_z;
  std::vector<double>                _tr_end_x;
  std::vector<double>                _tr_end_y;
  std::vector<double>                _tr_end_z;
  std::vector<double>                _tr_dir_x;
  std::vector<double>                _tr_dir_y;
  std::vector<double>                _tr_dir_z;
  std::vector<double>                _tr_ts0;
  std::vector<double>                _tr_ets0;
  std::vector<double>                _tr_ts1;
  std::vector<double>                _tr_ets1;
  std::vector<double>                _tr_pe;
  std::vector<double>                _tr_length;
  std::vector<double>                _tr_tof;
  std::vector<double>                _tr_theta;
  std::vector<double>                _tr_phi;
  std::vector<bool>                  _tr_triple;
  std::vector<int16_t>               _tr_tagger1;
  std::vector<int16_t>               _tr_tagger2;
  std::vector<int16_t>               _tr_tagger3;

  // Slice variables
  int fNSlices;

  //Space Point Variables
  std::vector<double> fSpacePointX;
  std::vector<double> fSpacePointY;
  std::vector<double> fSpacePointZ;
  std::vector<double> fSpacePointIntegral;



  //Waveforms
  std::vector<std::vector<double>> fRawChannelADC;
  std::vector<int> fRawChannelID;
  std::vector<double> fRawChannelPedestal;

  //Recob Wires
  unsigned int fNROIs;
  std::vector<std::vector<float>> fWireADC;
  std::vector<unsigned int> fWireID;
  std::vector<int> fWireStampTime;

  //Reconstructed vertex
  double fRecoVx;
  double fRecoVy;
  double fRecoVz;
  int fRecoVU;
  int fRecoVV;
  int fRecoVC;
  int fRecoVTimeTick;
  bool fRecoInFV;
  
  // Charge barycenter variables
  std::vector<double> fChargeWeightX;
  std::vector<double> fChargeWeightY;
  std::vector<double> fChargeWeightZ;
  std::vector<double> fChargeTotalWeight;

  std::vector<double> fChargeBarycenterX;
  std::vector<double> fChargeBarycenterY;
  std::vector<double> fChargeBarycenterZ;


  // Reco track start/end points
  std::vector<std::vector<double>> fPFTrackStart;
  std::vector<std::vector<double>> fPFTrackEnd;
  std::vector<double> fPFPDGCode;
  std::vector<double> fPFPT0;

  int fNAnalyzedEvents;

  //const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  art::ServiceHandle<geo::Geometry> fGeom;
  geo::WireReadoutGeom const& channelMapAlg = art::ServiceHandle<geo::WireReadout const>()->Get();

  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;
  // Wire orientation
  double fCos60;
  double fSin60;
  double fWirePitch;

  opdet::sbndPDMapAlg fPDSMap;

};


void test::FlashMatchAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree=tfs->make<TTree>("AnaTPCTree", "Analysis Output Tree");
  fTree->Branch("RunID", &fRunID, "RunID/I");
  fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");

  if(fSaveTruth){
    fTree->Branch("TruePrimariesPDG", &fTruePrimariesPDG);
    fTree->Branch("TruePrimariesE", &fTruePrimariesE);
    fTree->Branch("TruePrimariesStartP", &fTruePrimariesStartP);
    fTree->Branch("TrueVx", &fTrueVx, "TrueVx/D");
    fTree->Branch("TrueVy", &fTrueVy, "TrueVy/D");
    fTree->Branch("TrueVz", &fTrueVz, "TrueVz/D");
    fTree->Branch("TrueVt", &fTrueVt, "TrueVt/D");
    fTree->Branch("TrueVU", &fTrueVU, "TrueVU/I");
    fTree->Branch("TrueVV", &fTrueVV, "TrueVV/I");
    fTree->Branch("TrueVC", &fTrueVC, "TrueVC/I");
    fTree->Branch("TrueVTimeTick", &fTrueVTimeTick, "TrueVC/I");
    fTree->Branch("TrueVEnergy", &fTrueVEnergy, "TrueVEnergy/D");
    fTree->Branch("IntMode", &fIntMode, "IntMode/I");
    fTree->Branch("IntCCNC", &fIntCCNC, "IntCCMC/I");
    fTree->Branch("IntNProtons", &fIntNProtons, "IntNProtons/I");
    fTree->Branch("IntNNeutrons", &fIntNNeutrons, "IntNNeutrons/I");
    fTree->Branch("IntNPi0", &fIntNPi0, "IntNPi0/I");
    fTree->Branch("IntNPip", &fIntNPip, "IntNPip/I");
    fTree->Branch("IntNPim", &fIntNPim, "IntNPim/I");
    fTree->Branch("IntNMuonP", &fIntNMuonP, "IntNMuonP/I");
    fTree->Branch("IntNMuonM", &fIntNMuonM, "IntNMuonM/I");
    fTree->Branch("IntNElectronP", &fIntNElectronP, "IntNElectronP/I");
    fTree->Branch("IntNElectronM", &fIntNElectronM, "IntNElectronM/I");
    fTree->Branch("IntNLambda", &fIntNLambda, "IntNLambda/I");
    fTree->Branch("IntInFV", &fIntInFV, "IntInFV/O");
  }

  if(fSaveSimED){
    fTree->Branch("EnDepE", &fEnDepE);
    fTree->Branch("EnDepX", &fEnDepX);
    fTree->Branch("EnDepY", &fEnDepY);
    fTree->Branch("EnDepZ", &fEnDepZ);
    fTree->Branch("EnDepU", &fEnDepU);
    fTree->Branch("EnDepV", &fEnDepV);
    fTree->Branch("EnDepC", &fEnDepC);
    fTree->Branch("EnDepT", &fEnDepT);
    fTree->Branch("EnDepPDG", &fEnDepPDG);
  }

  if(fSaveSimEDOut){
    fTree->Branch("EnDepEOut", &fEnDepEOut);
    fTree->Branch("EnDepXOut", &fEnDepXOut);
    fTree->Branch("EnDepYOut", &fEnDepYOut);
    fTree->Branch("EnDepZOut", &fEnDepZOut);
    fTree->Branch("EnDepTOut", &fEnDepTOut);
  }

  if(fSaveWaveforms){
    fTree->Branch("RawChannelADC", &fRawChannelADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  if(fSaveWires){
    fTree->Branch("NROIs", &fNROIs);
    fTree->Branch("WireID", &fWireID);
    fTree->Branch("WireStampTime", &fWireStampTime);
    fTree->Branch("WireADC", &fWireADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  if(fSaveWires){
    fTree->Branch("RawChannelADC", &fRawChannelADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  fTree->Branch("HitsViewReco1", &fHitsViewReco1);
  fTree->Branch("HitsIntegralReco1", &fHitsIntegralReco1);
  fTree->Branch("HitsSummedADCReco1", &fHitsSummedADCReco1);
  fTree->Branch("HitsPeakTimeReco1", &fHitsPeakTimeReco1);
  fTree->Branch("HitsChannelReco1", &fHitsChannelReco1);
  fTree->Branch("HitsAmplitudeReco1", &fHitsAmplitudeReco1);
  fTree->Branch("HitsRMSReco1", &fHitsRMSReco1);
  fTree->Branch("HitsStartTReco1", &fHitsStartTReco1);
  fTree->Branch("HitsEndTReco1", &fHitsEndTReco1);
  fTree->Branch("HitsWidthReco1", &fHitsWidthReco1);
  fTree->Branch("HitsChi2Reco1", &fHitsChi2Reco1);
  fTree->Branch("HitsNDFReco1", &fHitsNDFReco1);

  fTree->Branch("HitsView", &fHitsViewReco2);
  fTree->Branch("HitsSummedADC", &fHitsSummedADCReco2);
  fTree->Branch("HitsPeakTime", &fHitsPeakTimeReco2);
  fTree->Branch("HitsChannel", &fHitsChannelReco2);
  fTree->Branch("HitsAmplitude", &fHitsAmplitudeReco2);
  fTree->Branch("HitsRMS", &fHitsRMSReco2);
  fTree->Branch("HitsStartT", &fHitsStartTReco2);
  fTree->Branch("HitsEndT", &fHitsEndTReco2);
  fTree->Branch("HitsWidth", &fHitsWidthReco2);
  fTree->Branch("HitsChi2", &fHitsChi2Reco2);
  fTree->Branch("HitsNDF", &fHitsNDFReco2);
  fTree->Branch("HitsIntegral", &fHitsIntegralReco2);
  fTree->Branch("HitsX", &fHitsXReco2);
  fTree->Branch("HitsY", &fHitsYReco2);
  fTree->Branch("HitsZ", &fHitsZReco2);

  fTree->Branch("SpacePointX", &fSpacePointX);
  fTree->Branch("SpacePointY", &fSpacePointY);
  fTree->Branch("SpacePointZ", &fSpacePointZ);
  fTree->Branch("SpacePointIntegral", &fSpacePointIntegral);

  fTree->Branch("RecoVx", &fRecoVx, "RecoVx/D");
  fTree->Branch("RecoVy", &fRecoVy, "RecoVy/D");
  fTree->Branch("RecoVz", &fRecoVz, "RecoVz/D");
  fTree->Branch("RecoVU", &fRecoVU, "RecoVU/I");
  fTree->Branch("RecoVV", &fRecoVV, "RecoVV/I");
  fTree->Branch("RecoVC", &fRecoVC, "RecoVC/I");
  fTree->Branch("RecoVTimeTick", &fRecoVTimeTick, "RecoVC/I");
  fTree->Branch("RecoInFV", &fRecoInFV, "RecoInFV/O");

  fTree->Branch("ChargeBarycenterX", &fChargeBarycenterX);
  fTree->Branch("ChargeBarycenterY", &fChargeBarycenterY);
  fTree->Branch("ChargeBarycenterZ", &fChargeBarycenterZ);

  fTree->Branch("PFTrackStart", &fPFTrackStart);
  fTree->Branch("PFTrackEnd", &fPFTrackEnd);
  fTree->Branch("PFPDGCode", &fPFPDGCode);
  fTree->Branch("PFPT0", &fPFPT0);

  fOpAnaTree = tfs->make<TTree>("OpAnaTree", "PDS Analysis Output Tree");
  fOpAnaTree->Branch("RunID", &fRunID, "RunID/I");
  fOpAnaTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fOpAnaTree->Branch("EventID", &fEventID, "EventID/I");

  // OpFlashes
  if(fSaveOpFlashes){
    fOpAnaTree->Branch("nopflash", &_nopflash, "nopflash/I");
    fOpAnaTree->Branch("flash_id","std::vector<int>", &_flash_id);
    fOpAnaTree->Branch("flash_time","std::vector<double>", &_flash_time);
    fOpAnaTree->Branch("flash_total_pe", "std::vector<double>", &_flash_total_pe);
    fOpAnaTree->Branch("flash_pe_v","std::vector<std::vector<double>>", &_flash_pe_v);
    fOpAnaTree->Branch("flash_tpc", "std::vector<int>", &_flash_tpc);
    fOpAnaTree->Branch("flash_y","std::vector<double>", &_flash_y);
    fOpAnaTree->Branch("flash_yerr", "std::vector<double>", &_flash_yerr);
    fOpAnaTree->Branch("flash_z","std::vector<double>", &_flash_z);
    fOpAnaTree->Branch("flash_zerr", "std::vector<double>", &_flash_zerr);
    fOpAnaTree->Branch("flash_x","std::vector<double>", &_flash_x);
    fOpAnaTree->Branch("flash_xerr", "std::vector<double>", &_flash_xerr);
    fOpAnaTree->Branch("flash_time_distribution", "std::vector<TH1D*>" ,&_flash_time_distribution);
    fOpAnaTree->Branch("flash_ophit_time", "std::vector<std::vector<double>>", &_flash_ophit_time);
    fOpAnaTree->Branch("flash_ophit_risetime", "std::vector<std::vector<double>>", &_flash_ophit_risetime);
    fOpAnaTree->Branch("flash_ophit_starttime", "std::vector<std::vector<double>>", &_flash_ophit_starttime);
    fOpAnaTree->Branch("flash_ophit_amp", "std::vector<std::vector<double>>", &_flash_ophit_amp);
    fOpAnaTree->Branch("flash_ophit_area", "std::vector<std::vector<double>>", &_flash_ophit_area);
    fOpAnaTree->Branch("flash_ophit_width", "std::vector<std::vector<double>>", &_flash_ophit_width);
    fOpAnaTree->Branch("flash_ophit_pe", "std::vector<std::vector<double>>", &_flash_ophit_pe);
    fOpAnaTree->Branch("flash_ophit_ch", "std::vector<std::vector<int>>", &_flash_ophit_ch);
    fOpAnaTree->Branch("PECoated", "std::vector<double>", &fPECoated);
    fOpAnaTree->Branch("PEUncoated", "std::vector<double>", &fPEUncoated);
    fOpAnaTree->Branch("PMTRatioPE","std::vector<double>", &fPMTRatioPE);
    fOpAnaTree->Branch("PMTRatioPEPerBox","std::vector<double>", &fPMTRatioPEPerBox);
  }

  fCRTTree = tfs->make<TTree>("CRTTree", "CRT Analysis Output Tree");
  fCRTTree->Branch("RunID", &fRunID, "RunID/I");
  fCRTTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fCRTTree->Branch("EventID", &fEventID, "EventID/I");
  // CRT info
  if(fSaveCRTTracks)
  {
    fCRTTree->Branch("tr_start_x", "std::vector<double>", &_tr_start_x);
    fCRTTree->Branch("tr_start_y", "std::vector<double>", &_tr_start_y);
    fCRTTree->Branch("tr_start_z", "std::vector<double>", &_tr_start_z);
    fCRTTree->Branch("tr_end_x", "std::vector<double>", &_tr_end_x);
    fCRTTree->Branch("tr_end_y", "std::vector<double>", &_tr_end_y);
    fCRTTree->Branch("tr_end_z", "std::vector<double>", &_tr_end_z);
    fCRTTree->Branch("tr_dir_x", "std::vector<double>", &_tr_dir_x);
    fCRTTree->Branch("tr_dir_y", "std::vector<double>", &_tr_dir_y);
    fCRTTree->Branch("tr_dir_z", "std::vector<double>", &_tr_dir_z);
    fCRTTree->Branch("tr_ts0", "std::vector<double>", &_tr_ts0);
    fCRTTree->Branch("tr_ets0", "std::vector<double>", &_tr_ets0);
    fCRTTree->Branch("_trts1", "std::vector<double>", &_tr_ts1);
    fCRTTree->Branch("tr_ets1", "std::vector<double>", &_tr_ets1);
    fCRTTree->Branch("tr_pe", "std::vector<double>", &_tr_pe);
    fCRTTree->Branch("tr_length", "std::vector<double>", &_tr_length);
    fCRTTree->Branch("tr_tof", "std::vector<double>", &_tr_tof);
    fCRTTree->Branch("tr_theta", "std::vector<double>", &_tr_theta);
    fCRTTree->Branch("tr_phi", "std::vector<double>", &_tr_phi);
    fCRTTree->Branch("tr_triple", "std::vector<bool>", &_tr_triple);
    fCRTTree->Branch("tr_tagger1", "std::vector<int16_t>", &_tr_tagger1);
    fCRTTree->Branch("tr_tagger2", "std::vector<int16_t>", &_tr_tagger2);
    fCRTTree->Branch("tr_tagger3", "std::vector<int16_t>", &_tr_tagger3);
  }

  fNAnalyzedEvents=0;
}

void test::FlashMatchAnalyzer::endJob(){

  if(fCreateTPCMap){
    std::ofstream fileout("TPCMapping.txt");
    std::ofstream fileoutXYZ("TPCMappingXYZ.txt");
    if(fileout.is_open()){
        double xyz_start[3], xyz_end[3];;
        for(unsigned int ch=0; ch<channelMapAlg.Nchannels(); ch++){
          std::vector<geo::WireID> wireV = channelMapAlg.ChannelToWire(ch);
          for(size_t w=0; w<wireV.size(); w++){
            channelMapAlg.WireEndPoints(wireV[w], xyz_start, xyz_end);
            fileout<<ch<<" "<<wireV[w].Plane<<" "<<wireV[w].TPC<<std::endl;
            fileoutXYZ<<ch<<" "<<wireV[w].Plane<<" "<<wireV[w].TPC<<" ";
            fileoutXYZ<<xyz_start[0]<<" "<<xyz_start[1]<<" "<<xyz_start[2]<<" ";
            fileoutXYZ<<xyz_end[0]<<" "<<xyz_end[1]<<" "<<xyz_end[2]<<std::endl;
            std::cout<<ch<<":"<<w<<" ID="<<wireV[w].Wire<<" Plane="<<wireV[w].Plane<<" TPC="<<wireV[w].TPC<<std::endl;
          }
        }
    }
    fileout.close();
    fileoutXYZ.close();
  }

}

DEFINE_ART_MODULE(test::FlashMatchAnalyzer)