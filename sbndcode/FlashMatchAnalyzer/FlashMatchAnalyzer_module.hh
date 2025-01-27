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


#include "TTree.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

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
  std::string fCalorimetryLabel;
  std::string fParticleIDLabel;  
  std::vector<std::string> fOpFlashesModuleLabel;
  bool fSaveReco2;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fSaveSimEDOut;
  bool fSaveWaveforms;
  bool fSaveWires;
  bool fSaveHits;
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


  TTree* fTree;
  TTree* fOpAnaTree;
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

  //Hit variables
  std::vector<int> fHitsView;
  std::vector<double> fHitsPeakTime;
  std::vector<double> fHitsIntegral;
  std::vector<double> fHitsSummedADC;
  std::vector<double> fHitsChannel;
  std::vector<double> fHitsAmplitude;
  std::vector<double> fHitsRMS;
  std::vector<double> fHitsStartT;
  std::vector<double> fHitsEndT;
  std::vector<double> fHitsWidth;
  std::vector<double> fHitsChi2;
  std::vector<double> fHitsNDF;
  std::vector<int> fHitsClusterID;
  std::vector<double> fHitsX;
  std::vector<double> fHitsY;
  std::vector<double> fHitsZ;


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

  // Reco track start/end points
  std::vector<std::vector<double>> fPFTrackStart;
  std::vector<std::vector<double>> fPFTrackEnd;
  std::vector<double> fPFPDGCode;

  int fNAnalyzedEvents;

  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

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

  if(fSaveHits){
    fTree->Branch("HitsView", &fHitsView);
    fTree->Branch("HitsIntegral", &fHitsIntegral);
    fTree->Branch("HitsSummedADC", &fHitsSummedADC);
    fTree->Branch("HitsPeakTime", &fHitsPeakTime);
    fTree->Branch("HitsChannel", &fHitsChannel);
    fTree->Branch("HitsAmplitude", &fHitsAmplitude);
    fTree->Branch("HitsRMS", &fHitsRMS);
    fTree->Branch("HitsStartT", &fHitsStartT);
    fTree->Branch("HitsEndT", &fHitsEndT);
    fTree->Branch("HitsWidth", &fHitsWidth);
    fTree->Branch("HitsChi2", &fHitsChi2);
    fTree->Branch("HitsNDF", &fHitsNDF);
    fTree->Branch("HitsClusterID", &fHitsClusterID);
    fTree->Branch("HitsX", &fHitsX);
    fTree->Branch("HitsY", &fHitsY);
    fTree->Branch("HitsZ", &fHitsZ);
    
  }

  if(fSaveSpacePoints){
    fTree->Branch("SpacePointX", &fSpacePointX);
    fTree->Branch("SpacePointY", &fSpacePointY);
    fTree->Branch("SpacePointZ", &fSpacePointZ);
    fTree->Branch("SpacePointIntegral", &fSpacePointIntegral);
  }

  if(fSaveVertex){
    fTree->Branch("RecoVx", &fRecoVx, "RecoVx/D");
    fTree->Branch("RecoVy", &fRecoVy, "RecoVy/D");
    fTree->Branch("RecoVz", &fRecoVz, "RecoVz/D");
    fTree->Branch("RecoVU", &fRecoVU, "RecoVU/I");
    fTree->Branch("RecoVV", &fRecoVV, "RecoVV/I");
    fTree->Branch("RecoVC", &fRecoVC, "RecoVC/I");
    fTree->Branch("RecoVTimeTick", &fRecoVTimeTick, "RecoVC/I");
    fTree->Branch("RecoInFV", &fRecoInFV, "RecoInFV/O");
  }

  if(fSaveReco2){
    fTree->Branch("PFTrackStart", &fPFTrackStart);
    fTree->Branch("PFTrackEnd", &fPFTrackEnd);
    fTree->Branch("PFPDGCode", &fPFPDGCode);
  }

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
    fOpAnaTree->Branch("flash_ophit_time", "std::vector<std::vector<double>>", &_flash_ophit_time);
    fOpAnaTree->Branch("flash_ophit_risetime", "std::vector<std::vector<double>>", &_flash_ophit_risetime);
    fOpAnaTree->Branch("flash_ophit_starttime", "std::vector<std::vector<double>>", &_flash_ophit_starttime);
    fOpAnaTree->Branch("flash_ophit_amp", "std::vector<std::vector<double>>", &_flash_ophit_amp);
    fOpAnaTree->Branch("flash_ophit_area", "std::vector<std::vector<double>>", &_flash_ophit_area);
    fOpAnaTree->Branch("flash_ophit_width", "std::vector<std::vector<double>>", &_flash_ophit_width);
    fOpAnaTree->Branch("flash_ophit_pe", "std::vector<std::vector<double>>", &_flash_ophit_pe);
    fOpAnaTree->Branch("flash_ophit_ch", "std::vector<std::vector<int>>", &_flash_ophit_ch);
  }

  fNAnalyzedEvents=0;
}

void test::FlashMatchAnalyzer::endJob(){

  if(fCreateTPCMap){
    std::ofstream fileout("TPCMapping.txt");
    std::ofstream fileoutXYZ("TPCMappingXYZ.txt");
    if(fileout.is_open()){
        double xyz_start[3], xyz_end[3];;
        for(unsigned int ch=0; ch<fGeom->Nchannels(); ch++){
          std::vector<geo::WireID> wireV = fGeom->ChannelToWire(ch);
          for(size_t w=0; w<wireV.size(); w++){
            fGeom->WireEndPoints(wireV[w], xyz_start, xyz_end);
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