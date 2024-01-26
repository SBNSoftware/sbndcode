////////////////////////////////////////////////////////////////////////
// Class:       ToFAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        ToFAnalyzer_module.cc
//
// Generated at Tue Feb  1 09:34:47 2022 by Varuna Crishan Meddage using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "sbndcode/RecoUtils/RecoUtils.h"
#include "lardataobj/MCBase/MCShower.h"

// sbndcode includes
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"

#include "lardataobj/RecoBase/OpHit.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>


using namespace std;

namespace sbnd{

class ToFAnalyzer : public art::EDAnalyzer {
public:
  using CRTSpacePoint = sbnd::crt::CRTSpacePoint;
  explicit ToFAnalyzer(fhicl::ParameterSet const& p);
  ToFAnalyzer(ToFAnalyzer const&) = delete;
  ToFAnalyzer(ToFAnalyzer&&) = delete;
  ToFAnalyzer& operator=(ToFAnalyzer const&) = delete;
  ToFAnalyzer& operator=(ToFAnalyzer&&) = delete;

  void beginJob() override;
  void analyze(art::Event const& evt) override;
  bool SpacePointCompare(const art::Ptr<CRTSpacePoint>& sp1, const art::Ptr<CRTSpacePoint>& sp2);
  void ClearVecs();
  double length(const simb::MCParticle& part, TVector3& start, TVector3& end);

private:
  const double LAR_PROP_DELAY = 1.0/1.35e1; //[ns/cm] - from fcl dump on G4 fcl file (vuv_vgroup_mean: 1.35e1)
  
  art::InputTag fGenLabel;
  art::InputTag fSimLabel;
  art::InputTag fOpHitModuleLabel;
  art::InputTag fOpFlashModuleLabel0;
  art::InputTag fOpFlashModuleLabel1;
  art::InputTag fCrtSpacePointModuleLabel;
  art::InputTag fCrtTrackModuleLabel;

  double fCoinWindow;
  double fOpDelay;
  double fCRTSpacePointThresh;
  double fFlashPeThresh;
  double fHitPeThresh;
  double fBeamLow;
  double fBeamUp; 
  bool fLFlash;
  bool fLFlash_hit;
  bool fCFlash;
  bool fCFlash_hit;
  bool fLhit;
  bool fChit;
  bool fSaveNuInfo;
  bool fSaveG4Info;
  double fG4timeUp;
  double fG4timeLow;
  bool fkeeponlytracks;
  bool fSaveTrueToFInfo;
  
  sbnd::crt::CRTBackTrackerAlg* bt;
  map<int,art::InputTag> fFlashLabels;
  geo::GeometryCore const* fGeometryService;
  
  TTree* fMatchTree; 
  
  int    frun;                  
  int    fsubrun;               
  int    fevent;
  
  vector<int> fnu_pdg;
  vector<double> fnu_E;
  vector<int> fnu_mode; //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  vector<int> fnu_CCNC; //0=CC 1=NC
  vector<double> fnu_posX;
  vector<double> fnu_posY; 
  vector<double> fnu_posZ; 
  vector<double> fnu_T;
  vector<bool> fnu_CRYO;
  vector<bool> fnu_TPC; 
  
  vector<int> fg4_pdg;
  vector<int> fg4_trkid;
  vector<bool> fg4st_TPC;
  vector<bool> fg4en_TPC;
  vector<bool> fg4st_CRYO;
  vector<bool> fg4en_CRYO;
  vector<bool> fg4is_prim;
  vector<double> fg4_T0;
  vector<double> fg4_T0_CRYO;
  vector<double> fg4_T0_TPC;
  vector<int> fg4_org;
  vector<double> fg4_stX; 
  vector<double> fg4_stY;
  vector<double> fg4_stZ;
  vector<double> fg4_enX; 
  vector<double> fg4_enY;
  vector<double> fg4_enZ;
  vector<double> fg4_tlen;
  vector<double> fg4_E;
  vector<bool> fg4_is_in_TPC;
  vector<bool> fg4_is_in_CRYO; 
  
  vector<double> fLhit_tof_vec;
  vector<bool> fLhit_frmtrk_vec;
  vector<bool> fLhit_frmhit_vec;
  vector<int> fLhit_crtspkey_vec;
  vector<int> fLhit_crttrkkey_vec;
  vector<double> fLhit_crttime_t1_vec;
  //  vector<double> fLhit_crttime_t0_vec;
  vector<double> fLhit_crtpos_X_vec;
  vector<double> fLhit_crtpos_Y_vec;
  vector<double> fLhit_crtpos_Z_vec;
  vector<double> fLhit_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fLhit_crttgr_vec;
  vector<int> fLhit_pmthitkey_vec;
  vector<double> fLhit_pmthitT_vec;
  vector<double> fLhit_pmthitX_vec;
  vector<double> fLhit_pmthitY_vec;
  vector<double> fLhit_pmthitZ_vec;
  vector<double> fLhit_pmthitpe_vec;
  
  vector<double> fChit_tof_vec;
  vector<bool> fChit_frmtrk_vec;
  vector<bool> fChit_frmhit_vec;
  vector<int> fChit_crtspkey_vec;
  vector<int> fChit_crttrkkey_vec;
  vector<double> fChit_crttime_t1_vec;
  //  vector<double> fChit_crttime_t0_vec;
  vector<double> fChit_crtpos_X_vec;
  vector<double> fChit_crtpos_Y_vec;
  vector<double> fChit_crtpos_Z_vec;
  vector<double> fChit_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fChit_crttgr_vec;
  vector<int> fChit_pmthitkey_vec;
  vector<double> fChit_pmthitT_vec;
  vector<double> fChit_pmthitX_vec;
  vector<double> fChit_pmthitY_vec;
  vector<double> fChit_pmthitZ_vec;
  vector<double> fChit_pmthitpe_vec;
  
  vector<double> fLflsh_tof_vec;
  vector<bool> fLflsh_frmtrk_vec;
  vector<bool> fLflsh_frmhit_vec;
  vector<int> fLflsh_crtspkey_vec;
  vector<int> fLflsh_crttrkkey_vec;
  vector<double> fLflsh_crttime_t1_vec;
  //  vector<double> fLflsh_crttime_t0_vec;
  vector<double> fLflsh_crtpos_X_vec;
  vector<double> fLflsh_crtpos_Y_vec;
  vector<double> fLflsh_crtpos_Z_vec;
  vector<double> fLflsh_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fLflsh_crttgr_vec;
  vector<int> fLflsh_pmtflshtpcID_vec;
  vector<int> fLflsh_pmtflshkey_vec;
  vector<double> fLflsh_pmtflshT_vec;
  vector<double> fLflsh_pmtflshY_vec;
  vector<double> fLflsh_pmtflshZ_vec;
  vector<double> fLflsh_pmtflshpe_vec;
  
  vector<double> fCflsh_tof_vec;
  vector<bool> fCflsh_frmtrk_vec;
  vector<bool> fCflsh_frmhit_vec;
  vector<int> fCflsh_crtspkey_vec;
  vector<int> fCflsh_crttrkkey_vec;
  vector<double> fCflsh_crttime_t1_vec;
  //  vector<double> fCflsh_crttime_t0_vec;
  vector<double> fCflsh_crtpos_X_vec;
  vector<double> fCflsh_crtpos_Y_vec;
  vector<double> fCflsh_crtpos_Z_vec;
  vector<double> fCflsh_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fCflsh_crttgr_vec;
  vector<int> fCflsh_pmtflshtpcID_vec;
  vector<int> fCflsh_pmtflshkey_vec;
  vector<double> fCflsh_pmtflshT_vec;
  vector<double> fCflsh_pmtflshY_vec;
  vector<double> fCflsh_pmtflshZ_vec;
  vector<double> fCflsh_pmtflshpe_vec;
  
  vector<double> fLflshhit_tof_vec;
  vector<bool> fLflshhit_frmtrk_vec;
  vector<bool> fLflshhit_frmhit_vec;
  vector<int> fLflshhit_crtspkey_vec;
  vector<int> fLflshhit_crttrkkey_vec;
  vector<double> fLflshhit_crttime_t1_vec;
  //  vector<double> fLflshhit_crttime_t0_vec;
  vector<double> fLflshhit_crtpos_X_vec;
  vector<double> fLflshhit_crtpos_Y_vec;
  vector<double> fLflshhit_crtpos_Z_vec;
  vector<double> fLflshhit_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fLflshhit_crttgr_vec;
  vector<int> fLflshhit_pmtflshtpcID_vec;
  vector<int> fLflshhit_pmtflshkey_vec;
  vector<int> fLflshhit_pmtkey_vec;
  vector<double> fLflshhit_pmtflshT_vec;
  vector<double> fLflshhit_pmtflshX_vec;
  vector<double> fLflshhit_pmtflshY_vec;
  vector<double> fLflshhit_pmtflshZ_vec;
  vector<double> fLflshhit_pmtflshpe_vec;
  
  vector<double> fCflshhit_tof_vec;
  vector<bool> fCflshhit_frmtrk_vec;
  vector<bool> fCflshhit_frmhit_vec;
  vector<int> fCflshhit_crtspkey_vec;
  vector<int> fCflshhit_crttrkkey_vec;
  vector<double> fCflshhit_crttime_t1_vec;
  //  vector<double> fCflshhit_crttime_t0_vec;
  vector<double> fCflshhit_crtpos_X_vec;
  vector<double> fCflshhit_crtpos_Y_vec;
  vector<double> fCflshhit_crtpos_Z_vec;
  vector<double> fCflshhit_crtpe_vec;
  vector<sbnd::crt::CRTTagger> fCflshhit_crttgr_vec;
  vector<int> fCflshhit_pmtflshtpcID_vec;
  vector<int> fCflshhit_pmtflshkey_vec;
  vector<int> fCflshhit_pmtkey_vec;
  vector<double> fCflshhit_pmtflshT_vec;
  vector<double> fCflshhit_pmtflshX_vec;
  vector<double> fCflshhit_pmtflshY_vec;
  vector<double> fCflshhit_pmtflshZ_vec;
  vector<double> fCflshhit_pmtflshpe_vec;
  
  vector<double> fTrue_TOF;
  vector<bool> fTrue_TOF_hit; // from hit == true, from track == alse
  vector<int> fTrue_TOF_pdg;
  vector<int> fTrue_TOF_part_ID;
  vector<bool> fTrue_TOF_part_prim;
  vector<bool> fTrue_TOF_part_org;
  vector<double> fTrue_TOF_traj_T;
  vector<double> fTrue_TOF_traj_X;
  vector<double> fTrue_TOF_traj_Y;
  vector<double> fTrue_TOF_traj_Z;
  vector<bool> fTrue_TOF_traj_in_TPC; // if trajectory point is in TPC = true, else false
};

ToFAnalyzer::ToFAnalyzer(fhicl::ParameterSet const& pset): 
EDAnalyzer{pset},
fGenLabel(pset.get<art::InputTag>("GenLabel")),
fSimLabel(pset.get<art::InputTag>("SimLabel")),
fOpHitModuleLabel(pset.get<art::InputTag>("OpHitModuleLabel")),				
fOpFlashModuleLabel0(pset.get<art::InputTag>("OpFlashModuleLabel0")),
fOpFlashModuleLabel1(pset.get<art::InputTag>("OpFlashModuleLabel1")),
fCrtSpacePointModuleLabel(pset.get<art::InputTag>("CrtSpacePointModuleLabel")),
fCrtTrackModuleLabel(pset.get<art::InputTag>("CrtTrackModuleLabel")),
fCoinWindow(pset.get<double>("CoincidenceWindow")),
fOpDelay(pset.get<double>("OpDelay")), // the cable time delay (135 ns) + PMTTransit Time delay (55 ns) + average TPB emission time delay (~4.5 ns)
fCRTSpacePointThresh(pset.get<double>("CRTSpacePointThresh")),
fFlashPeThresh(pset.get<int>("FlashPeThresh")),
fHitPeThresh(pset.get<int>("HitPeThresh")),
fBeamLow(pset.get<double>("BeamLow")),
fBeamUp(pset.get<double>("BeamUp")),
fLFlash(pset.get<bool>("LFlash")),
fLFlash_hit(pset.get<bool>("LFlash_hit")),
fCFlash(pset.get<bool>("CFlash")),
fCFlash_hit(pset.get<bool>("CFlash_hit")),
fLhit(pset.get<bool>("Lhit")),
fChit(pset.get<bool>("Chit")),
fSaveNuInfo(pset.get<bool>("SaveNuInfo")),
fSaveG4Info(pset.get<bool>("SaveG4Info")),
fG4timeUp(pset.get<double>("G4timeUp")),
fG4timeLow(pset.get<double>("G4timeLow")),
fkeeponlytracks(pset.get<bool>("keeponlytracks")),
fSaveTrueToFInfo(pset.get<bool>("SaveTrueToFInfo")),
  bt(new sbnd::crt::CRTBackTrackerAlg(pset.get<fhicl::ParameterSet>("CRTBackTrackerAlg")))
{
fFlashLabels[0] = fOpFlashModuleLabel0;
fFlashLabels[1] = fOpFlashModuleLabel1; 
}

//========================================================================
void ToFAnalyzer::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  
  fGeometryService = lar::providerFrom<geo::Geometry>();  
  
  art::ServiceHandle<art::TFileService> tfs;
  fMatchTree = tfs->make<TTree>("matchTree","CRTSpacePoint - OpHit/Flash matching analysis");
  
  fMatchTree->Branch("run", &frun, "run/I");
  fMatchTree->Branch("subrun", &fsubrun, "subrun/I");
  fMatchTree->Branch("event", &fevent, "event/I");
  
  fMatchTree->Branch("nu_pdg", &fnu_pdg);
  fMatchTree->Branch("nu_E", &fnu_E);
  fMatchTree->Branch("nu_mode", &fnu_mode);
  fMatchTree->Branch("nu_CCNC", &fnu_CCNC);
  fMatchTree->Branch("nu_posX", &fnu_posX);
  fMatchTree->Branch("nu_posY", &fnu_posY);
  fMatchTree->Branch("nu_posZ", &fnu_posZ);
  fMatchTree->Branch("nu_T", &fnu_T);
  fMatchTree->Branch("nu_CRYO", &fnu_CRYO);
  fMatchTree->Branch("nu_TPC", &fnu_TPC);
  
  fMatchTree->Branch("g4_pdg", &fg4_pdg);
  fMatchTree->Branch("g4_trkid", &fg4_trkid);
  fMatchTree->Branch("g4st_TPC", &fg4st_TPC);
  fMatchTree->Branch("g4en_TPC", &fg4en_TPC);
  fMatchTree->Branch("g4st_CRYO", &fg4st_CRYO);
  fMatchTree->Branch("g4en_CRYO", &fg4en_CRYO);
  fMatchTree->Branch("g4is_prim", &fg4is_prim);
  fMatchTree->Branch("g4_T0", &fg4_T0);
  fMatchTree->Branch("g4_T0_CRYO", &fg4_T0_CRYO);
  fMatchTree->Branch("g4_T0_TPC", &fg4_T0_TPC);
  fMatchTree->Branch("g4_org", &fg4_org);
  fMatchTree->Branch("g4_stX", &fg4_stX);
  fMatchTree->Branch("g4_stY", &fg4_stY);
  fMatchTree->Branch("g4_stZ", &fg4_stZ);
  fMatchTree->Branch("g4_enX", &fg4_enX);
  fMatchTree->Branch("g4_enY", &fg4_enY);
  fMatchTree->Branch("g4_enZ", &fg4_enZ);
  fMatchTree->Branch("g4_tlen", &fg4_tlen);
  fMatchTree->Branch("g4_E", &fg4_E);
  fMatchTree->Branch("g4_is_in_TPC", &fg4_is_in_TPC);
  fMatchTree->Branch("g4_is_in_CRYO", &fg4_is_in_CRYO);
  
  fMatchTree->Branch("Lhit_tof_vec", &fLhit_tof_vec);
  fMatchTree->Branch("Lhit_frmtrk_vec", &fLhit_frmtrk_vec);
  fMatchTree->Branch("Lhit_frmhit_vec", &fLhit_frmhit_vec);
  fMatchTree->Branch("Lhit_crtspkey_vec", &fLhit_crtspkey_vec);
  fMatchTree->Branch("Lhit_crttrkkey_vec", &fLhit_crttrkkey_vec);
  fMatchTree->Branch("Lhit_crttime_t1_vec", &fLhit_crttime_t1_vec);
  //  fMatchTree->Branch("Lhit_crttime_t0_vec", &fLhit_crttime_t0_vec);
  fMatchTree->Branch("Lhit_crtpos_X_vec", &fLhit_crtpos_X_vec);
  fMatchTree->Branch("Lhit_crtpos_Y_vec", &fLhit_crtpos_Y_vec);
  fMatchTree->Branch("Lhit_crtpos_Z_vec", &fLhit_crtpos_Z_vec);
  fMatchTree->Branch("Lhit_crtpe_vec", &fLhit_crtpe_vec);
  fMatchTree->Branch("Lhit_crttgr_vec", &fLhit_crttgr_vec);
  fMatchTree->Branch("Lhit_pmthitkey_vec", &fLhit_pmthitkey_vec);
  fMatchTree->Branch("Lhit_pmthitT_vec", &fLhit_pmthitT_vec);
  fMatchTree->Branch("Lhit_pmthitX_vec", &fLhit_pmthitX_vec);
  fMatchTree->Branch("Lhit_pmthitY_vec", &fLhit_pmthitY_vec);
  fMatchTree->Branch("Lhit_pmthitZ_vec", &fLhit_pmthitZ_vec);
  fMatchTree->Branch("Lhit_pmthitpe_vec", &fLhit_pmthitpe_vec);
  
  fMatchTree->Branch("Chit_tof_vec", &fChit_tof_vec);
  fMatchTree->Branch("Chit_frmtrk_vec", &fChit_frmtrk_vec);
  fMatchTree->Branch("Chit_frmhit_vec", &fChit_frmhit_vec);
  fMatchTree->Branch("Chit_crtspkey_vec", &fChit_crtspkey_vec);
  fMatchTree->Branch("Chit_crttrkkey_vec", &fChit_crttrkkey_vec);
  fMatchTree->Branch("Chit_crttime_t1_vec", &fChit_crttime_t1_vec);
  //  fMatchTree->Branch("Chit_crttime_t0_vec", &fChit_crttime_t0_vec);
  fMatchTree->Branch("Chit_crtpos_X_vec", &fChit_crtpos_X_vec);
  fMatchTree->Branch("Chit_crtpos_Y_vec", &fChit_crtpos_Y_vec);
  fMatchTree->Branch("Chit_crtpos_Z_vec", &fChit_crtpos_Z_vec);
  fMatchTree->Branch("Chit_crtpe_vec", &fChit_crtpe_vec);
  fMatchTree->Branch("Chit_crttgr_vec", &fChit_crttgr_vec);
  fMatchTree->Branch("Chit_pmthitkey_vec", &fChit_pmthitkey_vec);
  fMatchTree->Branch("Chit_pmthitT_vec", &fChit_pmthitT_vec);
  fMatchTree->Branch("Chit_pmthitX_vec", &fChit_pmthitX_vec);
  fMatchTree->Branch("Chit_pmthitY_vec", &fChit_pmthitY_vec);
  fMatchTree->Branch("Chit_pmthitZ_vec", &fChit_pmthitZ_vec);
  fMatchTree->Branch("Chit_pmthitpe_vec", &fChit_pmthitpe_vec);
  
  fMatchTree->Branch("Lflsh_tof_vec", &fLflsh_tof_vec);
  fMatchTree->Branch("Lflsh_frmtrk_vec", &fLflsh_frmtrk_vec);
  fMatchTree->Branch("Lflsh_frmhit_vec", &fLflsh_frmhit_vec);
  fMatchTree->Branch("Lflsh_crtspkey_vec", &fLflsh_crtspkey_vec);
  fMatchTree->Branch("Lflsh_crttrkkey_vec", &fLflsh_crttrkkey_vec);
  fMatchTree->Branch("Lflsh_crttime_t1_vec", &fLflsh_crttime_t1_vec);
  //  fMatchTree->Branch("Lflsh_crttime_t0_vec", &fLflsh_crttime_t0_vec);
  fMatchTree->Branch("Lflsh_crtpos_X_vec", &fLflsh_crtpos_X_vec);
  fMatchTree->Branch("Lflsh_crtpos_Y_vec", &fLflsh_crtpos_Y_vec);
  fMatchTree->Branch("Lflsh_crtpos_Z_vec", &fLflsh_crtpos_Z_vec);
  fMatchTree->Branch("Lflsh_crtpe_vec", &fLflsh_crtpe_vec);
  fMatchTree->Branch("Lflsh_crttgr_vec", &fLflsh_crttgr_vec);
  fMatchTree->Branch("Lflsh_pmtflshtpcID_vec", &fLflsh_pmtflshtpcID_vec);
  fMatchTree->Branch("Lflsh_pmtflshkey_vec", &fLflsh_pmtflshkey_vec);
  fMatchTree->Branch("Lflsh_pmtflshT_vec", &fLflsh_pmtflshT_vec);
  fMatchTree->Branch("Lflsh_pmtflshY_vec", &fLflsh_pmtflshY_vec);
  fMatchTree->Branch("Lflsh_pmtflshZ_vec", &fLflsh_pmtflshZ_vec);
  fMatchTree->Branch("Lflsh_pmtflshpe_vec", &fLflsh_pmtflshpe_vec);
  
  fMatchTree->Branch("Cflsh_tof_vec", &fCflsh_tof_vec);
  fMatchTree->Branch("Cflsh_frmtrk_vec", &fCflsh_frmtrk_vec);
  fMatchTree->Branch("Cflsh_frmhit_vec", &fCflsh_frmhit_vec);
  fMatchTree->Branch("Cflsh_crtspkey_vec", &fCflsh_crtspkey_vec);
  fMatchTree->Branch("Cflsh_crttrkkey_vec", &fCflsh_crttrkkey_vec);
  fMatchTree->Branch("Cflsh_crttime_t1_vec", &fCflsh_crttime_t1_vec);
  //  fMatchTree->Branch("Cflsh_crttime_t0_vec", &fCflsh_crttime_t0_vec);
  fMatchTree->Branch("Cflsh_crtpos_X_vec", &fCflsh_crtpos_X_vec);
  fMatchTree->Branch("Cflsh_crtpos_Y_vec", &fCflsh_crtpos_Y_vec);
  fMatchTree->Branch("Cflsh_crtpos_Z_vec", &fCflsh_crtpos_Z_vec);
  fMatchTree->Branch("Cflsh_crtpe_vec", &fCflsh_crtpe_vec);
  fMatchTree->Branch("Cflsh_crttgr_vec", &fCflsh_crttgr_vec);
  fMatchTree->Branch("Cflsh_pmtflshtpcID_vec", &fCflsh_pmtflshtpcID_vec);
  fMatchTree->Branch("Cflsh_pmtflshkey_vec", &fCflsh_pmtflshkey_vec);
  fMatchTree->Branch("Cflsh_pmtflshT_vec", &fCflsh_pmtflshT_vec);
  fMatchTree->Branch("Cflsh_pmtflshY_vec", &fCflsh_pmtflshY_vec);
  fMatchTree->Branch("Cflsh_pmtflshZ_vec", &fCflsh_pmtflshZ_vec);
  fMatchTree->Branch("Cflsh_pmtflshpe_vec", &fCflsh_pmtflshpe_vec);
  
  fMatchTree->Branch("Lflshhit_tof_vec", &fLflshhit_tof_vec);
  fMatchTree->Branch("Lflshhit_frmtrk_vec", &fLflshhit_frmtrk_vec);
  fMatchTree->Branch("Lflshhit_frmhit_vec", &fLflshhit_frmhit_vec);
  fMatchTree->Branch("Lflshhit_crtspkey_vec", &fLflshhit_crtspkey_vec);
  fMatchTree->Branch("Lflshhit_crttrkkey_vec", &fLflshhit_crttrkkey_vec);
  fMatchTree->Branch("Lflshhit_crttime_t1_vec", &fLflshhit_crttime_t1_vec);
  //  fMatchTree->Branch("Lflshhit_crttime_t0_vec", &fLflshhit_crttime_t0_vec);
  fMatchTree->Branch("Lflshhit_crtpos_X_vec", &fLflshhit_crtpos_X_vec);
  fMatchTree->Branch("Lflshhit_crtpos_Y_vec", &fLflshhit_crtpos_Y_vec);
  fMatchTree->Branch("Lflshhit_crtpos_Z_vec", &fLflshhit_crtpos_Z_vec);
  fMatchTree->Branch("Lflshhit_crtpe_vec", &fLflshhit_crtpe_vec);
  fMatchTree->Branch("Lflshhit_crttgr_vec", &fLflshhit_crttgr_vec);
  fMatchTree->Branch("Lflshhit_pmtflshtpcID_vec", &fLflshhit_pmtflshtpcID_vec);
  fMatchTree->Branch("Lflshhit_pmtflshkey_vec", &fLflshhit_pmtflshkey_vec);
  fMatchTree->Branch("Lflshhit_pmtkey_vec", &fLflshhit_pmtkey_vec);
  fMatchTree->Branch("Lflshhit_pmtflshT_vec", &fLflshhit_pmtflshT_vec);
  fMatchTree->Branch("Lflshhit_pmtflshX_vec", &fLflshhit_pmtflshX_vec);
  fMatchTree->Branch("Lflshhit_pmtflshY_vec", &fLflshhit_pmtflshY_vec);
  fMatchTree->Branch("Lflshhit_pmtflshZ_vec", &fLflshhit_pmtflshZ_vec);
  fMatchTree->Branch("Lflshhit_pmtflshpe_vec", &fLflshhit_pmtflshpe_vec);
  
  fMatchTree->Branch("Cflshhit_tof_vec", &fCflshhit_tof_vec);
  fMatchTree->Branch("Cflshhit_frmtrk_vec", &fCflshhit_frmtrk_vec);
  fMatchTree->Branch("Cflshhit_frmhit_vec", &fCflshhit_frmhit_vec);
  fMatchTree->Branch("Cflshhit_crtspkey_vec", &fCflshhit_crtspkey_vec);
  fMatchTree->Branch("Cflshhit_crttrkkey_vec", &fCflshhit_crttrkkey_vec);
  fMatchTree->Branch("Cflshhit_crttime_t1_vec", &fCflshhit_crttime_t1_vec);
  //  fMatchTree->Branch("Cflshhit_crttime_t0_vec", &fCflshhit_crttime_t0_vec);
  fMatchTree->Branch("Cflshhit_crtpos_X_vec", &fCflshhit_crtpos_X_vec);
  fMatchTree->Branch("Cflshhit_crtpos_Y_vec", &fCflshhit_crtpos_Y_vec);
  fMatchTree->Branch("Cflshhit_crtpos_Z_vec", &fCflshhit_crtpos_Z_vec);
  fMatchTree->Branch("Cflshhit_crtpe_vec", &fCflshhit_crtpe_vec);
  fMatchTree->Branch("Cflshhit_crttgr_vec", &fCflshhit_crttgr_vec);
  fMatchTree->Branch("Cflshhit_pmtflshtpcID_vec", &fCflshhit_pmtflshtpcID_vec);
  fMatchTree->Branch("Cflshhit_pmtflshkey_vec", &fCflshhit_pmtflshkey_vec);
  fMatchTree->Branch("Cflshhit_pmtkey_vec", &fCflshhit_pmtkey_vec);
  fMatchTree->Branch("Cflshhit_pmtflshT_vec", &fCflshhit_pmtflshT_vec);
  fMatchTree->Branch("Cflshhit_pmtflshX_vec", &fCflshhit_pmtflshX_vec);
  fMatchTree->Branch("Cflshhit_pmtflshY_vec", &fCflshhit_pmtflshY_vec);
  fMatchTree->Branch("Cflshhit_pmtflshZ_vec", &fCflshhit_pmtflshZ_vec);
  fMatchTree->Branch("Cflshhit_pmtflshpe_vec", &fCflshhit_pmtflshpe_vec);
  
  fMatchTree->Branch("True_TOF", &fTrue_TOF);
  fMatchTree->Branch("True_TOF_hit", &fTrue_TOF_hit);
  fMatchTree->Branch("True_TOF_pdg", &fTrue_TOF_pdg);
  fMatchTree->Branch("True_TOF_part_ID", &fTrue_TOF_part_ID);
  fMatchTree->Branch("True_TOF_part_prim", &fTrue_TOF_part_prim);
  fMatchTree->Branch("True_TOF_part_org", &fTrue_TOF_part_org);
  fMatchTree->Branch("True_TOF_traj_T", &fTrue_TOF_traj_T);
  fMatchTree->Branch("True_TOF_traj_X", &fTrue_TOF_traj_X);
  fMatchTree->Branch("True_TOF_traj_Y", &fTrue_TOF_traj_Y);
  fMatchTree->Branch("True_TOF_traj_Z", &fTrue_TOF_traj_Z);
  fMatchTree->Branch("True_TOF_traj_in_TPC", &fTrue_TOF_traj_in_TPC);
}



void ToFAnalyzer::analyze(art::Event const& evt)
{
 ClearVecs();
 
 geo::CryostatGeo const& cryo0 = fGeometryService->Cryostat();
 geo::TPCGeo const& tpc00 = cryo0.TPC(0);
 geo::TPCGeo const& tpc01 = cryo0.TPC(1);
 
 frun=evt.run();
 fsubrun=evt.subRun();
 fevent=evt.id().event(); 
 
 //std::cout << "******************  " <<  frun << "  " << fsubrun << "  " << fevent << "  *****************\n";
 
 //===================Neutrino information =========================
 
 if(fSaveNuInfo){
    art::Handle< std::vector<simb::MCTruth> > MCTruthListHandle;
    std::vector< art::Ptr<simb::MCTruth> >    MCTruthkList;
    if( evt.getByLabel(fGenLabel,MCTruthListHandle))
     art::fill_ptr_vector(MCTruthkList, MCTruthListHandle);
    
    for(auto const& mctruth : MCTruthkList){
        if (!mctruth->NeutrinoSet()) continue;
	auto const& nu = mctruth->GetNeutrino();
	const TLorentzVector xyzt = nu.Nu().Position(0);
	fnu_pdg.push_back(nu.Nu().PdgCode());
	fnu_E.push_back(nu.Nu().E());
	fnu_mode.push_back(nu.Mode());
	fnu_CCNC.push_back(nu.CCNC());
	fnu_posX.push_back(xyzt.X());
	fnu_posY.push_back(xyzt.Y());
	fnu_posZ.push_back(xyzt.Z());
	fnu_T.push_back(xyzt.T());
	double point[3] = {fnu_posX.back(), fnu_posY.back(), fnu_posZ.back()};
	if(cryo0.ContainsPosition(point)) fnu_CRYO.push_back(true);
	else fnu_CRYO.push_back(false);
	if(tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point)) fnu_TPC.push_back(true);
	else fnu_TPC.push_back(false);
    }
 }
	 
//===================================================================
 
//=============G4 level particle information ========================

 if(fSaveG4Info){
    art::Handle< std::vector<simb::MCParticle> > SimPartListHandle;
    std::vector< art::Ptr<simb::MCParticle> >    SimPartList;
    if( evt.getByLabel(fSimLabel,SimPartListHandle))
     art::fill_ptr_vector(SimPartList, SimPartListHandle);
    
    for(auto const& part : SimPartList){
        const cheat::ParticleInventory *inventory_service=lar::providerFrom<cheat::ParticleInventoryService>();
	if(!(part->T()<=fG4timeUp && part->T()>=fG4timeLow)) continue;
	if(fkeeponlytracks){
	   if(!(TMath::Abs(part->PdgCode())==13 || TMath::Abs(part->PdgCode())==211 || TMath::Abs(part->PdgCode())==2212 || TMath::Abs(part->PdgCode())==321)) continue;
        }
	
	fg4_pdg.push_back(part->PdgCode());
	fg4_trkid.push_back(part->TrackId());
        if(part->Process()=="primary") fg4is_prim.push_back(true);
	else fg4is_prim.push_back(false);
	fg4_T0.push_back(part->T());
	art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(part->TrackId());
	fg4_org.push_back(truth->NeutrinoSet());
	//std::cout << "******* Total number of trajectory points : " << part->NumberTrajectoryPoints() <<"\n";
	const TLorentzVector st_xyzt = part->Position(0);
	double st_point[3] = {st_xyzt.X(), st_xyzt.Y(), st_xyzt.Z()};
	fg4_stX.push_back(st_point[0]); 
        fg4_stY.push_back(st_point[1]);
        fg4_stZ.push_back(st_point[2]);
	if(cryo0.ContainsPosition(st_point)) fg4st_CRYO.push_back(true);
	else fg4st_CRYO.push_back(false);
	if(tpc00.ContainsPosition(st_point) || tpc01.ContainsPosition(st_point)) fg4st_TPC.push_back(true);
	else fg4st_TPC.push_back(false);
	const TLorentzVector en_xyzt = part->EndPosition();
	double en_point[3] = {en_xyzt.X(), en_xyzt.Y(), en_xyzt.Z()};
	fg4_enX.push_back(en_point[0]); 
        fg4_enY.push_back(en_point[1]);
        fg4_enZ.push_back(en_point[2]);
	if(cryo0.ContainsPosition(en_point)) fg4en_CRYO.push_back(true);
	else fg4en_CRYO.push_back(false);
	if(tpc00.ContainsPosition(en_point) || tpc01.ContainsPosition(en_point)) fg4en_TPC.push_back(true);
	else fg4en_TPC.push_back(false);
        fg4_E.push_back(part->E());
	TVector3 mcstart, mcend;
        fg4_tlen.push_back(length(*part, mcstart, mcend));
	bool is_in_TPC=false;
	bool is_in_CRYO=false;
	bool is_filled=false;
	fg4_T0_CRYO.push_back(-9999);
	fg4_T0_TPC.push_back(-9999);
        for(size_t i=0; i<part->NumberTrajectoryPoints(); i++){
	    const TLorentzVector& pos = part->Position(i);
            double point[3] = {pos.X(),pos.Y(),pos.Z()};
	    //std::cout << "Trjectory point : " << i << "  " << pos.X() << "  " << pos.Y() << "  " << pos.Z() << "\n";
            if(cryo0.ContainsPosition(point) && !is_filled){ 
	       is_filled=true;	   
	       is_in_CRYO=true;
	       fg4_T0_CRYO.back()=pos.T();
            }
		   
	    if(tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point)){ 
	       is_in_TPC=true; 
	       fg4_T0_TPC.back()=pos.T(); 
	       break;
	    }
        }
	fg4_is_in_TPC.push_back(is_in_TPC);
	fg4_is_in_CRYO.push_back(is_in_CRYO);
     }
 } 
 
 //==================================================================
 
 art::Handle< std::vector<sbnd::crt::CRTTrack> > crtTrackListHandle;
 std::vector< art::Ptr<sbnd::crt::CRTTrack> >    crtTrackList;
 if( evt.getByLabel(fCrtTrackModuleLabel,crtTrackListHandle))
     art::fill_ptr_vector(crtTrackList, crtTrackListHandle);
 
 art::FindManyP<CRTSpacePoint> findManySPs(crtTrackListHandle, evt, fCrtTrackModuleLabel); 
 std::vector<std::vector<art::Ptr<CRTSpacePoint>>> tracksps;
 
 //================================================================
 
 for(size_t itrk=0; itrk<crtTrackList.size(); itrk++){
     std::vector<art::Ptr<CRTSpacePoint>> trksps = findManySPs.at(itrk);
     std::sort(trksps.begin(),trksps.end(),
     [](const art::Ptr<CRTSpacePoint>& a, const art::Ptr<CRTSpacePoint>& b)->bool
     { 
       return a->Time() < b->Time(); 
     });
     tracksps.push_back(trksps);
 } // Crt space points coming from are ordered on ascending order by looking into ts1_ns variable
 
 //==================================================================
 
 art::Handle< std::vector<CRTSpacePoint> > crtSPListHandle;
 std::vector< art::Ptr<CRTSpacePoint> >    crtSPList;
 if( evt.getByLabel(fCrtSpacePointModuleLabel,crtSPListHandle))
     art::fill_ptr_vector(crtSPList, crtSPListHandle);

 art::FindOneP<sbnd::crt::CRTCluster> spToCluster(crtSPListHandle, evt, fCrtSpacePointModuleLabel);
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Lhit_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpHit>> > Lhit_tof_op_hits;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Chit_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpHit>> > Chit_tof_op_hits;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Lflsh_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpFlash>> > Lflsh_tof_op_flashes;
 map<int, std::vector<int>> Lflsh_tof_op_tpc;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Cflsh_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpFlash>> > Cflsh_tof_op_flashes;
 map<int, std::vector<int>> Cflsh_tof_op_tpc;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Lflshhit_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpHit>> > Lflshhit_tof_op_hits;
 map<int, std::vector<art::Ptr<recob::OpFlash>> > Lflshhit_tof_op_flashes;
 map<int, std::vector<int>> Lflshhit_tof_op_tpc;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > Cflshhit_tof_crt_sps;
 map<int, std::vector<art::Ptr<recob::OpHit>> > Cflshhit_tof_op_hits;
 map<int, std::vector<art::Ptr<recob::OpFlash>> > Cflshhit_tof_op_flashes;
 map<int, std::vector<int>> Cflshhit_tof_op_tpc;
 
 map<int, std::vector<art::Ptr<CRTSpacePoint>> > True_tof_crt_sps;
 map<int,std::vector<const simb::MCParticle*>> True_tof_sim_particles;
 
 for(auto const& crt : crtSPList){

   art::Ptr<sbnd::crt::CRTCluster> cluster = spToCluster.at(crt.key());

     if(!(crt->Time() >= fBeamLow &&  crt->Time()<= fBeamUp)) continue;
     if(crt->PE() < fCRTSpacePointThresh) continue;
     
     bool frm_trk=false;
     int index=0;
     
     for(auto const& trksps: tracksps){
	 for(size_t isp=0; isp<trksps.size(); isp++){
	     if(SpacePointCompare(trksps[isp],crt)){
		frm_trk=true;
		break;
	     }
	  }
	  if(frm_trk) break;
	  index++;
     }
     
     //======================== Doing a truth level study of ToF======================================
     
     if(fSaveTrueToFInfo){
	const cheat::ParticleInventory *inventory_service=lar::providerFrom<cheat::ParticleInventoryService>();
	sbnd::crt::CRTBackTrackerAlg::TruthMatchMetrics truthMatch=bt->TruthMatching(evt, cluster);
	int trackID = truthMatch.trackid;

	auto const& simparticles = *evt.getValidHandle<vector<simb::MCParticle>>(fSimLabel);
        map<int,const simb::MCParticle*> particleMap;
	for(auto const& particle : simparticles) particleMap[particle.TrackId()] = &particle;
	
	if(particleMap.find(abs(trackID))!=particleMap.end()){
	   if(frm_trk){
	      True_tof_crt_sps[index].push_back(crt);
	      True_tof_sim_particles[index].push_back(particleMap[abs(trackID)]);
	   }
	   else{
	       auto const& particle=particleMap[abs(trackID)];
	       bool found_tpc_traj_point=false;
	       for(size_t i=0; i<particle->NumberTrajectoryPoints(); i++){
	           const TLorentzVector& pos = particle->Position(i);
                   geo::Point_t const point{pos.X(),pos.Y(),pos.Z()};
		   if(tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point)){
                     auto const opDetPos = cryo0.OpDet(cryo0.GetClosestOpDet(point)).GetCenter();
                     double dprop=(opDetPos - point).R();
		     double tprop=pos.T() + dprop*LAR_PROP_DELAY;
		     fTrue_TOF.push_back(crt->Time()-tprop);
		     fTrue_TOF_hit.push_back(true);
		     fTrue_TOF_pdg.push_back(particle->PdgCode());
		     fTrue_TOF_part_ID.push_back(particle->TrackId());
		     if(particle->Process()=="primary") fTrue_TOF_part_prim.push_back(true);
		     else fTrue_TOF_part_prim.push_back(false);
		     art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(particle->TrackId());
		     fTrue_TOF_part_org.push_back(truth->NeutrinoSet());
		     fTrue_TOF_traj_T.push_back(pos.T());
                     fTrue_TOF_traj_X.push_back(point.X());
                     fTrue_TOF_traj_Y.push_back(point.Y());
                     fTrue_TOF_traj_Z.push_back(point.Z());
		     fTrue_TOF_traj_in_TPC.push_back(true);
		     found_tpc_traj_point=true;
		     break;
		   } // trajectory point is in TPC
	       } // loop over trajectory points
	       
	       if(!found_tpc_traj_point){
	          for(size_t i=0; i<particle->NumberTrajectoryPoints(); i++){
		      const TLorentzVector& pos = particle->Position(i);
                      geo::Point_t const point{pos.X(),pos.Y(),pos.Z()};
		      if(cryo0.ContainsPosition(point)){
                        auto const opDetPos = cryo0.OpDet(cryo0.GetClosestOpDet(point)).GetCenter();
                        double dprop=(opDetPos - point).R();
		         double tprop=pos.T() + dprop*LAR_PROP_DELAY;
		         fTrue_TOF.push_back(crt->Time()-tprop);
		         fTrue_TOF_hit.push_back(true);
		         fTrue_TOF_pdg.push_back(particle->PdgCode());
		         fTrue_TOF_part_ID.push_back(particle->TrackId());
		         if(particle->Process()=="primary") fTrue_TOF_part_prim.push_back(true);
		         else fTrue_TOF_part_prim.push_back(false);
		         art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(particle->TrackId());
		         fTrue_TOF_part_org.push_back(truth->NeutrinoSet());
		         fTrue_TOF_traj_T.push_back(pos.T());
                         fTrue_TOF_traj_X.push_back(point.X());
                         fTrue_TOF_traj_Y.push_back(point.Y());
                         fTrue_TOF_traj_Z.push_back(point.Z());
		         fTrue_TOF_traj_in_TPC.push_back(false);
		         break;
		      } // found a trajectory point inside CRYO
		  } // loop over trajectory points
	       } // did not find a trajectory point inside TPC
	       
	   } // lonely CRT hit
	   
	} // CRT hit is truth matched
     } // Save true ToF information
     
     //===============================================================================================
     
     //============================== Calculating ToF using largest optical hit ======================
     
     if(fLhit){
        double pehit_max=0;
	bool found_tof = false;
	int ophit_index = -1;
	
	art::Handle< std::vector<recob::OpHit> > opHitListHandle;
        std::vector< art::Ptr<recob::OpHit> >    opHitList;
        if( evt.getByLabel(fOpHitModuleLabel,opHitListHandle) )
            art::fill_ptr_vector(opHitList, opHitListHandle);
	
	for(auto const& hit : opHitList){
	    if(hit->PE()<fHitPeThresh) continue;
	    double thit = hit->PeakTime()*1e3-fOpDelay;
	    
	    if(abs(crt->Time()-thit)<fCoinWindow && hit->PE()>pehit_max){
	       pehit_max = hit->PE();
	       ophit_index = hit.key();
	       found_tof = true;
	    }
	} // loop over optical hit list
	
	if(found_tof){
	   if(frm_trk){
	      Lhit_tof_crt_sps[index].push_back(crt);
	      Lhit_tof_op_hits[index].push_back(opHitList[ophit_index]);
	   }
	   
	   else{
               fLhit_tof_vec.push_back(crt->Time() - (opHitList[ophit_index]->PeakTime()*1e3-fOpDelay));
	       fLhit_frmtrk_vec.push_back(false);
	       fLhit_frmhit_vec.push_back(true);
	       fLhit_crtspkey_vec.push_back(crt.key());
	       fLhit_crttrkkey_vec.push_back(-9999);
	       fLhit_crttime_t1_vec.push_back(crt->Time());
	       //	       fLhit_crttime_t0_vec.push_back(crt->ts0_ns);
	       fLhit_crtpos_X_vec.push_back(crt->X());
	       fLhit_crtpos_Y_vec.push_back(crt->Y());
	       fLhit_crtpos_Z_vec.push_back(crt->Z());
	       fLhit_crtpe_vec.push_back(crt->PE());
	       fLhit_crttgr_vec.push_back(cluster->Tagger());
	       fLhit_pmthitkey_vec.push_back(ophit_index);
	       fLhit_pmthitT_vec.push_back(opHitList[ophit_index]->PeakTime()*1e3-fOpDelay);
               auto const pos = fGeometryService->OpDetGeoFromOpChannel(opHitList[ophit_index]->OpChannel()).GetCenter();
               fLhit_pmthitX_vec.push_back(pos.X());
               fLhit_pmthitY_vec.push_back(pos.Y());
               fLhit_pmthitZ_vec.push_back(pos.Z());
	       fLhit_pmthitpe_vec.push_back(opHitList[ophit_index]->PE());
	       //std::cout << "CRT time : " << crt->ts1_ns << "  CRT PE : " << crt->peshit << "   CRT Tagger : " << cluster->Tagger() << "\n"; 
	   }
	}
     }
     
     //===============================================================================================
     
     //============================== Calculating ToF using closest optical hit ======================
     
     if(fChit){
        double ophit_minTOF = DBL_MAX;
	bool found_tof = false;
	int ophit_index = -1;
	
	art::Handle< std::vector<recob::OpHit> > opHitListHandle;
        std::vector< art::Ptr<recob::OpHit> >    opHitList;
        if( evt.getByLabel(fOpHitModuleLabel,opHitListHandle) )
            art::fill_ptr_vector(opHitList, opHitListHandle);
	
	for(auto const& hit : opHitList){
	    if(hit->PE()<fHitPeThresh) continue;
	    double thit = hit->PeakTime()*1e3-fOpDelay;
	    
	    if(abs(crt->Time()-thit)<fCoinWindow && abs(crt->Time()-thit)<ophit_minTOF){
	       ophit_minTOF = abs(crt->Time()-thit);
	       ophit_index = hit.key();
	       found_tof = true;
	    }
	} // loop over optical hit list
	
	if(found_tof){
	   if(frm_trk){
	      Chit_tof_crt_sps[index].push_back(crt);
	      Chit_tof_op_hits[index].push_back(opHitList[ophit_index]);
	   }
	   
	   else{
	       fChit_tof_vec.push_back(crt->Time() - (opHitList[ophit_index]->PeakTime()*1e3-fOpDelay));
	       fChit_frmtrk_vec.push_back(false);
	       fChit_frmhit_vec.push_back(true);
	       fChit_crtspkey_vec.push_back(crt.key());
	       fChit_crttrkkey_vec.push_back(-9999);
	       fChit_crttime_t1_vec.push_back(crt->Time());
	       //	       fChit_crttime_t0_vec.push_back(crt->ts0_ns);
	       fChit_crtpos_X_vec.push_back(crt->X());
	       fChit_crtpos_Y_vec.push_back(crt->Y());
	       fChit_crtpos_Z_vec.push_back(crt->Z());
	       fChit_crtpe_vec.push_back(crt->PE());
	       fChit_crttgr_vec.push_back(cluster->Tagger());
	       fChit_pmthitkey_vec.push_back(ophit_index);
	       fChit_pmthitT_vec.push_back(opHitList[ophit_index]->PeakTime()*1e3-fOpDelay);
               auto const pos = fGeometryService->OpDetGeoFromOpChannel(opHitList[ophit_index]->OpChannel()).GetCenter();
               fChit_pmthitX_vec.push_back(pos.X());
               fChit_pmthitY_vec.push_back(pos.Y());
               fChit_pmthitZ_vec.push_back(pos.Z());
	       fChit_pmthitpe_vec.push_back(opHitList[ophit_index]->PE());
	   }
	}
     }
     
     //===============================================================================================
     
     //============================== Calculating ToF using largest flash ============================
     
     if(fLFlash){
        double peflash_max=0;
	bool found_tof = false;
	int opflash_index = -1;
	int flash_tpc = -1;
	
	std::map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
        std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
	
	for(int i=0; i<2; i++) {
            if( evt.getByLabel(fFlashLabels[i],flashHandles[i]) )
                art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
        }
	
	for(auto const& flashList : opFlashLists){
	    for(size_t iflash=0; iflash<flashList.second.size(); iflash++){
	        auto const& flash = flashList.second[iflash];
		if(flash->TotalPE()<fFlashPeThresh) continue;
		double tflash = flash->AbsTime()*1e3-fOpDelay;
		if(abs(crt->Time()-tflash)<fCoinWindow && flash->TotalPE()>peflash_max){
		   peflash_max=flash->TotalPE();
		   opflash_index = flash.key();	
		   found_tof = true;
		   flash_tpc = flashList.first;
		} 
	    } 
	}
	
	if(found_tof){
	   if(frm_trk){
	      Lflsh_tof_crt_sps[index].push_back(crt);
	      Lflsh_tof_op_flashes[index].push_back(opFlashLists[flash_tpc][opflash_index]);
	      Lflsh_tof_op_tpc[index].push_back(flash_tpc);
	   }
	   else{
	       fLflsh_tof_vec.push_back(crt->Time() - (opFlashLists[flash_tpc][opflash_index]->Time()*1e3-fOpDelay));
	       fLflsh_frmtrk_vec.push_back(false);
	       fLflsh_frmhit_vec.push_back(true);
	       fLflsh_crtspkey_vec.push_back(crt.key());
	       fLflsh_crttrkkey_vec.push_back(-9999);
	       fLflsh_crttime_t1_vec.push_back(crt->Time());
	       //	       fLflsh_crttime_t0_vec.push_back(crt->ts0_ns);
	       fLflsh_crtpos_X_vec.push_back(crt->X());
	       fLflsh_crtpos_Y_vec.push_back(crt->Y());
	       fLflsh_crtpos_Z_vec.push_back(crt->Z());
	       fLflsh_crtpe_vec.push_back(crt->PE());
	       fLflsh_crttgr_vec.push_back(cluster->Tagger());
	       fLflsh_pmtflshtpcID_vec.push_back(flash_tpc);
	       fLflsh_pmtflshkey_vec.push_back(opflash_index);
	       fLflsh_pmtflshT_vec.push_back(opFlashLists[flash_tpc][opflash_index]->Time()*1e3-fOpDelay);
	       fLflsh_pmtflshY_vec.push_back(opFlashLists[flash_tpc][opflash_index]->YCenter());
	       fLflsh_pmtflshZ_vec.push_back(opFlashLists[flash_tpc][opflash_index]->ZCenter());
	       fLflsh_pmtflshpe_vec.push_back(opFlashLists[flash_tpc][opflash_index]->TotalPE());
	   }
	} 
     }
     
     //===============================================================================================
     
     //============================== Calculating ToF using closest flash ============================
     
     if(fCFlash){
        double flash_minTOF = DBL_MAX;
	bool found_tof = false;
	int opflash_index = -1;
	int flash_tpc = -1;
	
	std::map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
        std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
	
	for(int i=0; i<2; i++) {
            if( evt.getByLabel(fFlashLabels[i],flashHandles[i]) )
                art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
        }
	
	for(auto const& flashList : opFlashLists){
	    for(size_t iflash=0; iflash<flashList.second.size(); iflash++){
	        auto const& flash = flashList.second[iflash];
		if(flash->TotalPE()<fFlashPeThresh) continue;
		double tflash = flash->Time()*1e3-fOpDelay;
		if(abs(crt->Time()-tflash)<fCoinWindow && abs(crt->Time()-tflash)<flash_minTOF){
		   flash_minTOF= abs(crt->Time()-tflash);
		   opflash_index = flash.key();	
		   found_tof = true;
		   flash_tpc = flashList.first;
		} 
	    } 
	}
	
	if(found_tof){
	   if(frm_trk){
	      Cflsh_tof_crt_sps[index].push_back(crt);
	      Cflsh_tof_op_flashes[index].push_back(opFlashLists[flash_tpc][opflash_index]);
	      Cflsh_tof_op_tpc[index].push_back(flash_tpc);
	   }
	   else{
	       fCflsh_tof_vec.push_back(crt->Time() - (opFlashLists[flash_tpc][opflash_index]->Time()*1e3-fOpDelay));
	       fCflsh_frmtrk_vec.push_back(false);
	       fCflsh_frmhit_vec.push_back(true);
	       fCflsh_crtspkey_vec.push_back(crt.key());
	       fCflsh_crttrkkey_vec.push_back(-9999);
	       fCflsh_crttime_t1_vec.push_back(crt->Time());
	       //	       fCflsh_crttime_t0_vec.push_back(crt->ts0_ns);
	       fCflsh_crtpos_X_vec.push_back(crt->X());
	       fCflsh_crtpos_Y_vec.push_back(crt->Y());
	       fCflsh_crtpos_Z_vec.push_back(crt->Z());
	       fCflsh_crtpe_vec.push_back(crt->PE());
	       fCflsh_crttgr_vec.push_back(cluster->Tagger());
	       fCflsh_pmtflshtpcID_vec.push_back(flash_tpc);
	       fCflsh_pmtflshkey_vec.push_back(opflash_index);
	       fCflsh_pmtflshT_vec.push_back(opFlashLists[flash_tpc][opflash_index]->Time()*1e3-fOpDelay);
	       fCflsh_pmtflshY_vec.push_back(opFlashLists[flash_tpc][opflash_index]->YCenter());
	       fCflsh_pmtflshZ_vec.push_back(opFlashLists[flash_tpc][opflash_index]->ZCenter());
	       fCflsh_pmtflshpe_vec.push_back(opFlashLists[flash_tpc][opflash_index]->TotalPE());
	   }
	} 
     }
     
     //===============================================================================================
     
     //================ Calculating ToF using earliest hit of the largest flash ======================
     
     if(fLFlash_hit){
        double peflash_max=0;
	bool found_tof = false;
	int opflash_index = -1;
	int flash_tpc = -1;
	int ophit_index = -1;
	
	std::map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
        std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
	
	for(int i=0; i<2; i++) {
            if( evt.getByLabel(fFlashLabels[i],flashHandles[i]) )
                art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
        }
	
	art::Handle< std::vector<recob::OpHit> > opHitListHandle;
        std::vector< art::Ptr<recob::OpHit> >    opHitList;
        if( evt.getByLabel(fOpHitModuleLabel,opHitListHandle) )
            art::fill_ptr_vector(opHitList, opHitListHandle);
	
	for(auto const& flashList : opFlashLists){	
	    art::FindManyP<recob::OpHit> findManyOpHits(flashHandles[flashList.first], evt, fFlashLabels[flashList.first]);	
	    for(size_t iflash=0; iflash<flashList.second.size(); iflash++){
	        auto const& flash = flashList.second[iflash];
		if(flash->TotalPE()<fFlashPeThresh) continue;
		double tflash = flash->AbsTime()*1e3-fOpDelay;
		if(abs(crt->Time()-tflash)<fCoinWindow && flash->TotalPE()>peflash_max){
		   peflash_max=flash->TotalPE();
		   opflash_index = flash.key();	
		   found_tof = true;
		   flash_tpc = flashList.first;
		   vector<art::Ptr<recob::OpHit>> hits = findManyOpHits.at(flash.key());
		   double flashMinHitT = DBL_MAX;
		   for(auto const& hit : hits){
		       double tPmt = hit->PeakTime()*1e3-fOpDelay; 
		       if(tPmt < flashMinHitT){
			  flashMinHitT = tPmt;
			  ophit_index =  hit.key();    
		       } // getting the earliest hit
		   } // loop over associated ophits of the flash
		} // with in conincidence window and getting the largest flash
	    } // loop over flash list
	} // loop over flash list map
	
	if(found_tof){
	   if(frm_trk){	    
	       Lflshhit_tof_crt_sps[index].push_back(crt);
	       Lflshhit_tof_op_flashes[index].push_back(opFlashLists[flash_tpc][opflash_index]);
	       Lflshhit_tof_op_tpc[index].push_back(flash_tpc);
	       Lflshhit_tof_op_hits[index].push_back(opHitList[ophit_index]);
	    }
	    
	    else{
	         fLflshhit_tof_vec.push_back(crt->Time() - (opHitList[ophit_index]->PeakTime()*1e3-fOpDelay));
		 fLflshhit_frmtrk_vec.push_back(false);
                 fLflshhit_frmhit_vec.push_back(true);   
		 fLflshhit_crtspkey_vec.push_back(crt.key());
                 fLflshhit_crttrkkey_vec.push_back(-9999);
                 fLflshhit_crttime_t1_vec.push_back(crt->Time());
		 //                 fLflshhit_crttime_t0_vec.push_back(crt->ts0_ns);
                 fLflshhit_crtpos_X_vec.push_back(crt->X());
                 fLflshhit_crtpos_Y_vec.push_back(crt->Y());
                 fLflshhit_crtpos_Z_vec.push_back(crt->Z());
                 fLflshhit_crtpe_vec.push_back(crt->PE());
		 fLflshhit_crttgr_vec.push_back(cluster->Tagger());
		 fLflshhit_pmtflshtpcID_vec.push_back(flash_tpc);
                 fLflshhit_pmtflshkey_vec.push_back(opflash_index);
                 fLflshhit_pmtkey_vec.push_back(ophit_index);
                 fLflshhit_pmtflshT_vec.push_back(opHitList[ophit_index]->PeakTime()*1e3-fOpDelay);
                 auto const pos = fGeometryService->OpDetGeoFromOpChannel(opHitList[ophit_index]->OpChannel()).GetCenter();
                 fLflshhit_pmtflshX_vec.push_back(pos.X());
                 fLflshhit_pmtflshY_vec.push_back(pos.Y());
                 fLflshhit_pmtflshZ_vec.push_back(pos.Z());
                 fLflshhit_pmtflshpe_vec.push_back(opHitList[ophit_index]->PE());
            }
	}
     }
     
     //===============================================================================================
     
     //================ Calculating ToF using earliest hit of the closest flash ======================
     
     if(fCFlash_hit){
        double flash_minTOF = DBL_MAX;
	bool found_tof = false;
	int opflash_index = -1;
	int flash_tpc = -1;
	int ophit_index = -1;
	
	std::map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
        std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
	
	for(int i=0; i<2; i++) {
            if( evt.getByLabel(fFlashLabels[i],flashHandles[i]) )
                art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
        }
	
	art::Handle< std::vector<recob::OpHit> > opHitListHandle;
        std::vector< art::Ptr<recob::OpHit> >    opHitList;
        if( evt.getByLabel(fOpHitModuleLabel,opHitListHandle) )
            art::fill_ptr_vector(opHitList, opHitListHandle);
	
	for(auto const& flashList : opFlashLists){	
	    art::FindManyP<recob::OpHit> findManyOpHits(flashHandles[flashList.first], evt, fFlashLabels[flashList.first]);	
	    for(size_t iflash=0; iflash<flashList.second.size(); iflash++){
	        auto const& flash = flashList.second[iflash];
		if(flash->TotalPE()<fFlashPeThresh) continue;
		double tflash = flash->AbsTime()*1e3-fOpDelay;
		if(abs(crt->Time()-tflash)<fCoinWindow && abs(crt->Time()-tflash)<flash_minTOF){
		   flash_minTOF= abs(crt->Time()-tflash);
		   opflash_index = flash.key();	
		   found_tof = true;
		   flash_tpc = flashList.first;
		   vector<art::Ptr<recob::OpHit>> hits = findManyOpHits.at(flash.key());
		   double flashMinHitT = DBL_MAX;
		   for(auto const& hit : hits){
		       double tPmt = hit->PeakTime()*1e3-fOpDelay; 
		       if(tPmt < flashMinHitT){
			  flashMinHitT = tPmt;
			  ophit_index =  hit.key();    
		       } // getting the earliest hit
		   } // loop over associated ophits of the flash
		} // with in conincidence window and getting the largest flash
	    } // loop over flash list
	} // loop over flash list map
	
	if(found_tof){
	   if(frm_trk){	    
	       Cflshhit_tof_crt_sps[index].push_back(crt);
	       Cflshhit_tof_op_flashes[index].push_back(opFlashLists[flash_tpc][opflash_index]);
	       Cflshhit_tof_op_tpc[index].push_back(flash_tpc);
	       Cflshhit_tof_op_hits[index].push_back(opHitList[ophit_index]);
	    }
	    
	    else{
	         fCflshhit_tof_vec.push_back(crt->Time() - (opHitList[ophit_index]->PeakTime()*1e3-fOpDelay));
		 fCflshhit_frmtrk_vec.push_back(false);
                 fCflshhit_frmhit_vec.push_back(true);   
		 fCflshhit_crtspkey_vec.push_back(crt.key());
                 fCflshhit_crttrkkey_vec.push_back(-9999);
                 fCflshhit_crttime_t1_vec.push_back(crt->Time());
		 //                 fCflshhit_crttime_t0_vec.push_back(crt->ts0_ns);
                 fCflshhit_crtpos_X_vec.push_back(crt->X());
                 fCflshhit_crtpos_Y_vec.push_back(crt->Y());
                 fCflshhit_crtpos_Z_vec.push_back(crt->Z());
                 fCflshhit_crtpe_vec.push_back(crt->PE());
		 fCflshhit_crttgr_vec.push_back(cluster->Tagger());
		 fCflshhit_pmtflshtpcID_vec.push_back(flash_tpc);
                 fCflshhit_pmtflshkey_vec.push_back(opflash_index);
                 fCflshhit_pmtkey_vec.push_back(ophit_index);
                 fCflshhit_pmtflshT_vec.push_back(opHitList[ophit_index]->PeakTime()*1e3-fOpDelay);
                 auto const pos = fGeometryService->OpDetGeoFromOpChannel(opHitList[ophit_index]->OpChannel()).GetCenter();
                 fCflshhit_pmtflshX_vec.push_back(pos.X());
                 fCflshhit_pmtflshY_vec.push_back(pos.Y());
                 fCflshhit_pmtflshZ_vec.push_back(pos.Z());
                 fCflshhit_pmtflshpe_vec.push_back(opHitList[ophit_index]->PE());
            }
	}
     }
     
     //===============================================================================================
     
 } // loop over CRT hits
 
 
 //============================== Calculating ToF using largest optical hit =========================
 
 if(fLhit){
    if(!Lhit_tof_crt_sps.empty()){
       for (auto& ele: Lhit_tof_crt_sps){
	    double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
            for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    } 
	    
	    fLhit_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Lhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay));
	    fLhit_frmtrk_vec.push_back(true);
	    fLhit_frmhit_vec.push_back(false);
	    fLhit_crtspkey_vec.push_back(tracksps[ele.first].front().key());
	    fLhit_crttrkkey_vec.push_back(ele.first);
	    fLhit_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //	    fLhit_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
	    fLhit_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
	    fLhit_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
	    fLhit_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
	    fLhit_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fLhit_crttgr_vec.push_back(clstr->Tagger());
	    fLhit_pmthitkey_vec.push_back(Lhit_tof_op_hits[ele.first][min_index].key());
	    fLhit_pmthitT_vec.push_back(Lhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay);
            auto const pos = fGeometryService->OpDetGeoFromOpChannel(Lhit_tof_op_hits[ele.first][min_index]->OpChannel()).GetCenter();
            fLhit_pmthitX_vec.push_back(pos.X());
            fLhit_pmthitY_vec.push_back(pos.Y());
            fLhit_pmthitZ_vec.push_back(pos.Z());
	    fLhit_pmthitpe_vec.push_back(Lhit_tof_op_hits[ele.first][min_index]->PE());
	}
    }
 }
 
 //==================================================================================================
 
 //============================== Calculating ToF using closest optical hit =========================
 
 if(fChit){
    if(!Chit_tof_crt_sps.empty()){
       for (auto& ele: Chit_tof_crt_sps){
	    double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
            for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    } 
	    
	    fChit_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Chit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay));
	    fChit_frmtrk_vec.push_back(true);
	    fChit_frmhit_vec.push_back(false);
	    fChit_crtspkey_vec.push_back(tracksps[ele.first].front().key());
	    fChit_crttrkkey_vec.push_back(ele.first);
	    fChit_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //	    fChit_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
	    fChit_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
	    fChit_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
	    fChit_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
	    fChit_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fChit_crttgr_vec.push_back(clstr->Tagger());
	    fChit_pmthitkey_vec.push_back(Chit_tof_op_hits[ele.first][min_index].key());
	    fChit_pmthitT_vec.push_back(Chit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay);
            auto const pos = fGeometryService->OpDetGeoFromOpChannel(Chit_tof_op_hits[ele.first][min_index]->OpChannel()).GetCenter();
            fChit_pmthitX_vec.push_back(pos.X());
            fChit_pmthitY_vec.push_back(pos.Y());
            fChit_pmthitZ_vec.push_back(pos.Z());
	    fChit_pmthitpe_vec.push_back(Chit_tof_op_hits[ele.first][min_index]->PE());
	}
    }
 }
 
 //==================================================================================================
 
 //============================== Calculating ToF using largest optical flash =======================
 
 if(fLFlash){
    if(!Lflsh_tof_crt_sps.empty()){
       for (auto& ele: Lflsh_tof_crt_sps){
            double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
	    for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    }
	    
	    fLflsh_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Lflsh_tof_op_flashes[ele.first][min_index]->Time()*1e3-fOpDelay));
	    fLflsh_frmtrk_vec.push_back(true);
	    fLflsh_frmhit_vec.push_back(false);
	    fLflsh_crtspkey_vec.push_back(tracksps[ele.first].front().key());
	    fLflsh_crttrkkey_vec.push_back(ele.first);
	    fLflsh_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //	    fLflsh_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
	    fLflsh_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
	    fLflsh_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
	    fLflsh_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
	    fLflsh_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fLflsh_crttgr_vec.push_back(clstr->Tagger());
	    fLflsh_pmtflshtpcID_vec.push_back(Lflsh_tof_op_tpc[ele.first][min_index]);
	    fLflsh_pmtflshkey_vec.push_back(Lflsh_tof_op_flashes[ele.first][min_index].key());
	    fLflsh_pmtflshT_vec.push_back(Lflsh_tof_op_flashes[ele.first][min_index]->Time()*1e3-fOpDelay);
	    fLflsh_pmtflshY_vec.push_back(Lflsh_tof_op_flashes[ele.first][min_index]->YCenter());
	    fLflsh_pmtflshZ_vec.push_back(Lflsh_tof_op_flashes[ele.first][min_index]->ZCenter());   
	    fLflsh_pmtflshpe_vec.push_back(Lflsh_tof_op_flashes[ele.first][min_index]->TotalPE());
       }
    }
 }
 
 //==================================================================================================
 
 //============================== Calculating ToF using closest optical flash =======================
 
 if(fCFlash){
    if(!Cflsh_tof_crt_sps.empty()){
       for (auto& ele: Cflsh_tof_crt_sps){
            double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
	    for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    }
	    
	    fCflsh_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Cflsh_tof_op_flashes[ele.first][min_index]->Time()*1e3-fOpDelay));
	    fCflsh_frmtrk_vec.push_back(true);
	    fCflsh_frmhit_vec.push_back(false);
	    fCflsh_crtspkey_vec.push_back(tracksps[ele.first].front().key());
	    fCflsh_crttrkkey_vec.push_back(ele.first);
	    fCflsh_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //	    fCflsh_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
	    fCflsh_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
	    fCflsh_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
	    fCflsh_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
	    fCflsh_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fCflsh_crttgr_vec.push_back(clstr->Tagger());
	    fCflsh_pmtflshtpcID_vec.push_back(Cflsh_tof_op_tpc[ele.first][min_index]);
	    fCflsh_pmtflshkey_vec.push_back(Cflsh_tof_op_flashes[ele.first][min_index].key());
	    fCflsh_pmtflshT_vec.push_back(Cflsh_tof_op_flashes[ele.first][min_index]->Time()*1e3-fOpDelay);
	    fCflsh_pmtflshY_vec.push_back(Cflsh_tof_op_flashes[ele.first][min_index]->YCenter());
	    fCflsh_pmtflshZ_vec.push_back(Cflsh_tof_op_flashes[ele.first][min_index]->ZCenter());   
	    fCflsh_pmtflshpe_vec.push_back(Cflsh_tof_op_flashes[ele.first][min_index]->TotalPE());
       }
    }
 }
 
 //==================================================================================================
 
 //=============Calculation ToF values using Earliest hit of the Largest flash =====================
 
 if(fLFlash_hit){
    if(!Lflshhit_tof_crt_sps.empty()){
       for (auto& ele: Lflshhit_tof_crt_sps){
            double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
	    for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    }
	    
	    fLflshhit_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Lflshhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay));
            fLflshhit_frmtrk_vec.push_back(true);
            fLflshhit_frmhit_vec.push_back(false);   
            fLflshhit_crtspkey_vec.push_back(tracksps[ele.first].front().key());
            fLflshhit_crttrkkey_vec.push_back(ele.first);
	    fLflshhit_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //            fLflshhit_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
            fLflshhit_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
            fLflshhit_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
            fLflshhit_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
            fLflshhit_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fLflshhit_crttgr_vec.push_back(clstr->Tagger());
            fLflshhit_pmtflshtpcID_vec.push_back(Lflshhit_tof_op_tpc[ele.first][min_index]);
            fLflshhit_pmtflshkey_vec.push_back(Lflshhit_tof_op_flashes[ele.first][min_index].key());
            fLflshhit_pmtkey_vec.push_back(Lflshhit_tof_op_hits[ele.first][min_index].key());
            fLflshhit_pmtflshT_vec.push_back(Lflshhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay);
            auto const pos = fGeometryService->OpDetGeoFromOpChannel(Lflshhit_tof_op_hits[ele.first][min_index]->OpChannel()).GetCenter();
            fLflshhit_pmtflshX_vec.push_back(pos.X());
            fLflshhit_pmtflshY_vec.push_back(pos.Y());
            fLflshhit_pmtflshZ_vec.push_back(pos.Z());
            fLflshhit_pmtflshpe_vec.push_back(Lflshhit_tof_op_hits[ele.first][min_index]->PE());
       }
    }
 } 
 
 //===================================================================================================
 
 //=============Calculation ToF values using Earliest hit of the Closest flash =======================
 
 if(fCFlash_hit){
    if(!Cflshhit_tof_crt_sps.empty()){
       for (auto& ele: Cflshhit_tof_crt_sps){
            double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
	    for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    }
	    
	    fCflshhit_tof_vec.push_back(tracksps[ele.first].front()->Time() - (Cflshhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay));
            fCflshhit_frmtrk_vec.push_back(true);
            fCflshhit_frmhit_vec.push_back(false);   
            fCflshhit_crtspkey_vec.push_back(tracksps[ele.first].front().key());
            fCflshhit_crttrkkey_vec.push_back(ele.first);
            fCflshhit_crttime_t1_vec.push_back(tracksps[ele.first].front()->Time());
	    //            fCflshhit_crttime_t0_vec.push_back(tracksps[ele.first].front()->ts0_ns);
            fCflshhit_crtpos_X_vec.push_back(tracksps[ele.first].front()->X());
            fCflshhit_crtpos_Y_vec.push_back(tracksps[ele.first].front()->Y());
            fCflshhit_crtpos_Z_vec.push_back(tracksps[ele.first].front()->Z());
            fCflshhit_crtpe_vec.push_back(tracksps[ele.first].front()->PE());
	    const art::Ptr<sbnd::crt::CRTCluster> clstr = spToCluster.at(tracksps[ele.first].front().key());
	    fCflshhit_crttgr_vec.push_back(clstr->Tagger());
            fCflshhit_pmtflshtpcID_vec.push_back(Cflshhit_tof_op_tpc[ele.first][min_index]);
            fCflshhit_pmtflshkey_vec.push_back(Cflshhit_tof_op_flashes[ele.first][min_index].key());
            fCflshhit_pmtkey_vec.push_back(Cflshhit_tof_op_hits[ele.first][min_index].key());
            fCflshhit_pmtflshT_vec.push_back(Cflshhit_tof_op_hits[ele.first][min_index]->PeakTime()*1e3-fOpDelay);
            auto const pos = fGeometryService->OpDetGeoFromOpChannel(Cflshhit_tof_op_hits[ele.first][min_index]->OpChannel()).GetCenter();
            fCflshhit_pmtflshX_vec.push_back(pos.X());
            fCflshhit_pmtflshY_vec.push_back(pos.Y());
            fCflshhit_pmtflshZ_vec.push_back(pos.Z());
            fCflshhit_pmtflshpe_vec.push_back(Cflshhit_tof_op_hits[ele.first][min_index]->PE());
       }
    }
 } 
 
 //===================================================================================================
 
 //=============Calculation ToF values uisng truth level information ================================
 
 if(fSaveTrueToFInfo){
    if(!True_tof_crt_sps.empty()){
       const cheat::ParticleInventory *inventory_service=lar::providerFrom<cheat::ParticleInventoryService>();
       for (auto& ele: True_tof_crt_sps){
            double min_time = DBL_MAX; 
	    int all_index = 0;
	    int min_index = 0;  
	    for (auto const& hit:  ele.second){
		  if(hit->Time() < min_time){ 
	            min_time = hit->Time();
		    min_index = all_index;
		 }
		 all_index++;
	    }
	    
	    bool found_tpc_traj_point=false;
	    for(size_t i=0; i<True_tof_sim_particles[ele.first][min_index]->NumberTrajectoryPoints(); i++){
	        const TLorentzVector& pos = True_tof_sim_particles[ele.first][min_index]->Position(i);
                geo::Point_t const point{pos.X(),pos.Y(),pos.Z()};
		if(tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point)){
                   auto const opDetPos = cryo0.OpDet(cryo0.GetClosestOpDet(point)).GetCenter();
                   double dprop=(opDetPos - point).R();
		   double tprop=pos.T() + dprop*LAR_PROP_DELAY;
		   fTrue_TOF.push_back(tracksps[ele.first].front()->Time()-tprop);
		   fTrue_TOF_hit.push_back(false);
		   fTrue_TOF_pdg.push_back(True_tof_sim_particles[ele.first][min_index]->PdgCode());
		   fTrue_TOF_part_ID.push_back(True_tof_sim_particles[ele.first][min_index]->TrackId());
		   if(True_tof_sim_particles[ele.first][min_index]->Process()=="primary") fTrue_TOF_part_prim.push_back(true);
		   else fTrue_TOF_part_prim.push_back(false);
		   art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(True_tof_sim_particles[ele.first][min_index]->TrackId());
		   fTrue_TOF_part_org.push_back(truth->NeutrinoSet());
		   fTrue_TOF_traj_T.push_back(pos.T());
                   fTrue_TOF_traj_X.push_back(point.X());
                   fTrue_TOF_traj_Y.push_back(point.Y());
                   fTrue_TOF_traj_Z.push_back(point.Z());
		   fTrue_TOF_traj_in_TPC.push_back(true);
		   found_tpc_traj_point=true;
		   break;
		} // trajectory point is in TPC
	     } // loop over trajectory points
	     
	     if(!found_tpc_traj_point){
	          for(size_t i=0; i<True_tof_sim_particles[ele.first][min_index]->NumberTrajectoryPoints(); i++){
		      const TLorentzVector& pos = True_tof_sim_particles[ele.first][min_index]->Position(i);
                      const geo::Point_t point{pos.X(),pos.Y(),pos.Z()};
		      if(cryo0.ContainsPosition(point)){
                         auto const opDetPos = cryo0.OpDet(cryo0.GetClosestOpDet(point)).GetCenter();
                         double dprop=(opDetPos - point).R();
		         double tprop=pos.T() + dprop*LAR_PROP_DELAY;
		         fTrue_TOF.push_back(tracksps[ele.first].front()->Time()-tprop);
		         fTrue_TOF_hit.push_back(false);
		         fTrue_TOF_pdg.push_back(True_tof_sim_particles[ele.first][min_index]->PdgCode());
		         fTrue_TOF_part_ID.push_back(True_tof_sim_particles[ele.first][min_index]->TrackId());
		         if(True_tof_sim_particles[ele.first][min_index]->Process()=="primary") fTrue_TOF_part_prim.push_back(true);
		         else fTrue_TOF_part_prim.push_back(false);
		         art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(True_tof_sim_particles[ele.first][min_index]->TrackId());
		         fTrue_TOF_part_org.push_back(truth->NeutrinoSet());
		         fTrue_TOF_traj_T.push_back(pos.T());
                         fTrue_TOF_traj_X.push_back(point.X());
                         fTrue_TOF_traj_Y.push_back(point.Y());
                         fTrue_TOF_traj_Z.push_back(point.Z());
		         fTrue_TOF_traj_in_TPC.push_back(false);
		         break;
		      } // found a trajectory point inside CRYO
		  } // loop over trajectory points
	      } // did not find a trajectory point inside TPC
	 }
    }
 }
 
 //==================================================================================================
 
 fMatchTree->Fill();
} // End of Analyze function

//==========================================================================================
bool ToFAnalyzer::SpacePointCompare(const art::Ptr<CRTSpacePoint>& sp1, const art::Ptr<CRTSpacePoint>& sp2) {
     if(sp1->Time()     != sp2->Time())     return false;
     if(sp1->Pos()      != sp2->Pos())      return false;
     if(sp1->Err()      != sp2->Err())      return false;
     if(sp1->PE()       != sp2->PE())       return false;
     if(sp1->TimeErr()  != sp2->TimeErr())  return false;
     if(sp1->Complete() != sp2->Complete()) return false;

     return true;
}
//===========================================================================================

//===========================================================================================
double ToFAnalyzer::length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  art::ServiceHandle<geo::Geometry> geom;
  double xmin = -2.0 * geom->DetHalfWidth();
  double xmax = 2.0 * geom->DetHalfWidth();
  double ymin = -geom->DetHalfHeight();
  double ymax = geom->DetHalfHeight();
  double zmin = 0.;
  double zmax = geom->DetLength();
  
  double result = 0.;
  TVector3 disp;
  int n = part.NumberTrajectoryPoints();
  bool first = true;
  
 for(int i = 0; i < n; ++i) {
    double mypos[3] = {part.Vx(i), part.Vy(i), part.Vz(i)};
    if(mypos[0] >= xmin && mypos[0] <= xmax && mypos[1] >= ymin && mypos[1] <= ymax && mypos[2] >= zmin && mypos[2] <= zmax){
       double xGen   = part.Vx(i);
       double newX = xGen;
       
       TVector3 pos(newX,part.Vy(i),part.Vz(i));
       if(first){
          start = pos;
       }
       else {
          disp -= pos;
          result += disp.Mag();
       }
       first=false;
       disp = pos;
       end = pos;
   }
  }
  return result;
}	
//================================================================================================

//================================================================================================
void ToFAnalyzer::ClearVecs()
{
   frun = -9999;
   fsubrun = -9999;
   fevent = -9999;
   
   fnu_pdg.clear();
   fnu_E.clear();
   fnu_mode.clear(); 
   fnu_CCNC.clear();
   fnu_posX.clear();
   fnu_posY.clear(); 
   fnu_posZ.clear(); 
   fnu_T.clear();
   fnu_CRYO.clear();
   fnu_TPC.clear(); 
   
   fg4_pdg.clear();
   fg4_trkid.clear();
   fg4st_TPC.clear();
   fg4en_TPC.clear();
   fg4st_CRYO.clear();
   fg4en_CRYO.clear();
   fg4is_prim.clear();
   fg4_T0.clear();
   fg4_T0_CRYO.clear();
   fg4_T0_TPC.clear();
   fg4_org.clear();
   fg4_stX.clear(); 
   fg4_stY.clear();
   fg4_stZ.clear();
   fg4_enX.clear(); 
   fg4_enY.clear();
   fg4_enZ.clear();
   fg4_tlen.clear();
   fg4_E.clear();
   fg4_is_in_TPC.clear();
   fg4_is_in_CRYO.clear(); 
   
   fLhit_tof_vec.clear();
   fLhit_frmtrk_vec.clear();
   fLhit_frmhit_vec.clear();
   fLhit_crtspkey_vec.clear();
   fLhit_crttrkkey_vec.clear();
   fLhit_crttime_t1_vec.clear();
   //   fLhit_crttime_t0_vec.clear();
   fLhit_crtpos_X_vec.clear();
   fLhit_crtpos_Y_vec.clear();
   fLhit_crtpos_Z_vec.clear();
   fLhit_crtpe_vec.clear();
   fLhit_crttgr_vec.clear();
   fLhit_pmthitkey_vec.clear();
   fLhit_pmthitT_vec.clear();
   fLhit_pmthitX_vec.clear();
   fLhit_pmthitY_vec.clear();
   fLhit_pmthitZ_vec.clear();
   fLhit_pmthitpe_vec.clear();
   
   fChit_tof_vec.clear();
   fChit_frmtrk_vec.clear();
   fChit_frmhit_vec.clear();
   fChit_crtspkey_vec.clear();
   fChit_crttrkkey_vec.clear();
   fChit_crttime_t1_vec.clear();
   //   fChit_crttime_t0_vec.clear();
   fChit_crtpos_X_vec.clear();
   fChit_crtpos_Y_vec.clear();
   fChit_crtpos_Z_vec.clear();
   fChit_crtpe_vec.clear();
   fChit_crttgr_vec.clear();
   fChit_pmthitkey_vec.clear();
   fChit_pmthitT_vec.clear();
   fChit_pmthitX_vec.clear();
   fChit_pmthitY_vec.clear();
   fChit_pmthitZ_vec.clear();
   fChit_pmthitpe_vec.clear();
   
   fLflsh_tof_vec.clear();
   fLflsh_frmtrk_vec.clear();
   fLflsh_frmhit_vec.clear();
   fLflsh_crtspkey_vec.clear();
   fLflsh_crttrkkey_vec.clear();
   fLflsh_crttime_t1_vec.clear();
   //   fLflsh_crttime_t0_vec.clear();
   fLflsh_crtpos_X_vec.clear();
   fLflsh_crtpos_Y_vec.clear();
   fLflsh_crtpos_Z_vec.clear();
   fLflsh_crtpe_vec.clear();
   fLflsh_crttgr_vec.clear();
   fLflsh_pmtflshtpcID_vec.clear();
   fLflsh_pmtflshkey_vec.clear();
   fLflsh_pmtflshT_vec.clear();
   fLflsh_pmtflshY_vec.clear();
   fLflsh_pmtflshZ_vec.clear();
   fLflsh_pmtflshpe_vec.clear();
   
   fCflsh_tof_vec.clear();
   fCflsh_frmtrk_vec.clear();
   fCflsh_frmhit_vec.clear();
   fCflsh_crtspkey_vec.clear();
   fCflsh_crttrkkey_vec.clear();
   fCflsh_crttime_t1_vec.clear();
   //   fCflsh_crttime_t0_vec.clear();
   fCflsh_crtpos_X_vec.clear();
   fCflsh_crtpos_Y_vec.clear();
   fCflsh_crtpos_Z_vec.clear();
   fCflsh_crtpe_vec.clear();
   fCflsh_crttgr_vec.clear();
   fCflsh_pmtflshtpcID_vec.clear();
   fCflsh_pmtflshkey_vec.clear();
   fCflsh_pmtflshT_vec.clear();
   fCflsh_pmtflshY_vec.clear();
   fCflsh_pmtflshZ_vec.clear();
   fCflsh_pmtflshpe_vec.clear();
   
   fLflshhit_tof_vec.clear();
   fLflshhit_frmtrk_vec.clear();
   fLflshhit_frmhit_vec.clear();
   fLflshhit_crtspkey_vec.clear();
   fLflshhit_crttrkkey_vec.clear();
   fLflshhit_crttime_t1_vec.clear();
   //   fLflshhit_crttime_t0_vec.clear();
   fLflshhit_crtpos_X_vec.clear();
   fLflshhit_crtpos_Y_vec.clear();
   fLflshhit_crtpos_Z_vec.clear();
   fLflshhit_crtpe_vec.clear();
   fLflshhit_crttgr_vec.clear();
   fLflshhit_pmtflshtpcID_vec.clear();
   fLflshhit_pmtflshkey_vec.clear();
   fLflshhit_pmtkey_vec.clear();
   fLflshhit_pmtflshT_vec.clear();
   fLflshhit_pmtflshX_vec.clear();
   fLflshhit_pmtflshY_vec.clear();
   fLflshhit_pmtflshZ_vec.clear();
   fLflshhit_pmtflshpe_vec.clear();
   
   fCflshhit_tof_vec.clear();
   fCflshhit_frmtrk_vec.clear();
   fCflshhit_frmhit_vec.clear();
   fCflshhit_crtspkey_vec.clear();
   fCflshhit_crttrkkey_vec.clear();
   fCflshhit_crttime_t1_vec.clear();
   //   fCflshhit_crttime_t0_vec.clear();
   fCflshhit_crtpos_X_vec.clear();
   fCflshhit_crtpos_Y_vec.clear();
   fCflshhit_crtpos_Z_vec.clear();
   fCflshhit_crtpe_vec.clear();
   fCflshhit_crttgr_vec.clear();
   fCflshhit_pmtflshtpcID_vec.clear();
   fCflshhit_pmtflshkey_vec.clear();
   fCflshhit_pmtkey_vec.clear();
   fCflshhit_pmtflshT_vec.clear();
   fCflshhit_pmtflshX_vec.clear();
   fCflshhit_pmtflshY_vec.clear();
   fCflshhit_pmtflshZ_vec.clear();
   fCflshhit_pmtflshpe_vec.clear();
   
   fTrue_TOF.clear();
   fTrue_TOF_hit.clear();
   fTrue_TOF_pdg.clear();
   fTrue_TOF_part_ID.clear();
   fTrue_TOF_part_prim.clear();
   fTrue_TOF_part_org.clear();
   fTrue_TOF_traj_T.clear();
   fTrue_TOF_traj_X.clear();
   fTrue_TOF_traj_Y.clear();
   fTrue_TOF_traj_Z.clear();
   fTrue_TOF_traj_in_TPC.clear();
}
//================================================================================================
DEFINE_ART_MODULE(ToFAnalyzer)
}
