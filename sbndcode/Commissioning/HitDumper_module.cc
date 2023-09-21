////////////////////////////////////////////////////////////////////////
// Class:       Hitdumper
// Module Type: analyzer
// File:        Hitdumper_module.cc
//
////////////////////////////////////////////////////////////////////////
#ifndef Hitdumper_Module
#define Hitdumper_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// SBN/SBND includes
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Commissioning/MuonTrack.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
#include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"
#include "sbnobj/SBND/Trigger/CRTmetric.hh"

// Truth includes
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

const int MAX_INT = std::numeric_limits<int>::max();
const long int TIME_CORRECTION = (long int) std::numeric_limits<int>::max() * 2;
const int DEFAULT_VALUE = -9999;

enum CRTOrientation {
  kCRTNotDefined = -1,   ///< Not defined
  kCRTHorizontal = 0,    ///< Horizontal Tagger
  kCRTVertical = 1,      ///< VertivalTagger
};

enum PhotoDetectorType {
  kPDNotDefined = -1,   ///< Not defined
  kPMTCoated = 0,       ///< Coated PMT
  kPMTUnCoated = 1,     ///< Uncoated PMT
  kXArapucaVis,         ///< Arapuca Vis
  kXArapucaVuv,         ///< Arapuca VUV
};

class Hitdumper : public art::EDAnalyzer {
public:
  explicit Hitdumper(fhicl::ParameterSet const & p);
  virtual ~Hitdumper();

  // This method is called once, at the start of the job. In this
  // example, it will define the histograms and n-tuples we'll write.
  void beginJob() override;

  // This method is called once, at the start of each run. It's a
  // good place to read databases or files that may have
  // run-dependent information.
  // void beginRun(const art::Run& run);

  // This method reads in any parameters from the .fcl files. This
  // method is called 'reconfigure' because it might be called in the
  // middle of a job; e.g., if the user changes parameter values in an
  // interactive event display.
  void reconfigure(fhicl::ParameterSet const& pset);

  // The analysis routine, called once per event.
  void analyze (const art::Event& evt) override;

  // Called at the beginning of every subrun
  virtual void beginSubRun(art::SubRun const& sr) override;
private:

  /// Resets the variables that are saved to the TTree
  void ResetVars();
  /// Resets wire hits tree variables
  void ResetWireHitsVars(int n);
  /// Resets crt strips tree variables
  void ResetCRTStripsVars(int n);
  /// Resets custom crt tracks tree variables
  void ResetCRTCustomTracksVars(int n);
  /// Resets crt tracks tree variables
  void ResetCRTTracksVars(int n);
  /// Resets crt hits tree variables
  void ResetCRTHitsVars(int n);
  /// Resets optical hits tree variables
  void ResetOpHitsVars(int n);
  /// Resets pmt hardware trigger variables
  void ResetPmtTriggerVars(int n);
  /// Rests pmt software trigger variables
  void ResetPmtSoftTriggerVars();
  /// Resets crt software trigger variables
  void ResetCrtSoftTriggerVars();
  /// Resets crossing muon tracks tree variables
  void ResetMuonTracksVars(int n);
  /// Resets crossing muon hit tree variables 
  void ResetMuonHitVars(int n); 
  /// Resize the data structure for MCNeutrino particles
  void ResizeMCNeutrino(int nNeutrinos);
  /// Resize the data structure for Genie primaries
  void ResizeGenie(int nPrimaries);
  /// Resize the data structure for MCParticles
  void ResizeMCParticle(int nParticles);
  /// Resize the data structure for MCTracks
  void ResizeMCTrack(int nTracks);
  /// Resize the data structure for MCShowers
  void ResizeMCShower(int nShowers);

  opdet::sbndPDMapAlg _pd_map;

  TTree* fTree;
  //run information
  int _run;        ///< The run number
  int _subrun;     ///< The subrun number
  int _event;      ///< The event number
  double _evttime; ///< The event time
  int _t0;         ///< The t0

  // Wire hits variables
  int                 _nhits;                   ///< Number of reco hits in the event
  std::vector<int>    _hit_cryostat;            ///< Cryostat where the hit belongs to
  std::vector<int>    _hit_tpc;                 ///< TPC where the hit belongs to
  std::vector<int>    _hit_plane;               ///< Plane where the hit belongs to
  std::vector<int>    _hit_wire;                ///< Wire where the hit belongs to
  std::vector<int>    _hit_channel;             ///< Channel where the hit belongs to
  std::vector<double> _hit_peakT;               ///< Hit peak time
  std::vector<double> _hit_charge;              ///< Hit charge
  std::vector<double> _hit_ph;                  ///< Hit ph?
  std::vector<double> _hit_width;               ///< Hit width
  std::vector<double> _hit_full_integral;       ///< Hit charge integral
  std::vector<int>    _waveform_number;         ///< Number for each waveform, to allow for searching
  std::vector<double> _adc_on_wire;             ///< ADC on wire to draw waveform
  std::vector<int>    _time_for_waveform;       ///<Time for waveform to plot
  int                 _adc_count;               ///<Used for plotting waveforms
  std::vector<int>    _waveform_integral;       ///<Used to see progression of the waveform integral
  std::vector<int>    _adc_count_in_waveform;   ///<Used to view all waveforms on a hitplane together


  // CRT strip variables
  int _nstrips;                          ///< Number of CRT strips
  std::vector<int> _crt_plane;           ///< CRT plane
  std::vector<int> _crt_module;          ///< CRT module
  std::vector<int> _crt_strip;           ///< CRT strip
  std::vector<int> _crt_orient;          ///< CRT orientation (0 for y (horizontal) and 1 for x (vertical))
  std::vector<double> _crt_time;         ///< CRT time
  std::vector<double> _crt_adc;          ///< CRT adc
  std::vector<double> _crt_pos_x;          ///< CRT position X
  std::vector<double> _crt_pos_y;          ///< CRT position Y
  std::vector<double> _crt_pos_z;          ///< CRT position Z

  // CRT track variables
  int _nctrks;                          ///< Number of created CRT tracks
  std::vector<double> _ctrk_x1;         ///< CRT track x1
  std::vector<double> _ctrk_y1;         ///< CRT track y1
  std::vector<double> _ctrk_z1;         ///< CRT track z1
  std::vector<double> _ctrk_t1;         ///< CRT track t1
  std::vector<double> _ctrk_adc1;       ///< CRT track adc1
  std::vector<int> _ctrk_mod1x;         ///< CRT track mod2x
  std::vector<double> _ctrk_x2;         ///< CRT track x2
  std::vector<double> _ctrk_y2;         ///< CRT track y2
  std::vector<double> _ctrk_z2;         ///< CRT track z2
  std::vector<double> _ctrk_t2;         ///< CRT track t2
  std::vector<double> _ctrk_adc2;       ///< CRT track adc2
  std::vector<int> _ctrk_mod2x;         ///< CRT track mod2x

  // CRT hits variables
  int _nchits;                           ///< Number of CRT hits
  std::vector<double> _chit_x;           ///< CRT hit x
  std::vector<double> _chit_y;           ///< CRT hit y
  std::vector<double> _chit_z;           ///< CRT hit z
  std::vector<double> _chit_time;        ///< CRT hit time
  // std::vector<double> _chit_adc;         ///< CRT hit adc
  std::vector<int> _chit_plane;          ///< CRT hit plane

  // CRT track variables
  int _ncts;                            ///< Number of CRT tracks
  std::vector<double> _ct_time;         ///< CRT track time
  std::vector<double> _ct_pes;          ///< CRT track PEs
  std::vector<double> _ct_x1;           ///< CRT track x1
  std::vector<double> _ct_y1;           ///< CRT track y1
  std::vector<double> _ct_z1;           ///< CRT track z1
  std::vector<double> _ct_x2;           ///< CRT track x2
  std::vector<double> _ct_y2;           ///< CRT track y2
  std::vector<double> _ct_z2;           ///< CRT track z2

  int _nophits;                               ///< Number of Optical Hits
  std::vector<int> _ophit_opch;               ///< OpChannel of the optical hit
  std::vector<int> _ophit_opdet;              ///< OpDet of the optical hit
  std::vector<double> _ophit_peakT;           ///< Peak time of the optical hit [us]
  std::vector<double> _ophit_startT;          ///< Start time of the optical hit [us]
  std::vector<double> _ophit_riseT;           ///< Rise time of the optical hit [ns]
  std::vector<double> _ophit_width;           ///< Width of the optical hit
  std::vector<double> _ophit_area;            ///< Area of the optical hit
  std::vector<double> _ophit_amplitude;       ///< Amplitude of the optical hit
  std::vector<double> _ophit_pe;              ///< PEs of the optical hit
  std::vector<double> _ophit_opdet_x;         ///< OpDet X coordinate of the optical hit
  std::vector<double> _ophit_opdet_y;         ///< OpDet Y coordinate of the optical hit
  std::vector<double> _ophit_opdet_z;         ///< OpDet Z coordinate of the optical hit
  std::vector<int> _ophit_opdet_type;         ///< OpDet tyoe of the optical hit

  //pmt hardware trigger variables
  std::vector<int> _pmtTrigger_npmtshigh;    ///< number of pmt pairs above threshold, index = time during trigger window (usually beam spill)
  int _pmtTrigger_maxpassed;    ///< maximum number of pmt pairs above threshold during trigger window (usually beam spill)

  // PMT software trigger variables 
  bool   _pmtSoftTrigger_foundBeamTrigger;   /// Whether the beam spill was found or not 
  int    _pmtSoftTrigger_tts;                /// Trigger Time Stamp (TTS), ns (relative to start of beam spill)
  double _pmtSoftTrigger_promptPE;           /// Total photoelectron count 100 ns after the TTS
  double _pmtSoftTrigger_prelimPE;           /// Total photoelectron count before the TTS, during the beam spill             
  int    _pmtSoftTrigger_nAboveThreshold;    /// number of individual PMTs above ADC threshold (fcl) during the beam spill
  // std::vector<sbnd::trigger::pmtInfo> _pmtSoftTrigger_pmtInfoVec; /// vector of PMT information 

  // CRT software trigger variables
  int    _crtSoftTrigger_hitsperplane[7];       ///< Number of (very low level) CRT hits per plane

  // Muon track variables 
  int _nmuontrks;                            ///< number of muon tracks
  std::vector<double> _muontrk_t0;           ///< t0 (time of interaction)
  std::vector<float>  _muontrk_x1;           ///< x coordinate closer to anode
  std::vector<float>  _muontrk_y1;           ///< y coordinate closer to anode 
  std::vector<float>  _muontrk_z1;           ///< z coordinate closer to anode 
  std::vector<float>  _muontrk_x2;           ///< x coordinate closer to cathode 
  std::vector<float>  _muontrk_y2;           ///< y coordinate closer to cathode 
  std::vector<float>  _muontrk_z2;           ///< z coordinate closer to cathode
  std::vector<float>  _muontrk_theta_xz;     ///< theta_xz trajectory angle 
  std::vector<float>  _muontrk_theta_yz;     ///< theta_yz trajectory angle 
  std::vector<int>    _muontrk_tpc;          ///< tpc that muon is located in 
  std::vector<int>    _muontrk_type;         ///< type of muon track

  // Muon Hit variables
  int                 _nmhits;               ///< Number of muon collection hits per track
  std::vector<int>    _mhit_trk;             ///< Track number that the hit belongs to
  std::vector<int>    _mhit_tpc;             ///< TPC where the hit belongs to
  std::vector<int>    _mhit_wire;            ///< Wire where the hit belongs to
  std::vector<int>    _mhit_channel;         ///< Channel where the hit belongs to
  std::vector<double> _mhit_peakT;           ///< Hit peak time
  std::vector<double> _mhit_charge;          ///< Hit charge
  
  //mctruth information
  size_t MaxMCNeutrinos;     ///! The number of MCNeutrinos there is currently room for
  Int_t     mcevts_truth;                     ///< number of neutrino Int_teractions in the spill
  std::vector<Int_t>     nuScatterCode_truth; ///< Scattering code given by Genie for each neutrino
  std::vector<Int_t>     nuID_truth;          ///< Unique ID of each true neutrino
  std::vector<Int_t>     nuPDG_truth;         ///< neutrino PDG code
  std::vector<Int_t>     ccnc_truth;          ///< 0=CC 1=NC
  std::vector<Int_t>     mode_truth;          ///< 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  std::vector<Float_t>   enu_truth;           ///< true neutrino energy
  std::vector<Float_t>   Q2_truth;            ///< Momentum transfer squared
  std::vector<Float_t>   W_truth;             ///< hadronic invariant mass
  std::vector<Int_t>     hitnuc_truth;        ///< hit nucleon
  std::vector<Float_t>   nuvtxx_truth;        ///< neutrino vertex x
  std::vector<Float_t>   nuvtxy_truth;        ///< neutrino vertex y
  std::vector<Float_t>   nuvtxz_truth;        ///< neutrino vertex z
  std::vector<Float_t>   nu_dcosx_truth;      ///< neutrino dcos x
  std::vector<Float_t>   nu_dcosy_truth;      ///< neutrino dcos y
  std::vector<Float_t>   nu_dcosz_truth;      ///< neutrino dcos z
  std::vector<Float_t>   lep_mom_truth;       ///< lepton momentum
  std::vector<Float_t>   lep_dcosx_truth;     ///< lepton dcos x
  std::vector<Float_t>   lep_dcosy_truth;     ///< lepton dcos y
  std::vector<Float_t>   lep_dcosz_truth;     ///< lepton dcos z

  //flux information
  std::vector<Float_t>  tpx_flux;             ///< Px of parent particle leaving BNB target
  std::vector<Float_t>  tpy_flux;             ///< Py of parent particle leaving BNB target
  std::vector<Float_t>  tpz_flux;             ///< Pz of parent particle leaving BNB target
  std::vector<Int_t>    tptype_flux;         ///< Type of parent particle leaving BNB target

  //genie information
  size_t MaxGeniePrimaries = 0;
  Int_t     genie_no_primaries;
  std::vector<Int_t>     genie_primaries_pdg;
  std::vector<Float_t>   genie_Eng;
  std::vector<Float_t>   genie_Px;
  std::vector<Float_t>   genie_Py;
  std::vector<Float_t>   genie_Pz;
  std::vector<Float_t>   genie_P;
  std::vector<Int_t>     genie_status_code;
  std::vector<Float_t>   genie_mass;
  std::vector<Int_t>     genie_trackID;
  std::vector<Int_t>     genie_ND;
  std::vector<Int_t>     genie_mother;

  //MCParticle Info
  size_t MaxMCParticles = 0;
  Int_t     mcpart_no_primaries;                 
  std::vector<Int_t>    mcpart_pdg;              
  std::vector<Int_t>    mcpart_status;           
  std::vector<std::string>    mcpart_process;
  std::vector<std::string>    mcpart_endprocess;
  std::vector<Float_t>  mcpart_Eng;              
  std::vector<Float_t>  mcpart_EndE;
  std::vector<Float_t>  mcpart_Mass;
  std::vector<Float_t>  mcpart_Px;
  std::vector<Float_t>  mcpart_Py;
  std::vector<Float_t>  mcpart_Pz;
  std::vector<Float_t>  mcpart_P;
  std::vector<Float_t>  mcpart_StartPointx;
  std::vector<Float_t>  mcpart_StartPointy;
  std::vector<Float_t>  mcpart_StartPointz;
  std::vector<Float_t>  mcpart_StartT;  
  std::vector<Float_t>  mcpart_EndT;          
  std::vector<Float_t>  mcpart_EndPointx;
  std::vector<Float_t>  mcpart_EndPointy;
  std::vector<Float_t>  mcpart_EndPointz;
  std::vector<Float_t>  mcpart_theta_xz;    
  std::vector<Float_t>  mcpart_theta_yz;    
  std::vector<Int_t>    mcpart_NumberDaughters;
  std::vector<Int_t>    mcpart_TrackId;
  std::vector<Int_t>    mcpart_Mother;

  //MCTrack info
  size_t MaxMCTracks = 0;
  Int_t mctrack_no_primaries;
  std::vector<Int_t>    mctrack_pdg;                      
  std::vector<Int_t>    mctrack_TrackId;

  //MCShower info
  size_t MaxMCShowers = 0;
  Int_t mcshower_no_primaries;
  std::vector<Int_t>    mcshower_pdg;                       
  std::vector<Int_t>    mcshower_TrackId;


  TTree* _sr_tree; ///< A tree filled per subrun (for POT accounting)
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< POTs in each subrun


  int _max_hits;                    ///< maximum number of hits (to be set via fcl)
  int _max_ophits;                  ///< maximum number of hits (to be set via fcl)
  int _max_samples;                 ///< maximum number of samples (to be set via fcl)
  int _max_chits;                   ///< maximum number of CRT hits (to be set via fcl)
  int _max_nctrks;                  ///< maximum number of CRT tracks (to be set via fcl)

  std::string fHitsModuleLabel;     ///< Label for Hit dataproduct (to be set via fcl)
  std::string fLArG4ModuleLabel;    ///< Label for LArG4 dataproduct (to be set via fcl)
  std::string fCRTStripModuleLabel; ///< Label for CRTStrip dataproduct (to be set via fcl)
  std::string fCRTHitModuleLabel;   ///< Label for CRTHit dataproduct (to be set via fcl)
  std::string fCRTTrackModuleLabel; ///< Label for CRTTrack dataproduct (to be set via fcl)
  std::string fpmtTriggerModuleLabel; ///< Label for pmtTrigger dataproduct (to be set vis fcl)
  std::string fpmtSoftTriggerModuleLabel; ///< Label for pmt software trigger data product (to be set via fcl)
  std::string fcrtSoftTriggerModuleLabel; ///< Label for crt software trigger data product (to be set via fcl)
  std::string fMuonTrackModuleLabel;  ///< Label for MuonTrack dataproduct (to be set via fcl)
  std::string fDigitModuleLabel;    ///< Label for digitizer (to be set via fcl)
  std::string fGenieGenModuleLabel; ///< Label for Genie dataproduct (to be set via fcl)
  std::string fMCParticleModuleLabel; ///< Label for MCParticle dataproduct (to be set via fcl)
  std::string fMCTrackModuleLabel; ///< Label for MCTrack dataproduct (to be set via fcl)
  std::string fMCShowerModuleLabel; ///< Label for MCShower dataproduct (to be set via fcl)
  std::vector<std::string> fOpHitsModuleLabels; ///< Labels for OpHit dataproducts (to be set via fcl)

  // double fSelectedPDG;

  bool fkeepCRThits;       ///< Keep the CRT hits (to be set via fcl)
  bool fkeepCRTstrips;     ///< Keep the CRT strips (to be set via fcl)
  bool fmakeCRTtracks;     ///< Make the CRT tracks (to be set via fcl)
  bool freadCRTtracks;     ///< Keep the CRT tracks (to be set via fcl)
  bool freadOpHits;        ///< Add OpHits to output (to be set via fcl)
  bool freadMuonTracks;    ///< Add MuonTracks to output (to be set via fcl)
  bool freadMuonHits;      ///< Add MuonTrack hits to output(to be set via fcl)
  bool freadTruth;         ///< Add Truth info to output (to be set via fcl)
  bool freadMCParticle;    ///< Add MCParticle info to output (to be set via fcl)
  bool freadpmtTrigger;    ///< Add pmt hardware trigger info to output (to be set via fcl)
  bool freadpmtSoftTrigger;///< Add pmt software trigger info to output (to be set via fcl)
  bool freadcrtSoftTrigger;///< Add crt software trigger info to output (to be set via fcl)
  bool fsavePOTInfo;       ///< Add POT info to output (to be set via fcl)
  bool fcheckTransparency; ///< Checks for wire transprency (to be set via fcl)
  bool fUncompressWithPed; ///< Uncompresses the waveforms if true (to be set via fcl)
  int fWindow;
  bool fSkipInd;           ///< If true, induction planes are not saved (to be set via fcl)
  // double fSelectedPDG;

  std::vector<int> fKeepTaggerTypes = {0, 1, 2, 3, 4, 5, 6}; ///< Taggers to keep (to be set via fcl)

  sbnd::CRTHitRecoAlg hitAlg;

  geo::GeometryCore const* fGeometryService;
  // detinfo::ElecClock fTrigClock;
  art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
  const geo::AuxDetGeometry* fAuxDetGeo;
  const geo::AuxDetGeometryCore* fAuxDetGeoCore;
};


Hitdumper::Hitdumper(fhicl::ParameterSet const& pset)
: EDAnalyzer(pset)
{

  fGeometryService = lar::providerFrom<geo::Geometry>();
  // fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  // fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // fTrigClock = fDetectorClocks->TriggerClock();
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

  // Read in the parameters from the .fcl file.
  this->reconfigure(pset);
}

void Hitdumper::reconfigure(fhicl::ParameterSet const& p)
{

  _max_hits = p.get<int>("MaxHits", 50000);
  _max_ophits = p.get<int>("MaxOpHits", 50000);
  _max_samples = p.get<int>("MaxSamples", 5001);
  _max_chits = p.get<int>("MaxCRTHits", 5000);
  _max_nctrks = p.get<int>("MaxCRTTracks", 10);

  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
  fDigitModuleLabel    = p.get<std::string>("DigitModuleLabel", "daq");
  fLArG4ModuleLabel    = p.get<std::string>("LArG4ModuleLabel", "largeant");
  fCRTStripModuleLabel = p.get<std::string>("CRTStripModuleLabel", "crt");
  fCRTHitModuleLabel   = p.get<std::string>("CRTHitModuleLabel", "crthit");
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel", "crttrack");
  fOpHitsModuleLabels  = p.get<std::vector<std::string>>("OpHitsModuleLabel");
  fpmtTriggerModuleLabel = p.get<std::string>("pmtTriggerModuleLabel", "pmttriggerproducer");
  fpmtSoftTriggerModuleLabel = p.get<std::string>("pmtSoftTriggerModuleLabel", "pmtSoftwareTrigger");
  fcrtSoftTriggerModuleLabel = p.get<std::string>("crtSoftTriggerModuleLabel", "MetricProducer");
  fMuonTrackModuleLabel  = p.get<std::string>("MuonTrackModuleLabel", "MuonTrackProducer");
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel", "generator");
  fMCParticleModuleLabel    = p.get<std::string>("MCParticleModuleLabel ", "largeant");
  fMCTrackModuleLabel    = p.get<std::string>("MCTrackModuleLabel ", "mcreco");
  fMCShowerModuleLabel    = p.get<std::string>("MCShowerModuleLabel ", "mcreco");

  fkeepCRThits       = p.get<bool>("keepCRThits",true);
  fkeepCRTstrips     = p.get<bool>("keepCRTstrips",false);
  fmakeCRTtracks     = p.get<bool>("makeCRTtracks",true);
  freadCRTtracks     = p.get<bool>("readCRTtracks",true);
  freadOpHits        = p.get<bool>("readOpHits",true);
  freadpmtTrigger    = p.get<bool>("readpmtTrigger",true);
  freadpmtSoftTrigger= p.get<bool>("readpmtSoftTrigger",true);
  freadcrtSoftTrigger= p.get<bool>("readcrtSoftTrigger",false);
  freadMuonTracks    = p.get<bool>("readMuonTracks",true);
  freadMuonHits      = p.get<bool>("readMuonHits",false);
  fcheckTransparency = p.get<bool>("checkTransparency",false);
  freadTruth         = p.get<bool>("readTruth",true);
  freadMCParticle    = p.get<bool>("readMCParticle",false);
  fsavePOTInfo       = p.get<bool>("savePOTinfo",true);
  fUncompressWithPed = p.get<bool>("UncompressWithPed",false);

  fWindow            = p.get<int>("window",100);
  fKeepTaggerTypes   = p.get<std::vector<int>>("KeepTaggerTypes");

  fSkipInd           = p.get<bool>("SkipInduction",false);
}


Hitdumper::~Hitdumper()
{
  // Clean up dynamic memory and other resources here.
}


void Hitdumper::analyze(const art::Event& evt)
{
  // Reset the TTree variables
  ResetVars();

  _run = evt.run();
  _subrun = evt.subRun();
  _event = evt.id().event();

  _t0 = 0.;
  // t0 = detprop->TriggerOffset();  // units of TPC ticks

  //
  // Hits
  //
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)) {
    art::fill_ptr_vector(hitlist, hitListHandle);
    _nhits = hitlist.size();

    // Calculate how many hits we will save if skipping the induction planes
    if (fSkipInd) {
      _nhits = 0;
      for (auto h : hitlist) {
        if (h->WireID().Plane == 2) {
          _nhits++;
        }
      }
    }
  }
  else {
    std::cout << "Failed to get recob::Hit data product." << std::endl;
    _nhits = 0;
  }

  if (_nhits > _max_hits) {
    std::cout << "Available hits are " << _nhits
              << ", which is above the maximum number allowed to store." << std::endl;
    std::cout << "Will only store " << _max_hits << "hits." << std::endl;
    _nhits = _max_hits;
  }

  ResetWireHitsVars(_nhits);

  size_t counter = 0;
  for (size_t i = 0; i < hitlist.size(); ++i) {
    geo::WireID wireid = hitlist[i]->WireID();
    if (fSkipInd && wireid.Plane != 2) {
      continue;
    }

    _hit_cryostat[counter] = wireid.Cryostat;
    _hit_tpc[counter] = wireid.TPC;
    _hit_plane[counter] = wireid.Plane;
    _hit_wire[counter] = wireid.Wire;
    _hit_channel[counter] = hitlist[i]->Channel();
    // peak time needs plane dependent offset correction applied.
    _hit_peakT[counter] = hitlist[i]->PeakTime();
    _hit_charge[counter] = hitlist[i]->Integral();
    _hit_ph[counter] = hitlist[i]->PeakAmplitude();
    _hit_width[counter] = hitlist[i]->RMS();
    counter ++;
  }

  //
  // CRT strips
  //
  int _nstr = 0;
  art::Handle<std::vector<sbnd::crt::CRTData> > crtStripListHandle;
  std::vector<art::Ptr<sbnd::crt::CRTData> > striplist;
  // art::Handle< std::vector<sbnd::crt::CRTData> > crtStripListHandle;
  // std::vector< art::Ptr<sbnd::crt::CRTData> > striplist;
  if (evt.getByLabel(fCRTStripModuleLabel, crtStripListHandle))  {
    art::fill_ptr_vector(striplist, crtStripListHandle);
    _nstr = striplist.size();
  } else {
    std::cout << "Failed to get sbnd::crt::CRTData data product." << std::endl;
  }

  int ns = 0;
  if (_nstr > _max_chits) _nstr = _max_chits;
  // strips are always in pairs, one entry for each sipm (2 sipms per strip)

  ResetCRTStripsVars(_nstr);

  for (int i = 0; i < _nstr; i += 2){
    uint32_t chan = striplist[i]->Channel();

    //    std::pair<std::string,unsigned> tagger = CRTHitRecoAlg::ChannelToTagger(chan);
    std::pair<std::string,unsigned> tagger = hitAlg.ChannelToTagger(chan);
    sbnd::CRTPlane ip = sbnd::CRTCommonUtils::GetPlaneIndex(tagger.first);

    bool keep_tagger = false;
    for (auto t : fKeepTaggerTypes) {
      if (ip == t) {
        keep_tagger = true;
      }
    }
    // std::cout << "Tagger name " << tagger.first << ", ip " << ip << ", kept? " << (keep_tagger ? "yes" : "no") << std::endl;

    if (ip != sbnd::kCRTNotDefined && keep_tagger) {

      uint32_t ttime = striplist[i]->T0();
      float ctime = (int)ttime * 0.001; // convert form ns to us
      // recover simulation bug where neg times are sotred as unsigned integers
      // if (ttime > MAX_INT) ctime = 0.001 * (ttime - MAX_INT * 2);
      if (ctime < 1600. && ctime > -1400.) {
        uint32_t adc1 = striplist[i]->ADC();
        uint32_t adc2 = striplist[i+1]->ADC();
        if (adc1 > 4095) adc1 = 4095;
        if (adc2 > 4095) adc2 = 4095;
        //    std::cout << tagger.first << " " << tagger.second << std::endl;
        //    int sipm = chan & 1;  // 0 or 1
        int strip = (chan >> 1) & 15;
        int module = (chan>> 5);
        //
        std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
        auto const center = fAuxDetGeoCore->AuxDetChannelToPosition(name, 2*strip);
        _crt_plane.push_back(ip);
        _crt_module.push_back(module);
        _crt_strip.push_back(strip);
        _crt_orient.push_back(tagger.second);
        _crt_time.push_back(ctime);
        _crt_adc.push_back(adc1 + adc2 - 127.2); // -127.2/131.9 correct for gain and 2*ped to get pe
        _crt_pos_x.push_back(center.X());
        _crt_pos_y.push_back(center.Y());
        _crt_pos_z.push_back(center.Z());
        ns++;
      }
    }
  }
  _nstrips = ns;



  //
  // CRT Custom Tracks
  //
  ResetCRTCustomTracksVars(_nstrips);
  _nctrks = 0;
  if (fmakeCRTtracks) {
    int ntr = 0;
    int iflag[1000] = {0};
    for (int i = 0; i < (ns - 1); ++i) {
      if (iflag[i] == 0) {
        iflag[i] = 1;
        float plane1x = 0, plane2x = 0;
        float plane1y = 0, plane2y = 0;
        float plane1tx = 0, plane2tx = 0;
        float plane1ty = 0, plane2ty = 0;
        float plane1xm = -1, plane2xm = -1;
        float plane1ym = -1, plane2ym = -1;
        int  nh1x = 0, nh2x = 0;
        int  nh1y = 0, nh2y = 0;
        float adc1x = 0, adc2x = 0;
        float adc1y = 0, adc2y = 0;
        if (_crt_plane[i] == sbnd::kCRTFaceSouth) { // 1
          if (_crt_orient[i] == kCRTVertical && _crt_adc[i] > 500) { // < 500 hardcoded
            if (nh1x == 0 || (_crt_module[i] == plane1xm)) {
              nh1x++;
              if (_crt_adc[i] > adc1x) {
                plane1tx = _crt_time[i];
                adc1x += _crt_adc[i];
                plane1x = _crt_pos_x[i];
                plane1xm = _crt_module[i];
              }
            }
          }
          else if (_crt_orient[i]==kCRTHorizontal && _crt_adc[i]>500) { // < 500 hardcoded
            if (nh1y==0 || (_crt_module[i]==plane1ym)) {
              nh1y++;
              if (_crt_adc[i]>adc1y) {
                plane1ty=_crt_time[i];
                adc1y+=_crt_adc[i];
                plane1y=_crt_pos_y[i];
                plane1ym=_crt_module[i];
              }
            }
          }
        }
        else {
          if (_crt_orient[i]==kCRTVertical && _crt_adc[i]>500) { // < 500 hardcoded
            if (nh2x==0 ||  (_crt_module[i]==plane2xm)) {
              nh2x++;
              if (_crt_adc[i]>adc2x) {
                plane2tx=_crt_time[i];
                adc2x+=_crt_adc[i];
                plane2x=_crt_pos_x[i];
                plane2xm=_crt_module[i];
              }
            }
          }
          else if (_crt_orient[i]==kCRTHorizontal && _crt_adc[i]>500) { // < 500 hardcoded
            if (nh2y==0 ||  (_crt_module[i]==plane2ym)) {
              nh2y++;
              if (_crt_adc[i]>adc2y) {
                plane2ty=_crt_time[i];
                adc2y+=_crt_adc[i];
                plane2y=_crt_pos_y[i];
                plane2ym=_crt_module[i];
              }
            }
          }
        }
        for (int j=i+1;j<ns;++j) {
          float tdiff = fabs(_crt_time[i]-_crt_time[j]);
          if (tdiff<0.1) {
            iflag[j]=1;
            if (_crt_plane[j]==sbnd::kCRTFaceSouth) {
              if (_crt_orient[j]==kCRTVertical && _crt_adc[j]>1000) {
                if (nh1x==0 ||  (_crt_module[j]==plane1xm)) {
                  nh1x++;
                  if (_crt_adc[j]>adc1x) {
                    plane1tx=_crt_time[j];
                    adc1x+=_crt_adc[j];
                    plane1x=_crt_pos_x[j];
                    plane1xm=_crt_module[j];
                  }
                }
              }
              else if (_crt_orient[j]==kCRTHorizontal && _crt_adc[j]>1000) {
                if (nh1y==0 ||  (_crt_module[j]==plane1ym)) {
                  nh1y++;
                  if (_crt_adc[j]>adc1y) {
                    plane1ty=_crt_time[j];
                    adc1y+=_crt_adc[j];
                    plane1y=_crt_pos_y[j];
                    plane1ym=_crt_module[j];
                  }
                }
              }
            }
            else {
              if (_crt_orient[j]==kCRTVertical && _crt_adc[j]>1000) {
                if (nh2x==0 ||  (_crt_module[j]==plane2xm)) {
                  nh2x++;
                  if (_crt_adc[j]>adc2x) {
                    plane2tx=_crt_time[j];
                    adc2x+=_crt_adc[j];
                    plane2x=_crt_pos_x[j];
                    plane2xm=_crt_module[j];
                  }
                }
              }
              else if (_crt_orient[j]==kCRTHorizontal && _crt_adc[j]>1000) {
                if (nh2y==0 ||  (_crt_module[j]==plane2ym)) {
                  nh2y++;
                  if (_crt_adc[j]>adc2y) {
                    plane2ty=_crt_time[j];
                    adc2y+=_crt_adc[j];
                    plane2y=_crt_pos_y[j];
                    plane2ym=_crt_module[j];
                  }
                }
              }
            }
          }
	      } // look for hits at the same time as hit i
	      if (nh1x>0 && nh1y>0 && nh2x>0 && nh2y>0 && adc1x<9000 && adc1y<9000 && adc2x<9000 && adc2y<9000) {
	      // make a track!
          _ctrk_x1.push_back(plane1x);
          _ctrk_y1.push_back(plane1y);
          _ctrk_z1.push_back(-239.95);
          _ctrk_t1.push_back(0.5*(plane1tx+plane1ty));
          _ctrk_adc1.push_back(adc1x+adc1y);
          _ctrk_mod1x.push_back(plane1xm);
          _ctrk_x2.push_back(plane2x);
          _ctrk_y2.push_back(plane2y);
          _ctrk_z2.push_back(656.25);
          _ctrk_t2.push_back(0.5*(plane2tx+plane2ty));
          _ctrk_adc2.push_back(adc2x+adc2y);
          _ctrk_mod2x.push_back(plane2xm);
          ntr++;
	        // std::cout << "track " << ntr << std::endl;
	        // std::cout <<  "x y t adc: plane 1 " << plane1x << " " << plane1y << " " <<
	        //   0.5*(plane1tx+plane1ty) << " " << adc1x << " " << adc1y << std::endl;
	        // std::cout <<  "         : plane 2 " << plane2x << " " << plane2y << " " <<
	        //   0.5*(plane2tx+plane2ty) << " " << adc2x << " " << adc2y << std::endl;
        }
      } // i is the first hit with this time
    } // loop over hits
    _nctrks = ntr;
  }  // end if make tracks

  //
  // CRT hits
  //
  if (fkeepCRThits) {
    art::Handle<std::vector<sbn::crt::CRTHit> > crtHitListHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit> > chitlist;
    // art::Handle< std::vector<sbnd::crt::CRTData> > crtStripListHandle;
    // std::vector< art::Ptr<sbnd::crt::CRTData> > striplist;
    if (evt.getByLabel(fCRTHitModuleLabel, crtHitListHandle))  {
      art::fill_ptr_vector(chitlist, crtHitListHandle);
      _nchits = chitlist.size();
    }
    else {
      std::cout << "Failed to get sbn::crt::CRTHit data product." << std::endl;
      _nchits = 0;
    }

    if (_nchits > _max_chits) {
      std::cout << "Available CRT hits are " << _nchits
                << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << _max_chits << "CRT hits." << std::endl;
      _nchits = _max_chits;
    }

    ResetCRTHitsVars(_nchits);

    for (int i = 0; i < _nchits; ++i){
      // int ip = kNotDefined;
      sbnd::CRTPlane ip = sbnd::CRTCommonUtils::GetPlaneIndex(chitlist[i]->tagger);

      _chit_time[i]=chitlist[i]->ts1_ns*0.001;
      if (chitlist[i]->ts1_ns > MAX_INT) {
        _chit_time[i] = 0.001 * (chitlist[i]->ts1_ns - TIME_CORRECTION);
      }

      _chit_x[i] = chitlist[i]->x_pos;
      _chit_y[i] = chitlist[i]->y_pos;
      _chit_z[i] = chitlist[i]->z_pos;
      _chit_plane[i] = ip;
    }
  }

  //
  // CRT tracks
  //
  _ncts = 0;
  if (freadCRTtracks) {
    art::Handle<std::vector<sbn::crt::CRTTrack> > crtTrackListHandle;
    std::vector<art::Ptr<sbn::crt::CRTTrack> > ctrklist;
    if (evt.getByLabel(fCRTTrackModuleLabel, crtTrackListHandle))  {
      art::fill_ptr_vector(ctrklist, crtTrackListHandle);
      _ncts = ctrklist.size();
      if (_ncts > _max_nctrks) _ncts = _max_nctrks;

      ResetCRTTracksVars(_ncts);

      for (int i = 0; i < _ncts; ++i){
        _ct_pes[i] = ctrklist[i]->peshit;
        _ct_time[i] = ctrklist[i]->ts1_ns*0.001;
        if (ctrklist[i]->ts1_ns > MAX_INT) {
          _ct_time[i] = 0.001 * (ctrklist[i]->ts1_ns - TIME_CORRECTION);
        }
        _ct_x1[i] = ctrklist[i]->x1_pos;
        _ct_y1[i] = ctrklist[i]->y1_pos;
        _ct_z1[i] = ctrklist[i]->z1_pos;
        _ct_x2[i] = ctrklist[i]->x2_pos;
        _ct_y2[i] = ctrklist[i]->y2_pos;
        _ct_z2[i] = ctrklist[i]->z2_pos;
      }
    } else {
      std::cout << "Failed to get sbn::crt::CRTTrack data product." << std::endl;
    }
  }


  //
  // Optical Hits
  //
  if (freadOpHits) {
    _nophits = 0;
    size_t previous_nophits = 0;

    // Loop over all the ophits labels
    for (auto ophit_label : fOpHitsModuleLabels) {

      art::Handle<std::vector<recob::OpHit>> ophitListHandle;
      std::vector<art::Ptr<recob::OpHit>> ophitlist;
      if (evt.getByLabel(ophit_label, ophitListHandle)) {
        art::fill_ptr_vector(ophitlist, ophitListHandle);
        _nophits += ophitlist.size();
      }
      else {
        std::cout << "Failed to get recob::OpHit data product." << std::endl;
      }

      if (_nophits > _max_ophits) {
        std::cout << "Available optical hits are " << _nophits << ", which is above the maximum number allowed to store." << std::endl;
        std::cout << "Will only store " << _max_ophits << " optical hits." << std::endl;
        _nophits = _max_ophits;
      }

      ResetOpHitsVars(_nophits);

      for (size_t i = 0; i < ophitlist.size(); ++i) {
        size_t index = previous_nophits + i;
        _ophit_opch[index] = ophitlist.at(i)->OpChannel();
        _ophit_opdet[index] = fGeometryService->OpDetFromOpChannel(ophitlist.at(i)->OpChannel());
        _ophit_peakT[index] = ophitlist.at(i)->PeakTime();
        _ophit_startT[index] = ophitlist.at(i)->StartTime();
        _ophit_riseT[index] = ophitlist.at(i)->RiseTime();
        _ophit_width[index] = ophitlist.at(i)->Width();
        _ophit_area[index] = ophitlist.at(i)->Area();
        _ophit_amplitude[index] = ophitlist.at(i)->Amplitude();
        _ophit_pe[index] = ophitlist.at(i)->PE();
        auto opdet_center = fGeometryService->OpDetGeoFromOpChannel(ophitlist.at(i)->OpChannel()).GetCenter();
        _ophit_opdet_x[index] = opdet_center.X();
        _ophit_opdet_y[index] = opdet_center.Y();
        _ophit_opdet_z[index] = opdet_center.Z();
        auto pd_type = _pd_map.pdType(ophitlist.at(i)->OpChannel());
        if (pd_type == "pmt_coated") {_ophit_opdet_type[index] = kPMTCoated;}
        else if (pd_type == "pmt_uncoated") {_ophit_opdet_type[index] = kPMTUnCoated;}
        else if (pd_type == "xarapuca_vis") {_ophit_opdet_type[index] = kXArapucaVis;}
        else if (pd_type == "xarapuca_vuv") {_ophit_opdet_type[index] = kXArapucaVuv;}
        else {_ophit_opdet_type[index] = kPDNotDefined;}
      }
      previous_nophits = _nophits;
    }
  }

  //
  // pmt hardware trigger
  //

  if (freadpmtTrigger){
    art::Handle<std::vector<sbnd::comm::pmtTrigger> > pmtTriggerListHandle;
    std::vector<art::Ptr<sbnd::comm::pmtTrigger> > pmttriggerlist;

    if (evt.getByLabel(fpmtTriggerModuleLabel, pmtTriggerListHandle)){
      art::fill_ptr_vector(pmttriggerlist, pmtTriggerListHandle);
      ResetPmtTriggerVars( (int)pmttriggerlist[0]->numPassed.size());

      for (int i=0; i < (int)pmttriggerlist[0]->numPassed.size(); i++){
        _pmtTrigger_npmtshigh[i] = pmttriggerlist[0]->numPassed[i];
      }
      _pmtTrigger_maxpassed = pmttriggerlist[0]->maxPMTs;

    }
    else{
      std::cout << "Failed to get sbnd::comm::pmtTrigger data product" << std::endl;
    }
  }

  //
  // PMT Software Trigger
  //
  if (freadpmtSoftTrigger){
    art::Handle<std::vector<sbnd::trigger::pmtSoftwareTrigger>> pmtSoftTriggerListHandle;
    std::vector<art::Ptr<sbnd::trigger::pmtSoftwareTrigger>> pmtsofttriggerlist;

    if (evt.getByLabel(fpmtSoftTriggerModuleLabel, pmtSoftTriggerListHandle)){
      art::fill_ptr_vector(pmtsofttriggerlist,pmtSoftTriggerListHandle);
      ResetPmtSoftTriggerVars();

      auto pmtSoftTriggerMetrics = pmtsofttriggerlist[0];
      _pmtSoftTrigger_foundBeamTrigger = pmtSoftTriggerMetrics->foundBeamTrigger;
      _pmtSoftTrigger_tts = pmtSoftTriggerMetrics->triggerTimestamp;
      _pmtSoftTrigger_promptPE = pmtSoftTriggerMetrics->promptPE;
      _pmtSoftTrigger_prelimPE = pmtSoftTriggerMetrics->prelimPE;
      _pmtSoftTrigger_nAboveThreshold = pmtSoftTriggerMetrics->nAboveThreshold;
      // _pmtSoftTrigger_pmtInfoVec = pmtsofttriggerlist[0]->pmtInfoVec;
    }
    else{
      std::cout << "Failed to get sbnd::trigger::pmtSoftwareTrigger data product" << std::endl;
    }
  }

  //
  // CRT Software Trigger
  //
  if (freadcrtSoftTrigger){
    art::Handle<std::vector<sbndaq::CRTmetric>> crtSoftTriggerListHandle;
    std::vector<art::Ptr<sbndaq::CRTmetric>>    crtsofttriggerlist;
    if (evt.getByLabel(fcrtSoftTriggerModuleLabel, crtSoftTriggerListHandle)){
      art::fill_ptr_vector(crtsofttriggerlist, crtSoftTriggerListHandle);
      ResetCrtSoftTriggerVars();
      auto crtSoftTriggerMetrics = crtsofttriggerlist[0];

      for (int i=0; i<7; i++){
	      _crtSoftTrigger_hitsperplane[i] = crtSoftTriggerMetrics->hitsperplane[i];
      }
    }
    else{
      std::cout << "Failed to get sbndaq::crtMetric data product" << std::endl;
    }

  }

  //
  // Muon tracks 
  //
  _nmuontrks = 0; 
  if (freadMuonTracks){
    art::Handle<std::vector<sbnd::comm::MuonTrack> > muonTrackListHandle;
    std::vector<art::Ptr<sbnd::comm::MuonTrack> > muontrklist;

    if (evt.getByLabel(fMuonTrackModuleLabel, muonTrackListHandle)){
      art::fill_ptr_vector(muontrklist, muonTrackListHandle); 
      _nmuontrks = muontrklist.size();
      ResetMuonTracksVars(_nmuontrks);

      for (int i=0; i < _nmuontrks; i++){ 
        
        _muontrk_t0[i] = muontrklist[i]->t0_us; 
        _muontrk_x1[i] = muontrklist[i]->x1_pos;
        _muontrk_y1[i] = muontrklist[i]->y1_pos;
        _muontrk_z1[i] = muontrklist[i]->z1_pos;
        _muontrk_x2[i] = muontrklist[i]->x2_pos;
        _muontrk_y2[i] = muontrklist[i]->y2_pos;
        _muontrk_z2[i] = muontrklist[i]->z2_pos; 
        _muontrk_theta_xz[i] = muontrklist[i]->theta_xz; 
        _muontrk_theta_yz[i] = muontrklist[i]->theta_yz;
        _muontrk_tpc[i] = muontrklist[i]->tpc; 
        _muontrk_type[i] = muontrklist[i]->type;
      }
      if (freadMuonHits){
        art::FindMany<recob::Hit> muontrkassn(muonTrackListHandle, evt, fMuonTrackModuleLabel);
        ResetMuonHitVars(3000); //estimate of maximum collection hits
        _nmhits = 0;
        for (int i=0; i < _nmuontrks; i++){ 
        std::vector< const recob::Hit*> muonhitsVec = muontrkassn.at(i);
          _nmhits += (muonhitsVec.size()); 
          for (size_t j=0; j<muonhitsVec.size(); j++){
            auto muonhit = muonhitsVec.at(j);
            geo::WireID wireid = muonhit->WireID();
            _mhit_trk.push_back(i);
            _mhit_tpc.push_back(wireid.TPC);
            _mhit_wire.push_back(wireid.Wire);
            _mhit_channel.push_back(muonhit->Channel());
            _mhit_peakT.push_back(muonhit->PeakTime());
            _mhit_charge.push_back(muonhit->Integral());
          }
        }
      }
    }
    else{
      std::cout << "Failed to get sbnd::comm::MuonTrack data product" << std::endl;
    }
  }

  if (fcheckTransparency) {

    _waveform_number.resize(_max_hits*_max_samples, -9999.);
    _adc_on_wire.resize(_max_hits*_max_samples, -9999.);
    _time_for_waveform.resize(_max_hits*_max_samples, -9999.);
    _waveform_integral.resize(_max_hits*_max_samples, -9999.);
    _adc_count_in_waveform.resize(_max_hits*_max_samples, -9999.);

    art::Handle<std::vector<raw::RawDigit>> digitVecHandle;

    bool retVal = evt.getByLabel(fDigitModuleLabel, digitVecHandle);
    if(retVal == true) {
      mf::LogInfo("HitDumper")    << "I got fDigitModuleLabel: "         << fDigitModuleLabel << std::endl;
    } else {
      mf::LogWarning("HitDumper") << "Could not get fDigitModuleLabel: " << fDigitModuleLabel << std::endl;
    }

    int waveform_number_tracker = 0;
    int adc_counter = 1;
    _adc_count = _nhits * (fWindow * 2 + 1);

    // loop over waveforms
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter) {

      //GET THE REFERENCE TO THE CURRENT raw::RawDigit.
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      int channel   = digitVec->Channel();
      auto fDataSize = digitVec->Samples();
      std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.
      rawadc.resize(fDataSize);

      // see if there is a hit on this channel
      for (int ihit = 0; ihit < _nhits; ++ihit) {
        if (_hit_channel[ihit] == channel) {

          int pedestal = (int)digitVec->GetPedestal();
          //UNCOMPRESS THE DATA.
          if (fUncompressWithPed) {
            raw::Uncompress(digitVec->ADCs(), rawadc, pedestal, digitVec->Compression());
          }
          else {
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
          }

          unsigned int bin = _hit_peakT[ihit];
          unsigned int low_edge,high_edge;
          if((int)bin > fWindow and _hit_plane[ihit] == 0) {
            low_edge = bin - (2*fWindow);
          }
          else if ((int)bin>fWindow) {
            low_edge = bin-fWindow;
          }
          else {
            low_edge = 0;
          }
          high_edge = bin + fWindow;
          if (high_edge > (fDataSize-1)) {
            high_edge = fDataSize - 1;
          }
          double integral = 0.0;
          waveform_number_tracker++;
          int counter_for_adc_in_waveform = 0;
          for (size_t ibin = low_edge; ibin <= high_edge; ++ibin) {
            _adc_count_in_waveform[adc_counter] = counter_for_adc_in_waveform;
            counter_for_adc_in_waveform++;
            _waveform_number[adc_counter] = waveform_number_tracker;
            _adc_on_wire[adc_counter] = rawadc[ibin]-pedestal;
            _time_for_waveform[adc_counter] = ibin;
            //std::cout << "DUMP: " << _waveform_number[adc_counter] << " " << _adc_count << " " << _hit_plane[ihit] << " " << _hit_wire[ihit] << " " <<ibin << " " << (rawadc[ibin]-pedestal) << " " << _time_for_waveform[adc_counter] << " " << _adc_on_wire[adc_counter] << std::endl;
            integral+=_adc_on_wire[adc_counter];
            _waveform_integral[adc_counter] = integral;
            adc_counter++;
          }
          std::cout << "DUMP SUM: " << _hit_tpc[ihit] << " " << _hit_plane[ihit] << " " << _hit_wire[ihit] << " " <<  integral << " " << waveform_number_tracker << std::endl;
          _hit_full_integral[ihit] = integral;
        } // if hit channel matches waveform channel
      } //end loop over hits
    }// end loop over waveforms
  }// end if fCheckTrasparency


  if (freadMCParticle){
    //MCParticle
    art::Handle<std::vector<simb::MCParticle>> MCParticleListHandle;
    std::vector<art::Ptr<simb::MCParticle>> MCParticleList;
    if (evt.getByLabel(fMCParticleModuleLabel,MCParticleListHandle)){
      art::fill_ptr_vector(MCParticleList,MCParticleListHandle);
      mcpart_no_primaries = MCParticleList.size();
      ResizeMCParticle(MCParticleList.size()); //Set vectors
      for (size_t iMCPart = 0; iMCPart < MCParticleList.size(); iMCPart++){
        art::Ptr<simb::MCParticle> pPart = MCParticleList[iMCPart]; //get particle pointer
        //Geant info
        mcpart_Mother[iMCPart] = pPart->Mother();
        mcpart_TrackId[iMCPart] = pPart->TrackId();
        mcpart_pdg[iMCPart] = pPart->PdgCode();
        mcpart_status[iMCPart] =  pPart->StatusCode();
        mcpart_process[iMCPart] =  pPart->Process();
        mcpart_endprocess[iMCPart] =  pPart->EndProcess();
        mcpart_Eng[iMCPart] = pPart->E();
        mcpart_EndE[iMCPart] = pPart->EndE();
        mcpart_Mass[iMCPart] = pPart->Mass();
        mcpart_Px[iMCPart] = pPart->Px();
        mcpart_Py[iMCPart] = pPart->Py();
        mcpart_Pz[iMCPart] = pPart->Pz();
        mcpart_P[iMCPart] = pPart->Momentum().Vect().Mag();
        mcpart_StartPointx[iMCPart] = pPart->Vx();
        mcpart_StartPointy[iMCPart] = pPart->Vy();
        mcpart_StartPointz[iMCPart] = pPart->Vz();
        mcpart_StartT[iMCPart] = pPart->T();
        mcpart_EndPointx[iMCPart] = pPart->EndPosition()[0];
        mcpart_EndPointy[iMCPart] = pPart->EndPosition()[1];
        mcpart_EndPointz[iMCPart] = pPart->EndPosition()[2];
        mcpart_EndT[iMCPart] = pPart->EndT();
        mcpart_theta_xz[iMCPart] =  std::atan2(pPart->Px(), pPart->Pz());
        mcpart_theta_yz[iMCPart] =  std::atan2(pPart->Py(), pPart->Pz());
        mcpart_NumberDaughters[iMCPart] = pPart->NumberDaughters();
      }
    }//endif get label
    else {
      std::cout << "Failed to get MCParticle data product." << std::endl;
    }
    //MCTracks
    art::Handle<std::vector<sim::MCTrack>> MCTrackListHandle;
    std::vector<art::Ptr<sim::MCTrack>> MCTrackList;
    if (evt.getByLabel(fMCTrackModuleLabel,MCTrackListHandle)){
      art::fill_ptr_vector(MCTrackList,MCTrackListHandle);
      mctrack_no_primaries = MCTrackList.size();
      ResizeMCTrack(MCTrackList.size());
      for (size_t iMCTrack = 0; iMCTrack < MCTrackList.size(); iMCTrack++){
        art::Ptr<sim::MCTrack> pTrack = MCTrackList[iMCTrack]; //get particle pointer
        mctrack_pdg[iMCTrack] = pTrack->PdgCode();
        mctrack_TrackId[iMCTrack] = pTrack->TrackID();
      }
    }
    else{
      std::cout << "Failed to get MCTrack data product." << std::endl;
    }

    //MCShowers
    art::Handle<std::vector<sim::MCShower>> MCShowerListHandle;
    std::vector<art::Ptr<sim::MCShower>> MCShowerList;
    if (evt.getByLabel(fMCShowerModuleLabel,MCShowerListHandle)){
      art::fill_ptr_vector(MCShowerList,MCShowerListHandle);
      mcshower_no_primaries = MCShowerList.size();
      ResizeMCShower(MCShowerList.size());
      for (size_t iMCShower = 0; iMCShower < MCShowerList.size(); iMCShower++){
        art::Ptr<sim::MCShower> pShower = MCShowerList[iMCShower]; //get particle pointer
        mcshower_pdg[iMCShower] = pShower->PdgCode();
        mcshower_TrackId[iMCShower] = pShower->TrackID();
      }
    }
    else{
      std::cout << "Failed to get MCShower data product." << std::endl;
    }

  } // end read mcparticle
  if (freadTruth){
    //Genie
    int nGeniePrimaries = 0, nMCNeutrinos = 0;
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){
      art::fill_ptr_vector(mclist, mctruthListHandle);
    }else {
      std::cout << "Failed to get Genie data product." << std::endl;
    }

    art::Ptr<simb::MCTruth> mctruth;

      if (!mclist.empty()) {//at least one mc record

        mctruth = mclist[0];

        if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();

    } // if have MC truth
      MF_LOG_DEBUG("HitDumper") << "Expected " << nGeniePrimaries << " GENIE particles";

    //Initially call the number of neutrinos to be stored the number of MCTruth objects.  This is not strictly true i.e. BNB + cosmic overlay but we will count the number of neutrinos later
    nMCNeutrinos = mclist.size();

    ResizeGenie(nGeniePrimaries);
    ResizeMCNeutrino(nMCNeutrinos);

    mcevts_truth = mclist.size();
    //Brailsford 2017/10/16
    //Issue 17917
    //To keep a 1:1 between neutrinos and 'flux' we need the assns
    art::FindManyP<simb::MCFlux> fmFluxNeutrino(mctruthListHandle, evt, fGenieGenModuleLabel);
    // Get GTruth information for scattering code
    art::FindManyP< simb::GTruth > fmgt( mctruthListHandle, evt, fGenieGenModuleLabel );

    if (mcevts_truth > 0){//at least one mc record

      //Brailsford 2017/10/16
      //Issue 17917
      //Loop over every truth in the spill rather than just looking at the first one.
      //Because MCTruth could be a neutrino OR something else (e.g. cosmics) we are going to have to count up how many neutrinos there are
      mcevts_truth = 0;
      for (unsigned int i_mctruth = 0; i_mctruth < mclist.size(); i_mctruth++){
        //fetch an mctruth
        art::Ptr<simb::MCTruth> curr_mctruth = mclist[i_mctruth];
        //Check if it's a neutrino
        if (!curr_mctruth->NeutrinoSet()) continue;

        // Genie Truth association only for the neutrino
        if (fmgt.size()>i_mctruth) {
          std::vector< art::Ptr<simb::GTruth> > mcgtAssn = fmgt.at(i_mctruth);

          nuScatterCode_truth[i_mctruth] = mcgtAssn[0]->fGscatter;
        } else {
          nuScatterCode_truth[i_mctruth] = -1.;
        }

        nuPDG_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().PdgCode();
        ccnc_truth[i_mctruth] = curr_mctruth->GetNeutrino().CCNC();
        mode_truth[i_mctruth] = curr_mctruth->GetNeutrino().Mode();
        Q2_truth[i_mctruth] = curr_mctruth->GetNeutrino().QSqr();
        W_truth[i_mctruth] = curr_mctruth->GetNeutrino().W();
        hitnuc_truth[i_mctruth] = curr_mctruth->GetNeutrino().HitNuc();
        enu_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().E();
        nuvtxx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vx();
        nuvtxy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vy();
        nuvtxz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vz();
        if (curr_mctruth->GetNeutrino().Nu().P()){
          nu_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Px()/curr_mctruth->GetNeutrino().Nu().P();
          nu_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Py()/curr_mctruth->GetNeutrino().Nu().P();
          nu_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Pz()/curr_mctruth->GetNeutrino().Nu().P();
        }
        lep_mom_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().P();
        if (curr_mctruth->GetNeutrino().Lepton().P()){
          lep_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Px()/curr_mctruth->GetNeutrino().Lepton().P();
          lep_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Py()/curr_mctruth->GetNeutrino().Lepton().P();
          lep_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Pz()/curr_mctruth->GetNeutrino().Lepton().P();
        }
        //Brailsford
        //2017/10/17
        //Issue 12918
        //Use the art::Ptr key as the neutrino's unique ID
        nuID_truth[i_mctruth] = curr_mctruth.key();
        //We need to also store N 'flux' neutrinos per event so now check that the FindOneP is valid and, if so, use it!
        if (fmFluxNeutrino.isValid()){
          if (fmFluxNeutrino.at(0).size()>i_mctruth){
          art::Ptr<simb::MCFlux> curr_mcflux = fmFluxNeutrino.at(0).at(i_mctruth);
          tpx_flux[i_mctruth] = curr_mcflux->ftpx;
          tpy_flux[i_mctruth] = curr_mcflux->ftpy;
          tpz_flux[i_mctruth] = curr_mcflux->ftpz;
          tptype_flux[i_mctruth] = curr_mcflux->ftptype;
          }
        }

        //Let's increase the neutrino count
        mcevts_truth++;
      }

      if (mctruth->NeutrinoSet()){
        //genie particles information
        genie_no_primaries = mctruth->NParticles();

        size_t StoreParticles = std::min((size_t) genie_no_primaries, MaxGeniePrimaries);
        if (genie_no_primaries > (int) StoreParticles) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("HitDumper") << "event has "
            << genie_no_primaries << " MC particles, only "
            << StoreParticles << " stored in tree";
        }
        for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
          const simb::MCParticle& part(mctruth->GetParticle(iPart));
          genie_primaries_pdg[iPart]=part.PdgCode();
          genie_Eng[iPart]=part.E();
          genie_Px[iPart]=part.Px();
          genie_Py[iPart]=part.Py();
          genie_Pz[iPart]=part.Pz();
          genie_P[iPart]=part.P();
          genie_status_code[iPart]=part.StatusCode();
          genie_mass[iPart]=part.Mass();
          genie_trackID[iPart]=part.TrackId();
          genie_ND[iPart]=part.NumberDaughters();
          genie_mother[iPart]=part.Mother();
        } // for particle
      } //if neutrino set
    }//if (mcevts_truth)



  }//if (fReadTruth){



  fTree->Fill();

}

 void Hitdumper::beginJob()
 {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("hitdumpertree","analysis tree");
  fTree->Branch("run",&_run,"run/I");
  fTree->Branch("subrun",&_subrun,"subrun/I");
  fTree->Branch("event",&_event,"event/I");
  fTree->Branch("evttime",&_evttime,"evttime/D");
  fTree->Branch("t0",&_t0,"t0/I");

  fTree->Branch("nhits", &_nhits, "nhits/I");
  fTree->Branch("hit_cryostat", &_hit_cryostat);
  fTree->Branch("hit_tpc", &_hit_tpc);
  fTree->Branch("hit_plane", &_hit_plane);
  fTree->Branch("hit_wire", &_hit_wire);
  fTree->Branch("hit_channel", &_hit_channel);
  fTree->Branch("hit_peakT", &_hit_peakT);
  fTree->Branch("hit_charge", &_hit_charge);
  fTree->Branch("hit_ph", &_hit_ph);
  fTree->Branch("hit_width", &_hit_width);
  fTree->Branch("hit_full_integral", &_hit_full_integral);

  if (fcheckTransparency) {
    fTree->Branch("adc_count", &_adc_count,"adc_count/I");
    fTree->Branch("waveform_number", &_waveform_number);
    fTree->Branch("time_for_waveform",&_time_for_waveform);
    fTree->Branch("adc_on_wire", &_adc_on_wire);
    fTree->Branch("waveform_integral", &_waveform_integral);
    fTree->Branch("adc_count_in_waveform", &_adc_count_in_waveform);
  }

  if (fkeepCRTstrips) {
    fTree->Branch("nstrips", &_nstrips, "nstrips/I");
    fTree->Branch("crt_plane", &_crt_plane);
    fTree->Branch("crt_module", &_crt_module);
    fTree->Branch("crt_strip", &_crt_strip);
    fTree->Branch("crt_orient", &_crt_orient);
    fTree->Branch("crt_time", &_crt_time);
    fTree->Branch("crt_adc", &_crt_adc);
    fTree->Branch("crt_pos_x", &_crt_pos_x);
    fTree->Branch("crt_pos_y", &_crt_pos_y);
    fTree->Branch("crt_pos_z", &_crt_pos_z);
  }

  if (fmakeCRTtracks) {
    fTree->Branch("nctrks", &_nctrks, "nctrks/I");
    fTree->Branch("ctrk_x1", &_ctrk_x1);
    fTree->Branch("ctrk_y1", &_ctrk_y1);
    fTree->Branch("ctrk_z1", &_ctrk_z1);
    fTree->Branch("ctrk_t1", &_ctrk_t1);
    fTree->Branch("ctrk_adc1", &_ctrk_adc1);
    fTree->Branch("ctrk_mod1x", &_ctrk_mod1x);
    fTree->Branch("ctrk_x2", &_ctrk_x2);
    fTree->Branch("ctrk_y2", &_ctrk_y2);
    fTree->Branch("ctrk_z2", &_ctrk_z2);
    fTree->Branch("ctrk_t2", &_ctrk_t2);
    fTree->Branch("ctrk_adc2", &_ctrk_adc2);
    fTree->Branch("ctrk_mod2x", &_ctrk_mod2x);
  }
  if (fkeepCRThits) {
    fTree->Branch("nchits", &_nchits, "nchits/I");
    fTree->Branch("chit_x", &_chit_x);
    fTree->Branch("chit_y", &_chit_y);
    fTree->Branch("chit_z", &_chit_z);
    fTree->Branch("chit_time", &_chit_time);
    fTree->Branch("chit_plane", &_chit_plane);
  }
  if (freadCRTtracks) {
    fTree->Branch("ncts", &_ncts, "ncts/I");
    fTree->Branch("ct_x1", &_ct_x1);
    fTree->Branch("ct_y1", &_ct_y1);
    fTree->Branch("ct_z1", &_ct_z1);
    fTree->Branch("ct_x2", &_ct_x2);
    fTree->Branch("ct_y2", &_ct_y2);
    fTree->Branch("ct_z2", &_ct_z2);
    fTree->Branch("ct_time", &_ct_time);
    fTree->Branch("ct_pes", &_ct_pes);
  }

  if (freadOpHits) {
    fTree->Branch("nophits", &_nophits, "nophits/I");
    fTree->Branch("ophit_opch", &_ophit_opch);
    fTree->Branch("ophit_opdet", &_ophit_opdet);
    fTree->Branch("ophit_peakT", &_ophit_peakT);
    fTree->Branch("ophit_startT", &_ophit_startT);
    fTree->Branch("ophit_riseT", &_ophit_riseT);
    fTree->Branch("ophit_width", &_ophit_width);
    fTree->Branch("ophit_area", &_ophit_area);
    fTree->Branch("ophit_amplitude", &_ophit_amplitude);
    fTree->Branch("ophit_pe", &_ophit_pe);
    fTree->Branch("ophit_opdet_x", &_ophit_opdet_x);
    fTree->Branch("ophit_opdet_y", &_ophit_opdet_y);
    fTree->Branch("ophit_opdet_z", &_ophit_opdet_z);
    fTree->Branch("ophit_opdet_type", &_ophit_opdet_type);
  }

  if (freadpmtTrigger){
    fTree->Branch("pmtTrigger_npmtshigh", &_pmtTrigger_npmtshigh);
    fTree->Branch("pmtTrigger_maxpassed", &_pmtTrigger_maxpassed, "pmtTrigger_maxpassed/I");
  }

  if (freadpmtSoftTrigger){
    fTree->Branch("pmtSoftTrigger_foundBeamTrigger", &_pmtSoftTrigger_foundBeamTrigger);
    fTree->Branch("pmtSoftTrigger_tts", &_pmtSoftTrigger_tts);
    fTree->Branch("pmtSoftTrigger_promptPE", &_pmtSoftTrigger_promptPE);
    fTree->Branch("pmtSoftTrigger_prelimPE", &_pmtSoftTrigger_prelimPE);
    fTree->Branch("pmtSoftTrigger_nAboveThreshold", &_pmtSoftTrigger_nAboveThreshold);
    // fTree->Branch("pmtSoftTrigger_", &_pmtSoftTigger_)
  }

  if (freadcrtSoftTrigger){
    fTree->Branch("crtSoftTrigger_hitsperplane", &_crtSoftTrigger_hitsperplane,"crtSoftTrigger_hitsperplane[7]/I");
  }

  if (freadMuonTracks) {
    fTree->Branch("nmuontrks", &_nmuontrks, "nmuontrks/I");
    fTree->Branch("muontrk_t0", &_muontrk_t0);
    fTree->Branch("muontrk_x1", &_muontrk_x1);
    fTree->Branch("muontrk_y1", &_muontrk_y1);
    fTree->Branch("muontrk_z1", &_muontrk_z1);
    fTree->Branch("muontrk_x2", &_muontrk_x2);
    fTree->Branch("muontrk_y2", &_muontrk_y2);
    fTree->Branch("muontrk_z2", &_muontrk_z2);
    fTree->Branch("muontrk_theta_xz", &_muontrk_theta_xz);
    fTree->Branch("muontrk_theta_yz", &_muontrk_theta_yz);
    fTree->Branch("muontrk_tpc", &_muontrk_tpc); 
    fTree->Branch("muontrk_type", &_muontrk_type); 
  }

    if (freadMuonHits) {
    fTree->Branch("nmhits", &_nmhits, "nmhits/I");
    fTree->Branch("mhit_trk", &_mhit_trk);
    fTree->Branch("mhit_tpc", &_mhit_tpc);
    fTree->Branch("mhit_wire", &_mhit_wire); 
    fTree->Branch("mhit_channel", &_mhit_channel);
    fTree->Branch("mhit_peakT", &_mhit_peakT);
    fTree->Branch("mhit_charge", &_mhit_charge); 
  }

  if (freadTruth) {
    fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    fTree->Branch("nuScatterCode_truth",&nuScatterCode_truth);
    fTree->Branch("nuID_truth",&nuID_truth);
    fTree->Branch("nuPDG_truth",&nuPDG_truth);
    fTree->Branch("ccnc_truth",&ccnc_truth);
    fTree->Branch("mode_truth",&mode_truth);
    fTree->Branch("enu_truth",&enu_truth);
    fTree->Branch("Q2_truth",&Q2_truth);
    fTree->Branch("W_truth",&W_truth);
    fTree->Branch("hitnuc_truth",&hitnuc_truth);
    fTree->Branch("nuvtxx_truth",&nuvtxx_truth);
    fTree->Branch("nuvtxy_truth",&nuvtxy_truth);
    fTree->Branch("nuvtxz_truth",&nuvtxz_truth);
    fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth);
    fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth);
    fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth);
    fTree->Branch("lep_mom_truth",&lep_mom_truth);
    fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth);
    fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth);
    fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth);

    fTree->Branch("tpx_flux",&tpx_flux);
    fTree->Branch("tpy_flux",&tpy_flux);
    fTree->Branch("tpz_flux",&tpz_flux);
    fTree->Branch("tptype_flux",&tptype_flux);

    fTree->Branch("genie_no_primaries",&genie_no_primaries);
    fTree->Branch("genie_primaries_pdg",&genie_primaries_pdg);
    fTree->Branch("genie_Eng",&genie_Eng);
    fTree->Branch("genie_Px",&genie_Px);
    fTree->Branch("genie_Py",&genie_Py);
    fTree->Branch("genie_Pz",&genie_Pz);
    fTree->Branch("genie_P",&genie_P);
    fTree->Branch("genie_status_code",&genie_status_code);
    fTree->Branch("genie_mass",&genie_mass);
    fTree->Branch("genie_trackID",&genie_trackID);
    fTree->Branch("genie_ND",&genie_ND);
    fTree->Branch("genie_mother",&genie_mother);
  }
  if (freadMCParticle){
    //MCParticle
    fTree->Branch("mcpart_pdg",&mcpart_pdg);
    fTree->Branch("mcpart_status",&mcpart_status);    
    fTree->Branch("mcpart_process",&mcpart_process);  
    fTree->Branch("mcpart_endprocess",&mcpart_endprocess);  
    fTree->Branch("mcpart_Eng",&mcpart_Eng);
    fTree->Branch("mcpart_EndE",&mcpart_EndE);
    fTree->Branch("mcpart_Mass",&mcpart_Mass);
    fTree->Branch("mcpart_Px",&mcpart_Px);
    fTree->Branch("mcpart_Py",&mcpart_Py);
    fTree->Branch("mcpart_Pz",&mcpart_Pz);
    fTree->Branch("mcpart_P",&mcpart_P);
    fTree->Branch("mcpart_StartPointx",&mcpart_StartPointx);
    fTree->Branch("mcpart_StartPointy",&mcpart_StartPointy);
    fTree->Branch("mcpart_StartPointz",&mcpart_StartPointz);
    fTree->Branch("mcpart_StartT",&mcpart_StartT);
    fTree->Branch("mcpart_EndT",&mcpart_EndT);
    fTree->Branch("mcpart_EndPointx",&mcpart_EndPointx);
    fTree->Branch("mcpart_EndPointy",&mcpart_EndPointy);
    fTree->Branch("mcpart_EndPointz",&mcpart_EndPointz);         
    fTree->Branch("mcpart_theta_xz",&mcpart_theta_xz);
    fTree->Branch("mcpart_theta_yz",&mcpart_theta_yz);   
    fTree->Branch("mcpart_NumberDaughters",&mcpart_NumberDaughters);
    fTree->Branch("mcpart_TrackId",&mcpart_TrackId);
    fTree->Branch("mcpart_Mother",&mcpart_Mother);

    //MCTrack info
    fTree->Branch("mctrack_no_primaries",&mctrack_no_primaries);
    fTree->Branch("mctrack_pdg",&mctrack_pdg);                        
    fTree->Branch("mctrack_TrackId",&mctrack_TrackId);

    //MCShower info
    fTree->Branch("mcshower_no_primaries",&mcshower_no_primaries);
    fTree->Branch("mcshower_pdg",&mcshower_pdg);                        
    fTree->Branch("mcshower_TrackId",&mcshower_TrackId);
  }

  if (fsavePOTInfo) {
    _sr_tree = tfs->make<TTree>("pottree","");
    _sr_tree->Branch("run", &_sr_run, "run/I");
    _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
    _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
    _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
    _sr_tree->Branch("pot", &_sr_pot, "pot/D");
  }

}


void Hitdumper::ResetWireHitsVars(int n) {
  _hit_cryostat.assign(n, DEFAULT_VALUE);
  _hit_tpc.assign(n, DEFAULT_VALUE);
  _hit_plane.assign(n, DEFAULT_VALUE);
  _hit_wire.assign(n, DEFAULT_VALUE);
  _hit_channel.assign(n, DEFAULT_VALUE);
  _hit_peakT.assign(n, DEFAULT_VALUE);
  _hit_charge.assign(n, DEFAULT_VALUE);
  _hit_ph.assign(n, DEFAULT_VALUE);
  _hit_width.assign(n, DEFAULT_VALUE);
  _hit_full_integral.assign(n, DEFAULT_VALUE);
}

void Hitdumper::ResetCRTStripsVars(int n) {
  _crt_plane.clear();
  _crt_module.clear();
  _crt_strip.clear();
  _crt_orient.clear();
  _crt_time.clear();
  _crt_adc.clear();
  _crt_pos_x.clear();
  _crt_pos_y.clear();
  _crt_pos_z.clear();

  _crt_plane.reserve(n);
  _crt_module.reserve(n);
  _crt_strip.reserve(n);
  _crt_orient.reserve(n);
  _crt_time.reserve(n);
  _crt_adc.reserve(n);
  _crt_pos_x.reserve(n);
  _crt_pos_y.reserve(n);
  _crt_pos_z.reserve(n);
}

void Hitdumper::ResetCRTCustomTracksVars(int n) {
  _ctrk_x1.clear();
  _ctrk_y1.clear();
  _ctrk_z1.clear();
  _ctrk_t1.clear();
  _ctrk_adc1.clear();
  _ctrk_mod1x.clear();
  _ctrk_x2.clear();
  _ctrk_y2.clear();
  _ctrk_z2.clear();
  _ctrk_t2.clear();
  _ctrk_adc2.clear();
  _ctrk_mod2x.clear();

  _ctrk_x1.reserve(n);
  _ctrk_y1.reserve(n);
  _ctrk_z1.reserve(n);
  _ctrk_t1.reserve(n);
  _ctrk_adc1.reserve(n);
  _ctrk_mod1x.reserve(n);
  _ctrk_x2.reserve(n);
  _ctrk_y2.reserve(n);
  _ctrk_z2.reserve(n);
  _ctrk_t2.reserve(n);
  _ctrk_adc2.reserve(n);
  _ctrk_mod2x.reserve(n);
}

void Hitdumper::ResetCRTTracksVars(int n) {
  _ct_x1.assign(n, DEFAULT_VALUE);
  _ct_y1.assign(n, DEFAULT_VALUE);
  _ct_z1.assign(n, DEFAULT_VALUE);
  _ct_time.assign(n, DEFAULT_VALUE);
  _ct_pes.assign(n, DEFAULT_VALUE);
  _ct_x2.assign(n, DEFAULT_VALUE);
  _ct_y2.assign(n, DEFAULT_VALUE);
  _ct_z2.assign(n, DEFAULT_VALUE);
}

void Hitdumper::ResetCRTHitsVars(int n) {
  _chit_plane.assign(n, DEFAULT_VALUE);
  _chit_time.assign(n, DEFAULT_VALUE);
  _chit_x.assign(n, DEFAULT_VALUE);
  _chit_y.assign(n, DEFAULT_VALUE);
  _chit_z.assign(n, DEFAULT_VALUE);
}

void Hitdumper::ResetOpHitsVars(int n) {
  _ophit_opch.resize(n, DEFAULT_VALUE);
  _ophit_opdet.resize(n, DEFAULT_VALUE);
  _ophit_peakT.resize(n, DEFAULT_VALUE);
  _ophit_startT.resize(n, DEFAULT_VALUE);
  _ophit_riseT.resize(n, DEFAULT_VALUE);
  _ophit_width.resize(n, DEFAULT_VALUE);
  _ophit_area.resize(n, DEFAULT_VALUE);
  _ophit_amplitude.resize(n, DEFAULT_VALUE);
  _ophit_pe.resize(n, DEFAULT_VALUE);
  _ophit_opdet_x.resize(n, DEFAULT_VALUE);
  _ophit_opdet_y.resize(n, DEFAULT_VALUE);
  _ophit_opdet_z.resize(n, DEFAULT_VALUE);
  _ophit_opdet_type.resize(n, DEFAULT_VALUE);
}

void Hitdumper::ResetPmtTriggerVars(int n){
  _pmtTrigger_npmtshigh.assign(n, DEFAULT_VALUE);
  _pmtTrigger_maxpassed = 0;
}

void Hitdumper::ResetPmtSoftTriggerVars(){
  _pmtSoftTrigger_foundBeamTrigger = false;
  _pmtSoftTrigger_tts = 0;
  _pmtSoftTrigger_promptPE = 0;
  _pmtSoftTrigger_prelimPE = 0;
  _pmtSoftTrigger_nAboveThreshold = 0;
}

void Hitdumper::ResetCrtSoftTriggerVars(){
  for (int i=0; i<7; i++){
    _crtSoftTrigger_hitsperplane[i] = 0;
  }
}

void Hitdumper::ResetMuonTracksVars(int n){
  _muontrk_t0.assign(n, DEFAULT_VALUE);
  _muontrk_x1.assign(n, DEFAULT_VALUE);
  _muontrk_y1.assign(n, DEFAULT_VALUE);
  _muontrk_z1.assign(n, DEFAULT_VALUE);
  _muontrk_x2.assign(n, DEFAULT_VALUE);
  _muontrk_y2.assign(n, DEFAULT_VALUE);
  _muontrk_z2.assign(n, DEFAULT_VALUE); 
  _muontrk_theta_xz.assign(n, DEFAULT_VALUE); 
  _muontrk_theta_yz.assign(n, DEFAULT_VALUE);
  _muontrk_tpc.assign(n, DEFAULT_VALUE);
  _muontrk_type.assign(n, DEFAULT_VALUE);
}

void Hitdumper::ResetMuonHitVars(int n){
  _mhit_trk.clear(); 
  _mhit_tpc.clear();
  _mhit_wire.clear();
  _mhit_channel.clear();
  _mhit_peakT.clear();
  _mhit_charge.clear();

  _mhit_trk.reserve(n);
  _mhit_tpc.reserve(n);
  _mhit_wire.reserve(n);
  _mhit_channel.reserve(n);
  _mhit_peakT.reserve(n);
  _mhit_charge.reserve(n);
}

void Hitdumper::ResetVars() {


  _run = -99999;
  _subrun = -99999;
  _event = -99999;
  _evttime = -99999;
  _t0 = -99999;

  mcevts_truth = 0;
  genie_no_primaries = 0;
  mcpart_no_primaries = 0;
  mctrack_no_primaries = 0;
  mcshower_no_primaries = 0;

}

void Hitdumper::ResizeMCNeutrino(int nNeutrinos) {

  //min size is 1, to guarantee an address
  MaxMCNeutrinos = (size_t) std::max(nNeutrinos, 1);

  nuScatterCode_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuID_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuPDG_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  ccnc_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  mode_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  enu_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  Q2_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  W_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  hitnuc_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_mom_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  //Also resize the flux information here as it's a 1:1 with the MCNeutrino
  tpx_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tpy_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tpz_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tptype_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);

}

void Hitdumper::ResizeGenie(int nPrimaries) {

  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);

  genie_primaries_pdg.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Eng.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Px.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Py.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Pz.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_P.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_status_code.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_mass.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_trackID.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_ND.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_mother.assign(MaxGeniePrimaries, DEFAULT_VALUE);

}

void Hitdumper::ResizeMCParticle(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxMCParticles = (size_t) std::max(nParticles, 1);
              
  mcpart_pdg.assign(MaxMCParticles,DEFAULT_VALUE);              
  mcpart_status.assign(MaxMCParticles,DEFAULT_VALUE); 
  mcpart_process.assign(MaxMCParticles,"Dummy"); 
  mcpart_endprocess.assign(MaxMCParticles,"Dummy");           
  mcpart_Eng.assign(MaxMCParticles,DEFAULT_VALUE);              
  mcpart_EndE.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Mass.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Px.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Py.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Pz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_P.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointx.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointy.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartT.assign(MaxMCParticles,DEFAULT_VALUE);  
  mcpart_EndT.assign(MaxMCParticles,DEFAULT_VALUE);          
  mcpart_EndPointx.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_EndPointy.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_EndPointz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_theta_xz.assign(MaxMCParticles,DEFAULT_VALUE);    
  mcpart_theta_yz.assign(MaxMCParticles,DEFAULT_VALUE);    
  mcpart_NumberDaughters.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_TrackId.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Mother.assign(MaxMCParticles,DEFAULT_VALUE);

}

void Hitdumper::ResizeMCTrack(int nTracks) {
  MaxMCTracks = (size_t) std::max(nTracks,1);

  mctrack_pdg.assign(MaxMCTracks,DEFAULT_VALUE);
  mctrack_TrackId.assign(MaxMCTracks,DEFAULT_VALUE);
}

void Hitdumper::ResizeMCShower(int nShowers) {
  MaxMCShowers = (size_t) std::max(nShowers,1);

  mcshower_pdg.assign(MaxMCShowers,DEFAULT_VALUE);
  mcshower_TrackId.assign(MaxMCShowers,DEFAULT_VALUE);
}


void Hitdumper::beginSubRun(art::SubRun const& sr) {

  if (!fsavePOTInfo) {
    return;
  }

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(fGenieGenModuleLabel, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }

  _sr_tree->Fill();

}

DEFINE_ART_MODULE(Hitdumper)

#endif // Hitdumper_Module
