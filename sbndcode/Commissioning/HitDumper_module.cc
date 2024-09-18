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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
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
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Commissioning/MuonTrack.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/CRTmetric.hh"

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
  /// Resets crt strip hit tree variables
  void ResetCRTStripHitVars();
  /// Resets crt tracks tree variables
  void ResetCRTTracksVars();
  /// Resets crt spacepoint tree variables
  void ResetCRTSpacePointVars();
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
  std::vector<double> _hit_ph;                  ///< Hit pulse height
  std::vector<double> _hit_width;               ///< Hit width
  std::vector<double> _hit_full_integral;       ///< Hit charge integral
  std::vector<int>    _waveform_number;         ///< Number for each waveform, to allow for searching
  std::vector<double> _adc_on_wire;             ///< ADC on wire to draw waveform
  std::vector<int>    _time_for_waveform;       ///<Time for waveform to plot
  int                 _adc_count;               ///<Used for plotting waveforms
  std::vector<int>    _waveform_integral;       ///<Used to see progression of the waveform integral
  std::vector<int>    _adc_count_in_waveform;   ///<Used to view all waveforms on a hitplane together


  // CRT strip variables
  uint _n_crt_strip_hits;                          ///< Number of CRT strip hits
  std::vector<int> _crt_strip_hit_tagger;          ///< CRT strip hit tagger enum
  std::vector<int> _crt_strip_hit_module;          ///< CRT strip hit module
  std::vector<int> _crt_strip_hit_channel;         ///< CRT strip hit channel
  std::vector<int> _crt_strip_hit_orient;          ///< CRT strip hit orientation (0 for y (horizontal) and 1 for x (vertical))
  std::vector<uint> _crt_strip_hit_t0;             ///< CRT strip hit t0
  std::vector<uint> _crt_strip_hit_t1;             ///< CRT strip hit t1
  std::vector<uint> _crt_strip_hit_adc1;           ///< CRT strip hit adc1
  std::vector<uint> _crt_strip_hit_adc2;           ///< CRT strip hit adc2
  std::vector<double> _crt_strip_hit_pos;          ///< CRT strip hit position
  std::vector<double> _crt_strip_hit_pos_err;      ///< CRT strip hit position error

  // CRT space point variables (3D hits)
  uint _n_crt_space_points;                         ///< Number of CRT SpacePoints
  std::vector<double> _crt_space_point_x;           ///< CRT SpacePoint x
  std::vector<double> _crt_space_point_y;           ///< CRT SpacePoint y
  std::vector<double> _crt_space_point_z;           ///< CRT SpacePoint z
  std::vector<double> _crt_space_point_time;        ///< CRT SpacePoint time
  std::vector<double> _crt_space_point_x_err;       ///< CRT SpacePoint x error
  std::vector<double> _crt_space_point_y_err;       ///< CRT SpacePoint y error
  std::vector<double> _crt_space_point_z_err;       ///< CRT SpacePoint z error
  std::vector<double> _crt_space_point_time_err;    ///< CRT SpacePoint time error
  std::vector<double> _crt_space_point_pe;          ///< CRT SpacePoint PE
  std::vector<int> _crt_space_point_tagger;         ///< CRT SpacePoint tagger enum
  std::vector<int> _crt_space_point_nhits;          ///< CRT SpacePoint nhits

  // CRT track variables
  uint _n_crt_tracks;                            ///< Number of CRT tracks
  std::vector<double> _crt_track_time;         ///< CRT track time
  std::vector<double> _crt_track_pes;          ///< CRT track PEs
  std::vector<double> _crt_track_x1;           ///< CRT track x1
  std::vector<double> _crt_track_y1;           ///< CRT track y1
  std::vector<double> _crt_track_z1;           ///< CRT track z1
  std::vector<double> _crt_track_x2;           ///< CRT track x2
  std::vector<double> _crt_track_y2;           ///< CRT track y2
  std::vector<double> _crt_track_z2;           ///< CRT track z2
  std::vector<int> _crt_track_tagger1;         ///< Hit 1 tagger
  std::vector<int> _crt_track_tagger2;         ///< Hit 2 tagger
  std::vector<double> _crt_track_theta;        ///< CRT track theta
  std::vector<double> _crt_track_phi;          ///< CRT track phi
  std::vector<double> _crt_track_length;       ///< CRT track length

  // Optical hit variables
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

  // PMT hardware trigger variables
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
  std::vector<int>    _mhit_plane;           ///< Wire plane where the hit belongs to
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
  uint _max_crt_strip_hits;         ///< maximum number of CRT strip hits (to be set via fcl)
  uint _max_crt_space_points;       ///< maximum number of CRT space points (to be set via fcl)
  uint _max_crt_tracks;                 ///< maximum number of CRT tracks (to be set via fcl)

  std::string fHitsModuleLabel;     ///< Label for Hit dataproduct (to be set via fcl)
  std::string fLArG4ModuleLabel;    ///< Label for LArG4 dataproduct (to be set via fcl)
  std::string fCRTStripHitModuleLabel; ///< Label for CRTStrip dataproduct (to be set via fcl)
  std::string fCRTSpacePointModuleLabel;   ///< Label for CRTSpacePoint dataproduct (to be set via fcl)
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

  bool fKeepCRTSpacePoints;///< Keep the CRT spacepoints (to be set via fcl)
  bool fKeepCRTStripHits;  ///< Keep the CRT strips (to be set via fcl)
  bool fKeepCRTTracks;     ///< Keep the CRT tracks (to be set via fcl)
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

  sbnd::crt::CRTGeoAlg fCRTGeoAlg;

  geo::GeometryCore const* fGeometryService;
  // detinfo::ElecClock fTrigClock;
  art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
  const geo::AuxDetGeometry* fAuxDetGeo;
  const geo::AuxDetGeometryCore* fAuxDetGeoCore;

  const int MAX_INT = std::numeric_limits<int>::max();
  const long int TIME_CORRECTION = (long int) std::numeric_limits<int>::max() * 2;
};


Hitdumper::Hitdumper(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
  , fCRTGeoAlg(pset.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
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
  _max_crt_strip_hits = p.get<uint>("MaxCRTStripHits", 5000);
  _max_crt_space_points = p.get<uint>("MaxCRTSpacePoints", 500);
  _max_crt_tracks = p.get<int>("MaxCRTTracks", 100);

  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
  fDigitModuleLabel    = p.get<std::string>("DigitModuleLabel", "daq");
  fLArG4ModuleLabel    = p.get<std::string>("LArG4ModuleLabel", "largeant");
  fCRTStripHitModuleLabel = p.get<std::string>("CRTStripHitModuleLabel", "crtstrips");
  fCRTSpacePointModuleLabel   = p.get<std::string>("CRTSpacePointModuleLabel", "crtspacepoints");
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel", "crttracks");
  fOpHitsModuleLabels  = p.get<std::vector<std::string>>("OpHitsModuleLabel");
  fpmtTriggerModuleLabel = p.get<std::string>("pmtTriggerModuleLabel", "pmttriggerproducer");
  fpmtSoftTriggerModuleLabel = p.get<std::string>("pmtSoftTriggerModuleLabel", "pmtSoftwareTrigger");
  fcrtSoftTriggerModuleLabel = p.get<std::string>("crtSoftTriggerModuleLabel", "MetricProducer");
  fMuonTrackModuleLabel  = p.get<std::string>("MuonTrackModuleLabel", "MuonTrackProducer");
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel", "generator");
  fMCParticleModuleLabel    = p.get<std::string>("MCParticleModuleLabel ", "largeant");
  fMCTrackModuleLabel    = p.get<std::string>("MCTrackModuleLabel ", "mcreco");
  fMCShowerModuleLabel    = p.get<std::string>("MCShowerModuleLabel ", "mcreco");

  fKeepCRTSpacePoints = p.get<bool>("KeepCRTSpacePoints",true);
  fKeepCRTStripHits   = p.get<bool>("KeepCRTStripHits", true);
  fKeepCRTTracks      = p.get<bool>("KeepCRTTracks",true);
  freadOpHits         = p.get<bool>("readOpHits",true);
  freadpmtTrigger     = p.get<bool>("readpmtTrigger",true);
  freadpmtSoftTrigger = p.get<bool>("readpmtSoftTrigger",true);
  freadcrtSoftTrigger = p.get<bool>("readcrtSoftTrigger",false);
  freadMuonTracks     = p.get<bool>("readMuonTracks",true);
  freadMuonHits       = p.get<bool>("readMuonHits",false);
  fcheckTransparency  = p.get<bool>("checkTransparency",false);
  freadTruth          = p.get<bool>("readTruth",true);
  freadMCParticle     = p.get<bool>("readMCParticle",false);
  fsavePOTInfo        = p.get<bool>("savePOTinfo",true);
  fUncompressWithPed  = p.get<bool>("UncompressWithPed",false);

  fWindow             = p.get<int>("window",100);
  fKeepTaggerTypes    = p.get<std::vector<int>>("KeepTaggerTypes");

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
  art::Handle<std::vector<sbnd::crt::CRTStripHit>> crtStripHitHandle;
  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> crtStripHitVector;
  if (evt.getByLabel(fCRTStripHitModuleLabel, crtStripHitHandle))  {
    art::fill_ptr_vector(crtStripHitVector, crtStripHitHandle);
    _n_crt_strip_hits = crtStripHitVector.size();
  }
  else
    std::cout << "Failed to get sbnd::crt::CRTStripHit data product ("<<fCRTStripHitModuleLabel<<")." << std::endl;

  if (_n_crt_strip_hits > _max_crt_strip_hits) _n_crt_strip_hits = _max_crt_strip_hits;

  ResetCRTStripHitVars();

  for(uint i = 0; i < _n_crt_strip_hits; ++i)
    {
      const art::Ptr<sbnd::crt::CRTStripHit> stripHit = crtStripHitVector[i];

      _crt_strip_hit_tagger.push_back(fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel()));
      _crt_strip_hit_module.push_back(fCRTGeoAlg.GetModule(stripHit->Channel()).adID);
      _crt_strip_hit_channel.push_back(stripHit->Channel());
      _crt_strip_hit_orient.push_back(fCRTGeoAlg.ChannelToOrientation(stripHit->Channel()));
      _crt_strip_hit_t0.push_back(stripHit->Ts0());
      _crt_strip_hit_t1.push_back(stripHit->Ts1());
      _crt_strip_hit_adc1.push_back(stripHit->ADC1());
      _crt_strip_hit_adc2.push_back(stripHit->ADC2());
      _crt_strip_hit_pos.push_back(stripHit->Pos());
      _crt_strip_hit_pos_err.push_back(stripHit->Error());
    }

  //
  // CRT SpacePoints (3D hits)
  //
  if (fKeepCRTSpacePoints) {
    art::Handle<std::vector<sbnd::crt::CRTSpacePoint> > crtSPHandle;
    std::vector<art::Ptr<sbnd::crt::CRTSpacePoint> > crtSPVector;
    if (evt.getByLabel(fCRTSpacePointModuleLabel, crtSPHandle))  {
      art::fill_ptr_vector(crtSPVector, crtSPHandle);
      _n_crt_space_points = crtSPVector.size();
    }
    else {
      std::cout << "Failed to get sbnd::crt::CRTSpacePoint data product ("<<fCRTSpacePointModuleLabel<<")." << std::endl;
      _n_crt_space_points = 0;
    }

    art::FindOneP<sbnd::crt::CRTCluster> spToCluster(crtSPHandle, evt, fCRTSpacePointModuleLabel);

    if (_n_crt_space_points > _max_crt_space_points) {
      std::cout << "Available CRT space points are " << _n_crt_space_points
                << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << _max_crt_space_points << " CRT space points." << std::endl;
      _n_crt_space_points = _max_crt_space_points;
    }

    ResetCRTSpacePointVars();

    for (uint i = 0; i < _n_crt_space_points; ++i){
      const art::Ptr<sbnd::crt::CRTSpacePoint> spacePoint = crtSPVector[i];

      if(!spacePoint->Complete())
	continue;

      const art::Ptr<sbnd::crt::CRTCluster> cluster = spToCluster.at(spacePoint.key());

      _crt_space_point_x.push_back(spacePoint->X());
      _crt_space_point_y.push_back(spacePoint->Y());
      _crt_space_point_z.push_back(spacePoint->Z());
      _crt_space_point_time.push_back(spacePoint->Time());
      _crt_space_point_x_err.push_back(spacePoint->XErr());
      _crt_space_point_y_err.push_back(spacePoint->YErr());
      _crt_space_point_z_err.push_back(spacePoint->ZErr());
      _crt_space_point_time_err.push_back(spacePoint->TimeErr());
      _crt_space_point_pe.push_back(spacePoint->PE());
      _crt_space_point_tagger.push_back(cluster->Tagger());
      _crt_space_point_nhits.push_back(cluster->NHits());
    }
  }

  //
  // CRT tracks
  //

  _n_crt_tracks = 0;
  if (fKeepCRTTracks) {

    art::Handle<std::vector<sbnd::crt::CRTTrack>> crtTrackListHandle;
    std::vector<art::Ptr<sbnd::crt::CRTTrack>> ctrklist;

    if (evt.getByLabel(fCRTTrackModuleLabel, crtTrackListHandle))  {
      art::fill_ptr_vector(ctrklist, crtTrackListHandle);
      _n_crt_tracks = ctrklist.size();
      if (_n_crt_tracks > _max_crt_tracks) _n_crt_tracks = _max_crt_tracks;

      ResetCRTTracksVars();

      for (uint i = 0; i < _n_crt_tracks; ++i){
        const art::Ptr<sbnd::crt::CRTTrack> crttrack=ctrklist[i];
        _crt_track_pes.push_back(crttrack->PE());
        _crt_track_time.push_back(crttrack->Time());

	const geo::Point_t start = crttrack->Start();
	const geo::Point_t end   = crttrack->End();

        _crt_track_x1.push_back(start.X());
        _crt_track_y1.push_back(start.Y());
        _crt_track_z1.push_back(start.Z());
        _crt_track_x2.push_back(end.X());
        _crt_track_y2.push_back(end.Y());
        _crt_track_z2.push_back(end.Z());

  	_crt_track_tagger1.push_back(fCRTGeoAlg.WhichTagger(start.X(),start.Y(),start.Z()));
	_crt_track_tagger2.push_back(fCRTGeoAlg.WhichTagger(end.X(),end.Y(),end.Z()));
     
	_crt_track_theta.push_back(crttrack->Theta()*(180.0/M_PI));
	_crt_track_phi.push_back(crttrack->Phi()*(180.0/M_PI));
	_crt_track_length.push_back(crttrack->Length());
      }
    } else {
      std::cout << "Failed to get sbnd::crt::CRTTrack data product ("<<fCRTTrackModuleLabel<<")." << std::endl;
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
        std::cout << "Failed to get recob::OpHit data product: " << ophit_label << std::endl;
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
        else if (pd_type == "xarapuca_vuv") {
          _ophit_opdet_type[index] = kXArapucaVuv;
          //std::cout<<"XA VUV: "<< pd_type <<std::endl;
        }
        else {
          _ophit_opdet_type[index] = kPDNotDefined;
        }
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
      _pmtSoftTrigger_tts = pmtSoftTriggerMetrics->trig_ts;
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
            _mhit_plane.push_back(wireid.Plane);
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

  if (fKeepCRTStripHits) {
    fTree->Branch("n_crt_strip_hits", &_n_crt_strip_hits);
    fTree->Branch("crt_strip_hit_tagger", &_crt_strip_hit_tagger);
    fTree->Branch("crt_strip_hit_module", &_crt_strip_hit_module);
    fTree->Branch("crt_strip_hit_channel", &_crt_strip_hit_channel);
    fTree->Branch("crt_strip_hit_orient", &_crt_strip_hit_orient);
    fTree->Branch("crt_strip_hit_t0", &_crt_strip_hit_t0);
    fTree->Branch("crt_strip_hit_t1", &_crt_strip_hit_t1);
    fTree->Branch("crt_strip_hit_adc1", &_crt_strip_hit_adc1);
    fTree->Branch("crt_strip_hit_adc2", &_crt_strip_hit_adc2);
    fTree->Branch("crt_strip_hit_pos", &_crt_strip_hit_pos);
    fTree->Branch("crt_strip_hit_pos_err", &_crt_strip_hit_pos_err);
  }

  if (fKeepCRTSpacePoints) {
    fTree->Branch("n_crt_space_points", &_n_crt_space_points);
    fTree->Branch("crt_space_point_x", &_crt_space_point_x);
    fTree->Branch("crt_space_point_y", &_crt_space_point_y);
    fTree->Branch("crt_space_point_z", &_crt_space_point_z);
    fTree->Branch("crt_space_point_time", &_crt_space_point_time);
    fTree->Branch("crt_space_point_x_err", &_crt_space_point_x_err);
    fTree->Branch("crt_space_point_y_err", &_crt_space_point_y_err);
    fTree->Branch("crt_space_point_z_err", &_crt_space_point_z_err);
    fTree->Branch("crt_space_point_time_err", &_crt_space_point_time_err);
    fTree->Branch("crt_space_point_pe", &_crt_space_point_pe);
    fTree->Branch("crt_space_point_tagger", &_crt_space_point_tagger);
    fTree->Branch("crt_space_point_nhits", &_crt_space_point_nhits);    
  }
  if (fKeepCRTTracks) {
    fTree->Branch("n_crt_tracks", &_n_crt_tracks, "n_crt_tracks/I");
    fTree->Branch("crt_track_x1", &_crt_track_x1);
    fTree->Branch("crt_track_y1", &_crt_track_y1);
    fTree->Branch("crt_track_z1", &_crt_track_z1);
    fTree->Branch("crt_track_x2", &_crt_track_x2);
    fTree->Branch("crt_track_y2", &_crt_track_y2);
    fTree->Branch("crt_track_z2", &_crt_track_z2);
    fTree->Branch("crt_track_time", &_crt_track_time);
    fTree->Branch("crt_track_pes", &_crt_track_pes);
    fTree->Branch("crt_track_tagger1", &_crt_track_tagger1);
    fTree->Branch("crt_track_tagger2", &_crt_track_tagger2);
    fTree->Branch("crt_track_theta", &_crt_track_theta);
    fTree->Branch("crt_track_phi", &_crt_track_phi);
    fTree->Branch("crt_track_length", &_crt_track_length);
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
    fTree->Branch("mhit_plane", &_mhit_plane);
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

void Hitdumper::ResetCRTStripHitVars() {
  _crt_strip_hit_tagger.clear();
  _crt_strip_hit_module.clear();
  _crt_strip_hit_channel.clear();
  _crt_strip_hit_orient.clear();
  _crt_strip_hit_t0.clear();
  _crt_strip_hit_t1.clear();
  _crt_strip_hit_adc1.clear();
  _crt_strip_hit_adc2.clear();
  _crt_strip_hit_pos.clear();
  _crt_strip_hit_pos_err.clear();
}

void Hitdumper::ResetCRTTracksVars() {
  _crt_track_x1.clear();
  _crt_track_y1.clear();
  _crt_track_z1.clear();
  _crt_track_time.clear();
  _crt_track_pes.clear();
  _crt_track_x2.clear();
  _crt_track_y2.clear();
  _crt_track_z2.clear();
  _crt_track_tagger1.clear();
  _crt_track_tagger2.clear();
  _crt_track_theta.clear();
  _crt_track_phi.clear();
  _crt_track_length.clear();
}

void Hitdumper::ResetCRTSpacePointVars() {
  _crt_space_point_x.clear();
  _crt_space_point_y.clear();
  _crt_space_point_z.clear();
  _crt_space_point_time.clear();
  _crt_space_point_x_err.clear();
  _crt_space_point_y_err.clear();
  _crt_space_point_z_err.clear();
  _crt_space_point_time_err.clear();
  _crt_space_point_pe.clear();
  _crt_space_point_tagger.clear();
  _crt_space_point_nhits.clear();
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
  _mhit_plane.clear();
  _mhit_wire.clear();
  _mhit_channel.clear();
  _mhit_peakT.clear();
  _mhit_charge.clear();

  _mhit_trk.reserve(n);
  _mhit_tpc.reserve(n);
  _mhit_plane.reserve(n);
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
