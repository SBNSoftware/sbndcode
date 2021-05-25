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
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
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
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// Truth includes
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

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

enum CRTPos {
  kNotDefined = -1,   ///< Not defined
  kBot = 0,           ///< Bot
  kFaceFront = 1,     ///< FaceFront
  kFaceBack,          ///< FaceBack
  kSideLeft,          ///< SideLeft
  kSideRight,         ///< SideRight
  kTopLow,            ///< TopLow
  kTopHigh,           ///< TopHigh
  kCRTPosMax
};

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
  /// Resize the data structure for MCNeutrino particles
  void ResizeMCNeutrino(int nNeutrinos);
  /// Resize the data strutcure for Genie primaries
  void ResizeGenie(int nPrimaries);

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
  std::vector<double> _crt_pos;          ///< CRT position

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
  std::vector<double> _ophit_peakT;           ///< Peak time of the optical hit
  std::vector<double> _ophit_width;           ///< Width of the optical hit
  std::vector<double> _ophit_area;            ///< Area of the optical hit
  std::vector<double> _ophit_amplitude;       ///< Amplitude of the optical hit
  std::vector<double> _ophit_pe;              ///< PEs of the optical hit
  std::vector<double> _ophit_opdet_x;         ///< OpDet X coordinate of the optical hit
  std::vector<double> _ophit_opdet_y;         ///< OpDet Y coordinate of the optical hit
  std::vector<double> _ophit_opdet_z;         ///< OpDet Z coordinate of the optical hit
  std::vector<int> _ophit_opdet_type;         ///< OpDet tyoe of the optical hit

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
  std::vector<Int_t>     tptype_flux;         ///< Type of parent particle leaving BNB target

  //genie information
  size_t MaxGeniePrimaries = 0;
  Int_t     genie_no_primaries;
  std::vector<Int_t>     genie_primaries_pdg;
  std::vector<Float_t>  genie_Eng;
  std::vector<Float_t>  genie_Px;
  std::vector<Float_t>  genie_Py;
  std::vector<Float_t>  genie_Pz;
  std::vector<Float_t>  genie_P;
  std::vector<Int_t>     genie_status_code;
  std::vector<Float_t>  genie_mass;
  std::vector<Int_t>     genie_trackID;
  std::vector<Int_t>     genie_ND;
  std::vector<Int_t>     genie_mother;


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
  std::string fOpHitsModuleLabel;   ///< Label for OpHit dataproduct (to be set via fcl)
  std::string fDigitModuleLabel;    ///< Label for digitizer (to be set via fcl)
  std::string fGenieGenModuleLabel; ///< Label for Genie dataproduct (to be set via fcl)

  // double fSelectedPDG;

  bool fkeepCRThits;       ///< Keep the CRT hits (to be set via fcl)
  bool fkeepCRTstrips;     ///< Keep the CRT strips (to be set via fcl)
  bool fmakeCRTtracks;     ///< Make the CRT tracks (to be set via fcl)
  bool freadCRTtracks;     ///< Keep the CRT tracks (to be set via fcl)
  bool freadOpHits;        ///< Add OpHits to output (to be set via fcl)
  bool freadTruth;         ///< Add Truth info to output (to be set via fcl)
  bool fsavePOTInfo;       ///< Add POT info to output (to be set via fcl)
  bool fcheckTransparency; ///< Checks for wire transprency (to be set via fcl)
  bool fUncompressWithPed; ///< Uncompresses the waveforms if true (to be set via fcl)
  int fWindow;
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
  fDigitModuleLabel    = p.get<std::string>("DigitModuleLabel",    "daq");
  fLArG4ModuleLabel    = p.get<std::string>("LArG4ModuleLabel",    "largeant");
  fCRTStripModuleLabel = p.get<std::string>("CRTStripModuleLabel", "crt");
  fCRTHitModuleLabel   = p.get<std::string>("CRTHitModuleLabel",   "crthit");
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel", "crttrack");
  fOpHitsModuleLabel   = p.get<std::string>("OpHitsModuleLabel");
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel", "generator");

  fkeepCRThits       = p.get<bool>("keepCRThits",true);
  fkeepCRTstrips     = p.get<bool>("keepCRTstrips",false);
  fmakeCRTtracks     = p.get<bool>("makeCRTtracks",true);
  freadCRTtracks     = p.get<bool>("readCRTtracks",true);
  freadOpHits        = p.get<bool>("readOpHits",true);
  fcheckTransparency = p.get<bool>("checkTransparency",false);
  freadTruth         = p.get<bool>("readTruth",true);
  fsavePOTInfo       = p.get<bool>("savePOTinfo",true);
  fUncompressWithPed = p.get<bool>("UncompressWithPed",false);

  fWindow            = p.get<int>("window",100);
  fKeepTaggerTypes   = p.get<std::vector<int>>("KeepTaggerTypes");

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
  }
  else {
    std::cout << "Failed to get recob::Hit data product." << std::endl;
    _nhits = 0;
  }

  if (_nhits > _max_hits) {
    std::cout << "Available hits are " << _nhits << ", which is above the maximum number allowed to store." << std::endl;
    std::cout << "Will only store " << _max_hits << "hits." << std::endl;
    _nhits = _max_hits;
  }
  for (int i = 0; i<_nhits; ++i){
    geo::WireID wireid = hitlist[i]->WireID();
    _hit_cryostat.push_back(wireid.Cryostat);
    _hit_tpc.push_back(wireid.TPC);
    _hit_plane.push_back(wireid.Plane);
    _hit_wire.push_back(wireid.Wire);
    _hit_channel.push_back(hitlist[i]->Channel());
    // peak time needs plane dependent offset correction applied.
    _hit_peakT.push_back(hitlist[i]->PeakTime());
    _hit_charge.push_back(hitlist[i]->Integral());
    _hit_ph.push_back(hitlist[i]->PeakAmplitude());
    _hit_width.push_back(hitlist[i]->RMS());
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
  for (int i = 0; i < _nstr; i += 2){
    uint32_t chan = striplist[i]->Channel();

    //    std::pair<std::string,unsigned> tagger = CRTHitRecoAlg::ChannelToTagger(chan);
    std::pair<std::string,unsigned> tagger = hitAlg.ChannelToTagger(chan);
    CRTPos ip = kNotDefined;
    if  (tagger.first=="volTaggerFaceFront_0" )    ip = kFaceFront;
    else if (tagger.first=="volTaggerFaceBack_0")  ip = kFaceBack;
    else if (tagger.first=="volTaggerSideLeft_0")  ip = kSideLeft;
    else if (tagger.first=="volTaggerSideRight_0") ip = kSideRight;
    else if (tagger.first=="volTaggerTopLow_0")    ip = kTopLow;
    else if (tagger.first=="volTaggerTopHigh_0")   ip = kTopHigh;
    else if (tagger.first=="volTaggerBot_0")       ip = kBot;

    bool keep_tagger = false;
    for (auto t : fKeepTaggerTypes) {
      if (ip == t) {
        keep_tagger = true;
      }
    }
    // std::cout << "Tagger name " << tagger.first << ", ip " << ip << ", kept? " << (keep_tagger ? "yes" : "no") << std::endl;

    if (ip != kNotDefined && keep_tagger) {

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
        TVector3 center = fAuxDetGeoCore->AuxDetChannelToPosition(2*strip, name);
        _crt_plane.push_back(ip);
        _crt_module.push_back(module);
        _crt_strip.push_back(strip);
        _crt_orient.push_back(tagger.second);
        _crt_time.push_back(ctime);
        _crt_adc.push_back(adc1 + adc2 - 127.2); // -127.2/131.9 correct for gain and 2*ped to get pe
        _crt_pos.push_back(center.Y());
        if (tagger.second == kCRTVertical) {
          _crt_pos.push_back(center.X());
        };
        ns++;
      }
    }
  }
  _nstrips = ns;

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
        if (_crt_plane[i] == kFaceFront) { // 1
          if (_crt_orient[i] == kCRTVertical && _crt_adc[i] > 500) { // < 500 hardcoded
            if (nh1x == 0 || (_crt_module[i] == plane1xm)) {
              nh1x++;
              if (_crt_adc[i] > adc1x) {
                plane1tx = _crt_time[i];
                adc1x += _crt_adc[i];
                plane1x = _crt_pos[i];
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
                plane1y=_crt_pos[i];
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
                plane2x=_crt_pos[i];
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
                plane2y=_crt_pos[i];
                plane2ym=_crt_module[i];
              }
            }
          }
        }
        for (int j=i+1;j<ns;++j) {
          float tdiff = fabs(_crt_time[i]-_crt_time[j]);
          if (tdiff<0.1) {
            iflag[j]=1;
            if (_crt_plane[j]==kFaceFront) {
              if (_crt_orient[j]==kCRTVertical && _crt_adc[j]>1000) {
                if (nh1x==0 ||  (_crt_module[j]==plane1xm)) {
                  nh1x++;
                  if (_crt_adc[j]>adc1x) {
                    plane1tx=_crt_time[j];
                    adc1x+=_crt_adc[j];
                    plane1x=_crt_pos[j];
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
                    plane1y=_crt_pos[j];
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
                    plane2x=_crt_pos[j];
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
                    plane2y=_crt_pos[j];
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
      _nchits = hitlist.size();
    }
    else {
      std::cout << "Failed to get sbn::crt::CRTHit data product." << std::endl;
      _nchits = 0;
    }

    if (_nchits > _max_chits) {
      std::cout << "Available CRT hits are " << _nchits << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << _max_chits << "CRT hits." << std::endl;
      _nchits = _max_chits;
    }
    //  std::cout << " number CRT hits " << nchits << std::endl;
    for (int i = 0; i < _nchits; ++i){
      int ip = kNotDefined;
      if  (chitlist[i]->tagger=="volTaggerFaceFront_0" )    ip = kFaceFront;
      else if (chitlist[i]->tagger=="volTaggerFaceBack_0")  ip = kFaceBack;
      else if (chitlist[i]->tagger=="volTaggerSideLeft_0")  ip = kSideLeft;
      else if (chitlist[i]->tagger=="volTaggerSideRight_0") ip = kSideRight;
      else if (chitlist[i]->tagger=="volTaggerTopLow_0")    ip = kTopLow;
      else if (chitlist[i]->tagger=="volTaggerTopHigh_0")   ip = kTopHigh;
      else if (chitlist[i]->tagger=="volTaggerBot_0")       ip = kBot;
      else {
        mf::LogWarning("HitDumper") << "Cannot identify tagger of type " << chitlist[i]->tagger << std::endl;
      }

      _chit_time[i]=chitlist[i]->ts1_ns*0.001;
      if (chitlist[i]->ts1_ns > MAX_INT) {
        _chit_time[i] = 0.001 * (chitlist[i]->ts1_ns - TIME_CORRECTION);
      }

      _chit_x.push_back(chitlist[i]->x_pos);
      _chit_y.push_back(chitlist[i]->y_pos);
      _chit_z.push_back(chitlist[i]->z_pos);
      _chit_plane.push_back(ip);
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
      for (int i = 0; i < _ncts; ++i){
        _ct_pes.push_back(ctrklist[i]->peshit);
        _ct_time.push_back(ctrklist[i]->ts1_ns*0.001);
        if (ctrklist[i]->ts1_ns > MAX_INT) {
          _ct_time.push_back(0.001 * (ctrklist[i]->ts1_ns - TIME_CORRECTION));
        }
        _ct_x1.push_back(ctrklist[i]->x1_pos);
        _ct_y1.push_back(ctrklist[i]->y1_pos);
        _ct_z1.push_back(ctrklist[i]->z1_pos);
        _ct_x2.push_back(ctrklist[i]->x2_pos);
        _ct_y2.push_back(ctrklist[i]->y2_pos);
        _ct_z2.push_back(ctrklist[i]->z2_pos);
      }
    } else {
      std::cout << "Failed to get sbn::crt::CRTTrack data product." << std::endl;
    }
  }


  //
  // Optical Hits
  //
  if (freadOpHits) {
    art::Handle< std::vector<recob::OpHit> > ophitListHandle;
    std::vector<art::Ptr<recob::OpHit> > ophitlist;
    if (evt.getByLabel(fOpHitsModuleLabel, ophitListHandle)) {
      art::fill_ptr_vector(ophitlist, ophitListHandle);
      _nophits = ophitlist.size();
    }
    else {
      std::cout << "Failed to get recob::OpHit data product." << std::endl;
      _nophits = 0;
    }

    if (_nophits > _max_ophits) {
      std::cout << "Available optical hits are " << _nophits << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << _max_ophits << " optical hits." << std::endl;
      _nophits = _max_ophits;
    }
    int counter = 0;
    for (int i = 0; i < _nophits; ++i) {
      _ophit_opch.push_back(ophitlist.at(i)->OpChannel());
      _ophit_opdet.push_back(fGeometryService->OpDetFromOpChannel(ophitlist.at(i)->OpChannel()));
      _ophit_peakT.push_back(ophitlist.at(i)->PeakTime());
      _ophit_width.push_back(ophitlist.at(i)->Width());
      _ophit_area.push_back(ophitlist.at(i)->Area());
      _ophit_amplitude.push_back(ophitlist.at(i)->Amplitude());
      _ophit_pe.push_back(ophitlist.at(i)->PE());
      auto opdet_center = fGeometryService->OpDetGeoFromOpChannel(ophitlist.at(i)->OpChannel()).GetCenter();
      _ophit_opdet_x.push_back(opdet_center.X());
      _ophit_opdet_y.push_back(opdet_center.Y());
      _ophit_opdet_z.push_back(opdet_center.Z());
      auto pd_type = _pd_map.pdType(ophitlist.at(i)->OpChannel());
      if (pd_type == "pmt_coated") {_ophit_opdet_type.push_back(kPMTCoated);}
      else if (pd_type == "pmt_uncoated") {_ophit_opdet_type.push_back(kPMTUnCoated);}
      else if (pd_type == "xarapuca_vis") {_ophit_opdet_type.push_back(kXArapucaVis);}
      else if (pd_type == "xarapuca_vuv") {_ophit_opdet_type.push_back(kXArapucaVuv);}
      else {_ophit_opdet_type.push_back(kPDNotDefined);}
      counter++;
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


  if (freadTruth){

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
    fTree->Branch("crt_pos", &_crt_pos);
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
    fTree->Branch("ophit_width", &_ophit_width);
    fTree->Branch("ophit_area", &_ophit_area);
    fTree->Branch("ophit_amplitude", &_ophit_amplitude);
    fTree->Branch("ophit_pe", &_ophit_pe);
    fTree->Branch("ophit_opdet_x", &_ophit_opdet_x);
    fTree->Branch("ophit_opdet_y", &_ophit_opdet_y);
    fTree->Branch("ophit_opdet_z", &_ophit_opdet_z);
    fTree->Branch("ophit_opdet_type", &_ophit_opdet_type);
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

  if (fsavePOTInfo) {
    _sr_tree = tfs->make<TTree>("pottree","");
    _sr_tree->Branch("run", &_sr_run, "run/I");
    _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
    _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
    _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
    _sr_tree->Branch("pot", &_sr_pot, "pot/D");
  }

}

void Hitdumper::ResetVars(){

  _run = -99999;
  _subrun = -99999;
  _event = -99999;
  _evttime = -99999;
  _t0 = -99999;
  // _adc_count = 0;

  _nhits = 0;
  _hit_cryostat.clear();
  _hit_tpc.clear();
  _hit_plane.clear();
  _hit_wire.clear();
  _hit_channel.clear();
  _hit_peakT.clear();
  _hit_charge.clear();
  _hit_ph.clear();
  _hit_width.clear();
  _hit_full_integral.clear();

  _waveform_number.clear();
  _adc_on_wire.clear();
  _time_for_waveform.clear();
  _adc_count_in_waveform.clear();

  _nstrips=0;
  _crt_plane.clear();
  _crt_module.clear();
  _crt_strip.clear();
  _crt_orient.clear();
  _crt_time.clear();
  _crt_adc.clear();
  _crt_pos.clear();

  _nctrks=0;
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

  _ncts=0;
  _ct_x1.clear();
  _ct_y1.clear();
  _ct_z1.clear();
  _ct_time.clear();
  _ct_pes.clear();
  _ct_x2.clear();
  _ct_y2.clear();
  _ct_z2.clear();

  _nchits=0;
  _chit_plane.clear();
  _chit_time.clear();
  _chit_x.clear();
  _chit_y.clear();
  _chit_z.clear();

  _nophits = 0;
  _ophit_opch.clear();
  _ophit_opdet.clear();
  _ophit_peakT.clear();
  _ophit_width.clear();
  _ophit_area.clear();
  _ophit_amplitude.clear();
  _ophit_pe.clear();
  _ophit_opdet_x.clear();
  _ophit_opdet_y.clear();
  _ophit_opdet_z.clear();
  _ophit_opdet_type.clear();

  mcevts_truth = 0;
  nuScatterCode_truth.clear();
  nuID_truth.clear();
  nuPDG_truth.clear();
  ccnc_truth.clear();
  mode_truth.clear();
  enu_truth.clear();
  Q2_truth.clear();
  W_truth.clear();
  hitnuc_truth.clear();
  nuvtxx_truth.clear();
  nuvtxy_truth.clear();
  nuvtxz_truth.clear();
  nu_dcosx_truth.clear();
  nu_dcosy_truth.clear();
  nu_dcosz_truth.clear();
  lep_mom_truth.clear();
  lep_dcosx_truth.clear();
  lep_dcosy_truth.clear();
  lep_dcosz_truth.clear();
  tpx_flux.clear();
  tpy_flux.clear();
  tpz_flux.clear();
  tptype_flux.clear();

  genie_no_primaries = 0;

}

void Hitdumper::ResizeMCNeutrino(int nNeutrinos){

  //min size is 1, to guarantee an address
  MaxMCNeutrinos = (size_t) std::max(nNeutrinos, 1);
  nuScatterCode_truth.resize(MaxMCNeutrinos);
  nuID_truth.resize(MaxMCNeutrinos);
  nuPDG_truth.resize(MaxMCNeutrinos);
  ccnc_truth.resize(MaxMCNeutrinos);
  mode_truth.resize(MaxMCNeutrinos);
  enu_truth.resize(MaxMCNeutrinos);
  Q2_truth.resize(MaxMCNeutrinos);
  W_truth.resize(MaxMCNeutrinos);
  hitnuc_truth.resize(MaxMCNeutrinos);
  nuvtxx_truth.resize(MaxMCNeutrinos);
  nuvtxy_truth.resize(MaxMCNeutrinos);
  nuvtxz_truth.resize(MaxMCNeutrinos);
  nu_dcosx_truth.resize(MaxMCNeutrinos);
  nu_dcosy_truth.resize(MaxMCNeutrinos);
  nu_dcosz_truth.resize(MaxMCNeutrinos);
  lep_mom_truth.resize(MaxMCNeutrinos);
  lep_dcosx_truth.resize(MaxMCNeutrinos);
  lep_dcosy_truth.resize(MaxMCNeutrinos);
  lep_dcosz_truth.resize(MaxMCNeutrinos);
  //Also resize the flux information here as it's a 1:1 with the MCNeutrino
  tpx_flux.resize(MaxMCNeutrinos);
  tpy_flux.resize(MaxMCNeutrinos);
  tpz_flux.resize(MaxMCNeutrinos);
  tptype_flux.resize(MaxMCNeutrinos);

  return;
} // sbnd::AnalysisTreeDataStruct::ResizeMCNeutrino()

void Hitdumper::ResizeGenie(int nPrimaries) {

  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);
  genie_primaries_pdg.resize(MaxGeniePrimaries);
  genie_Eng.resize(MaxGeniePrimaries);
  genie_Px.resize(MaxGeniePrimaries);
  genie_Py.resize(MaxGeniePrimaries);
  genie_Pz.resize(MaxGeniePrimaries);
  genie_P.resize(MaxGeniePrimaries);
  genie_status_code.resize(MaxGeniePrimaries);
  genie_mass.resize(MaxGeniePrimaries);
  genie_trackID.resize(MaxGeniePrimaries);
  genie_ND.resize(MaxGeniePrimaries);
  genie_mother.resize(MaxGeniePrimaries);

} // sbnd::AnalysisTreeDataStruct::ResizeGenie()


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
