////////////////////////////////////////////////////////////////////////
// Class:       EventAna
// Module Type: analyzer
// File:        EventAna_module.cc
// Description: Prints out information about each event.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/WhiteRabbitFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTB_dataFormat.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPCFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPC/NevisTPCTypes.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/NevisTPC/NevisTPCUtilities.hh"

#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/CRTmetric.hh"

//#include "art/Framework/Services/Optional/TFileService.h" //before art_root_io transition
#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <bitset>

namespace sbndaq {
  class EventAna;
}

/**************************************************************************************************/

class sbndaq::EventAna : public art::EDAnalyzer {

public:
  struct Config {
    //--one atom for each parameter
    fhicl::Atom<bool> include_caen {
      fhicl::Name("include_caen"),
	fhicl::Comment("look for caen 1730 fragments true/false"),
	true
	};
    fhicl::Atom<bool> caen_keepwaveforms {
      fhicl::Name("caen_keepwaveforms"),
	fhicl::Comment("put waveforms into tree true/false"),
	false
	};
    fhicl::Atom<int> Shift {
      fhicl::Name("shift_fragment_id"),
	fhicl::Comment("Number to subtract to the fragment_id"),
	0
	};
    fhicl::Atom<bool> include_wr {
      fhicl::Name("include_wr"),
	fhicl::Comment("look for wr dio fragments true/false"),
	false
	};
    fhicl::Atom<int> window_wr {
      fhicl::Name("window_wr"),
	fhicl::Comment("window around CAEN TTT to keep WR timestamps.  integer values from 0 to 10^9"),
	10000
	};
    fhicl::Atom<bool> include_berncrt {
      fhicl::Name("include_berncrt"),
	fhicl::Comment("look for bern CRT V2 fragments true/false"),
	false
	};
    fhicl::Atom<bool> crt_keepall {
      fhicl::Name("crt_keepall"),
	fhicl::Comment("put all crt fluff into tree true/false"),
	false
	};
    fhicl::Atom<bool> verbose {
      fhicl::Name("verbose"),
	fhicl::Comment("lots of text output if set to true"),
	false
	};
    fhicl::Atom<bool> include_ptb {
      fhicl::Name("include_ptb"),
      fhicl::Comment("look for ptb fragments")
    };
    fhicl::Atom<bool> include_ntb {
      fhicl::Name("include_ntb"),
      fhicl::Comment("look for ntb fragments")
    };
    fhicl::Atom<bool> include_tdc {
      fhicl::Name("include_tdc"),
	fhicl::Comment("look for spec tdc fragments (only save the nanoseconds fraction) true/false"),
	false
	};
    fhicl::Atom<bool> tdc_utc {
      fhicl::Name("tdc_utc"),
	fhicl::Comment("also save the full timestamp in utc format"),
	false
	};
  fhicl::Atom<bool> include_crtsoft {
    fhicl::Name("include_crtsoft"),
  fhicl::Comment("save the crt software trigger metric variables"),
  false
  };
  fhicl::Atom<bool> include_pmtsoft {
    fhicl::Name("include_pmtsoft"),
  fhicl::Comment("save the pmt software trigger metric variables"),
  false
  };
  fhicl::Atom<std::string> CRTInputModule {
    fhicl::Name("CRTInputModule"),
  fhicl::Comment("crt software trigger module label name"),
  "crttriggerproducer"
  };
  fhicl::Atom<std::string> PMTInputModule {
    fhicl::Name("PMTInputModule"),
  fhicl::Comment("pmt software trigger module label name"),
  "pmttriggerproducer"
  };


  }; //--configuration
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit EventAna(Parameters const & pset);
  virtual ~EventAna();

  void analyze(const art::Event& evt) override;
  void beginJob() override;
  void endJob() override;

private:

  void analyze_caen_fragment(artdaq::Fragment & frag);
  void analyze_wr_fragment_dio(artdaq::Fragment & frag);
  void analyze_bern_fragment(artdaq::Fragment & frag);
  void analyze_tdc_fragment(artdaq::Fragment & frag);
  void analyze_ntb_fragment(artdaq::Fragment & frag);
  // include ptb private classes
  void extract_triggers(artdaq::Fragment & frag);
  void reset_ptb_variables();
  void reset_caen_variables();

  //--default values
  uint32_t nChannels;//    = 16;
  uint32_t Ttt_DownSamp;// =  4;
 /* the waveforms are sampled at 500MHz sampling. The trigger timestamp is
                               * sampled 4 times slower than input channels*/

  TNtuple* nt_header;

  TH1F*    hEventCounter;
  TH1D*    hTriggerTimeTag;
  TH1F*    h_wvfm_ev0_ch0;

  TTree* events;
  int fRun;
  art::EventNumber_t fEvent;
  
  //CAEN Fragment

  std::vector<int> TTT;  // will be set to value in CAEN fragement header
  std::vector<int> TTT_ns;
  std::vector<uint64_t>  caen_frag_ts;
  std::vector<uint64_t>   ntb_frag_ts; //this is never a vector. if the ntb is in it it has to be pushing, but I don't want to fix it right now for reasons
  std::vector<uint64_t>  fTicksVec;
  std::vector< std::vector<uint16_t> >  fWvfmsVec;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch0;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch1;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch2;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch3;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch4;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch5;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch6;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch7;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch8;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch9;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch10;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch11;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch12;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch13;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch14;
  std::vector< std::vector<uint16_t> >  fWvfmsVec_ch15;
  std::vector<int> fPMT_ch0;
  std::vector<int> fPMT_ch1;
  std::vector<int> fPMT_ch2;
  std::vector<int> fPMT_ch3;
  std::vector<int> fPMT_ch4;
  std::vector<int> fPMT_ch5;
  std::vector<int> fPMT_ch6;
  std::vector<int> fPMT_ch7;
  std::vector<int> fPMT_ch8;
  std::vector<int> fPMT_ch9;
  std::vector<int> fPMT_ch10;
  std::vector<int> fPMT_ch11;
  std::vector<int> fPMT_ch12;
  std::vector<int> fPMT_ch13;
  std::vector<int> fPMT_ch14;
  std::vector<int> fPMT_ch15;
  std::vector<int> ffragID;

  //WR DIO  
  int fnstamps0;
  int fnstamps1;
  int fnstamps2;
  int fnstamps3;
  int fnstamps4;
  std::vector<float> fWR_ch0;
  std::vector<float> fWR_ch1;
  std::vector<float> fWR_ch2;
  std::vector<float> fWR_ch3;
  std::vector<float> fWR_ch4;
  bool firstEvt = true;

  int first_wr_ch0;
  int first_wr_ch1;
  int first_wr_ch2;
  int first_wr_ch3;
  int first_wr_ch4;

  //BernCRTV2 Information

  std::vector<int> flags;
  std::vector<int> lostcpu;
  std::vector<int> lostfpga;
  std::vector<uint> ts0;
  std::vector<uint> ts1;
  std::vector<int> adc0;
  std::vector<int> adc1;
  std::vector<int> adc2;
  std::vector<int> adc3;
  std::vector<int> adc4;
  std::vector<int> adc5;
  std::vector<int> adc6;
  std::vector<int> adc7;
  std::vector<int> adc8;
  std::vector<int> adc9;
  std::vector<int> adc10;
  std::vector<int> adc11;
  std::vector<int> adc12;
  std::vector<int> adc13;
  std::vector<int> adc14;
  std::vector<int> adc15;
  std::vector<int> adc16;
  std::vector<int> adc17;
  std::vector<int> adc18;
  std::vector<int> adc19;
  std::vector<int> adc20;
  std::vector<int> adc21;
  std::vector<int> adc22;
  std::vector<int> adc23;
  std::vector<int> adc24;
  std::vector<int> adc25;
  std::vector<int> adc26;
  std::vector<int> adc27;
  std::vector<int> adc28;
  std::vector<int> adc29;
  std::vector<int> adc30;
  std::vector<int> adc31;
  std::vector<int> coinc;

  std::vector<std::pair<int, int>> max_adc;
  std::vector<int>                 max_chan;

  std::vector<int>  feb_hit_number          ; //hit counter for individual FEB, including hits lost in FEB or fragment generator
  std::vector<uint64_t>  timestamp               ; //absolute timestamp
  std::vector<uint>  last_accepted_timestamp ; //timestamp of previous accepted hit
  std::vector<int>  lost_hits               ; //number of lost hits from the previous one

  // CRT metadata
  std::vector<int>  mac5; //last 8 bits of FEB mac5 address
  std::vector<uint>  run_start_time;
  std::vector<uint>  this_poll_start;
  std::vector<uint>  this_poll_end;
  std::vector<uint>  last_poll_start;
  std::vector<uint>  last_poll_end;
  std::vector<int>   system_clock_deviation;
  std::vector<int>  feb_hits_in_poll;
  std::vector<int>  feb_hits_in_fragment;
  
  //WR spectdc data
  std::vector<uint64_t> ftdc_ch0; //fractional part of the timestamp
  std::vector<uint64_t> ftdc_ch1;
  std::vector<uint64_t> ftdc_ch2;
  std::vector<uint64_t> ftdc_ch3;
  std::vector<uint64_t> ftdc_ch4;
  std::vector<uint64_t> ftdc_ch0_utc; //timestamp in UTC format
  std::vector<uint64_t> ftdc_ch1_utc;
  std::vector<uint64_t> ftdc_ch2_utc;
  std::vector<uint64_t> ftdc_ch3_utc;
  std::vector<uint64_t> ftdc_ch4_utc;

  //information from fragment header
  std::vector<int>  sequence_id;

  bool finclude_tdc;
  bool ftdc_utc;
  bool finclude_caen;
  bool fcaen_keepwaveforms;
  int fShift;
  bool finclude_wr;
  int fWindow;
  bool finclude_berncrt;
  bool fcrt_keepall;
  bool fverbose;
  bool finclude_ptb;
  bool finclude_ntb;
  bool finclude_crtsoft;
  bool finclude_pmtsoft;
  std::string fcrtSoftTriggerModuleLabel;
  std::string fpmtSoftTriggerModuleLabel;


  // including ptb information on the tree
  bool unknown_or_error_word; // flag to indicate the event has
  int ts_word_count;
  int hlt_word_count;
  std::vector<uint64_t> ptb_frag_ts;
  std::vector<uint64_t> llt_trigger;
  std::vector<uint64_t> llt_ts;
  std::vector<uint64_t> hlt_trigger;
  std::vector<uint64_t> hlt_ts;
  std::vector<uint64_t> hlt_gateCount;
  std::vector<uint16_t> crt_status;
  std::vector<uint16_t> beam_status;
  std::vector<uint16_t> mtca_status;
  std::vector<uint16_t> nim_status;
  std::vector<uint32_t> auxpds_status;
  std::vector<uint64_t> chan_stat_ts;

  // PMT software trigger variables
  bool   _pmtSoftTrigger_foundBeamTrigger;   /// Whether the beam spill was found or not
  int    _pmtSoftTrigger_tts;                /// Trigger Time Stamp (TTS), ns (relative to start of beam spill)
  double _pmtSoftTrigger_promptPE;           /// Total photoelectron count 100 ns after the TTS
  double _pmtSoftTrigger_prelimPE;           /// Total photoelectron count before the TTS, during the beam spill
  int    _pmtSoftTrigger_nAboveThreshold;    /// number of individual PMTs above ADC threshold (fcl) during the beam spill
  // std::vector<sbnd::trigger::pmtInfo> _pmtSoftTrigger_pmtInfoVec; /// vector of PMT information

  // CRT software trigger variables
  int    _crtSoftTrigger_hitsperplane[7];       ///< Number of (very low level) CRT hits per plane


}; //--class EventAna


sbndaq::EventAna::EventAna(EventAna::Parameters const& pset): art::EDAnalyzer(pset)
{
  finclude_caen = pset().include_caen();
  fcaen_keepwaveforms = pset().caen_keepwaveforms();
  fShift = pset().Shift();
  finclude_wr = pset().include_wr();
  finclude_tdc = pset().include_tdc();
  finclude_ntb= pset().include_ntb();
  ftdc_utc = pset().tdc_utc();
  fWindow = pset().window_wr();
  if (fWindow<0 || fWindow>1000000000) {
    fWindow=1000000000;
    std::cout << "Bad value for fcl parameter window_wr=" << fWindow << "  setting to default value 1000000000 nsec = 1 sec"
	    << std::endl;
  }
  fverbose = pset().verbose();
  finclude_berncrt = pset().include_berncrt();
  fcrt_keepall = pset().crt_keepall();
  finclude_ptb = pset().include_ptb();

  finclude_crtsoft = pset().include_crtsoft();
  finclude_pmtsoft = pset().include_pmtsoft();
  fcrtSoftTriggerModuleLabel = pset().CRTInputModule();
  fpmtSoftTriggerModuleLabel = pset().PMTInputModule();
}

void sbndaq::EventAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  nt_header       = tfs->make<TNtuple>("nt_header","Multi Header Ntuple","art_ev:caen_ev:caenv_ev_tts");
  /************************************************************************************************/
  hEventCounter   = tfs->make<TH1F>("hEventCounter","Event Counter Histogram",10000,0,10000);
  hTriggerTimeTag = tfs->make<TH1D>("hTriggerTimeTag","Trigger Time Tag Histogram",12500,0,125000000);
  h_wvfm_ev0_ch0  = tfs->make<TH1F>("h_wvfm_ev0_ch0","Waveform",2000,0,2000);
  /************************************************************************************************/
  //--make tree to store the channel waveform info:
  events = tfs->make<TTree>("events","waveform tree");
  events->Branch("fRun",&fRun,"fRun/I");
  events->Branch("fEvent",&fEvent,"fEvent/I");
  if (finclude_caen) {
    events->Branch("TTT_ns",&TTT_ns);
    events->Branch("caen_frag_ts",&caen_frag_ts);
    if (fcaen_keepwaveforms) {
      events->Branch("fTicksVec",&fTicksVec);
      events->Branch("fWvfmsVec",&fWvfmsVec);
      events->Branch("fWvfmsVec_ch0",&fWvfmsVec_ch0);
      events->Branch("fWvfmsVec_ch1",&fWvfmsVec_ch1);
      events->Branch("fWvfmsVec_ch2",&fWvfmsVec_ch2);
      events->Branch("fWvfmsVec_ch3",&fWvfmsVec_ch3);
      events->Branch("fWvfmsVec_ch4",&fWvfmsVec_ch4);
      events->Branch("fWvfmsVec_ch5",&fWvfmsVec_ch5);
      events->Branch("fWvfmsVec_ch6",&fWvfmsVec_ch6);
      events->Branch("fWvfmsVec_ch7",&fWvfmsVec_ch7);
      events->Branch("fWvfmsVec_ch8",&fWvfmsVec_ch8);
      events->Branch("fWvfmsVec_ch9",&fWvfmsVec_ch9);
      events->Branch("fWvfmsVec_ch10",&fWvfmsVec_ch10);
      events->Branch("fWvfmsVec_ch11",&fWvfmsVec_ch11);
      events->Branch("fWvfmsVec_ch12",&fWvfmsVec_ch12);
      events->Branch("fWvfmsVec_ch13",&fWvfmsVec_ch13);
      events->Branch("fWvfmsVec_ch14",&fWvfmsVec_ch14);
      events->Branch("fWvfmsVec_ch15",&fWvfmsVec_ch15);
    }
    events->Branch("fPMT_ch0",&fPMT_ch0);
    events->Branch("fPMT_ch1",&fPMT_ch1);
    events->Branch("fPMT_ch2",&fPMT_ch2);
    events->Branch("fPMT_ch3",&fPMT_ch3);
    events->Branch("fPMT_ch4",&fPMT_ch4);
    events->Branch("fPMT_ch5",&fPMT_ch5);
    events->Branch("fPMT_ch6",&fPMT_ch6);
    events->Branch("fPMT_ch7",&fPMT_ch7);
    events->Branch("fPMT_ch8",&fPMT_ch8);
    events->Branch("fPMT_ch9",&fPMT_ch9);
    events->Branch("fPMT_ch10",&fPMT_ch10);
    events->Branch("fPMT_ch11",&fPMT_ch11);
    events->Branch("fPMT_ch12",&fPMT_ch12);
    events->Branch("fPMT_ch13",&fPMT_ch13);
    events->Branch("fPMT_ch14",&fPMT_ch14);
    events->Branch("fPMT_ch15",&fPMT_ch15);
    events->Branch("ffragID",&ffragID);
  }
  if (finclude_wr) {
    events->Branch("fWR_ch0",&fWR_ch0);
    events->Branch("fWR_ch1",&fWR_ch1);
    events->Branch("fWR_ch2",&fWR_ch2);
    events->Branch("fWR_ch3",&fWR_ch3);
    events->Branch("fWR_ch4",&fWR_ch4);
  }
  if (finclude_tdc) {
    events->Branch("ftdc_ch0",&ftdc_ch0);
    events->Branch("ftdc_ch1",&ftdc_ch1);
    events->Branch("ftdc_ch2",&ftdc_ch2);
    events->Branch("ftdc_ch3",&ftdc_ch3);
    events->Branch("ftdc_ch4",&ftdc_ch4);
    if(ftdc_utc){
      events->Branch("ftdc_ch0_utc",&ftdc_ch0_utc);
      events->Branch("ftdc_ch1_utc",&ftdc_ch1_utc);
      events->Branch("ftdc_ch2_utc",&ftdc_ch2_utc);
      events->Branch("ftdc_ch3_utc",&ftdc_ch3_utc);
      events->Branch("ftdc_ch4_utc",&ftdc_ch4_utc);
    }
  }
  if (finclude_ntb){
    events->Branch("ntb_frag_ts",&ntb_frag_ts);
  }
  if (finclude_berncrt) {
    events->Branch("flags",         &flags);
    events->Branch("coinc",         &coinc);
    events->Branch("mac5",          &mac5);
    events->Branch("ts0",           &ts0);
    events->Branch("ts1",           &ts1);
    events->Branch("adc0",           &adc0);
    events->Branch("adc1",           &adc1);
    events->Branch("adc2",           &adc2);
    events->Branch("adc3",           &adc3);
    events->Branch("adc4",           &adc4);
    events->Branch("adc5",           &adc5);
    events->Branch("adc6",           &adc6);
    events->Branch("adc7",           &adc7);
    events->Branch("adc8",           &adc8);
    events->Branch("adc9",           &adc9);
    events->Branch("adc10",           &adc10);
    events->Branch("adc11",           &adc11);
    events->Branch("adc12",           &adc12);
    events->Branch("adc13",           &adc13);
    events->Branch("adc14",           &adc14);
    events->Branch("adc15",           &adc15);
    events->Branch("adc16",           &adc16);
    events->Branch("adc17",           &adc17);
    events->Branch("adc18",           &adc18);
    events->Branch("adc19",           &adc19);
    events->Branch("adc20",           &adc20);
    events->Branch("adc21",           &adc21);
    events->Branch("adc22",           &adc22);
    events->Branch("adc23",           &adc23);
    events->Branch("adc24",           &adc24);
    events->Branch("adc25",           &adc25);
    events->Branch("adc26",           &adc26);
    events->Branch("adc27",           &adc27);
    events->Branch("adc28",           &adc28);
    events->Branch("adc29",           &adc29);
    events->Branch("adc30",           &adc30);
    events->Branch("adc31",           &adc31);
    events->Branch("max_adc",         &max_adc);
    events->Branch("max_chan",        &max_chan);
    if (fcrt_keepall) {
      events->Branch("timestamp",     &timestamp);
      events->Branch("lostcpu",       &lostcpu);
      events->Branch("lostfpga",      &lostfpga);
      events->Branch("feb_hit_number",&feb_hit_number);
      events->Branch("last_accepted_timestamp",&last_accepted_timestamp);
      events->Branch("lost_hits",     &lost_hits);
      events->Branch("run_start_time",            &run_start_time);
      events->Branch("this_poll_start",           &this_poll_start);
      events->Branch("this_poll_end",             &this_poll_end);
      events->Branch("last_poll_start",           &last_poll_start);
      events->Branch("last_poll_end",             &last_poll_end);
      events->Branch("system_clock_deviation",    &system_clock_deviation);
      events->Branch("feb_hits_in_poll",          &feb_hits_in_poll);
      events->Branch("feb_hits_in_fragment",      &feb_hits_in_fragment);
      events->Branch("sequence_id",               &sequence_id);
    }
  }
  // include ptb branches
  if(finclude_ptb){
    events->Branch("unknown_or_error_word", &unknown_or_error_word);
    events->Branch("ts_word_count", &ts_word_count);
    events->Branch("hlt_word_count", &hlt_word_count);
    events->Branch("ptb_frag_ts", &ptb_frag_ts);
    // Trigger words and TS
    events->Branch("hlt_trigger", &hlt_trigger);
    events->Branch("hlt_ts",      &hlt_ts);
    events->Branch("hlt_gateCount", &hlt_gateCount);
    events->Branch("llt_trigger", &llt_trigger);
    events->Branch("llt_ts",      &llt_ts);
    // Channel status words & TS
    events->Branch("beam_status",   &beam_status);
    events->Branch("crt_status",    &crt_status);
    events->Branch("mtca_status",   &mtca_status);
    events->Branch("nim_status",    &nim_status);
    events->Branch("auxpds_status", &auxpds_status);
    events->Branch("chan_stat_ts",  &chan_stat_ts);
  }

  if(finclude_crtsoft){
    events->Branch("crtSoftTrigger_hitsperplane", &_crtSoftTrigger_hitsperplane,"crtSoftTrigger_hitsperplane[7]/I");
  }
  if(finclude_pmtsoft){
    events->Branch("pmtSoftTrigger_foundBeamTrigger", &_pmtSoftTrigger_foundBeamTrigger);
    events->Branch("pmtSoftTrigger_tts", &_pmtSoftTrigger_tts);
    events->Branch("pmtSoftTrigger_promptPE", &_pmtSoftTrigger_promptPE);
    events->Branch("pmtSoftTrigger_prelimPE", &_pmtSoftTrigger_prelimPE);
    events->Branch("pmtSoftTrigger_nAboveThreshold", &_pmtSoftTrigger_nAboveThreshold);
  }

}

void sbndaq::EventAna::endJob()
{
  if (fverbose)  std::cout << "Ending EventAna...\n";
}


sbndaq::EventAna::~EventAna()
{
}


void sbndaq::EventAna::analyze(const art::Event& evt)
{
  fRun = evt.run();
  fEvent = evt.event();
  if (fverbose)
    std::cout << "\n\nRun " << fRun << " event " << fEvent << std::endl;

  /************************************************************************************************/
  // need to clear tree variables at the beginning of the event

  //WR DIO
  fnstamps0=0;   fnstamps1=0;   fnstamps2=0;   fnstamps3=0;   fnstamps4=0;
  fWR_ch0.clear();   fWR_ch1.clear();   fWR_ch2.clear();   fWR_ch3.clear();   fWR_ch4.clear();
  first_wr_ch0=0;  first_wr_ch1=0;  first_wr_ch2=0;  first_wr_ch3=0;  first_wr_ch4=0;

  //WR SPEC TDC
  ftdc_ch0.clear();   ftdc_ch1.clear();   ftdc_ch2.clear();   ftdc_ch3.clear();   ftdc_ch4.clear();
  ftdc_ch0_utc.clear();   ftdc_ch1_utc.clear();   ftdc_ch2_utc.clear();   ftdc_ch3_utc.clear();   ftdc_ch4_utc.clear();
  /************************************************************************************************/
  //BERN CRT 
  mac5.clear();    flags.clear();   lostcpu.clear();   lostfpga.clear();   ts0.clear();     ts1.clear();
  adc0.clear();    adc1.clear();    adc2.clear();      adc3.clear();       adc4.clear();    adc5.clear();     adc6.clear();
  adc7.clear();    adc8.clear();    adc9.clear();      adc10.clear();      adc11.clear();   adc12.clear();    adc13.clear();
  adc14.clear();   adc15.clear();   adc16.clear();     adc17.clear();      adc18.clear();   adc19.clear();    adc20.clear();
  adc21.clear();   adc22.clear();   adc23.clear();     adc24.clear();      adc25.clear();   adc26.clear();    adc27.clear();
  adc28.clear();   adc29.clear();   adc30.clear();     adc31.clear();      coinc.clear();   max_adc.clear();  max_chan.clear();


  feb_hit_number.clear()       ;   timestamp.clear()      ;    last_accepted_timestamp.clear();
  lost_hits.clear()            ;   run_start_time.clear() ;    this_poll_start.clear()        ;   this_poll_end.clear();
  last_poll_start.clear()      ;   last_poll_end.clear()  ;    system_clock_deviation.clear();    feb_hits_in_poll.clear();
  feb_hits_in_fragment.clear() ;   sequence_id.clear();

  /************************************************************************************************/
  ntb_frag_ts.clear();

  // Reset PTB variables
  reset_ptb_variables();

  // Reset 1730 variables
  reset_caen_variables();

  /************************************************************************************************/

  for (int i=0; i<7; i++){
    _crtSoftTrigger_hitsperplane[i] = 0;
  }

  _pmtSoftTrigger_foundBeamTrigger = false;
  _pmtSoftTrigger_tts = 0;
  _pmtSoftTrigger_promptPE = 0;
  _pmtSoftTrigger_prelimPE = 0;
  _pmtSoftTrigger_nAboveThreshold = 0;
 /************************************************************************************************/


  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;

  fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();

  /************************************************************************************************/
  if (finclude_caen) { 
    for (auto handle : fragmentHandles) {
      if (!handle.isValid() || handle->size() == 0) continue;

      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
	//Container fragment
	for (auto cont : *handle) {
	  artdaq::ContainerFragment contf(cont);
	  if (contf.fragment_type()==sbndaq::detail::FragmentType::CAENV1730) {
	    if (fverbose) 	  std::cout << "    Found " << contf.block_count() << " CAEN Fragments in container " << std::endl;
 	    fWvfmsVec.resize(16*contf.block_count());
	    for (size_t ii = 0; ii < contf.block_count(); ++ii)
	      analyze_caen_fragment(*contf[ii].get());
	  }
	}
      }
      else {
	//normal fragment
	if (handle->front().type()==sbndaq::detail::FragmentType::CAENV1730) {
	  if (fverbose)	std::cout << "   found normal caen fragments " << handle->size() << std::endl;
	  fWvfmsVec.resize(16*handle->size());
 	  for (auto frag : *handle)
	    analyze_caen_fragment(frag);
	}
      }
    } // loop over frag handles
  }  // if include caen

  /************************************************************************************************/
  if (finclude_wr) {
    for (auto handle : fragmentHandles) {
      if (!handle.isValid() || handle->size() == 0) continue;
      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
	//Container fragment
	for (auto cont : *handle) {
	  artdaq::ContainerFragment contf(cont);
	  if (contf.fragment_type()==sbndaq::detail::FragmentType::WhiteRabbit) {
	    if (fverbose) 	  std::cout << "    Found " << contf.block_count() << " WR Fragments in container " << std::endl;
	    for (size_t ii = 0; ii < contf.block_count(); ++ii)
	      analyze_wr_fragment_dio(*contf[ii].get());
	  }
	}
      }
      else {
	//normal fragment
	if (handle->front().type()==sbndaq::detail::FragmentType::WhiteRabbit) {
	  for (auto frag : *handle)
	    analyze_wr_fragment_dio(frag);
	}
      }
    } // loop over frag handles

    if (fverbose){ std::cout << " WR ch 0 " << fnstamps0 << " WR ch 1 " << fnstamps1 << " WR ch 2 " << fnstamps2 << " WR ch 3 " <<
		    fnstamps3 << " WR ch 4 " << fnstamps4 << std::endl;
    }
  } // if (include_wr)


  /************************************************************************************************/
  if (finclude_berncrt){

    std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;

    fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();

    for (auto handle : fragmentHandles) {
      if (!handle.isValid() || handle->size() == 0)
	continue;

      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
        //Container fragment
        for (auto cont : *handle) {
          artdaq::ContainerFragment contf(cont);
          if (contf.fragment_type() != sbndaq::detail::FragmentType::BERNCRTV2)
            continue;
          for (size_t ii = 0; ii < contf.block_count(); ++ii)
            analyze_bern_fragment(*contf[ii].get());
        }
      }
      else {
        //normal fragment
        if (handle->front().type() != sbndaq::detail::FragmentType::BERNCRTV2) continue;
        for (auto frag : *handle)
          analyze_bern_fragment(frag);
      }
    }  //  loop over frag handles
  }// if include_berncrt

  /************************************************************************************************/
  if (finclude_tdc) {
    for (auto handle : fragmentHandles) {
      if (!handle.isValid() || handle->size() == 0) continue;
      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
	//Container fragment
	for (auto cont : *handle) {
	  artdaq::ContainerFragment contf(cont);
	  if (contf.fragment_type()==sbndaq::detail::FragmentType::TDCTIMESTAMP) {
	    if (fverbose) 	  std::cout << "    Found " << contf.block_count() << " TDC Timestamp Fragments in container " << std::endl;
	    for (size_t ii = 0; ii < contf.block_count(); ++ii)
	      analyze_tdc_fragment(*contf[ii].get());
	  }
	}
      }
      else {
	//normal fragment
	if (handle->front().type()==sbndaq::detail::FragmentType::TDCTIMESTAMP) {
	  for (auto frag : *handle)
	    analyze_tdc_fragment(frag);
	}
      }
    } // loop over frag handles

  } // if (include_tdc)

   /************************************************************************************************/
   // Save PTB data in tree
   if(finclude_ptb) {
     std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;

         fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();

     for (auto handle : fragmentHandles) {
       if (!handle.isValid() || handle->size() == 0) continue;

       if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
         //Container fragment
         for (auto cont : *handle) {
           artdaq::ContainerFragment contf(cont);
           if (contf.fragment_type() != sbndaq::detail::FragmentType::PTB) continue;
           for (size_t ii = 0; ii < contf.block_count(); ++ii) {
             extract_triggers(*contf[ii].get());
           }
         }
       }
       else {
         //normal fragment
         if (handle->front().type() != sbndaq::detail::FragmentType::PTB) continue;
         for (auto frag : *handle) {
           extract_triggers(frag);
         }
       }
     }
   } // if includes ptb

   /************************************************************************************************/
   // Save Software Trigger Metrics

   if (finclude_crtsoft){
     art::Handle<sbndaq::CRTmetric> crtSoftTriggerHandle;
     if (evt.getByLabel(fcrtSoftTriggerModuleLabel, crtSoftTriggerHandle)){
      const sbndaq::CRTmetric &crtSoftTriggerMetrics = (*crtSoftTriggerHandle);
      for (int i=0; i<7; i++){
	       _crtSoftTrigger_hitsperplane[i] = crtSoftTriggerMetrics.hitsperplane[i];
      }
    }
    else{
      std::cout << "Failed to get sbndaq::crtMetric data product" << std::endl;
    }

  }//if include crt soft trig

  if (finclude_pmtsoft){
    art::Handle<sbnd::trigger::pmtSoftwareTrigger> pmtSoftTriggerHandle;
    if (evt.getByLabel(fpmtSoftTriggerModuleLabel, pmtSoftTriggerHandle)){
      const sbnd::trigger::pmtSoftwareTrigger &pmtSoftTriggerMetrics = (*pmtSoftTriggerHandle);
      _pmtSoftTrigger_foundBeamTrigger = pmtSoftTriggerMetrics.foundBeamTrigger;
      _pmtSoftTrigger_tts = pmtSoftTriggerMetrics.trig_ts;
      _pmtSoftTrigger_promptPE = pmtSoftTriggerMetrics.promptPE;
      _pmtSoftTrigger_prelimPE = pmtSoftTriggerMetrics.prelimPE;
      _pmtSoftTrigger_nAboveThreshold = pmtSoftTriggerMetrics.nAboveThreshold;
      // _pmtSoftTrigger_pmtInfoVec = pmtsofttriggerlist[0]->pmtInfoVec;
    }
    else{
      std::cout << "Failed to get sbnd::trigger::pmtSoftwareTrigger data product" << std::endl;
    }
  }//if include pmt soft trig


   // Fill the tree
   events->Fill();

}



void sbndaq::EventAna::analyze_wr_fragment_dio(artdaq::Fragment & frag)  {



  if (fverbose)       std::cout <<  "     timestamp is  " << frag.timestamp() << std::endl;
  if (fverbose)       std::cout <<  "     seq ID is " << frag.sequenceID() << std::endl;


  const WhiteRabbitEvent *event_ptr = reinterpret_cast<WhiteRabbitEvent const*>(frag.dataBeginBytes());
  timespec sysTime=event_ptr->systemTime;
  if (fverbose)      std::cout << "systime  " << sysTime.tv_sec << " " << sysTime.tv_nsec << std::endl;
  WhiteRabbitData fragdata=event_ptr->data;

  // each WR fragment has data from only one channel. The fragments are not always in time order
  if (fverbose)       std::cout << "WR: command " << fragdata.command << std::endl;
  if (fverbose)       std::cout << "WR: channel " << fragdata.channel << std::endl;
  if (fverbose)       std::cout << "WR: value " << fragdata.value << std::endl;
  if (fverbose)       std::cout << "WR: flags " << fragdata.flags << std::endl;
  if (fverbose)       std::cout << "WR: number of time stamps in this fragment " << fragdata.nstamp << std::endl;
  if (fragdata.channel==0)  fnstamps0+=fragdata.nstamp;
  if (fragdata.channel==1)	fnstamps1+=fragdata.nstamp;
  if (fragdata.channel==2)	fnstamps2+=fragdata.nstamp;
  if (fragdata.channel==3)  fnstamps3+=fragdata.nstamp;
  if (fragdata.channel==4)	fnstamps4+=fragdata.nstamp;

  std::cout << "WR timestamp count " << fnstamps0 << " + " << fnstamps1 << " + " << fnstamps2 << " + " << fnstamps3 << " + " << fnstamps4 << std::endl;

  for (int i=0;i<(int)fragdata.nstamp;++i) {

    int caenTTT = TTT_ns[0]; //check the TTT of the first caen fragment, assuming that a caen fragment is already filled
    int diff = 0;
    uint this_time_sec = fragdata.timeStamp[i].tv_sec;
    uint this_time_ns = fragdata.timeStamp[i].tv_nsec;
    // if (caenTTT_ns>this_time) {
    //    std::cout << this_time_sec << std::endl;
    diff = caenTTT-this_time_ns;
    if (diff>500000000) diff = 1000000000-diff;
    else if (diff<-500000000) diff = 1000000000+diff;
    // }
    // else {

    //   diff = this_time-caenTTT_ns;
    //   if (diff>500000000) diff = 1000000000-diff;
    // }
    //	if (diff>500000000) std::cout<< "diff " << diff << "this_time_ns " << this_time_ns << " caenTTT " << caenTTT << std::endl;



    if (fabs(diff)<fWindow && fragdata.channel==1) 	 {
      if (first_wr_ch1==0) first_wr_ch1=fragdata.timeStamp[i].tv_sec;
      int secdiff=0;
      if (fragdata.timeStamp[i].tv_sec>first_wr_ch1) secdiff=this_time_sec-first_wr_ch1;
      fWR_ch1.emplace_back((int)this_time_ns+1e9*secdiff);
      if (fverbose){ 	  std::cout << " Event " << fEvent << " PMT" <<
			    " Timestamp " << i << "  : " << std::setw(16) << fragdata.timeStamp[i].tv_sec <<
			    " " << std::setw(9) << fragdata.timeStamp[i].tv_nsec <<
			    " TTT " << std::setw(9) << caenTTT <<
			    " TTT diff  " << std::setw(9)  << diff << std::endl;
      }
    }
    if (fabs(diff)< 50000000 && fragdata.channel==2) 	{
      if (first_wr_ch2==0) first_wr_ch2=fragdata.timeStamp[i].tv_sec;
      int secdiff=0;
      if (fragdata.timeStamp[i].tv_sec>first_wr_ch2) secdiff=fragdata.timeStamp[i].tv_sec-first_wr_ch2;
      fWR_ch2.emplace_back((int)fragdata.timeStamp[i].tv_nsec+1e9*secdiff);
      if (fverbose){ 	  std::cout << " Event " << fEvent << " RWM" <<
			    " Timestamp " << i << "  : " << std::setw(16) << fragdata.timeStamp[i].tv_sec <<
			    " " << std::setw(9) << fragdata.timeStamp[i].tv_nsec <<
			    " TTT " << std::setw(9) << caenTTT <<
			    " TTT diff  " << std::setw(9)  << diff << std::endl;
      }
    }
    if (fragdata.channel==0 ) 	{
      if (first_wr_ch0==0) first_wr_ch0=fragdata.timeStamp[i].tv_sec;
      int secdiff=0;
      if (fragdata.timeStamp[i].tv_sec>first_wr_ch0) secdiff=fragdata.timeStamp[i].tv_sec-first_wr_ch0;
      fWR_ch0.emplace_back((int)fragdata.timeStamp[i].tv_nsec+1e9*secdiff);
      if (fverbose){ 	  std::cout << " Event " << fEvent << " PPS" <<
			    " Timestamp " << i << "  : " << std::setw(16) << fragdata.timeStamp[i].tv_sec <<
			    " " << std::setw(9) << fragdata.timeStamp[i].tv_nsec <<
			    " TTT " << std::setw(9) << caenTTT <<
			    " TTT diff  " << std::setw(9)  << diff << std::endl;
      }
    }
    // if (diff<(uint)fWindow && fragdata.channel==3 )
    if ( fabs(diff)<5000000 && fragdata.channel==3) 	 {
      if (first_wr_ch3==0) first_wr_ch3=fragdata.timeStamp[i].tv_sec;
      int secdiff=0;
      if (fragdata.timeStamp[i].tv_sec>first_wr_ch3) secdiff=fragdata.timeStamp[i].tv_sec-first_wr_ch3;
      fWR_ch3.emplace_back((int)fragdata.timeStamp[i].tv_nsec+1e9*secdiff);
      if (fverbose){ 	  std::cout << " Event " << fEvent << " TRIG" <<
			    " Timestamp " << i << "  : " << std::setw(16) << fragdata.timeStamp[i].tv_sec <<
			    " " << std::setw(10) << fragdata.timeStamp[i].tv_nsec <<
			    " TTT " << std::setw(10) << caenTTT <<
			    " TTT diff  " << std::setw(10)  << diff << std::endl;
      }
    }
  }
  std::cout << " ----------------------- " << std::endl;

}

void sbndaq::EventAna::analyze_caen_fragment(artdaq::Fragment & frag)  {

  if (fverbose) std::cout <<  "     timestamp is  " << frag.timestamp() << std::endl;
  if (fverbose) std::cout <<  "     seq ID is " << frag.sequenceID() << std::endl;

  caen_frag_ts.push_back(frag.timestamp());

  CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();
  CAENV1730Event const* event_ptr = bb.Event();
  CAENV1730EventHeader header = event_ptr->Header;

  int fragId = static_cast<int>(frag.fragmentID());
  ffragID.push_back(fragId);
  
  if (fverbose)      std::cout << "\tFrom CAEN header, event counter is "  << header.eventCounter   << "\n";
  if (fverbose)      std::cout << "\tFrom CAEN header, triggerTimeTag is " << header.triggerTimeTag << "\n";
  if (fverbose)       std::cout << "\tFrom CAEN header, board id is "       << header.boardID       << "\n";
  if (fverbose)       std::cout << "\tFrom CAEN fragment, fragment id is "  << fragId << "\n";
  if (fverbose)       std::cout << "\tFragment counter for this event "  << fShift << "\n";

  uint32_t t0 = header.triggerTimeTag;
//  TTT = (int)t0;
//  TTT_ns = t0*8;
  TTT.push_back((int)t0);
  TTT_ns.push_back((int)t0*8);

  if (fverbose)       std::cout << "\n\tTriggerTimeTag in ns is " << (int)t0*8 << "\n";  // 500 MHz is 2 ns per tick
  hEventCounter->Fill(header.eventCounter);
  hTriggerTimeTag->Fill((int)t0);
  nt_header->Fill(fEvent,header.eventCounter,t0);
  nChannels = md->nChannels;
  if (fverbose){      std::cout << "\tNumber of channels: " << nChannels << "\n";
  }

  //--get the number of 32-bit words (quad_bytes) from the header
  uint32_t ev_size_quad_bytes = header.eventSize;
  if (fverbose){       std::cout << "Event size in quad bytes is: " << ev_size_quad_bytes << "\n";
  }
  uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
  uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
  uint32_t wfm_length = data_size_double_bytes/nChannels;
  //--note, needs to take into account channel mask
  if (fverbose) std::cout << "Channel waveform length = " << wfm_length << "\n";

  //--store the tick value for each acquisition
  fTicksVec.resize(wfm_length);

  const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes()
						      	 + sizeof(CAENV1730EventHeader));
  const uint16_t* value_ptr =  data_begin;
  uint16_t value = 0;
  size_t ch_offset = 0;
  
  fShift = 0; //temporary fix for having multiple caen fragments...
  std::vector<uint16_t> temp;
    
  for (size_t i_ch=0; i_ch<nChannels; ++i_ch){

    fWvfmsVec[i_ch+nChannels*fShift].resize(wfm_length);
    ch_offset = (size_t)(i_ch * wfm_length);
    
    temp.resize(wfm_length);

    //--loop over waveform samples
    for(size_t i_t=0; i_t<wfm_length; ++i_t){
      value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
      value = *(value_ptr);

      temp[i_t] = value;

      fTicksVec[i_t] = t0*Ttt_DownSamp + i_t;   /*timestamps, event level, TODO: Ttt_DownSamp is not initiated*/
      fWvfmsVec[i_ch+nChannels*fShift][i_t] = value;
      if (i_ch == 0 && firstEvt) h_wvfm_ev0_ch0->SetBinContent(i_t,value);
     
    } //--end waveforms loop
  
    if(i_ch == 0) fWvfmsVec_ch0.push_back(temp);
    if(i_ch == 1) fWvfmsVec_ch1.push_back(temp);
    if(i_ch == 2) fWvfmsVec_ch2.push_back(temp);
    if(i_ch == 3) fWvfmsVec_ch3.push_back(temp);
    if(i_ch == 4) fWvfmsVec_ch4.push_back(temp);
    if(i_ch == 5) fWvfmsVec_ch5.push_back(temp);
    if(i_ch == 6) fWvfmsVec_ch6.push_back(temp);
    if(i_ch == 7) fWvfmsVec_ch7.push_back(temp);
    if(i_ch == 8) fWvfmsVec_ch8.push_back(temp);
    if(i_ch == 9) fWvfmsVec_ch9.push_back(temp);
    if(i_ch == 10) fWvfmsVec_ch10.push_back(temp);
    if(i_ch == 11) fWvfmsVec_ch11.push_back(temp);
    if(i_ch == 12) fWvfmsVec_ch12.push_back(temp);
    if(i_ch == 13) fWvfmsVec_ch13.push_back(temp);
    if(i_ch == 14) fWvfmsVec_ch14.push_back(temp);
    if(i_ch == 15) fWvfmsVec_ch15.push_back(temp);
 
    temp.clear();
    firstEvt = false;
  } //--end channels loop 

  // threshold values and fragID are hardcoded, should be fcl params instead.
  int threshold[]= { 9000, 9000, 9000, 9000,9000, 9000, 9000, 9000,9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000};

  int thisind = ffragID.size();
  thisind--;
  if (ffragID.at(thisind)==9) {
  // find leading edges in waveforms
  int toggle=0;
  int i_ch =0 ;
  auto this_wf = fWvfmsVec_ch0.at(thisind);
  wfm_length=this_wf.size();
  auto this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch0.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 1; toggle=0;
  this_wf = fWvfmsVec_ch1.at(thisind);
  this_value =this_wf[0];
  //  std::cout << "starting : this+value " << this_value << " thresh " << threshold[i_ch] << std::endl;
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    //    std::cout << " i_t " << i_t << " this_value " << this_value << " toggle " << toggle << std::endl;
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch1.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 2; toggle=0;
  this_wf = fWvfmsVec_ch2.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch2.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 3; toggle=0;
  this_wf = fWvfmsVec_ch2.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch3.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 4; toggle=0;
  this_wf = fWvfmsVec_ch4.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch4.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 5; toggle=0;
  this_wf = fWvfmsVec_ch5.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch5.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 6; toggle=0;
  this_wf = fWvfmsVec_ch6.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch6.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 7; toggle=0;
  this_wf = fWvfmsVec_ch7.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch7.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 8; toggle=0;
  this_wf = fWvfmsVec_ch8.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch8.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 9; toggle=0;
  this_wf = fWvfmsVec_ch9.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch9.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 10; toggle=0;
  this_wf = fWvfmsVec_ch10.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch10.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 11; toggle=0;
  this_wf = fWvfmsVec_ch11.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch11.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 12; toggle=0;
  this_wf = fWvfmsVec_ch12.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch12.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 13; toggle=0;
  this_wf = fWvfmsVec_ch13.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch13.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 14; toggle=0;
  this_wf = fWvfmsVec_ch14.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch14.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }

  // find leading edges in waveforms
  i_ch = 15; toggle=0;
  this_wf = fWvfmsVec_ch15.at(thisind);
  this_value =this_wf[0];
  if (this_value>threshold[i_ch]) toggle=1;
  for(size_t i_t=1; i_t<wfm_length; ++i_t){
    this_value = this_wf[i_t];
    if (toggle==0 && this_value>threshold[i_ch]) {
      toggle=1;
      fPMT_ch15.emplace_back(i_t);
    }
    if (toggle==1 && this_value<threshold[i_ch]) toggle=0;
  }
  }// end if fragID==9
  fShift++;

}


void sbndaq::EventAna::reset_caen_variables() {
  fShift=0;   // fragment counter
  
  TTT.clear();
  TTT_ns.clear();// =0;  // will be set to value in CAEN fragement header
  caen_frag_ts.clear();
  ffragID.clear();  

  fPMT_ch0.clear();   fPMT_ch1.clear();   fPMT_ch2.clear(); fPMT_ch3.clear(); fPMT_ch4.clear();
  fPMT_ch5.clear();   fPMT_ch6.clear();   fPMT_ch7.clear(); fPMT_ch8.clear(); fPMT_ch9.clear();
  fPMT_ch10.clear();   fPMT_ch11.clear();   fPMT_ch12.clear(); fPMT_ch13.clear(); fPMT_ch14.clear();
  fPMT_ch15.clear();

  fWvfmsVec_ch0.clear(); fWvfmsVec_ch1.clear(); fWvfmsVec_ch2.clear(); fWvfmsVec_ch3.clear();
  fWvfmsVec_ch4.clear(); fWvfmsVec_ch5.clear(); fWvfmsVec_ch6.clear(); fWvfmsVec_ch7.clear();
  fWvfmsVec_ch8.clear(); fWvfmsVec_ch9.clear(); fWvfmsVec_ch10.clear(); fWvfmsVec_ch11.clear();
  fWvfmsVec_ch12.clear(); fWvfmsVec_ch13.clear(); fWvfmsVec_ch14.clear(); fWvfmsVec_ch15.clear();

  fTicksVec.clear();
}// end reset caen 1730

void sbndaq::EventAna::analyze_bern_fragment(artdaq::Fragment & frag)  {

  BernCRTFragmentV2 bern_fragment(frag);
  const BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);

    if (fcrt_keepall) {
      //metadata
      sequence_id.push_back(frag.sequenceID());
      run_start_time.push_back(            md->run_start_time());
      this_poll_start.push_back(           md->this_poll_start());
      this_poll_end.push_back(             md->this_poll_end());
      last_poll_start.push_back(           md->last_poll_start());
      last_poll_end.push_back(             md->last_poll_end());
      system_clock_deviation.push_back(    md->system_clock_deviation());
      feb_hits_in_poll.push_back(          md->hits_in_poll());
      feb_hits_in_fragment.push_back(      md->hits_in_fragment());
      timestamp.push_back(                 bevt->timestamp);
      //event info
      lostcpu.push_back(                   bevt->lostcpu);
      lostfpga.push_back(                  bevt->lostfpga);
      feb_hit_number.push_back(            bevt->feb_hit_number);
      last_accepted_timestamp.push_back(   bevt->last_accepted_timestamp);
      lost_hits.push_back(                 bevt->lost_hits);
    }

    mac5.push_back(                      md->MAC5());
    flags.push_back(                     bevt->flags);
    ts0.push_back(                       bevt->ts0);
    ts1.push_back(                       bevt->ts1);
    coinc.push_back(                     bevt->coinc);
    adc0.push_back(                      bevt->adc[0]);
    adc1.push_back(                      bevt->adc[1]);
    adc2.push_back(                      bevt->adc[2]);
    adc3.push_back(                      bevt->adc[3]);
    adc4.push_back(                      bevt->adc[4]);
    adc5.push_back(                      bevt->adc[5]);
    adc6.push_back(                      bevt->adc[6]);
    adc7.push_back(                      bevt->adc[7]);
    adc8.push_back(                      bevt->adc[8]);
    adc9.push_back(                      bevt->adc[9]);
    adc10.push_back(                      bevt->adc[10]);
    adc11.push_back(                      bevt->adc[11]);
    adc12.push_back(                      bevt->adc[12]);
    adc13.push_back(                      bevt->adc[13]);
    adc14.push_back(                      bevt->adc[14]);
    adc15.push_back(                      bevt->adc[15]);
    adc16.push_back(                      bevt->adc[16]);
    adc17.push_back(                      bevt->adc[17]);
    adc18.push_back(                      bevt->adc[18]);
    adc19.push_back(                      bevt->adc[19]);
    adc20.push_back(                      bevt->adc[20]);
    adc21.push_back(                      bevt->adc[21]);
    adc22.push_back(                      bevt->adc[22]);
    adc23.push_back(                      bevt->adc[23]);
    adc24.push_back(                      bevt->adc[24]);
    adc25.push_back(                      bevt->adc[25]);
    adc26.push_back(                      bevt->adc[26]);
    adc27.push_back(                      bevt->adc[27]);
    adc28.push_back(                      bevt->adc[28]);
    adc29.push_back(                      bevt->adc[29]);
    adc30.push_back(                      bevt->adc[30]);
    adc31.push_back(                      bevt->adc[31]);
    
    int max     = -std::numeric_limits<int>::max();
    int max_id  = std::numeric_limits<int>::max();
    int counter = 0;

    for(auto const &adc : bevt->adc)
      {
	if(adc > max)
	  {
	    max    = adc;
	    max_id = counter;
	  }
	++counter;
      }

    int max_other_sipm_id = max_id % 2 ? max_id - 1 : max_id + 1;
    
    max_adc.push_back({max, bevt->adc[max_other_sipm_id]});
    max_chan.push_back(max_id);

    
    
    if (fverbose) {
      std::cout << "  mac5                "     <<    (int)(md->MAC5())                   << std::endl;
      std::cout << "  run_start_time      "     <<  md->run_start_time()              << std::endl;
      std::cout << "  this_poll_start     "     << md->this_poll_start()            << std::endl;
      std::cout << "  this_poll_end        "    << md->this_poll_end()              << std::endl;
      std::cout << "  last_poll_start      "    << md->last_poll_start()            << std::endl;
      std::cout << "  last_poll_end          "  << md->last_poll_end()              << std::endl;
      std::cout << "  system_clock_deviation "  << md->system_clock_deviation()     << std::endl;
      std::cout << "  feb_hits_in_poll       "  << md->hits_in_poll()           << std::endl;
      std::cout << "  feb_hits_in_fragment  "   << md->hits_in_fragment()       << std::endl;
      std::cout << "  flags  "                  << (int)(bevt->flags)                   << std::endl;
      std::cout << "  lostcpu         "         << bevt->lostcpu                 << std::endl;
      std::cout << "  lostfpga    "             <<     bevt->lostfpga     << std::endl;
      std::cout << "  ts0         "             <<     bevt->ts0          << std::endl;
      std::cout << "  ts1                    "  <<     bevt->ts1             << std::endl;
      std::cout << "  coinc                  "  <<     bevt->coinc                   << std::endl;
      std::cout << "  feb_hit_number         "  <<     bevt->feb_hit_number          << std::endl;
      std::cout << "  timestamp              "  <<     bevt->timestamp               << std::endl;
      std::cout << "  last_accepted_timestamp"  <<     bevt->last_accepted_timestamp << std::endl;
      std::cout << "  lost_hits              "  <<     bevt->lost_hits               << std::endl;

      for(int ch=0; ch<32; ch++)
      	std::cout << "channel " << ch << " has adc value " << bevt->adc[ch] << std::endl;
    }// if verbose
    

  }// end loop over fragments

}//analyze_bern_fragment


void print_fragment_words(artdaq::Fragment& frag, size_t wordNum, size_t bitsPerWord ) {
  
  const __uint8_t* data_ptr = reinterpret_cast<const __uint8_t*>(frag.dataBegin());
  size_t wordCount = frag.dataSizeBytes() / (bitsPerWord / 8);
  for (size_t w = 0; w < wordCount; w++) {
    // Check if the current word index matches the specified word index
    if (w == wordNum) {
        // Print the bits for the specified n-bit word
      for (int i = (bitsPerWord-1); i >= 0; --i) {
            // Calculate the byte index and bit index
	size_t byteIndex = w * (bitsPerWord / 8) + (i / 8);
	int bitIndex = (i % 8); // Get the correct bit position in the byte
	// Print the bit
	std::cout << static_cast<int>((data_ptr[byteIndex] >> bitIndex) & 1);
	if (i % 8 == 0) {
	  std::cout << " "; 
	}
      }
      std::cout << std::endl; // New line after printing the word
    }
  }      
}

// Extract the PTB words/data from the artDAQ fragments
void sbndaq::EventAna::extract_triggers(artdaq::Fragment & frag) {

  // Construct PTB fragment overlay class giving us access to all the helpful decoder functions
  CTBFragment ptb_fragment(frag);

  if(fverbose){
    std::cout << "PTB Fragment ID: " << frag.sequenceID() << " TS: " << frag.timestamp()
	      << " Containing " << ptb_fragment.NWords() << " words" << std::endl;
  }

  /*********************
  *
  * Note: Below for the Timestamp conversion,
  *       The PTB TS is in UTC seconds since the Unix epoch in 50MHz clock ticks.
  *       This means to recover seconds since the Unix epoch use: sec_since_epoch = TS / 50e6
  *       and of course nanosec_since_epoch = (TS / 50e6) * 1e9 = TS * 20
  *
  * Note: The `CTBFragment` constructor grabs the chunk of memory the artDAQ fragment occupies.
  *       Since we know a priori the number of PTB words (from the TCP header) and the size (128b for all words)
  *       we can loop over the memory addresses, casting each 128b into the correct word type given by each word's `type`.
  *
  * There are 5 word types with the following bit arrangement:
  *       Feeback        = | 3b Word Type | 61b Payload | 64b Timestamp | (here Payload is split into code, source, payload1, payload2)
  *       LLT, HLT       = | 3b Word Type | 61b Payload | 64b Timestamp |
  *       Channel Status = | 3b Word Type | 64b Payload | 61b Timestamp | (Larger Payload to fit all input channels)
  *
  * Word Type:
  *       - 0x0 = Feedback/Error Word -> Errors from the firmware, should abort the run
  *       - 0x1 = Low Level Trigger (LLT) -> Holds a record of any asserted LLTs
  *       - 0x2 = High Level Trigger (HLT) -> Holds a record of any asserted HLTs
  *       - 0x3 = Channel status -> Holds a record of any asserted inputs
  *       - 0x7 = Timestamp Word -> No payload just a timestamp, these are periodic
  *
  * Note: Payload AND'd with mask the size of the expected number of bits just
  *       to make sure we don't get any unexpected garbage bits since we aren't using
  *       the full bits of the variable type e.g. uint64_t, unint16_t..
  *
  **********************/

  // Loop through all the PTB words in the fragment, casting to
  // one of the 5 word types. The 3 Msb hold the word type
  for ( size_t i = 0; i < ptb_fragment.NWords(); i++ ) {
    if (fverbose) std::cout << "PTB Word type [" << ptb_fragment.Word(i)->word_type << "]" << std::endl;
    //std::cout << "PTB Word type [" << ptb_fragment.Word(i)->word_type << "]" ;
    switch ( ptb_fragment.Word(i)->word_type ) {
      case 0x0 : // Feedback (errors) Word
        // Only get this word if something goes wrong at the firmware level requires expert knowledge
        // to interpret. The appearance of this word should have crashed the run.
        unknown_or_error_word = true;
        std::cout << "Feedback Word! Code: " << ptb_fragment.Feedback(i)->code
                  << " Source: "             << ptb_fragment.Feedback(i)->source
                  << " Payload: "            <<  ptb_fragment.Feedback(i)->payload
                  << " Timestamp: "          << ptb_fragment.TimeStamp(i) << std::endl;
        break;
      case 0x1 : // LL Trigger
        if (fverbose) std::cout << "LLT Payload: " << ptb_fragment.Trigger(i)->trigger_word << std::endl;
        llt_trigger.emplace_back( ptb_fragment.Trigger(i)->trigger_word & 0x1FFFFFFFFFFFFFFF ); // bit map of asserted LLTs
        llt_ts.emplace_back( ptb_fragment.TimeStamp(i) * 20 ); // Timestamp of the word
	//std::cout <<ptb_fragment.Trigger(i)->trigger_word << "  LLT Timestamp: " <<  ptb_fragment.TimeStamp(i) << "  " <<std::bitset<64>( ptb_fragment.TimeStamp(i) ) <<  std::endl;
        break;
      case 0x2 : // HL Trigger
        if (fverbose) std::cout << "HLT Payload: " << ptb_fragment.Trigger(i)->trigger_word << std::endl;
	if (fverbose) std::cout << "HLT ts: " << ptb_fragment.TimeStamp(i) << std::endl;
        hlt_gateCount.emplace_back(ptb_fragment.Trigger(i)->gate_counter & 0xFF);
        hlt_trigger.emplace_back( ptb_fragment.Trigger(i)->trigger_word & 0x1FFFFFFFFFFFFFFF );
        hlt_ts.emplace_back( ptb_fragment.TimeStamp(i) * 20 );
        ptb_frag_ts.push_back(frag.timestamp());
        hlt_word_count++;
	//std::cout << "HLT Word: " << ptb_fragment.Trigger(i)->trigger_word << " HLT GateCount: " << ptb_fragment.Trigger(i)->gate_counter  << " TS: " << ptb_fragment.TimeStamp(i) * 20<< " Prev TS: " <<ptb_fragment.PTBWord(i)->prevTS *20<< std::endl;
	break;
      case 0x3 : // Channel Status
        // Each PTB input gets a bit map e.g. CRT has 14 inputs and is 14b
        // (1 is channel asserted 0 otherwise)
        // TODO add MTCA and NIM channel status words
        auxpds_status.emplace_back( ptb_fragment.ChStatus(i)->pds & 0x3FF );
        crt_status.emplace_back( ptb_fragment.ChStatus(i)->crt & 0x3FFF );
	//std::cout <<std::bitset<14>(ptb_fragment.ChStatus(i)->crt) << "                        CRT Timestamp: " <<  ptb_fragment.TimeStamp(i) << "  " <<std::bitset<64>( ptb_fragment.TimeStamp(i)) <<  std::endl;
        beam_status.emplace_back( ptb_fragment.ChStatus(i)->beam & 0x3 );
        chan_stat_ts.emplace_back( ptb_fragment.TimeStamp(i) * 20 );
        break;
      case 0x7 : // Timestamp Word
        // We don't care about this word, it only has a TS and is sent periodically.
        ptb_fragment.TimeStamp(i);
        ts_word_count++;
        break;
      default : // Unknown, should never happen!
        unknown_or_error_word = true;
        std::cout << "Unknown PTB word type = " << ptb_fragment.Word(i)->word_type << std::endl;
    }
  }

}  // extract trigger fragments for the PTB

void sbndaq::EventAna::reset_ptb_variables() {

  // Initialize/reset the variables
  unknown_or_error_word = false;
  ts_word_count = 0;
  hlt_word_count = 0;
  ptb_frag_ts.clear();
  llt_trigger.clear();
  llt_ts.clear();
  hlt_trigger.clear();
  hlt_gateCount.clear();
  hlt_ts.clear();
  crt_status.clear();
  beam_status.clear();
  mtca_status.clear();
  nim_status.clear();
  auxpds_status.clear();
  chan_stat_ts.clear();

} // reset for ptb parameters

void sbndaq::EventAna::analyze_tdc_fragment(artdaq::Fragment & frag)  {

  const auto tsfrag = TDCTimestampFragment(frag);
  const auto ts = tsfrag.getTDCTimestamp();

  // each TDCTimstamp fragment has data from only one channel. The fragments are not always in time order

  if (ts->vals.channel==0)  ftdc_ch0.emplace_back(ts->nanoseconds());
  if (ts->vals.channel==1)  ftdc_ch1.emplace_back(ts->nanoseconds());
  if (ts->vals.channel==2)  ftdc_ch2.emplace_back(ts->nanoseconds());
  if (ts->vals.channel==3)  ftdc_ch3.emplace_back(ts->nanoseconds());
  if (ts->vals.channel==4)  ftdc_ch4.emplace_back(ts->nanoseconds());

  if(ftdc_utc){
    if (ts->vals.channel==0)  ftdc_ch0_utc.emplace_back(ts->timestamp_ns());
    if (ts->vals.channel==1)  ftdc_ch1_utc.emplace_back(ts->timestamp_ns());
    if (ts->vals.channel==2)  ftdc_ch2_utc.emplace_back(ts->timestamp_ns());
    if (ts->vals.channel==3)  ftdc_ch3_utc.emplace_back(ts->timestamp_ns());
    if (ts->vals.channel==4)  ftdc_ch4_utc.emplace_back(ts->timestamp_ns());
  }

  if(fverbose){
    std::cout << "=====================================" << std::endl;
    std::cout << "seq ID: " << frag.sequenceID() << std::endl;
    std::cout << "channel: " << ts->vals.channel << std::endl;
    std::cout << "name: " << ts->vals.name[0]
                          << ts->vals.name[1]
    			  << ts->vals.name[2]
    			  << ts->vals.name[3]
                          << ts->vals.name[4]
                          << ts->vals.name[5]
                          << ts->vals.name[6]
                          << ts->vals.name[7]
                          << std::endl;
    std::cout << "seconds: " << ts->vals.seconds << " s" << std::endl;
    std::cout << "coarse: " << ts->vals.coarse << " tick for 8ns/tick" << std::endl;
    std::cout << "frac: " << ts->vals.frac << " bit for 8ns/4096bit" << std::endl;
    std::cout << std::endl;
    std::cout << "fragment ts:  " << frag.timestamp() << std::endl;
    std::cout << "channel ts:   " << ts->timestamp_ns() << std::endl;
    std::cout << "channel frac: " << ts->nanoseconds() << std::endl;
    std::cout << "=====================================" << std::endl;
  };

}//end of analyze tdc timstamp fragment
void sbndaq::EventAna::analyze_ntb_fragment(artdaq::Fragment & frag)  {
  ntb_frag_ts.push_back(frag.timestamp());
}//end of analyze ntb fragment

DEFINE_ART_MODULE(sbndaq::EventAna)
