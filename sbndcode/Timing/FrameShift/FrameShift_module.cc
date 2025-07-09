////////////////////////////////////////////////////////////////////////
// Class:       FrameShift
// Plugin Type: Producer
// File:        FrameShift_module.cc
//
// Author: Lan Nguyen (vclnguyen@ucsb.edu)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "TTree.h"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

#include <bitset>

namespace sbnd {
  class FrameShift;
}


class sbnd::FrameShift : public art::EDProducer {
public:
  explicit FrameShift(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FrameShift(FrameShift const&) = delete;
  FrameShift(FrameShift&&) = delete;
  FrameShift& operator=(FrameShift const&) = delete;
  FrameShift& operator=(FrameShift&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void ResetEventVars();
  double SubtractUTCTimestmap(const uint64_t& ts1, const uint64_t& ts2);

private:

  // Declare member data here.

  //PTB
  art::InputTag fPtbDecodeLabel;
  
  std::vector<int> fPtbEtrigHlts; 
  std::vector<int> fPtbGateHlts; 
  std::vector<int> fPtbCrtT1Hlts; 
  
  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;
  std::vector<uint64_t> _ptb_hlt_unmask_timestamp;
  std::vector<int> _ptb_hlt_trunmask;

  int _hlt_etrig;
  uint64_t _hlt_etrig_ts;
  uint64_t _hlt_gate_ts;
  uint64_t _hlt_crtt1_ts;
 
  //TDC
  art::InputTag fTdcDecodeLabel;

  std::vector<uint64_t> _tdc_ch0;
  std::vector<uint64_t> _tdc_ch1;
  std::vector<uint64_t> _tdc_ch2;
  std::vector<uint64_t> _tdc_ch4;

  uint64_t _tdc_crtt1_ts;
  uint64_t _tdc_bes_ts;
  uint64_t _tdc_rwm_ts;
  uint64_t _tdc_etrig_ts;

  //Timing Reference
  art::InputTag fTimingRefPmtLabel;
  art::InputTag fTimingRefCrtLabel;

  uint16_t _pmt_timing_type; // e.g. SPECTDC = 0; PTB HLT = 1; CAEN-only = 3
  uint16_t _pmt_timing_ch;
  uint16_t _crt_timing_type; // e.g. SPECTDC = 0; PTB HLT = 1;
  uint16_t _crt_timing_ch;

  //Frame Shift
  uint16_t _global_timing_type; // e.g. SPECTDC = 0; PTB HLT = 1; 
  uint64_t _global_frame;

  double _frame_tdc_crtt1;
  double _frame_tdc_bes;
  double _frame_tdc_rwm;
  double _frame_hlt_crtt1;
  double _frame_hlt_gate;

  //Debug
  bool fDebugPtb;
  bool fDebugTdc;
  bool fDebugTimingRef;
  bool fDebugFrame;

  //---TREE PARAMETERS
  TTree *fTree;
  art::ServiceHandle<art::TFileService> tfs;
  int _run, _subrun, _event;
};


sbnd::FrameShift::FrameShift(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");
  fPtbDecodeLabel = p.get<art::InputTag>("PtbDecodeLabel", "ptbdecoder");
  fTimingRefPmtLabel = p.get<art::InputTag>("TimingRefPmtLabel", "pmtdecoder");
  fTimingRefCrtLabel = p.get<art::InputTag>("TimingRefCrtLabel", "crtstrips");

  fPtbEtrigHlts = p.get<std::vector<int>>("PtbEtrigHlts", {1, 2, 3, 4, 5, 14, 15});
  fPtbGateHlts = p.get<std::vector<int>>("PtbGateHlts", {26, 27});
  fPtbCrtT1Hlts = p.get<std::vector<int>>("PtbCrtT1Hlts", {20, 21});

  fDebugPtb = p.get<bool>("DebugPtb", false);
  fDebugTdc = p.get<bool>("DebugTdc", false);
  fDebugTimingRef = p.get<bool>("DebugTimingRef", false);
  fDebugFrame = p.get<bool>("DebugFrame", false);
  
  produces< raw::FrameShiftInfo >();

}

void sbnd::FrameShift::produce(art::Event& e)
{
  std::unique_ptr< raw::FrameShiftInfo > newFrameShiftInfo(new raw::FrameShiftInfo());

  ResetEventVars();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if (fDebugTdc | fDebugPtb | fDebugTimingRef | fDebugFrame)
        std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";

  //---------------------------TDC-----------------------------//
  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
  e.getByLabel(fTdcDecodeLabel, tdcHandle);

  if (!tdcHandle.isValid() || tdcHandle->size() == 0){
    if(fDebugTdc) {
      //TODO: Sometimes TDC products are not available, but this is not an error.
      throw cet::exception("FrameShift") << "No sbnd::timing::DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!";
    } 
  }
  else{
    std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
    for (size_t i=0; i<tdc_v.size(); i++){
      auto tdc = tdc_v[i];
      const uint32_t  ch = tdc.Channel();
      const uint64_t  ts = tdc.Timestamp();
      //Also TODO: Make use of picosecond timestamps
      //const uint64_t  ts_ps = tdc.TimestampPs();

      if (ch == 0) {
        _tdc_ch0.push_back(ts);
      }
      else if (ch == 1) {
        _tdc_ch1.push_back(ts);
      }
      else if (ch == 2) {
        _tdc_ch2.push_back(ts);
      }
      else if (ch == 4) {
        _tdc_ch4.push_back(ts);
      }
    }
  }
  tdcHandle.removeProduct();

  //TODO: What if saves more than one timestamp per channel?
  if (_tdc_ch0.size() > 0) _tdc_crtt1_ts = _tdc_ch0[0];
  if (_tdc_ch1.size() > 0) _tdc_bes_ts = _tdc_ch1[0];
  if (_tdc_ch2.size() > 0) _tdc_rwm_ts = _tdc_ch2[0];
  if (_tdc_ch4.size() > 0) _tdc_etrig_ts = _tdc_ch4[0];
  
  if (fDebugTdc){
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "TDC Channel 0 (CRTT1) Timestamp: " 
                      << " (s) = " << _tdc_crtt1_ts/uint64_t(1e9)
                      << ", (ns) = " << _tdc_crtt1_ts%uint64_t(1e9) 
                      << std::endl;
    std::cout << "TDC Channel 1 (BES) Timestamp: "
                      << " (s) = " << _tdc_bes_ts/uint64_t(1e9)
                      << ", (ns) = " << _tdc_bes_ts%uint64_t(1e9) 
                      << std::endl;       
    std::cout << "TDC Channel 2 (RWM) Timestamp: "
                      << " (s) = " << _tdc_rwm_ts/uint64_t(1e9)
                      << ", (ns) = " << _tdc_rwm_ts%uint64_t(1e9) 
                      << std::endl;
    std::cout << "TDC Channel 4 (ETRIG) Timestamp: "
                      << " (s) = " << _tdc_etrig_ts/uint64_t(1e9)
                      << ", (ns) = " << _tdc_etrig_ts%uint64_t(1e9)
                      << std::endl;
  } 
  //---------------------------PTB-----------------------------//
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  std::vector<art::Ptr<raw::ptb::sbndptb>> ptb_v;
  e.getByLabel(fPtbDecodeLabel, ptbHandle);
  
  if ((!ptbHandle.isValid() || ptbHandle->size() == 0)){
      //TODO: Sometimes PTB products are not available, but this is not an error.
      throw cet::exception("FrameShift") << "No raw::ptb::sbndptb found w/ tag " << fPtbDecodeLabel << ". Check data quality!";
  }
  else{
    art::fill_ptr_vector(ptb_v, ptbHandle);
    // HLT Words
    unsigned nHLTs = 0;

    for(auto const& ptb : ptb_v)
      nHLTs += ptb->GetNHLTriggers();
  
    _ptb_hlt_trigger.resize(nHLTs);
    _ptb_hlt_timestamp.resize(nHLTs);
    _ptb_hlt_trunmask.resize(nHLTs);
    _ptb_hlt_unmask_timestamp.resize(nHLTs);

    unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
    unsigned h_i = 0; //For trigger with bitmask

    for(auto const& ptb : ptb_v){
      for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i){
        _ptb_hlt_trigger[h_i] = ptb->GetHLTrigger(i).trigger_word;
        _ptb_hlt_timestamp[h_i] = ptb->GetHLTrigger(i).timestamp * uint64_t(20); //Units can be found in the Decoder Module 
        h_i++;
  
        int val = ptb->GetHLTrigger(i).trigger_word;
        int upBit[32];
    
        for (int u=0; u<32; u++){ //setting default values for maximum of 32 bits
          upBit[u]=-1;
        }
  
        int numOfTrig =0;
        for(int b=0; b<32;b++){
          if ((val & 0x01) ==1){
            upBit[numOfTrig] = b;
            numOfTrig++;
          }
          val = val >> 1;
        }
    
        if (numOfTrig ==1){
          _ptb_hlt_unmask_timestamp[hlt_i] = _ptb_hlt_timestamp[h_i-1];
          _ptb_hlt_trunmask[hlt_i] = upBit[0];
          hlt_i++;
        }//End of if statement for single upbit
        else if (numOfTrig > 1){
          nHLTs += (numOfTrig -1);
          _ptb_hlt_unmask_timestamp.resize(nHLTs);
          _ptb_hlt_trunmask.resize(nHLTs);
  
          for (int mult =0; mult < numOfTrig; mult++){ 
            _ptb_hlt_trunmask[hlt_i] = upBit[mult];
            _ptb_hlt_unmask_timestamp[hlt_i] = _ptb_hlt_timestamp[h_i-1];
            hlt_i++;
          } //End of loop over multiple upbits
        } //End of else statement for multiple triggers
      } //End of loop over nHLTriggers
    } //End of loop over ptb in ptb_v
    if (fDebugPtb){
      for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "HLT " << _ptb_hlt_trunmask[i] 
                            << " sec (s) = " << _ptb_hlt_unmask_timestamp[i]/uint64_t(1e9)
                            << ", ts (ns) = " << _ptb_hlt_unmask_timestamp[i]%uint64_t(1e9)
                            <<std::endl;
      }
    }
  }
  ptbHandle.removeProduct();

  for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
    for (size_t j = 0; j < fPtbEtrigHlts.size(); j++){
      if (_ptb_hlt_trunmask[i] == fPtbEtrigHlts[j]){
        _hlt_etrig = _ptb_hlt_trunmask[i];
        _hlt_etrig_ts = _ptb_hlt_unmask_timestamp[i];
        break;
      }
    }
  }
  
  for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
    for (size_t j = 0; j < fPtbGateHlts.size(); j++){
      if(_ptb_hlt_trunmask[i] == fPtbGateHlts[j]){
        _hlt_gate_ts = _ptb_hlt_unmask_timestamp[i];
        break;
      }
    }
  }
  
  for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
    for (size_t j = 0; j < fPtbCrtT1Hlts.size(); j++){
      if(_ptb_hlt_trunmask[i] == fPtbCrtT1Hlts[j]){
        _hlt_crtt1_ts = _ptb_hlt_unmask_timestamp[i];
        break;
      }
    }
  }

  //---------------------------TimingRef-----------------------------//
  art::Handle<raw::TimingReferenceInfo> timingRefPmtHandle;
  e.getByLabel(fTimingRefPmtLabel, timingRefPmtHandle);

  if (!timingRefPmtHandle.isValid()){
    throw cet::exception("FrameShift") << "No raw::TimingReferenceInfo found w/ tag " << fTimingRefPmtLabel << ". Check data quality!";
  }
  else{
    raw::TimingReferenceInfo const& pmt_timing(*timingRefPmtHandle);
    _pmt_timing_type = pmt_timing.timingType;
    _pmt_timing_ch = pmt_timing.timingChannel;
  }
  timingRefPmtHandle.removeProduct();

  art::Handle<raw::TimingReferenceInfo> timingRefCrtHandle;
  e.getByLabel(fTimingRefCrtLabel, timingRefCrtHandle);

  if (!timingRefCrtHandle.isValid()){
    throw cet::exception("FrameShift") << "No raw::TimingReferenceInfo found w/ tag " << fTimingRefCrtLabel << ". Check data quality!";
  }
  else{
    raw::TimingReferenceInfo const& crt_timing(*timingRefCrtHandle);
    _crt_timing_type = crt_timing.timingType;
    _crt_timing_ch = crt_timing.timingChannel;
  }
  timingRefCrtHandle.removeProduct();

  //Check what frame it is decoded to?
  if (_pmt_timing_type != _crt_timing_type){
    throw cet::exception("FrameShift") << "Timing Reference for PMT and CRT are not the same! PMT type = " 
                                       << _pmt_timing_type << ", CRT type = " << _crt_timing_type;
  }else{
    _global_timing_type = _pmt_timing_type;
  }

  if (fDebugTimingRef){
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Timing Reference For Decoding Pmt" << std::endl;
    std::cout << "   Type = " << _pmt_timing_type << " (SPECTDC = 0; PTB HLT = 1; CAEN-only = 3)." << std::endl;
    std::cout << "   Channel = " << _pmt_timing_ch << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Timing Reference For Decoding Crt" << std::endl;
    std::cout << "   Type = " << _crt_timing_type << " (SPECTDC = 0; PTB HLT = 1)." << std::endl;
    std::cout << "   Channel = " << _crt_timing_ch << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Global Timing Type = " << _global_timing_type << " (SPECTDC = 0; PTB HLT = 1;)." << std::endl;
  }

  //-----------------------Compute Frame Shift-----------------------//

  if (_global_timing_type == 0) _global_frame = _tdc_etrig_ts;
  else if (_global_timing_type == 1) _global_frame = _hlt_etrig_ts;

  //TODO: account for cable length
  //for beam stream
  if ((_hlt_etrig == 1) | (_hlt_etrig ==2)){
    if(_tdc_crtt1_ts != 0) _frame_tdc_crtt1 = SubtractUTCTimestmap(_global_frame, _tdc_crtt1_ts);
    if(_tdc_bes_ts != 0) _frame_tdc_bes = SubtractUTCTimestmap(_global_frame, _tdc_bes_ts);
    if(_tdc_rwm_ts != 0) _frame_tdc_rwm = SubtractUTCTimestmap(_global_frame, _tdc_rwm_ts);
    if(_hlt_crtt1_ts != 0) _frame_hlt_crtt1 = SubtractUTCTimestmap(_global_frame, _hlt_crtt1_ts);
    if(_hlt_gate_ts != 0) _frame_hlt_gate = SubtractUTCTimestmap(_global_frame, _hlt_gate_ts);
  }
  //for offbeam stream
  else if ((_hlt_etrig == 3) | (_hlt_etrig ==4)){
    if(_tdc_crtt1_ts != 0) _frame_tdc_crtt1 = SubtractUTCTimestmap(_global_frame, _tdc_crtt1_ts);
    if(_hlt_crtt1_ts != 0) _frame_hlt_crtt1 = SubtractUTCTimestmap(_global_frame, _hlt_crtt1_ts);
    if(_hlt_gate_ts != 0) _frame_hlt_gate = SubtractUTCTimestmap(_global_frame, _hlt_gate_ts);
  }
  //TODO: for cosmic stream

  if (fDebugFrame){
    std::cout << "--------------------------------------" << std::endl;
    if ((_hlt_etrig == 1) | (_hlt_etrig ==2)) std::cout << "This is beam stream!" << std::endl;
    if ((_hlt_etrig == 3) | (_hlt_etrig ==4)) std::cout << "This is offbeam stream!" << std::endl;

    std::cout << std::setprecision(9) << " Frame tdc crtt1 = " << _frame_tdc_crtt1 << std::endl;
    std::cout << std::setprecision(9) << " Frame tdc bes = " << _frame_tdc_bes << std::endl;
    std::cout << std::setprecision(9) << " Frame tdc rwm = " << _frame_tdc_rwm << std::endl;
    std::cout << std::setprecision(9) << " Frame hlt crtt1 = " << _frame_hlt_crtt1 << std::endl;
    std::cout << std::setprecision(9) << " Frame hlt gate = " << _frame_hlt_gate << std::endl;
  }

  //Put product in event
  newFrameShiftInfo->frameTdcCrtt1 = _frame_tdc_crtt1;
  newFrameShiftInfo->frameTdcBes = _frame_tdc_bes;
  newFrameShiftInfo->frameTdcRwm = _frame_tdc_rwm;
  newFrameShiftInfo->frameHltCrtt1 = _frame_hlt_crtt1;
  newFrameShiftInfo->frameHltBeamGate = _frame_hlt_gate;

  e.put(std::move(newFrameShiftInfo));

  //Fill the tree
  fTree->Fill();
}

void sbnd::FrameShift::beginJob()
{
  // Implementation of optional member function here.
  //Event Tree
  fTree = tfs->make<TTree>("events", "");
  
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);
  fTree->Branch("event", &_event);
  
  //TDC
  fTree->Branch("tdc_ch0", &_tdc_ch0);
  fTree->Branch("tdc_ch1", &_tdc_ch1);
  fTree->Branch("tdc_ch2", &_tdc_ch2);
  fTree->Branch("tdc_ch4", &_tdc_ch4);

  fTree->Branch("tdc_crtt1_ts", &_tdc_crtt1_ts);
  fTree->Branch("tdc_bes_ts", &_tdc_bes_ts);
  fTree->Branch("tdc_rwm_ts", &_tdc_rwm_ts);
  fTree->Branch("tdc_etrig_ts", &_tdc_etrig_ts);

  //PTB
  fTree->Branch("ptb_hlt_trigger", &_ptb_hlt_trigger);
  fTree->Branch("ptb_hlt_timestamp", &_ptb_hlt_timestamp);
  fTree->Branch("ptb_hlt_unmask_timestamp", &_ptb_hlt_unmask_timestamp);
  fTree->Branch("ptb_hlt_trunmask", &_ptb_hlt_trunmask);

  fTree->Branch("hlt_etrig", &_hlt_etrig);
  fTree->Branch("hlt_etrig_ts", &_hlt_etrig_ts);
  fTree->Branch("hlt_gate_ts", &_hlt_gate_ts);
  fTree->Branch("hlt_crtt1_ts", &_hlt_crtt1_ts);

  //Timing Reference
  fTree->Branch("crt_timing_ch", &_crt_timing_ch);
  fTree->Branch("crt_timing_type", &_crt_timing_type);
  fTree->Branch("pmt_timing_ch", &_pmt_timing_ch);
  fTree->Branch("pmt_timing_type", &_pmt_timing_type);

  //Frame Shift
  fTree->Branch("global_timing_type", &_global_timing_type);
  fTree->Branch("global_frame", &_global_frame);
  fTree->Branch("frame_tdc_crtt1", &_frame_tdc_crtt1);
  fTree->Branch("frame_tdc_bes", &_frame_tdc_bes);
  fTree->Branch("frame_tdc_rwm", &_frame_tdc_rwm);
  fTree->Branch("frame_hlt_crtt1", &_frame_hlt_crtt1);
  fTree->Branch("frame_hlt_gate", &_frame_hlt_gate);
}

void sbnd::FrameShift::ResetEventVars()
{
  _run = -1;
  _subrun = -1;
  _event = -1;

  _tdc_ch0.clear();
  _tdc_ch1.clear();
  _tdc_ch2.clear();
  _tdc_ch4.clear();

  _tdc_crtt1_ts = 0;
  _tdc_bes_ts = 0;
  _tdc_rwm_ts = 0;
  _tdc_etrig_ts = 0;

  _ptb_hlt_trigger.clear();
  _ptb_hlt_timestamp.clear();
  _ptb_hlt_unmask_timestamp.clear();
  _ptb_hlt_trunmask.clear();

  _hlt_etrig = -1;
  _hlt_etrig_ts = 0;
  _hlt_gate_ts = 0;
  _hlt_crtt1_ts = 0;

  _crt_timing_ch = 99;
  _crt_timing_type = 99;
  _pmt_timing_ch = 99;
  _pmt_timing_type = 99;
  
  _global_timing_type = 99;
  _global_frame = 0;

  _frame_tdc_crtt1 = 0;
  _frame_tdc_bes = 0;
  _frame_tdc_rwm = 0;
  _frame_hlt_crtt1 = 0;
  _frame_hlt_gate = 0;
}

double sbnd::FrameShift::SubtractUTCTimestmap(const uint64_t& ts1, const uint64_t& ts2)
{
  //Subtract two timestamps in UTC format
  //ts1 and ts2 are in nanoseconds in uint64_t format
  //Return the difference in nanoseconds as double

  double ts1_s = ts1 / uint64_t(1e9);
  double ts1_ns = ts1 % uint64_t(1e9);
  double ts2_s = ts2 / uint64_t(1e9);
  double ts2_ns = ts2 % uint64_t(1e9);

  //std::cout << std::setprecision(15) << "ts1_s = " << ts1_s << std::setprecision(9)<< ", ts1_ns = " << ts1_ns << std::endl;
  //std::cout << std::setprecision(15) << "ts2_s = " << ts2_s << std::setprecision(9)<< ", ts2_ns = " << ts2_ns << std::endl;
  
  double diff_s = 0;
  double diff_ns = 0;
  
  //If the same PPS, just subtract the nanoseconds
  if(ts1_s == ts2_s){
    diff_ns = ts1_ns - ts2_ns;
  }
  //If ts1 is later than ts2, then subtract the seconds and add the nanoseconds
  else{
    diff_s = ts1_s - ts2_s;
    diff_ns = diff_s * (double)1e9 + ts1_ns - ts2_ns;
  }

  return diff_ns;
}

DEFINE_ART_MODULE(sbnd::FrameShift)
