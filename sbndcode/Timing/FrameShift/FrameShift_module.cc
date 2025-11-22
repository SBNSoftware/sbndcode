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

#include "artdaq-core/Data/RawEvent.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/Timing/TimingInfo.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"

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
    
  static constexpr uint64_t InvalidTimestamp = std::numeric_limits<uint64_t>::max(); ///< Invalid frame.    
  static constexpr uint16_t InvalidTimingChannel = std::numeric_limits<uint16_t>::max(); ///< Invalid frame.    
  static constexpr uint16_t InvalidTimingType = std::numeric_limits<uint16_t>::max(); ///< Invalid frame.    

private:

  // Declare member data here.
  //DAQ Header
  std::string fDAQHeaderModuleLabel;
  std::string fDAQHeaderInstanceLabel;
  uint64_t fRawTSCorrection; //Correction for the Event timestamp to account for NTB
  uint64_t fMaxAllowedRefTimeDiff;
  uint64_t _raw_ts; //ns
  
  //Timing Reference
  art::InputTag fTimingRefPmtLabel;
  art::InputTag fTimingRefCrtLabel;

  //PTB
  art::InputTag fPtbDecodeLabel;
  
  std::vector<int> fPtbEtrigHlts; 
  
  std::vector<int> fBeamEtrigHlt;
  std::vector<int> fOffbeamEtrigHlt;
  std::vector<int> fXmuonEtrigHlt;
  
  int fBeamCrtT1Hlt;
  int fOffbeamCrtT1Hlt;
  
  int fBeamGateHlt;
  int fOffbeamGateHlt;
  
  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;
  std::vector<uint64_t> _ptb_hlt_unmask_timestamp;
  std::vector<int> _ptb_hlt_trunmask;

  std::vector<int> _ptb_hlt_etrig;
  std::vector<uint64_t> _ptb_hlt_etrig_ts;

  bool _isBeam;
  bool _isOffbeam;
  bool _isXmuon;  

  int _hlt_etrig;
  uint64_t _hlt_etrig_ts;
  int _hlt_gate;
  uint64_t _hlt_gate_ts;
  int _hlt_crtt1;
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

  //Frame Shift
  uint64_t _frame_crtt1; //ns
  uint16_t _timing_type_crtt1;
  uint16_t _timing_channel_crtt1;
  uint64_t _frame_gate; //ns
  uint16_t _timing_type_gate;
  uint16_t _timing_channel_gate;
  uint64_t _frame_etrig; //ns
  uint16_t _timing_type_etrig;
  uint16_t _timing_channel_etrig;

  //Value to apply at downstream modules depending on which stream
  uint64_t _frame_default;
  uint16_t _timing_type_default;
  uint16_t _timing_channel_default;

  //Value to shift Data to MC -- so that data agree with MC [ns]
  //TODO: Derive this value and verify if it is consistent across pmt/crt
  //TODO: Get this value from database instead of fhicl parameter
  uint64_t fShiftData2MC; //ns

  //Value to move RWM frame to agree with HLT Gate Frame
  //This is derived by subtracting: TDC RWM - HLT Gate.
  //Using the MC2025B dataset, this distribution has a mean of 1738 ns and std of 9 ns.
  //TODO: Get this value from database instead of fhicl parameter
  uint64_t fShiftRWM2Gate; //ns

  //Value to move TDC values to PTB HLT values
  uint64_t fShiftTDC2PTB; //ns

  //Debug
  bool fDebugDAQHeader;
  bool fDebugPtb;
  bool fDebugTdc;
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
  fDAQHeaderInstanceLabel = p.get<std::string>("DAQHeaderInstanceLabel");
  fDAQHeaderModuleLabel = p.get<std::string>("DAQHeaderModuleLabel");
  fRawTSCorrection = p.get<uint64_t>("RawTSCorrection");
  fMaxAllowedRefTimeDiff = p.get<uint64_t>("MaxAllowedRefTimeDiff");
  
  fTimingRefPmtLabel = p.get<art::InputTag>("TimingRefPmtLabel");
  fTimingRefCrtLabel = p.get<art::InputTag>("TimingRefCrtLabel");

  fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel");
  fPtbDecodeLabel = p.get<art::InputTag>("PtbDecodeLabel");

  fPtbEtrigHlts = p.get<std::vector<int>>("PtbEtrigHlts");
  
  fBeamEtrigHlt = p.get<std::vector<int>>("BeamEtrigHlt");
  fOffbeamEtrigHlt = p.get<std::vector<int>>("OffbeamEtrigHlt");
  fXmuonEtrigHlt = p.get<std::vector<int>>("XmuonEtrigHlt");
  
  fBeamCrtT1Hlt = p.get<int>("BeamCrtT1Hlt");
  fOffbeamCrtT1Hlt = p.get<int>("OffbeamCrtT1Hlt");

  fBeamGateHlt = p.get<int>("BeamGateHlt");
  fOffbeamGateHlt = p.get<int>("OffbeamGateHlt");

  fDebugDAQHeader = p.get<bool>("DebugDAQHeader", false);
  fDebugPtb = p.get<bool>("DebugPtb", false);
  fDebugTdc = p.get<bool>("DebugTdc", false);
  fDebugFrame = p.get<bool>("DebugFrame", false);
  
  //TODO: Get from database instead of fhicl parameters
  fShiftData2MC = p.get<uint64_t>("ShiftData2MC");
  fShiftRWM2Gate = p.get<uint64_t>("ShiftRWM2Gate"); 
  fShiftTDC2PTB = p.get<uint64_t>("ShiftTDC2PTB");
  
  produces< sbnd::timing::FrameShiftInfo >();
  produces< sbnd::timing::TimingInfo >();
}

void sbnd::FrameShift::produce(art::Event& e)
{

  ResetEventVars();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if (fDebugTdc | fDebugPtb | fDebugFrame)
        std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";

  //----------------------------DAQ Header-----------------------------//
  art::Handle<artdaq::detail::RawEventHeader> DAQHeaderHandle;
  e.getByLabel(fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel, DAQHeaderHandle);

  if (!DAQHeaderHandle.isValid()){
    throw cet::exception("FrameShift") << "No artdaq::detail::RawEventHeader found w/ tag " << fDAQHeaderModuleLabel << ". Check data quality!";
  }
  else{
      artdaq::RawEvent rawHeaderEvent = artdaq::RawEvent(*DAQHeaderHandle);
      _raw_ts = rawHeaderEvent.timestamp() - fRawTSCorrection;
  }
  if (fDebugDAQHeader){
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "DAQ Header Timestamp: "
                      << " (s) = " << _raw_ts/uint64_t(1e9)
                      << ", (ns) = " << _raw_ts%uint64_t(1e9) 
                      << std::endl;
  }

  //---------------------------TDC-----------------------------//
  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
  e.getByLabel(fTdcDecodeLabel, tdcHandle);

  if (!tdcHandle.isValid() || tdcHandle->size() == 0){
    mf::LogInfo("FrameShift") << "No sbnd::timing::DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!\n";
    //throw cet::exception("FrameShift") << "No sbnd::timing::DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!";
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

  //---------------------------PTB-----------------------------//
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  std::vector<art::Ptr<raw::ptb::sbndptb>> ptb_v;
  e.getByLabel(fPtbDecodeLabel, ptbHandle);
  
  if ((!ptbHandle.isValid() || ptbHandle->size() == 0)){
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
  }
  ptbHandle.removeProduct();

  if (fDebugPtb){
    for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << "HLT " << _ptb_hlt_trunmask[i] 
                          << " sec (s) = " << _ptb_hlt_unmask_timestamp[i]/uint64_t(1e9)
                          << ", ts (ns) = " << _ptb_hlt_unmask_timestamp[i]%uint64_t(1e9)
                          <<std::endl;
    }
  }

  //---------------------------Get ETRIG From TDC and PTB-----------------------------//

  // TDC ch4: ETRIG
  if (_tdc_ch4.size() == 0){
    throw cet::exception("FrameShift") << "No TDC ETRIG timestamps found! Check data quality!";
  }
  else if (_tdc_ch4.size() == 1){
    _tdc_etrig_ts = _tdc_ch4[0];
  }
  else if(_tdc_ch4.size() > 1){
    uint64_t min_diff = std::numeric_limits<uint64_t>::max();
    for(auto ts : _tdc_ch4){
      uint64_t diff = _raw_ts > ts ? _raw_ts - ts : ts - _raw_ts; //raw_ts must be found for every event else throw exception
      if(diff < min_diff)
      {
        min_diff    = diff;
        _tdc_etrig_ts = ts;
      }
    }
  }

  //HLT ETRIG 
  //Grab all the ETRIG HLTS -- there might be more than 1
  for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
    for (size_t j = 0; j < fPtbEtrigHlts.size(); j++){
      if (_ptb_hlt_trunmask[i] == fPtbEtrigHlts[j]){
        _ptb_hlt_etrig.push_back(_ptb_hlt_trunmask[i]);
        _ptb_hlt_etrig_ts.push_back(_ptb_hlt_unmask_timestamp[i]);
      }
    }
  }


  if (_ptb_hlt_etrig.size() == 0){
    throw cet::exception("FrameShift") << "No HLT ETRIG timestamps found! Check data quality!";
  }
  else if (_ptb_hlt_etrig.size() == 1){
    _hlt_etrig = _ptb_hlt_etrig[0];
    _hlt_etrig_ts = _ptb_hlt_etrig_ts[0];
  }
  else if(_ptb_hlt_etrig.size() > 1){
    uint64_t min_diff = std::numeric_limits<uint64_t>::max();
    for(size_t i = 0; i < _ptb_hlt_etrig.size(); i++){
      uint64_t diff = _raw_ts > _ptb_hlt_etrig_ts[i] ? _raw_ts - _ptb_hlt_etrig_ts[i] : _ptb_hlt_etrig_ts[i] - _raw_ts;
      if(diff < min_diff)
      {
        min_diff    = diff;
        _hlt_etrig = _ptb_hlt_etrig[i];
        _hlt_etrig_ts = _ptb_hlt_etrig_ts[i];
      }
    }
  }

  //Decide which global frame to use as reference 
  uint64_t _global_frame = InvalidTimestamp;
  if (_tdc_etrig_ts != InvalidTimestamp){ _global_frame = _tdc_etrig_ts; }
  else if (_hlt_etrig_ts != InvalidTimestamp){ _global_frame = _hlt_etrig_ts; }
  else { _global_frame = _raw_ts;}

  if (fDebugFrame){
    std::cout << "----------------------------------------------------" << std::endl;
    if (_tdc_etrig_ts != InvalidTimestamp){
      std::cout << "Using TDC ETRIG as Global Frame Reference" << std::endl;
    }
    else if (_hlt_etrig_ts != InvalidTimestamp){
      std::cout << "Using PTB HLT ETRIG as Global Frame Reference" << std::endl;
    }
    else {
      std::cout << "Using DAQ Header Timestamp as Global Frame Reference" << std::endl;
    }
    std::cout << "Global Frame Timestamp: "
                      << " (s) = " << _global_frame/uint64_t(1e9)
                      << ", (ns) = " << _global_frame%uint64_t(1e9) 
                      << std::endl;
  }

  //---------------------------TDC Frame-----------------------------//
  // ch0: CRT T1
  if (_tdc_ch0.size() == 1){
    _tdc_crtt1_ts = _tdc_ch0[0];
  }
  else if(_tdc_ch0.size() > 1){
    //Get the one closest to the global frame
    uint64_t min_diff = std::numeric_limits<uint64_t>::max();
    for(auto ts : _tdc_ch0){
      uint64_t diff = _global_frame > ts ? _global_frame - ts : ts - _global_frame;
      if(diff < min_diff)
      {
        min_diff    = diff;
        _tdc_crtt1_ts = ts;
      }
    }
  }

  // ch1: BES
  if (_tdc_ch1.size() == 1){
    _tdc_bes_ts = _tdc_ch1[0];
  }
  else if(_tdc_ch1.size() > 1){
    //Get the one closest to the global frame
    uint64_t min_diff = std::numeric_limits<uint64_t>::max();
    for(auto ts : _tdc_ch1){
      uint64_t diff = _global_frame > ts ? _global_frame - ts : ts - _global_frame;
      if(diff < min_diff)
      {
        min_diff    = diff;
        _tdc_bes_ts = ts;
      }
    }
  }

  // ch2: RWM
  if (_tdc_ch2.size() == 1){
    _tdc_rwm_ts = _tdc_ch2[0];
  }
  else if(_tdc_ch2.size() > 1){
    //Get the one closest to the global frame
    uint64_t min_diff = std::numeric_limits<uint64_t>::max();
    for(auto ts : _tdc_ch2){
      uint64_t diff = _global_frame > ts ? _global_frame - ts : ts - _global_frame;
      if(diff < min_diff)
      {
        min_diff    = diff;
        _tdc_rwm_ts = ts;
      }
    }
  }

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
    std::cout << "----------------------------------------------------" << std::endl;
  } 
  
  //---------------------------PTB Frame-----------------------------//

  //Check which Gate/CRT T1 HLT to use based on ETRIG HLT
  //Order to check: Beam -> Offbeam -> Xmuon
  for (size_t i = 0; i < fBeamEtrigHlt.size(); i++){
    if (_hlt_etrig == fBeamEtrigHlt[i]){
      _hlt_gate = fBeamGateHlt;
      _hlt_crtt1 = fBeamCrtT1Hlt;
      _isBeam = true;
      break;
    }
  }
  if (!_isBeam){
    for (size_t i = 0; i < fOffbeamEtrigHlt.size(); i++){
      if (_hlt_etrig == fOffbeamEtrigHlt[i]){
        _hlt_gate = fOffbeamGateHlt;
        _hlt_crtt1 = fOffbeamCrtT1Hlt;
        _isOffbeam = true;
        break;
      }
    }
  }
  if (!_isBeam & !_isOffbeam){
    for (size_t i = 0; i < fXmuonEtrigHlt.size(); i++){
      if (_hlt_etrig == fXmuonEtrigHlt[i]){
        _isXmuon = true;
        break;
      }
    }
  } 

  if( !_isBeam & !_isOffbeam & !_isXmuon){
    throw cet::exception("FrameShift") << "ETRIG HLT " << _hlt_etrig << " does not match any known Beam/Offbeam/Xmuon ETRIG HLT! Check data quality!";
  }

  //Get Gate and CRT T1 HLT timestamps 
  //TODO: What if there is no Gate or CRT T1 HLT?
  for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
    if(_ptb_hlt_trunmask[i] == _hlt_gate){
      _hlt_gate_ts = _ptb_hlt_unmask_timestamp[i];
    }
    if(_ptb_hlt_trunmask[i] == _hlt_crtt1){
      _hlt_crtt1_ts = _ptb_hlt_unmask_timestamp[i];
    }
  }

  if (fDebugPtb){

    std::cout << "----------------------------------------------------" << std::endl;
    if (_isBeam) std::cout << "This is Beam Stream!" << std::endl; 
    if (_isOffbeam) std::cout << "This is Offbeam Stream!" << std::endl; 
    std::cout << "HLT ETRIG = " << _hlt_etrig
                      << ", Timestamp: "
                      << " (s) = " << _hlt_etrig_ts/uint64_t(1e9)
                      << ", (ns) = " << _hlt_etrig_ts%uint64_t(1e9) 
                      << std::endl;
    std::cout << "HLT Gate = " << _hlt_gate
                      << ", Timestamp: "
                      << " (s) = " << _hlt_gate_ts/uint64_t(1e9)
                      << ", (ns) = " << _hlt_gate_ts%uint64_t(1e9) 
                      << std::endl;
    std::cout << "HLT CRT T1 = " << _hlt_crtt1
                      << ", Timestamp: "
                      << " (s) = " << _hlt_crtt1_ts/uint64_t(1e9)
                      << ", (ns) = " << _hlt_crtt1_ts%uint64_t(1e9) 
                      << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
  }

  //-----------------------Pick default frame-----------------------//
  // The follow picks which frame to apply at downstream stage and store it as frame_default, based on the stream
  //
  // 1. Beam Stream: recontruct the beam spill + porch relative to RWM record in TDC since it has better resolution. 
  // Then, apply a constant to shift the RWM fram 
  // i.e. frame_tdc_rwm + shift rwm to beam gate HLT 26
  //
  // 2. Offbeam Stream: reconstruct the porch relative to the beam gate opening frame
  // i.e. frame_hlt_gate equiv
  //
  // 3. Xmuon Stream: shift to etrig
  //
  // TODO: Align Data Beam and Offbeam with MC with frame_data2mc
  //    + MC: t = 0 = first proton in spill
  //    + Data: t = 0 = abitrary. All subsystem electronics time is reference to the last PPS
    
  if (_isBeam){

    //Frame CRT T1
    if(_tdc_crtt1_ts != InvalidTimestamp){
      _frame_crtt1 = _tdc_crtt1_ts; //TODO: Add shift from TDC to PTB
      _timing_type_crtt1 = 0; //SPECTDC
      _timing_channel_crtt1 = 0;
    }
    else if(_hlt_crtt1_ts != InvalidTimestamp) {
      _frame_crtt1 = _hlt_crtt1_ts;
      _timing_type_crtt1 = 1; //PTB HLT
      _timing_channel_crtt1 = _hlt_crtt1;
    }
    else{
      _frame_crtt1 = InvalidTimestamp;
      _timing_type_crtt1 = 2;
      _timing_channel_crtt1 = InvalidTimingChannel;
    }

    //Frame Beam Gate
    if(_tdc_rwm_ts != InvalidTimestamp){
      _frame_gate = _tdc_rwm_ts + fShiftRWM2Gate; //TODO: + fShiftData2MC;
      _timing_type_gate = 0; //SPECTDC
      _timing_channel_gate = 2;
    }
    else if(_hlt_gate_ts != InvalidTimestamp){
      _frame_gate = _hlt_gate_ts;
      _timing_type_gate = 1; //PTB HLT
      _timing_channel_gate = _hlt_gate;
    }else{
      _frame_gate = InvalidTimestamp;
      _timing_type_gate = 2;
      _timing_channel_gate = InvalidTimingChannel;
    }

    //Frame ETRIG
    if(_tdc_etrig_ts != InvalidTimestamp){
      _frame_etrig = _tdc_etrig_ts + fShiftTDC2PTB;
      _timing_type_etrig = 0; //SPECTDC
      _timing_channel_etrig = 4;
    }
    else if(_hlt_etrig_ts != InvalidTimestamp){
      _frame_etrig = _hlt_etrig_ts;
      _timing_type_etrig = 1; //PTB HLT
      _timing_channel_etrig = _hlt_etrig;
    }
    else {
      _frame_etrig = InvalidTimestamp;
      _timing_type_etrig = 2;
      _timing_channel_etrig = InvalidTimingChannel;
    }

    //Pick default stream -- beam gate
    _frame_default = _frame_gate;
    _timing_type_default = _timing_type_gate;
    _timing_channel_default = _timing_channel_gate;
  }
  else if (_isOffbeam){

    //Frame CRT T1
    if(_tdc_crtt1_ts != InvalidTimestamp){
      _frame_crtt1 = _tdc_crtt1_ts; //TODO: Add shift from TDC to PTB
      _timing_type_crtt1 = 0; //SPECTDC
      _timing_channel_crtt1 = 0;
    }
    else if(_hlt_crtt1_ts != InvalidTimestamp) {
      _frame_crtt1 = _hlt_crtt1_ts;
      _timing_type_crtt1 = 1; //PTB HLT
      _timing_channel_crtt1 = _hlt_crtt1;
    }
    else{
      _frame_crtt1 = InvalidTimestamp;
      _timing_type_crtt1 = 2;
      _timing_channel_crtt1 = InvalidTimingChannel;
    }

    //Frame Gate -- TODO: I think HLT Gate is recorded in TDC as FTRIG and can be found
    if(_hlt_gate_ts != InvalidTimestamp) {
      _frame_gate = _hlt_gate_ts; // TODO: + fShiftData2MC;
      _timing_type_gate = 1; //PTB HLT
      _timing_channel_gate = _hlt_gate;
    }
    else {
      _frame_gate = InvalidTimestamp;
      _timing_type_gate = 2;
      _timing_channel_gate = InvalidTimingChannel;
    }

    //Frame ETRIG
    if(_tdc_etrig_ts != InvalidTimestamp) {
      _frame_etrig = _tdc_etrig_ts + fShiftTDC2PTB;
      _timing_type_etrig = 0; //SPECTDC
      _timing_channel_etrig = 4;
    }
    else if(_hlt_etrig_ts !=0) {
      _frame_etrig = _hlt_etrig_ts;
      _timing_type_etrig = 1; //PTB HLT
      _timing_channel_etrig = _hlt_etrig;
    }
    else {
      _frame_etrig = InvalidTimestamp;
      _timing_type_etrig = 2;
      _timing_channel_etrig = InvalidTimingChannel;
    }

    //Pick default stream -- beam gate
    _frame_default = _frame_gate;
    _timing_type_default = _timing_type_gate;
    _timing_channel_default = _timing_channel_gate;
  }
  else if (_isXmuon){

    //Frame ETRIG
    if(_tdc_etrig_ts != InvalidTimestamp) {
      _frame_etrig = _tdc_etrig_ts + fShiftTDC2PTB;
      _timing_type_etrig = 0; //SPECTDC
      _timing_channel_etrig = 4;
    }
    else if(_hlt_etrig_ts != InvalidTimestamp) {
      _frame_etrig = _hlt_etrig_ts;
      _timing_type_etrig = 1; //PTB HLT
      _timing_channel_etrig = _hlt_etrig;
    }
    else {
      _frame_etrig = InvalidTimestamp;
      _timing_type_etrig = 2;
      _timing_channel_etrig = InvalidTimingChannel;
    }

    //Pick default stream -- ETRIG
    _frame_default = _frame_etrig;
    _timing_type_default = _timing_type_etrig;
    _timing_channel_etrig = _timing_channel_etrig;
  }

  if (fDebugFrame){
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Frame Shift Results:" << std::endl;
    std::cout << "Frame CRT T1 type " << _timing_type_crtt1
                      <<", channel " << _timing_channel_crtt1
                      << ": (s) = " << _frame_crtt1/uint64_t(1e9) 
                      << ", (ns) = " << _frame_crtt1%uint64_t(1e9) 
                      << std::endl;
    std::cout << "Frame Beam Gate  type " << _timing_type_gate
                      <<", channel " << _timing_channel_gate
                      << ": (s) = " << _frame_gate/uint64_t(1e9) 
                      << ", (ns) = " << _frame_gate%uint64_t(1e9) 
                      << std::endl;
    std::cout << "Frame ETRIG type " << _timing_type_etrig
                      <<", channel " << _timing_channel_etrig
                      << ": (s) = " << _frame_etrig/uint64_t(1e9) 
                      << ", (ns) = " << _frame_etrig%uint64_t(1e9) 
                      << std::endl;
    std::cout << "Default Frame type " << _timing_type_default
                      <<", channel " << _timing_channel_default
                      << " : (s) = " << _frame_default/uint64_t(1e9) 
                      << ", (ns) = " << _frame_default%uint64_t(1e9) 
                      << std::endl;
    std::cout << "--------------------------------------" << std::endl;
  }
  

  //Put product in event
  std::unique_ptr< sbnd::timing::FrameShiftInfo > newFrameShiftInfo(new sbnd::timing::FrameShiftInfo(_frame_crtt1, _timing_type_crtt1, _timing_channel_crtt1,
                                                                                      _frame_gate, _timing_type_gate, _timing_channel_gate,
                                                                                      _frame_etrig, _timing_type_etrig,_timing_channel_etrig,
                                                                                      _frame_default, _timing_type_default, _timing_channel_default));

  std::unique_ptr< sbnd::timing::TimingInfo > newTimingInfo(new sbnd::timing::TimingInfo(_raw_ts, _tdc_crtt1_ts, _tdc_bes_ts, _tdc_rwm_ts, _tdc_etrig_ts, _hlt_crtt1_ts, _hlt_etrig_ts, _hlt_gate_ts));

  e.put(std::move(newTimingInfo));
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
  fTree->Branch("ptb_hlt_unmask_timestamp", &_ptb_hlt_unmask_timestamp);
  fTree->Branch("ptb_hlt_trunmask", &_ptb_hlt_trunmask);

  fTree->Branch("hlt_etrig", &_hlt_etrig);
  fTree->Branch("hlt_etrig_ts", &_hlt_etrig_ts);
  fTree->Branch("hlt_gate", &_hlt_gate);
  fTree->Branch("hlt_gate_ts", &_hlt_gate_ts);
  fTree->Branch("hlt_crtt1", &_hlt_crtt1);
  fTree->Branch("hlt_crtt1_ts", &_hlt_crtt1_ts);

  //Frame Shift
  fTree->Branch("frame_crtt1", &_frame_crtt1);
  fTree->Branch("timing_type_crtt1", &_timing_type_crtt1);
  fTree->Branch("timing_channel_crtt1", &_timing_channel_crtt1);

  fTree->Branch("frame_gate", &_frame_gate);
  fTree->Branch("timing_type_gate", &_timing_type_gate);
  fTree->Branch("timing_channel_gate", &_timing_channel_gate);

  fTree->Branch("frame_etrig", &_frame_etrig);  
  fTree->Branch("timing_type_etrig", &_timing_type_etrig);
  fTree->Branch("timing_channel_etrig", &_timing_channel_etrig);

  fTree->Branch("frame_default", &_frame_default);
  fTree->Branch("timing_type_default", &_timing_type_default);
  fTree->Branch("timing_channel_default", &_timing_channel_default);
}

void sbnd::FrameShift::ResetEventVars()
{
  _run = -1;
  _subrun = -1;
  _event = -1;

  _raw_ts = InvalidTimestamp;

  _tdc_ch0.clear();
  _tdc_ch1.clear();
  _tdc_ch2.clear();
  _tdc_ch4.clear();

  _tdc_crtt1_ts = InvalidTimestamp;
  _tdc_bes_ts = InvalidTimestamp;
  _tdc_rwm_ts = InvalidTimestamp;
  _tdc_etrig_ts = InvalidTimestamp;

  _ptb_hlt_trigger.clear();
  _ptb_hlt_timestamp.clear();
  _ptb_hlt_unmask_timestamp.clear();
  _ptb_hlt_trunmask.clear();

  _ptb_hlt_etrig.clear();
  _ptb_hlt_etrig_ts.clear();

  _isBeam = false;
  _isOffbeam = false;
  _isXmuon = false;  

  _hlt_etrig = std::numeric_limits<int>::max();
  _hlt_etrig_ts = InvalidTimestamp;
  _hlt_gate = std::numeric_limits<int>::max();
  _hlt_gate_ts = InvalidTimestamp;
  _hlt_crtt1= std::numeric_limits<int>::max();
  _hlt_crtt1_ts = InvalidTimestamp;
  
  _frame_crtt1 = InvalidTimestamp;
  _timing_type_crtt1 = InvalidTimingType;
  _timing_channel_crtt1 = InvalidTimingChannel;
  
  _frame_gate = InvalidTimestamp;
  _timing_type_gate = InvalidTimingType;
  _timing_channel_gate = InvalidTimingChannel;
  
  _frame_etrig = InvalidTimestamp;
  _timing_type_etrig = InvalidTimingType;
  _timing_channel_etrig = InvalidTimingChannel;

  _frame_default = InvalidTimestamp;
  _timing_type_default = InvalidTimingType;
  _timing_channel_default = InvalidTimingChannel;
}

DEFINE_ART_MODULE(sbnd::FrameShift)
