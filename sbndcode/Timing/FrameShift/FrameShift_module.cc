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

namespace sbnd::timing {
  class FrameShift;
}


class sbnd::timing::FrameShift : public art::EDProducer {
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

  uint64_t FindClosest(const std::vector<uint64_t> &timestamps, const uint64_t &reference);
  std::string PrintFormatTimestamp(const uint64_t &timestamp);

  void GetRawTimestamp(const art::Event &e);
  void GetTDCTimestamps(const art::Event &e);
  void GetPTBTimestamps(const art::Event &e);
  void FindETRIGs();

  uint64_t DecideGlobalEtrigTimestamp();
  void DecideRelevantTDCTimestamps(const uint64_t &global_etrig_ts);
  void DecideRelevantPTBTimestamps(const uint64_t &global_etrig_ts);

private:

  // FCL Parameters
  std::string fDAQHeaderModuleLabel;   // Instance label for DAQ header product
  std::string fDAQHeaderInstanceLabel; // Module label for DAQ header product
  uint64_t fRawTSCorrection;           // Correction to the DAQ header timestamp to account for NTB offset
  art::InputTag fPtbDecodeLabel;       // Module label for decoded PTB data
  std::vector<int> fPtbEtrigHlts;      // All allowed PTB HLT codes for event triggers
  std::vector<int> fBeamEtrigHlt;      // PTB HLT codes associated with beam trigger
  std::vector<int> fOffbeamEtrigHlt;   // PTB HLT codes associated with offbeam trigger
  std::vector<int> fXmuonEtrigHlt;     // PTB HLT codes associated with crossing muon trigger
  int fBeamCrtT1Hlt;                   // PTB HLT code for beam CRT T1 reset
  int fOffbeamCrtT1Hlt;                // PTB HLT code for offbeam CRT T1 reset
  int fBeamGateHlt;                    // PTB HLT code for beam gate opening
  int fOffbeamGateHlt;                 // PTB HLT code for offbeam gate opening
  art::InputTag fTdcDecodeLabel;       // Module label for decoded TDC data
  uint64_t fShiftData2MC;              // Value to shift Data to MC -- so that data agree with MC [ns]
                                       // TODO: Derive this value and verify if it is consistent across pmt/crt
                                       // TODO: Get this value from database instead of fhicl parameter
  uint64_t fShiftTdcRwm2PtbGate;       // Value to move Tdc RWM (ch2) frame to agree with HLT Gate Frame
                                       // This is derived by subtracting: TDC RWM - HLT Gate.
                                       // Using the MC2025B dataset, this distribution has a mean of 1738 ns and std of 9 ns.
                                       // TODO: Get this value from database instead of fhicl parameter
  uint64_t fShiftTdcEtrig2PtbEtrig;    // Value to move Tdc Etrig (ch4) to agree with HLT ETRIG
                                       // TODO: Derive + edit fcl
                                       // TODO: Get from database
  uint64_t fShiftTdcCrtt12PtbCrtt1;    // Value to move Tdc Crt T1 Reset (ch0) to agree with HLT CRT T1
                                       // TODO: Derive + edit fcl
                                       // TODO: Get from database
  bool fMakeTree;                      // Whether to produce a TTree in the hist file
  bool fDebugDAQHeader;                // Whether to print debug statements relevant to DAQ header
  bool fDebugPtb;                      // Whether to print debug statements relevant to PTB
  bool fDebugTdc;                      // Whether to print debug statements relevant to TDC
  bool fDebugFrame;                    // Whether to print debug statements relevant to frame shifts

  // Global Variables, set in processing
  int _run, _subrun, _event;                       // Stores the unique run, subrun, event number combination for this event

  uint64_t _raw_ts;                                // Stores DAQ header timestamp
  std::vector<uint64_t> _ptb_hlt_trigger;          // Stores full trigger word from PTB HLT (indexed per HLT)
  std::vector<uint64_t> _ptb_hlt_timestamp;        // Stores timestamp associated with PTB HLT (indexed per HLT)
  std::vector<uint64_t> _ptb_hlt_unmask_timestamp; // Stores timestamp associated with PTB HLT for each element (indexed per HLT element once unmasked)
  std::vector<int> _ptb_hlt_trunmask;              // Stores PTB HLT code for each element (indexed per HLT element once unmasked)
  std::vector<int> _ptb_hlt_etrig;                 // Stores all PTB HLTs that are an 'allowed type' for this stream
  std::vector<uint64_t> _ptb_hlt_etrig_ts;         // Stores the associated timestamps for the above HLTs

  bool _isBeam;                                    // Event is a beam trigger (BNB or BNBLight)
  bool _isOffbeam;                                 // Event is an offbeam trigger (OffBeam or OffBeamLight)
  bool _isXmuon;                                   // Event is a crossing muon trigger

  int _hlt_etrig;                                  // Stores the HLT event trigger code, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _hlt_etrig_ts;                          // Stores the HLT timestamp, once we've decided which one was closest to the DAQ header timestamp
  int _hlt_gate;                                   // Stores the HLT gate opening code, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _hlt_gate_ts;                           // Stores the HLT gate opening timestamp, once we've decided which one was closest to the DAQ header timestamp
  int _hlt_crtt1;                                  // Stores the HLT CRT T1 code, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _hlt_crtt1_ts;                          // Stores the HLT CRT T1 timestamp, once we've decided which one was closest to the DAQ header timestamp
 
  std::vector<uint64_t> _tdc_ch0;                  // Stores all the timestamps recorded in the TDC channel 0 (CRT T1)
  std::vector<uint64_t> _tdc_ch1;                  // Stores all the timestamps recorded in the TDC channel 1 (BES)
  std::vector<uint64_t> _tdc_ch2;                  // Stores all the timestamps recorded in the TDC channel 2 (RWM)
  std::vector<uint64_t> _tdc_ch4;                  // Stores all the timestamps recorded in the TDC channel 4 (Event trigger)

  uint64_t _tdc_crtt1_ts;                          // Stores the TDC CRT T1 timestamp, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _tdc_bes_ts;                            // Stores the TDC BES timestamp, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _tdc_rwm_ts;                            // Stores the TDC RWM timestamp, once we've decided which one was closest to the DAQ header timestamp
  uint64_t _tdc_etrig_ts;                          // Stores the TDC event trigger timestamp, once we've decided which one was closest to the DAQ header timestamp

  uint64_t _frame_crtt1;                           // Stores the shift required to move from the PPS frame to the CRT T1 frame
  uint16_t _timing_type_crtt1;                     // Stores the type (TDC/PTB) used to produce the above shift
  uint16_t _timing_channel_crtt1;                  // Stores the channel (TDC) / code (PTB) used to produce the above shift
  uint64_t _frame_gate;                            // Stores the shift required to move from the PPS frame to the gate opening frame
  uint16_t _timing_type_gate;                      // Stores the type (TDC/PTB) used to produce the above shift
  uint16_t _timing_channel_gate;                   // Stores the channel (TDC) / code (PTB) used to produce the above shift
  uint64_t _frame_etrig;                           // Stores the shift required to move from the PPS frame to the event trigger frame
  uint16_t _timing_type_etrig;                     // Stores the type (TDC/PTB) used to produce the above shift
  uint16_t _timing_channel_etrig;                  // Stores the channel (TDC) / code (PTB) used to produce the above shift

  uint64_t _frame_default;                         // Stores the shift required to move from the PPS frame to the default frame for this stream
  uint16_t _timing_type_default;                   // Stores the type (TDC/PTB) used to produce the above shift
  uint16_t _timing_channel_default;                // Stores the channel (TDC) / code (PTB) used to produce the above shift

  // Tree production
  TTree *fTree;
  art::ServiceHandle<art::TFileService> tfs;

  // Useful value
  static constexpr uint64_t kSecondInNanoseconds = static_cast<uint64_t>(1e9);
};


sbnd::timing::FrameShift::FrameShift(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  fDAQHeaderInstanceLabel = p.get<std::string>("DAQHeaderInstanceLabel");
  fDAQHeaderModuleLabel = p.get<std::string>("DAQHeaderModuleLabel");
  fRawTSCorrection = p.get<uint64_t>("RawTSCorrection");
  fPtbDecodeLabel = p.get<art::InputTag>("PtbDecodeLabel");
  fPtbEtrigHlts = p.get<std::vector<int>>("PtbEtrigHlts");
  fBeamEtrigHlt = p.get<std::vector<int>>("BeamEtrigHlt");
  fOffbeamEtrigHlt = p.get<std::vector<int>>("OffbeamEtrigHlt");
  fXmuonEtrigHlt = p.get<std::vector<int>>("XmuonEtrigHlt");
  fBeamCrtT1Hlt = p.get<int>("BeamCrtT1Hlt");
  fOffbeamCrtT1Hlt = p.get<int>("OffbeamCrtT1Hlt");
  fBeamGateHlt = p.get<int>("BeamGateHlt");
  fOffbeamGateHlt = p.get<int>("OffbeamGateHlt");
  fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel");
  fShiftData2MC = p.get<uint64_t>("ShiftData2MC"); //TODO: Define this parameter + Get from database instead of fhicl parameters
  fShiftTdcRwm2PtbGate = p.get<uint64_t>("ShiftTdcRwm2PtbGate"); //TODO: Get from database instead of fhicl parameters
  fShiftTdcEtrig2PtbEtrig = p.get<uint64_t>("ShiftTdcEtrig2PtbEtrig"); //TODO: Get from database instead of fhicl parameters
  fShiftTdcCrtt12PtbCrtt1 = p.get<uint64_t>("ShiftTdcCrtt12PtbCrtt1"); //TODO: Get from database instead of fhicl parameters
  fMakeTree = p.get<bool>("MakeTree", false);
  fDebugDAQHeader = p.get<bool>("DebugDAQHeader", false);
  fDebugPtb = p.get<bool>("DebugPtb", false);
  fDebugTdc = p.get<bool>("DebugTdc", false);
  fDebugFrame = p.get<bool>("DebugFrame", false);
  
  produces<FrameShiftInfo>();
  produces<TimingInfo>();
}

void sbnd::timing::FrameShift::produce(art::Event& e)
{
  ResetEventVars();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if(fDebugTdc | fDebugPtb | fDebugFrame)
    std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";

  GetRawTimestamp(e);
  GetTDCTimestamps(e);
  GetPTBTimestamps(e);
  FindETRIGs();

  uint64_t global_etrig_ts = DecideGlobalEtrigTimestamp();
  DecideRelevantTDCTimestamps(global_etrig_ts);
  DecideRelevantPTBTimestamps(global_etrig_ts);

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
    
  //Frame CRT T1
  if(_isBeam || _isOffbeam)
    {
      if(_tdc_crtt1_ts != kInvalidTimestamp)
        {
          _frame_crtt1          = _tdc_crtt1_ts; //TODO: Add shift from TDC to PTB: +fShiftTdcCrtt12PtbCrtt1
          _timing_type_crtt1    = kSPECTDCType;
          _timing_channel_crtt1 = 0;
        }
      else if(_hlt_crtt1_ts != kInvalidTimestamp)
        {
          _frame_crtt1          = _hlt_crtt1_ts;
          _timing_type_crtt1    = kPTBHLTType;
          _timing_channel_crtt1 = _hlt_crtt1;
        }
      else
        {
          _frame_crtt1          = 0;
          _timing_type_crtt1    = kNoShiftType;
          _timing_channel_crtt1 = kInvalidChannel;
        }

      //Frame Beam Gate
      if(_isBeam && _tdc_rwm_ts != kInvalidTimestamp) // TODO: For Offbeam, I think HLT Gate is recorded in TDC as FTRIG and can be found
        {
          _frame_gate          = _tdc_rwm_ts + fShiftTdcRwm2PtbGate; //TODO: Add  shift from Data to MD: + fShiftData2MC
          _timing_type_gate    = kSPECTDCType;
          _timing_channel_gate = 2;
        }
      else if(_hlt_gate_ts != kInvalidTimestamp)
        {
          _frame_gate          = _hlt_gate_ts;
          _timing_type_gate    = kPTBHLTType;
          _timing_channel_gate = _hlt_gate;
        }
      else
        {
          _frame_gate          = 0;
          _timing_type_gate    = kNoShiftType;
          _timing_channel_gate = kInvalidChannel;
        }
    }

  //Frame ETRIG
  if(_tdc_etrig_ts != kInvalidTimestamp)
    {
      _frame_etrig          = _tdc_etrig_ts + fShiftTdcEtrig2PtbEtrig;
      _timing_type_etrig    = kSPECTDCType;
      _timing_channel_etrig = 4;
    }
  else if(_hlt_etrig_ts != kInvalidTimestamp)
    {
      _frame_etrig          = _hlt_etrig_ts;
      _timing_type_etrig    = kPTBHLTType;
      _timing_channel_etrig = _hlt_etrig;
    }
  else
    {
      _frame_etrig          = 0;
      _timing_type_etrig    = kNoShiftType;
      _timing_channel_etrig = kInvalidChannel;
    }

  if(_isBeam || _isOffbeam)
    {
      //Pick default stream -- beam gate
      _frame_default          = _frame_gate;
      _timing_type_default    = _timing_type_gate;
      _timing_channel_default = _timing_channel_gate;
    }
  else if(_isXmuon)
    {
      //Pick default stream -- ETRIG
      _frame_default        = _frame_etrig;
      _timing_type_default  = _timing_type_etrig;
      _timing_channel_etrig = _timing_channel_etrig;
    }

  if(fDebugFrame)
    {
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "Frame Shift Results:" << std::endl;
      std::cout << "Frame CRT T1 type " << _timing_type_crtt1 <<", channel " << _timing_channel_crtt1 << ": "
                << PrintFormatTimestamp(_frame_crtt1) << std::endl;
      std::cout << "Frame Beam Gate  type " << _timing_type_gate <<", channel " << _timing_channel_gate << ": "
                << PrintFormatTimestamp(_frame_gate) << std::endl;
      std::cout << "Frame ETRIG type " << _timing_type_etrig <<", channel " << _timing_channel_etrig << ": "
                << PrintFormatTimestamp(_frame_etrig) << std::endl;
      std::cout << "Default Frame type " << _timing_type_default <<", channel " << _timing_channel_default << ": "
                << PrintFormatTimestamp(_frame_default) << std::endl;
      std::cout << "--------------------------------------" << std::endl;
    }

  //Put product in event
  std::unique_ptr<FrameShiftInfo> newFrameShiftInfo(new FrameShiftInfo(_frame_crtt1, _timing_type_crtt1, _timing_channel_crtt1,
                                                                       _frame_gate, _timing_type_gate, _timing_channel_gate,
                                                                       _frame_etrig, _timing_type_etrig,_timing_channel_etrig,
                                                                       _frame_default, _timing_type_default, _timing_channel_default));

  std::unique_ptr<TimingInfo> newTimingInfo(new TimingInfo(_raw_ts, _tdc_crtt1_ts, _tdc_bes_ts, _tdc_rwm_ts, _tdc_etrig_ts,
                                                           _hlt_crtt1_ts, _hlt_etrig_ts, _hlt_gate_ts));

  e.put(std::move(newTimingInfo));
  e.put(std::move(newFrameShiftInfo));

  //Fill the tree
  if(fMakeTree)
    fTree->Fill();
}

uint64_t sbnd::timing::FrameShift::FindClosest(const std::vector<uint64_t> &timestamps, const uint64_t &reference)
{
  uint64_t min_diff = kInvalidTimestamp;
  uint64_t closest  = kInvalidTimestamp;

  for(const uint64_t &timestamp : timestamps)
    {
      uint64_t diff = reference > timestamp ? reference - timestamp : timestamp - reference;

      if(diff < min_diff)
        {
          min_diff = diff;
          closest  = timestamp;
        }
    }

  return closest;
}

std::string sbnd::timing::FrameShift::PrintFormatTimestamp(const uint64_t &timestamp)
{
  std::stringstream ss;
  ss << "(s) = " <<  timestamp / kSecondInNanoseconds << ", (ns) = " << timestamp % kSecondInNanoseconds;

  return ss.str();
}

void sbnd::timing::FrameShift::GetRawTimestamp(const art::Event &e)
{
  art::Handle<artdaq::detail::RawEventHeader> DAQHeaderHandle;
  e.getByLabel(fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel, DAQHeaderHandle);

  if(!DAQHeaderHandle.isValid())
    throw cet::exception("FrameShift") << "No artdaq::detail::RawEventHeader found w/ tag " << fDAQHeaderModuleLabel << ". Check data quality!";
  else
    {
      artdaq::RawEvent rawHeaderEvent = artdaq::RawEvent(*DAQHeaderHandle);
      _raw_ts = rawHeaderEvent.timestamp() - fRawTSCorrection;
    }

  if(fDebugDAQHeader)
    {
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << "DAQ Header Timestamp: " << PrintFormatTimestamp(_raw_ts) << std::endl;
    }
}

void sbnd::timing::FrameShift::GetTDCTimestamps(const art::Event &e)
{
  art::Handle<std::vector<DAQTimestamp>> tdcHandle;
  e.getByLabel(fTdcDecodeLabel, tdcHandle);
  
  if(!tdcHandle.isValid() || tdcHandle->size() == 0)
    mf::LogInfo("FrameShift") << "No DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!\n";
  else
    {
      for(auto const& tdc : *tdcHandle)
        {
          const uint32_t  ch = tdc.Channel();
          const uint64_t  ts = tdc.Timestamp();
          //Also TODO: Make use of picosecond timestamps
          //const uint64_t  ts_ps = tdc.TimestampPs();

          switch(ch)
            {
            case 0:
              _tdc_ch0.push_back(ts);
              break;
            case 1:
              _tdc_ch1.push_back(ts);
              break;
            case 2:
              _tdc_ch2.push_back(ts);
              break;
            case 4:
              _tdc_ch4.push_back(ts);
              break;
            }
        }
    }
}

void sbnd::timing::FrameShift::GetPTBTimestamps(const art::Event &e)
{
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  e.getByLabel(fPtbDecodeLabel, ptbHandle);
  
  if((!ptbHandle.isValid() || ptbHandle->size() == 0))
    throw cet::exception("FrameShift") << "No raw::ptb::sbndptb found w/ tag " << fPtbDecodeLabel << ". Check data quality!";
  else
    {
      // HLT Words
      unsigned nHLTs = 0;

      for(auto const& ptb : *ptbHandle)
        nHLTs += ptb.GetNHLTriggers();
  
      _ptb_hlt_trigger.resize(nHLTs);
      _ptb_hlt_timestamp.resize(nHLTs);
      _ptb_hlt_trunmask.resize(nHLTs);
      _ptb_hlt_unmask_timestamp.resize(nHLTs);

      unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
      unsigned h_i = 0; //For trigger with bitmask

      for(auto const& ptb : *ptbHandle)
        {
          for(unsigned i = 0; i < ptb.GetNHLTriggers(); ++i)
            {
              _ptb_hlt_trigger[h_i]   = ptb.GetHLTrigger(i).trigger_word;
              _ptb_hlt_timestamp[h_i] = ptb.GetHLTrigger(i).timestamp * uint64_t(20); //Units can be found in the Decoder Module 
              h_i++;
  
              int val = ptb.GetHLTrigger(i).trigger_word;
              int upBit[32];
    
              for(int u=0; u<32; u++)
                upBit[u] = -1; //setting default values for maximum of 32 bits
  
              int numOfTrig = 0;
              for(int b=0; b<32; b++)
                {
                  if((val & 0x01) == 1)
                    {
                      upBit[numOfTrig] = b;
                      numOfTrig++;
                    }

                  val = val >> 1;
                }
    
              if(numOfTrig == 1)
                {
                  _ptb_hlt_unmask_timestamp[hlt_i] = _ptb_hlt_timestamp[h_i-1];
                  _ptb_hlt_trunmask[hlt_i]         = upBit[0];
                  hlt_i++;
                }
              else if(numOfTrig > 1)
                {
                  nHLTs += (numOfTrig - 1);
                  _ptb_hlt_unmask_timestamp.resize(nHLTs);
                  _ptb_hlt_trunmask.resize(nHLTs);
  
                  for(int mult = 0; mult < numOfTrig; mult++)
                    {
                      _ptb_hlt_unmask_timestamp[hlt_i] = _ptb_hlt_timestamp[h_i-1];
                      _ptb_hlt_trunmask[hlt_i]         = upBit[mult];
                      hlt_i++;
                    }
                }
            } //End of loop over nHLTriggers
        } //End of loop over ptb in ptb_v
    }

  if(fDebugPtb)
    {
      for(size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++)
        {
          std::cout << "----------------------------------------------------" << std::endl;
          std::cout << "HLT " << _ptb_hlt_trunmask[i] << ": " << PrintFormatTimestamp(_ptb_hlt_unmask_timestamp[i]) << std::endl;
        }
    }
}

void sbnd::timing::FrameShift::FindETRIGs()
{
  if(_tdc_ch4.size() == 0)
    throw cet::exception("FrameShift") << "No TDC ETRIG timestamps found! Check data quality!";
  else
    _tdc_etrig_ts = FindClosest(_tdc_ch4, _raw_ts);

  //Grab all the ETRIG HLTS -- there might be more than 1
  for(size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++)
    {
      for(size_t j = 0; j < fPtbEtrigHlts.size(); j++)
        {
          if(_ptb_hlt_trunmask[i] == fPtbEtrigHlts[j])
            {
              _ptb_hlt_etrig.push_back(_ptb_hlt_trunmask[i]);
              _ptb_hlt_etrig_ts.push_back(_ptb_hlt_unmask_timestamp[i]);
            }
        }
    }

  if(_ptb_hlt_etrig.size() == 0)
    throw cet::exception("FrameShift") << "No HLT ETRIG timestamps found! Check data quality!";
  else
    {
      uint64_t min_diff = std::numeric_limits<uint64_t>::max();
      for(size_t i = 0; i < _ptb_hlt_etrig.size(); i++)
        {
          uint64_t diff = _raw_ts > _ptb_hlt_etrig_ts[i] ? _raw_ts - _ptb_hlt_etrig_ts[i] : _ptb_hlt_etrig_ts[i] - _raw_ts;

          if(diff < min_diff)
            {
              min_diff      = diff;
              _hlt_etrig    = _ptb_hlt_etrig[i];
              _hlt_etrig_ts = _ptb_hlt_etrig_ts[i];
            }
        }
    }
}

uint64_t sbnd::timing::FrameShift::DecideGlobalEtrigTimestamp()
{
  //Decide which global frame to use as reference
  // Prioritise TDC ETRIG then PTB ETRIG then the raw header.
  uint64_t global_etrig_ts = kInvalidTimestamp;

  if(_tdc_etrig_ts != kInvalidTimestamp)
    global_etrig_ts = _tdc_etrig_ts;
  else if(_hlt_etrig_ts != kInvalidTimestamp)
    global_etrig_ts = _hlt_etrig_ts;
  else
    global_etrig_ts = _raw_ts;

  if(fDebugFrame)
    {
      std::cout << "----------------------------------------------------" << std::endl;
      if(_tdc_etrig_ts != kInvalidTimestamp)
        std::cout << "Using TDC ETRIG as Global Frame Reference" << std::endl;
      else if(_hlt_etrig_ts != kInvalidTimestamp)
        std::cout << "Using PTB HLT ETRIG as Global Frame Reference" << std::endl;
      else
        std::cout << "Using DAQ Header Timestamp as Global Frame Reference" << std::endl;

      std::cout << "Global Frame Timestamp: " << PrintFormatTimestamp(global_etrig_ts) << std::endl;
    }

  return global_etrig_ts;
}

void sbnd::timing::FrameShift::DecideRelevantTDCTimestamps(const uint64_t &global_etrig_ts)
{
  // ch0: CRT T1
  if(_tdc_ch0.size() != 0)
    _tdc_crtt1_ts = FindClosest(_tdc_ch0, global_etrig_ts);

  // ch1: BES
  if(_tdc_ch1.size() != 0)
    _tdc_bes_ts = FindClosest(_tdc_ch1, global_etrig_ts);

  // ch2: RWM
  if(_tdc_ch2.size() != 0)
    _tdc_rwm_ts = FindClosest(_tdc_ch2, global_etrig_ts);

  if(fDebugTdc)
    {
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << "TDC Channel 0 (CRTT1) Timestamp: " << PrintFormatTimestamp(_tdc_crtt1_ts) << std::endl;
      std::cout << "TDC Channel 1 (BES) Timestamp: " << PrintFormatTimestamp(_tdc_bes_ts) << std::endl;
      std::cout << "TDC Channel 2 (RWM) Timestamp: "<< PrintFormatTimestamp(_tdc_rwm_ts) << std::endl;
      std::cout << "TDC Channel 4 (ETRIG) Timestamp: "<< PrintFormatTimestamp(_tdc_etrig_ts) << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
    }
}

void sbnd::timing::FrameShift::DecideRelevantPTBTimestamps(const uint64_t &global_etrig_ts)
{
  //Check which Gate/CRT T1 HLT to use based on ETRIG HLT
  //Order to check: Beam -> Offbeam -> Xmuon

  for(const int &beam_etrig_hlt : fBeamEtrigHlt)
    {
      if(_hlt_etrig == beam_etrig_hlt)
        {
          _hlt_gate  = fBeamGateHlt;
          _hlt_crtt1 = fBeamCrtT1Hlt;
          _isBeam    = true;
          break;
        }
    }

  if(!_isBeam)
    {
      for(const int &offbeam_etrig_hlt : fOffbeamEtrigHlt)
        {
          if(_hlt_etrig == offbeam_etrig_hlt)
            {
              _hlt_gate  = fOffbeamGateHlt;
              _hlt_crtt1 = fOffbeamCrtT1Hlt;
              _isOffbeam = true;
              break;
            }
        }
    }

  if(!_isBeam & !_isOffbeam)
    {
      for(const int &xmuon_etrig_hlt : fXmuonEtrigHlt)
        {
          if(_hlt_etrig == xmuon_etrig_hlt)
            {
              _isXmuon = true;
              break;
            }
        }
    }

  if( !_isBeam & !_isOffbeam & !_isXmuon)
    throw cet::exception("FrameShift") << "ETRIG HLT " << _hlt_etrig << " does not match any known Beam/Offbeam/Xmuon ETRIG HLT! Check data quality!";

  //Get Gate and CRT T1 HLT timestamps
  //TODO: What if there is no Gate or CRT T1 HLT?
  for(size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++)
    {
      if(_ptb_hlt_trunmask[i] == _hlt_gate)
        _hlt_gate_ts = _ptb_hlt_unmask_timestamp[i];
      if(_ptb_hlt_trunmask[i] == _hlt_crtt1)
        _hlt_crtt1_ts = _ptb_hlt_unmask_timestamp[i];
    }

  if(fDebugPtb)
    {
      std::cout << "----------------------------------------------------" << std::endl;
      if(_isBeam) std::cout << "This is Beam Stream!" << std::endl;
      if(_isOffbeam) std::cout << "This is Offbeam Stream!" << std::endl;
      if(_isXmuon) std::cout << "This is Crossing Muon Stream!" << std::endl;
      std::cout << "HLT ETRIG = " << _hlt_etrig << ", Timestamp: " << PrintFormatTimestamp(_hlt_etrig_ts) << std::endl;
      std::cout << "HLT Gate = " << _hlt_gate << ", Timestamp: " << PrintFormatTimestamp(_hlt_gate_ts) << std::endl;
      std::cout << "HLT CRT T1 = " << _hlt_crtt1 << ", Timestamp: " << PrintFormatTimestamp(_hlt_crtt1_ts) << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
    }
}

void sbnd::timing::FrameShift::beginJob()
{
  if(fMakeTree)
    {
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
}

void sbnd::timing::FrameShift::ResetEventVars()
{
  _run = -1;
  _subrun = -1;
  _event = -1;

  _raw_ts = kInvalidTimestamp;

  _tdc_ch0.clear();
  _tdc_ch1.clear();
  _tdc_ch2.clear();
  _tdc_ch4.clear();

  _tdc_crtt1_ts = kInvalidTimestamp;
  _tdc_bes_ts = kInvalidTimestamp;
  _tdc_rwm_ts = kInvalidTimestamp;
  _tdc_etrig_ts = kInvalidTimestamp;

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
  _hlt_etrig_ts = kInvalidTimestamp;
  _hlt_gate = std::numeric_limits<int>::max();
  _hlt_gate_ts = kInvalidTimestamp;
  _hlt_crtt1= std::numeric_limits<int>::max();
  _hlt_crtt1_ts = kInvalidTimestamp;
  
  _frame_crtt1 = kInvalidTimestamp;
  _timing_type_crtt1 = kInvalidType;
  _timing_channel_crtt1 = kInvalidChannel;
  
  _frame_gate = kInvalidTimestamp;
  _timing_type_gate = kInvalidType;
  _timing_channel_gate = kInvalidChannel;
  
  _frame_etrig = kInvalidTimestamp;
  _timing_type_etrig = kInvalidType;
  _timing_channel_etrig = kInvalidChannel;

  _frame_default = kInvalidTimestamp;
  _timing_type_default = kInvalidType;
  _timing_channel_default = kInvalidChannel;
}

DEFINE_ART_MODULE(sbnd::timing::FrameShift)
