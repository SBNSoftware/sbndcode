////////////////////////////////////////////////////////////////////////
// Class:       TimeStampDumper
// Plugin Type: analyzer (Unknown Unknown)
// File:        TimeStampDumper_module.cc
//
// Generated at Tue Jan  7 16:11:21 2025 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"


#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "artdaq-core/Data/RawEvent.hh"


#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"


#include "TTree.h"

class TimeStampDumper;


class TimeStampDumper : public art::EDAnalyzer {
public:
  explicit TimeStampDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TimeStampDumper(TimeStampDumper const&) = delete;
  TimeStampDumper(TimeStampDumper&&) = delete;
  TimeStampDumper& operator=(TimeStampDumper const&) = delete;
  TimeStampDumper& operator=(TimeStampDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void ResetTree();

private:

  // Declare member data here.
  std::string fPTBLabel;
  std::string fPMTWaveformLabel;
  int fTotalCAENBoards;
  int fPMTPerCAEN;
  uint32_t                  fraw_ts_correction;
  std::string               fspectdc_product_name;
  uint32_t                  fspectdc_ftrig_ch;
  uint32_t                  fspectdc_etrig_ch;
  TTree* evtTree;
  int NLLT;
  int NHLT;
  int tree_event;
  int tree_run;
  int tree_subrun;
  int EventTime_s;
  int EventTime_ns;
  int TriggerPulseID;
  std::vector<int> LLT_Type;
  std::vector<int> HLT_Type;
  std::vector<long> LLT_Time;
  std::vector<long> HLT_Time;
  std::vector<long> Flash_Time;

};

void TimeStampDumper::ResetTree()
  {
    tree_event = -1;
    tree_run = -1;
    tree_subrun = -1;
    EventTime_s = -9999;
    EventTime_ns= -9999;
    TriggerPulseID=-1;
    NLLT=-1;
    NHLT=-1;
    int LLTSize=600;
    LLT_Type.clear();
    LLT_Type.resize(LLTSize);
    LLT_Time.clear();
    LLT_Time.resize(LLTSize);
    int HLTSize=25;
    HLT_Type.clear();
    HLT_Type.resize(HLTSize);
    HLT_Time.clear();
    HLT_Time.resize(HLTSize);
    int NFlashMax = 50;
    Flash_Time.clear();
    Flash_Time.resize(NFlashMax);
  }



TimeStampDumper::TimeStampDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fPMTWaveformLabel = p.get< std::string >("PMTLabel",  "pmtdecoder:PMTChannels:DECODE");
  fPTBLabel = p.get< std::string >("PTBLabel",  "ptbdecoder::DECODE");
  fTotalCAENBoards = p.get<int>("TotalCAENBoards", 8);
  fPMTPerCAEN = p.get<int>("PMTPerCAEN", 15);
  fspectdc_product_name = p.get<std::string>("spectdc_product_name","tdcdecoder");
  fspectdc_ftrig_ch = p.get<uint32_t>("spectdc_ftrig_ch",3);
  fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch",4);
  fraw_ts_correction = p.get<uint>("raw_ts_correction",367000); // ns
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  ResetTree();
  art::ServiceHandle<art::TFileService> tfs;
  evtTree = tfs->make<TTree>("TimeStamps","Time Stamp Dump Tree");
  evtTree->Branch("event",&tree_event);
  evtTree->Branch("run",&tree_run);
  evtTree->Branch("subrun",&tree_subrun);
  evtTree->Branch("NLLT",&NLLT);
  evtTree->Branch("NHLT",&NHLT);
  evtTree->Branch("EventTime_s",&EventTime_s);
  evtTree->Branch("EventTime_ns",&EventTime_ns);
  evtTree->Branch("TriggerPulseID",&TriggerPulseID);
  evtTree->Branch("LLT_Type",&LLT_Type);
  evtTree->Branch("HLT_Type",&HLT_Type);
  evtTree->Branch("LLT_Time",&LLT_Time);
  evtTree->Branch("HLT_Time",&HLT_Time);
  evtTree->Branch("Flash_Time",&Flash_Time);
}

void TimeStampDumper::analyze(art::Event const& e)
{
  ResetTree();
  tree_run = e.run();
  tree_subrun = e.subRun();
  tree_event =  e.id().event();;
  // Implementation of required member function here.
  EventTime_s=e.time().timeHigh();
  EventTime_ns=e.time().timeLow();

  //PTB Trigger times
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle_2;
  e.getByLabel(fPTBLabel,ptbHandle_2);
  for(int index=0; index<int(ptbHandle_2->size()); index++)
  {
    auto ptb = (*ptbHandle_2)[index];
    auto hltrigs = ptb.GetHLTriggers();
    auto lltrigs = ptb.GetLLTriggers();
    NHLT = NHLT+int(hltrigs.size());
    NLLT = NLLT+int(lltrigs.size());
    for(int HLT=0; HLT<int(hltrigs.size()); HLT++)
    {
      int Power=0;
      while(Power<64)
      {
        if(hltrigs[HLT].trigger_word & (0x1 << Power)) break;
        else Power=Power+1;
      }
      HLT_Type[HLT] = Power;
      HLT_Time[HLT] = hltrigs[HLT].timestamp*20;
    }
    for(int LLT=0; LLT<int(lltrigs.size()); LLT++)
    {
      int Power=0;
      while(Power<64)
      {
        if(lltrigs[LLT].trigger_word & (0x1 << Power)) break;
        else Power=Power+1;
      }
      LLT_Type[LLT] = Power;
      LLT_Time[LLT] = lltrigs[LLT].timestamp*20;
    }
  }
  //PMT timestamps
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
  e.getByLabel(fPMTWaveformLabel, waveHandle);
  int PMTPerCAEN=fPMTPerCAEN;
  int TotalFlash = waveHandle->size()/(fTotalCAENBoards*PMTPerCAEN);
  std::cout << "\t It also has " << TotalFlash << " opDetWaveforms " << std::endl;
  long MinTime = 9999999999;
  for(int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
  {
    int WaveIndex = FlashCounter*PMTPerCAEN;
    double currentTimeStamp = (*waveHandle)[WaveIndex].TimeStamp();
    Flash_Time[FlashCounter] = currentTimeStamp;
    if(TMath::Abs(currentTimeStamp)<MinTime)
    {
      MinTime=currentTimeStamp;
      TriggerPulseID = FlashCounter;
    }
  }
//TDC timestamps
art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
e.getByLabel(fspectdc_product_name,tdcHandle);
bool found_ett = false;
std::vector<uint64_t> tdc_etrig_v;
uint64_t min_raw_tdc_diff = uint64_t(1e12);
const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
art::Handle<artdaq::detail::RawEventHeader> header_handle;
uint64_t raw_timestamp = 0;
e.getByLabel("daq", "RawEventHeader", header_handle);
auto rawheader = artdaq::RawEvent(*header_handle); 
raw_timestamp = rawheader.timestamp() - fraw_ts_correction; // includes sec + ns portion
std::cout << "Raw timestamp (w/ correction) -> "  << "ts (ns): " << raw_timestamp % uint64_t(1e9) << ", sec (s): " << raw_timestamp / uint64_t(1e9) << std::endl;
std::cout << "vs event object itself " << EventTime_s << " sec " << EventTime_ns << std::endl;
for (size_t i=0; i<tdc_v.size(); i++){
    auto tdc = tdc_v[i];
    const uint32_t  ch = tdc.Channel();
    const uint64_t  ts = tdc.Timestamp();
    const uint64_t  offset = tdc.Offset();
    const std::string name  = tdc.Name();

    if (true){
        std::cout << "      TDC CH " << ch << " -> "
        << "name: " << name
        << ", ts (ns): " << ts%uint64_t(1e9)
        << ", sec (s): " << ts/uint64_t(1e9)
        << ", offset: " << offset 
        << std::endl;
    }
    if (ch==fspectdc_etrig_ch){
        found_ett = true;
        tdc_etrig_v.push_back(ts);
    }
}
uint64_t event_trigger_time = 0; // in ns
if (tdc_etrig_v.size()==1)
  {
      event_trigger_time = tdc_etrig_v.front()%uint64_t(1e9);
  }
else
  { // finding the closest ETRIG to the raw timestamp
      for (size_t i=0; i < dc_etrig_v.size(); i++){
          auto tdc_etrig = tdc_etrig_v[i];
          uint64_t diff;
          if (tdc_etrig < (raw_timestamp))
              diff = raw_timestamp - tdc_etrig;
          else
              diff = tdc_etrig - raw_timestamp;
          if (diff < min_raw_tdc_diff){
              event_trigger_time = tdc_etrig%uint64_t(1e9);
              min_raw_tdc_diff = diff;
          }
      }
  }


evtTree->Fill();
}

DEFINE_ART_MODULE(TimeStampDumper)
