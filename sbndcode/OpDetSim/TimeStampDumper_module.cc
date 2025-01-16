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
  std::vector<int> TDC_Channel;
  std::vector<std::string> TDC_Name;
  std::vector<uint64_t> TDC_TimeStamp;

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
    int LLTSize=2400;
    LLT_Type.clear();
    LLT_Type.resize(LLTSize);
    LLT_Time.clear();
    LLT_Time.resize(LLTSize);
    int HLTSize=200;
    HLT_Type.clear();
    HLT_Type.resize(HLTSize);
    HLT_Time.clear();
    HLT_Time.resize(HLTSize);
    int NFlashMax = 200;
    Flash_Time.clear();
    Flash_Time.resize(NFlashMax);
    int TDCSize= 200;
    TDC_Channel.clear();
    TDC_Channel.resize(TDCSize);
    TDC_Name.clear();
    TDC_Name.resize(TDCSize);
    TDC_TimeStamp.clear();
    TDC_TimeStamp.resize(TDCSize);

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
  evtTree->Branch("TDC_Channel",&TDC_Channel);
  evtTree->Branch("TDC_Time",&TDC_TimeStamp);
  evtTree->Branch("TDC_Name",&TDC_Name);
}

void TimeStampDumper::analyze(art::Event const& e)
{
  ResetTree();
  tree_run = e.run();
  tree_subrun = e.subRun();
  tree_event =  e.id().event();;
  // Implementation of required member function here.
  EventTime_s=e.time().timeHigh(); //DONT USE THIS USE EVENT HEADDER
  EventTime_ns=e.time().timeLow();

  //PTB Trigger times
  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle_2;
  e.getByLabel(fPTBLabel,ptbHandle_2);
  int HLTIndexer=0;
  int LLTIndexer=0;
  for(int index=0; index<int(ptbHandle_2->size()); index++)
  {
    auto ptb = (*ptbHandle_2)[index];
    auto hltrigs = ptb.GetHLTriggers();
    auto lltrigs = ptb.GetLLTriggers();
    for(int HLT=0; HLT<int(hltrigs.size()); HLT++)
    {
      int Power=0;
      std::vector<int> GrabbedHLTs;
      while(Power<32)
      {
        if(hltrigs[HLT].trigger_word & (0x1 << Power))
        {
          GrabbedHLTs.push_back(Power);
        }
        Power=Power+1;
      }
      //Do a loop over vector elements of HLT types
      for(int i=0; i<int(GrabbedHLTs.size()); i++ )
      {
        if(HLTIndexer>=int(HLT_Type.size())) return; //Error handling for weird events
        HLT_Type[HLTIndexer] = GrabbedHLTs[i];
        HLT_Time[HLTIndexer] = hltrigs[HLT].timestamp*20;
        HLTIndexer=HLTIndexer+1;
      }
    }
    for(int LLT=0; LLT<int(lltrigs.size()); LLT++)
    {
      int Power=0;
      std::vector<int> GrabbedLLTs;
      while(Power<32)
      {
        if(lltrigs[LLT].trigger_word & (0x1 << Power))
        {
          GrabbedLLTs.push_back(Power);
        }
        Power=Power+1;
      }
      for(int i=0; i<int(GrabbedLLTs.size()); i++ )
      {
        if(LLTIndexer>=int(LLT_Type.size())) return; //Error handling for weird events
        LLT_Type[LLTIndexer] = GrabbedLLTs[i];
        LLT_Time[LLTIndexer] = lltrigs[LLT].timestamp*20;
        LLTIndexer=LLTIndexer+1;
      }
    }
  }
  NHLT = HLTIndexer+1;
  NLLT = LLTIndexer+1;
  if(NHLT>100 || NLLT >2000) std::cout << NHLT << "  " << NLLT << "   something overflowed " << std::endl; 
  //PMT timestamps
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;
  e.getByLabel(fPMTWaveformLabel, waveHandle);
  int PMTPerCAEN=fPMTPerCAEN;
  int TotalFlash = waveHandle->size()/(fTotalCAENBoards*PMTPerCAEN);
  long MinTime = 9999999999;
  if( TotalFlash>int(Flash_Time.size()) ) return; //Error handling for weird events
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
std::vector<uint64_t> tdc_etrig_v;
const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
art::Handle<artdaq::detail::RawEventHeader> header_handle;
uint64_t raw_timestamp = 0;
e.getByLabel("daq", "RawEventHeader", header_handle);
auto rawheader = artdaq::RawEvent(*header_handle); 
raw_timestamp = rawheader.timestamp() - fraw_ts_correction; // includes sec + ns portion
EventTime_s = raw_timestamp/ uint64_t(1e9);
EventTime_ns = raw_timestamp % uint64_t(1e9);
if( tdc_v.size()> TDC_TimeStamp.size() ) return; //Error handling for weird events
for (size_t i=0; i<tdc_v.size(); i++){
    auto tdc = tdc_v[i];
    const uint32_t  ch = tdc.Channel();
    const uint64_t  ts = tdc.Timestamp();
    const uint64_t  offset = tdc.Offset();
    const std::string name  = tdc.Name();
    TDC_TimeStamp[i] = ts;
    TDC_Name[i]=name;
    TDC_Channel[i]=ch;
    if (ch==fspectdc_etrig_ch){
        tdc_etrig_v.push_back(ts-offset);
    }
}
//Save to tree

evtTree->Fill();
}

DEFINE_ART_MODULE(TimeStampDumper)