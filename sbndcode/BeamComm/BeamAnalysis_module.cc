////////////////////////////////////////////////////////////////////////
// Class:       BeamAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        BeamAnalysis_module.cc
//
// Generated at Tue Oct 29 09:11:32 2024 by Vu Chi Lan Nguyen using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "art_root_io/TFileService.h"

#include "TTree.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndcode/Decoders/PMT/sbndpmt.h"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

namespace sbnd {
  class BeamAnalysis;
}


class sbnd::BeamAnalysis : public art::EDAnalyzer {
public:
    explicit BeamAnalysis(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    BeamAnalysis(BeamAnalysis const&) = delete;
    BeamAnalysis(BeamAnalysis&&) = delete;
    BeamAnalysis& operator=(BeamAnalysis const&) = delete;
    BeamAnalysis& operator=(BeamAnalysis&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    void ResetEventVars();

private:

    // Event Tree
    TTree *_tree;

    //---Tree Variables
    int _run, _subrun, _event;

    // TDC stuff
    std::vector<uint64_t> tdc_ch0;
    std::vector<uint64_t> tdc_ch1;
    std::vector<uint64_t> tdc_ch2;
    std::vector<uint64_t> tdc_ch3;
    std::vector<uint64_t> tdc_ch4;

    std::vector<uint64_t> tdc_ch0_utc;
    std::vector<uint64_t> tdc_ch1_utc;
    std::vector<uint64_t> tdc_ch2_utc;
    std::vector<uint64_t> tdc_ch3_utc;
    std::vector<uint64_t> tdc_ch4_utc;

    // CRT stuff
    std::vector<double> crt_x;
    std::vector<double> crt_y;
    std::vector<double> crt_z;
    std::vector<double> crt_xe;
    std::vector<double> crt_ye;
    std::vector<double> crt_ze;
    std::vector<double> crt_ts0;
    std::vector<double> crt_ts1;
    std::vector<double> crt_ts0e;
    std::vector<double> crt_ts1e;
    std::vector<int> crt_tagger;

    // PMT Timing
    uint16_t pmt_timing_type;
    uint16_t pmt_timing_ch;

    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fCrtSpacePointLabel;
    art::InputTag fPmtFtrigDecodeLabel;
    art::InputTag fPmtFtrigBoardLabel;
    art::InputTag fPmtTimingLabel;

    // Debug
    bool fDebugTdc;
    bool fDebugCrt;
    bool fDebugPmt;
    
    // Which data products to include
    bool fIncludeCrt;
    bool fIncludePmt;
};


sbnd::BeamAnalysis::BeamAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // 
    // More initializers here.
{
    fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");

    fCrtSpacePointLabel = p.get<art::InputTag>("CrtSpacePointLabel", "crtspacepoints");

    fPmtTimingLabel = p.get<art::InputTag>("fPmtTimingLabel", "pmtdecoder");

    fPmtFtrigDecodeLabel = p.get<art::InputTag>("PmtFtrigDecodeLabel", "pmtdecoder:FTrigChannels");
    fPmtFtrigBoardLabel = p.get<art::InputTag>("PmtFtrigBoardLabel", "pmtdecoder:FTrigTiming");

    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugCrt = p.get<bool>("DebugCrt", false);
    fDebugPmt = p.get<bool>("DebugPmt", false);

    fIncludeCrt = p.get<bool>("IncludeCrt", true);
    fIncludePmt = p.get<bool>("IncludePmt", true);
}

void sbnd::BeamAnalysis::analyze(art::Event const& e)
{
    
    ResetEventVars();

    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  =  e.id().event();

    if (fDebugTdc | fDebugCrt | fDebugPmt)
        std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";

    //---------------------------TDC-----------------------------//
    art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
    e.getByLabel(fTdcDecodeLabel, tdcHandle);

    if (!tdcHandle.isValid() || tdcHandle->size() == 0){
        if (fDebugTdc) std::cout << "No SPECTDC products found. Skip this event." << std::endl;
        return;
    }
    else{
        
        const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
        
        for (size_t i=0; i<tdc_v.size(); i++){
            auto tdc = tdc_v[i];
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();
            //const uint64_t  offset = tdc.Offset();
            const std::string name  = tdc.Name();

            if (fDebugTdc) std::cout << "Ch " << name
                                << ", ts (ns) = " << ts%uint64_t(1e9)
                                << ", sec (s) = " << ts/uint64_t(1e9)
                                << std::endl;

            if(ch == 0){
                tdc_ch0_utc.push_back(ts);
                tdc_ch0.push_back(ts%uint64_t(1e9));
            }
            if(ch == 1){
                tdc_ch1_utc.push_back(ts);
                tdc_ch1.push_back(ts%uint64_t(1e9));
            }
            if(ch == 2){
                tdc_ch2_utc.push_back(ts);
                tdc_ch2.push_back(ts%uint64_t(1e9));
            }
            if(ch == 3){
                tdc_ch3_utc.push_back(ts);
                tdc_ch3.push_back(ts%uint64_t(1e9));
            }
            if(ch == 4){
                tdc_ch4_utc.push_back(ts);
                tdc_ch4.push_back(ts%uint64_t(1e9));
            }
        } 
    }
    //------------------------------------------------------------//

    //---------------------------CRT-----------------------------//
   
    if (fIncludeCrt){

        art::Handle<std::vector<sbnd::crt::CRTSpacePoint>> crtSpacePointHandle;
        std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_sp_v;
        e.getByLabel(fCrtSpacePointLabel, crtSpacePointHandle);

        if (!crtSpacePointHandle.isValid() || crtSpacePointHandle->size() == 0){
            if (fDebugCrt) std::cout << "No CRT Space Point products found." << std::endl;
        }
        else{

            art::fill_ptr_vector(crt_sp_v, crtSpacePointHandle);
            art::FindManyP<sbnd::crt::CRTCluster> crtSPClusterAssoc(crt_sp_v, e, fCrtSpacePointLabel);

            for (auto const& crt_sp: crt_sp_v){
                
                if (!crt_sp->Complete()) continue;

                const std::vector<art::Ptr<sbnd::crt::CRTCluster>> crt_cluster_v(crtSPClusterAssoc.at(crt_sp.key()));

                if(crt_cluster_v.size() != 1 ) continue;

                const art::Ptr<sbnd::crt::CRTCluster>& crt_cluster(crt_cluster_v.front());

                crt_x.push_back(crt_sp->X());
                crt_y.push_back(crt_sp->Y());
                crt_z.push_back(crt_sp->Z());
                crt_xe.push_back(crt_sp->XErr());
                crt_ye.push_back(crt_sp->YErr());
                crt_ze.push_back(crt_sp->ZErr());
                crt_ts0.push_back(crt_sp->Ts0());
                crt_ts1.push_back(crt_sp->Ts1());
                crt_ts0e.push_back(crt_sp->Ts0Err());
                crt_ts1e.push_back(crt_sp->Ts1Err());
                crt_tagger.push_back(crt_cluster->Tagger());

                if (fDebugCrt){
                    std::cout << "CRT Space Point------------------------------------" << std::endl;
                    std::cout << "   x = " << crt_x.back() << ", y = " << crt_y.back() << ", z = " << crt_z.back() << std::endl;
                    std::cout << "   ts0 = " << crt_ts0.back() << ", ts1 = " << crt_ts1.back() << std::endl;
                    std::cout << "   tagger = " << crt_tagger.back() << std::endl;
                }
            }
        }
    }

    //-----------------------------------------------------------------//

    if (fIncludePmt){

        //------------------------PMT Timing--------------------------//
        art::Handle<raw::pmt::eventTimingInfo> pmtTimingHandle;
        e.getByLabel(fPmtTimingLabel, pmtTimingHandle);

        if (!pmtTimingHandle.isValid()){
            if (fDebugPmt) std::cout << "No PMT Timing products found." << std::endl;
        }
        else{
            raw::pmt::eventTimingInfo const& pmt_timing(*pmtTimingHandle);
            
            pmt_timing_type = pmt_timing.timingType;
            pmt_timing_ch = pmt_timing.timingChannel;

            if (fDebugPmt){
                std::cout << "Timing Reference For Decoding PMT" << std::endl;
                std::cout << "   Type = " << pmt_timing_type << " (SPECTDC = 0; PTB HLT = 1; CAEN-only = 3)." << std::endl;
                std::cout << "   Channel = " << pmt_timing_ch << " (TDC ETRIG = 4; PTB BNB Beam+Light = 2)." << std::endl;
            }
        }

        //---------------------------PMT-----------------------------//
        art::Handle<std::vector<raw::OpDetWaveform>> pmtFtrigHandle;
        std::vector<art::Ptr<raw::OpDetWaveform>> pmt_ftrig_v;
        e.getByLabel(fPmtFtrigDecodeLabel, pmtFtrigHandle);
        
        if (!pmtFtrigHandle.isValid() || pmtFtrigHandle->size() == 0){
            if (fDebugPmt) std::cout << "No PMT Decoder FTRIG products found." << std::endl;
        }
        else{

            art::fill_ptr_vector(pmt_ftrig_v, pmtFtrigHandle);
            art::FindManyP<raw::pmt::boardTimingInfo> pmtBoardAssoc(pmt_ftrig_v, e, fPmtFtrigBoardLabel);

            std::cout << "Found OpDetWaveform FTRIG size = " << pmt_ftrig_v.size() << std::endl;

            for (auto const& pmt: pmt_ftrig_v){

                const std::vector<art::Ptr<raw::pmt::boardTimingInfo>> pmt_board_v(pmtBoardAssoc.at(pmt.key()));

                if(pmt_board_v.size() != 1 ) continue;

                const art::Ptr<raw::pmt::boardTimingInfo> pmt_board(pmt_board_v.front());


                std::cout << "   channel id = " << pmt->ChannelNumber() << std::endl;
                std::cout << "   board postpercent = " << pmt_board->postPercent << std::endl;
                std::cout << "   board ttt = " <<  pmt_board->triggerTimeTag.front() << std::endl;

            }
        }
    }

    //-----------------------------------------------------------//
    
    if (fDebugTdc | fDebugCrt | fDebugPmt) std::cout <<"#--------------------------------------------------------#" << std::endl;

    //Fill once every event
    _tree->Fill();
}

void sbnd::BeamAnalysis::beginJob()
{
 
    art::ServiceHandle<art::TFileService> tfs;

    //Event Tree
    _tree = tfs->make<TTree>("events", "");
  
    _tree->Branch("run", &_run);
    _tree->Branch("subrun", &_subrun);
    _tree->Branch("event", &_event);
    _tree->Branch("tdc_ch0", &tdc_ch0);
    _tree->Branch("tdc_ch1", &tdc_ch1);
    _tree->Branch("tdc_ch2", &tdc_ch2);
    _tree->Branch("tdc_ch3", &tdc_ch3);
    _tree->Branch("tdc_ch4", &tdc_ch4);
    _tree->Branch("tdc_ch0_utc", &tdc_ch0_utc);
    _tree->Branch("tdc_ch1_utc", &tdc_ch1_utc);
    _tree->Branch("tdc_ch2_utc", &tdc_ch2_utc);
    _tree->Branch("tdc_ch3_utc", &tdc_ch3_utc);
    _tree->Branch("tdc_ch4_utc", &tdc_ch4_utc);

    _tree->Branch("crt_x", &crt_x);
    _tree->Branch("crt_y", &crt_y);
    _tree->Branch("crt_z", &crt_z);
    _tree->Branch("crt_xe", &crt_xe);
    _tree->Branch("crt_ye", &crt_ye);
    _tree->Branch("crt_ze", &crt_ze);
    _tree->Branch("crt_ts0", &crt_ts0);
    _tree->Branch("crt_ts1", &crt_ts1);
    _tree->Branch("crt_ts0e", &crt_ts0e);
    _tree->Branch("crt_ts1e", &crt_ts1e);
    _tree->Branch("crt_tagger", &crt_tagger);

    _tree->Branch("pmt_timing_type", &pmt_timing_type);
    _tree->Branch("pmt_timing_ch", &pmt_timing_ch);
}

void sbnd::BeamAnalysis::endJob()
{
  // Implementation of optional member function here.
}

void sbnd::BeamAnalysis::ResetEventVars()
{
    _run = -1; _subrun = -1; _event = -1;
   
    tdc_ch0.clear();
    tdc_ch1.clear();
    tdc_ch2.clear();
    tdc_ch3.clear();
    tdc_ch4.clear();

    tdc_ch0_utc.clear();
    tdc_ch1_utc.clear();
    tdc_ch2_utc.clear();
    tdc_ch3_utc.clear();
    tdc_ch4_utc.clear();

    crt_x.clear();
    crt_y.clear();
    crt_z.clear();
    crt_xe.clear();
    crt_ye.clear();
    crt_ze.clear();
    crt_ts0.clear();
    crt_ts1.clear();
    crt_ts0e.clear();
    crt_ts1e.clear();
    crt_tagger.clear();

    pmt_timing_type = -1;
    pmt_timing_ch = -1;
}

DEFINE_ART_MODULE(sbnd::BeamAnalysis)
