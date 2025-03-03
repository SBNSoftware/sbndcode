////////////////////////////////////////////////////////////////////////
// Class:       TimingCAEN
// Plugin Type: analyzer (Unknown Unknown)
// File:        TimingCAEN_module.cc
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

#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TSystem.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbndcode/ChannelMaps/PMT/PMTChannelMapService.h"

#include "sbndcode/Decoders/PTB/sbndptb.h"

#include <set>
#include <vector>


namespace sbnd {
  class TimingCAEN;
}


class sbnd::TimingCAEN : public art::EDAnalyzer {
public:
    explicit TimingCAEN(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TimingCAEN(TimingCAEN const&) = delete;
    TimingCAEN(TimingCAEN&&) = delete;
    TimingCAEN& operator=(TimingCAEN const&) = delete;
    TimingCAEN& operator=(TimingCAEN&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    void ResetEventVars();

private:

    // Event Tree
    TTree *fTree;
    std::stringstream _histName; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;

    std::vector<std::string> crtTagger{"Bottom","South","North","West","East","TopLow","TopHigh"};
    size_t nCrt = crtTagger.size(); 
    TH1D* _hTopHatCRTT0[8];
    
    std::vector<std::string> pmtBoard{"0","1","2","3","4","5","6","7","8"};
    size_t nPmt = pmtBoard.size(); 

    //---Tree Variables
    int _run, _subrun, _event;

    // TDC stuff
    std::vector<uint64_t> _tdc_ch0;
    std::vector<uint64_t> _tdc_ch1;
    std::vector<uint64_t> _tdc_ch2;
    std::vector<uint64_t> _tdc_ch3;
    std::vector<uint64_t> _tdc_ch4;

    std::vector<uint64_t> _tdc_ch0_utc;
    std::vector<uint64_t> _tdc_ch1_utc;
    std::vector<uint64_t> _tdc_ch2_utc;
    std::vector<uint64_t> _tdc_ch3_utc;
    std::vector<uint64_t> _tdc_ch4_utc;

    // Ptb stuff
    std::vector<uint64_t> _ptb_hlt_trigger;
    std::vector<uint64_t> _ptb_hlt_timestamp;
    std::vector<int> _ptb_hlt_trunmask;
      
    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;

    //Service
    art::ServiceHandle<SBND::PMTChannelMapService> fPMTChannelMapService;
    
    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fPtbDecodeLabel;
    art::InputTag fTimingRefLabel;
    art::InputTag fTimingDecodeLabel;
    art::InputTag fTimingBoardLabel;

    // Include
    bool fIncludePtb;
    bool fIncludePmt;

    // Debug
    bool fDebugTdc;
    bool fDebugPtb;
    bool fDebugTiming;
};


sbnd::TimingCAEN::TimingCAEN(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // 
    // More initializers here.
{
    fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");
    fTimingRefLabel = p.get<art::InputTag>("fTimingRefLabel", "pmtdecoder");
    fTimingDecodeLabel = p.get<art::InputTag>("fTimingDecodeLabel", "pmtdecoder:TimingChannels");
    fTimingBoardLabel = p.get<art::InputTag>("fTimingBoardLabel", "pmtdecoder:TimingTiming");

    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugPtb = p.get<bool>("DebugPtb", false);
    fDebugTiming = p.get<bool>("DebugTiming", false);

    fIncludePtb = p.get<bool>("IncludePtb", true);
    fIncludePmt = p.get<bool>("IncludePmt", true);
}

void sbnd::TimingCAEN::analyze(art::Event const& e)
{
    
    ResetEventVars();
  
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  =  e.id().event();

    if (fDebugTdc | fDebugPtb| fDebugTiming)
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

            if (fDebugTdc) std::cout << "Ch " << ch << " " << name
                                << ", ts (ns) = " << ts%uint64_t(1e9)
                                << ", sec (s) = " << ts/uint64_t(1e9)
                                << std::endl;

            if(ch == 0){
                _tdc_ch0_utc.push_back(ts);
                _tdc_ch0.push_back(ts%uint64_t(1e9));
            }
            if(ch == 1){
                _tdc_ch1_utc.push_back(ts);
                _tdc_ch1.push_back(ts%uint64_t(1e9));
            }
            if(ch == 2){
                _tdc_ch2_utc.push_back(ts);
                _tdc_ch2.push_back(ts%uint64_t(1e9));
            }
            if(ch == 3){
                _tdc_ch3_utc.push_back(ts);
                _tdc_ch3.push_back(ts%uint64_t(1e9));
            }
            if(ch == 4){
                _tdc_ch4_utc.push_back(ts);
                _tdc_ch4.push_back(ts%uint64_t(1e9));
            }
        } 
    }
    //------------------------------------------------------------//
    if (fIncludePtb){
        art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
        std::vector<art::Ptr<raw::ptb::sbndptb>> ptb_v;
        e.getByLabel(fPtbDecodeLabel,ptbHandle);
                    
        if ((!ptbHandle.isValid() || ptbHandle->size() == 0)){
            if (fDebugPtb) std::cout << "No Ptb products found." << std::endl;
        }
        else{
            art::fill_ptr_vector(ptb_v, ptbHandle);
            
            unsigned nHLTs = 0;
            for(auto const& ptb : ptb_v)
                nHLTs += ptb->GetNHLTriggers();
   
            if (fDebugPtb) std::cout << "--------------Found nHLTs = " << nHLTs << std::endl;

            _ptb_hlt_trigger.resize(nHLTs);
            _ptb_hlt_timestamp.resize(nHLTs);
            _ptb_hlt_trunmask.resize(nHLTs);

            unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
            unsigned h_i = 0; //For trigger with bitmask
            for(auto const& ptb : ptb_v){
                for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i){

                    _ptb_hlt_trigger[h_i] = ptb->GetHLTrigger(i).trigger_word;
                    _ptb_hlt_timestamp[h_i] = ptb->GetHLTrigger(i).timestamp * 20; //Units can be found in the Decoder Module 
                    
                    h_i++;

                    int val = ptb->GetHLTrigger(i).trigger_word;
                    int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
                    int numOfTrig =0;
                    for(int i=0; i<32;i++){
                        if ((val & 0x01) ==1){
                            upBit[numOfTrig] = i;
                            numOfTrig++;
                        }
                        val = val >> 1;
                    }
                    
                    if (numOfTrig ==1){
                        _ptb_hlt_trunmask[hlt_i] = upBit[0];
                        hlt_i++;
                    }//End of if statement for single upbit
                    
                    else if (numOfTrig > 1){
                        nHLTs += (numOfTrig -1);
                        _ptb_hlt_timestamp.resize(nHLTs);
                        _ptb_hlt_trunmask.resize(nHLTs);

                        for (int mult =0; mult < numOfTrig; mult++){ 
                          _ptb_hlt_trunmask[hlt_i] = upBit[mult];
                          hlt_i++;
                        } //End of loop over multiple upbits
                    } //End of else statement for multiple triggers
                    

                } //End of loop over nHLTriggers
            } //End of loop over ptb in PtbVec
            
            if (fDebugPtb){
                for (size_t i = 0; i <_ptb_hlt_trunmask.size(); i++){
                    std::cout << " HLT " << _ptb_hlt_trunmask[i]
                            << ", ts = " << _ptb_hlt_timestamp[i] << std::endl;
                }
            }
        }
    }

    //-----------------------------------------------------------------//

    if (fIncludePmt){

        //------------------------PMT Timing--------------------------//
        art::Handle<raw::TimingReferenceInfo> timingRefHandle;
        e.getByLabel(fTimingRefLabel, timingRefHandle);

        if (!timingRefHandle.isValid()){
            if (fDebugTiming) std::cout << "No Timing Reference products found." << std::endl;
        }
        else{
            raw::TimingReferenceInfo const& pmt_timing(*timingRefHandle);
            
            _pmt_timing_type = pmt_timing.timingType;
            _pmt_timing_ch = pmt_timing.timingChannel;

            if (fDebugTiming){
                std::cout << "Timing Reference For Decoding PMT" << std::endl;
                std::cout << "   Type = " << _pmt_timing_type << " (SPECTDC = 0; Ptb HLT = 1; CAEN-only = 3)." << std::endl;
                std::cout << "   Channel = " << _pmt_timing_ch << " (TDC ETRIG = 4; Ptb BNB Beam+Light = 2)." << std::endl;
            }
        }
        
        //Timing waveforms
        art::Handle<std::vector<raw::OpDetWaveform>> wfTimingHandle;
        std::vector<art::Ptr<raw::OpDetWaveform>> wf_timing_v;
        e.getByLabel(fTimingDecodeLabel, wfTimingHandle);
        
        if (!wfTimingHandle.isValid() || wfTimingHandle->size() == 0){
            throw cet::exception("TimingCAEN") << "No raw::OpDetWaveform found w/ tag " << fTimingDecodeLabel << ". Check data quality!";
        }
        else{

            art::fill_ptr_vector(wf_timing_v, wfTimingHandle);
            art::FindManyP<raw::pmt::BoardTimingInfo> timingBoardAssn(wf_timing_v, e, fTimingBoardLabel);

            //There is only 1 timing CAEN in the hardware and not every channel is saved

            if (fDebugTiming) std::cout << std::endl << "Found OpDetWaveform Timing size = " << wf_timing_v.size() << std::endl;

            for (size_t i = 0; i < wf_timing_v.size(); i++){
                art::Ptr<raw::OpDetWaveform> wf(wf_timing_v.at(i));
                //Get assn
                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(timingBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 )
                    throw cet::exception("TimingCAEN") << "No raw::wf::BoardTimingInfo found w/ tag" << fTimingBoardLabel <<". Check data quality!";

                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
                
                if (wf->ChannelNumber() != 902) continue;
                std::cout << wf->ChannelNumber() << std::endl;
                std::cout << wf_board->triggerTimeTag[0] << std::endl;
            }
        }
        wfTimingHandle.removeProduct();
        wf_timing_v.clear();

    }
    //-----------------------------------------------------------//
    
    if (fDebugTdc | fDebugPtb| fDebugTiming)
        std::cout <<"#--------------------------------------------------------#" << std::endl;

    //Fill once every event
    fTree->Fill();
}

void sbnd::TimingCAEN::beginJob()
{
    for(size_t i = 0; i < nCrt; i++){
        _hTopHatCRTT0[i] = tfs->make<TH1D>(Form("hCRTT0_Tagger_%s", crtTagger[i].c_str()), "", 100, -5, 5);
    }

    //Event Tree
    fTree = tfs->make<TTree>("events", "");
  
    fTree->Branch("run", &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event", &_event);
    fTree->Branch("tdc_ch0", &_tdc_ch0);
    fTree->Branch("tdc_ch1", &_tdc_ch1);
    fTree->Branch("tdc_ch2", &_tdc_ch2);
    fTree->Branch("tdc_ch3", &_tdc_ch3);
    fTree->Branch("tdc_ch4", &_tdc_ch4);
    fTree->Branch("tdc_ch0_utc", &_tdc_ch0_utc);
    fTree->Branch("tdc_ch1_utc", &_tdc_ch1_utc);
    fTree->Branch("tdc_ch2_utc", &_tdc_ch2_utc);
    fTree->Branch("tdc_ch3_utc", &_tdc_ch3_utc);
    fTree->Branch("tdc_ch4_utc", &_tdc_ch4_utc);

    if (fIncludePtb){
        fTree->Branch("ptb_hlt_trigger", &_ptb_hlt_trigger);
        fTree->Branch("ptb_hlt_timestamp", &_ptb_hlt_timestamp);
        fTree->Branch("ptb_hlt_trunmask", &_ptb_hlt_trunmask);
    }

    if (fIncludePmt){
        fTree->Branch("_pmt_timing_type", &_pmt_timing_type);
        fTree->Branch("_pmt_timing_ch", &_pmt_timing_ch);
    }
}

void sbnd::TimingCAEN::endJob()
{
}

void sbnd::TimingCAEN::ResetEventVars()
{
    _run = -1; _subrun = -1; _event = -1;
   
    _tdc_ch0.clear();
    _tdc_ch1.clear();
    _tdc_ch2.clear();
    _tdc_ch3.clear();
    _tdc_ch4.clear();

    _tdc_ch0_utc.clear();
    _tdc_ch1_utc.clear();
    _tdc_ch2_utc.clear();
    _tdc_ch3_utc.clear();
    _tdc_ch4_utc.clear();

    if (fIncludePtb){
        _ptb_hlt_trigger.clear();
        _ptb_hlt_timestamp.clear();
        _ptb_hlt_trunmask.clear();
    }

    if (fIncludePmt){
        _pmt_timing_type = -1;
        _pmt_timing_ch = -1;
    }
}

DEFINE_ART_MODULE(sbnd::TimingCAEN)
