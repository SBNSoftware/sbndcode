////////////////////////////////////////////////////////////////////////
// Class:       WaveformAlignment
// Plugin Type: Producer
// File:        WaveformAlignment_module.cc
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

#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TSystem.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

#include <bitset>

namespace sbnd {
    class WaveformAlignment;
}


class sbnd::WaveformAlignment : public art::EDProducer {
public:
    explicit WaveformAlignment(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WaveformAlignment(WaveformAlignment const&) = delete;
    WaveformAlignment(WaveformAlignment&&) = delete;
    WaveformAlignment& operator=(WaveformAlignment const&) = delete;
    WaveformAlignment& operator=(WaveformAlignment&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    void ResetEventVars();

    std::pair<double, double> FitFtrig(art::Ptr<raw::OpDetWaveform> wf, const int flashId, const int boardId);
    uint64_t FindNearestFrame(const uint64_t tsFtrig, const std::vector<uint64_t> frameVec, const double frameOffset);
    bool IsSecondRollover(const uint64_t tsFtrig, const uint64_t tsFrame);
    double ComputeShift(const uint64_t decodeFrame, const double decodeTs,  const double tickFtrig
                                            , const std::vector<uint64_t> frameVec
                                            , const double frameOffset
                                            , const bool isFrameEarly
                                            );
    double ComputeShiftPs(const uint64_t decodeFrame, const double decodeTs,  const double tickFtrig
                                            , const std::vector<uint64_t> frameVec
                                            , const std::vector<uint64_t> frameVecPs
                                            , const double frameOffset
                                            , const bool isFrameEarly
                                            );
    bool CheckShift(const double shift, const int boardId);

    void PlotFtrigCompare(const int flashId);
    void PlotFtrigFit(TGraph *g, TF1 *fitf, const bool converged, const int boardId, const int flashId);
    void ResetPlotVars();
private:

    //---GLOBAL PARAMETERS
    
    enum TimingStatus {
        kUndefined = -99, //Something went wrong
        kFitFailure = -1, //Fit does not converge - default 0
        kMissingTDC = -11, //Missing hardware information from TDC
        kOutOfBoundTDC = -12, //Cannot find correct frame from TDC
        kMissingPTB = -21, //Missing hardware information from PTB
        kOutOfBoundPTB = -22, //Cannot find correct frame from PTB
        kOutOfBoundCAEN = -3, //Cannot align with other CAEN digitisers since 
        kGoodTDC = 1, //Align using TDC O(2ns)
        kGoodPTB = 2, //Align using PTB O(2ns)
        kGoodCAEN = 3 //Align with CAEN earlier in the chain O(8ns)
    };
    
    // Plotting
    std::stringstream _plotName; //raw waveform hist name
    std::map<uint16_t, double> tickVec;
    std::map<uint16_t, std::vector<double>> xVec;
    std::map<uint16_t, std::vector<double>> yVec;
    std::map<uint16_t, double> boardMidX; 
    std::map<uint16_t, double> boardMidY; 

    // Shifting
    std::vector<double> boardJitter[9];
    std::vector<double> boardStatus[9];
    std::vector<double> boardFrame;

    // Board counter
    int nTotalBoard;
    int nTimingBoard;
    int nPmtBoard;
    std::vector<int> boardId_v;

    // Flash counter
    int nFtrigFlash;

    //---TREE PARAMETERS
    TTree *fTree;
    art::ServiceHandle<art::TFileService> tfs;

    int _run, _subrun, _event;

    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;

    // TDC stuff
    std::vector<uint64_t> _tdc_ch3;
    std::vector<uint64_t> _tdc_ch3_ps;
    std::vector<uint64_t> _tdc_ch4;

    //PTB stuff
    std::vector<uint64_t> _ptb_hlt_trigger;
    std::vector<uint64_t> _ptb_hlt_timestamp;
    std::vector<uint64_t> _ptb_hlt_unmask_timestamp;
    std::vector<int> _ptb_hlt_trunmask;

    //---SERVICE
    sbndDB::PMTCalibrationDatabase const* fPMTCalibrationDatabaseService;
    
    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTimingRefLabel;
    art::InputTag fTdcDecodeLabel;
    art::InputTag fPtbDecodeLabel;
    art::InputTag fPmtDecodeLabel;
    art::InputTag fPmtBoardLabel;
    art::InputTag fTimingDecodeLabel;
    art::InputTag fTimingBoardLabel;
    art::InputTag fFtrigDecodeLabel;
    art::InputTag fFtrigBoardLabel;

    //PTB
    std::vector<int> fPtbAllowedHlts; 

    // DAQ
    double fWfLength;
    int fnChperBoard;
    std::vector<int> fPmtBoard;
    std::vector<int> fPmtCol;

    // Shift configuration
    int fTickLbPmt;
    int fTickUbPmt;
    
    int fTickLbTiming;
    int fTickUbTiming;
    
    int fFitBound;
    int fFitTries;
    std::vector<double> fGradientInitialGuess;

    double fPmtJitterLowerBound;
    double fPmtJitterUpperBound;
    double fTimingJitterLowerBound;
    double fTimingJitterUpperBound;
    
    //Cable length
    std::vector<double> fTdc3CaenOffset;
    std::vector<double> fPtbCaenOffset;
    
    bool fCorrectCableOnly;
    
    // Debug
    bool fDebugTimeRef;
    bool fDebugFtrig;
    bool fDebugPmt;
    bool fDebugTiming;

    // New product labels
    std::string fFtrigNewLabel;
    std::string fFtrigBoardNewLabel;

    std::string fPmtNewLabel;
    std::string fPmtBoardNewLabel;
    std::string fPmtAlignNewLabel;

    std::string fTimingNewLabel;
    std::string fTimingBoardNewLabel;
    std::string fTimingAlignNewLabel;

    // Plotting
    bool fSaveGoodFit;
    bool fSaveBadFit;
    bool fSaveCompare;
    std::string fSavePath;
};


sbnd::WaveformAlignment::WaveformAlignment(fhicl::ParameterSet const& p)
    : EDProducer{p}  // 
    // More initializers here.
{
    fTimingRefLabel = p.get<art::InputTag>("fTimingRefLabel", "pmtdecoder");
    fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");
    fPtbDecodeLabel = p.get<art::InputTag>("PtbDecodeLabel", "ptbdecoder");

    fPmtDecodeLabel = p.get<art::InputTag>("PmtDecodeLabel", "pmtdecoder:PMTChannels");
    fPmtBoardLabel = p.get<art::InputTag>("PmtBoardLabel", "pmtdecoder:PMTTiming");
    
    fTimingDecodeLabel = p.get<art::InputTag>("TimingDecodeLabel", "pmtdecoder:TimingChannels");
    fTimingBoardLabel = p.get<art::InputTag>("TimingBoardLabel", "pmtdecoder:TimingTiming");

    fFtrigDecodeLabel = p.get<art::InputTag>("FtrigDecodeLabel", "pmtdecoder:FTrigChannels");
    fFtrigBoardLabel = p.get<art::InputTag>("FtrigBoardLabel", "pmtdecoder:FTrigTiming");
    
    fPtbAllowedHlts = p.get<std::vector<int>>("PtbAllowedHlts", {1, 2, 3, 4, 5, 14, 15});

    fTdc3CaenOffset = p.get<std::vector<double>>("Tdc3CaenOffset", {101, 101, 101, 101, 101, 101, 101, 101, 101});
    fPtbCaenOffset = p.get<std::vector<double>>("PtbCaenOffset", {166, 166, 166, 166, 166, 166, 166, 166, 166});
  
    fWfLength = p.get<double>("WfLength", 5000);
    fnChperBoard = p.get<int>("nChperBoard", 15);
    
    //Pmt CAEN 0-7, Timing CAEN 8
    fPmtBoard = p.get<std::vector<int>>("PmtBoard", {0,1,2,3,4,5,6,7,8});
    fPmtCol = p.get<std::vector<int>>("PmtCol", {1, 2, 3, 4, 5, 6, 7, 28, 30});

    fTickLbPmt = p.get<int>("TickLbPmt", 1005);
    fTickUbPmt = p.get<int>("TickUbPmt", 1035);

    fTickLbTiming = p.get<int>("TickLbTiming", 975);
    fTickUbTiming = p.get<int>("TickUbTiming", 1005);

    fFitBound = p.get<int>("FitBound", 12);
    fFitTries = p.get<int>("FitTries", 100);
    fGradientInitialGuess = p.get<std::vector<double>>("GradientInitialGuess", {0.3, 0.25, 0.25, 0.25, 0.3, 0.3, 0.3, 0.3, 0.08});

    fPmtJitterLowerBound = p.get<int>("PmtJitterBound", -5);
    fPmtJitterUpperBound = p.get<int>("PmtJitterBound", 3);
    fTimingJitterLowerBound = p.get<int>("TimingJitterBound", -56);
    fTimingJitterUpperBound = p.get<int>("TimingJitterBound", -48);

    fCorrectCableOnly = p.get<bool>("CorrectCableOnly", false);

    fDebugTimeRef = p.get<bool>("DebugTimeRef", false);
    fDebugFtrig = p.get<bool>("DebugFtrig", false);
    fDebugPmt = p.get<bool>("DebugPmt", false);
    fDebugTiming = p.get<bool>("DebugTiming", false);
    
    fSaveGoodFit = p.get<bool>("SaveGoodFit", false);
    fSaveBadFit = p.get<bool>("SaveBadFit", false);
    fSaveCompare = p.get<bool>("SaveCompare", false);
    fSavePath = p.get<std::string>("SavePath", "");
    
    fFtrigNewLabel = p.get<std::string>("FtrigNewLabel","FTrigChannels");
    fFtrigBoardNewLabel = p.get<std::string>("FtrigBoardNewLabel","FTrigTiming");
    
    fPmtNewLabel = p.get<std::string>("PmtNewLabel","PMTChannels");
    fPmtBoardNewLabel = p.get<std::string>("PmtBoardNewLabel","PMTTiming");
    fPmtAlignNewLabel = p.get<std::string>("PmtAlignNewLabel","PMTAlign");
    
    fTimingNewLabel = p.get<std::string>("TimingNewLabel","TimingChannels");
    fTimingBoardNewLabel = p.get<std::string>("TimingBoardNewLabel","TimingTiming");
    fTimingAlignNewLabel = p.get<std::string>("TimingAlignNewLabel","TimingAlign");

    produces< std::vector< raw::OpDetWaveform > >(fFtrigNewLabel); 
    produces< art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fFtrigBoardNewLabel);
    
    produces< std::vector< raw::OpDetWaveform > >(fPmtNewLabel); 
    produces< art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fPmtBoardNewLabel);
    produces< art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform > >(fPmtAlignNewLabel);
    
    produces< std::vector< raw::OpDetWaveform > >(fTimingNewLabel); 
    produces< art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fTimingBoardNewLabel);
    produces< art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform > >(fTimingAlignNewLabel);
    
    produces< std::vector< raw::pmt::BoardAlignment > >();
}

void sbnd::WaveformAlignment::produce(art::Event& e)
{
    // Output data products
    
    // FTRIG
    std::unique_ptr< std::vector< raw::OpDetWaveform > > newFtrigWf (new std::vector< raw::OpDetWaveform >);
    art::PtrMaker<raw::OpDetWaveform> make_ftrigwf_ptr{e, fFtrigNewLabel}; 
    
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > newFtrigBoardAssn (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);

    // PMT
    std::unique_ptr< std::vector< raw::OpDetWaveform > > newPmtWf (new std::vector< raw::OpDetWaveform >);
    art::PtrMaker<raw::OpDetWaveform> make_pmtwf_ptr{e, fPmtNewLabel}; 
    
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > newPmtBoardAssn (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);
    std::unique_ptr< art::Assns<  raw::pmt::BoardAlignment, raw::OpDetWaveform> > newPmtAlignAssn (new art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform >);
   
    // TIMING
    std::unique_ptr< std::vector< raw::OpDetWaveform > > newTimingWf (new std::vector< raw::OpDetWaveform >);
    art::PtrMaker<raw::OpDetWaveform> make_timingwf_ptr{e, fTimingNewLabel}; 
    
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > newTimingBoardAssn (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);
    std::unique_ptr< art::Assns<  raw::pmt::BoardAlignment, raw::OpDetWaveform> > newTimingAlignAssn (new art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform >);
   
    // Board Alignment
    std::unique_ptr< std::vector< raw::pmt::BoardAlignment >> newBoardAlign (new std::vector< raw::pmt::BoardAlignment >);
    art::PtrMaker<raw::pmt::BoardAlignment> make_align_ptr{e}; 
    
    ResetEventVars();
    
    if(fSaveGoodFit | fSaveBadFit | fSaveCompare){
        gSystem->Exec(Form("mkdir -p %s", fSavePath.c_str()));
        
        gStyle->SetFrameLineWidth(2);
        gStyle->SetTextFont(62);
        gStyle->SetTextSize(0.07);
        gStyle->SetLabelFont(62, "xyz");
        gStyle->SetLabelSize(0.04, "xyz");
        gStyle->SetTitleSize(0.04, "xyz");
        gStyle->SetTitleFont(62, "xyz");
        gStyle->SetTitleX(0.55); //title X location
        gStyle->SetTitleY(1.1); //title Y location
        gStyle->SetTitleW(0.75); //title width
        gStyle->SetTitleH(0.3); //title height
        gStyle->SetTitleFontSize(0.04);
    }

    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  =  e.id().event();

    if (fDebugTimeRef | fDebugFtrig | fDebugPmt | fDebugTiming)
        std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";

    //------------------------PMT Timing--------------------------//
    art::Handle<raw::TimingReferenceInfo> timingRefHandle;
    e.getByLabel(fTimingRefLabel, timingRefHandle);

    if (!timingRefHandle.isValid()){
        throw cet::exception("WaveformAlignment") << "No raw::TimingReferenceInfo found w/ tag " << fTimingRefLabel << ". Check data quality!";
    }
    else{
        raw::TimingReferenceInfo const& pmt_timing(*timingRefHandle);
        
        _pmt_timing_type = pmt_timing.timingType;
        _pmt_timing_ch = pmt_timing.timingChannel;

        if (fDebugTimeRef){
            std::cout << "Timing Reference For Decoding PMT" << std::endl;
            std::cout << "   Type = " << _pmt_timing_type << " (SPECTDC = 0; PTB HLT = 1; CAEN-only = 3)." << std::endl;
            std::cout << "   Channel = " << _pmt_timing_ch << " (TDC ETRIG = 4; PTB BNB Beam+Light = 2)." << std::endl;
        }
    }
    timingRefHandle.removeProduct();
    
    //---------------------------TDC-----------------------------//
    art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
    e.getByLabel(fTdcDecodeLabel, tdcHandle);

    if (!tdcHandle.isValid() || tdcHandle->size() == 0){
        if (_pmt_timing_type == 0){
            throw cet::exception("WaveformAlignment") << "No sbnd::timing::DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!";
        }
    }
    else{

        std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);

        for (size_t i=0; i<tdc_v.size(); i++){
            auto tdc = tdc_v[i];
            const std::string name  = tdc.Name();
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();

            if ((ch == 3) || (ch ==4)){
                if (ch == 3) {
                    _tdc_ch3.push_back(ts);
                    _tdc_ch3_ps.push_back(tdc.TimestampPs());
                }
                if (ch == 4) _tdc_ch4.push_back(ts);

                if (fDebugTimeRef) std::cout << "Ch " << name
                                    << " sec (s) = " << ts/uint64_t(1e9)
                                    << ", ts (ns) = " << ts%uint64_t(1e9)
                                    << " picosecond = " << tdc.TimestampPs()
                                    << std::endl;
            }
        }
    }
    tdcHandle.removeProduct();
    //---------------------------PTB-----------------------------//
    art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
    std::vector<art::Ptr<raw::ptb::sbndptb>> ptb_v;
    e.getByLabel(fPtbDecodeLabel, ptbHandle);
    
    if ((!ptbHandle.isValid() || ptbHandle->size() == 0)){
        if (_pmt_timing_type == 1)
            throw cet::exception("WaveformAlignment") << "No raw::ptb::sbndptb found w/ tag " << fPtbDecodeLabel << ". Check data quality!";
    }
    else{
        
        art::fill_ptr_vector(ptb_v, ptbHandle);
    
        // HLT Words
        unsigned nHLTs = 0;

        for(auto const& ptb : ptb_v)
            nHLTs += ptb->GetNHLTriggers();
  
        _ptb_hlt_timestamp.resize(nHLTs);
        _ptb_hlt_trunmask.resize(nHLTs);
        _ptb_hlt_unmask_timestamp.resize(nHLTs);
  
        unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
        unsigned h_i = 0; //For trigger with bitmask
        
        for(auto const& ptb : ptb_v){
            for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i){
	            _ptb_hlt_timestamp[h_i] = ptb->GetHLTrigger(i).timestamp * 20; //20 ns/tick 
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
        } //End of loop over ptb in PTBVec

        if (fDebugTimeRef){
            for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
                std::cout << "HLT " << _ptb_hlt_trunmask[i] 
                                    << " sec (s) = " << _ptb_hlt_unmask_timestamp[i]/uint64_t(1e9)
                                    << ", ts (ns) = " << _ptb_hlt_unmask_timestamp[i]%uint64_t(1e9)
                                    <<std::endl;
            }
        }
    }

    //---------------------------FTRIG-----------------------------//
    art::Handle<std::vector<raw::OpDetWaveform>> wfFtrigHandle;
    std::vector<art::Ptr<raw::OpDetWaveform>> wf_ftrig_v;
    e.getByLabel(fFtrigDecodeLabel, wfFtrigHandle);

    if (!wfFtrigHandle.isValid() || wfFtrigHandle->size() == 0){
        throw cet::exception("WaveformAlignment") << "No raw::OpDetWaveform found w/ tag " << fFtrigDecodeLabel << ". Check data quality!";
    }
    else{

        art::fill_ptr_vector(wf_ftrig_v, wfFtrigHandle);
        art::FindManyP<raw::pmt::BoardTimingInfo> ftrigBoardAssn(wf_ftrig_v, e, fFtrigBoardLabel);
       
        //Check the board Id available in data, assuming the data structure is awalsy PMT boardId first, then Timing boardId is last
        for (size_t i = 0; i < wf_ftrig_v.size(); i++){
            art::Ptr<raw::OpDetWaveform> wf(wf_ftrig_v.at(i));
            int cnt = count(boardId_v.begin(), boardId_v.end(), wf->ChannelNumber());
            if(cnt == 0){
                if ((int)wf->ChannelNumber() == fPmtBoard.back()){
                    nTimingBoard++;
                }else{
                    nPmtBoard++;
                }
                boardId_v.push_back(wf->ChannelNumber());
            }
        }
        nTotalBoard = boardId_v.size();
        nFtrigFlash = wf_ftrig_v.size() / nTotalBoard;

        if (fDebugFtrig) std::cout << std::endl << "Found OpDetWaveform FTRIG size = " << wf_ftrig_v.size() << ", nFTRIG per board = " << nFtrigFlash << std::endl;

        for (int flashIdx = 0; flashIdx < nFtrigFlash; flashIdx++){
            if (fDebugFtrig) std::cout << "   Looping over OpDet id " << flashIdx << "..." << std::endl;

            for (int boardIdx = 0; boardIdx < nTotalBoard; boardIdx++){ //boardIdx and wf ChannelNumber are the same

                int wfIdx = flashIdx + boardIdx * nFtrigFlash;
                art::Ptr<raw::OpDetWaveform> wf(wf_ftrig_v.at(wfIdx));

                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(ftrigBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 )
                    throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fFtrigBoardLabel <<". Check data quality!";
                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
        
                //Fit FTRIG rising edge
                std::pair<double, bool> midPoint = FitFtrig(wf, flashIdx, wf->ChannelNumber());
                double tickFtrig = midPoint.first; 

                int status = kUndefined;
                double shift = 0; //Default as 0

                //Check if fit converged
                if (midPoint.second == false){
                    status = kFitFailure;
                    if (fDebugFtrig)  std::cout << "    board id = " << wf->ChannelNumber() << " has fit failure, status = " << status << ", shift value = " << shift << " ns." << std::endl;
                }else{
                    //TDC as reference frame if timing type is 0 
                    if (_pmt_timing_type == 0){
                        if (_tdc_ch4.size() == 0 || _tdc_ch3.size() == 0 || _tdc_ch3_ps.size() == 0){
                            status = kMissingTDC;
                        }else{
                            double tempShift = 0;
                            if (_tdc_ch3_ps.front() == 0){ //old decode version does not have picosecond
                                tempShift = ComputeShift( _tdc_ch4.front()
                                                          , wf->TimeStamp()
                                                          , tickFtrig
                                                          , _tdc_ch3
                                                          , fTdc3CaenOffset[wf->ChannelNumber()]
                                                          , false
                                                          );
                            }else{
                                tempShift = ComputeShiftPs( _tdc_ch4.front()
                                                            , wf->TimeStamp()
                                                            , tickFtrig
                                                            , _tdc_ch3
                                                            , _tdc_ch3_ps
                                                            , fTdc3CaenOffset[wf->ChannelNumber()]
                                                            , false
                                                            );
                            }
                            if (CheckShift(tempShift, wf->ChannelNumber())){
                                status = kGoodTDC;
                                shift = std::round(tempShift * 1000.0) / 1000.0; //round to 3 decimal place of ns
                                if(boardFrame.size() == 0){
                                    if((int)wf->ChannelNumber() != fPmtBoard.back()) boardFrame.push_back(wf->TimeStamp()*1e3 + tickFtrig*2.0 - shift);
                                }
                            }else{
                                status = kOutOfBoundTDC;
                            }
                        }
                    }
                    //PTB as reference frame if timing type is 1 
                    //OR timing type is 0 but cannot find the TDC frame
                    //TDC status error will be replaced with PTB status error
                    if((_pmt_timing_type == 1) | (status == kMissingTDC) | (status == kOutOfBoundTDC)){
                        //Find the HLT of ETRIG
                        uint64_t tsFrame = std::numeric_limits<uint64_t>::min();
                        for (size_t i = 0; i < _ptb_hlt_unmask_timestamp.size(); i++){
                            for (size_t j = 0; j < fPtbAllowedHlts.size(); j++){
                                if (_ptb_hlt_trunmask[i] == fPtbAllowedHlts[j]){
                                    tsFrame = _ptb_hlt_unmask_timestamp[i];
                                    break;
                                }
                            }
                        }
                        if(tsFrame == std::numeric_limits<uint64_t>::min() || _ptb_hlt_unmask_timestamp.size() == 0){   
                            status = kMissingPTB;
                        }else{
                            double tempShift = ComputeShift( tsFrame
                                                            , wf->TimeStamp()
                                                            , tickFtrig
                                                            , _ptb_hlt_unmask_timestamp
                                                            , fPtbCaenOffset[wf->ChannelNumber()]
                                                            , true
                                                            );
                            if (CheckShift(tempShift, wf->ChannelNumber())){
                                status = kGoodPTB;
                                shift = std::round(tempShift * 1000.0) / 1000.0; //round to 3 decimal place of ns
                                if(boardFrame.size() == 0){
                                    if((int)wf->ChannelNumber() != fPmtBoard.back()) boardFrame.push_back(wf->TimeStamp()*1e3 + tickFtrig*2.0 - shift);
                                }
                            }else{
                                status = kOutOfBoundPTB;
                            }
                        }
                    }
                    //CAEN as reference frame if timing type is 2
                    //OR timing type is 1 but cannot find the PTB frame
                    //PTB status error will be replaced with CAEN status error
                    if ((_pmt_timing_type == 2) | (status == kMissingPTB) | (status == kOutOfBoundPTB)){
                        if (boardFrame.size() == 0){
                            //Cannot be timing CAEN
                            if ((int)wf->ChannelNumber() != fPmtBoard.back()){
                                boardFrame.push_back(wf->TimeStamp()*1e3 + tickFtrig*2.0);
                                status = kGoodCAEN;
                                shift = 0;
                            }
                        }
                        else if (boardFrame.size() > 0){ 
                            double tsFrame = boardFrame.front();
                            double tsFtrig = wf->TimeStamp()*1e3 + tickFtrig*2.0;
                            double tempShift = tsFtrig - tsFrame;

                            if (CheckShift(tempShift, wf->ChannelNumber())){
                                status = kGoodCAEN;
                                shift = std::round(tempShift * 1000.0) / 1000.0; //round to 3 decimal place of ns
                            }else{
                                status = kOutOfBoundCAEN;
                            }
                        }
                    }
                    if (fDebugFtrig)  std::cout << "    board id = " << wf->ChannelNumber() << " status = " << status << ", rising tick = " << tickFtrig << " and shift value = " << shift << " ns." << std::endl;
                }
                
                //Save to vector 
                boardJitter[wf->ChannelNumber()].push_back(shift);
                boardStatus[wf->ChannelNumber()].push_back(status);

                //Save a fraction of FTRIG waveform
                int saveLb = std::numeric_limits<int>::min();
                int saveUb = std::numeric_limits<int>::min();
                if (wf->ChannelNumber() == 8){
                    saveLb = fTickLbTiming - 25;
                    saveUb = fTickLbTiming + 100;
                }else{
                    saveLb = fTickLbPmt - 25;
                    saveUb = fTickLbPmt + 100;
                }
                double saveTs = wf->TimeStamp() + saveLb*2.0; 
                std::vector<uint16_t> adc_vec(saveUb - saveLb, 0);

                for (int i = 0; i < (saveUb - saveLb); i++){
                    adc_vec[i] = wf->Waveform()[i+saveLb];
                }

                raw::OpDetWaveform new_wf(saveTs, wf->ChannelNumber(), adc_vec);
                newFtrigWf->push_back(new_wf);
                art::Ptr<raw::OpDetWaveform> wfPtr = make_ftrigwf_ptr(newFtrigWf->size()-1);
                newFtrigBoardAssn->addSingle(wf_board, wfPtr);

                if (fSaveCompare){
                    tickVec[wf->ChannelNumber()] = tickFtrig;
                    for (size_t i = 0; i < wf->Waveform().size(); i++){
                        yVec[wf->ChannelNumber()].push_back((double)wf->Waveform()[i]);
                        xVec[wf->ChannelNumber()].push_back(wf_board->triggerTimeTag[0] - (fWfLength - i)*2.0);
                    }
                }

            } //Done looping over boards

            if (fSaveCompare){
                if (yVec.size() == 0) continue;
                PlotFtrigCompare(flashIdx);
                ResetPlotVars();
            }
            boardFrame.clear();//clear the frame for next flash
        } //Done looping over OpDetWf

        // Save board alignment product
        for (size_t i = 0; i < std::size(boardId_v); i++){
            std::vector<double> shift_vec;
            std::vector<int> status_vec;
            int boardId = boardId_v[i];

            for (size_t j = 0; j < boardJitter[boardId].size(); j++){
                shift_vec.push_back(boardJitter[boardId][j]);
                status_vec.push_back(boardStatus[boardId][j]);
            }
            raw::pmt::BoardAlignment board_align(shift_vec, status_vec);
            newBoardAlign->push_back(board_align);
        }
    }
    wfFtrigHandle.removeProduct();
    wf_ftrig_v.clear();

    //---------------------------ALIGN PMT WAVEFORMS-----------------------------//
    art::Handle<std::vector<raw::OpDetWaveform>> wfPmtHandle;
    std::vector<art::Ptr<raw::OpDetWaveform>> wf_pmt_v;
    e.getByLabel(fPmtDecodeLabel, wfPmtHandle);
    
    if (!wfPmtHandle.isValid() || wfPmtHandle->size() == 0){
        throw cet::exception("WaveformAlignment") << "No raw::OpDetWaveform found w/ tag " << fPmtDecodeLabel << ". Check data quality!";
    }
    else{

        art::fill_ptr_vector(wf_pmt_v, wfPmtHandle);
        art::FindManyP<raw::pmt::BoardTimingInfo> pmtBoardAssn(wf_pmt_v, e, fPmtBoardLabel);
        
        if (fDebugPmt) std::cout << std::endl << "Found OpDetWaveform PMT size = " << wf_pmt_v.size() << ", nFtrigFlash per board = " << nFtrigFlash << std::endl;

        for (int boardIdx = 0 ; boardIdx < nPmtBoard; boardIdx++){
            if (fDebugPmt) std::cout << "   Looping over board " << boardId_v[boardIdx] << "..." << std::endl;
            
            for (int flashIdx = 0; flashIdx < nFtrigFlash; flashIdx++){
                if (fDebugPmt) std::cout << "     flash " << flashIdx << " has shift value " << boardJitter[boardId_v[boardIdx]][flashIdx] << std::endl;
                
                for (int chIdx = 0; chIdx < fnChperBoard; chIdx++){

                    int wfIdx = boardIdx * nFtrigFlash * fnChperBoard + flashIdx * fnChperBoard + chIdx;
                    art::Ptr<raw::OpDetWaveform> wf(wf_pmt_v.at(wfIdx));

                    //Get assn
                    std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(pmtBoardAssn.at(wf.key()));
                    if(wf_board_v.size() != 1 )
                        throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fPmtBoardLabel <<". Check data quality!";

                    art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
                    art::Ptr<raw::pmt::BoardAlignment> wf_align = make_align_ptr(boardIdx);

                    fPMTCalibrationDatabaseService = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();
                    double total_transit = fPMTCalibrationDatabaseService->getTotalTransitTime(wf->ChannelNumber());

                    double correction = 0;
                    correction -= total_transit; //total transit = propagation time of photon from PMT to digitiser
                    if(fCorrectCableOnly == false) correction -= boardJitter[boardId_v[boardIdx]][flashIdx];
                    correction /= 1000; //ns to us conversion

                    double new_ts = wf->TimeStamp() + correction; //us

                    std::vector<uint16_t> adc_vec(wf->Waveform().size(), 0);
                    for (size_t i = 0; i < wf->Waveform().size(); i++){
                        adc_vec[i] = wf->Waveform()[i];
                    }

                    raw::OpDetWaveform new_wf(new_ts, wf->ChannelNumber(), adc_vec);

                    newPmtWf->push_back(new_wf);
                    art::Ptr<raw::OpDetWaveform> wfPtr = make_pmtwf_ptr(newPmtWf->size()-1);

                    newPmtBoardAssn->addSingle(wf_board, wfPtr);
                    newPmtAlignAssn->addSingle(wf_align, wfPtr);
                }
            }
        }
    }
    wfPmtHandle.removeProduct();
    wf_pmt_v.clear();

    //---------------------------ALIGN TIMING WAVEFORMS-----------------------------//
    art::Handle<std::vector<raw::OpDetWaveform>> wfTimingHandle;
    std::vector<art::Ptr<raw::OpDetWaveform>> wf_timing_v;
    e.getByLabel(fTimingDecodeLabel, wfTimingHandle);
    
    if (!wfTimingHandle.isValid() || wfTimingHandle->size() == 0){
        throw cet::exception("WaveformAlignment") << "No raw::OpDetWaveform found w/ tag " << fTimingDecodeLabel << ". Check data quality!";
    }
    else{

        art::fill_ptr_vector(wf_timing_v, wfTimingHandle);
        art::FindManyP<raw::pmt::BoardTimingInfo> timingBoardAssn(wf_timing_v, e, fTimingBoardLabel);

        //There is only 1 timing CAEN in the hardware and not every channel is saved
        int nChTiming = wf_timing_v.size()/nFtrigFlash;
        if (fDebugTiming) std::cout << std::endl << "Found OpDetWaveform Timing size = " << wf_timing_v.size() << ", nFlash per board = " << nFtrigFlash << std::endl;

        if (fDebugTiming) std::cout << "    Looping over board " << fPmtBoard.back() << "..." << std::endl;
        for (int flashIdx = 0; flashIdx < nFtrigFlash; flashIdx++){

            if (fDebugTiming) std::cout << "     flash " << flashIdx << " has shift value " << boardJitter[fPmtBoard.back()][flashIdx] << std::endl;

            for (int chIdx = 0; chIdx < nChTiming; chIdx++){

                int wfIdx = flashIdx * nChTiming + chIdx;
                art::Ptr<raw::OpDetWaveform> wf(wf_timing_v.at(wfIdx));

                //Get assn
                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(timingBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 )
                    throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fTimingBoardLabel <<". Check data quality!";

                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
                art::Ptr<raw::pmt::BoardAlignment> wf_align = make_align_ptr(boardId_v.size()-1);

                //Timing CAEN is not connected to PMT
                double correction= 0;
                correction -= boardJitter[fPmtBoard.back()][flashIdx]; 
                correction /= 1000; //ns to us
             
                double new_ts = wf->TimeStamp() + correction; //ns to us
                std::vector<uint16_t> adc_vec(wf->Waveform().size(), 0);
                
                for (size_t i = 0; i < wf->Waveform().size(); i++){
                    adc_vec[i] = wf->Waveform()[i];
                }

                raw::OpDetWaveform new_wf(new_ts, wf->ChannelNumber(), adc_vec);

                newTimingWf->push_back(new_wf);
                art::Ptr<raw::OpDetWaveform> wfPtr = make_timingwf_ptr(newTimingWf->size()-1);

                newTimingBoardAssn->addSingle(wf_board, wfPtr);
                newTimingAlignAssn->addSingle(wf_align, wfPtr);
            }
        }
    }
    wfTimingHandle.removeProduct();
    wf_timing_v.clear();

    if (fDebugTimeRef | fDebugFtrig | fDebugPmt | fDebugTiming)
      std::cout <<"#--------------------------------------------------------#" << std::endl;
    
    //Put product in event
    e.put(std::move(newBoardAlign));

    e.put(std::move(newFtrigWf), fFtrigNewLabel);
    e.put(std::move(newFtrigBoardAssn), fFtrigBoardNewLabel);
    
    e.put(std::move(newPmtWf), fPmtNewLabel);
    e.put(std::move(newPmtBoardAssn), fPmtBoardNewLabel);
    e.put(std::move(newPmtAlignAssn), fPmtAlignNewLabel);
    
    e.put(std::move(newTimingWf), fTimingNewLabel);
    e.put(std::move(newTimingBoardAssn), fTimingBoardNewLabel);
    e.put(std::move(newTimingAlignAssn), fTimingAlignNewLabel);

    //Fill once every event
    fTree->Fill();
}

void sbnd::WaveformAlignment::beginJob()
{
    //Event Tree
    fTree = tfs->make<TTree>("events", "");
  
    fTree->Branch("run", &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event", &_event);
    fTree->Branch("tdc_ch3", &_tdc_ch3);
    fTree->Branch("tdc_ch3_ps", &_tdc_ch3_ps);
    fTree->Branch("tdc_ch4", &_tdc_ch4);

    fTree->Branch("pmt_timing_type", &_pmt_timing_type);
    fTree->Branch("pmt_timing_ch", &_pmt_timing_ch);
    
    fTree->Branch("shift_board_0", &boardJitter[0]);
    fTree->Branch("shift_board_1", &boardJitter[1]);
    fTree->Branch("shift_board_2", &boardJitter[2]);
    fTree->Branch("shift_board_3", &boardJitter[3]);
    fTree->Branch("shift_board_4", &boardJitter[4]);
    fTree->Branch("shift_board_5", &boardJitter[5]);
    fTree->Branch("shift_board_6", &boardJitter[6]);
    fTree->Branch("shift_board_7", &boardJitter[7]);
    fTree->Branch("shift_board_8", &boardJitter[8]);
    
    fTree->Branch("status_board_0", &boardStatus[0]);
    fTree->Branch("status_board_1", &boardStatus[1]);
    fTree->Branch("status_board_2", &boardStatus[2]);
    fTree->Branch("status_board_3", &boardStatus[3]);
    fTree->Branch("status_board_4", &boardStatus[4]);
    fTree->Branch("status_board_5", &boardStatus[5]);
    fTree->Branch("status_board_6", &boardStatus[6]);
    fTree->Branch("status_board_7", &boardStatus[7]);
    fTree->Branch("status_board_8", &boardStatus[8]);
}

void sbnd::WaveformAlignment::endJob()
{
}

void sbnd::WaveformAlignment::ResetEventVars()
{
    _run = -1; _subrun = -1; _event = -1;

    _pmt_timing_type = std::numeric_limits<uint16_t>::max();
    _pmt_timing_ch = std::numeric_limits<uint16_t>::max();
   
    _tdc_ch3.clear();
    _tdc_ch3_ps.clear();
    _tdc_ch4.clear();
    
    _ptb_hlt_trigger.clear();
    _ptb_hlt_timestamp.clear();
    _ptb_hlt_unmask_timestamp.clear();
    _ptb_hlt_trunmask.clear();

    nPmtBoard = 0;
    nTimingBoard = 0;
    nTotalBoard = -1;
    boardId_v.clear();

    nFtrigFlash = -1;

    boardFrame.clear();

    for (size_t i = 0; i < std::size(boardJitter); ++i){
        boardJitter[i].clear();
        boardStatus[i].clear();
    }
}

std::pair<double, double> sbnd::WaveformAlignment::FitFtrig(art::Ptr<raw::OpDetWaveform> wf, const int flashId, const int boardId)
{
    //--------Make TGraph
    size_t wfLength = wf->Waveform().size();
    std::vector<double> x(wfLength), y(wfLength);
    
    for (size_t i = 0; i < wfLength; i++){
        y[i] = (double)wf->Waveform()[i];
        x[i] = i;
    }

    TGraph *g = new TGraph(wfLength, &x[0], &y[0]);

    //---------Find baseline value
    double baseline = 0.;
    for(unsigned int i = 0 ; i < 800; i++) {
      baseline += y[i];
    }
    baseline /= 800;

    //-------Guess rising tick
    double tickLb = std::numeric_limits<double>::min();
    double tickUb = std::numeric_limits<double>::min();

    if (boardId == 8){
        tickLb = fTickLbTiming;
        tickUb = fTickUbTiming;
    }else{
        tickLb = fTickLbPmt;
        tickUb = fTickUbPmt;
    }

    double largestRise = std::numeric_limits<double>::min();
    int risingTickGuess = std::numeric_limits<int>::min();
    for(int i = tickLb; i < tickUb; i++){
        
        double diff = y[i + 1] -  y[i];
        if(diff > largestRise){
            largestRise = diff;
            risingTickGuess = i;
        }
    }

    //-------Find amplitude
    double amplitude = 0;
    for(int i = risingTickGuess; i < risingTickGuess+50; i++){
        amplitude += y[i];
    }
    amplitude /= 50;
    amplitude -= baseline;

    //---------Define Rising Edge
   
    int fitLb = risingTickGuess - fFitBound;
    int fitUb = risingTickGuess + fFitBound;

    g->GetXaxis()->SetRangeUser(x[fitLb], x[fitUb]);

    std::string rising_edge = "1/(1+exp(([0]-x)/[1]))*[2]+[3]";
    
    TF1 *fitf = new TF1("fitf", rising_edge.c_str(), x[fitLb], x[fitUb]);
    fitf->SetParameter(0, x[risingTickGuess]);
    fitf->SetParameter(1, fGradientInitialGuess[boardId]);
    fitf->SetParameter(2, amplitude);
    fitf->SetParameter(3, baseline);

    //---------Fit
    TFitResultPtr fp = g->Fit(fitf,"SRQ","", x[fitLb], x[fitUb]);
    bool converged = !(bool)(int(fp));
    
    //--------Force fitting a few more times if does not converge
    int tries = 1;
    while(!converged){

        fitf->FixParameter(0, fitf->GetParameter(0));
        fitf->SetParameter(1, fGradientInitialGuess[boardId]);
        fitf->FixParameter(2, fitf->GetParameter(2));
        fitf->FixParameter(3, fitf->GetParameter(3));

        TFitResultPtr fp = g->Fit(fitf,"SRQ","", x[fitLb], x[fitUb]);
        converged = !(bool)(int(fp));
        tries++;

        if((converged) | (tries == fFitTries)) {
            break;
        }
    }

    //-------Check if fit converges
    std::pair<double, bool> midPoint(std::numeric_limits<double>::min(), converged);
    if(converged) {
        double parA = fitf->GetParameter(0);
        double parMag = fitf->GetParameter(2);
        double parOffset = fitf->GetParameter(3);

        double midPointY = parMag/2 +  parOffset;
        double midPointX = fitf->GetX(midPointY, parA - 20, parA + 20);

        midPoint.first = midPointX;
    }

    //-------Save fits to plots
    if(converged && fSaveGoodFit) PlotFtrigFit(g, fitf, converged, boardId, flashId);
    if(!converged && fSaveBadFit) PlotFtrigFit(g, fitf, converged, boardId, flashId);

    return midPoint; 
}

uint64_t  sbnd::WaveformAlignment::FindNearestFrame(const uint64_t tsFtrig, const std::vector<uint64_t> frameVec, const double frameOffset)
{
    double smallestDiff = std::numeric_limits<double>::max();
    uint64_t nearestFrame = std::numeric_limits<uint64_t>::min();

    for (size_t i = 0; i < frameVec.size(); i++){
        uint64_t frame = (uint64_t)frameVec[i];
        uint64_t diff = std::numeric_limits<uint64_t>::min();
        
        if (frame > tsFtrig){
            diff = frame - tsFtrig;
        }
        else if (frame < tsFtrig){
            diff = tsFtrig - frame;
        }
        
        if (diff > frameOffset) diff -= frameOffset;
        else if (frameOffset > diff) diff = frameOffset - diff;

        if (diff < smallestDiff){
             smallestDiff = diff;
             nearestFrame = frame;
        }
    }
    return nearestFrame;
}

bool sbnd::WaveformAlignment::IsSecondRollover(const uint64_t tsFtrig, const uint64_t tsFrame){
    uint64_t tsFtrigSec = tsFtrig / uint64_t(1e9);
    uint64_t tsFrameSec = tsFrame / uint64_t(1e9);

    if(tsFtrigSec == tsFrameSec){
        return false;
    }
    return true;
}

double sbnd::WaveformAlignment::ComputeShiftPs(const uint64_t decodeFrame, const double decodeTs,  const double tickFtrig
                                            , const std::vector<uint64_t> frameVec
                                            , const std::vector<uint64_t> frameVecPs
                                            , const double frameOffset
                                            , const bool isFrameEarly
                                            ){
    //acquire UTC timestamp of the start of the wavef
    uint64_t wfTimestamp = 0;
    wfTimestamp += decodeFrame;
    uint64_t decodeTs_ns = decodeTs * 1e3;  //convert us to ns
    wfTimestamp += decodeTs_ns;

    //acquired UTC timestamp of ftrig
    uint64_t ftrigShift = tickFtrig*2.0;
    uint64_t tsFtrig = wfTimestamp + ftrigShift; 
    
    //find nearest frame using UTC timestamp
    uint64_t tsFrame = FindNearestFrame(tsFtrig, frameVec, frameOffset);

    //find index of frame in frameVecPs
    size_t frameIdx = std::numeric_limits<size_t>::min();
    for(size_t i = 0; i < frameVec.size(); i++){
        if (frameVec[i] == tsFrame){
            frameIdx = i;
            break;
        }
    }
    double tsFrame_ns = frameVecPs[frameIdx]/1e3; //convert ps to ns

    double tsFtrig_ns = wfTimestamp%uint64_t(1e9) + tickFtrig*2.0; //double to handle decimal place

    //Check if second rollover happened
    bool isRollover = IsSecondRollover(tsFtrig, tsFrame);
    if(isRollover){
        if (tsFtrig_ns < tsFrame_ns) tsFtrig_ns += 1e9; //add 1 second to tsFtrig
        else tsFrame_ns += 1e9; //add 1 second to tsFrame
    }

    double shift = 0;
    if (!isFrameEarly) shift =+ frameOffset - std::abs(tsFtrig_ns - tsFrame_ns); 
    else shift =+ std::abs(tsFtrig_ns - tsFrame_ns) - frameOffset;

    return shift;
}

double sbnd::WaveformAlignment::ComputeShift(const uint64_t decodeFrame, const double decodeTs,  const double tickFtrig
                                            , const std::vector<uint64_t> frameVec
                                            , const double frameOffset
                                            , const bool isFrameEarly
                                            ){
    //acquire UTC timestamp of the start of the wavef
    uint64_t wfTimestamp = 0;
    wfTimestamp += decodeFrame;
    uint64_t decodeTs_ns = decodeTs * 1e3;  //convert us to ns
    wfTimestamp += decodeTs_ns;

    //acquired UTC timestamp of ftrig
    uint64_t ftrigShift = tickFtrig*2.0;
    uint64_t tsFtrig = wfTimestamp + ftrigShift; 
    
    //find nearest frame using UTC timestamp
    uint64_t tsFrame = FindNearestFrame(tsFtrig, frameVec, frameOffset);

    double tsFtrig_ns = wfTimestamp%uint64_t(1e9) + tickFtrig*2.0; //double to handle decimal place
    double tsFrame_ns = tsFrame%uint64_t(1e9);

    //Check if second rollover happened
    bool isRollover = IsSecondRollover(tsFtrig, tsFrame);
    if(isRollover){
        if (tsFtrig_ns < tsFrame_ns) tsFtrig_ns += 1e9; //add 1 second to tsFtrig
        else tsFrame_ns += 1e9; //add 1 second to tsFrame
    }

    double shift = 0;
    if (!isFrameEarly) shift =+ frameOffset - std::abs(tsFtrig_ns - tsFrame_ns); 
    else shift =+ std::abs(tsFtrig_ns - tsFrame_ns) - frameOffset;

    return shift;
}

bool sbnd::WaveformAlignment::CheckShift(const double shift, const int boardId)
{
    if (boardId == 8){
        if ((shift > fTimingJitterUpperBound) | (shift < fTimingJitterLowerBound)){
            if (fDebugFtrig)  std::cout << "    board id = " << boardId 
                                        << " has jittering " << shift 
                                        << " out of bound [" << fTimingJitterLowerBound << ", " << fTimingJitterUpperBound << "]"
                                        <<  std::endl;
            return false;
        }
    }else{
        if ((shift > fPmtJitterUpperBound) | (shift < fPmtJitterLowerBound)){
            if (fDebugFtrig)  std::cout << "    board id = " << boardId 
                                        << " has jittering " << shift 
                                        << " out of bound [" << fPmtJitterLowerBound << ", " << fPmtJitterUpperBound << "]"
                                        <<  std::endl;
            return false;
        }
    }
    
    return true;
}

void sbnd::WaveformAlignment::PlotFtrigFit(TGraph *g, TF1 *fitf, const bool converged, const int boardId, const int flashId)
{

    _plotName.str(std::string());
    _plotName << "run" << _run  << "_subrun" << _subrun << "_event" << _event << "_board" << boardId << "_flash" << flashId;
    
    int fitCol = 8;

    if (!converged){
        fitCol = 2;
        _plotName << "_BADFIT";
    }
    
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.55);
    gStyle->SetOptStat(0); 
    gStyle->SetOptFit(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas *c = new TCanvas(_plotName.str().c_str(), _plotName.str().c_str(), 600, 500);
    c->cd();
    
    g->SetTitle(_plotName.str().c_str());
    g->GetXaxis()->SetTitle("Tick");
    g->GetYaxis()->SetTitle("ADC");
    
    g->SetMarkerStyle(5);
    g->SetMarkerSize(2);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(1);
    
    g->Draw("ALP");
    
    fitf->SetLineColor(fitCol);
    fitf->Draw("SAME");
    
    if (converged){
        TGraph *p = new TGraph();
        double halfPointY = fitf->GetParameter(2)/2 +  fitf->GetParameter(3);
        double risingEdge = fitf->GetX(halfPointY, fitf->GetParameter(0)  - 20, fitf->GetParameter(0) + 20);
        
        p->SetPoint(0, risingEdge, halfPointY);
        p->SetMarkerStyle(5);
        p->SetMarkerSize(4);
        p->SetMarkerColor(kRed);
        p->Draw("Psame");
    }
    
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.1);
    c->SetBottomMargin(0.25);
    
    g->SetMinimum(1800); g->SetMaximum(6500);
    g->GetYaxis()->SetNdivisions(505, true);
    g->GetXaxis()->SetNdivisions(515, true);
    g->GetXaxis()->SetTitleOffset(2.1);
    g->GetXaxis()->CenterTitle();
    g->GetXaxis()->SetNoExponent(1);
    
    for (int i=1;  i <= g->GetXaxis()->GetNlabels() + 2; i++){
        g->GetXaxis()->ChangeLabel(i, 60, 0.05, 31, -1, -1);
    }
    
    c->SaveAs(Form("%s/%s.png", fSavePath.c_str(), _plotName.str().c_str()));
}

void sbnd::WaveformAlignment::PlotFtrigCompare(const int flashId)
{
    auto *mg1 = new TMultiGraph();
    auto *mg2 = new TMultiGraph();
    auto *lg1 = new TLegend(0.2,0.5,0.35,0.9); 
    auto *lg2 = new TLegend(0.2,0.5,0.35,0.9); 
    
    int x_lb = 0;
    int x_ub = 0;

    for (auto& pair : tickVec) {
        x_lb = pair.second - 2;
        x_ub = pair.second + 4;
        break;
    }

    //Find the first filled boardId 
    int refBoardId = -1;
    for (const auto& pair : yVec) {
        refBoardId = pair.first;
        break;
    }

    //---------Find baseline value
    double baseline = 0.;
    for(unsigned int i = 0 ; i < 800; i++) {
      baseline += yVec[refBoardId][i];
    }
    baseline /= 800;

    for (const auto& pair : yVec) {
        int boardId = pair.first;
    
        std::vector<double> x1 = xVec[boardId];
        std::vector<double> y = yVec[boardId];

        if ((x1.size() == 0) | (y.size() == 0)) continue;

        //-------Apply shift
        double shift = boardJitter[boardId][flashId];
        std::vector<double> x2;
        for (size_t i = 0; i < x1.size(); i++){
            x2.push_back(x1[i] - shift);
        }
      
        //------Shift baseline
        double current_baseline = 0;
        for(size_t i = 0 ; i < 800; i++) {
          current_baseline += y[i];
        }
        current_baseline /= 800;

        double baseline_shift = baseline - current_baseline;
        for (size_t i = 0; i < y.size(); i++){
            y[i] += baseline_shift;
        }
        
        auto *g1 = new TGraph(fWfLength, &x1[0], &y[0]);
        g1->SetTitle(Form("Board_%s", std::to_string(fPmtBoard[boardId]).c_str()));
        g1->SetMarkerStyle(5);
        g1->SetMarkerSize(1);
        g1->SetMarkerColor(fPmtCol[boardId]);
        g1->SetLineColor(fPmtCol[boardId]);

        auto *g2 = new TGraph(fWfLength, &x2[0], &y[0]);
        g2->SetTitle(Form("Board_%s", std::to_string(fPmtBoard[boardId]).c_str()));
        g2->SetMarkerStyle(5);
        g2->SetMarkerSize(1);
        g2->SetMarkerColor(fPmtCol[boardId]);
        g2->SetLineColor(fPmtCol[boardId]);

        mg1->Add(g1);
        mg2->Add(g2);

        lg1->AddEntry(g1, g1->GetTitle(), "lp");
        lg2->AddEntry(g2, g2->GetTitle(), "lp");
    }
    _plotName.str(std::string());
    _plotName << "run" << _run  << "_subrun" << _subrun << "_event" << _event << "_flash" << flashId;

    TCanvas *c = new TCanvas(_plotName.str().c_str(), _plotName.str().c_str(), 600, 1000);
    c->Divide(1, 2, 0, 0);

    c->cd(1);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.3);

    mg1->Draw("ALP");

    mg1->SetTitle(_plotName.str().c_str());
    mg1->GetXaxis()->SetTitle("Time [ns]");
    mg1->GetYaxis()->SetTitle("ADC");

    mg1->SetMinimum(1800); mg1->SetMaximum(7000);
    mg1->GetYaxis()->SetNdivisions(505, true);

    mg1->GetXaxis()->SetLimits(xVec[refBoardId][x_lb], xVec[refBoardId][x_ub]);
    mg1->GetXaxis()->SetNdivisions(310, true);
    mg1->GetXaxis()->SetLabelOffset(-0.01);
    mg1->GetXaxis()->SetTitleOffset(3.1);
    mg1->GetXaxis()->SetNoExponent(1);

    for (int i=1;  i <= mg1->GetXaxis()->GetNlabels() + 2; i++){
        mg1->GetXaxis()->ChangeLabel(i, 60, 0.04, 31, -1, -1);
    }

    lg1->SetTextSize(0.03);
    lg1->Draw();
    
    c->Modified();

    c->cd(2);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.3);
    _plotName.str(std::string());
    _plotName << "run" << _run  << "_subrun" << _subrun << "_event" << _event << "_flash" << flashId << "_aligned";

    mg2->Draw("ALP");

    mg2->SetTitle(_plotName.str().c_str());
    mg2->GetXaxis()->SetTitle("Time [ns]");
    mg2->GetYaxis()->SetTitle("ADC");

    mg2->SetMinimum(1800); mg2->SetMaximum(7000);
    mg2->GetYaxis()->SetNdivisions(505, true);

    mg2->GetXaxis()->SetLimits(xVec[refBoardId][x_lb], xVec[refBoardId][x_ub]);
    mg2->GetXaxis()->SetNdivisions(310, true);
    mg2->GetXaxis()->SetLabelOffset(-0.01);
    mg2->GetXaxis()->SetTitleOffset(3.1);
    mg2->GetXaxis()->SetNoExponent(1);

    for (int i=1;  i <= mg2->GetXaxis()->GetNlabels() + 2; i++){
        mg2->GetXaxis()->ChangeLabel(i, 60, 0.04, 31, -1, -1);
    }

    lg2->SetTextSize(0.03);
    lg2->Draw();
    
    _plotName.str(std::string());
    _plotName << "run" << _run  << "_subrun" << _subrun << "_event" << _event << "_flash" << flashId;

    c->Modified();

    c->SaveAs(Form("%s/%s.png", fSavePath.c_str(), _plotName.str().c_str()));
}

void sbnd::WaveformAlignment::ResetPlotVars()
{
    tickVec.clear();
    xVec.clear();
    yVec.clear();
    boardMidX.clear();
    boardMidY.clear();
}
DEFINE_ART_MODULE(sbnd::WaveformAlignment)
