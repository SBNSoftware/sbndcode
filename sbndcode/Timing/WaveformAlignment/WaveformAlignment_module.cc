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
//#include "sbndcode/ChannelMaps/PMT/PMTChannelMapService.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbndcode/Timing/SBNDRawTimingObj.h"


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
    double FindNearestTdc(const double tsFtrig, const std::vector<uint64_t> tdcVec);
    bool CheckShift(const double shift, const int boardId);

    void PlotFtrigCompare(const int flashId);
    void PlotFtrigFit(TGraph *g, TF1 *fitf, const bool converged, const int boardId, const int flashId);
    void ResetPlotVars();
private:

    //---GLOBAL PARAMETERS
    
    enum TimingStatus {
        kUndefined = -1,
        kGood,
        kFitFailure,
        kOutOfBound
    };
    
    // Plotting
    std::stringstream _plotName; //raw waveform hist name
    std::map<uint16_t, double> tickVec;
    std::map<uint16_t, std::vector<double>> xVec;
    std::map<uint16_t, std::vector<double>> yVec;
    std::map<uint16_t, double> boardMidX; 
    std::map<uint16_t, double> boardMidY; 

    // Shifting
    //std::map<uint16_t, std::vector<double>> boardJitter; 
    //std::map<uint16_t, std::vector<int>> boardStatus; 

    std::vector<double> boardJitter[9];
    std::vector<double> boardStatus[9];

    // Flash counter
    int nFtrigFlash;
    int nPmtFlash;

    //---TREE PARAMETERS
    TTree *fTree;
    art::ServiceHandle<art::TFileService> tfs;
    
    std::vector<std::string> pmtBoard{"0","1","2","3","4","5","6","7","8"};
    std::vector<int> pmtCol{1, 2, 3, 4, 5, 6, 7, 28, 30};

    int _run, _subrun, _event;

    // TDC stuff
    std::vector<uint64_t> _tdc_ch0;
    std::vector<uint64_t> _tdc_ch1;
    std::vector<uint64_t> _tdc_ch2;
    std::vector<uint64_t> _tdc_ch3;
    std::vector<uint64_t> _tdc_ch4;

    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;

    //---SERVICE
    //art::ServiceHandle<SBND::PMTChannelMapService> fPMTChannelMapService;
    
    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fTimingRefLabel;
    art::InputTag fPmtDecodeLabel;
    art::InputTag fPmtBoardLabel;
    art::InputTag fTimingDecodeLabel;
    art::InputTag fTimingBoardLabel;
    art::InputTag fFtrigDecodeLabel;
    art::InputTag fFtrigBoardLabel;

    // DAQ
    double fWfLength;
    std::vector<double> fTdc3CaenOffset;
    int fnPmtBoard;
    int fnTimingBoard;
    int fnChperBoard;

    // Shift configuration
    int fTickLbPmt;
    int fTickUbPmt;

    int fTickLbTiming;
    int fTickUbTiming;

    int fPmtFitBound;
    int fTimingFitBound;

    double fPmtJitterBound;
    double fTimingJitterBound;

    // Debug
    bool fDebugTdc;
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
    fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");
    fTimingRefLabel = p.get<art::InputTag>("fTimingRefLabel", "pmtdecoder");

    fPmtDecodeLabel = p.get<art::InputTag>("PmtDecodeLabel", "pmtdecoder:PMTChannels");
    fPmtBoardLabel = p.get<art::InputTag>("PmtBoardLabel", "pmtdecoder:PMTTiming");
    
    fTimingDecodeLabel = p.get<art::InputTag>("TimingDecodeLabel", "pmtdecoder:TimingChannels");
    fTimingBoardLabel = p.get<art::InputTag>("TimingBoardLabel", "pmtdecoder:TimingTiming");

    fFtrigDecodeLabel = p.get<art::InputTag>("FtrigDecodeLabel", "pmtdecoder:FTrigChannels");
    fFtrigBoardLabel = p.get<art::InputTag>("FtrigBoardLabel", "pmtdecoder:FTrigTiming");

    fWfLength = p.get<double>("WfLength", 5000);
    fTdc3CaenOffset = p.get<std::vector<double>>("Tdc3CaenOffset", {});
    
    fnPmtBoard = p.get<int>("nPmtBoard", 8);
    fnTimingBoard = p.get<int>("nTimingBoard", 1);

    fnChperBoard = p.get<int>("nChperBoard", 15);

    fTickLbPmt = p.get<int>("TickLbPmt", 1005);
    fTickUbPmt = p.get<int>("TickUbPmt", 1035);

    fTickLbTiming = p.get<int>("TickLbTiming", 975);
    fTickUbTiming = p.get<int>("TickUbTiming", 1005);

    fPmtFitBound = p.get<int>("PmtFitBound", 35);
    fTimingFitBound = p.get<int>("TimingFitBound", 35);

    fPmtJitterBound = p.get<int>("PmtJitterBound", 16);
    fTimingJitterBound = p.get<int>("TimingJitterBound", 68);

    fDebugTdc = p.get<bool>("DebugTdc", false);
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
   
    // Board Timing
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

    if (fDebugTdc | fDebugTimeRef | fDebugFtrig | fDebugPmt | fDebugTiming)
        std::cout <<"#----------RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event <<"----------#\n";
    
    //---------------------------TDC-----------------------------//
    art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
    e.getByLabel(fTdcDecodeLabel, tdcHandle);

    if (!tdcHandle.isValid() || tdcHandle->size() == 0){
        throw cet::exception("WaveformAlignment") << "No sbnd::timing::DAQTimestamp found w/ tag " << fTdcDecodeLabel << ". Check data quality!";
    }
    else{
        
        std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);
        
        for (size_t i=0; i<tdc_v.size(); i++){
            auto tdc = tdc_v[i];
            const std::string name  = tdc.Name();
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();

            if (fDebugTdc) std::cout << "Ch " << name
                                << ", ts (ns) = " << ts%uint64_t(1e9)
                                << ", sec (s) = " << ts/uint64_t(1e9)
                                << std::endl;
            if(ch == 0){
                _tdc_ch0.push_back(ts%uint64_t(1e9));
            }
            if(ch == 1){
                _tdc_ch1.push_back(ts%uint64_t(1e9));
            }
            if(ch == 2){
                _tdc_ch2.push_back(ts%uint64_t(1e9));
            }
            if(ch == 3){
                _tdc_ch3.push_back(ts%uint64_t(1e9));
            }
            if(ch == 4){
                _tdc_ch4.push_back(ts%uint64_t(1e9));
            }
        }
    }
    tdcHandle.removeProduct();

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
        
        int nTotalBoard = fnPmtBoard + fnTimingBoard;
    
        if (std::fmod(wf_ftrig_v.size(), nTotalBoard) != 0)
            throw cet::exception("WaveformAlignment") << "raw::OpDetWaveform found w/ tag " << fFtrigDecodeLabel << " does not have the same number of flashes per board. Check data quality!";

        nFtrigFlash = wf_ftrig_v.size() / nTotalBoard;

        if (fDebugFtrig)
            std::cout << std::endl << "Found OpDetWaveform FTRIG size = " << wf_ftrig_v.size() << ", nFTRIG per board = " << nFtrigFlash << std::endl;

        for (int flashIdx = 0; flashIdx < nFtrigFlash; flashIdx++){

            if (fDebugFtrig) std::cout << "   Looping over OpDet id " << flashIdx << "..." << std::endl;

            for (int boardIdx = 0; boardIdx < nTotalBoard; boardIdx++){ //wf->ChannelNumber() and j are the same

                int status = kUndefined;

                int wfIdx = flashIdx + boardIdx * nFtrigFlash;
                art::Ptr<raw::OpDetWaveform> wf(wf_ftrig_v.at(wfIdx));

                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(ftrigBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 )
                    throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fFtrigBoardLabel <<". Check data quality!";

                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
        
                //Fit FTRIG rising edge
                std::pair<double, double> midPoint = FitFtrig(wf, flashIdx, wf->ChannelNumber());
                double tickFtrig = midPoint.first; 
                double shift = 0; //Default as 0

                //Check if fit converged
                if (tickFtrig == std::numeric_limits<double>::min()){
                    status = kFitFailure;
                    if (fDebugFtrig)  std::cout << "    board id = " << wf->ChannelNumber() << " has FTRIG fit failure." <<  std::endl;
                }else{
                    
                    //convert tick to ts
                    double tsFtrig = wf_board->triggerTimeTag[0] - (fWfLength - tickFtrig)*2.0; 

                    //TDC as reference frame TODO: add picosecond correction
                    if ((_pmt_timing_type == 0) & (_pmt_timing_ch == 4)){
                        
                        double tsFrame = FindNearestTdc(tsFtrig, _tdc_ch3); 
                        
                        double tempShift = std::abs(tsFtrig - tsFrame) - fTdc3CaenOffset[wf->ChannelNumber()]; 

                        if (CheckShift(tempShift, wf->ChannelNumber())){
                            status = kGood;
                            //if (fDebugFtrig)  std::cout << "    board id = " << wf->ChannelNumber() << " has rising tick value = " << tickFtrig << " and shift value = " << shift << " ns." << std::endl;
                            shift = tempShift;
                        }else{
                            status = kOutOfBound; //either jitter too much or no reference frame 
                        }

                    }

                    //TODO: Add reference for PTB and CAEN only
                    
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
                    //boardMidX[wf->ChannelNumber()] = tsFtrig;
                    //boardMidY[wf->ChannelNumber()] = midPoint.second;
                    tickVec[wf->ChannelNumber()] = tickFtrig;
                    for (size_t i = 0; i < wf->Waveform().size(); i++){
                        yVec[wf->ChannelNumber()].push_back((double)wf->Waveform()[i]);
                        xVec[wf->ChannelNumber()].push_back(wf_board->triggerTimeTag[0] - (fWfLength - i)*2.0);
                    }
                }
            } //Done looping over boards

            if (fSaveCompare){
                PlotFtrigCompare(flashIdx);
                ResetPlotVars();
            }
        } //Done looping over OpDetWf

        for (size_t i = 0; i < std::size(boardJitter); i++){
            raw::pmt::BoardAlignment board_align;
            board_align.boardId = i;
            
            std::vector<double> shift_vec;
            for (auto j: boardJitter[i]){
                shift_vec.push_back(j);
            }
            board_align.shift = shift_vec;
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
        
        if (std::fmod((wf_pmt_v.size() / fnChperBoard), fnPmtBoard) != 0)
            throw cet::exception("WaveformAlignment") << "raw::OpDetWaveform found w/ tag " << fPmtDecodeLabel << " does not have the same number of flashes per board. Check data quality!";
        nPmtFlash = (wf_pmt_v.size() / fnChperBoard) / fnPmtBoard;
        
        if (fDebugPmt) std::cout << std::endl << "Found OpDetWaveform PMT size = " << wf_pmt_v.size() << ", nPmtFlash per board = " << nPmtFlash << std::endl;

        for (int boardIdx = 0 ; boardIdx < fnPmtBoard; boardIdx++){
            
            if (fDebugPmt) std::cout << "   Looping over board " << boardIdx << "..." << std::endl;
            
            for (int flashIdx = 0; flashIdx < nPmtFlash; flashIdx++){

                if (fDebugPmt) std::cout << "     flash " << flashIdx << " has shift value " << boardJitter[boardIdx][flashIdx] << std::endl;

                for (int chIdx = 0; chIdx < fnChperBoard; chIdx++){

                    int wfIdx = boardIdx * nPmtFlash * fnChperBoard + flashIdx * fnChperBoard + chIdx;
                    art::Ptr<raw::OpDetWaveform> wf(wf_pmt_v.at(wfIdx));

                    //Get assn
                    std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(pmtBoardAssn.at(wf.key()));
                    if(wf_board_v.size() != 1 )
                        throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fPmtBoardLabel <<". Check data quality!";

                    art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
                    art::Ptr<raw::pmt::BoardAlignment> wf_align = make_align_ptr(boardIdx);

                    //TODO: turn this to database
                    //SBND::PMTChannelMapService::PMTInfo_t pmtInfo = fPMTChannelMapService->GetPMTInfoFromChannelID(wf->ChannelNumber());
                    //if(!pmtInfo.valid)
                    //    throw cet::exception("WaveformAlignment") << "No PMTChannelMapService found for ch " << wf->ChannelNumber() << std::endl;
                    //double total_transit = pmtInfo.TotalTransit;

                    //TODO: Get cable length + PMT Response from Channel Map Service
                    //double correction = boardJitter[boardIdx][flashIdx] - total_transit;// - pmt response
                    double correction = boardJitter[boardIdx][flashIdx];// - total_transit;// - pmt response
                    //double correction = total_transit;// - pmt response
                 
                    double new_ts = wf->TimeStamp() + correction; //ns to us
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
        int nChTiming = wf_timing_v.size()/nPmtFlash;

        if (fDebugTiming) std::cout << std::endl << "Found OpDetWaveform Timing size = " << wf_timing_v.size() << ", nFlash per board = " << nPmtFlash << std::endl;

        //TODO: fix hardcoded
        int boardIdx = 8;

        if (fDebugTiming) std::cout << "    Looping over board " << boardIdx << "..." << std::endl;

        for (int flashIdx = 0; flashIdx < nPmtFlash; flashIdx++){

            if (fDebugTiming) std::cout << "     flash " << flashIdx << " has shift value " << boardJitter[boardIdx][flashIdx] << std::endl;

            for (int chIdx = 0; chIdx < nChTiming; chIdx++){

                int wfIdx = flashIdx * nChTiming + chIdx;
                art::Ptr<raw::OpDetWaveform> wf(wf_timing_v.at(wfIdx));

                //Get assn
                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(timingBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 )
                    throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found w/ tag" << fTimingBoardLabel <<". Check data quality!";

                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
                art::Ptr<raw::pmt::BoardAlignment> wf_align = make_align_ptr(boardIdx);

                //Timing CAEN is not connected to PMT
                double correction = boardJitter[boardIdx][flashIdx]; 
             
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

    if (fDebugTdc | fDebugFtrig) 
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
    fTree->Branch("tdc_ch0", &_tdc_ch0);
    fTree->Branch("tdc_ch1", &_tdc_ch1);
    fTree->Branch("tdc_ch2", &_tdc_ch2);
    fTree->Branch("tdc_ch3", &_tdc_ch3);
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
   
    _tdc_ch0.clear();
    _tdc_ch1.clear();
    _tdc_ch2.clear();
    _tdc_ch3.clear();
    _tdc_ch4.clear();

    _pmt_timing_type = -1;
    _pmt_timing_ch = -1;

    nFtrigFlash = -1;
    nPmtFlash = -1;

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
   
    int fitLb = std::numeric_limits<int>::min();
    int fitUb = std::numeric_limits<int>::min();
    
    if (boardId == 8){
        //TODO: Maybe need to change the equation for Timing CAEN?
        fitLb = risingTickGuess - 85;
        fitUb = risingTickGuess + 55;
    }else{
        fitLb = risingTickGuess - fPmtFitBound; //25;
        fitUb = risingTickGuess + fPmtFitBound; //25;
    }

    g->GetXaxis()->SetRangeUser(x[fitLb], x[fitUb]);

    std::string rising_edge = "1/(1+exp(([0]-x)/[1]))*[2]+[3]";

    TF1 *fitf = new TF1("fitf", rising_edge.c_str(), x[fitLb], x[fitUb]);
    fitf->SetParameter(0, x[risingTickGuess]);
    fitf->SetParameter(1, 0.3);
    fitf->SetParameter(2, amplitude);
    fitf->SetParameter(3, baseline);

    //---------Fit
    TFitResultPtr fp = g->Fit(fitf,"SRQ","", x[fitLb], x[fitUb]);
    bool converged = !(bool)(int(fp));

    //-------Check if fit converges
    std::pair<double, double> midPoint(std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
    if(converged) {
        double parA = fitf->GetParameter(0);
        double parMag = fitf->GetParameter(2);
        double parOffset = fitf->GetParameter(3);

        double midPointY = parMag/2 +  parOffset;
        double midPointX = fitf->GetX(midPointY, parA - 20, parA + 20);

        midPoint.first = midPointX;
        midPoint.second = midPointY;
    }
    
    //-------Save fits to plots
    if(converged && fSaveGoodFit) PlotFtrigFit(g, fitf, converged, boardId, flashId);
    if(!converged && fSaveBadFit) PlotFtrigFit(g, fitf, converged, boardId, flashId);

    return midPoint; 
}

double  sbnd::WaveformAlignment::FindNearestTdc(const double tsFtrig, const std::vector<uint64_t> tdcVec)
{
    double smallestDiff = std::numeric_limits<double>::max();
    double nearestTdc = std::numeric_limits<double>::min();

    for (size_t i = 0; i < tdcVec.size(); i++){
       double tdc = (double)tdcVec[i];
       double diff = std::abs(tsFtrig - tdc);
       
       if (diff < smallestDiff){
            smallestDiff = diff;
            nearestTdc = tdc;
       }
    }

    return nearestTdc;
}

bool sbnd::WaveformAlignment::CheckShift(const double shift, const int boardId)
{
    if (boardId == 8){
        if (std::abs(shift) > fTimingJitterBound){
            if (fDebugFtrig)  std::cout << "    board id = " << boardId << " has jittering " << shift << " > " << fTimingJitterBound <<  std::endl;
            return false;
        }
    }else{
        if (std::abs(shift) > fPmtJitterBound){
            if (fDebugFtrig)  std::cout << "    board id = " << boardId << " has jittering " << shift << " > " << fPmtJitterBound <<  std::endl;
            return false;
        }
    }

    //TODO: add more check?
    
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
    int x_lb = std::abs(tickVec[0] - 2);
    int x_ub = std::abs(tickVec[0] + 4);

    //---------Find baseline value
    double baseline = 0.;
    for(unsigned int i = 0 ; i < 800; i++) {
      baseline += yVec[0][i];
    }
    baseline /= 800;

    for(int k = 0; k < (fnPmtBoard + fnTimingBoard); k++){
    
        std::vector<double> x1 = xVec[k];
        std::vector<double> y = yVec[k];

        if ((x1.size() == 0) | (y.size() == 0)) continue;

        //-------Apply shift
        double shift = boardJitter[k][flashId];
        std::vector<double> x2;
        for (size_t i = 0; i < x1.size(); i++){
            x2.push_back(x1[i] + shift);
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
        g1->SetTitle(Form("Board_%s", pmtBoard[k].c_str()));
        g1->SetMarkerStyle(5);
        g1->SetMarkerSize(1);
        g1->SetMarkerColor(pmtCol[k]);
        g1->SetLineColor(pmtCol[k]);

        auto *g2 = new TGraph(fWfLength, &x2[0], &y[0]);
        g2->SetTitle(Form("Board_%s", pmtBoard[k].c_str()));
        g2->SetMarkerStyle(5);
        g2->SetMarkerSize(1);
        g2->SetMarkerColor(pmtCol[k]);
        g2->SetLineColor(pmtCol[k]);

        //auto *p = new TGraph();
        //p->SetPoint(0, boardMidX[k] + shift, boardMidY[k] + baseline_shift);
        //p->SetMarkerStyle(50);
        //p->SetMarkerSize(2);
        //p->SetMarkerColor(pmtCol[k]);
        //mg2->Add(p);

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

    //mg1->GetXaxis()->SetLimits(xVec[0][x_lb+68], xVec[0][x_ub+68]);
    mg1->GetXaxis()->SetLimits(xVec[0][x_lb], xVec[0][x_ub]);
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

    //mg2->GetXaxis()->SetLimits(xVec[0][x_lb+68], xVec[0][x_ub+68]);
    mg2->GetXaxis()->SetLimits(xVec[0][x_lb], xVec[0][x_ub]);
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
