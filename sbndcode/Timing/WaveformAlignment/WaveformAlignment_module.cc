////////////////////////////////////////////////////////////////////////
// Class:       WaveformAlignment
// Plugin Type: Producer
// File:        WaveformAlignment_module.cc
//
// Generated at Tue Oct 29 09:11:32 2024 by Vu Chi Lan Nguyen using cetskelgen
// from cetlib version 3.18.02.
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
#include "sbndcode/ChannelMaps/PMT/PMTChannelMapService.h"
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

    double GetFtrigRisingEdge(art::Ptr<raw::OpDetWaveform> wf, const int flashId, const int tickLb, const int tickUb);
    double GetFtrigOffsetFromTdc(const double tickFtrig, art::Ptr<raw::pmt::BoardTimingInfo> rawTime, const std::vector<uint64_t> tdcVec, const int boardId);
    bool CheckTick(const double tick);
    double CheckShift(const double shift);

    void PlotFtrigCompare(std::map<uint16_t, std::vector<double>> xVec,
                        std::map<uint16_t, std::vector<double>> xVec_shift,
                        std::map<uint16_t, std::vector<double>> yVec,
                        std::map<uint16_t, double> tickVec,
                        const int flashId);
    void PlotFtrigFit(TGraph *g, TF1 *fitf, const bool converged, const int boardId, const int flashId);
    void ResetPlotVars();
private:

    //---GLOBAL PARAMETERS
    
    // Plotting
    std::stringstream _plotName; //raw waveform hist name
    std::map<uint16_t, double> tickVec;
    std::map<uint16_t, std::vector<double>> xVec;
    std::map<uint16_t, std::vector<double>> xVec_shift;
    std::map<uint16_t, std::vector<double>> yVec;

    // Shifting
    std::map<uint16_t, std::vector<double>> boardJitter; 

    //---TREE PARAMETERS
    TTree *fTree;
    art::ServiceHandle<art::TFileService> tfs;
    
    std::vector<std::string> pmtBoard{"0","1","2","3","4","5","6","7","8"};
    std::vector<int> pmtCol{1, 2, 3, 4, 5, 6, 7, 28, 30};
    size_t nCAEN = pmtBoard.size(); 
    TH1D* _hFTRIG[9];

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

    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;

    //---SERVICE
    art::ServiceHandle<SBND::PMTChannelMapService> fPMTChannelMapService;
    
    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fTimingLabel;
    art::InputTag fFtrigDecodeLabel;
    art::InputTag fFtrigBoardLabel;
    art::InputTag fPmtDecodeLabel;
    art::InputTag fPmtBoardLabel;

    //DAQ parameters
    double fWfLength;
    std::vector<double> fTdc3CaenOffset;
    int fnCAEN;
    int fnChperBoard;

    int fTickLbPmt;
    int fTickUbPmt;

    int fTickLbTiming;
    int fTickUbTiming;

    // Debug
    bool fDebugTdc;
    bool fDebugTimeRef;
    bool fDebugFtrig;
    bool fDebugPmt;

    // Products
    bool fApplyShift;
    std::vector<int> fExcludeFtrigBoard;
    std::string fPmtNewLabel;
    std::string fPmtBoardNewLabel;
    std::string fPmtAlignNewLabel;

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
    fTimingLabel = p.get<art::InputTag>("fTimingLabel", "pmtdecoder");

    fFtrigDecodeLabel = p.get<art::InputTag>("FtrigDecodeLabel", "pmtdecoder:FTrigChannels");
    fFtrigBoardLabel = p.get<art::InputTag>("FtrigBoardLabel", "pmtdecoder:FTrigTiming");

    fPmtDecodeLabel = p.get<art::InputTag>("PmtDecodeLabel", "pmtdecoder:PMTChannels");
    fPmtBoardLabel = p.get<art::InputTag>("PmtBoardLabel", "pmtdecoder:PMTTiming");

    fWfLength = p.get<double>("WfLength", 5000);
    fTdc3CaenOffset = p.get<std::vector<double>>("Tdc3CaenOffset", {});
    fnCAEN = p.get<int>("nCAEN", 9);
    fnChperBoard = p.get<int>("nChperBoard", 15);

    fTickLbPmt = p.get<int>("TickLbPmt", 1005);
    fTickUbPmt = p.get<int>("TickUbPmt", 1035);

    fTickLbTiming = p.get<int>("TickLbTiming", 975);
    fTickUbTiming = p.get<int>("TickUbTiming", 1005);

    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugTimeRef = p.get<bool>("DebugTimeRef", false);
    fDebugFtrig = p.get<bool>("DebugFtrig", false);
    fDebugPmt = p.get<bool>("DebugPmt", false);

    fApplyShift = p.get<bool>("ApplyShift", true);
    fExcludeFtrigBoard = p.get<std::vector<int>>("ExcludeFtrigBoard", {});

    fSaveGoodFit = p.get<bool>("SaveGoodFit", false);
    fSaveBadFit = p.get<bool>("SaveBadFit", false);
    fSaveCompare = p.get<bool>("SaveCompare", false);
    fSavePath = p.get<std::string>("SavePath", "");

    fPmtNewLabel = p.get<std::string>("PmtNewLabel","PMTChannels");
    fPmtBoardNewLabel = p.get<std::string>("PmtBoardNewLabel","PMTTiming");
    fPmtAlignNewLabel = p.get<std::string>("PmtAlignNewLabel","PMTAlign");

    produces< std::vector< raw::OpDetWaveform > >(fPmtNewLabel); 
    produces< art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fPmtBoardNewLabel);
    produces< std::vector< raw::pmt::BoardAlignment > >();
    produces< art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform > >(fPmtAlignNewLabel);
}

void sbnd::WaveformAlignment::produce(art::Event& e)
{

    // Output data products
    std::unique_ptr< std::vector< raw::OpDetWaveform > > newPmtWf (new std::vector< raw::OpDetWaveform >);
    std::unique_ptr< std::vector< raw::pmt::BoardAlignment >> newPmtAlign (new std::vector< raw::pmt::BoardAlignment >);

    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > newPmtBoardAssn (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);
    std::unique_ptr< art::Assns<  raw::pmt::BoardAlignment, raw::OpDetWaveform> > newPmtAlignAssn (new art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform >);

    art::PtrMaker<raw::OpDetWaveform> make_pmtwf_ptr{e, fPmtNewLabel}; 
    art::PtrMaker<raw::pmt::BoardAlignment> make_align_ptr{e}; 
    

    //art::PtrMaker<raw::pmt::BoardAlignment> make_wfshift_ptr{e};
    //std::unique_ptr< art::Assns<  raw::pmt::BoardAlignment, raw::OpDetWaveform> > pmtShiftAssn (new art::Assns< raw::pmt::BoardAlignment, raw::OpDetWaveform >);

    

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

    if (fDebugTdc | fDebugTimeRef | fDebugFtrig)
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
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();
            //const uint64_t  offset = tdc.Offset();
            const std::string name  = tdc.Name();

            if (fDebugTdc) std::cout << "Ch " << name
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

        tdcHandle.removeProduct();
        tdc_v.clear();
    }

    //------------------------PMT Timing--------------------------//
    art::Handle<raw::TimingReferenceInfo> timingRefHandle;
    e.getByLabel(fTimingLabel, timingRefHandle);

    if (!timingRefHandle.isValid()){
        throw cet::exception("WaveformAlignment") << "No raw::TimingReferenceInfo found w/ tag " << fTimingLabel << ". Check data quality!";
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
        timingRefHandle.removeProduct();
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
        
        if (fDebugFtrig) std::cout << "Found OpDetWaveform FTRIG size = " << wf_ftrig_v.size() << ", nFTRIG per board = " << wf_ftrig_v.size()/fnCAEN << std::endl << std::endl;

        int nWfperBoard = wf_ftrig_v.size() / fnCAEN;

        for (int i = 0; i < nWfperBoard; i++){

            if (fDebugFtrig) std::cout << "   Looping over OpDet id " << i << "..." << std::endl;

            for (int j = 0; j < fnCAEN; j++){

                int wfIdx = i + j * nWfperBoard;
                art::Ptr<raw::OpDetWaveform> wf(wf_ftrig_v.at(wfIdx));

                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(ftrigBoardAssn.at(wf.key()));
                if(wf_board_v.size() != 1 ) throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found associated with raw::OpDetWaveform!";
                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());

                //TODO: Check which reference frame the OpDetWaveform was decoded, assuming TDC for now
                //If it's not TDC then has to change how to get offset
                if (_pmt_timing_ch != 4) std::cout << "Timing Ref is not TDC  it is " << _pmt_timing_ch << std::endl; 
        
                //Check if certain CAEN board is excluded
                bool excludeThisBoard = false;
                for(auto ch: fExcludeFtrigBoard)
                    if (int(wf->ChannelNumber())==ch) excludeThisBoard = true;

                if (excludeThisBoard){
                    //FTRIG channel ID is actually digitiser board ID
                    if(fDebugFtrig) std::cout << "      board id = " << wf->ChannelNumber() << " is exluded." << std::endl;
                    continue;
                }

                //Determine FTRIG rising edge
                double tickFtrig = -9999;
                //TODO: Add CAEN09
                if (int(wf->ChannelNumber())==8){
                    tickFtrig = GetFtrigRisingEdge(wf, i, fTickLbTiming, fTickUbTiming);
                }else{
                    tickFtrig = GetFtrigRisingEdge(wf, i, fTickLbPmt, fTickUbPmt);
                }

                //Compute shift amount in ns
                double shift = 0;
                if (CheckTick(tickFtrig)) shift = GetFtrigOffsetFromTdc(tickFtrig, wf_board, _tdc_ch3, wf->ChannelNumber());
                //shift = CheckShift(shift);
                
                if (fDebugFtrig)  std::cout << "        board id = " << wf->ChannelNumber() << " has rising tick value = " << tickFtrig << " and shift value = " << shift << " ns." << std::endl;
                
                //Save to vector
                //TODO: manual add shift values 
                //boardJitter[wf->ChannelNumber()].push_back(i+j*nWfperBoard);
                boardJitter[wf->ChannelNumber()].push_back(shift);


                //Save to tree
                _hFTRIG[wf->ChannelNumber()]->Fill(shift);

                if (fSaveCompare){
                    tickVec[wf->ChannelNumber()] = tickFtrig;
                    for (size_t i = 0; i < wf->Waveform().size(); i++){
                        yVec[wf->ChannelNumber()].push_back((double)wf->Waveform()[i]);
                        xVec[wf->ChannelNumber()].push_back(wf_board->triggerTimeTag[0] - (fWfLength - i)*2.0);
                        xVec_shift[wf->ChannelNumber()].push_back(wf_board->triggerTimeTag[0] - (fWfLength - i)*2.0 + shift);
                    }
                }

            } //Done looping over boards

            if (fSaveCompare){
                PlotFtrigCompare(xVec, xVec_shift, yVec, tickVec, i);
                ResetPlotVars();
            }
        } //Done looping over OpDetWf

        wfFtrigHandle.removeProduct();
        wf_ftrig_v.clear();

        //TODO: Do I want to save just digitiser jittering or everything?
        for (auto i : boardJitter){
            raw::pmt::BoardAlignment board_shift;
            board_shift.boardId = i.first;
            
            std::vector<double> shift_vec;
            for (auto j: i.second){
                shift_vec.push_back(j);
            }
            board_shift.shift = shift_vec;

            newPmtAlign->push_back(board_shift);
        }
    }

    //---------------------------ALIGN PMT WAVEFORMS-----------------------------//
    // 
    //// Grab pmt waveforms
    //art::Handle<std::vector<raw::OpDetWaveform>> wfPmtHandle;
    //std::vector<art::Ptr<raw::OpDetWaveform>> wf_pmt_v;
    //e.getByLabel(fPmtDecodeLabel, wfPmtHandle);
    
    //if (!wfPmtHandle.isValid() || wfPmtHandle->size() == 0){
    //    throw cet::exception("WaveformAlignment") << "No raw::OpDetWaveform found w/ tag " << fPmtDecodeLabel << ". Check data quality!";
    //}
    //else{

    //    art::fill_ptr_vector(wf_pmt_v, wfPmtHandle);
    //    art::FindManyP<raw::pmt::BoardTimingInfo> pmtBoardAssn(wf_pmt_v, e, fPmtBoardLabel);
    //    
    //    int fnPmtCAEN = fnCAEN - 1; //remove 1 Timing CAEN

    //    if (fDebugPmt) std::cout << "Found OpDetWaveform PMT size = " << wf_pmt_v.size() << ", nWf per board = " << (wf_pmt_v.size()/fnPmtCAEN)/fnChperBoard<< std::endl << std::endl;
    //   
    //    int nWfperBoard = (wf_pmt_v.size() / fnPmtCAEN) / fnChperBoard;
    //    

    //    for (int i = 0 ; i < fnPmtCAEN; i++){
    //        for (int j = 0; j < nWfperBoard; j++){

    //            if (fDebugPmt) std::cout << "   Board " << i << " flash " << j << " has shift value " << boardJitter[i][j] << std::endl;

    //            for (int k = 0; k < fnChperBoard; k++){
    //                int wfIdx = i * nWfperBoard * fnChperBoard + j * fnChperBoard + k;
    //                art::Ptr<raw::OpDetWaveform> wf(wf_pmt_v.at(wfIdx));

    //                //Get assn
    //                std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> wf_board_v(pmtBoardAssn.at(wf.key()));
    //                if(wf_board_v.size() != 1 ) throw cet::exception("WaveformAlignment") << "No raw::wf::BoardTimingInfo found associated with raw::OpDetWaveform!";

    //                art::Ptr<raw::pmt::BoardTimingInfo> wf_board(wf_board_v.front());
    //                art::Ptr<raw::pmt::BoardAlignment> wf_align = make_align_ptr(i);

    //                //TODO: Get cable length + PMT Response from Channel Map Service
    //                double correction = boardJitter[i][j]; // - cable length - pmt response
    //             
    //                double new_ts = wf->TimeStamp() + correction; //ns to us
    //                std::vector<uint16_t> adc_vec(wf->Waveform().size(), 0);
    //                
    //                for (size_t i = 0; i < wf->Waveform().size(); i++){
    //                    adc_vec[i] = wf->Waveform()[i];
    //                }

    //                raw::OpDetWaveform new_wf(new_ts, wf->ChannelNumber(), adc_vec);

    //                newPmtWf->push_back(new_wf);
    //                art::Ptr<raw::OpDetWaveform> wfPtr = make_pmtwf_ptr(newPmtWf->size()-1);

    //                newPmtBoardAssn->addSingle(wf_board, wfPtr);
    //                newPmtAlignAssn->addSingle(wf_align, wfPtr);
    //            }
    //        }
    //    }

    //}
    //wfPmtHandle.removeProduct();

    if (fDebugTdc | fDebugFtrig) std::cout <<"#--------------------------------------------------------#" << std::endl;
    
    //Put product in event
    e.put(std::move(newPmtWf), fPmtNewLabel);
    e.put(std::move(newPmtBoardAssn), fPmtBoardNewLabel);
    e.put(std::move(newPmtAlign));
    e.put(std::move(newPmtAlignAssn), fPmtAlignNewLabel);

    //Fill once every event
    fTree->Fill();
}

bool sbnd::WaveformAlignment::CheckTick(const double tick)
{
    //TODO: I feel fishy about this, what does it mean when a fit fails?
    if (tick == 9999){
        if (fDebugFtrig) std::cout << "     Fail fitting FTRIG waveform. Check waveform fitting. Also check data quality" << std::endl;;
        return false;
    }
   
    //TODO: to include both PMT and Timing CAEN
    //if ((tick < fTickLbPmt) | (tick > fTickUbPmt)){
    //    std::stringstream warning;;
    //    warning << "     Rising tick position " << tick << " is out of tolerance "<< fTickLbPmt << " - " << fTickUbPmt << ". Check data quality." << std::endl;

    //    if (fDebugFtrig) std::cout << warning.str() << std::endl;
    //    throw cet::exception("WaveformAlignment") << warning.str();
    //}

    return true;
}

double sbnd::WaveformAlignment::CheckShift(const double shift)
{

    //TODO:
    //if (shift == 0){
    //}

    //TODO: Add rounding up 

    return shift;
}

void sbnd::WaveformAlignment::beginJob()
{
    for(size_t i = 0; i < nCAEN; i++){
        _hFTRIG[i] = tfs->make<TH1D>(Form("hFTRIG_RisingEdge_Board_%s", pmtBoard[i].c_str()), "", 32, -16, 16);
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

    fTree->Branch("_pmt_timing_type", &_pmt_timing_type);
    fTree->Branch("_pmt_timing_ch", &_pmt_timing_ch);
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

    _tdc_ch0_utc.clear();
    _tdc_ch1_utc.clear();
    _tdc_ch2_utc.clear();
    _tdc_ch3_utc.clear();
    _tdc_ch4_utc.clear();

    _pmt_timing_type = -1;
    _pmt_timing_ch = -1;

    boardJitter.clear();
}

double sbnd::WaveformAlignment::GetFtrigRisingEdge(art::Ptr<raw::OpDetWaveform> wf, const int flashId, const int tickLb, const int tickUb)
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
    double largestRise = -9999;
    int risingTickGuess = 0;
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
    int fitLb = risingTickGuess - 25;
    int fitUb = risingTickGuess + 25;
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
    double risingEdge = 9999;
    double halfPointY = 9999;
    if(converged) {
        double parA = fitf->GetParameter(0);
        double parMag = fitf->GetParameter(2);
        double parOffset = fitf->GetParameter(3);

        halfPointY = parMag/2 +  parOffset;
        risingEdge = fitf->GetX(halfPointY, parA - 20, parA + 20);
    }
    
    //-------Save fits to plots
    if(converged && fSaveGoodFit) PlotFtrigFit(g, fitf, converged, wf->ChannelNumber(), flashId);
    if(!converged && fSaveBadFit) PlotFtrigFit(g, fitf, converged, wf->ChannelNumber(), flashId);

    return risingEdge; 
}

double  sbnd::WaveformAlignment::GetFtrigOffsetFromTdc(const double tickFtrig, art::Ptr<raw::pmt::BoardTimingInfo> rawTime, const std::vector<uint64_t> tdcVec, const int boardId)
{
    //Get the first TTT of the waveform
    double tsFtrig = rawTime->triggerTimeTag[0] - (fWfLength - tickFtrig)*2.0; 

    //Find nearest TDC FTRIG
    double smallestDiff = 999999999;

    for (size_t i = 0; i < tdcVec.size(); i++){
       double tdc = (double)tdcVec[i];
       double diff = std::abs(tsFtrig - tdc);
       
       if (diff < smallestDiff){
            smallestDiff = diff;
       }
    }

    //Compute shift for alignment
    double shift = smallestDiff - fTdc3CaenOffset[boardId]; 

    return shift;
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
    //gStyle->SetOptFit(0);
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

void sbnd::WaveformAlignment::PlotFtrigCompare(std::map<uint16_t, std::vector<double>> xVec,
                                             std::map<uint16_t, std::vector<double>> xVec_shift,
                                             std::map<uint16_t, std::vector<double>> yVec,
                                             std::map<uint16_t, double> tickVec,
                                             const int flashId)
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

    for(int k = 0; k < fnCAEN; k++){
    
        std::vector<double> x1 = xVec[k];
        std::vector<double> x2 = xVec_shift[k];
        std::vector<double> y = yVec[k];

        if ((x1.size() == 0) | (y.size() == 0) | (x2.size() == 0)) continue;
       
        double current_baseline = 0;
        for(size_t i = 0 ; i < 800; i++) {
          current_baseline += y[i];
        }
        current_baseline /= 800;

        double baseline_shift = baseline - current_baseline;
        std::cout << baseline_shift << " " << baseline << " " << current_baseline << std::endl;
        for (size_t i = 0; i < y.size(); i++){
            y[i] += baseline_shift;
        }
        
        auto *g1 = new TGraph(fWfLength, &x1[0], &y[0]);
        g1->SetTitle(Form("Board_%s", pmtBoard[k].c_str()));
        g1->SetMarkerStyle(5);
        g1->SetMarkerSize(2);
        g1->SetMarkerColor(pmtCol[k]);
        g1->SetLineColor(pmtCol[k]);

        auto *g2 = new TGraph(fWfLength, &x2[0], &y[0]);
        g2->SetTitle(Form("Board_%s", pmtBoard[k].c_str()));
        g2->SetMarkerStyle(5);
        g2->SetMarkerSize(2);
        g2->SetMarkerColor(pmtCol[k]);
        g2->SetLineColor(pmtCol[k]);

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
    xVec_shift.clear();
    yVec.clear();
}
DEFINE_ART_MODULE(sbnd::WaveformAlignment)
