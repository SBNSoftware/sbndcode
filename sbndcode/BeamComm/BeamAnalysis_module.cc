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

    double GetFtrigRisingEdge(const art::Ptr<raw::OpDetWaveform> wf, uint16_t pp);

private:

    // Event Tree
    TTree *_tree;
    std::stringstream _histName; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;

    std::vector<std::string> crtTagger{"Bottom","South","North","West","East","TopLow","TopHigh"};
    size_t nCrt = crtTagger.size(); 
    TH1D* _hTopHatCRTT0[8];
    
    std::vector<std::string> pmtBoard{"0","1","2","3","4","5","6","7","8"};
    size_t nPmt = pmtBoard.size(); 
    TH1D* _hFTRIG[9];

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
    std::map<uint16_t, double> board_jitter;

    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fCrtSpacePointLabel;
    art::InputTag fPmtFtrigDecodeLabel;
    art::InputTag fPmtFtrigBoardLabel;
    art::InputTag fPmtTimingLabel;
    art::InputTag fOpHitLabel;

    // Debug
    bool fDebugTdc;
    bool fDebugCrt;
    bool fDebugPmt;
    
    // Which data products to include
    bool fIncludeCrt;
    bool fIncludePmt;

    std::vector<int> fExcludePmtCh;
    bool fSavePmt;
    std::string fSavePmtPath;
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

    fOpHitLabel = p.get<art::InputTag>("OpHitLabel","");

    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugCrt = p.get<bool>("DebugCrt", false);
    fDebugPmt = p.get<bool>("DebugPmt", true);

    fIncludeCrt = p.get<bool>("IncludeCrt", true);
    fIncludePmt = p.get<bool>("IncludePmt", true);

    fExcludePmtCh = p.get<std::vector<int>>("ExcludePmtCh", {1});
    fSavePmt = p.get<bool>("SavePmt", true);
    fSavePmtPath = p.get<std::string>("SavePmtPath", "/exp/sbnd/data/users/lnguyen/BeamComm/BeamAna/plots/");
}

void sbnd::BeamAnalysis::analyze(art::Event const& e)
{
    
    ResetEventVars();
  
    if(fSavePmt){
        gSystem->Exec(Form("mkdir -p %s", fSavePmtPath.c_str()));

        gStyle->SetFrameLineWidth(2);
        gStyle->SetTextFont(62);
        gStyle->SetTextSize(0.07);
        gStyle->SetLabelFont(62, "xyz");
        gStyle->SetLabelSize(0.05, "xyz");
        gStyle->SetTitleSize(0.05, "xyz");
        gStyle->SetTitleFont(62, "xyz");
    }

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
                
                double ts0_ns = crt_ts0.back() - (tdc_ch2.back() - tdc_ch4.back());
                double ts0_us = ts0_ns / 1'000;

                _hTopHatCRTT0[crt_tagger.back()]->Fill(ts0_us);
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

            if (fDebugPmt) std::cout << "Found OpDetWaveform FTRIG size = " << pmt_ftrig_v.size() << std::endl;

            for (auto const& pmt: pmt_ftrig_v){

                const std::vector<art::Ptr<raw::pmt::boardTimingInfo>> pmt_board_v(pmtBoardAssoc.at(pmt.key()));

                if(pmt_board_v.size() != 1 ) continue;

                const art::Ptr<raw::pmt::boardTimingInfo> pmt_board(pmt_board_v.front());
              
                bool excludeThisCh = false;
                for(auto const ch: fExcludePmtCh){
                    if (int(pmt->ChannelNumber())==ch) excludeThisCh = true;
                }

                if (excludeThisCh){
                    //FTRIG channel ID is actually digitiser board ID
                    if(fDebugPmt) std::cout << "   board id = " << pmt->ChannelNumber() << " is exluded." << std::endl << std::endl;; 
                    continue;
                }

                double tsFtrig = GetFtrigRisingEdge(pmt, pmt_board->postPercent);
                
                if (fDebugPmt){ 
                    std::cout << "   board id = " << pmt->ChannelNumber() << " has FTRIG rising edge = " << tsFtrig << " ns." << std::endl << std::endl;
                }

                board_jitter[pmt->ChannelNumber()] = tsFtrig;
                _hFTRIG[pmt->ChannelNumber()]->Fill(tsFtrig);
            }
        }

        //---------------------------OpHit-----------------------------//
        art::Handle<std::vector<recob::OpHit>> opHitHandle;
        std::vector<art::Ptr<recob::OpHit>> ophit_v;
        e.getByLabel(fOpHitLabel, opHitHandle);
        
        if (!opHitHandle.isValid() || opHitHandle->size() == 0){
            if (fDebugPmt) std::cout << "No OpHit products found." << std::endl;
        }
        else{
            art::fill_ptr_vector(ophit_v, opHitHandle);

            for (auto const& ophit: ophit_v){
                std::cout << "hello" << std::endl;

                //get channel ID
                uint64_t ch = ophit->OpChannel();
                std::cout << ch << std::endl;
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
    for(size_t i = 0; i < nCrt; i++){
        _hTopHatCRTT0[i] = tfs->make<TH1D>(Form("hCRTT0_Tagger_%s", crtTagger[i].c_str()), "", 100, -5, 5);
    }

    for(size_t i = 0; i < nPmt; i++){
        _hFTRIG[i] = tfs->make<TH1D>(Form("hFTRIG_RisingEdge_Board_%s", pmtBoard[i].c_str()), "", 100, -5, 5);
    }

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

double sbnd::BeamAnalysis::GetFtrigRisingEdge(const art::Ptr<raw::OpDetWaveform> wf, uint16_t pp){

    //---------Fill waveform hist
    _histName.str(std::string());
    _histName << "run" << _run  
            << "_subrun" << _subrun
            << "_event" << _event
            << "_pmtnum" << wf->ChannelNumber();
   
    unsigned int recordLength = wf->Waveform().size();

    TH1D *hWvfm = new TH1D(Form("run%s_subrun%s_event%s_pmtnum%s", std::to_string(_run).c_str(), std::to_string(_subrun).c_str(),std::to_string(_event).c_str(), std::to_string(wf->ChannelNumber()).c_str()), "", recordLength, 0, recordLength);

    for(unsigned int i = 0 ; i < recordLength; i++) {
      hWvfm->SetBinContent(i, (double)wf->Waveform()[i]);
    }

    //---------Find baseline value
    double baseline = 0.;
    unsigned int baselineNBin = recordLength * 0.4;
    for(unsigned int i = 0 ; i < baselineNBin; i++) {
      baseline += hWvfm->GetBinContent(i);
    }
    baseline /= baselineNBin;

    //-------Find rising/falling edge
    double postPercent = pp / 100.0;
    double largestRise = -9999;
    double largestFall = 9999;
    double risingTickGuess = 0;
    double fallingTickGuess = 0;
    for(unsigned int i = recordLength * postPercent; i < recordLength - 1; i++){
       double diff = hWvfm->GetBinContent(i + 1) -  hWvfm->GetBinContent(i);
       if(diff > largestRise){
            largestRise = diff;
            risingTickGuess = i;
       }
       if(diff < largestFall){
           largestFall = diff;
           fallingTickGuess = i;
       }
    }

    //-------Find amplitude
    double amplitude = 0;
    for(unsigned int i = risingTickGuess; i < fallingTickGuess; i++){
        amplitude += hWvfm->GetBinContent(i);
    }
    amplitude /= (fallingTickGuess - risingTickGuess);
    amplitude -= baseline;

    //---------Define Square Wave
    std::string square_wave = "(-1/(1+exp(([0]-x)/[4])) + 1/(1+exp(([1]-x)/[5])))*[2]+[3]";

    double wf_lb = recordLength * postPercent - 50;
    double wf_ub = recordLength * postPercent + 150;

    TF1 *fitf = new TF1("fitf", square_wave.c_str(), wf_lb, wf_ub);
    fitf->SetParameter(0, risingTickGuess);
    fitf->SetParameter(1, fallingTickGuess);
    fitf->SetParameter(2, amplitude);
    fitf->SetParameter(3, baseline);
    fitf->SetParameter(4, -1);
    fitf->SetParameter(5, -1);
    fitf->Draw();

    //---------Fit
    TFitResultPtr fp = hWvfm->Fit(fitf,"SRQ","", wf_lb, wf_ub);
    bool converged = !(bool)(int(fp));

    //-------Save fits to plots
    if(fSavePmt){
        
        gStyle->SetOptStat(0); 
        gStyle->SetOptFit(1);
    
        TCanvas *c = new TCanvas(_histName.str().c_str(), _histName.str().c_str());
        c->cd();

        hWvfm->GetXaxis()->SetTitle("Tick");
        hWvfm->GetYaxis()->SetTitle("ADC");
        hWvfm->GetXaxis()->SetRangeUser(wf_lb, wf_ub);
        
        hWvfm->SetLineColor(kBlack);
        hWvfm->SetLineWidth(2);
        hWvfm->Draw("histe");
        fitf->SetLineColor(kSpring-6);
        fitf->Draw("same");
        
        c->SaveAs(Form("%s/%s.png", fSavePmtPath.c_str(), _histName.str().c_str()));
        //c->SaveAs(Form("%s/%s.pdf", fSavePmtPath.c_str(), _histName.str().c_str()));
    }

    //-------Check if fit converges
    if(!converged) {
        if(fDebugPmt) std::cout << "Fitting FTRIG waveform does not converge!" << std::endl;
        return -99999;
    }

    double parA = fitf->GetParameter(0);
    double parOffset = fitf->GetParameter(2);
    double parMag = fitf->GetParameter(3);
    double halfPointY = parMag/2 + parOffset;

    double risingEdgeTick = fitf->GetX(halfPointY, parA - 10, parA + 10);
    double risingEdgeNs = std::round(wf->TimeStamp()*1'000 + risingEdgeTick*2.0);

    return risingEdgeNs;
}

DEFINE_ART_MODULE(sbnd::BeamAnalysis)
