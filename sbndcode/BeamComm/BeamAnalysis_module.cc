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
#include "lardataobj/RecoBase/OpFlash.h"

#include "sbndcode/Timing/SBNDRawTimingObj.h"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbndcode/ChannelMaps/PMT/PMTChannelMapService.h"


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
    TTree *fTree;
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

    // CRT stuff
    std::vector<double> _crt_x;
    std::vector<double> _crt_y;
    std::vector<double> _crt_z;
    std::vector<double> _crt_xe;
    std::vector<double> _crt_ye;
    std::vector<double> _crt_ze;
    std::vector<double> _crt_ts0;
    std::vector<double> _crt_ts1;
    std::vector<double> _crt_ts0e;
    std::vector<double> _crt_ts1e;
    std::vector<int> _crt_tagger;

    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;
    std::map<uint16_t, double> board_jitter;

    //OpHit
    int _nophits;
    std::vector<int>  _ophit_opch;           ///< OpChannel of the optical hit
    std::vector<double> _ophit_peakT;       ///< Peak time of the optical hit
    std::vector<double> _ophit_startT;       ///< Peak time of the optical hit
    std::vector<double> _ophit_riseT;       ///< Peak time of the optical hit
    std::vector<double> _ophit_width;       ///< Width of the optical hit
    std::vector<double> _ophit_area;        ///< Area of the optical hit
    std::vector<double> _ophit_amplitude;   ///< Amplitude of the optical hit
    std::vector<double> _ophit_pe;          ///< PEs of the optical hit
    std::vector<double> _ophit_totalTransit;  ///< Cable delay from PMT to digitiser
    std::vector<double> _ophit_boardJitter;  ///< Board jittering deduced from FTRIG

    // OpFlash
    int _nopflash;
    std::vector<int> _flash_id;
    std::vector<double> _flash_time;
    std::vector<double> _flash_total_pe;
    std::vector<std::vector<double>> _flash_pe_v;
    std::vector<double> _flash_y;
    std::vector<double> _flash_yerr ;
    std::vector<double> _flash_z;
    std::vector<double> _flash_zerr;
    std::vector<double> _flash_x;
    std::vector<double> _flash_xerr;
    std::vector<int> _flash_tpc;
    std::vector<std::vector<double>> _flash_ophit_time;
    std::vector<std::vector<double>> _flash_ophit_risetime;
    std::vector<std::vector<double>> _flash_ophit_starttime;
    std::vector<std::vector<double>>_flash_ophit_amp;
    std::vector<std::vector<double>> _flash_ophit_area;
    std::vector<std::vector<double>> _flash_ophit_width;
    std::vector<std::vector<double>> _flash_ophit_pe;
    std::vector<std::vector<int>> _flash_ophit_ch;
    std::vector<std::vector<double>> _flash_ophit_totalTransit;
    std::vector<std::vector<double>> _flash_ophit_boardJitter;

    //Service
    art::ServiceHandle<SBND::PMTChannelMapService> fPMTChannelMapService;
    
    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fCrtSpacePointLabel;
    art::InputTag fPmtFtrigDecodeLabel;
    art::InputTag fPmtFtrigBoardLabel;
    art::InputTag fPmtTimingLabel;
    art::InputTag fOpHitLabel;
    std::vector<art::InputTag> fOpFlashLabels;

    // Debug
    bool fDebugTdc;
    bool fDebugCrt;
    bool fDebugPmt;
    
    // Which data products to include
    bool fIncludeCrt;
    bool fIncludePmt;

    std::vector<int> fExcludeFtrigBoard;
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

    fOpHitLabel = p.get<art::InputTag>("OpHitLabel","ophitpmt");
    fOpFlashLabels = p.get<std::vector<art::InputTag>>("OpFlashLabel", {"opflashtpc0","opflashtpc1"});

    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugCrt = p.get<bool>("DebugCrt", false);
    fDebugPmt = p.get<bool>("DebugPmt", false);

    fIncludeCrt = p.get<bool>("IncludeCrt", true);
    fIncludePmt = p.get<bool>("IncludePmt", true);

    fExcludeFtrigBoard = p.get<std::vector<int>>("ExcludeFtrigBoard", {});
    fSavePmt = p.get<bool>("SavePmt", false);
    fSavePmtPath = p.get<std::string>("SavePmtPath", "");
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

                _crt_x.push_back(crt_sp->X());
                _crt_y.push_back(crt_sp->Y());
                _crt_z.push_back(crt_sp->Z());
                _crt_xe.push_back(crt_sp->XErr());
                _crt_ye.push_back(crt_sp->YErr());
                _crt_ze.push_back(crt_sp->ZErr());
                _crt_ts0.push_back(crt_sp->Ts0());
                _crt_ts1.push_back(crt_sp->Ts1());
                _crt_ts0e.push_back(crt_sp->Ts0Err());
                _crt_ts1e.push_back(crt_sp->Ts1Err());
                _crt_tagger.push_back(crt_cluster->Tagger());

                if (fDebugCrt){
                    std::cout << "CRT Space Point------------------------------------" << std::endl;
                    std::cout << "   x = " << _crt_x.back() << ", y = " << _crt_y.back() << ", z = " << _crt_z.back() << std::endl;
                    std::cout << "   ts0 = " << _crt_ts0.back() << ", ts1 = " << _crt_ts1.back() << std::endl;
                    std::cout << "   tagger = " << _crt_tagger.back() << std::endl;
                }
                
                double ts0_ns = _crt_ts0.back() - (_tdc_ch2.back() - _tdc_ch4.back());
                double ts0_us = ts0_ns / 1'000;

                _hTopHatCRTT0[_crt_tagger.back()]->Fill(ts0_us);
            }
        }
    }

    //-----------------------------------------------------------------//

    if (fIncludePmt){

        //------------------------PMT Timing--------------------------//
        art::Handle<raw::TimingReferenceInfo> timingRefHandle;
        e.getByLabel(fPmtTimingLabel, timingRefHandle);

        if (!timingRefHandle.isValid()){
            if (fDebugPmt) std::cout << "No Timing Reference products found." << std::endl;
        }
        else{
            raw::TimingReferenceInfo const& pmt_timing(*timingRefHandle);
            
            _pmt_timing_type = pmt_timing.timingType;
            _pmt_timing_ch = pmt_timing.timingChannel;

            if (fDebugPmt){
                std::cout << "Timing Reference For Decoding PMT" << std::endl;
                std::cout << "   Type = " << _pmt_timing_type << " (SPECTDC = 0; PTB HLT = 1; CAEN-only = 3)." << std::endl;
                std::cout << "   Channel = " << _pmt_timing_ch << " (TDC ETRIG = 4; PTB BNB Beam+Light = 2)." << std::endl;
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
            art::FindManyP<raw::pmt::BoardTimingInfo> pmtBoardAssoc(pmt_ftrig_v, e, fPmtFtrigBoardLabel);

            if (fDebugPmt) std::cout << "Found OpDetWaveform FTRIG size = " << pmt_ftrig_v.size() << std::endl;
                
            for (auto const& pmt: pmt_ftrig_v){

                const std::vector<art::Ptr<raw::pmt::BoardTimingInfo>> pmt_board_v(pmtBoardAssoc.at(pmt.key()));

                if(pmt_board_v.size() != 1 ) continue;

                const art::Ptr<raw::pmt::BoardTimingInfo> pmt_board(pmt_board_v.front());
              
                bool excludeThisBoard = false;
                for(auto const ch: fExcludeFtrigBoard){
                    if (int(pmt->ChannelNumber())==ch) excludeThisBoard = true;
                }

                if (excludeThisBoard){
                    //FTRIG channel ID is actually digitiser board ID
                    if(fDebugPmt) std::cout << "   board id = " << pmt->ChannelNumber() << " is exluded." << std::endl << std::endl;
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

            if (fDebugPmt) std::cout << "Found OpHit size = " << ophit_v.size() << std::endl;
      
            _nophits = ophit_v.size();

            for (auto const& ophit: ophit_v){

                 _ophit_opch.push_back( ophit->OpChannel() );
                 _ophit_peakT.push_back( ophit->PeakTimeAbs() );
                 _ophit_startT.push_back( ophit->StartTime() );
                 _ophit_riseT.push_back( ophit->RiseTime() );
                 _ophit_width.push_back( ophit->Width() );
                 _ophit_area.push_back( ophit->Area() );
                 _ophit_amplitude.push_back(  ophit->Amplitude() );
                 _ophit_pe.push_back( ophit->PE() );

                //get PMT Info
                SBND::PMTChannelMapService::PMTInfo_t pmtInfo = fPMTChannelMapService->GetPMTInfoFromChannelID(ophit->OpChannel());

                if(!pmtInfo.valid){
                    if (fDebugPmt) std::cout << "Cannot find PMT Info for channel ID " << ophit->OpChannel() << std::endl;
                    _ophit_totalTransit.push_back(-99999);
                    _ophit_boardJitter.push_back(-99999);
                }else{
                    _ophit_totalTransit.push_back(pmtInfo.TotalTransit);
                    _ophit_boardJitter.push_back(board_jitter[pmtInfo.digitiserBoardID]);
                }
            }
        }
        //---------------------------OpFlash-----------------------------//
        art::Handle<std::vector<recob::OpFlash>> opFlashHandle;
        std::vector<art::Ptr<recob::OpFlash>> opflash_v;

        // Loop over all the OpFlash labels
        for (size_t s = 0; s < fOpFlashLabels.size(); s++) {
            e.getByLabel(fOpFlashLabels[s], opFlashHandle);

            if(!opFlashHandle.isValid() || opFlashHandle->size() == 0){
                if (fDebugPmt) std::cout << "No OpFlash products found with label " << fOpFlashLabels[s] << std::endl;
            }else{
            
                art::fill_ptr_vector(opflash_v, opFlashHandle);
                art::FindManyP<recob::OpHit> flashOpHitAssoc(opflash_v, e, fOpFlashLabels[s]);

                if (fDebugPmt) std::cout << "Found OpFlash with label " << fOpFlashLabels[s] << " size = " << opflash_v.size() << std::endl;

                for (auto const& flash: opflash_v){

                    const std::vector<art::Ptr<recob::OpHit>> ophit_v(flashOpHitAssoc.at(flash.key()));

                    _flash_id.push_back( _nopflash );
                    _flash_time.push_back( flash->AbsTime() );
                    _flash_total_pe.push_back( flash->TotalPE() );
                    _flash_pe_v.push_back( flash->PEs() );
                    _flash_tpc.push_back( s );
                    _flash_y.push_back( flash->YCenter() );
                    _flash_yerr.push_back( flash->YWidth() );
                    _flash_x.push_back( flash->XCenter() );
                    _flash_xerr.push_back( flash->XWidth() );
                    _flash_z.push_back( flash->ZCenter() );
                    _flash_zerr.push_back( flash->ZWidth() );
          
                    _flash_ophit_time.push_back({});
                    _flash_ophit_risetime.push_back({});
                    _flash_ophit_starttime.push_back({});
                    _flash_ophit_amp.push_back({});
                    _flash_ophit_area.push_back({});
                    _flash_ophit_width.push_back({});
                    _flash_ophit_pe.push_back({});
                    _flash_ophit_ch.push_back({});
                    _flash_ophit_totalTransit.push_back({});
                    _flash_ophit_boardJitter.push_back({});

                    if (fDebugPmt) std::cout << "  OpHit associated with this flash size = " << ophit_v.size() << std::endl;
                    
                    for (auto ophit : ophit_v) {
                        _flash_ophit_time[_nopflash].push_back(ophit->PeakTimeAbs());
                        _flash_ophit_risetime[_nopflash].push_back(ophit->RiseTime());
                        _flash_ophit_starttime[_nopflash].push_back(ophit->StartTime());
                        _flash_ophit_amp[_nopflash].push_back(ophit->Amplitude());
                        _flash_ophit_area[_nopflash].push_back(ophit->Area());
                        _flash_ophit_width[_nopflash].push_back(ophit->Width());
                        _flash_ophit_pe[_nopflash].push_back(ophit->PE());
                        _flash_ophit_ch[_nopflash].push_back(ophit->OpChannel());

                        SBND::PMTChannelMapService::PMTInfo_t pmtInfo = fPMTChannelMapService->GetPMTInfoFromChannelID(ophit->OpChannel());

                        if(!pmtInfo.valid){
                            if (fDebugPmt) std::cout << "Cannot find PMT Info for channel ID " << ophit->OpChannel() << std::endl;
                            _flash_ophit_totalTransit[_nopflash].push_back(-99999);
                            _flash_ophit_boardJitter[_nopflash].push_back(-99999);
                        }else{
                            _flash_ophit_totalTransit[_nopflash].push_back(pmtInfo.TotalTransit);
                            _flash_ophit_boardJitter[_nopflash].push_back(board_jitter[pmtInfo.digitiserBoardID]);
                        }
                    }

                    _nopflash++;
                } 
            }
        }
    }
    //-----------------------------------------------------------//
    
    if (fDebugTdc | fDebugCrt | fDebugPmt) std::cout <<"#--------------------------------------------------------#" << std::endl;

    //Fill once every event
    fTree->Fill();
}

void sbnd::BeamAnalysis::beginJob()
{
    for(size_t i = 0; i < nCrt; i++){
        _hTopHatCRTT0[i] = tfs->make<TH1D>(Form("hCRTT0_Tagger_%s", crtTagger[i].c_str()), "", 100, -5, 5);
    }

    for(size_t i = 0; i < nPmt; i++){
        _hFTRIG[i] = tfs->make<TH1D>(Form("hFTRIG_RisingEdge_Board_%s", pmtBoard[i].c_str()), "", 40, -120, -80);
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

    if (fIncludeCrt){
        fTree->Branch("crt_x", &_crt_x);
        fTree->Branch("crt_y", &_crt_y);
        fTree->Branch("crt_z", &_crt_z);
        fTree->Branch("crt_xe", &_crt_xe);
        fTree->Branch("crt_ye", &_crt_ye);
        fTree->Branch("crt_ze", &_crt_ze);
        fTree->Branch("crt_ts0", &_crt_ts0);
        fTree->Branch("crt_ts1", &_crt_ts1);
        fTree->Branch("crt_ts0e", &_crt_ts0e);
        fTree->Branch("crt_ts1e", &_crt_ts1e);
        fTree->Branch("crt_tagger", &_crt_tagger);
    }
    
    if (fIncludePmt){
        fTree->Branch("_pmt_timing_type", &_pmt_timing_type);
        fTree->Branch("_pmt_timing_ch", &_pmt_timing_ch);

        fTree->Branch("nophits", &_nophits);
        fTree->Branch("ophit_opch", &_ophit_opch);
        fTree->Branch("ophit_peakT", &_ophit_peakT);
        fTree->Branch("ophit_startT", &_ophit_startT);
        fTree->Branch("ophit_riseT", &_ophit_riseT);
        fTree->Branch("ophit_width", &_ophit_width);
        fTree->Branch("ophit_area", &_ophit_area);
        fTree->Branch("ophit_amplitude", &_ophit_amplitude);
        fTree->Branch("ophit_pe", &_ophit_pe);
        fTree->Branch("ophit_totalTransit", &_ophit_totalTransit);
        fTree->Branch("ophit_boardJitter", &_ophit_boardJitter);

        fTree->Branch("nopflash", &_nopflash);
        fTree->Branch("flash_id", &_flash_id);
        fTree->Branch("flash_time", &_flash_time);
        fTree->Branch("flash_total_pe", &_flash_total_pe);
        fTree->Branch("flash_pe_v", &_flash_pe_v);
        fTree->Branch("flash_y",&_flash_y);
        fTree->Branch("flash_yerr",&_flash_yerr);
        fTree->Branch("flash_z",&_flash_z);
        fTree->Branch("flash_zerr",&_flash_zerr);
        fTree->Branch("flash_x",&_flash_x);
        fTree->Branch("flash_xerr",&_flash_xerr);
        fTree->Branch("flash_tpc",&_flash_tpc);
        fTree->Branch("flash_ophit_time",&_flash_ophit_time);
        fTree->Branch("flash_ophit_risetime",&_flash_ophit_risetime);
        fTree->Branch("flash_ophit_starttime",&_flash_ophit_starttime);
        fTree->Branch("flash_ophit_amp",&_flash_ophit_amp);
        fTree->Branch("flash_ophit_area",&_flash_ophit_area);
        fTree->Branch("flash_ophit_width",&_flash_ophit_width);
        fTree->Branch("flash_ophit_pe",&_flash_ophit_pe);
        fTree->Branch("flash_ophit_ch",&_flash_ophit_ch);
        fTree->Branch("flash_ophit_totalTransit",&_flash_ophit_totalTransit);
        fTree->Branch("flash_ophit_boardJitter",&_flash_ophit_boardJitter);
    }
}

void sbnd::BeamAnalysis::endJob()
{
}

void sbnd::BeamAnalysis::ResetEventVars()
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

    if (fIncludeCrt){
        _crt_x.clear();
        _crt_y.clear();
        _crt_z.clear();
        _crt_xe.clear();
        _crt_ye.clear();
        _crt_ze.clear();
        _crt_ts0.clear();
        _crt_ts1.clear();
        _crt_ts0e.clear();
        _crt_ts1e.clear();
        _crt_tagger.clear();
    }


    if (fIncludePmt){
        _pmt_timing_type = -1;
        _pmt_timing_ch = -1;

        for (unsigned int i = 0; i < nPmt; i++){
            board_jitter[i] = -99999;
        } 

        _nophits = 0;
        _ophit_opch.clear();
        _ophit_peakT.clear();
        _ophit_startT.clear();
        _ophit_riseT.clear();
        _ophit_width.clear();
        _ophit_area.clear();
        _ophit_amplitude.clear();
        _ophit_pe.clear();
        _ophit_totalTransit.clear();
        _ophit_boardJitter.clear();

        _nopflash = 0;
        _flash_id.clear();
        _flash_time.clear();
        _flash_total_pe.clear();
        _flash_pe_v.clear();
        _flash_y.clear();
        _flash_yerr.clear();
        _flash_z.clear();
        _flash_zerr.clear();
        _flash_x.clear();
        _flash_xerr.clear();
        _flash_tpc.clear();
        _flash_ophit_time.clear();
        _flash_ophit_risetime.clear();
        _flash_ophit_starttime.clear();
        _flash_ophit_amp.clear();
        _flash_ophit_area.clear();
        _flash_ophit_width.clear();
        _flash_ophit_pe.clear();
        _flash_ophit_ch.clear();
        _flash_ophit_totalTransit.clear();
        _flash_ophit_boardJitter.clear();

    }
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
    double postPercent = 1 - (pp / 100.0);
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

    //-------Check if fit converges
    if(!converged) {
        if(fDebugPmt) std::cout << "Fitting FTRIG waveform does not converge!" << std::endl;
        return -99999;
    }

    //-------Save fits to plots
    if(fSavePmt){
        
        gStyle->SetOptStat(0); 
        gStyle->SetOptFit(1);
    
        TCanvas *c = new TCanvas(_histName.str().c_str(), _histName.str().c_str());
        c->cd();

        hWvfm->GetXaxis()->SetTitle("Tick");
        hWvfm->GetYaxis()->SetTitle("ADC");
        //hWvfm->GetXaxis()->SetRangeUser(wf_lb, wf_ub);
        
        hWvfm->SetLineColor(kBlack);
        hWvfm->SetLineWidth(2);
        hWvfm->Draw("histe");
        fitf->SetLineColor(kSpring-6);
        fitf->Draw("same");
        
        c->SaveAs(Form("%s/%s.png", fSavePmtPath.c_str(), _histName.str().c_str()));
        //c->SaveAs(Form("%s/%s.pdf", fSavePmtPath.c_str(), _histName.str().c_str()));
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
