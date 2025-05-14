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
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"

#include <set>
#include <vector>


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
    TTree *fTree;
    std::stringstream _histName; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;

    std::vector<std::string> crtTagger{"Bottom","South","North","West","East","TopLow","TopHigh"};
    size_t nCrt = crtTagger.size(); 
    
    std::vector<std::string> pmtBoard{"0","1","2","3","4","5","6","7","8"};
    size_t nPmt = pmtBoard.size(); 

    std::vector<int> hlt{1, 2, 3, 4, 5, 14, 15};

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
    std::vector<uint64_t> _ptb_hlt_unmask_timestamp;
    std::vector<int> _ptb_hlt_trunmask;
      
    // CRT spacepoint
    std::vector<double> _crt_sp_x;
    std::vector<double> _crt_sp_y;
    std::vector<double> _crt_sp_z;
    std::vector<double> _crt_sp_xe;
    std::vector<double> _crt_sp_ye;
    std::vector<double> _crt_sp_ze;
    std::vector<double> _crt_sp_ts0;
    std::vector<double> _crt_sp_ts1;
    std::vector<double> _crt_sp_ts0e;
    std::vector<double> _crt_sp_ts1e;
    std::vector<int> _crt_sp_tagger;

    // CRT track
    int _ncrt_trk;
    std::vector<int>    _crt_trk_id;
    std::vector<double> _crt_trk_ts0;     // average time according to T0 clock [ns]
    std::vector<double> _crt_trk_ts0Err;  // average time error according to T0 clock [ns]
    std::vector<double> _crt_trk_ts1;     // average time according to T1 clock [ns]
    std::vector<double> _crt_trk_ts1Err;  // average time error according to T1 clock [ns]
    std::vector<double> _crt_trk_pe;      // total PE
    std::vector<double> _crt_trk_tof;     // time from first space point to last [ns]
    std::vector<std::vector<int>> _crt_trk_taggers; // which taggers were used to create the track
    std::vector<double> _crt_trk_theta;
    std::vector<double> _crt_trk_phi;
    std::vector<double> _crt_trk_length;
    std::vector<double> _crt_trk_start_x;
    std::vector<double> _crt_trk_start_y;
    std::vector<double> _crt_trk_start_z;
    std::vector<double> _crt_trk_end_x;
    std::vector<double> _crt_trk_end_y;
    std::vector<double> _crt_trk_end_z;
    std::vector<std::vector<double>> _crt_trk_sp_x;
    std::vector<std::vector<double>> _crt_trk_sp_y;
    std::vector<std::vector<double>> _crt_trk_sp_z;
    std::vector<std::vector<double>> _crt_trk_sp_xe;
    std::vector<std::vector<double>> _crt_trk_sp_ye;
    std::vector<std::vector<double>> _crt_trk_sp_ze;
    std::vector<std::vector<double>> _crt_trk_sp_ts0;
    std::vector<std::vector<double>> _crt_trk_sp_ts1;
    std::vector<std::vector<double>> _crt_trk_sp_ts0e;
    std::vector<std::vector<double>> _crt_trk_sp_ts1e;
    std::vector<std::vector<int>>    _crt_trk_sp_tagger;

    // PMT Timing
    uint16_t _pmt_timing_type;
    uint16_t _pmt_timing_ch;

    // PMT Metric
    float _metric_peakpe;
    float _metric_peaktime;
    float _metric_dt;

    //OpHit
    //int _nophits;
    //std::vector<int>  _ophit_opch;           ///< OpChannel of the optical hit
    //std::vector<double> _ophit_peakT;       ///< Peak time of the optical hit
    //std::vector<double> _ophit_startT;       ///< Peak time of the optical hit
    //std::vector<double> _ophit_riseT;       ///< Peak time of the optical hit
    //std::vector<double> _ophit_width;       ///< Width of the optical hit
    //std::vector<double> _ophit_area;        ///< Area of the optical hit
    //std::vector<double> _ophit_amplitude;   ///< Amplitude of the optical hit
    //std::vector<double> _ophit_pe;          ///< PEs of the optical hit

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

    //---FHICL CONFIG PARAMETERS
    
    // Product label
    art::InputTag fTdcDecodeLabel;
    art::InputTag fPtbDecodeLabel;
    art::InputTag fCrtSpacePointLabel;
    art::InputTag fCrtTrackLabel;
    art::InputTag fPmtTimingLabel;\
    art::InputTag fOpHitLabel;
    std::vector<art::InputTag> fOpFlashLabels;
    std::vector<art::InputTag> fPmtMetricLabels;

    // Debug
    bool fDebugTdc;
    bool fDebugPtb;
    bool fDebugCrtSP;
    bool fDebugCrtTrack;
    bool fDebugPmt;
    
    // Which data products to include
    bool fIncludeCrtSP;
    bool fIncludePtb;
    bool fIncludeCrtTrack;
    bool fIncludePmtMetric;
    bool fIncludePmt;
};


sbnd::BeamAnalysis::BeamAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // 
    // More initializers here.
{
    fTdcDecodeLabel = p.get<art::InputTag>("TdcDecodeLabel", "tdcdecoder");
    fPtbDecodeLabel = p.get<art::InputTag>("PtbDecodeLabel", "ptbdecoder");
    fCrtSpacePointLabel = p.get<art::InputTag>("CrtSpacePointLabel", "crtspacepoints");
    fCrtTrackLabel = p.get<art::InputTag>("CrtTrackLabel", "crttracks");
    fPmtTimingLabel = p.get<art::InputTag>("fPmtTimingLabel", "pmtdecoder");
    fOpHitLabel = p.get<art::InputTag>("OpHitLabel","ophitpmt");
    fOpFlashLabels = p.get<std::vector<art::InputTag>>("OpFlashLabel", {"opflashtpc0","opflashtpc1"});
    fPmtMetricLabels = p.get<std::vector<art::InputTag>>("PmtMetricLabels", {"pmtmetriccrossingmuon", "pmtmetricbnbzero", "pmtmetricbnblight", "pmtmetricoffbeamzero", "pmtmetricoffbeamlight"});
    fDebugTdc = p.get<bool>("DebugTdc", false);
    fDebugPtb = p.get<bool>("DebugPtb", false);
    fDebugCrtSP = p.get<bool>("DebugCrtSP", false);
    fDebugCrtTrack = p.get<bool>("DebugCrtTrack", false);
    fDebugPmt = p.get<bool>("DebugPmt", false);

    fIncludePtb = p.get<bool>("IncludePtb", true);
    fIncludeCrtSP = p.get<bool>("IncludeCrtSP", false);
    fIncludeCrtTrack = p.get<bool>("IncludeCrtTrack", false);
    fIncludePmtMetric = p.get<bool>("IncludePmtMetric", false);
    fIncludePmt = p.get<bool>("IncludePmt", false);
}

void sbnd::BeamAnalysis::analyze(art::Event const& e)
{
    
    ResetEventVars();
  
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  =  e.id().event();

    if (fDebugTdc | fDebugCrtSP | fDebugCrtTrack | fDebugPmt)
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
            _ptb_hlt_unmask_timestamp.resize(nHLTs);

            unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
            unsigned h_i = 0; //For trigger with bitmask
            for(auto const& ptb : ptb_v){
                for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i){
                    _ptb_hlt_trigger[h_i] = ptb->GetHLTrigger(i).trigger_word;
                    _ptb_hlt_timestamp[h_i] = ptb->GetHLTrigger(i).timestamp*20; //Units can be found in the Decoder Module 
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

            if (fDebugPtb){
                for (size_t i = 0; i <_ptb_hlt_trunmask.size(); i++){
                    for(auto h: hlt){
                        if (_ptb_hlt_trunmask[i] == h) std::cout << "This stream is hlt " << h << std::endl;
                    }

                    //std::cout << " HLT " << _ptb_hlt_trunmask[i]
                    //        << ", ts = " << _ptb_hlt_unmask_timestamp[i] << std::endl;
                }
            }
        }
    }

    //---------------------------CRT-----------------------------//
   
    if (fIncludeCrtSP){

        art::Handle<std::vector<sbnd::crt::CRTSpacePoint>> crtSpacePointHandle;
        std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_sp_v;
        e.getByLabel(fCrtSpacePointLabel, crtSpacePointHandle);

        if (!crtSpacePointHandle.isValid() || crtSpacePointHandle->size() == 0){
            if (fDebugCrtSP) std::cout << "No CRT Space Point products found." << std::endl;
        }
        else{

            art::fill_ptr_vector(crt_sp_v, crtSpacePointHandle);
            art::FindManyP<sbnd::crt::CRTCluster> crtSPClusterAssoc(crt_sp_v, e, fCrtSpacePointLabel);

            for (auto const& crt_sp: crt_sp_v){
                
                if (!crt_sp->Complete()) continue;

                const std::vector<art::Ptr<sbnd::crt::CRTCluster>> crt_cluster_v(crtSPClusterAssoc.at(crt_sp.key()));

                if(crt_cluster_v.size() != 1 ) continue;

                const art::Ptr<sbnd::crt::CRTCluster>& crt_cluster(crt_cluster_v.front());

                _crt_sp_x.push_back(crt_sp->X());
                _crt_sp_y.push_back(crt_sp->Y());
                _crt_sp_z.push_back(crt_sp->Z());
                _crt_sp_xe.push_back(crt_sp->XErr());
                _crt_sp_ye.push_back(crt_sp->YErr());
                _crt_sp_ze.push_back(crt_sp->ZErr());
                _crt_sp_ts0.push_back(crt_sp->Ts0());
                _crt_sp_ts1.push_back(crt_sp->Ts1());
                _crt_sp_ts0e.push_back(crt_sp->Ts0Err());
                _crt_sp_ts1e.push_back(crt_sp->Ts1Err());
                _crt_sp_tagger.push_back(crt_cluster->Tagger());

                if (fDebugCrtSP){
                    std::cout << "CRT Space Point------------------------------------" << std::endl;
                    std::cout << "   x = " << _crt_sp_x.back() << ", y = " << _crt_sp_y.back() << ", z = " << _crt_sp_z.back() << std::endl;
                    std::cout << "   ts0 = " << _crt_sp_ts0.back() << ", ts1 = " << _crt_sp_ts1.back() << std::endl;
                    std::cout << "   tagger = " << _crt_sp_tagger.back() << std::endl;
            
                }
            }
        }
    }

    if (fIncludeCrtTrack){

        art::Handle<std::vector<sbnd::crt::CRTTrack>> crtTrackHandle;
        std::vector<art::Ptr<sbnd::crt::CRTTrack>> crt_trk_v;
        e.getByLabel(fCrtTrackLabel, crtTrackHandle);

        art::Handle<std::vector<sbnd::crt::CRTSpacePoint>> crtSpacePointHandle;
        std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_sp_v;
        e.getByLabel(fCrtSpacePointLabel, crtSpacePointHandle);
        
        std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_trk_sp_v;

        if (!crtTrackHandle.isValid() || crtTrackHandle->size() == 0){
            if (fDebugCrtTrack) std::cout << "No CRT Track products found." << std::endl;
        }
        else{

            art::fill_ptr_vector(crt_trk_v, crtTrackHandle);
            art::FindManyP<sbnd::crt::CRTSpacePoint> crtTrkSPAssoc(crt_trk_v, e, fCrtTrackLabel);
            
            art::fill_ptr_vector(crt_sp_v, crtSpacePointHandle);
            art::FindManyP<sbnd::crt::CRTCluster> crtSPClusterAssoc(crt_sp_v, e, fCrtSpacePointLabel);

            for (auto const& crt_trk: crt_trk_v){

                _crt_trk_id.push_back( _ncrt_trk );
                _crt_trk_ts0.push_back(crt_trk->Ts0());
                _crt_trk_ts0Err.push_back(crt_trk->Ts0Err());
                _crt_trk_ts1.push_back(crt_trk->Ts1());  
                _crt_trk_ts1Err.push_back(crt_trk->Ts1Err()); 
                _crt_trk_pe.push_back(crt_trk->PE());     
                _crt_trk_tof.push_back(crt_trk->ToF());    
                _crt_trk_theta.push_back(crt_trk->Theta());
                _crt_trk_phi.push_back(crt_trk->Phi());
                std::vector<int> taggers;
                for(auto const t: crt_trk->Taggers())
                    taggers.push_back(t);

                _crt_trk_taggers.push_back(taggers);

                _crt_trk_length.push_back(crt_trk->Length());
                _crt_trk_start_x.push_back(crt_trk->Start().X());
                _crt_trk_start_y.push_back(crt_trk->Start().Y());
                _crt_trk_start_z.push_back(crt_trk->Start().Z());
                _crt_trk_end_x.push_back(crt_trk->End().X());
                _crt_trk_end_y.push_back(crt_trk->End().Y());
                _crt_trk_end_z.push_back(crt_trk->End().Z());

                _crt_trk_sp_x.push_back({});
                _crt_trk_sp_y.push_back({});
                _crt_trk_sp_z.push_back({});
                _crt_trk_sp_xe.push_back({});
                _crt_trk_sp_ye.push_back({});
                _crt_trk_sp_ze.push_back({});
                _crt_trk_sp_ts0.push_back({});
                _crt_trk_sp_ts1.push_back({});
                _crt_trk_sp_ts0e.push_back({});
                _crt_trk_sp_ts1e.push_back({});
                _crt_trk_sp_tagger.push_back({});

                const std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_trk_sp_v(crtTrkSPAssoc.at(crt_trk.key()));

                if (fDebugCrtTrack) std::cout << "Found associated CRT Space Point size = " << crt_sp_v.size() << std::endl;

                for (auto const& crt_sp: crt_trk_sp_v){

                    const std::vector<art::Ptr<sbnd::crt::CRTCluster>> crt_cluster_v(crtSPClusterAssoc.at(crt_sp.key()));
                    const art::Ptr<sbnd::crt::CRTCluster>& crt_cluster(crt_cluster_v.front());

                    _crt_trk_sp_x[_ncrt_trk].push_back(crt_sp->X());
                    _crt_trk_sp_y[_ncrt_trk].push_back(crt_sp->Y());
                    _crt_trk_sp_z[_ncrt_trk].push_back(crt_sp->Z());
                    _crt_trk_sp_xe[_ncrt_trk].push_back(crt_sp->XErr());
                    _crt_trk_sp_ye[_ncrt_trk].push_back(crt_sp->YErr());
                    _crt_trk_sp_ze[_ncrt_trk].push_back(crt_sp->ZErr());
                    _crt_trk_sp_ts0[_ncrt_trk].push_back(crt_sp->Ts0());
                    _crt_trk_sp_ts1[_ncrt_trk].push_back(crt_sp->Ts1());
                    _crt_trk_sp_ts0e[_ncrt_trk].push_back(crt_sp->Ts0Err());
                    _crt_trk_sp_ts1e[_ncrt_trk].push_back(crt_sp->Ts1Err());
                    _crt_trk_sp_tagger[_ncrt_trk].push_back(crt_cluster->Tagger());
                }

                _ncrt_trk++;

                if (fDebugCrtTrack){
                    std::cout << "CRT Track------------------------------------" << std::endl;
                    //std::cout << std::setprecision(9) << "   ts0 = " << _crt_trk_ts0.back() << ", ts1 = " << _crt_trk_ts1.back() << std::endl;
                    std::cout << "  theta = " << _crt_trk_theta.back() << ", phi = " << _crt_trk_phi.back() << std::endl;
                    //std::cout << "   tagger = ";
                    for (auto const t: _crt_trk_taggers.back()){
                        std::cout << crtTagger[t] << " ";
                    }
                    std::cout << std::endl;
                    
                    auto start = crt_trk->Start();
                    auto end = crt_trk->End();
                    std::cout << " Length = " << crt_trk->Length() << std::endl;
                    std::cout << " Start x = " << start.X()  << ", y = " << start.Y() << ", z = " << start.Z() << std::endl;
                    std::cout << " End x = " << end.X()  << ", y = " << end.Y() << ", z = " << end.Z() << std::endl;
                }
            }
        }
    }

    //-----------------------------------------------------------------//

    if (fIncludePmtMetric){

        //------------------------PMT Metric--------------------------//
        size_t nlabels = fPmtMetricLabels.size();
        for (size_t i=0; i<nlabels; i++){
            auto label = fPmtMetricLabels[i];
      
            art::Handle<std::vector<sbnd::trigger::pmtSoftwareTrigger>> metricHandle;
            e.getByLabel(label,metricHandle);
         
            if (!metricHandle.isValid() || metricHandle->size() == 0) continue;
        
            const std::vector<sbnd::trigger::pmtSoftwareTrigger> metric_v(*metricHandle);
            
            auto metric = metric_v[0];
            _metric_peakpe   = metric.peakPE;
            _metric_peaktime = metric.peaktime; // us
            _metric_dt = float(metric.trig_ts)/1e3;

            break;
        }
    }
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
                std::cout << "   Type = " << _pmt_timing_type << " (SPECTDC = 0; Ptb HLT = 1; CAEN-only = 3)." << std::endl;
                std::cout << "   Channel = " << _pmt_timing_ch << " (TDC ETRIG = 4; Ptb BNB Beam+Light = 2)." << std::endl;
            }
        }

        //---------------------------OpHit-----------------------------//
        //art::Handle<std::vector<recob::OpHit>> opHitHandle;
        //std::vector<art::Ptr<recob::OpHit>> ophit_v;
        //e.getByLabel(fOpHitLabel, opHitHandle);
        //
        //if (!opHitHandle.isValid() || opHitHandle->size() == 0){
        //    if (fDebugPmt) std::cout << "No OpHit products found." << std::endl;
        //}
        //else{
        //    art::fill_ptr_vector(ophit_v, opHitHandle);

        //    if (fDebugPmt) std::cout << "Found OpHit size = " << ophit_v.size() << std::endl;
      
        //    _nophits = ophit_v.size();

        //    for (auto const& ophit: ophit_v){
        //         _ophit_opch.push_back( ophit->OpChannel() );
        //         _ophit_peakT.push_back( ophit->PeakTimeAbs() );
        //         _ophit_startT.push_back( ophit->StartTime() );
        //         _ophit_riseT.push_back( ophit->RiseTime() );
        //         _ophit_width.push_back( ophit->Width() );
        //         _ophit_area.push_back( ophit->Area() );
        //         _ophit_amplitude.push_back(  ophit->Amplitude() );
        //         _ophit_pe.push_back( ophit->PE() );
        //    }
        //}
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
                    }

                    _nopflash++;
                } 
            }
        }
    }
    //-----------------------------------------------------------//
    
    if (fDebugTdc | fDebugCrtSP | fDebugCrtTrack | fDebugPmt) std::cout <<"#--------------------------------------------------------#" << std::endl;

    //Fill once every event
    fTree->Fill();
}

void sbnd::BeamAnalysis::beginJob()
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
    fTree->Branch("tdc_ch0_utc", &_tdc_ch0_utc);
    fTree->Branch("tdc_ch1_utc", &_tdc_ch1_utc);
    fTree->Branch("tdc_ch2_utc", &_tdc_ch2_utc);
    fTree->Branch("tdc_ch3_utc", &_tdc_ch3_utc);
    fTree->Branch("tdc_ch4_utc", &_tdc_ch4_utc);

    if (fIncludePtb){
        fTree->Branch("ptb_hlt_trigger", &_ptb_hlt_trigger);
        fTree->Branch("ptb_hlt_timestamp", &_ptb_hlt_timestamp);
        fTree->Branch("ptb_hlt_trunmask", &_ptb_hlt_trunmask);
        fTree->Branch("ptb_hlt_unmask_timestamp", &_ptb_hlt_unmask_timestamp);
    }

    if (fIncludeCrtSP){
        fTree->Branch("crt_sp_x", &_crt_sp_x);
        fTree->Branch("crt_sp_y", &_crt_sp_y);
        fTree->Branch("crt_sp_z", &_crt_sp_z);
        fTree->Branch("crt_sp_xe", &_crt_sp_xe);
        fTree->Branch("crt_sp_ye", &_crt_sp_ye);
        fTree->Branch("crt_sp_ze", &_crt_sp_ze);
        fTree->Branch("crt_sp_ts0", &_crt_sp_ts0);
        fTree->Branch("crt_sp_ts1", &_crt_sp_ts1);
        fTree->Branch("crt_sp_ts0e", &_crt_sp_ts0e);
        fTree->Branch("crt_sp_ts1e", &_crt_sp_ts1e);
        fTree->Branch("crt_sp_tagger", &_crt_sp_tagger);
    }
    
    if (fIncludeCrtTrack){
        fTree->Branch("ncrt_trk", &_ncrt_trk);
        fTree->Branch("crt_trk_id", &_crt_trk_id);
        fTree->Branch("crt_trk_ts0", &_crt_trk_ts0); 
        fTree->Branch("crt_trk_ts0Err", &_crt_trk_ts0Err); 
        fTree->Branch("crt_trk_ts1", &_crt_trk_ts1); 
        fTree->Branch("crt_trk_ts1Err", &_crt_trk_ts1Err); 
        fTree->Branch("crt_trk_pe", &_crt_trk_pe); 
        fTree->Branch("crt_trk_tof", &_crt_trk_tof); 
        fTree->Branch("crt_trk_taggers", &_crt_trk_taggers); 
        fTree->Branch("crt_trk_theta", &_crt_trk_theta); 
        fTree->Branch("crt_trk_phi", &_crt_trk_phi); 
        fTree->Branch("crt_trk_length", &_crt_trk_length);
        fTree->Branch("crt_trk_start_x", &_crt_trk_start_x);
        fTree->Branch("crt_trk_start_y", &_crt_trk_start_y);
        fTree->Branch("crt_trk_start_z", &_crt_trk_start_z);
        fTree->Branch("crt_trk_end_x", &_crt_trk_end_x);
        fTree->Branch("crt_trk_end_y", &_crt_trk_end_y);
        fTree->Branch("crt_trk_end_z", &_crt_trk_end_z);
        fTree->Branch("crt_trk_sp_x", &_crt_trk_sp_x);
        fTree->Branch("crt_trk_sp_y", &_crt_trk_sp_y);
        fTree->Branch("crt_trk_sp_z", &_crt_trk_sp_z);
        fTree->Branch("crt_trk_sp_xe", &_crt_trk_sp_xe);
        fTree->Branch("crt_trk_sp_ye", &_crt_trk_sp_ye);
        fTree->Branch("crt_trk_sp_ze", &_crt_trk_sp_ze);
        fTree->Branch("crt_trk_sp_ts0", &_crt_trk_sp_ts0);
        fTree->Branch("crt_trk_sp_ts1", &_crt_trk_sp_ts1);
        fTree->Branch("crt_trk_sp_ts0e", &_crt_trk_sp_ts0e);
        fTree->Branch("crt_trk_sp_ts1e", &_crt_trk_sp_ts1e);
        fTree->Branch("crt_trk_sp_tagger", &_crt_trk_sp_tagger);
    }

    if (fIncludePmtMetric){
        fTree->Branch("metric_peakpe", &_metric_peakpe);
        fTree->Branch("metric_peaktime", &_metric_peaktime);
        fTree->Branch("metric_dt", &_metric_dt);
    }

    if (fIncludePmt){
        fTree->Branch("pmt_timing_type", &_pmt_timing_type);
        fTree->Branch("pmt_timing_ch", &_pmt_timing_ch);

        //fTree->Branch("nophits", &_nophits);
        //fTree->Branch("ophit_opch", &_ophit_opch);
        //fTree->Branch("ophit_peakT", &_ophit_peakT);
        //fTree->Branch("ophit_startT", &_ophit_startT);
        //fTree->Branch("ophit_riseT", &_ophit_riseT);
        //fTree->Branch("ophit_width", &_ophit_width);
        //fTree->Branch("ophit_area", &_ophit_area);
        //fTree->Branch("ophit_amplitude", &_ophit_amplitude);
        //fTree->Branch("ophit_pe", &_ophit_pe);

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

    if (fIncludePtb){
        _ptb_hlt_trigger.clear();
        _ptb_hlt_timestamp.clear();
        _ptb_hlt_unmask_timestamp.clear();
        _ptb_hlt_trunmask.clear();
    }

    if (fIncludeCrtSP){
        _crt_sp_x.clear();
        _crt_sp_y.clear();
        _crt_sp_z.clear();
        _crt_sp_xe.clear();
        _crt_sp_ye.clear();
        _crt_sp_ze.clear();
        _crt_sp_ts0.clear();
        _crt_sp_ts1.clear();
        _crt_sp_ts0e.clear();
        _crt_sp_ts1e.clear();
        _crt_sp_tagger.clear();
    }

    if (fIncludeCrtTrack){
        _ncrt_trk = 0;
        _crt_trk_id.clear();
        _crt_trk_ts0.clear();   
        _crt_trk_ts0Err.clear();
        _crt_trk_ts1.clear();   
        _crt_trk_ts1Err.clear();
        _crt_trk_pe.clear();    
        _crt_trk_tof.clear();   
        _crt_trk_taggers.clear();   
        _crt_trk_theta.clear();
        _crt_trk_phi.clear();
        _crt_trk_length.clear();
        _crt_trk_start_x.clear();
        _crt_trk_start_y.clear();
        _crt_trk_start_z.clear();
        _crt_trk_end_x.clear();
        _crt_trk_end_y.clear();
        _crt_trk_end_z.clear();
        _crt_trk_sp_x.clear();
        _crt_trk_sp_y.clear();
        _crt_trk_sp_z.clear();
        _crt_trk_sp_xe.clear();
        _crt_trk_sp_ye.clear();
        _crt_trk_sp_ze.clear();
        _crt_trk_sp_ts0.clear();
        _crt_trk_sp_ts1.clear();
        _crt_trk_sp_ts0e.clear();
        _crt_trk_sp_ts1e.clear();
        _crt_trk_sp_tagger.clear();
    }

    if (fIncludePmtMetric){
        _metric_peakpe = -1;
        _metric_peaktime = -9999;
        _metric_dt = -9999;
    }

    if (fIncludePmt){
        _pmt_timing_type = -1;
        _pmt_timing_ch = -1;

        //_nophits = 0;
        //_ophit_opch.clear();
        //_ophit_peakT.clear();
        //_ophit_startT.clear();
        //_ophit_riseT.clear();
        //_ophit_width.clear();
        //_ophit_area.clear();
        //_ophit_amplitude.clear();
        //_ophit_pe.clear();

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
    }
}

DEFINE_ART_MODULE(sbnd::BeamAnalysis)
