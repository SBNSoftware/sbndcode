////////////////////////////////////////////////////////////////////////
// Class:       SBNDPMTDecoder
// Plugin Type: producer 
// File:        SBNDPMTDecoder_module.cc
// Author:      Lynn Tung (lynnt@uchicago.edu)
//
// Generated at Wed Aug 23 12:28:48 2023 by Lynn Tung using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "artdaq-core/Data/RawEvent.hh"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"
#include "sbndcode/OpDetSim/TriggerEmulationService.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "art_root_io/TFileService.h"
#include "TH1D.h"
#include "TNtuple.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <bitset>
#include <memory>
#include <string>

namespace sbndaq {
    class SBNDPMTDecoder;
}


class sbndaq::SBNDPMTDecoder : public art::EDProducer {
public:
    explicit SBNDPMTDecoder(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SBNDPMTDecoder(SBNDPMTDecoder const&) = delete;
    SBNDPMTDecoder(SBNDPMTDecoder&&) = delete;
    SBNDPMTDecoder& operator=(SBNDPMTDecoder const&) = delete;
    SBNDPMTDecoder& operator=(SBNDPMTDecoder&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;
    void get_fragments(artdaq::Fragment & frag, std::vector<std::vector<artdaq::Fragment>> & board_frag_v);
    void get_waveforms(artdaq::Fragment & frag, std::vector<std::vector<uint16_t>> & wvfm_v);

    uint32_t get_length(artdaq::Fragment & frag);
    uint32_t get_ttt(artdaq::Fragment & frag);
    uint32_t get_boardid(artdaq::Fragment & frag);
    void     get_timing(artdaq::Fragment & frag, uint16_t & postpercent, uint32_t & ttt, uint32_t & len, int & tick);

    std::vector<uint> fill_chmap(sbndDB::PMTCalibrationDatabase const* pmt_calib_db);

private:
    uint fdebug;

    std::vector<std::string>  fcaen_fragment_name;
    std::string               fcaen_module_label;
    std::string               fframeshift_module_label;
    
    std::vector<uint32_t>     fignore_fragid;
    uint32_t                  fnominal_length; 
    uint                      fallowed_time_diff;

    std::string fpmt_instance_name;
    std::string fflt_instance_name;
    std::string ftim_instance_name;

    std::string fpmt_timing_instance_name;
    std::string fflt_timing_instance_name;
    std::string ftim_timing_instance_name;

    bool foutput_ftrig_wvfm;
    bool foutput_timing_wvfm;
    std::vector<int> fignore_timing_ch;

    int fn_maxflashes;
    uint fn_caenboards;
    uint fn_caenchannels;
    uint ftiming_caen_offset;
    uint16_t fthreshold_ftrig;
    uint16_t fdefault_postpercent; // should be a number between 0 and 100

    uint ffragid_offset;
    uint fhist_evt;

    std::vector<uint> fset_fragid_map;
    bool fuse_set_map;
    int fmon_threshold;

    // histogram info  
    std::stringstream histname; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;
    uint evt_counter = 0;

    sbndDB::PMTCalibrationDatabase const* fpmt_calib_db;
    opdet::sbndPDMapAlg opdetmap; //map for photon detector types
};


sbndaq::SBNDPMTDecoder::SBNDPMTDecoder(fhicl::ParameterSet const& p)
    : EDProducer{p}  // ,
{
    fdebug          = p.get<uint>("debug",0);

    fcaen_fragment_name = p.get<std::vector<std::string>>("caen_fragment_name");
    fcaen_module_label  = p.get<std::string>("caen_module_label","daq");
    fframeshift_module_label = p.get<std::string>("frameshift_module_label");

    fignore_fragid = p.get<std::vector<uint32_t>>("ignore_fragid",{});
    fnominal_length = p.get<uint32_t>("nominal_length",5000);
    fallowed_time_diff = p.get<uint>("allowed_time_diff",3000); // us!!! 

    // configure output labels 
    fpmt_instance_name = p.get<std::string>("pmtInstanceName","PMTChannels");
    fflt_instance_name = p.get<std::string>("ftrigInstanceName","FTrigChannels");
    ftim_instance_name = p.get<std::string>("timingInstanceName","TimingChannels");
    fpmt_timing_instance_name = p.get<std::string>("pmtTimingInstanceName","PMTTiming");
    fflt_timing_instance_name = p.get<std::string>("ftrigTimingInstanceName","FTrigTiming");
    ftim_timing_instance_name = p.get<std::string>("timingTimingInstanceName","TimingTiming");

    foutput_ftrig_wvfm = p.get<bool>("output_ftrig_wvfm",true);
    foutput_timing_wvfm = p.get<bool>("output_timing_wvfm",true);
    fignore_timing_ch = p.get<std::vector<int>>("ignore_timing_ch",{});

    fn_maxflashes    = p.get<int>("n_maxflashes",30);
    fn_caenboards    = p.get<uint>("n_caenboards",8);
    fn_caenchannels  = p.get<uint>("n_caenchannels",15);
    ftiming_caen_offset = p.get<uint>("timing_caen_offset",900);
    fthreshold_ftrig = p.get<uint16_t>("threshold_ftrig",16350);
    fdefault_postpercent = p.get<uint16_t>("default_postpercent",80);
    ffragid_offset   = p.get<uint>("fragid_offset",40960);
    fhist_evt        = p.get<int>("hist_evt",1);

    fset_fragid_map = p.get<std::vector<uint>>("set_fragid_map",{});
    fuse_set_map    = p.get<bool>("use_set_map",false);

    fpmt_calib_db    = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();
    fmon_threshold   = p.get<int>("mon_threshold", 15);
 
    produces< std::vector< raw::OpDetWaveform > >(fpmt_instance_name); 
    produces< std::vector< raw::OpDetWaveform > >(fflt_instance_name);
    produces< std::vector< raw::OpDetWaveform > >(ftim_instance_name);

    produces< raw::TimingReferenceInfo >();
    produces< std::vector< raw::pmt::BoardTimingInfo > >();
    produces< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fpmt_timing_instance_name);
    produces< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(fflt_timing_instance_name);
    produces< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform > >(ftim_timing_instance_name);

    produces< std::vector<int> >("MonPulses");
    produces< std::vector<int> >("MonPulseSizes"); 
}

void sbndaq::SBNDPMTDecoder::produce(art::Event& evt)
{
    // output data products
    std::unique_ptr< std::vector< raw::OpDetWaveform > > pmtwvfmVec (new std::vector< raw::OpDetWaveform >);
    std::unique_ptr< std::vector< raw::OpDetWaveform > > fltwvfmVec (new std::vector< raw::OpDetWaveform >);
    std::unique_ptr< std::vector< raw::OpDetWaveform > > timwvfmVec (new std::vector< raw::OpDetWaveform >);

    std::unique_ptr< raw::TimingReferenceInfo > evtTimingInfo (new raw::TimingReferenceInfo());
    std::unique_ptr< std::vector< raw::pmt::BoardTimingInfo > > brdTimingInfoVec (new std::vector< raw::pmt::BoardTimingInfo >);

    // creating PtrMakers for the associations
    art::PtrMaker<raw::pmt::BoardTimingInfo> make_pmtbrd_ptr{evt};

    art::PtrMaker<raw::OpDetWaveform> make_pmtwvfm_ptr{evt,fpmt_instance_name};
    art::PtrMaker<raw::OpDetWaveform> make_fltwvfm_ptr{evt,fflt_instance_name};
    art::PtrMaker<raw::OpDetWaveform> make_timwvfm_ptr{evt,ftim_instance_name};

    // making the associations
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > pmtTimingAssns (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > fltTimingAssns (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);
    std::unique_ptr< art::Assns<  raw::pmt::BoardTimingInfo, raw::OpDetWaveform> > timTimingAssns (new art::Assns< raw::pmt::BoardTimingInfo, raw::OpDetWaveform >);

    evt_counter++;

    std::vector<std::vector<artdaq::Fragment>> board_frag_v(fn_caenboards);
    std::vector<uint> ch_map = fill_chmap(fpmt_calib_db);

    uint ncont = 0; // counter for number of containers

    bool found_caen = false;
    // get CAEN fragments 
    for (const std::string &caen_name : fcaen_fragment_name){
        art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
        evt.getByLabel(fcaen_module_label, caen_name, fragmentHandle);

        if (!fragmentHandle.isValid() || fragmentHandle->size() == 0){
            if (fdebug>0) std::cout << "No CAEN V1730 fragments with label " << caen_name << " found." << std::endl;
            continue;
        }
        else
            found_caen = true;

        if (fragmentHandle->front().type() == artdaq::Fragment::ContainerFragmentType){
            for (auto cont : *fragmentHandle){
                artdaq::ContainerFragment contf(cont);
                ncont++;

                if (contf.fragment_type() == sbndaq::detail::FragmentType::CAENV1730){
                    if (ncont==1 && fdebug>1)
                        std::cout << "Found " << contf.block_count() << " CAEN fragments (flash triggers) inside the container." << std::endl;
                    for (size_t ii = 0; ii < contf.block_count(); ++ii)
                        get_fragments(*contf[ii].get(),board_frag_v);
                }
            }
        }
        else if (fragmentHandle->front().type() == sbndaq::detail::FragmentType::CAENV1730){
        	if (fdebug>1) std::cout << "Found " << fragmentHandle->size() << " normal CAEN fragments " << std::endl;
            for (size_t ii = 0; ii < fragmentHandle->size(); ++ii){
                auto frag = fragmentHandle->at(ii);
                get_fragments(frag,board_frag_v);
            }
        }
        fragmentHandle.removeProduct();
    }

    if (found_caen==false){
        if (fdebug>0) std::cout << "No CAEN V1730 fragments of any type found, pushing empty waveforms." << std::endl;

        evt.put(std::move(pmtwvfmVec),fpmt_instance_name);  
        evt.put(std::move(fltwvfmVec),fflt_instance_name);
        evt.put(std::move(timwvfmVec),ftim_instance_name);

        evt.put(std::move(evtTimingInfo));
        evt.put(std::move(brdTimingInfoVec));

        evt.put(std::move(pmtTimingAssns),fpmt_timing_instance_name);
        evt.put(std::move(fltTimingAssns),fflt_timing_instance_name);
        evt.put(std::move(timTimingAssns),ftim_timing_instance_name);
 
        auto flatPtr = std::make_unique<std::vector<int>>();
        auto sizesPtr = std::make_unique<std::vector<int>>();
        evt.put(std::move(flatPtr), "MonPulses");
        evt.put(std::move(sizesPtr), "MonPulseSizes");
        return;
    }
    
    uint64_t event_trigger_time = 0; // in ns

    art::Handle<sbnd::timing::FrameShiftInfo> frameShiftHandle;    
    evt.getByLabel(fframeshift_module_label, frameShiftHandle);

    if(!frameShiftHandle.isValid())
        throw std::runtime_error("Frame Shift Info object is invalid, check data quality");
    
    event_trigger_time           = frameShiftHandle->FrameDefault();
    evtTimingInfo->timingType    = frameShiftHandle->TimingTypeDefault();
    evtTimingInfo->timingChannel = frameShiftHandle->TimingChannelDefault(); 

    // if frameshift has no shift, the timestamp will be equal to
    // the start of the waveform (according to CAEN TTT and wvfm length)
    if(evtTimingInfo->timingType == sbnd::timing::kNoShiftType)
        event_trigger_time = 0;

    bool extended_flag = false;

    for (size_t iboard=0; iboard < board_frag_v.size(); iboard++){
        // check to see if there are any extended triggers in this event 
        if (extended_flag==false){
            for (size_t itrig=0; itrig < board_frag_v.at(iboard).size(); itrig++){
                auto frag = board_frag_v.at(iboard).at(itrig);
                auto length = get_length(frag);
                if (length < fnominal_length){
                    extended_flag = true;
                    break;
                }
            }
        }
        if (std::find(fignore_fragid.begin(), fignore_fragid.end(), iboard) != fignore_fragid.end())
            continue;
    }

    // store the waveform itself
    // store the waveform start time (for OpDetWaveform timestamp)
    for (size_t iboard=0; iboard < board_frag_v.size(); iboard++){
        // the vector of fragments from a single digitizer 
        auto frag_v = board_frag_v.at(iboard);
        if (std::find(fignore_fragid.begin(), fignore_fragid.end(), iboard) != fignore_fragid.end())
            continue;
        if (frag_v.empty()) continue;
        auto fragid = iboard;
        // trigger number bookkeeping
        auto trig_counter = 0;
        auto full_found = 0;
        auto extensions_found = 0;
        auto extensions_used  = 0;

        for (size_t itrig=0; itrig < frag_v.size(); itrig++){
            raw::pmt::BoardTimingInfo board_info;
            std::vector<uint32_t> board_ttt_v; 
            std::vector<std::vector<uint16_t>> iwvfm_v;
            auto ifrag = frag_v.at(itrig);
            uint16_t ipostpercent = 0; // nanosecond portion of the timestamp (at flash trigger)
            uint32_t ittt = 0; // trigger time tag in ns (at end of the waveform)
            uint32_t ilen = 0;
            int      itick = 0;
            
            get_timing(ifrag, ipostpercent, ittt, ilen, itick);
            // if this waveform is short, skip it 
            if (ilen < fnominal_length){
                extensions_found++;
                continue;
            }
            full_found++;
            board_info.postPercent = ipostpercent;
            board_ttt_v.push_back(ittt);

            get_waveforms(ifrag, iwvfm_v);
            uint iwvfm_start = 0;
            auto iwvfm_end   = ittt;

            // need to check if there is rollover between ttt and start
            if (ittt < ilen*2)
                iwvfm_start = 1e9 - (ilen*2 - ittt);
            else iwvfm_start = ittt - ilen*2;
            
            if (fdebug>2){
                std::cout << "      Board: " << fragid
                          << ", ttt: " << ittt
                          << ", frag ts: " << ifrag.timestamp()%uint(1e9)
                          << ", length: " << ilen
                          << std::endl;
            }
            if (extended_flag){
                for (size_t jtrig=itrig+1; jtrig < frag_v.size(); jtrig++){
                    auto jfrag = frag_v.at(jtrig);
                    auto jlen = get_length(jfrag);
                    auto jttt = get_ttt(jfrag);
                    // to be an extended trigger, need to fulfill both conditions:
                    // 1. the length of the extension will be less than the nominal length 
                    // 2. the time between the first ttt and the extended ttt will be less than the nominal length*2 (2 ns tick)
                    bool length_check = (jlen < fnominal_length);
                    bool ttt_check=false;
                    if (iwvfm_end < uint(1e9-fnominal_length*2)){
                        if ((jttt-iwvfm_end) < fnominal_length*2)
                            ttt_check = true;
                    }
                    // if the end of wvfm is less than 10 us away from end of the second, may be second rollover
                    else if ((iwvfm_end > uint(1e9-fnominal_length*2)) && (jttt < fnominal_length*2)) {
                        if ((jttt + 1e9 - iwvfm_end) < fnominal_length*2){
                            ttt_check = true;
                        }
                    }
                    if (length_check && ttt_check){
                        if (fdebug>3){
                        std::cout << "      Board: " << fragid
                          << ", ttt: " << jttt
                          << ", frag ts: " << jfrag.timestamp()%uint(1e9)
                          << ", length: " << jlen
                          << std::endl;
                        }
                        extensions_used++;
                        board_ttt_v.push_back(jttt);
                        std::vector<std::vector<uint16_t>> jwvfm_v;
                        get_waveforms(jfrag, jwvfm_v);
                        for (size_t ich=0; ich < iwvfm_v.size(); ich++){
                            std::vector<uint16_t>& orig_wvfm = iwvfm_v.at(ich);
                            std::vector<uint16_t>  extd_wvfm = jwvfm_v.at(ich); // extended waveform 
                            orig_wvfm.insert( orig_wvfm.end(), extd_wvfm.begin(), extd_wvfm.end());
                        } // end ch loop
                    
                    } // end if extension condition
                    else break; // if the immediate next trigger is not an extension, break
                    iwvfm_end = jttt;
                } // end jtrig loop
                if (fdebug>2){
                    std::cout << "      Board: " << fragid
                              << " -> start time: " << int(iwvfm_start) - int(event_trigger_time)
                              << " , extended wvfm length: " << iwvfm_v.at(0).size()
                              << std::endl;
                }
            } // end extended flag
            board_info.triggerTimeTag = board_ttt_v;
            brdTimingInfoVec->push_back(board_info);
            art::Ptr<raw::pmt::BoardTimingInfo> brdTimingInfoPtr = make_pmtbrd_ptr(brdTimingInfoVec->size()-1);

            for (size_t i = 0; i < iwvfm_v.size(); i++){
                auto combined_wvfm = iwvfm_v[i];
                if (evt_counter==fhist_evt){                    
                    // histo: save waveforms section for combined waveforms
                    if ((fragid==8) && (std::find(fignore_timing_ch.begin(), fignore_timing_ch.end(), i) != fignore_timing_ch.end()))
                        continue;
                    histname.str(std::string());
                    histname << "evt" << evt.event() << "_trig" << trig_counter << "_frag" << fragid << "_ch" << i << "_combined_wvfm";

                    TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), combined_wvfm.size(), 0, combined_wvfm.size());
                    wvfmHist->GetXaxis()->SetTitle("ticks");
                    for(unsigned int n = 0; n < combined_wvfm.size(); n++) 
                        wvfmHist->SetBinContent(n + 1, (double)combined_wvfm[n]);
                }
                double time_diff = (int(iwvfm_start) - int(event_trigger_time))*1e-3; // us
                if ((std::abs(time_diff) > fallowed_time_diff) && (event_trigger_time!=0)){
                    // second rollover between reference time and waveform start 
                    if (std::abs(time_diff + 1e6) < fallowed_time_diff)
                        time_diff += 1e6; // us 
                    // second rollover between waveform start and reference time
                    else if (std::abs(time_diff - 1e6) < fallowed_time_diff)
                        time_diff -= 1e6; // us
                    else if (fdebug>1)
                        std::cout << "WARNING: TIME DIFFERENCE IS GREATER THAN " << fallowed_time_diff << " us. Event timestamp: " << event_trigger_time << ", waveform timestamp: " << iwvfm_start << std::endl;
                }
                if ((i == 15) && (foutput_ftrig_wvfm)){
                    // for the timing ch of each board, the opdetwaveform channel is equal to the board id
                    raw::OpDetWaveform waveform(time_diff, fragid, combined_wvfm);
                    fltwvfmVec->push_back(waveform);
                    art::Ptr<raw::OpDetWaveform> wvfmPtr = make_fltwvfm_ptr(fltwvfmVec->size()-1);
                    fltTimingAssns->addSingle(brdTimingInfoPtr, wvfmPtr);
                }
                else if ((fragid == 8)){ // fyi: this hardcodes the timing caen board index
                    // for the timing caen, the opdetwaveform channel is offset (900) + chidx 
                    raw::OpDetWaveform waveform(time_diff, ftiming_caen_offset + i, combined_wvfm);
                    if ((foutput_timing_wvfm) && (std::find(fignore_timing_ch.begin(), fignore_timing_ch.end(), i) == fignore_timing_ch.end())){
                        timwvfmVec->push_back(waveform);
                        art::Ptr<raw::OpDetWaveform> wvfmPtr = make_timwvfm_ptr(timwvfmVec->size()-1);
                        timTimingAssns->addSingle(brdTimingInfoPtr, wvfmPtr);
                    }
                }
                else{
                    // for normal pmt chs, the opdetwaveform channel is from the ch_map
                    raw::OpDetWaveform waveform(time_diff, ch_map.at(fragid*15+i), combined_wvfm);
                    pmtwvfmVec->push_back(waveform);
                    art::Ptr<raw::OpDetWaveform> wvfmPtr = make_pmtwvfm_ptr(pmtwvfmVec->size()-1);
                    pmtTimingAssns->addSingle(brdTimingInfoPtr, wvfmPtr);
                }
            }
            trig_counter++;       
        } // end itrig loop

        if (fdebug>1){
            bool trigger_bookkeeping = (full_found + extensions_used) == int(frag_v.size());
            bool found_used = (extensions_used == extensions_found);
            std::cout << "      Board: " << fragid
                      << ", passed trigger bookkeeping: " << (trigger_bookkeeping && found_used)
                      << std::endl;
            if (fdebug>2){
            // print number of triggers, number of extensions found, and number of extensions used
                std::cout << "      Board: " << fragid
                          << " -> # of triggers: " << frag_v.size()
                          << ", full found: " << full_found
                          << ", extensions found: " << extensions_found
                          << ", extensions used: " << extensions_used
                          << std::endl;
            }
        }  
    } // end board loop
    board_frag_v.clear();

    // loop through flashes
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<calib::TriggerEmulationService> fTriggerService;
    int PMTPerBoard = fTriggerService->getPMTPerBoard();
    std::vector<int> MonPulsesFlat;
    std::vector<int> pulseSizes;
    MonPulsesFlat.clear();
    pulseSizes.clear();
    int TotalFlash = pmtwvfmVec->size()/((int)fn_caenboards*PMTPerBoard); // pmtwvfmVec = waveHandle ???
    for (int FlashCounter=0; FlashCounter<TotalFlash; FlashCounter++)
    {
      int WaveIndex = FlashCounter*PMTPerBoard;
      int WaveformSize = (*pmtwvfmVec)[WaveIndex].size();
      std::vector<int> *MonPulse = new std::vector<int>(WaveformSize);
      fTriggerService->ConstructMonPulse(*pmtwvfmVec, fmon_threshold, MonPulse, FlashCounter);
      //MonPulsesAll.push_back(std::move(MonPulse));
      MonPulsesFlat.insert(MonPulsesFlat.end(), (*MonPulse).begin(), (*MonPulse).end());
      pulseSizes.push_back(MonPulse->size());
      delete MonPulse;
    }
    // make ptrs
    auto flatPtr = std::make_unique<std::vector<int>>(std::move(MonPulsesFlat));
    auto sizesPtr = std::make_unique<std::vector<int>>(std::move(pulseSizes));

    evt.put(std::move(pmtwvfmVec),fpmt_instance_name);  
    evt.put(std::move(fltwvfmVec),fflt_instance_name);
    evt.put(std::move(timwvfmVec),ftim_instance_name);

    evt.put(std::move(evtTimingInfo));
    evt.put(std::move(brdTimingInfoVec));

    evt.put(std::move(pmtTimingAssns),fpmt_timing_instance_name);
    evt.put(std::move(fltTimingAssns),fflt_timing_instance_name);
    evt.put(std::move(timTimingAssns),ftim_timing_instance_name);

    evt.put(std::move(flatPtr), "MonPulses");
    evt.put(std::move(sizesPtr), "MonPulseSizes");
}

void sbndaq::SBNDPMTDecoder::get_fragments(artdaq::Fragment & frag, std::vector<std::vector<artdaq::Fragment>> & board_frag_v){
    auto fragid = frag.fragmentID() - ffragid_offset;
    if (fuse_set_map){
        auto it = std::find(fset_fragid_map.begin(), fset_fragid_map.end(), fragid);
        if (it != fset_fragid_map.end()){
            auto idx = std::distance(fset_fragid_map.begin(), it);
            fragid = idx;
        }
    }
    // ignore boards that are in the list of boards to ignore
    if (std::find(fignore_fragid.begin(), fignore_fragid.end(), fragid) != fignore_fragid.end())
        return;
    else if ((fragid>=0) && (fragid < fn_caenboards )){
        board_frag_v.at(fragid).push_back(frag);
    }
    // if fragid is out of range
    else{
        std::cout << "Fragment ID " << fragid << " is out of range. FragID offset/FragID map may be misconfigured, or this FragID is not attributed to a PMT Digitizer. Skipping this fragment..." << std::endl;
    }
}

void sbndaq::SBNDPMTDecoder::get_waveforms(artdaq::Fragment & frag, std::vector<std::vector<uint16_t>> & wvfm_v){
    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    auto nch = md->nChannels;    

    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t ev_size_quad_bytes = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wvfm_length = data_size_double_bytes/nch;

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
                                   + sizeof(CAENV1730EventHeader));
    const uint16_t* value_ptr =  data_begin;
    uint16_t value = 0;

    for (size_t i_ch=0; i_ch<nch; ++i_ch){
        auto ch_offset = (size_t)(i_ch * wvfm_length);      
        std::vector<uint16_t> wvfm(wvfm_length,0);

        //--loop over waveform entries
        for(size_t i_t=0; i_t<wvfm_length; ++i_t){ 
            value_ptr = data_begin + ch_offset + i_t; 
            value = *(value_ptr);
            wvfm.at(i_t) = value;
        }
        wvfm_v.push_back(wvfm);
    }
}

void sbndaq::SBNDPMTDecoder::get_timing(artdaq::Fragment & frag,
                                       uint16_t & frag_postpercent,
                                       uint32_t & frag_ttt, 
                                       uint32_t & frag_len, 
                                       int & frag_tick){
    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    auto nch = md->nChannels;
    
    if (nch < 16){
        std::cout << "No Ch15 in the fragment, FTRIG channel missing? FTRIG value set to 0." << std::endl;
        frag_tick = 0;
        return;
    }

    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t ev_size_quad_bytes = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wvfm_length = data_size_double_bytes/nch;

    auto frag_ts = frag.timestamp()%uint(1e9);
    frag_ttt = header.triggerTimeTag*8;
    frag_len = wvfm_length;

    if (wvfm_length < fnominal_length){ // if an extension
        frag_postpercent = fdefault_postpercent;
    }
    else{
        if (frag_ts!=frag_ttt){
            // this assumes that if timestamp != ttt... then the time tag shift is enabled in the caen configuration
            if (frag_ttt>frag_ts)
                frag_postpercent = round((100.0*(frag_ttt-frag_ts)/(frag_len*2.0)));
            else // second rollover between ts and ttt
                frag_postpercent = round((100.0*(1e9+frag_ttt-frag_ts)/(frag_len*2.0))); 
            if (frag_postpercent > 100){
                std::cout << "WARNING: Postpercent is not correctly configured!! Saving default postpercent." << std::endl;
                frag_postpercent = fdefault_postpercent;
            }
        }
        else
            frag_postpercent = fdefault_postpercent;
    }

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
								   + sizeof(CAENV1730EventHeader));
    const uint16_t* value_ptr =  data_begin;
    uint16_t value = 0;

    auto ch15_offset = (size_t)(15 * wvfm_length);

    for(size_t i_t=0; i_t<wvfm_length; ++i_t){ 
        value_ptr = data_begin + ch15_offset + i_t; 
        value = *(value_ptr);
        if (value > fthreshold_ftrig){
            frag_tick = i_t;
            break;
        }
    }
}

uint32_t sbndaq::SBNDPMTDecoder::get_length(artdaq::Fragment & frag){
    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    auto nch = md->nChannels;    

    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t ev_size_quad_bytes = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wvfm_length = data_size_double_bytes/nch;

    return wvfm_length;
}

uint32_t sbndaq::SBNDPMTDecoder::get_ttt(artdaq::Fragment & frag){
    CAENV1730Fragment bb(frag);
    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t ttt = header.triggerTimeTag*8;
    return ttt;
}

uint32_t sbndaq::SBNDPMTDecoder::get_boardid(artdaq::Fragment & frag){
    CAENV1730Fragment bb(frag);
    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t boardid = header.boardID;
    return boardid;
}

std::vector<uint> sbndaq::SBNDPMTDecoder::fill_chmap(sbndDB::PMTCalibrationDatabase const* pmt_calib_db){
    std::vector<uint> ch_map(fn_caenboards*fn_caenchannels,9999);
    if (pmt_calib_db==nullptr){
        throw std::runtime_error("PMT Calibration Database pointer is null.");
    }
    auto nopdets = opdetmap.size();
    for (size_t idet=0; idet<nopdets; idet++){
        std::string pd_type = opdetmap.pdType(idet);
        if (pd_type.find("pmt") == std::string::npos) continue;
        // the boards in the database are 1-indexed, need to convert to 0-indexed
        int board = pmt_calib_db->getCAENDigitizer(idet) - 1;
        int channel = pmt_calib_db->getCAENDigitizerChannel(idet);
        ch_map[board*15 + channel] = idet;
    }
    return ch_map;
}

DEFINE_ART_MODULE(sbndaq::SBNDPMTDecoder)
