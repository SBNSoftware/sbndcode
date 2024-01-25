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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndcode/Decoders/PTB/sbndptb.h"

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
    bool check_fragments(std::vector<artdaq::Fragment> frag_v, uint32_t& trig_ttt, uint32_t& trig_len);
    void get_waveforms(artdaq::Fragment & frag, std::vector<std::vector<uint16_t>> & wvfm_v);

    uint32_t get_length(artdaq::Fragment & frag);
    uint32_t get_ttt(artdaq::Fragment & frag);
    uint32_t get_boardid(artdaq::Fragment & frag);

    void get_timing(artdaq::Fragment & frag, uint32_t & ttt, uint32_t & len, int & tick);
    
    std::vector<raw::OpDetWaveform> fWaveforms;
    std::vector<raw::OpDetWaveform> fTriggerWaveforms;

    std::vector<std::vector<artdaq::Fragment>>  trig_frag_v; // every entry should correspond to 1 [flash] trigger 
    uint64_t event_trigger_time = 0; 

private:
    bool fverbose;

    std::vector<std::string>  fcaen_fragment_name;
    std::string               fcaen_module_label;

    std::vector<uint32_t>     fignore_boards;
    uint32_t                  fnominal_length; 

    std::string               fspectdc_product_name;
    uint32_t                  fspectdc_ftrig_ch;
    uint32_t                  fspectdc_etrig_ch;

    std::string               fptb_product_name;

    std::string fch_instance_name;
    std::string ftr_instance_name;

    int fn_maxflashes;
    int fn_caenboards;
    uint16_t fthreshold_ftrig;
    uint fhist_evt;

    // histogram info  
    std::stringstream histname; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;
    uint evt_counter = 0;
};


sbndaq::SBNDPMTDecoder::SBNDPMTDecoder(fhicl::ParameterSet const& p)
    : EDProducer{p}  // ,
{
    fverbose          = p.get<bool>("Verbose",true);

    fcaen_fragment_name = p.get<std::vector<std::string>>("caen_fragment_name");
    fcaen_module_label  = p.get<std::string>("caen_module_label","daq");

    fignore_boards = p.get<std::vector<uint32_t>>("ignore_boards",{});
    fnominal_length = p.get<uint32_t>("nominal_length",5000);

    fspectdc_product_name = p.get<std::string>("spectdc_product_name","tdcdecoder");
    fspectdc_ftrig_ch = p.get<uint32_t>("spectdc_ftrig_ch",3);
    fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch",4);

    fptb_product_name = p.get<std::string>("ptb_product_name","ptbdecoder");

    fch_instance_name = p.get<std::string>("pmtInstanceName","PMTChannels");
    ftr_instance_name = p.get<std::string>("ftrigInstanceName","FTrigChannels");

    fn_maxflashes = p.get<int>("n_maxflashes",30);
    fn_caenboards = p.get<int>("n_caenboards",8);
    fthreshold_ftrig = p.get<uint16_t>("threshold_ftrig",16350);
    fhist_evt     = p.get<int>("hist_evt",1);

    produces< std::vector< raw::OpDetWaveform > >(fch_instance_name); 
    produces< std::vector< raw::OpDetWaveform > >(ftr_instance_name);
}

void sbndaq::SBNDPMTDecoder::produce(art::Event& evt)
{
    std::unique_ptr< std::vector< raw::OpDetWaveform > > wvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    std::unique_ptr< std::vector< raw::OpDetWaveform > > twvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());

    std::cout << "event: " << evt.id().event() << std::endl;
    evt_counter++;

    trig_frag_v.clear(); 
    trig_frag_v.resize(fn_maxflashes);
    for (auto &v : trig_frag_v)
        v.reserve(fn_caenboards);

    fWaveforms.clear();
    fTriggerWaveforms.clear();
    uint ncont = 0; // counter for number of containers

    bool found_caen = false;
    // get CAEN fragments 
    for (const std::string &caen_name : fcaen_fragment_name){
        std::cout << "Looking for CAEN V1730 fragments with label " << caen_name << "..." << std::endl;
        art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
        evt.getByLabel(fcaen_module_label, caen_name, fragmentHandle);

        if (!fragmentHandle.isValid() || fragmentHandle->size() == 0){
            std::cout << "No CAEN V1730 fragments found." << std::endl;
            continue;
        }
        else
            found_caen = true;

        if (fragmentHandle->front().type() == artdaq::Fragment::ContainerFragmentType){
            for (auto cont : *fragmentHandle){
                artdaq::ContainerFragment contf(cont);
                ncont++;

                if (contf.fragment_type() == sbndaq::detail::FragmentType::CAENV1730){
                    if (ncont==1){ 
                        // use the size of the first container to re-initialize 
                        trig_frag_v.resize(contf.block_count());
                        std::cout << "Found " << contf.block_count() << " CAEN fragments (flash triggers) inside the container." << std::endl;
                    }
                    for (size_t ii = 0; ii < contf.block_count(); ++ii){
                        auto frag = *contf[ii].get();                        
                        CAENV1730Fragment bb(frag);
                        CAENV1730Event const* event_ptr = bb.Event();
                        CAENV1730EventHeader header = event_ptr->Header;
                        auto boardID = header.boardID;                        
                    
                        // ignore boards that are not in the list of boards to ignore
                        if (std::find(fignore_boards.begin(), fignore_boards.end(), boardID) != fignore_boards.end())
                            continue;
                        if (int(ii) >= fn_maxflashes){
                            std::cout << "Warning: more than " << fn_maxflashes << " flash triggers found, update the fcl! Skipping the rest." << std::endl;
                            break;
                        }
                        trig_frag_v[ii].push_back(frag);

                    }
                }
            }
        }
        else if (fragmentHandle->front().type() == sbndaq::detail::FragmentType::CAENV1730){
        	std::cout << "Found " << fragmentHandle->size() << " normal CAEN fragments " << std::endl;
        }
    }

    if (found_caen==false){
        std::cout << "No Trigger information found, pushing empty waveforms." << std::endl;
        evt.put(std::move(wvfmVec),fch_instance_name);  
        evt.put(std::move(twvfmVec),ftr_instance_name);
        return;
    }
    
    // get spec tdc product 
    art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
    evt.getByLabel(fspectdc_product_name,tdcHandle);
    if (!tdcHandle.isValid() || tdcHandle->size() == 0){
        std::cout << "No SPECTDC products found." << std::endl;
    }
    else{
        std::cout << "SPECTDC (decoded) products found: " << std::endl;
        const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);

        for (size_t i=0; i<tdc_v.size(); i++){
            auto tdc = tdc_v[i];
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();
            const uint64_t  offset = tdc.Offset();
            const std::string name  = tdc.Name();
        
            if (ch==fspectdc_etrig_ch || ch==fspectdc_ftrig_ch)
            std::cout << "      TDC CH " << ch << " ->"
            << ", name: " << name
            << ", ts (ns): " << ts%uint64_t(1e9) 
            << ", offset: " << offset 
            << ", ts (s+ns): " << float(ts)*1e-9
            << std::endl;

            if (ch==fspectdc_etrig_ch)
                event_trigger_time = ts%uint64_t(1e9);
        }
    }

    // get ptb product
    art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
    evt.getByLabel(fptb_product_name,ptbHandle);
    if (!ptbHandle.isValid() || ptbHandle->size() == 0){
        std::cout << "No PTB products found." << std::endl;
    }
    else{
        std::cout << "PTB (decoded) products found: " << std::endl;
        const std::vector<raw::ptb::sbndptb> ptb_v(*ptbHandle);
        for (size_t i=0; i<ptb_v.size(); i++){
            auto ptb = ptb_v[i];
            auto hltrigs = ptb.GetHLTriggers();
            auto lltrigs = ptb.GetLLTriggers();

            for (size_t j=0; j < hltrigs.size(); j++){
                raw::ptb::Trigger trig = hltrigs.at(j);
                std::cout << "      PTB HLT " <<  j << "-> " 
                          << "ts (ns): " << (trig.timestamp * 20)%(uint(1e9))
                          << ", trigger word: " << std::hex << trig.trigger_word << std::dec
                          << ", word type: " << trig.word_type
                          << std::endl;
            }
            // for (size_t j=0; j < lltrigs.size(); j++){
            //     raw::ptb::Trigger trig = lltrigs.at(j);
            //     std::cout << "PTB LLT " << j << " -> " 
            //               << "ts (ns): " << (trig.timestamp * 20)%(uint(1e9))
            //               << ", trigger word: " << std::hex << trig.trigger_word << std::dec
            //               << ", word type: " << trig.word_type
            //               << std::endl;
            // }
        }
    }


    auto ntrig = trig_frag_v.size();

    std::vector<uint32_t> trig_ttt_v(ntrig,0);
    std::vector<uint32_t> trig_len_v(ntrig,0);

    for (size_t i=0; i < trig_frag_v.size(); i++){
        // frag_v contains all fragments for a single trigger 
        auto frag_v = trig_frag_v.at(i);
        std::cout << "CAEN Trigger " << i << " has " << frag_v.size() << " fragments" << std::endl;
        if (frag_v.empty()) continue;
        uint32_t trig_ttt = 0;
        uint32_t trig_len = 0;

        bool pass_check = check_fragments(frag_v, trig_ttt, trig_len);        
        if (pass_check == false)
            continue;
        // TODO: once we get the event trigger from other timing subystems,
        //       we can pass trigger based on the event trigger time
        //       instead of passing all triggers
        trig_ttt_v.at(i) = trig_ttt;
        trig_len_v.at(i) = trig_len;
    }

    // check to see if there are any extended triggers in this event 
    bool extended_flag = false;
    for (size_t itrig=0; itrig < ntrig; itrig++){
        auto ilen = trig_len_v.at(itrig);
        if (ilen < fnominal_length){
            extended_flag = true;
            break;
        }
    }

    for (size_t itrig=0; itrig < ntrig; itrig++){

        auto ittt = trig_ttt_v.at(itrig);
        auto ilen = trig_len_v.at(itrig);
        auto ifrag_v = trig_frag_v.at(itrig);

        std::vector<int> fragid_v(ifrag_v.size(),0);
        std::vector<uint32_t> id_board_v(ifrag_v.size(),0);

        // if this waveform is short, skip it 
        
        if (ilen < fnominal_length) continue;

        // store the waveform itself
        // store the waveform start time (for OpDetWaveform timestamp)
        // store the time at the end of the waveform (for extended trigger condition)
        std::vector<std::vector<uint16_t>> iwvfm_v;
        std::vector<int> iwvfm_start_v(ifrag_v.size(),0);
        auto iwvfm_end = ittt;

        for (size_t idx=0; idx<ifrag_v.size(); idx++){
            auto frag = ifrag_v.at(idx);
            get_waveforms(frag, iwvfm_v);
            fragid_v.at(idx) = frag.fragmentID();
            id_board_v.at(idx) = get_boardid(frag);
            iwvfm_start_v.at(idx) = get_ttt(frag) - 2*(get_length(frag));

            uint32_t frag_ttt = 0; 
            uint32_t frag_len = 0;
            int      frag_tick = 0;
            
            get_timing(frag, frag_ttt, frag_len, frag_tick);
            int frag_ts = frag_ttt - 2*(frag_len - frag_tick); // ns
            auto length = get_length(frag);
            std::cout << "      Board ID: " << get_boardid(frag)
                      << ", Frag ID: " << frag.fragmentID()
                      << ", ttt: " << frag_ttt
                      << ", tick: " << frag_tick
                      << ", trig ts: " << frag_ts
                      << ", start ts: " << frag_ttt - 2*(frag_len)
                      << ", length: " << length
                      << std::endl;
        }
        // if there exist fragments that may be part of the same trigger 
        if (extended_flag){
            for (size_t jtrig=itrig+1; jtrig < ntrig; jtrig++){
                bool pass_checks = true;
                auto jttt = trig_ttt_v.at(jtrig);
                auto jlen = trig_len_v.at(jtrig);
                auto jfrag_v = trig_frag_v.at(jtrig);

                // if the next trigger is more than 10 us away or is equal to the nominal length, stop looking
                if (((signed)(jttt - iwvfm_end) > 1e4) || (jlen >= fnominal_length)) break; 
                else if (((jttt - iwvfm_end) < 1e4) && (jlen < fnominal_length)){
                    std::vector<std::vector<uint16_t>> jwvfm_v;
                    for (size_t idx=0; idx<jfrag_v.size(); idx++){
                        auto frag = jfrag_v.at(idx);
                        if (fragid_v.at(idx) != frag.fragmentID()){
                            // check that the two fragments originated from the same board
                            std::cout << "Error: fragment IDs do not match between triggers " << itrig << " and " << jtrig << std::endl;
                            pass_checks=false;
                            break;
                        }
                        get_waveforms(frag, jwvfm_v);
                    }
                    // check that the number of waveforms are the same 
                    if (jwvfm_v.size() != iwvfm_v.size()){
                        std::cout << "Error: number of channels in extended fragment " << jtrig << " does not match number of channels in original fragment " << itrig << std::endl;
                        pass_checks=false;
                    }
                    if (pass_checks==false) continue;
                    // combine the waveforms
                    for (size_t ich=0; ich < iwvfm_v.size(); ich++){
                        std::vector<uint16_t>& orig_wvfm = iwvfm_v.at(ich);
                        std::vector<uint16_t>  extd_wvfm = jwvfm_v.at(ich); // extended waveform 
                        orig_wvfm.insert( orig_wvfm.end(), extd_wvfm.begin(), extd_wvfm.end());
                    }
                    // update the end time of the waveform
                    iwvfm_end = jttt;
                } // extended trigger found 
            } // loop over subsequent triggers
        }


        std::cout << "Obtained waveforms for TTT with TTT " << ittt << "\n" 
                    << "\tNumber of channels: " << iwvfm_v.size() << "\n"
                    << "\tNumber of entries: " << iwvfm_v.at(0).size() << "\n"
                    << "\tTimestamp (ns): " << int(event_trigger_time) - int(iwvfm_start_v.at(0)) << std::endl;

        for (size_t i = 0; i < iwvfm_v.size(); i++){
            auto combined_wvfm = iwvfm_v[i];
            int board_idx = int(i)/16; // integer division to get the board index

            if (evt_counter==fhist_evt){
            // histo: save waveforms section for combined waveforms
                histname.str(std::string());
                histname << "evt" << evt.event() << "_board" << id_board_v.at(board_idx) << "_ch" << i%16 << "_combined_wvfm";

                TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), combined_wvfm.size(), 0, combined_wvfm.size());
                wvfmHist->GetXaxis()->SetTitle("ticks");
                for(unsigned int n = 0; n < combined_wvfm.size(); n++) 
                    wvfmHist->SetBinContent(n + 1, (double)combined_wvfm[n]);
            }
            int time_diff = int(event_trigger_time) - int(iwvfm_start_v.at(board_idx));
            raw::OpDetWaveform waveform(time_diff, i%16, combined_wvfm);

            if (i%16 == 15)
                twvfmVec->push_back(waveform);
            else
                wvfmVec->push_back(waveform);
        }

    } // loop over triggers 

    
    evt.put(std::move(wvfmVec),fch_instance_name);  
    evt.put(std::move(twvfmVec),ftr_instance_name);
}

bool sbndaq::SBNDPMTDecoder::check_fragments(std::vector<artdaq::Fragment> frag_v, uint32_t& trig_ttt, uint32_t& trig_len){
    bool pass = true;
    uint32_t this_ttt = 0; 
    uint32_t this_len = 0;
    for (size_t ifrag=0; ifrag<frag_v.size(); ifrag++){
        CAENV1730Fragment bb(frag_v.at(ifrag));
        CAENV1730Event const* event_ptr = bb.Event();
        auto const * md = bb.Metadata();
        CAENV1730EventHeader header = event_ptr->Header;
        // uint32_t boardid = header.boardID;

        uint32_t ttt = header.triggerTimeTag*8;
        auto nch = md->nChannels;    
        uint32_t ev_size_quad_bytes = header.eventSize;
        uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
        uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
        uint32_t wvfm_length = data_size_double_bytes/nch; // downsampled trigger time tag; TTT points to end of wvfm w.r.t. pps

        if (ifrag==0){ 
            this_ttt = ttt; 
            this_len = wvfm_length; 
        }
        // if (this_ttt != ttt || this_len != wvfm_length){
        //     if (abs((signed)(this_ttt - ttt)) > 16 || abs((signed)(this_len - wvfm_length)) > 8 ){
        //         std::cout << "Mismatch is greater than maximum 16 ns jitter" << std::endl;
        //         pass=false;
        //         break;
        //     }
        // }
    }
    trig_ttt = this_ttt;
    trig_len = this_len;
    return pass;
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
                                       uint32_t & frag_ttt, 
                                       uint32_t & frag_len, 
                                       int & frag_tick){
    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    auto nch = md->nChannels;    

    CAENV1730Event const* event_ptr = bb.Event();
    CAENV1730EventHeader header = event_ptr->Header;
    uint32_t ev_size_quad_bytes = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wvfm_length = data_size_double_bytes/nch;

    frag_ttt = header.triggerTimeTag*8;
    frag_len = wvfm_length;

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

DEFINE_ART_MODULE(sbndaq::SBNDPMTDecoder)
