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
    bool check_fragments(std::vector<artdaq::Fragment> &frag_v, uint32_t& trig_ttt, uint32_t& trig_len);
    void get_waveforms(artdaq::Fragment & frag, std::vector<std::vector<uint16_t>> & wvfm_v);

    uint32_t get_length(artdaq::Fragment & frag);
    uint32_t get_ttt(artdaq::Fragment & frag);
    uint32_t get_boardid(artdaq::Fragment & frag);

    void get_timing(artdaq::Fragment & frag, uint32_t & ttt, uint32_t & len, int & tick);
    

private:
    uint fdebug;
    uint                      ftiming_type;

    std::vector<std::string>  fcaen_fragment_name;
    std::string               fcaen_module_label;

    std::vector<uint32_t>     fignore_fragid;
    uint32_t                  fnominal_length; 

    std::string               fspectdc_product_name;
    uint32_t                  fspectdc_ftrig_ch;
    uint32_t                  fspectdc_etrig_ch;

    std::string               fptb_product_name;
    ULong64_t                 fptb_etrig_trigword; 
    uint32_t                  fptb_etrig_wordtype;

    std::string fch_instance_name;
    std::string ftr_instance_name;
    bool foutput_ftrig_wvfm;

    int fn_maxflashes;
    int fn_caenboards;
    uint16_t fthreshold_ftrig;

    uint ffragid_offset;
    uint fhist_evt;
    std::vector<uint> fch_map;

    uint32_t fnch; // number of channels
    // histogram info  
    std::stringstream histname; //raw waveform hist name
    art::ServiceHandle<art::TFileService> tfs;
    uint evt_counter = 0;
};


sbndaq::SBNDPMTDecoder::SBNDPMTDecoder(fhicl::ParameterSet const& p)
    : EDProducer{p}  // ,
{
    fdebug          = p.get<uint>("debug",0);

    fcaen_fragment_name = p.get<std::vector<std::string>>("caen_fragment_name");
    fcaen_module_label  = p.get<std::string>("caen_module_label","daq");

    fignore_fragid = p.get<std::vector<uint32_t>>("ignore_fragid",{});
    fnominal_length = p.get<uint32_t>("nominal_length",5000);

    ftiming_type = p.get<uint>("timing_type",0);

    fspectdc_product_name = p.get<std::string>("spectdc_product_name","tdcdecoder");
    fspectdc_ftrig_ch = p.get<uint32_t>("spectdc_ftrig_ch",3);
    fspectdc_etrig_ch = p.get<uint32_t>("spectdc_etrig_ch",4);

    fptb_product_name = p.get<std::string>("ptb_product_name","ptbdecoder");
    fptb_etrig_trigword = p.get<ULong64_t>("ptb_etrig_trigword",0x0000000000000000);
    fptb_etrig_wordtype = p.get<uint32_t>("ptb_etrig_wordtype",2);

    fch_instance_name = p.get<std::string>("pmtInstanceName","PMTChannels");
    ftr_instance_name = p.get<std::string>("ftrigInstanceName","FTrigChannels");
    foutput_ftrig_wvfm = p.get<bool>("output_ftrig_wvfm",true);

    fn_maxflashes    = p.get<int>("n_maxflashes",30);
    fn_caenboards    = p.get<int>("n_caenboards",8);
    fthreshold_ftrig = p.get<uint16_t>("threshold_ftrig",16350);
    ffragid_offset   = p.get<uint>("fragid_offset",40960);
    fhist_evt        = p.get<int>("hist_evt",1);
    fch_map          = p.get<std::vector<uint>>("ch_map",{});

    produces< std::vector< raw::OpDetWaveform > >(fch_instance_name); 
    produces< std::vector< raw::OpDetWaveform > >(ftr_instance_name);
}

void sbndaq::SBNDPMTDecoder::produce(art::Event& evt)
{
    std::unique_ptr< std::vector< raw::OpDetWaveform > > wvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    std::unique_ptr< std::vector< raw::OpDetWaveform > > twvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());

    evt_counter++;

    std::vector<std::vector<artdaq::Fragment>>  trig_frag_v; // every entry should correspond to 1 [flash] trigger 
    trig_frag_v.resize(fn_maxflashes);
    for (auto &v : trig_frag_v)
        v.reserve(fn_caenboards);

    uint64_t event_trigger_time = 0; 

    uint ncont = 0; // counter for number of containers

    bool found_caen = false;
    // get CAEN fragments 
    for (const std::string &caen_name : fcaen_fragment_name){
        art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
        if (fdebug>0) std::cout << "Looking for CAEN V1730 fragments with label " << caen_name << "..." << std::endl;
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
                    if (ncont==1){ 
                        // use the size of the first container to re-initialize 
                        trig_frag_v.resize(contf.block_count());
                        if (fdebug>1) std::cout << "Found " << contf.block_count() << " CAEN fragments (flash triggers) inside the container." << std::endl;
                    }
                    else{
                        if (trig_frag_v.size() != contf.block_count())
                            std::cout << "WARNING! number of CAEN fragments in container " << ncont << " does not match the first container. Missing board fragments?" << std::endl;
                    }
                    for (size_t ii = 0; ii < contf.block_count(); ++ii){
                        auto frag = *contf[ii].get();                        
                        auto fragid = frag.fragmentID() - ffragid_offset;
                        // ignore boards that are not in the list of boards to ignore
                        if (std::find(fignore_fragid.begin(), fignore_fragid.end(), fragid) != fignore_fragid.end())
                            continue;
                        if (int(ii) >= fn_maxflashes){
                            if (fdebug>0) std::cout << "Warning: more than " << fn_maxflashes << " flash triggers found, update the fcl! Skipping the rest." << std::endl;
                            break;
                        }
                        trig_frag_v[ii].push_back(frag);

                    }
                }
            }
        }
        else if (fragmentHandle->front().type() == sbndaq::detail::FragmentType::CAENV1730){
        	if (fdebug>1) std::cout << "Found " << fragmentHandle->size() << " normal CAEN fragments " << std::endl;
        }
        fragmentHandle.removeProduct();
    }

    if (found_caen==false){
        if (fdebug>0) std::cout << "No CAEN V1730 fragments of any type found, pushing empty waveforms." << std::endl;

        evt.put(std::move(wvfmVec),fch_instance_name);  
        evt.put(std::move(twvfmVec),ftr_instance_name);
        return;
    }
    
    // create a timing type per event so the default doesn't get overwritten
    auto timing_type = ftiming_type; 

    // get spec tdc product
    if (timing_type==0){
        art::Handle<std::vector<sbnd::timing::DAQTimestamp>> tdcHandle;
        evt.getByLabel(fspectdc_product_name,tdcHandle);
        if (!tdcHandle.isValid() || tdcHandle->size() == 0){
            if (fdebug>0) std::cout << "No SPECTDC products found." << std::endl;
            timing_type++;
        }
        else{
            if (fdebug>1) std::cout << "SPECTDC (decoded) products found: " << std::endl;
            const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*tdcHandle);

            for (size_t i=0; i<tdc_v.size(); i++){
                auto tdc = tdc_v[i];
                const uint32_t  ch = tdc.Channel();
                const uint64_t  ts = tdc.Timestamp();
                const uint64_t  offset = tdc.Offset();
                const std::string name  = tdc.Name();
            
                // if (ch==fspectdc_etrig_ch || ch==fspectdc_ftrig_ch)
                if (fdebug>1){
                    std::cout << "      TDC CH " << ch << " -> "
                    << "name: " << name
                    << ", ts (ns): " << ts%uint64_t(1e9)
                    << ", sec (s): " << ts/uint64_t(1e9)
                    << ", offset: " << offset 
                    << std::endl;
                }
                if (ch==fspectdc_etrig_ch)
                    event_trigger_time = ts%uint64_t(1e9);
            }
        }
    }
    if (timing_type==1){
        // get ptb product
        art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
        evt.getByLabel(fptb_product_name,ptbHandle);
        if ((!ptbHandle.isValid() || ptbHandle->size() == 0)){
            if (fdebug>0) std::cout << "No PTB products found." << std::endl;
            timing_type++;
        }
        else{
            const std::vector<raw::ptb::sbndptb> ptb_v(*ptbHandle);
            for (size_t i=0; i<ptb_v.size(); i++){
                auto ptb = ptb_v[i];
                auto hltrigs = ptb.GetHLTriggers();

                if (!hltrigs.empty()){
                    if (fdebug>1) std::cout << "PTB (decoded) HLTs found: " << std::endl;
                    for (size_t j=0; j < hltrigs.size(); j++){
                        raw::ptb::Trigger trig = hltrigs.at(j);
                        if (fdebug>1){
                            std::cout << "      PTB HLT " <<  j << "-> " 
                                      << "ts (ns): " << (trig.timestamp * 20)%(uint(1e9)) 
                                      << ", trigger word: " << std::hex << trig.trigger_word << std::dec
                                      << ", word type: " << trig.word_type
                                      << std::endl;
                        }
                        if (trig.trigger_word == fptb_etrig_trigword && trig.word_type == fptb_etrig_wordtype){
                            event_trigger_time = (trig.timestamp * 20)%(uint64_t(1e9));
                            break;
                        }
                    } 
                } // end hlt check 
                else timing_type++;
            } // end of loop over ptb products
        } // end handle validity check
    } // end of timing type 1
    if (timing_type==2){
        // CAEN only configuration
        // if neither PTB or SPEC products are found, the timestamp will be equal to
        // the start of the waveform (according to CAEN TTT and wvfm length)
        event_trigger_time = 0;
    }

    auto ntrig = trig_frag_v.size();

    std::vector<uint32_t> trig_ttt_v(ntrig,0);
    std::vector<uint32_t> trig_len_v(ntrig,0);

    for (size_t i=0; i < trig_frag_v.size(); i++){
        // frag_v contains all fragments for a single trigger 
        auto frag_v = trig_frag_v.at(i);
        if (fdebug >1) std::cout << "CAEN Trigger " << i << " has " << frag_v.size() << " fragments" << std::endl;
        if (frag_v.empty()) continue;
        uint32_t trig_ttt = 0;
        uint32_t trig_len = 0;

        bool pass_check = check_fragments(frag_v, trig_ttt, trig_len);        
        if (pass_check == false)
            continue;
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

        std::vector<uint> fragid_v(ifrag_v.size(),0);

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
            fragid_v.at(idx) = frag.fragmentID() - ffragid_offset;
            iwvfm_start_v.at(idx) = get_ttt(frag) - 2*(get_length(frag));

            uint32_t frag_ttt = 0; 
            uint32_t frag_len = 0;
            int      frag_tick = 0;
            
            get_timing(frag, frag_ttt, frag_len, frag_tick);
            int frag_ts = frag_ttt - 2*(frag_len - frag_tick); // ns
            auto length = get_length(frag);
            if (fdebug>1){
                std::cout << "      Frag ID: " << frag.fragmentID() - ffragid_offset
                          << ", ttt: " << frag_ttt
                          << ", tick: " << frag_tick
                          << ", trig ts: " << frag_ts
                          << ", start ts: " << frag_ttt - 2*(frag_len)
                          << ", length: " << length
                          << std::endl;
            }
        }
        // if there exist fragments that may be part of the same trigger 
        if (extended_flag){
            for (size_t jtrig=itrig+1; jtrig < ntrig; jtrig++){
                bool pass_checks = true;
                auto jttt = trig_ttt_v.at(jtrig);
                auto jlen = trig_len_v.at(jtrig);
                auto jfrag_v = trig_frag_v.at(jtrig);

                // if the next trigger is more than 10 us away than the end of the wvfm or is equal to the nominal length, stop looking
                if (((signed)(jttt - iwvfm_end) > 1e4) || (jlen >= fnominal_length)) break; 
                else if (((jttt - iwvfm_end) < 1e4) && (jlen < fnominal_length)){
                    std::vector<std::vector<uint16_t>> jwvfm_v;
                    for (size_t idx=0; idx<jfrag_v.size(); idx++){
                        auto frag = jfrag_v.at(idx);
                        if (fragid_v.at(idx) != (uint)(frag.fragmentID() - ffragid_offset)){
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

        if (fdebug>0){
            std::cout << "Obtained waveforms for TTT with TTT " << ittt << "\n" 
                        << "\tNumber of channels: " << iwvfm_v.size() << "\n"
                        << "\tNumber of entries: " << iwvfm_v.at(0).size() << "\n"
                        << "\tTimestamp (ns): " << int(iwvfm_start_v.at(0)) - int(event_trigger_time) << std::endl;
        }
        for (size_t i = 0; i < iwvfm_v.size(); i++){
            auto combined_wvfm = iwvfm_v[i];
            int board_idx = int(i)/fnch; // integer division to get the board index

            if (evt_counter==fhist_evt){
                if (fdebug>2) std::cout << "Creating histograms for event " << evt_counter << std::endl;
                // histo: save waveforms section for combined waveforms
                histname.str(std::string());
                histname << "evt" << evt.event() << "_frag" << fragid_v.at(board_idx) << "_ch" << i%fnch << "_combined_wvfm";

                TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), combined_wvfm.size(), 0, combined_wvfm.size());
                wvfmHist->GetXaxis()->SetTitle("ticks");
                for(unsigned int n = 0; n < combined_wvfm.size(); n++) 
                    wvfmHist->SetBinContent(n + 1, (double)combined_wvfm[n]);
            }
            int time_diff = int(iwvfm_start_v.at(board_idx)) - int(event_trigger_time);
            uint ch;
            if (i%fnch == 15)
                ch = fragid_v.at(board_idx);
            else
                ch = fch_map.at(fragid_v.at(board_idx)*15 + i%fnch);
            
            raw::OpDetWaveform waveform(time_diff, ch, combined_wvfm);

            if (i%fnch == 15){
                if (foutput_ftrig_wvfm) twvfmVec->push_back(waveform);
            }
            else 
                wvfmVec->push_back(waveform);
        }

    } // loop over triggers 
    trig_frag_v.clear(); 

    evt.put(std::move(wvfmVec),fch_instance_name);  
    evt.put(std::move(twvfmVec),ftr_instance_name);
}

bool sbndaq::SBNDPMTDecoder::check_fragments(std::vector<artdaq::Fragment> &frag_v, uint32_t& trig_ttt, uint32_t& trig_len){
    bool pass = true;
    uint32_t this_ttt = 0; 
    uint32_t this_len = 0;
    int frag_counter = 0;
    for (size_t ifrag=0; ifrag<frag_v.size(); ifrag++){

        // assuming that the timing CAEN has fragmentId==8, skip it 
        if ((frag_v.at(ifrag).fragmentID() - ffragid_offset) == 8 )
            continue;
        CAENV1730Fragment bb(frag_v.at(ifrag));
        CAENV1730Event const* event_ptr = bb.Event();
        auto const * md = bb.Metadata();
        CAENV1730EventHeader header = event_ptr->Header;

        uint32_t ttt = header.triggerTimeTag*8;
        auto nch = md->nChannels;    
        fnch = nch; // set the # of channels for the event
        uint32_t ev_size_quad_bytes = header.eventSize;
        uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
        uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
        uint32_t wvfm_length = data_size_double_bytes/nch; // downsampled trigger time tag; TTT points to end of wvfm w.r.t. pps

        if (frag_counter==0){ 
            this_ttt = ttt; 
            this_len = wvfm_length;
        }
        frag_counter++;
        if (this_ttt != ttt || this_len != wvfm_length){
            if (abs((signed)(this_ttt - ttt)) > 16 || abs((signed)(this_len - wvfm_length)) > 8 ){
                std::cout << "Mismatch is greater than maximum 16 ns jitter" << std::endl;
                pass=false;
                break;
            }
        }
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
