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
    void get_fragments(artdaq::Fragment & frag, std::vector<std::vector<artdaq::Fragment>> & board_frag_v);
    void get_waveforms(artdaq::Fragment & frag, std::vector<std::vector<uint16_t>> & wvfm_v);

    uint32_t get_length(artdaq::Fragment & frag);
    uint32_t get_ttt(artdaq::Fragment & frag);
    uint32_t get_boardid(artdaq::Fragment & frag);
    void     get_timing(artdaq::Fragment & frag, uint32_t & ttt, uint32_t & len, int & tick);
    

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
    uint fn_caenboards;
    uint16_t fthreshold_ftrig;

    uint ffragid_offset;
    uint fhist_evt;

    std::vector<uint> fset_fragid_map;
    bool fuse_set_map;


    std::vector<uint> fch_map;

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
    fn_caenboards    = p.get<uint>("n_caenboards",8);
    fthreshold_ftrig = p.get<uint16_t>("threshold_ftrig",16350);
    ffragid_offset   = p.get<uint>("fragid_offset",40960);
    fhist_evt        = p.get<int>("hist_evt",1);

    fset_fragid_map = p.get<std::vector<uint>>("set_fragid_map",{});
    fuse_set_map    = p.get<bool>("use_set_map",false);

    fch_map          = p.get<std::vector<uint>>("ch_map",{});

    produces< std::vector< raw::OpDetWaveform > >(fch_instance_name); 
    produces< std::vector< raw::OpDetWaveform > >(ftr_instance_name);
}

void sbndaq::SBNDPMTDecoder::produce(art::Event& evt)
{
    std::unique_ptr< std::vector< raw::OpDetWaveform > > wvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    std::unique_ptr< std::vector< raw::OpDetWaveform > > twvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());

    evt_counter++;

    std::vector<std::vector<artdaq::Fragment>> board_frag_v(fn_caenboards);
    uint64_t event_trigger_time = 0; 

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
            
                if ((ch==fspectdc_etrig_ch || ch==fspectdc_ftrig_ch) && (fdebug>1)){
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
        if (fdebug>1) std::cout << "Board " << iboard << " has " << board_frag_v.at(iboard).size() << " trigger(s) || " ;
    }
    if (fdebug>1) std::cout << std::endl;

    // store the waveform itself
    // store the waveform start time (for OpDetWaveform timestamp)
    for (size_t iboard=0; iboard < board_frag_v.size(); iboard++){
        // the vector of fragments from a single digitizer 
        auto frag_v = board_frag_v.at(iboard);
        if (std::find(fignore_fragid.begin(), fignore_fragid.end(), iboard) != fignore_fragid.end())
            continue;
        if (frag_v.empty()) continue;
        auto fragid = iboard;
        auto trig_counter = 0;

        for (size_t itrig=0; itrig < frag_v.size(); itrig++){
            std::vector<std::vector<uint16_t>> iwvfm_v;
            auto ifrag = frag_v.at(itrig);
            uint32_t ittt = 0; 
            uint32_t ilen = 0;
            int      itick = 0;
            
            get_timing(ifrag, ittt, ilen, itick);
            // if this waveform is short, skip it 
            if (ilen < fnominal_length) continue;

            get_waveforms(ifrag, iwvfm_v);
            auto iwvfm_start = ittt - 2*ilen;
            auto iwvfm_end   = ittt;

            if (fdebug>1){
                std::cout << "      Frag ID: " << fragid
                          << ", ttt: " << ittt
                          << ", tick: " << itick
                          << ", frag ts: " << ifrag.timestamp()%uint(1e9)
                          << ", length: " << ilen
                          << std::endl;
            }
            if (extended_flag){
                for (size_t jtrig=itrig+1; jtrig < frag_v.size(); jtrig++){
                    auto jfrag = frag_v.at(jtrig);
                    auto jlen = get_length(jfrag);
                    auto jttt = get_ttt(jfrag);
                    // if the next trigger is more than 10 us away than the end of the wvfm or is equal to the nominal length, stop looking
                    if (((signed)(jttt - iwvfm_end) > (signed)(fnominal_length*2)) || (jlen >= fnominal_length)) break; 
                    else if (((jttt - iwvfm_end) < 1e4) && (jlen < fnominal_length)){
                        std::vector<std::vector<uint16_t>> jwvfm_v;
                        get_waveforms(jfrag, jwvfm_v);
                        for (size_t ich=0; ich < iwvfm_v.size(); ich++){
                            std::vector<uint16_t>& orig_wvfm = iwvfm_v.at(ich);
                            std::vector<uint16_t>  extd_wvfm = jwvfm_v.at(ich); // extended waveform 
                            orig_wvfm.insert( orig_wvfm.end(), extd_wvfm.begin(), extd_wvfm.end());
                        } // end ch loop 
                    } // end if extension condition
                    iwvfm_end = jttt;
                } // end jtrig loop
                if (fdebug>2){
                    std::cout << "      Frag ID: " << fragid
                              << " -> start time: " << int(iwvfm_start) - int(event_trigger_time)
                              << " , wvfm length: " << iwvfm_v.at(0).size()
                              << std::endl; 
                }
            } // end extended flag         
            for (size_t i = 0; i < iwvfm_v.size(); i++){
                auto combined_wvfm = iwvfm_v[i];
                if (evt_counter==fhist_evt){                    // histo: save waveforms section for combined waveforms
                    histname.str(std::string());
                    histname << "evt" << evt.event() << "_trig" << trig_counter << "_frag" << fragid << "_ch" << i << "_combined_wvfm";

                    TH1D *wvfmHist = tfs->make< TH1D >(histname.str().c_str(), histname.str().c_str(), combined_wvfm.size(), 0, combined_wvfm.size());
                    wvfmHist->GetXaxis()->SetTitle("ticks");
                    for(unsigned int n = 0; n < combined_wvfm.size(); n++) 
                        wvfmHist->SetBinContent(n + 1, (double)combined_wvfm[n]);
                }
                double time_diff = (int(iwvfm_start) - int(event_trigger_time))*1e-3; // us
                uint ch;
                if (i == 15){
                    ch = fragid;
                }
                else
                    ch = fch_map.at(fragid*15 + i);
                raw::OpDetWaveform waveform(time_diff, ch, combined_wvfm);

                if (i == 15){
                    if (foutput_ftrig_wvfm) twvfmVec->push_back(waveform);
                }
                else 
                    wvfmVec->push_back(waveform);
            }
            trig_counter++;

        } // end itrig loop 
    } // end board loop
    board_frag_v.clear();

    evt.put(std::move(wvfmVec),fch_instance_name);  
    evt.put(std::move(twvfmVec),ftr_instance_name);
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
    // ignore boards that are not in the list of boards to ignore
    if (std::find(fignore_fragid.begin(), fignore_fragid.end(), fragid) != fignore_fragid.end())
        return;
    // check our fragid, is it reasonable?
    if (fragid<0 || fragid >= fn_caenboards){
        std::cout << "Fragment ID " << fragid << " is out of range. FragID offset/FragID map may be misconfigured, or this FragID is not attributed to a PMT Digitizer. Skipping this fragment..." << std::endl;
        return;
    }
    board_frag_v.at(fragid).push_back(frag);
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
