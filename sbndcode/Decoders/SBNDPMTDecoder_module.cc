////////////////////////////////////////////////////////////////////////
// Class:       SBNDPMTDecoder
// Plugin Type: producer (Unknown Unknown)
// File:        SBNDPMTDecoder_module.cc
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

#include "art_root_io/TFileService.h"
#include "TH1F.h"
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
    struct Config {
        fhicl::Atom<double> ttt_post_percent {
            fhicl::Name("ttt_post_percent"),
            fhicl::Comment("the percent of the waveform kept after the trigger"),
            0.9
        };
        fhicl::Atom<int> ttt_downsample {
            fhicl::Name("ttt_downsample"),
            fhicl::Comment("the downsampling factor from clock counts to ticks, default=4"),
            4
        };
        fhicl::Atom<int> waveform_length {
            fhicl::Name("waveform_length"),
            fhicl::Comment("the expected length of the waveform, default=5000 ticks. Used to differentiate shortened waveforms"),
            5000
        };
        fhicl::Atom<bool> verbose {
            fhicl::Name("verbose"),
            fhicl::Comment("toggle for additional text printout"),
            true
        };
        fhicl::Atom<std::string> ch_instance_name {
            fhicl::Name("ch_instance_name"),
            fhicl::Comment("instance name for pmt channel wvfm data product"),
            "PMTChannels"
        };
        fhicl::Atom<std::string> tr_instance_name {
            fhicl::Name("tr_instance_name"),
            fhicl::Comment("instance name for trigger channels wvfm data product"),
            "TriggerChannels"
        };
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit SBNDPMTDecoder(Parameters const& config);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SBNDPMTDecoder(SBNDPMTDecoder const&) = delete;
    SBNDPMTDecoder(SBNDPMTDecoder&&) = delete;
    SBNDPMTDecoder& operator=(SBNDPMTDecoder const&) = delete;
    SBNDPMTDecoder& operator=(SBNDPMTDecoder&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;
    void analyze_caen_fragment(artdaq::Fragment &frag);
    
    std::vector<raw::OpDetWaveform> fWaveforms;
    std::vector<raw::OpDetWaveform> fTriggerWaveforms;


private:
    double fttt_post_percent;
    int fttt_downsample;
    int fwaveform_length;
    bool fverbose;
    std::string fch_instance_name;
    std::string ftr_instance_name;
};


sbndaq::SBNDPMTDecoder::SBNDPMTDecoder(Parameters const& config)
    : EDProducer{config}  // ,
{
    fttt_post_percent = config().ttt_post_percent();
    fttt_downsample = config().ttt_downsample();
    fwaveform_length = config().waveform_length();
    fverbose = config().verbose();
    fch_instance_name = config().ch_instance_name();
    ftr_instance_name = config().tr_instance_name();

    produces< std::vector< raw::OpDetWaveform > >(fch_instance_name); 
    produces< std::vector< raw::OpDetWaveform > >(ftr_instance_name);
}

void sbndaq::SBNDPMTDecoder::produce(art::Event& evt)
{
    std::unique_ptr< std::vector< raw::OpDetWaveform > > wvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    std::unique_ptr< std::vector< raw::OpDetWaveform > > twvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());

    // auto run = evt.run();
    // auto sub = evt.subRun();
    // auto evt = evt.event();

    std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
    fragmentHandles = evt.getMany<std::vector<artdaq::Fragment>>();
    
    fWaveforms.clear();
    fTriggerWaveforms.clear();

    for (auto handle : fragmentHandles){
        if (!handle.isValid() || handle->size() == 0){
            std::cout << "No CAEN V1730 fragments found." << std::endl;
            continue;
        }
        if (handle->front().type() == artdaq::Fragment::ContainerFragmentType){
            for (auto cont : *handle){
            artdaq::ContainerFragment contf(cont);
            if (contf.fragment_type() == sbndaq::detail::FragmentType::CAENV1730){
                std::cout << "Found " << contf.block_count() << " CAEN Fragments in container." << std::endl;
                for (size_t ii = 0; ii < contf.block_count(); ++ii){
                analyze_caen_fragment(*contf[ii].get());
                }
            }
            }
        }
        else{
            if (handle->front().type() == sbndaq::detail::FragmentType::CAENV1730){
                std::cout << "Found normal CAEN fragments." << std::endl;
                for (auto frag : *handle){
                    analyze_caen_fragment(frag);
                }
            }
        }
    }
    if (fWaveforms.empty() == false){
        std::cout << "Number of PMT waveforms: " << fWaveforms.size() << std::endl;
        std::cout << "Number of trigger waveforms: " << fTriggerWaveforms.size() << std::endl;
        for (const raw::OpDetWaveform &waveform : fWaveforms)
            wvfmVec->push_back(waveform);
        for (const raw::OpDetWaveform &waveform : fTriggerWaveforms)
            twvfmVec->push_back(waveform);
    }
    else 
        std::cout << "pushing empty waveform vec!" << std::endl;
    
    evt.put(std::move(wvfmVec),fch_instance_name);  
    evt.put(std::move(twvfmVec),ftr_instance_name);

}


void sbndaq::SBNDPMTDecoder::analyze_caen_fragment(artdaq::Fragment & frag)  {


    CAENV1730Fragment bb(frag);
    auto const* md = bb.Metadata();
    CAENV1730Event const* event_ptr = bb.Event();

    CAENV1730EventHeader header = event_ptr->Header;
    
    // fragID = static_cast<int>(frag.fragmentID()); 
    
    std::cout << "\tFrom header, event counter is "  << header.eventCounter   << "\n";
    // std::cout << "\tFrom header, triggerTimeTag is " << header.triggerTimeTag << "\n";
    // std::cout << "\tFrom fragment, fragment id is "  << fragID << "\n";
    
    auto boardID = header.boardID;

    uint32_t ttt = header.triggerTimeTag*4; // downsampled trigger time tag; TTT points to end of wvfm w.r.t. pps
    auto nch = md->nChannels;

    std::cout << "From header, downsampled TTT is " << ttt << "\n";
    std::cout << "\tNumber of channels: " << nch << "\n";
    
    //--get the number of 32-bit words (quad_bytes) from the header
    uint32_t ev_size_quad_bytes = header.eventSize;
    std::cout << "\tEvent size in quad bytes is: " << ev_size_quad_bytes << "\n";
    uint32_t evt_header_size_quad_bytes = sizeof(CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t wvfm_length = data_size_double_bytes/nch;
    std::cout << "\tChannel waveform length = " << wvfm_length << "\n";

    auto trigger_position = ttt - fttt_post_percent*wvfm_length; // ticks, position of trigger in waveform
    auto trigger_time = trigger_position*2*0.001; // (2 ns per tick) * (1 us per 1000 ns)
    std::cout << "\tTrigger time is " << trigger_time << " us" << std::endl; 

    // --store the tick value for each acquisition 
    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
								   + sizeof(CAENV1730EventHeader));
    const uint16_t* value_ptr =  data_begin;
    uint16_t value = 0;
    size_t ch_offset = 0;
    //--loop over channels
    for (size_t i_ch=0; i_ch<nch; ++i_ch){
        raw::OpDetWaveform wvfm;
        wvfm.Waveform() = std::vector<short>(wvfm_length,0);
        wvfm.SetChannelNumber(boardID*100 + (uint)i_ch); // dummy channel number: boardID*100 + channel number
        wvfm.SetTimeStamp(trigger_time);
        ch_offset = (size_t)(i_ch * wvfm_length);      
        //--loop over waveform samples
        for(size_t i_t=0; i_t<wvfm_length; ++i_t){ 
            value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
            value = *(value_ptr);

            wvfm.at(i_t) = value;
        }
        if (i_ch == 15)
            fTriggerWaveforms.push_back(wvfm);
        else
        fWaveforms.push_back(wvfm);
    }
}

DEFINE_ART_MODULE(sbndaq::SBNDPMTDecoder)
