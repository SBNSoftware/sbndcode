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

#include "TTree.h"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

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

    // Declare member data here.

    // Event Tree
    TTree *fTree;
    int _run, _subrun, _event;

    // TDC stuff
    std::vector<uint64_t> fTDC_ch0;
    std::vector<uint64_t> fTDC_ch1;
    std::vector<uint64_t> fTDC_ch2;
    std::vector<uint64_t> fTDC_ch3;
    std::vector<uint64_t> fTDC_ch4;

    std::vector<uint64_t> fTDC_ch0_utc;
    std::vector<uint64_t> fTDC_ch1_utc;
    std::vector<uint64_t> fTDC_ch2_utc;
    std::vector<uint64_t> fTDC_ch3_utc;
    std::vector<uint64_t> fTDC_ch4_utc;

    // Product label
    art::InputTag fTDCDecodeLabel;

    bool fDebug;
};


sbnd::BeamAnalysis::BeamAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // 
    // More initializers here.
{
    fTDCDecodeLabel = p.get<art::InputTag>("TDCDecodeLabel", "tdcdecoder");
    fDebug = p.get<bool>("Debug", true);
  
}

void sbnd::BeamAnalysis::analyze(art::Event const& e)
{
    
    ResetEventVars();

    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  =  e.id().event();

    art::Handle<std::vector<sbnd::timing::DAQTimestamp>> TDCHandle;
    e.getByLabel(fTDCDecodeLabel, TDCHandle);

    if (!TDCHandle.isValid() || TDCHandle->size() == 0){

        if (fDebug) std::cout << "No SPECTDC products found." << std::endl;
    }
    else{
        
        const std::vector<sbnd::timing::DAQTimestamp> tdc_v(*TDCHandle);
        
        for (size_t i=0; i<tdc_v.size(); i++){
            auto tdc = tdc_v[i];
            const uint32_t  ch = tdc.Channel();
            const uint64_t  ts = tdc.Timestamp();
            //const uint64_t  offset = tdc.Offset();
            const std::string name  = tdc.Name();

            if (fDebug) std::cout << "Event " << _event 
                                << ": ch " << name
                                << ", ts (ns) = " << ts%uint64_t(1e9)
                                << ", sec (s) = " << ts/uint64_t(1e9)
                                << std::endl;

            if(ch == 0){
                fTDC_ch0_utc.push_back(ts);
                fTDC_ch0.push_back(ts%uint64_t(1e9));
            }
            if(ch == 1){
                fTDC_ch1_utc.push_back(ts);
                fTDC_ch1.push_back(ts%uint64_t(1e9));
            }
            if(ch == 2){
                fTDC_ch2_utc.push_back(ts);
                fTDC_ch2.push_back(ts%uint64_t(1e9));
            }
            if(ch == 3){
                fTDC_ch3_utc.push_back(ts);
                fTDC_ch3.push_back(ts%uint64_t(1e9));
            }
            if(ch == 4){
                fTDC_ch4_utc.push_back(ts);
                fTDC_ch4.push_back(ts%uint64_t(1e9));
            }
        } 
    }

    //Fill once every event
    fTree->Fill();
}

void sbnd::BeamAnalysis::beginJob()
{
 
    art::ServiceHandle<art::TFileService> tfs;

    //Event Tree
    fTree = tfs->make<TTree>("events", "");
  
    fTree->Branch("run", &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event", &_event);
    fTree->Branch("TDC_ch0", &fTDC_ch0);
    fTree->Branch("TDC_ch1", &fTDC_ch1);
    fTree->Branch("TDC_ch2", &fTDC_ch2);
    fTree->Branch("TDC_ch3", &fTDC_ch3);
    fTree->Branch("TDC_ch4", &fTDC_ch4);
    fTree->Branch("TDC_ch0_utc", &fTDC_ch0_utc);
    fTree->Branch("TDC_ch1_utc", &fTDC_ch1_utc);
    fTree->Branch("TDC_ch2_utc", &fTDC_ch2_utc);
    fTree->Branch("TDC_ch3_utc", &fTDC_ch3_utc);
    fTree->Branch("TDC_ch4_utc", &fTDC_ch4_utc);
}

void sbnd::BeamAnalysis::endJob()
{
  // Implementation of optional member function here.
}

void sbnd::BeamAnalysis::ResetEventVars()
{
    _run = -1; _subrun = -1; _event = -1;
   
    fTDC_ch0.clear();
    fTDC_ch1.clear();
    fTDC_ch2.clear();
    fTDC_ch3.clear();
    fTDC_ch4.clear();

    fTDC_ch0_utc.clear();
    fTDC_ch1_utc.clear();
    fTDC_ch2_utc.clear();
    fTDC_ch3_utc.clear();
    fTDC_ch4_utc.clear();
}

DEFINE_ART_MODULE(sbnd::BeamAnalysis)
