////////////////////////////////////////////////////////////////////////
// Class:       ArtdaqFragmentProducer
// Module Type: producer
// File:        ArtdaqFragmentProducer_module.cc
//
// Producer module for creating DAQ-readable fragments from simulation
//
//
// Erin Yandel (eyandel@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"


// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// artdaq includes
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h"

// random numbers
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

// ROOT includes.
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TString.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <random>
#include <iomanip>

namespace sbnd {
  namespace trigger {
    class ArtdaqFragmentProducer;
  }
}

class sbnd::trigger::ArtdaqFragmentProducer : public art::EDProducer {
public:
  explicit ArtdaqFragmentProducer(fhicl::ParameterSet const& p);

  ArtdaqFragmentProducer(ArtdaqFragmentProducer const&) = delete;
  ArtdaqFragmentProducer(ArtdaqFragmentProducer&&) = delete;
  ArtdaqFragmentProducer & operator=(ArtdaqFragmentProducer const&) = delete;
  ArtdaqFragmentProducer & operator=(ArtdaqFragmentProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:

  int fRun;
  int fSubrun;
  int fEvent;

  //CRT

  // fcl file parameters
  art::InputTag fSimModuleLabel;      ///< name of detsim producer
  art::InputTag fCRTSimLabel;         ///< name of CRT producer
  std::string fFEBDataLabel;          ///< name of FEBData producer
  bool          fVerbose;             ///< print information about what's going on
  double        fClockSpeedCRT;
  size_t        fFirstFEBMac5;        ///< lowest mac5 address for CRT FEBs

  //limits for array sizes
  enum LIMITS{
    max_hits_in_fragment = 50,
    num_febs = 255,
    num_channels = 32
  };

  //arrays to store FEBinfo by module (want one fragment per module per event)
  std::string taggers[num_febs];
  uint16_t feb_hits_in_fragments[num_febs];
  uint16_t ADCs[num_febs][max_hits_in_fragment][num_channels];
  uint32_t T0s[num_febs][max_hits_in_fragment];
  uint32_t T1s[num_febs][max_hits_in_fragment];
  bool empty_fragment[num_febs];


  // Other variables shared between different methods.
  geo::GeometryCore const* fGeometryService;
  sbnd::crt::CRTGeoAlg fCrtGeo;

  //PMT

  // fhicl parameters
  std::string fInputModuleNameWvfm;
  std::string fInputModuleNameTrigger;
  int fBaseline; // baseline in simulation, default 8000 ADC (for expanding waveforms only, when not fully simulated)
  int fMultiplicityThreshold; // number of PMT pairs in hardware trigger to pass
  double fBeamWindowLength;
  uint32_t nChannelsFrag;
  uint32_t wfm_length; // ~10us, 2ns tick

  // PD information
  opdet::sbndPDMapAlg pdMap; // photon detector map
  std::vector<unsigned int> channelList;

  // waveforms
  std::vector<std::vector<short>> wvf_channel;

  // sampling rate
  double fSampling;

  // services
  art::ServiceHandle<art::TFileService> tfs;

  // random numbers
  CLHEP::HepRandomEngine&       fTriggerTimeEngine;

}; // class ArtdaqFragmentProducer


sbnd::trigger::ArtdaqFragmentProducer::ArtdaqFragmentProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fSimModuleLabel(p.get<std::string>("SimModuleLabel", "largeant")),
  fCRTSimLabel(p.get<std::string>("CRTSimLabel", "crt")),
  fFEBDataLabel(p.get<std::string>("FEBDataLabel", "crtsim")),
  fVerbose(p.get<bool>("Verbose", false)),
  fClockSpeedCRT(p.get<double>("ClockSpeedCRT")),
  fFirstFEBMac5(p.get<size_t>("FirstFEBMac5", 0)),
  fInputModuleNameWvfm(p.get<std::string>("InputModuleNameWvfm")),
  fInputModuleNameTrigger(p.get<std::string>("InputModuleNameTrigger")),
  fBaseline(p.get<int>("Baseline",8000)),
  fMultiplicityThreshold(p.get<int>("MultiplicityThreshold")),
  fBeamWindowLength(p.get<double>("BeamWindowLength", 1.6)),
  nChannelsFrag(p.get<double>("nChannelsFrag", 15)),
  wfm_length(p.get<double>("WfmLength", 5120)),
    fTriggerTimeEngine(art::ServiceHandle<rndm::NuRandomService>{}->registerAndSeedEngine(
                         createEngine(0, "HepJamesRandom", "trigger"), "HepJamesRandom", "trigger", p, "SeedTriggerTime"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<artdaq::Fragment> >();

  // Get a pointer to the fGeometryServiceetry service provider
  fGeometryService = lar::providerFrom<geo::Geometry>();

  // get clock
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSampling = clockData.OpticalClock().Frequency(); // MHz

  // build PD map and channel list
  auto subsetCondition = [](auto const& i)->bool { return i["pd_type"] == "pmt_coated" || i["pd_type"] == "pmt_uncoated"; };
  auto pmtMap = pdMap.getCollectionFromCondition(subsetCondition);
  if (fVerbose) std::cout << "Number of PDs selected: \t" << pmtMap.size() << "\n";
  for(auto const& i:pmtMap){
    channelList.push_back(i["channel"]);
  }

}


void sbnd::trigger::ArtdaqFragmentProducer::produce(art::Event& e)
{

    // Fetch basic event info
    // event information
    fRun = e.run();
    fSubrun = e.subRun();
    fEvent = e.id().event();
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<fRun<<", SubRun = "<<fSubrun<<", Event = "<<fEvent<<std::endl
               <<"============================================"<<std::endl;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          CRT
    //----------------------------------------------------------------------------------------------------------

    for(int i=0; i<num_febs; i++){empty_fragment[i]=true;}

    int num_module = fCrtGeo.NumModules();

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get FEB data from the event
    art::Handle<std::vector<sbnd::crt::FEBData>> feb_data_h;
    e.getByLabel(fFEBDataLabel, feb_data_h);

    // make sure hits look good
    if (!feb_data_h.isValid()) {
      throw art::Exception(art::errors::Configuration) << "could not locate FEBData." << std::endl;;
    }

    std::vector<art::Ptr<sbnd::crt::FEBData>> feb_data_v;
    art::fill_ptr_vector(feb_data_v, feb_data_h);

    art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febdata_to_ides (feb_data_h, e, fFEBDataLabel);

    //Fill arrays with default values
    std::fill(feb_hits_in_fragments, feb_hits_in_fragments + sizeof(feb_hits_in_fragments), 0);

    // fragments vector
    std::unique_ptr<std::vector<artdaq::Fragment>> vecFrag = std::make_unique<std::vector<artdaq::Fragment>>();

    // set properties of fragment that are common to event
    //quantities in fragment
    uint16_t lostcpu = 0;
    uint16_t lostfpga = 0;
    uint32_t coinc = 0;
    uint16_t lost_hits = 0; //number of lost hits from the previous one


  //metadata
    uint64_t  run_start_time = 0;
    uint64_t  this_poll_start = 0;
    uint64_t  this_poll_end = 0;
    int32_t   system_clock_deviation = 0;
    uint32_t  feb_hits_in_poll = 0;

  //information from fragment header
    uint64_t  sequence_id = fEvent;

    uint64_t temp_last_time = 0;

    //random number generator to bu used to simulate pedestal
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(175.0,23.0);

    //loop through FEBData product and sort by module/mac5
    for (size_t feb_i = 0; feb_i < feb_data_v.size(); feb_i++) {

        auto const feb_data = feb_data_v[feb_i];

        if(fVerbose){std::cout << "FEB " << feb_i << " with mac " << feb_data->Mac5() << std::endl;}

        empty_fragment[feb_data->Mac5()] = false;

        int channel = feb_data->Mac5() * 32;
        std::string stripName = fCrtGeo.ChannelToStripName(channel);
        std::string tagger = fCrtGeo.GetTaggerName(stripName);
        taggers[feb_data->Mac5()] = tagger;

        T0s[feb_data->Mac5()][feb_hits_in_fragments[feb_data->Mac5()]] = feb_data->Ts0();
        T1s[feb_data->Mac5()][feb_hits_in_fragments[feb_data->Mac5()]] = feb_data->Ts1();
        for (int i_adc = 0; i_adc<32; i_adc++){
          uint16_t adc = 0;
          if (feb_data->ADC()[i_adc] > 4089){
            adc = 4089;
          }else if (feb_data->ADC()[i_adc] == 0){
            //pull from a normal distribution to simulate pedestal
            adc = distribution(generator);
          }else{
            adc = feb_data->ADC()[i_adc];
          }
          ADCs[feb_data->Mac5()][feb_hits_in_fragments[feb_data->Mac5()]][i_adc] = adc;
        }

        feb_hits_in_fragments[feb_data->Mac5()]++;

    }//FEBData loop

    //make one fragment for every module
    for (size_t feb_i = fFirstFEBMac5; feb_i < num_module+fFirstFEBMac5; feb_i++){
        //quantities in fragment
        if (empty_fragment[feb_i]){
          //if no hits for a module, make a simulated "T1 reset" event to avoid missing fragments
          feb_hits_in_fragments[feb_i] = 1;
          int channel = feb_i * 32;
          std::string stripName = fCrtGeo.ChannelToStripName(channel);
          std::string tagger = fCrtGeo.GetTaggerName(stripName);
          taggers[feb_i] = tagger;

          T0s[feb_i][0] = 0;
          T1s[feb_i][0] = 0;
          for (int i_adc = 0; i_adc<32; i_adc++){
            uint16_t adc = distribution(generator);
            ADCs[feb_i][0][i_adc] = adc;
          }
        }

        if (feb_hits_in_fragments[feb_i]==0) continue;

      //metadata
        uint8_t  mac5 = (uint8_t)feb_i; //last 8 bits of FEB mac5 address
        uint16_t feb_hits_in_fragment = feb_hits_in_fragments[feb_i];

        sbndaq::BernCRTFragmentMetadataV2 metadata;

        metadata.set_mac5(mac5);
        metadata.set_clock_deviation(system_clock_deviation);
        metadata.set_hits_in_poll(feb_hits_in_poll);
        metadata.set_hits_in_fragment(feb_hits_in_fragment);
        metadata.set_run_start_time(run_start_time);
        metadata.update_poll_time(this_poll_start, this_poll_end);

        uint64_t  timestamp = (uint64_t)(T0s[feb_i][feb_hits_in_fragments[feb_i]]/fClockSpeedCRT); //absolute timestamp

        // create fragment
        sbnd::crt::CRTTagger tagger_num = sbnd::crt::CRTCommonUtils::GetTaggerEnum(taggers[feb_i]);
        uint16_t fragmentIDVal = 32768 + 12288 + (tagger_num * 256) + (uint16_t)mac5;
        if(fVerbose){std::cout<<"fragmentID: "<<std::bitset<16>{fragmentIDVal}<<std::endl;}
        auto fragment_uptr = artdaq::Fragment::FragmentBytes(sizeof(sbndaq::BernCRTHitV2)*metadata.hits_in_fragment(), //payload_size
            sequence_id,
            fragmentIDVal,
            sbndaq::detail::FragmentType::BERNCRTV2,
            metadata,
            timestamp
          ); // unique pointer

        // populate fragment
        std::vector<sbndaq::BernCRTHitV2> data_vec;
        for (int i_frag = 0; i_frag<feb_hits_in_fragments[feb_i]; i_frag++){
          sbndaq::BernCRTHitV2 hit;

          uint8_t flags = 3;
          uint32_t ts0 = T0s[feb_i][i_frag];
          uint32_t ts1 = T1s[feb_i][i_frag];

          if (ts0==0){flags=7;}else if (ts1==0){flags=11;}

          uint64_t  feb_hit_number = feb_hits_in_fragments[feb_i]; //hit counter for individual FEB, including hits lost in FEB or fragment generator
          timestamp = (uint64_t)ts0/fClockSpeedCRT; //absolute timestamp
          uint64_t  last_accepted_timestamp = temp_last_time; //timestamp of previous accepted hit
          temp_last_time = timestamp;

          hit.flags = (uint16_t)flags;
          hit.lostcpu = (uint16_t)lostcpu;
          hit.lostfpga = (uint16_t)lostfpga;
          hit.ts0 = (uint32_t)ts0;
          hit.ts1 = (uint32_t)ts1;
          for (int i_adc = 0; i_adc<32; i_adc++){
            hit.adc[i_adc] = ADCs[feb_i][i_frag][i_adc];
          }
          hit.coinc = (uint32_t)coinc;
          hit.feb_hit_number = (uint64_t)feb_hit_number;
          hit.timestamp = (uint64_t)timestamp;
          hit.last_accepted_timestamp = (uint64_t)last_accepted_timestamp;
          hit.lost_hits = (uint16_t)lost_hits;

          data_vec.emplace_back(hit);
        }//bern crt hit vector

        //copy data vector of hits into the fragment we made
        std::memcpy(fragment_uptr->dataBeginBytes(), data_vec.data(), sizeof(sbndaq::BernCRTHitV2)*metadata.hits_in_fragment());

        // add fragment to vector
        vecFrag->push_back(*fragment_uptr);

    }//module (mac5) loop

    int num_crt_frags = vecFrag->size();
    if(fVerbose) std::cout << "CRT Fragments written: " << num_crt_frags << std::endl;

    //----------------------------------------------------------------------------------------------------------
    //                                          PMT
    //----------------------------------------------------------------------------------------------------------

    // access PMT waveforms and hardware trigger information
  art::Handle< std::vector< raw::OpDetWaveform > > wvfmHandle;
  art::Handle< std::vector< sbnd::comm::pmtTrigger > > triggerHandle;
  e.getByLabel(fInputModuleNameWvfm, wvfmHandle);
  e.getByLabel(fInputModuleNameTrigger, triggerHandle);

  if(!wvfmHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any waveforms, input must contain OpDetWaveforms." << "\n";
  }
  if(!triggerHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any PMT hardware trigger object, must run PMT trigger producer before running this module." << "\n";
  }

  // create empty vectors to hold waveforms for each channel
  double fMinStartTime = -1510.0;//in us
  double fMaxEndTime = 1510.0;//in us

  for(auto const& wvf : (*wvfmHandle)) {
    double fChNumber = wvf.ChannelNumber();
    std::string opdetType = pdMap.pdType(fChNumber);
    // only look at pmts
    if (opdetType != "pmt_coated" && opdetType != "pmt_uncoated") continue;
      if (wvf.TimeStamp() < fMinStartTime){ fMinStartTime = wvf.TimeStamp(); }
      if ((double(wvf.size()) / fSampling + wvf.TimeStamp()) > fMaxEndTime){ fMaxEndTime = double(wvf.size()) / fSampling + wvf.TimeStamp();}
  }

  if (fVerbose){std::cout<<"MinStartTime: "<<fMinStartTime<<" MaxEndTime: "<<fMaxEndTime<<std::endl;}

  std::vector<short> wvf_0; wvf_0.reserve((int)(3020*1e6/2));
  for (double i = fMinStartTime; i<fMaxEndTime+(1./fSampling); i+=(1./fSampling)){
    wvf_0.push_back(fBaseline);
  }
  for (size_t i = 0; i < channelList.size(); i++){
    wvf_channel.push_back(wvf_0);
  }

  // counters
  int num_pmt_wvf = 0;

  // loop through waveform handles
  size_t wvf_id = -1;
  size_t hist_id = -1;
  for(auto const& wvf : (*wvfmHandle)) {
    wvf_id++;
    double fChNumber = wvf.ChannelNumber();
    std::string opdetType = pdMap.pdType(fChNumber);

    // only look at pmts
    if (opdetType != "pmt_coated" && opdetType != "pmt_uncoated") continue;

    num_pmt_wvf++;

    double fStartTime = wvf.TimeStamp(); // in us
    double fEndTime = double(wvf.size()) / fSampling + fStartTime; // in us

    // create full waveform
    std::vector<short> wvf_full; wvf_full.reserve((short)(3020*1e6/2));

    if (fStartTime > fMinStartTime){
      for (double i = fStartTime-fMinStartTime; i>0.; i-=(1./fSampling)){
        wvf_full.push_back(fBaseline);
      }
    }

    for(unsigned int i = 0; i < wvf.size(); i++) {
      wvf_full.push_back(wvf[i]);
    }

    if (fEndTime < fMaxEndTime){
      for (double i = fMaxEndTime-fEndTime; i>0.; i-=(1./fSampling)){
        wvf_full.push_back(fBaseline);
      }
    }

    // combine waveform with any other waveforms from same channel
    int i_ch = -1.;
    auto ich = std::find(channelList.begin(), channelList.end(), fChNumber);
    if (ich != channelList.end()){
      i_ch = ich - channelList.begin();
    }
    if (wvf_channel.at(i_ch).size() < wvf_full.size()){
      if (fVerbose) std::cout<<"Full waveform -- Previous Channel" << fChNumber <<" Size: "<<wvf_channel.at(i_ch).size()<<"New Channel" << fChNumber <<" Size: "<<wvf_full.size()<<std::endl;
      for(unsigned int i = wvf_channel.at(i_ch).size(); i < wvf_full.size(); i++) {
        wvf_channel.at(i_ch).push_back(fBaseline);
      }
    }
    for(unsigned int i = 0; i < wvf_full.size(); i++) {
       wvf_channel.at(i_ch)[i] += (wvf_full[i] - fBaseline);
    }

    hist_id++;

    wvfmHandle.clear();
  } // waveform handle loop

  // access hardware trigger information
  std::vector<size_t> triggerIndex;
  for(auto const& trigger : (*triggerHandle)) {
    for (size_t idx = 0; idx < trigger.numPassed.size(); idx++) {
      if (trigger.numPassed[idx] >= fMultiplicityThreshold) {
        // save index
        triggerIndex.push_back(idx);
        // skip ahead by 5120 samples (downsampled by 4) before checking again
        idx += 1280;
      }
    }
  } // trigger handle loop

  if (fVerbose) std::cout << "Number of PMT hardware triggers found: " << triggerIndex.size() << std::endl;

  // set properties of fragment that are common to event
  uint32_t nChannelsTotal = channelList.size();
  //uint32_t nChannelsFrag = 15;
  uint32_t nFrag = (uint32_t)nChannelsTotal/nChannelsFrag;
  //uint32_t wfm_length = 5120; // ~10us, 2ns tick

  // fragment properties
  uint32_t sequenceIDVal = fEvent;

  // create and populate common metadata
  sbndaq::CAENV1730FragmentMetadata metadata;
  metadata.nChannels = nChannelsFrag + 1; // 15 PMT channels + final channel to store beam window / trigger information
  metadata.nSamples = wfm_length;

  // fragment handle properties
  uint32_t eventCounterVal = fEvent;
  uint32_t boardIDVal = 0;
  uint32_t triggerTimeTagVal = (uint32_t)CLHEP::RandFlat::shoot(&fTriggerTimeEngine, 0, 1e9);
  uint32_t eventSizeVal = ((wfm_length * (nChannelsFrag+1)) * sizeof(uint16_t) + sizeof(sbndaq::CAENV1730EventHeader)) / sizeof(uint32_t);

  // loop over PMT hardware triggers
  for (auto wvfIdx : triggerIndex) {

    // index in full waveform, 2ns tick
    size_t trigIdx = wvfIdx*4;
    // size_t startIdx = trigIdx-500; // -1us
    size_t startIdx = abs(fMinStartTime)*1000/2 + trigIdx-500;

    // determine and set timestamp for particular trigger
    // double triggerTime = fMinStartTime + wvfIdx*0.008; // in us
    double triggerTime = wvfIdx*0.008; // in us
    double timestampVal = 0.5 + (triggerTime*1e-6); // in seconds // std::time(nullptr); // current time
    metadata.timeStampSec = (uint32_t)timestampVal;
    metadata.timeStampNSec = (uint32_t)(timestampVal*1e9);

    // create fragments to hold waveforms, set properties and populate
    // 15 PMTs stored per fragment, 120/15 = 8 fragments per trigger
    for (size_t counter = 0; counter < nFrag; counter++) {

      // create fragment
      const auto fragment_datasize_bytes = metadata.ExpectedDataSize();
      uint32_t fragmentIDVal = counter;
      auto fragment_uptr = artdaq::Fragment::FragmentBytes(fragment_datasize_bytes, sequenceIDVal, fragmentIDVal, sbndaq::detail::FragmentType::CAENV1730, metadata); // unique pointer
      fragment_uptr->setTimestamp(timestampVal);

      // populate fragment header
      auto header_ptr = reinterpret_cast<sbndaq::CAENV1730EventHeader*>(fragment_uptr->dataBeginBytes());

      header_ptr->eventCounter = eventCounterVal;
      header_ptr->boardID = boardIDVal;
      header_ptr->triggerTimeTag = triggerTimeTagVal;  // ns // set timetag as random value for event
      header_ptr->eventSize = eventSizeVal;

      // populate waveforms
      // populate fragment with waveform
      uint16_t* data_begin = reinterpret_cast<uint16_t*>(fragment_uptr->dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
      uint16_t* value_ptr = data_begin;
      uint16_t value = 0;
      size_t ch_offset = 0;

      // loop over channels
      for (size_t i_ch = 0; i_ch < nChannelsFrag; i_ch++) {
        ch_offset = (size_t)(i_ch*wfm_length);
        // loop over waveform
        for (size_t i_t = 0; i_t < wfm_length; i_t++) {
          // set value
          value = wvf_channel[counter*nChannelsFrag + i_ch][startIdx+i_t];
          value_ptr = data_begin + ch_offset + i_t;
          *value_ptr = value;
        }
      }

      // create add beam window trigger waveform
      ch_offset = (size_t)(nChannelsFrag*wfm_length);
      size_t beamStartIdx = abs(fMinStartTime)*1000/2;
      size_t beamEndIdx = beamStartIdx + fBeamWindowLength*1000/2;
      // loop over waveform
      for (size_t i_t = 0; i_t < wfm_length; i_t++) {
        // set value
        if (startIdx + i_t >= beamStartIdx && startIdx + i_t <= beamEndIdx) value = 1;
        else value = 0;
        value_ptr = data_begin + ch_offset + i_t;
        *value_ptr = value;
      }

      // add fragment to vector
      vecFrag->push_back(*fragment_uptr);
    }
  }

  if(fVerbose) std::cout << "PMT Fragments written: " << vecFrag->size() - num_crt_frags << std::endl;


    //produce fragment vector
    e.put(std::move(vecFrag));

    // clear variables
    wvf_channel.clear();
    wvf_channel.shrink_to_fit();

} // ArtdaqFragmentProducer::produce()



DEFINE_ART_MODULE(sbnd::trigger::ArtdaqFragmentProducer)
