////////////////////////////////////////////////////////////////////////
// Class:       CRTArtdaqFragmentProducer
// Module Type: producer
// File:        CRTArtdaqFragmentProducer_module.cc
//
// Producer module for creating CRT DAQ-readable fragment from simulation
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


// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// artdaq includes
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
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

// ROOT includes.
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TString.h"

// C++ includes
#include <bitset>
#include <map>
#include <vector>
#include <string>
#include <random>
#include <iomanip>

namespace sbnd {
  namespace trigger {
    class CRTArtdaqFragmentProducer;
  }
}

class sbnd::trigger::CRTArtdaqFragmentProducer : public art::EDProducer {
public:
  explicit CRTArtdaqFragmentProducer(fhicl::ParameterSet const& p);

  CRTArtdaqFragmentProducer(CRTArtdaqFragmentProducer const&) = delete;
  CRTArtdaqFragmentProducer(CRTArtdaqFragmentProducer&&) = delete;
  CRTArtdaqFragmentProducer & operator=(CRTArtdaqFragmentProducer const&) = delete;
  CRTArtdaqFragmentProducer & operator=(CRTArtdaqFragmentProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;


private:

  int fRun;
  int fSubrun;
  int fEvent;

  // fcl file parameters
  art::InputTag fSimModuleLabel;      ///< name of detsim producer
  art::InputTag fCRTSimLabel;         ///< name of CRT producer
  std::string   fFEBDataLabel;        ///< name of FEBData producer
  bool          fVerbose;             ///< print information about what's going on
  double        fClockSpeedCRT;
  size_t           fFirstFEBMac5;        ///< lowest mac5 address for CRT FEBs

  // Other variables shared between different methods.
  geo::GeometryCore const* fGeometryService;
  sbnd::crt::CRTGeoAlg fCrtGeo;

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




}; // class CRTArtdaqFragmentProducer


sbnd::trigger::CRTArtdaqFragmentProducer::CRTArtdaqFragmentProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fSimModuleLabel(p.get<std::string>("SimModuleLabel", "largeant")),
  fCRTSimLabel(p.get<std::string>("CRTSimLabel", "crt")),
  fFEBDataLabel(p.get<std::string>("FEBDataLabel", "crtsim")),
  fVerbose(p.get<bool>("Verbose", false)),
  fClockSpeedCRT(p.get<double>("ClockSpeedCRT")),
  fFirstFEBMac5(p.get<size_t>("FirstFEBMac5", 0))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<artdaq::Fragment> >();

  // Get a pointer to the fGeometryServiceetry service provider
  fGeometryService = lar::providerFrom<geo::Geometry>();

}


void sbnd::trigger::CRTArtdaqFragmentProducer::produce(art::Event& e)
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

    //random number generator to be used to simulate pedestal
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(175.0,23.0);

    //list all fragment ids - uncomment if needed for daq running
  /*if (fVerbose){std::cout<<"Num Modules: "<<fCrtGeo.NumModules()<<std::endl;}
    for (size_t mod_i = 0; mod_i<fCrtGeo.NumModules(); mod_i++){
      int plane = sbnd::CRTCommonUtils::GetPlaneIndex(fCrtGeo.GetTaggerName(fCrtGeo.ChannelToStripName(mod_i * 32)));
      if (fVerbose){std::cout << std::hex << (32768 + 12288 + (plane * 256) + (int)mod_i) << std::endl;}
    } */

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

          if (ts1==0){flags=11;}else if (ts0==0){flags=7;}

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

    if(fVerbose) std::cout << "Fragments written: " << vecFrag->size() << std::endl;

    e.put(std::move(vecFrag));


} // CRTArtdaqFragmentProducer::produce()



DEFINE_ART_MODULE(sbnd::trigger::CRTArtdaqFragmentProducer)
