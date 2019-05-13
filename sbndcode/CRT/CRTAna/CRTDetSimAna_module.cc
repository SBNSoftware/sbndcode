////////////////////////////////////////////////////////////////////////
// Class:       CRTDetSimAna
// Module Type: analyzer
// File:        CRTDetSimAna_module.cc
//
// Analysis module for evaluating CRT reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TString.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace sbnd {

  class CRTDetSimAna : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> CRTSimLabel {
        Name("CRTSimLabel"),
        Comment("tag of CRT simulation data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTDetSimAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTSimLabel;         ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on
    
    // Histograms
    // ADC as a function of distance from distance along width
    TH2D *hEndDistanceADC;
    TH2D *hSipmDistanceADC;
    // ADC as a function of distance from end
    // ADCs of crt data
    TH1D *hADC;
    // Time of CRT data
    TH1D *hTime;
    TH1D *hChannel;

    // Other variables shared between different methods.
    detinfo::DetectorClocks const* fDetectorClocks;
    detinfo::ElecClock fTrigClock;
    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

  }; // class CRTDetSimAna

  // Constructor
  CRTDetSimAna::CRTDetSimAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTSimLabel          (config().CRTSimLabel())
    , fVerbose              (config().Verbose())
  {
    // Get a pointer to the fGeometryServiceetry service provider
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fTrigClock = fDetectorClocks->TriggerClock();
  }

  void CRTDetSimAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    hADC          = tfs->make<TH1D>("ADC","",100, 0, 20000);
    hTime         = tfs->make<TH1D>("Time","",100, -5000, 5000);
    hChannel      = tfs->make<TH1D>("Channel", "", 400, 0, 10000);

    hEndDistanceADC = tfs->make<TH2D>("EndDistanceADC", "", 20, 0, 500, 20, 0, 10000);
    hSipmDistanceADC = tfs->make<TH2D>("SipmDistanceADC", "", 20, 0, 12, 20, 0, 10000);

    // Initial output
    std::cout<<"----------------- Full CRT Detector Simulation Analysis Module -------------------"<<std::endl;

  }// CRTDetSimAna::beginJob()

  void CRTDetSimAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }


    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get g4 particles
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get CRT data from the event
    art::Handle< std::vector<crt::CRTData>> crtDataHandle;
    std::vector<art::Ptr<crt::CRTData> > crtDataList;
    if (event.getByLabel(fCRTSimLabel, crtDataHandle))
      art::fill_ptr_vector(crtDataList, crtDataHandle);

    art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTSimLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
     
    std::map<int, simb::MCParticle> particles;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;

    }

    //----------------------------------------------------------------------------------------------------------
    //                                            ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    for(size_t i  = 0; i < crtDataList.size(); i++){
      // Get the IDEs associated with the crtData
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
      // Calculate the average position
      double x = 0, y = 0, z = 0;
      for(auto const& ide : ides){
        x += (ide->entryX + ide->exitX)/2.;
        y += (ide->entryY + ide->exitY)/2.;
        z += (ide->entryZ + ide->exitZ)/2.;
      }
      x = x/ides.size();
      y = y/ides.size();
      z = z/ides.size();
      geo::Point_t cross {x, y, z};
      int channel = crtDataList[i]->Channel();
      std::string stripName = fCrtGeo.ChannelToStripName(channel);
      // Calculate the distance to the end
      double distBetweenSipms = std::abs(fCrtGeo.DistanceBetweenSipms(cross, channel));
      hSipmDistanceADC->Fill(distBetweenSipms, crtDataList[i]->ADC());
      double distToEnd = std::abs(fCrtGeo.DistanceDownStrip(cross, stripName));
      hEndDistanceADC->Fill(distToEnd, crtDataList[i]->ADC());
      fTrigClock.SetTime(crtDataList[i]->T0());
      double time = fTrigClock.Time(); // [us]
      hADC->Fill(crtDataList[i]->ADC());
      hTime->Fill(time);
      hChannel->Fill(crtDataList[i]->Channel());
    }


  } // CRTDetSimAna::analyze()

  void CRTDetSimAna::endJob(){

  } // CRTDetSimAna::endJob()
  
  
  DEFINE_ART_MODULE(CRTDetSimAna)
} // namespace sbnd


