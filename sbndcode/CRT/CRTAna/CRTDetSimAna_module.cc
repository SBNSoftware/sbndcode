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
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
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

      fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
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
    // ADCs of crt data
    std::map<std::string, TH1D*> hADC;
    // Time of CRT data
    std::map<std::string, TH1D*> hTime;
    std::map<std::string, TH1D*> hStripDist;
    std::map<std::string, TH1D*> hSipmDist;

    // Efficiency plots
    std::map<std::string, TH1D*> hEffSipmDistTotal;
    std::map<std::string, TH1D*> hEffSipmDistReco;
    std::map<std::string, TH1D*> hEffStripDistTotal;
    std::map<std::string, TH1D*> hEffStripDistReco;
    std::map<std::string, TH1D*> hEffEDepTotal;
    std::map<std::string, TH1D*> hEffEDepReco;
    std::map<std::string, TH1D*> hEffLengthTotal;
    std::map<std::string, TH1D*> hEffLengthReco;

    // ADC as a function of distance from distance along width
    std::map<std::string, TH2D*> hStripDistADC;
    std::map<std::string, TH2D*> hSipmDistADC;

    // Other variables shared between different methods.
    detinfo::DetectorClocks const* fDetectorClocks;
    detinfo::ElecClock fTrigClock;
    geo::GeometryCore const* fGeometryService;
    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;
    CRTBackTracker fCrtBackTrack;

  }; // class CRTDetSimAna

  // Constructor
  CRTDetSimAna::CRTDetSimAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTSimLabel          (config().CRTSimLabel())
    , fVerbose              (config().Verbose())
    , fCrtBackTrack         (config().CrtBackTrack())
  {
    // Get a pointer to the fGeometryServiceetry service provider
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    fGeometryService = lar::providerFrom<geo::Geometry>();
  }

  void CRTDetSimAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
      std::string tagger = fCrtGeo.GetTagger(i).name;
      hADC[tagger]          = tfs->make<TH1D>(Form("ADC_%s", tagger.c_str()),       "", 40, 0, 10000);
      hTime[tagger]         = tfs->make<TH1D>(Form("Time_%s", tagger.c_str()),      "", 40, -5000, 5000);
      hStripDist[tagger]    = tfs->make<TH1D>(Form("StripDist_%s", tagger.c_str()), "", 40, 0,     450);
      hSipmDist[tagger]     = tfs->make<TH1D>(Form("SipmDist_%s", tagger.c_str()),  "", 40, 0,     11.2);

      hEffSipmDistTotal[tagger]  = tfs->make<TH1D>(Form("EffSipmDistTotal_%s", tagger.c_str()),  "", 20, 0, 11.2);
      hEffSipmDistReco[tagger]   = tfs->make<TH1D>(Form("EffSipmDistReco_%s", tagger.c_str()),   "", 20, 0, 11.2);
      hEffStripDistTotal[tagger] = tfs->make<TH1D>(Form("EffStripDistTotal_%s", tagger.c_str()), "", 20, 0, 450);
      hEffStripDistReco[tagger]  = tfs->make<TH1D>(Form("EffStripDistReco_%s", tagger.c_str()),  "", 20, 0, 450);
      hEffEDepTotal[tagger]      = tfs->make<TH1D>(Form("EffEdepTotal_%s", tagger.c_str()),      "", 20, 0, 0.01);
      hEffEDepReco[tagger]       = tfs->make<TH1D>(Form("EffEdepReco_%s", tagger.c_str()),       "", 20, 0, 0.01);
      hEffLengthTotal[tagger]    = tfs->make<TH1D>(Form("EffLengthTotal_%s", tagger.c_str()),    "", 20, 0, 10);
      hEffLengthReco[tagger]     = tfs->make<TH1D>(Form("EffLengthReco_%s", tagger.c_str()),     "", 20, 0, 10);

      hStripDistADC[tagger] = tfs->make<TH2D>(Form("StripDistADC_%s", tagger.c_str()), "", 20, 0, 450,  20, 0, 10000);
      hSipmDistADC[tagger]  = tfs->make<TH2D>(Form("SipmDistADC_%s", tagger.c_str()),  "", 20, 0, 11.2, 20, 0, 10000);
    }
    

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

    // Get all the auxdet IDEs from the event
    art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
    event.getByLabel(fSimModuleLabel, channels);

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
    //                                        CRT DATA ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    std::map<int, std::vector<crt::CRTData>> crtData;
    for(size_t i  = 0; i < crtDataList.size(); i++){

      // Truth matching
      std::vector<int> trueIDs = fCrtBackTrack.AllTrueIds(event, *crtDataList[i]);
      for(auto const& id : trueIDs){
        crtData[id].push_back(*crtDataList[i]);
      }

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
      std::string tagger = fCrtGeo.GetTaggerName(stripName);

      // Calculate distance to the sipm
      double sipmDist = fCrtGeo.DistanceBetweenSipms(cross, channel);
      hSipmDist[tagger]->Fill(sipmDist);
      hSipmDistADC[tagger]->Fill(sipmDist, crtDataList[i]->ADC());

      // Calculate the distance to the end
      double stripDist = std::abs(fCrtGeo.DistanceDownStrip(cross, stripName));
      hStripDist[tagger]->Fill(stripDist);
      hStripDistADC[tagger]->Fill(stripDist, crtDataList[i]->ADC());

      hADC[tagger]->Fill(crtDataList[i]->ADC());

      fTrigClock.SetTime(crtDataList[i]->T0());
      double time = fTrigClock.Time(); // [us]
      hTime[tagger]->Fill(time);
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          EFFICIENCIES
    //----------------------------------------------------------------------------------------------------------

    for (auto& adsc : *channels) {
      // Get the IDEs from the Aux Det channels
      const geo::AuxDetGeo& adGeo = fGeometryService->AuxDet(adsc.AuxDetID());
      const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());
      std::string stripName = adsGeo.TotalVolume()->GetName();
      stripName.append("_0");
      std::string tagger = fCrtGeo.GetTaggerName(stripName);
      std::vector<sim::AuxDetIDE> ides = adsc.AuxDetIDEs();

      for(auto const& ide : ides){
        geo::Point_t pos {(ide.entryX + ide.exitX)/2, 
                          (ide.entryY + ide.exitY)/2, 
                          (ide.entryZ + ide.exitZ)/2};
        int trueID = ide.trackID;
        // Calculate length in strip
        double length = std::sqrt(std::pow(ide.entryX - ide.exitX, 2)
                                  + std::pow(ide.entryY - ide.exitY, 2)
                                  + std::pow(ide.entryZ - ide.exitZ, 2));
        // Calculate distances of average point to sipms
        double sipmDistIde = fCrtGeo.DistanceBetweenSipms(pos, stripName);
        double stripDistIde = fCrtGeo.DistanceDownStrip(pos, stripName);

        // Only consider primary muons for efficiency
        if(particles.find(trueID) == particles.end()) continue;
        if(!(std::abs(particles[trueID].PdgCode()) == 13 && particles[trueID].Mother() == 0)) continue;

        // Fill total histograms
        hEffSipmDistTotal[tagger]->Fill(sipmDistIde);
        hEffStripDistTotal[tagger]->Fill(stripDistIde);
        hEffEDepTotal[tagger]->Fill(ide.energyDeposited);
        hEffLengthTotal[tagger]->Fill(length);

        // Get all the CRT data matched to the same true ID
        if(crtData.find(trueID) == crtData.end()) continue;
        bool match = false;
        // If there is CRT data produced on the same channel then it is matched
        for(auto const& data : crtData[trueID]){
          std::string stripNameData = fCrtGeo.ChannelToStripName(data.Channel());
          if(stripName == stripNameData) match = true;
        }
        if(!match) continue;

        // Fill reconstructed histograms
        hEffSipmDistReco[tagger]->Fill(sipmDistIde);
        hEffStripDistReco[tagger]->Fill(stripDistIde);
        hEffEDepReco[tagger]->Fill(ide.energyDeposited);
        hEffLengthReco[tagger]->Fill(length);
          
      }
    }


  } // CRTDetSimAna::analyze()

  void CRTDetSimAna::endJob(){

  } // CRTDetSimAna::endJob()
  
  
  DEFINE_ART_MODULE(CRTDetSimAna)
} // namespace sbnd


