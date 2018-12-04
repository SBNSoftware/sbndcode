////////////////////////////////////////////////////////////////////////
// Class:       CRTTest
// Module Type: analyzer
// File:        CRTTestAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTruthMatchUtils.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"


// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace sbnd {


  class CRTTest : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> CrtTrackModuleLabel {
        Name("CrtTrackModuleLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> CrtHitModuleLabel {
        Name("CrtHitModuleLabel"),
        Comment("tag of CRT track producer data product")
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTTest(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

  private:

    // fcl file parameters
    art::InputTag fCrtTrackModuleLabel; ///< name of CRT track producer
    art::InputTag fCrtHitModuleLabel; ///< name of CRT track producer

  }; // class CRTTest


  // Constructor
  CRTTest::CRTTest(Parameters const& config)
    : EDAnalyzer(config)
    , fCrtTrackModuleLabel  (config().CrtTrackModuleLabel())
    , fCrtHitModuleLabel    (config().CrtHitModuleLabel())
  {

  } // CRTTest()


  void CRTTest::beginJob()
  {
  } // CRTTest::beginJob()


  void CRTTest::analyze(const art::Event& event)
  {

    // Get CRT tracks
    art::Handle< std::vector<crt::CRTTrack> > crtTrackListHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCrtTrackModuleLabel, crtTrackListHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackListHandle);

    // Get track to hit associations
    art::FindManyP<crt::CRTHit> findManyHits(crtTrackListHandle, event, fCrtTrackModuleLabel);

    for(size_t i = 0; i < 4; i++){ 
      std::vector<art::Ptr<crt::CRTHit>> hits = findManyHits.at(i);
      std::cout<<"Hits size = "<<hits.size()<<"\n";
      art::FindManyP<crt::CRTData> findManyData(hits, event, "crthit");

      std::vector<int> ids = CRTTruthMatchUtils::AllTrueIds(crtTrackListHandle, event, fCrtTrackModuleLabel, fCrtHitModuleLabel, (int)i);
      for(size_t l = 0; l < ids.size(); l++){
        std::cout<<"ID = "<<ids[l]<<"\n";
      }
      for(size_t j = 0; j < hits.size(); j++){
        std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(j);
        std::cout<<"Data size = "<<data.size()<<"\n";
        art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, "crt");
        for(size_t k = 0; k < data.size(); k++){
          std::vector<art::Ptr<sim::AuxDetIDE>> ide = findManyIdes.at(k);
          std::cout<<"ide size = "<<ide.size()<<"\n";
          for(size_t m = 0; m < ide.size(); m++){
            std::cout<<ide[m]->trackID<<"\n";
          }
        }
      }
    }


    
  } // CRTTest::analyze()


  void CRTTest::endJob(){

  
  } // CRTTest::endJob()


  DEFINE_ART_MODULE(CRTTest)
} // namespace sbnd

