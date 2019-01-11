#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN


///////////////////////////////////////////////
// CRTHitRecoAlg.h
//
// Functions for CRT hit reconstruction
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"


namespace sbnd{

  struct CRTStrip {
    double t0; // [us]
    uint32_t channel; 
    double x; // [cm]
    double ex; // [cm]
    double pes;
    std::pair<std::string, unsigned> tagger;
    size_t dataID;
  };


  class CRTHitRecoAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<bool> UseReadoutWindow {
        Name("UseReadoutWindow"),
        Comment("Only reconstruct hits within readout window")
      };

      fhicl::Atom<double> QPed {
        Name("QPed"),
        Comment("Pedestal offset [ADC]")
      };

      fhicl::Atom<double> QSlope {
        Name("QSlope"),
        Comment("Pedestal slope [ADC/photon]")
      };

      fhicl::Atom<double> TimeCoincidenceLimit {
        Name("TimeCoincidenceLimit"),
        Comment("Minimum time between two overlapping hit crt strips")
      };

    };

    CRTHitRecoAlg(const Config& config);

    CRTHitRecoAlg(const fhicl::ParameterSet& pset) :
      CRTHitRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTHitRecoAlg();

    ~CRTHitRecoAlg();

    void reconfigure(const Config& config);

    std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> CreateTaggerStrips(std::vector<art::Ptr<crt::CRTData>> data);

    CRTStrip CreateCRTStrip(art::Ptr<crt::CRTData> sipm1, art::Ptr<crt::CRTData> sipm2, size_t ind);
    
    std::vector<std::pair<crt::CRTHit, std::vector<int>>> CreateCRTHits(std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> taggerStrips);

    // Function to calculate the strip position limits in real space from channel
    std::vector<double> ChannelToLimits(CRTStrip strip);
 
    // Function to calculate the overlap between two crt strips
    std::vector<double> CrtOverlap(std::vector<double> strip1, std::vector<double> strip2);
 
    // Function to return the CRT tagger name and module position from the channel ID
    std::pair<std::string,unsigned> ChannelToTagger(uint32_t channel);
 
    // Function to check if a CRT strip overlaps with a perpendicular module
    bool CheckModuleOverlap(uint32_t channel);
 
    // Function to make filling a CRTHit a bit faster
    sbnd::crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                           std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, int plane, 
                           double x, double ex, double y, double ey, double z, double ez, std::string tagger); 

  private:

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorClocks const* fDetectorClocks;
    detinfo::DetectorProperties const* fDetectorProperties;
    detinfo::ElecClock fTrigClock;
    art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    const geo::AuxDetGeometry* fAuxDetGeo;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

    bool fUseReadoutWindow;
    double fQPed;
    double fQSlope;
    double fTimeCoincidenceLimit;

  };

}

#endif
