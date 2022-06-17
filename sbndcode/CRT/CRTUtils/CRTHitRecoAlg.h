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
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

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

      fhicl::Atom<double> NpeScaleShift {
        Name("NpeScaleShift"),
        Comment("Parameter for correcting for distance along strip")
      };

      fhicl::Atom<double> TimeCoincidenceLimit {
        Name("TimeCoincidenceLimit"),
        Comment("Minimum time between two overlapping hit crt strips")
      };

      fhicl::Atom<double> ClockSpeedCRT {
        Name("ClockSpeedCRT"),
        Comment("Clock speed of the CRT system [MHz]")
      };

      fhicl::Atom<double> TimeOffset {
        Name("TimeOffset"),
        Comment("The global time offset to use for the CRT times")
      };
      fhicl::Atom<bool> UseG4RefTimeOffset {
        Name("UseG4RefTimeOffset"),
        Comment("Wheater or not to use the G4RefTime as GlobalT0Offset"),
      };
      fhicl::Atom<double> PropDelay {
        Name("PropDelay"),
        Comment("Delay in pulse arrival time [ns/m]"),
      };
      fhicl::Atom<double> TDelayNorm {
        Name("TDelayNorm"),
        Comment("Time delay fit: Gaussian normalization"),
      };
      fhicl::Atom<double> TDelayShift {
        Name("TDelayShift"),
        Comment("Time delay fit: Gaussian x shift"),
      };
      fhicl::Atom<double> TDelaySigma {
        Name("TDelaySigma"),
        Comment("Time delay fit: Gaussian width"),
      };
      fhicl::Atom<double> TDelayOffset {
        Name("TDelayOffset"),
        Comment("Time delay fit: Gaussian baseline offset"),
      };


    };

    CRTHitRecoAlg(const Config& config, double g4RefTime);

    CRTHitRecoAlg(const fhicl::ParameterSet& pset,  double g4RefTime) :
      CRTHitRecoAlg(fhicl::Table<Config>(pset, {})(), g4RefTime) {}

    CRTHitRecoAlg();

    ~CRTHitRecoAlg();

    void reconfigure(const Config& config);

    std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> CreateTaggerStrips(detinfo::DetectorClocksData const& clockData,
                                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                                         std::vector<art::Ptr<sbnd::crt::CRTData>> data);

    CRTStrip CreateCRTStrip(art::Ptr<sbnd::crt::CRTData> sipm1, art::Ptr<sbnd::crt::CRTData> sipm2, size_t ind);

    std::pair<double, double> DistanceBetweenSipms(art::Ptr<sbnd::crt::CRTData> sipm1, art::Ptr<sbnd::crt::CRTData> sipm2);
    
    std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> CreateCRTHits(std::map<std::pair<std::string, unsigned>, std::vector<CRTStrip>> taggerStrips);

    // Function to calculate the strip position limits in real space from channel
    std::vector<double> ChannelToLimits(CRTStrip strip);
 
    // Function to calculate the overlap between two crt strips
    std::vector<double> CrtOverlap(std::vector<double> strip1, std::vector<double> strip2);
 
    // Function to return the CRT tagger name and module position from the channel ID
    std::pair<std::string,unsigned> ChannelToTagger(uint32_t channel);
 
    // Function to check if a CRT strip overlaps with a perpendicular module
    bool CheckModuleOverlap(uint32_t channel);
 
    // Function to make filling a CRTHit a bit faster
    sbn::crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                           std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, int plane, 
                           double x, double ex, double y, double ey, double z, double ez, std::string tagger); 

    // Function to correct number of photoelectrons by distance down strip
    double CorrectNpe(CRTStrip strip1, CRTStrip strip2, TVector3 position);

    // Function to correct the time by distance down strip
    double CorrectTime(CRTStrip strip1, CRTStrip strip2, TVector3 position);

  private:

    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

    bool fUseReadoutWindow;
    double fQPed;
    double fQSlope;
    double fNpeScaleShift;
    double fTimeCoincidenceLimit;
    double fClockSpeedCRT;
    double fTimeOffset;
    bool fUseG4RefTimeOffset;
    double fPropDelay;
    double fTDelayNorm;
    double fTDelayShift;
    double fTDelaySigma;
    double fTDelayOffset;

  };

}

#endif
