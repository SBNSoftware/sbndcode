#ifndef SBND_CRTDETSIMALG_H
#define SBND_CRTDETSIMALG_H

// art includes
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// larsoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// CLHEP includes
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

//C++ includes
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>

// ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"

// CRT includes
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::string;
using sim::AuxDetIDE;

namespace sbnd {
    namespace crt {
        class CRTDetSimAlg;
    }
}

// struct Tagger {
//   std::vector<std::pair<unsigned, uint32_t> > planesHit;
//   std::vector<sbnd::crt::CRTData> data;
//   std::vector<std::vector<sim::AuxDetIDE>> ides;
// };

struct SiPMData {
    // uint16_t mac5;
    // unsigned planeID;
    int sipmID;
    int channel;
    uint64_t t0;
    uint64_t t1;
    uint16_t adc;
    SiPMData(int _sipmID, int _channel, uint64_t _t0, uint64_t _t1, uint16_t _adc) :
        sipmID(_sipmID)
      , channel(_channel)
      , t0(_t0)
      , t1(_t1)
      , adc(_adc)
    {};
    // sim::AuxDetIDE ide;
};

struct StripData {
    uint16_t mac5;
    unsigned planeID;
    SiPMData sipm0;
    SiPMData sipm1;
    // int channel;
    // uint64_t t0;
    // uint64_t t1;
    // uint16_t adc;
    sim::AuxDetIDE ide;

    StripData(uint16_t _mac5, unsigned _planeID, SiPMData _sipm0, SiPMData _sipm1, sim::AuxDetIDE _ide) :
        mac5(_mac5)
      , planeID(_planeID)
      , sipm0(_sipm0)
      , sipm1(_sipm1)
      , ide(_ide)
    {};
};

// struct ModuleData {
//     uint16_t mac5;
//     std::vector<uint16_t> adcs;
//     std::vector<uint64_t> times;
// };

struct Tagger {
  // std::vector<uint32_t> mac5;
  // std::vector<unsigned> planeID;
  // std::vector<uint32_t> time;
  std::vector<StripData> data;
  // std::vector<sim::AuxDetIDE> ide;
  int size() {return data.size();}
};



class sbnd::crt::CRTDetSimAlg {

public:

    CRTDetSimAlg(fhicl::ParameterSet const & p, CLHEP::HepRandomEngine& fRandEngine);
    void reconfigure(fhicl::ParameterSet const & p);


    /**
       * Function to clear member data at beginning of each art::event
       */
    void ClearTaggers();


    /**
       * Filles CRT taggers from AuxDetSimChannels.
       *
       * Intented to be called within loop over AuxDetChannels and provided the
       * AuxDetChannelID, AuxDetSensitiveChannelID, vector of AuxDetIDEs and
       * the number of ides from the end of the vector to include in the detector sim.
       * It handles deposited energy ti light output at SiPM (including attenuation)
       * and to PEs from the SiPM with associated time stamps.
       *
       * @param adid The AuxDetChannelID
       * @param adsid The AuxDetSensitiveChannelID
       * @param ides The vector of AuxDetIDE
       */
    void FillTaggers(const uint32_t adid, const uint32_t adsid, vector<AuxDetIDE> ides);

    /**
       * Returns FEBData objects.
       *
       * This function is called after loop over AuxDetSimChannels where FillTaggers
       * was used to perform first detsim step. This function applies trigger logic,
       * deadtime, and close-in-time signal biasing effects. it produces the
       * "triggered data" products which make it into the event.
       *
       * @return Vector of pairs (FEBData, vector of AuxDetIDE)
       */
    std::vector<std::pair<sbnd::crt::CRTData, std::vector<AuxDetIDE>>> CreateData();


private:

    double fGlobalT0Offset;  //!< Time delay fit: Gaussian normalization
    double fTDelayNorm;  //!< Time delay fit: Gaussian normalization
    double fTDelayShift;  //!< Time delay fit: Gaussian x shift
    double fTDelaySigma;  //!< Time delay fit: Gaussian width
    double fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
    double fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
    double fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
    double fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
    double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
    double fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
    double fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
    double fClockSpeedCRT; //!< Clock speed for the CRT system [MHz]
    double fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
    double fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
    double fQ0;  //!< Average energy deposited for mips, for charge scaling [GeV]
    double fQPed;  //!< ADC offset for the single-peak peak mean [ADC]
    double fQSlope;  //!< Slope in mean ADC / Npe [ADC]
    double fQRMS;  //!< ADC single-pe spectrum width [ADC]
    double fQThreshold;  //!< ADC charge threshold [ADC]
    double fTResInterpolator;  //!< Interpolator time resolution [ns]
    double fPropDelay;  //!< Delay in pulse arrival time [ns/m]
    double fPropDelayError;  //!< Delay in pulse arrival time, uncertainty [ns/m]
    double fStripCoincidenceWindow;  //!< Time window for two-fiber coincidence [ns]
    double fTaggerPlaneCoincidenceWindow;  //!< Time window for two-plane coincidence [ticks]
    double fAbsLenEff;  //!< Effective abs. length for transverse Npe scaling [cm]
    bool fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
    double fSipmTimeResponse; //!< Minimum time to resolve separate energy deposits [ns]
    uint32_t fAdcSaturation; //!< Saturation limit per SiPM in ADC counts
    double fDeadTime; //!< Saturation limit per SiPM in ADC counts

    CLHEP::HepRandomEngine& fEngine;

    // A list of hit taggers, before any coincidence requirement (mac5 -> tagger)
    std::map<std::string, Tagger> fTaggers;

    std::vector<sbnd::crt::FEBData> fFEBDatas;

    /**
       * Get the channel trigger time relative to the start of the MC event.
       *
       * @param engine The random number generator engine
       * @param clock The clock to count ticks on
       * @param t0 The starting time (which delay is added to)
       * @param npe Number of observed photoelectrons
       * @param r Distance between the energy deposit and strip readout end [mm]
       * @return Trigger clock ticks at this true hit time
       */
    uint32_t getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                    /*detinfo::ElecClock& clock,*/
                                    float t0, float npeMean, float r);

    void ProcessStrips(uint32_t trigger_time, uint32_t coinc, std::vector<StripData> strips);
    uint16_t WaveformEmulation(uint32_t time_delay, uint16_t adc);
    void AddADC(sbnd::crt::FEBData & feb_data, int sipmID, uint16_t adc);

};

#endif