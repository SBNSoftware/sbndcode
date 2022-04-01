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
#include "messagefacility/MessageLogger/MessageLogger.h"

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
#include "Math/Interpolator.h"

// CRT includes
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "CRTDetSimParams.h"

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::string;
using sim::AuxDetIDE;

namespace sbnd {
    namespace crt {
        class CRTDetSimAlg;
        struct SiPMData;
        struct StripData;
        struct Tagger;
    }
}

/** An enum type. 
 *  The documentation block cannot be put after the enum! 
 */
struct sbnd::crt::SiPMData {
    SiPMData(int _sipmID, int _channel, uint64_t _t0, uint64_t _t1, double _adc) :
        sipmID(_sipmID)
      , channel(_channel)
      , t0(_t0)
      , t1(_t1)
      , adc(_adc)
    {};

    int sipmID; ///< The SiPM ID in the strip (0-31)
    int channel; ///< The SiPM global channel number, calculated as 32 x (module) + 2 x (strip) + j
    uint64_t t0; ///< The SiPM Ts0
    uint64_t t1; ///< The SiPM Ts1
    double adc; ///< The SiPM simulated ADC value in double format
};

struct sbnd::crt::StripData {
    StripData(uint16_t _mac5, unsigned _planeID, SiPMData _sipm0, SiPMData _sipm1,
              bool _sipm_coinc, sim::AuxDetIDE _ide) :
        mac5(_mac5)
      , planeID(_planeID)
      , sipm0(_sipm0)
      , sipm1(_sipm1)
      , sipm_coinc(_sipm_coinc)
      , ide(_ide)
    {};

    uint16_t mac5; ///< The FEB number this strip is in
    unsigned planeID; ///< The plane ID (0 or 1 for horizontal or vertical)
    SiPMData sipm0; ///< One SiPM (the one that sees signal first) 
    SiPMData sipm1; ///< One SiPM
    bool sipm_coinc; ///< Stores the ID of the SiPM that sees signal first (0-31)
    sim::AuxDetIDE ide; ///< The AuxDetIDE associated with this strip
};

struct sbnd::crt::Tagger {
  std::vector<StripData> data; ///< A vector of strips that have energy deposits
  /** \brief Returns the number of strips with energy deposits for this tagger */
  int size() {return data.size();}
};



class sbnd::crt::CRTDetSimAlg {

public:

    using Parameters = fhicl::Table<CRTDetSimParams>;
    CRTDetSimAlg(const Parameters & p, CLHEP::HepRandomEngine& fRandEngine, double g4RefTime);

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
    void FillTaggers(const uint32_t adid, const uint32_t adsid, std::vector<AuxDetIDE> ides);

    /**
     * Returns FEBData objects.
     *
     * This function is called after loop over AuxDetSimChannels where FillTaggers
     * was used to perform first detsim step. This function applies trigger logic,
     * deadtime, and close-in-time signal biasing effects. it produces the
     * "triggered data" products which make it into the event. Use "GetData"
     * to retrieve the result.
     *
     * @return Vector of pairs (FEBData, vector of AuxDetIDE)
     */
    void CreateData();

    /**
     * Returns the simulated FEBData and AuxDetIDEs
     *
     * @return Vector of pairs (FEBData, vector of AuxDetIDE)
     */
    std::vector<std::pair<sbnd::crt::FEBData, std::vector<AuxDetIDE>>> GetData();



private:

    CRTDetSimParams fParams; //!< The table of CRT simulation parameters

    CLHEP::HepRandomEngine& fEngine; //!< The random-number engine

    double fG4RefTime;

    double fTimeOffset;

    std::unique_ptr<ROOT::Math::Interpolator> fInterpolator; //!< The interpolator used to estimate the CRT waveform

    std::map<std::string, Tagger> fTaggers; // A list of hit taggers, before any coincidence requirement (name -> tagger)

    // std::vector<sbnd::crt::FEBData> fFEBDatas;
    std::vector<std::pair<sbnd::crt::FEBData, std::vector<AuxDetIDE>>> fData;

    void ConfigureWaveform();
    void ConfigureTimeOffset();

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

    void ProcessStrips(/*const uint32_t & trigger_ts1,
                       const uint32_t & trigger_ts0,*/
                       const uint32_t & coinc,
                       const std::vector<StripData> & strips,
                       const std::string & tagger_name);

    uint16_t WaveformEmulation(const uint32_t & time_delay, const double & adc);

    void AddADC(sbnd::crt::FEBData & feb_data, const int & sipmID, const uint16_t & adc);

};

#endif