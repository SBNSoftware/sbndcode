/**
 * \brief Class for SBND CRT detector simulation.
 *
 * \details This class performs the SBND CRT detector simulation starting from AuxDetSimChannels.
 * Each AuxDetSimChannel can be passed to this class by calling method FillTaggers for each
 * AuxDetSimChannel. Then, the CreateData method performs the CRT detector simulation
 * (time and charge responde, and triggering).
 *
 * \author Andy Mastbaum
 * \author Marco Del Tutto
 */

#ifndef SBND_CRTDETSIMALG_H
#define SBND_CRTDETSIMALG_H

// art includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

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
#include <chrono>

// ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "Math/Interpolator.h"

// CRT includes
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "CRTDetSimParams.h"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

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
        struct Trigger;
    }
}

/** A struct to temporarily store information on a single SiPM.
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

/** A struct to temporarily store information on a single CRT Strip.
 */
struct sbnd::crt::StripData {
    StripData(uint16_t _mac5, uint16_t _flags, unsigned _orientation, SiPMData _sipm0, SiPMData _sipm1,
              uint32_t _unixs, bool _sipm_coinc, sim::AuxDetIDE _ide) :
        mac5(_mac5)
      , flags(_flags)
      , orientation(_orientation)
      , sipm0(_sipm0)
      , sipm1(_sipm1)
      , unixs(_unixs)
      , sipm_coinc(_sipm_coinc)
      , ide(_ide)
    {};

    uint16_t mac5; ///< The FEB number this strip is in
    uint16_t flags; ///< The flags given by the data of this strip
    unsigned orientation; ///< The orientation of the module (0 or 1 for horizontal or vertical)
    SiPMData sipm0; ///< One SiPM (the one that sees signal first)
    SiPMData sipm1; ///< One SiPM
    uint32_t unixs; ///< The unix time of the readout
    bool sipm_coinc; ///< Stores the ID of the SiPM that sees signal first (0-31)
    sim::AuxDetIDE ide; ///< The AuxDetIDE associated with this strip
};

/** A struct to temporarily store information on a CRT Tagger.
 */
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
    std::vector<std::pair<FEBData, std::vector<AuxDetIDE>>> GetData();

    /**
     * Returns the indeces of SiPMs associated to the AuxDetIDEs
     *
     * @return Vector of vector (1: FEBs, 2: SiPMs indeces per AuxDetIDE)
     */
    std::vector<std::vector<int>> GetAuxData();


    /**
     * Get the channel trigger time relative to the start of the MC event.
     *
     * @param clock The clock to count ticks on
     * @param t0 The starting time (which delay is added to)
     * @param npe Number of observed photoelectrons
     * @param r Distance between the energy deposit and strip readout end [mm]
     * @return Trigger clock ticks at this true hit time
     */
    uint32_t getChannelTriggerTicks(/*detinfo::ElecClock& clock,*/
                                    float t0, float npeMean, float r);


    /**
     * Simulated the CRT charge response. Does not include waveform emulation.
     *
     * @param eDep The simulated energy deposited on the strip from G4
     * @param d0 The distance to the optical fiber connected to SiPM 0
     * @param d1 The distance to the optical fiber connected to SiPM 1
     * @param distToReadout The distance to the readout
     * @param npe0 The output number of PEs for SiPM 0
     * @param npe1 The output number of PEs for SiPM 1
     * @param q0 The output ADC simulated value (double) for SiPM 0
     * @param q1 The output ADC simulated value (double) for SiPM 1
     */
    void ChargeResponse(double eDep, double d0, double d1, double distToReadout,
                        long & npe0, long & npe1, double & q0, double & q1);


    /**
     * Emulated the CRT slow-shaped waveform. This is not a full waveform simulation,
     * but it tries to encapsulate the main effect of the waveform, which is that the ADC
     * counts are reduced if the signal arrives with a time delay w.r.t. the primary event
     * trigger.
     *
     * @param time_delay The time delay of this SiPM signal w.r.t. the primary event trigger.
     * @param adc The simulated ADC counts (double) of this SiPM.
     * @return The simulated ADC counts after waveform emulation (in uint16_t format).
     */
    uint16_t WaveformEmulation(const uint32_t & time_delay, const double & adc);


    /**
     * Returns the CRT simulation parmeters
     *
     * @return The CRT simulation parameters
     */
    const CRTDetSimParams & Params() {return fParams;}



private:

    CRTDetSimParams fParams; //!< The table of CRT simulation parameters

    CLHEP::HepRandomEngine& fEngine; //!< The random-number engine

    double fG4RefTime; //!< The G4 reference time that can be used as a time offset
    double fTimeOffset; //!< The time that will be used in the simulation

    std::unique_ptr<ROOT::Math::Interpolator> fInterpolator; //!< The interpolator used to estimate the CRT waveform

    std::map<std::string, Tagger> fTaggers; //!< A list of hit taggers, before any coincidence requirement (name -> tagger)

    std::vector<std::pair<FEBData, std::vector<AuxDetIDE>>> fData; //!< This member stores the final FEBData for the CRT simulation

    std::vector<std::vector<int>> fAuxData; //!< This member stores the indeces of SiPM per AuxDetIDE

    CRTGeoAlg fCRTGeoAlg;

    /**
     * Configures the waveform by reading waveform points from configuration and
     * setting up the interpolator.
     */
    void ConfigureWaveform();

    /**
     * Configures the time offset to use, either a custom number,
     * of the G4RefTime, depending on user configuration.
     */
    void ConfigureTimeOffset();

    /**
     * Proccesses a set of CRT strips that belong to the same trigger. This method
     * takes as input all the strips that belong to a single CRT tagger-level trigger
     * and constructs FEBData objects from them.
     *
     * @param strips The set of strips that belong to the same trigger
     */
    void ProcessStrips(const std::vector<StripData> & strips);

    /**
     * Adds ADCs to a certain SiPM in a FEBData object
     *
     * @param feb_data The FEBData object.
     * @param sipmID The SiPM index (0-31).
     * @param adc ADC value to be added.
     */
    void AddADC(FEBData & feb_data, const int & sipmID, const uint16_t & adc);

};

#endif
