/**
 *  Service for the PMT timing corrections.
 *  Andrea Scarpelli (ascarpell@bnl.gov), Matteo Vicenzi (mvicenzi@bnl.gov)
 */
// Ported to SBND by Alejandro Sanchez-Castillo, Jan. 2025

#ifndef SBNDCODE_DATABASEINTERFACE_PMTCALIBRATIONDATABASEPROVIDER_H
#define SBNDCODE_DATABASEINTERFACE_PMTCALIBRATIONDATABASEPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"

// Local
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"

// C/C++ standard libraries
#include <string>
#include <map>
#include <stdint.h>


namespace sbndDB{ class PMTCalibrationDatabaseProvider; }
/**
 * @brief 
 * 
 * This module reads the PMT timing corrections from the database.
 *
 */

//--------------------------------------------------------------------------------

// PMT Scalings for systematic variations. All scale fields optional
struct PMTScaleConfig {
    fhicl::Atom<int> channel {
        fhicl::Name("Channel"),
        fhicl::Comment("Which channel to apply the scaling to (-1 for all channels)")
    };

    fhicl::OptionalAtom<float> spe_response {
        fhicl::Name("SPEAmplitude"),
        fhicl::Comment("Multiplicative scaling applied to SPE response")
    };
};


class sbndDB::PMTCalibrationDatabaseProvider : public PMTCalibrationDatabase {

    public:

        PMTCalibrationDatabaseProvider(const fhicl::ParameterSet& pset);

        void readPMTCalibrationDatabase(const art::Run& run);

        int getBreakoutBox( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).breakoutBox;
        };
        int getCAENDigitizer( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).caenDigitizer;
        };
        int getCAENDigitizerChannel( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).caenDigitizerChannel;
        };
        bool getOnPMT(unsigned int channelID) const override {
            return getChannelCorrOrDefault(channelID).onPMT;
        };
        bool getReconstructChannel(unsigned int channelID) const override {
            return getChannelCorrOrDefault(channelID).reconstructChannel;
        };
        double getTotalTransitTime( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).totalTransitTime;
        };
        double getCosmicTimeCorrection(unsigned int channelID) const override {
            return getChannelCorrOrDefault(channelID).cosmicTimeCorrection;
        };
        double getSPEAmplitude( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).spe_amplitude;
        };
        double getSPEAmplitudeStd( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).spe_amplitude_std;
        };
        double getGaussFilterPower( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).gauss_wc_power;
        };
        double getGaussFilterWC( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).gauss_wc;
        };
        double getNonLineatiryPESat( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).nonlinearity_pesat;
        };
        double getNonLineatiryAlpha( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).nonlinearity_alpha;
        };
        std::vector<double> getSER( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).ser;
        };
    private:
        
        bool fVerbose = false; ///< Whether to print the configuration we read.
        bool fApplyScales = false; // Apply scalings to DB parameters for systematic variations
        std::string fLogCategory; ///< Category tag for messages.
        std::string fPMTCalibrationDatabaseTag;  
        long fDatabaseTimeStamp;
        std::string fTableName;
        size_t fSERLength;

        /// Structure for single channel corrections
        struct PMTCalibrationDB { 

            size_t breakoutBox=0;
            size_t caenDigitizer=0;
            size_t caenDigitizerChannel=0;
            bool onPMT=true;
            bool reconstructChannel=false;
            double totalTransitTime=0.;
            double cosmicTimeCorrection=0.;
            double spe_amplitude=0.;
            double spe_amplitude_std=0.;
            double gauss_wc_power=0.;
            double gauss_wc=0.;
            double nonlinearity_pesat=0.;
            double nonlinearity_alpha=0.;
            std::vector<double> ser={};
            };

        /// Structure for scaling factors to DB values (for syst. studies)
        struct PMTCalibrationDBScales {
            // bools are OR'd with nominal value
            std::optional<bool> onPMT;
            std::optional<bool> reconstructChannel;

            // otherwise multiplicative
            std::optional<double> totalTransitTime;
            std::optional<double> cosmicTimeCorrection;
            std::optional<double> spe_amplitude;
            std::optional<double> spe_amplitude_std;
            std::optional<double> gauss_wc_power;
            std::optional<double> gauss_wc;
            std::optional<double> nonlinearity_pesat;
            std::optional<double> nonlinearity_alpha;

            // SER not implemented
        };
            
        const PMTCalibrationDB CorrectionDefaults = {0, 0, 0, true, false, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}};

	    /// Map of corrections by channel
        std::map<unsigned int, PMTCalibrationDB> fPMTCalibrationData;

        /// per-channel scales
        std::map<unsigned int, PMTCalibrationDBScales> fPMTCalibrationScales;
        
        /// Internal access to the channel correction record; returns defaults if not present.
        PMTCalibrationDB const& getChannelCorrOrDefault
            (unsigned int channelID) const
            {
                auto const it = fPMTCalibrationData.find(channelID);
                return (it == fPMTCalibrationData.end())? CorrectionDefaults: it->second;
            }

        /// Convert run number to internal database
        uint64_t RunToDatabaseTimestamp(uint32_t run) const;

        void ReadPMTCalibration(uint32_t run);
        PMTCalibrationDB scalePMTCalibrationData(const PMTCalibrationDB&, const PMTCalibrationDBScales&) const;
        void modifyScale(PMTCalibrationDBScales&, const PMTCalibrationDBScales&) const;
        std::string calibDataStr(const PMTCalibrationDB&) const;


}; // services class

#endif
