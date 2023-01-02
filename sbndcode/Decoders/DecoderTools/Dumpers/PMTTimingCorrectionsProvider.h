/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrectionsProvider.h
 * @brief  Service for the PMT timing corrections.
 * @author Andrea Scarpelli (apapadopoulou@anl.gov)
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTIMINGCORRECTIONSPROVIDER_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTIMINGCORRECTIONSPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrections.h"

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace sbndDB::details {
    
  struct PMTTimeCorrectionsDB {

    double triggerCableDelay=0; ///< Expect nanoseconds.

    double resetCableDelay=0; ///< Expect nanoseconds.

    double laserCableDelay=0; ///< Expect nanoseconds.

    double cosmicsCorrections=0; ///< Expect nanoseconds.
    
  };
  
} // sbndDB::details

namespace sbndDB{ class PMTTimingCorrectionsProvider; }
/**
 * @brief 
 * 
 * This module reads 
 * 
 * Input
 * ------
 * 
 * 
 * Output
 * -------
 * 
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * 
 */
class sbndDB::PMTTimingCorrectionsProvider : public PMTTimingCorrections {

    public: 

        PMTTimingCorrectionsProvider(const fhicl::ParameterSet& pset);

        void readTimeCorrectionDatabase(const art::Run& run);

        double getTriggerCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).triggerCableDelay;
        };

        double getResetCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).resetCableDelay;
        };

        double getLaserCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).laserCableDelay;
        };

        double getCosmicsCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).cosmicsCorrections;
        };

    private:
        
        using PMTTimeCorrectionsDB = details::PMTTimeCorrectionsDB;

        std::string fUrl;

        unsigned int fTimeout;

        bool fCorrectCablesDelay;

        bool fVerbose = false; ///< Whether to print the configuration we read.
  
        std::string fLogCategory; ///< Category tag for messages.

        /// Interface to LArSoft configuration for detector timing.
        detinfo::DetectorClocksData const fClocksData;

        static constexpr PMTTimeCorrectionsDB CorrectionDefaults {};
        

        std::map<unsigned int, PMTTimeCorrectionsDB> fDatabaseTimingCorrections;
        
        
        /// Internal access to the channel correction record; returns defaults if not present.
        PMTTimeCorrectionsDB const& getChannelCorrOrDefault
            (unsigned int channelID) const
            {
                auto const it = fDatabaseTimingCorrections.find(channelID);
                return (it == fDatabaseTimingCorrections.end())? CorrectionDefaults: it->second;
            }

        int ConnectToDataset(const std::string& name, 
            uint32_t run, Dataset& dataset ) const;

        void ReadPMTCablesCorrections(uint32_t run);

        void ReadLaserCorrections(uint32_t run);

        void ReadCosmicsCorrections(uint32_t run);

}; // services class

#endif 
