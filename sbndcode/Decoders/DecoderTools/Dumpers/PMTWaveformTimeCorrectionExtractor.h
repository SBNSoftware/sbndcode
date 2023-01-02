/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTWaveformTimeCorrectionExtractor.h
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H

// SBND/SBN libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/SBNDChannelMap.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrections.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTWaveformTimeCorrection.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace sbnd::timing { class PMTWaveformTimeCorrectionExtractor; }

/**
 * @brief Extracts timing corrections from a waveform.
 * 
 * This algorithm extracts a time correction from a reference signal, and
 * associates that correction with all the channels in the same PMT readout
 * crate.
 * 
 * The correction extraction is performed by `findWaveformTimeCorrections()`.
 * That method analyzes a reference waveform searching for a reference signal
 * assumed to have been generated in time with the global trigger, and emits
 * a correction so that this signal would appear at exactly the time of the
 * global trigger (`detinfo::DetectorClocks::TriggerTime()`).
 * 
 * 
 * Signal timing extraction
 * -------------------------
 * 
 * The detection algorithm is currently quite unsophisticated.
 * 
 * The reference signal is expected to be a sharp square wave in negative
 * polarity. The time of the correction is based on the front side of that wave:
 * 
 *  * the absolute minimum of the waveform is found
 *  * an interval starting 20 ticks before that minimum is considered
 *  * the baseline level is defined as the value at the start of that interval
 *  * the start time is set to the exact tick with an amplitude exceeding 20%
 *    of the maximum of the signal from the baseline
 * 
 */
class sbnd::timing::PMTWaveformTimeCorrectionExtractor {
	
	public: 

        // --- BEGIN -- Exceptions ---------------------------------------------

        /// Exception thrown when trying to overwrite a correction.
        struct Error: cet::exception { Error(std::string const& msg = ""); };

        /// Exception thrown when trying to overwrite a correction.
        struct MultipleCorrectionsForChannel: Error {
            MultipleCorrectionsForChannel(unsigned int existing, unsigned int additional);
        };

        /// Exception thrown when correction requested from a non-special channel.
        struct NotASpecialChannel: Error {
            NotASpecialChannel(unsigned int channel);
            private: static Error makeBaseException(unsigned int channel);
        };

        /// Exception thrown when PMT readout crate not recognised.
        struct UnknownCrate: Error {
            UnknownCrate(unsigned int channel);
            private: static Error makeBaseException(unsigned int channel);
        };

        // --- END ---- Exceptions ---------------------------------------------


        PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            sbndDB::SBNDChannelMap const & channelMapService,
            sbndDB::PMTTimingCorrections const* pmtTimingCorrectionsService, 
            bool verbose );

        /**
         * @brief Extracts a correction from `wave` and assigns it to channels.
         * @param wave the reference waveform to extract the correction from
         * @param correctCableDelay whether to apply the correction for cable delays
         * @param[in,out] corrections where to add the newly extracted corrections
         * @throw NotASpecialChannel if `wave` is not a special channel
         * @throw UnknownCrate if `wave` is from an unexpected readout crate
         * @throw MultipleCorrectionsForChannel if `corrections` already contains
         *   a correction for any of the channels we are associating to the new correction
         * 
         * This function performs the analysis of the reference waveform `wave`,
         * extracts the correction out of it and associates it in `corrections`
         * for all the channels that are in the same readout crate as the
         * reference waveform itself.
         * 
         * The `corrections` list is extended as needed. Since `corrections`
         * is a dense structure with as index the channel ID, it may happen that
         * the extension produces correction information for channels not in
         * this readout crate: in this case, these correction values are invalid
         * (as in `PMTTimingCorrections::isValid()` returning `false`).
         * Any attempt to overwrite a valid correction already in `corrections`
         * will cause throwing an exception.
         * 
         * The channel ID in each stored correction is the one of the reference
         * waveform the correction was extracted from, and likewise the `sample`
         * data member is where the reference signal can be found within that
         * waveform. The `startTime` correction is an offset to be _added_ to
         * the waveform timestamps to correct them.
         * 
         */
        void findWaveformTimeCorrections(   
            raw::OpDetWaveform const & wave,
            bool correctCableDelay,
            std::vector<PMTWaveformTimeCorrection> & corrections ) const;


	private:

        detinfo::DetectorClocksData const fClocksData;

        sbndDB::SBNDChannelMap const & fChannelMap;

        sbndDB::PMTTimingCorrections const* fPMTTimingCorrectionsService = nullptr;

        bool const fVerbose;

        std::map<unsigned int, std::vector<unsigned int>> const
            fCrateFragmentMap {
                {0x0070 , { 0 , 1 , 2  }}, 
                {0x0060 , { 3 , 4 , 5  }}, 
                {0x0050 , { 6 , 7 , 8  }}, 
                {0x0040 , { 9 , 10, 11 }},
                {0x0030 , { 12, 13, 14 }}, 
                {0x0020 , { 15, 16, 17 }}, 
                {0x0010 , { 18, 19, 20 }}, 
                {0x0000 , { 21, 22, 23 }}, 
            }; 

        template<typename T>
            static size_t getMaxBin( 
                std::vector<T> const& vv,
                size_t startElement, 
                size_t endElement);

        template<typename T>
            static size_t getMinBin( 
                std::vector<T> const& vv,
                size_t startElement, 
                size_t endElement);

        template<typename T>
            static size_t getStartSample( std::vector<T> const& vv );

};


#endif //SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H
