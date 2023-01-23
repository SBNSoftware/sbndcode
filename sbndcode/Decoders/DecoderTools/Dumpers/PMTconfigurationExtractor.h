/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTconfigurationExtractor.h
 * @brief  Utility to extract PMT readout configuration from data.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 23, 2023
 * @file   sbndcode/Decoders/DecoderTools/Dumper/PMTconfigurationExtractor.cxx
 * 
 * The design of this thing is still a big mess, winking to a templated
 * class hierarchy to be used for other components as well, but not getting
 * quite there.
 * 
 * Extra linking is required for `extractPMTreadoutConfiguration(Principal)`
 * to pick up the class used as `Principal`.
 * 
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTCONFIGURATIONEXTRACTOR_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTCONFIGURATIONEXTRACTOR_H


// SBND libraries
#include "sbndcode/Decoders/ChannelMapping/ISBNDChannelMap.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/ReadArtConfiguration.h" // util::readConfigurationFromArtPrincipal()
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"
#include "sbnobj/Common/PMT/Data/V1730Configuration.h"

// framework libraries
#include "art/Framework/Principal/DataViewImpl.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TFile.h"

// C/C++ standard libraries
#include <regex>
#include <string>
#include <optional>
#include <utility> // std::move(), std::pair<>
#include <initializer_list>



// -----------------------------------------------------------------------------
namespace sbnd {
  
  class PMTconfigurationExtractorBase;
  class PMTconfigurationExtractor;
  
  sbn::PMTconfiguration extractPMTreadoutConfiguration
    (std::string const& srcFileName, sbnd::PMTconfigurationExtractor extractor);
  sbn::PMTconfiguration extractPMTreadoutConfiguration
    (TFile& srcFile, sbnd::PMTconfigurationExtractor extractor);
  template <typename Principal>
  sbn::PMTconfiguration extractPMTreadoutConfiguration
    (Principal const& data, sbnd::PMTconfigurationExtractor extractor);
  
} // namespace sbnd


// -----------------------------------------------------------------------------
class sbnd::PMTconfigurationExtractorBase {
  
    public:
  
  using ConfigurationData_t = sbn::PMTconfiguration;
  
  // --- BEGIN -- Interface ----------------------------------------------------
  /// @name Interface
  /// @{
  
  /// Returns whether `pset` may contain the needed configuration.
  static bool mayHaveConfiguration(fhicl::ParameterSet const& pset)
    { return pset.has_key("configuration_documents"); }
  
  
  /**
   * @brief Extracts all supported PMT configuration from `config`.
   * @param config a FHiCL parameter set with component configuration
   * @return an object with the supported PMT configuration
   * 
   * All PMT-related configuration that is known to this code is extracted and
   * returned.
   * 
   * This function is undefined here: it must be overridden.
   */
  ConfigurationData_t extract(fhicl::ParameterSet const& config) const;
  
  
  /// Finalizes the content of `config` and returns it.
  ConfigurationData_t finalize(ConfigurationData_t config) const
    { return config; }
  
  
  /// @}
  // --- END ---- Interface ----------------------------------------------------
  
  
  // --- BEGIN -- Utility ------------------------------------------------------
  /// @name Utility
  /// @{
  
  /**
   * @brief Returns a parameter set with the content of
   *        `configuration_documents` key from `container`.
   * @param container parameter set including a table with key `configListKey`
   * @param configListKey name of the key in `container` with the configuration
   * @param components keys to be converted (as regular expressions)
   * @return a parameter set
   * 
   * The `configuration_documents` element of `container` is processed: for
   * each of its keys which match at least one of the `components` regular
   * expression patterns (`std::regex_match()`), the associated string value
   * is parsed with FHiCL parser, and the result is set as a FHiCL table in
   * the output parameter set.
   * For example, if the `components` are
   * `{ std::regex{".*pmt.*"}, std::regex{".*trigger.*"} }`, the returned
   * value is a parameter set that may have keys like `sbndpmtee01`,
   * `sbndpmtew02`, `sbndtrigger` etc., each one with a FHiCL table as
   * `value.
   */
  static fhicl::ParameterSet convertConfigurationDocuments(
    fhicl::ParameterSet const& container,
    std::string const& configListKey,
    std::initializer_list<std::regex const> components
    );
  
  /// Returns whether `key` matches at least one of the regular expressions
  /// in the [ `rbegin`, `rend` [ range.
  template <typename RBegin, typename REnd>
  static bool matchKey(std::string const& key, RBegin rbegin, REnd rend);
  
  /// @}
  // --- END ---- Utility ------------------------------------------------------
  
}; // sbnd::PMTconfigurationExtractorBase


// -----------------------------------------------------------------------------
/**
 * @brief Class to extract PMT readout board configuration.
 * 
 * This is an example of PMT readout board configuration taken from SBND run
 * 4774 (`destinations` and `metrics` have been omitted):
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * sbndpmteetop01: {
 *   daq: {
 *     fragment_receiver: {
 *       BaselineCh1: 14768
 *       BaselineCh10: 14843
 *       BaselineCh11: 14846
 *       BaselineCh12: 14840
 *       BaselineCh13: 14981
 *       BaselineCh14: 14896
 *       BaselineCh15: 14966
 *       BaselineCh16: 14800
 *       BaselineCh2: 14836
 *       BaselineCh3: 14806
 *       BaselineCh4: 14954
 *       BaselineCh5: 14766
 *       BaselineCh6: 14835
 *       BaselineCh7: 14861
 *       BaselineCh8: 14857
 *       BaselineCh9: 14949
 *       BoardChainNumber: 0
 *       CalibrateOnConfig: true
 *       ChargePedstalBitCh1: 14400
 *       CircularBufferSize: 5e8
 *       CombineReadoutWindows: false
 *       GetNextFragmentBunchSize: 20
 *       GetNextSleep: 100
 *       InterruptEventNumber: 1
 *       InterruptLevel: 0
 *       InterruptStatusID: 0
 *       LVDSLogicValueG1: 3
 *       LVDSLogicValueG2: 3
 *       LVDSLogicValueG3: 3
 *       LVDSLogicValueG4: 3
 *       LVDSLogicValueG5: 3
 *       LVDSLogicValueG6: 3
 *       LVDSLogicValueG7: 3
 *       LVDSLogicValueG8: 3
 *       LVDSOutWidthC1: 20
 *       LVDSOutWidthC10: 20
 *       LVDSOutWidthC11: 20
 *       LVDSOutWidthC12: 20
 *       LVDSOutWidthC13: 20
 *       LVDSOutWidthC14: 20
 *       LVDSOutWidthC15: 20
 *       LVDSOutWidthC16: 20
 *       LVDSOutWidthC2: 20
 *       LVDSOutWidthC3: 20
 *       LVDSOutWidthC4: 20
 *       LVDSOutWidthC5: 20
 *       LVDSOutWidthC6: 20
 *       LVDSOutWidthC7: 20
 *       LVDSOutWidthC8: 20
 *       LVDSOutWidthC9: 20
 *       LockTempCalibration: false
 *       MajorityCoincidenceWindow: 1
 *       MajorityLevel: 0
 *       MajorityTimeWindow: 4
 *       ModeLVDS: 1
 *       SWTrigger: false
 *       SelfTrigBit: 16
 *       SelfTriggerMask: 255
 *       SelfTriggerMode: 0
 *       TimeOffsetNanoSec: 0
 *       UseTimeTagForTimeStamp: false
 *       Verbosity: 1
 *       acqMode: 0
 *       allowTriggerOverlap: false
 *       analogMode: 1
 *       boardId: 0
 *       board_id: 18
 *       channelEnable0: true
 *       channelEnable1: true
 *       channelEnable10: true
 *       channelEnable11: true
 *       channelEnable12: true
 *       channelEnable13: true
 *       channelEnable14: true
 *       channelEnable15: true
 *       channelEnable2: true
 *       channelEnable3: true
 *       channelEnable4: true
 *       channelEnable5: true
 *       channelEnable6: true
 *       channelEnable7: true
 *       channelEnable8: true
 *       channelEnable9: true
 *       channelPedestal0: 6554
 *       channelPedestal1: 6554
 *       channelPedestal10: 6554
 *       channelPedestal11: 6554
 *       channelPedestal12: 6554
 *       channelPedestal13: 6554
 *       channelPedestal14: 6554
 *       channelPedestal15: 6554
 *       channelPedestal2: 6554
 *       channelPedestal3: 6554
 *       channelPedestal4: 6554
 *       channelPedestal5: 6554
 *       channelPedestal6: 6554
 *       channelPedestal7: 6554
 *       channelPedestal8: 6554
 *       channelPedestal9: 6554
 *       dacValue: 32768
 *       data_buffer_depth_fragments: 10000
 *       data_buffer_depth_mb: 2000
 *       debugLevel: 7
 *       destinations: {} # ...
 *       dynamicRange: 0
 *       enableReadout: 1
 *       eventCounterWarning: 1
 *       eventsPerInterrupt: 1
 *       extTrgMode: 3
 *       fragment_id: 18
 *       fragment_type: "CAENV1730"
 *       generator: "CAENV1730Readout"
 *       ioLevel: 1
 *       irqWaitTime: 1
 *       link: 2
 *       maxEventsPerTransfer: 1
 *       max_fragment_size_bytes: 1e7
 *       memoryAlmostFull: 2
 *       multicast_interface_ip: "192.168.191.0"
 *       nChannels: 16
 *       outputSignalMode: 0
 *       postPercent: 70
 *       readoutMode: 0
 *       receive_requests: true
 *       recordLength: 25000
 *       request_address: "227.128.1.129"
 *       request_mode: "sequence"
 *       request_port: 3502
 *       routing_table_config: { use_routing_master: false }
 *       runSyncMode: 0
 *       separate_data_thread: true
 *       swTrgMode: 0
 *       testPattern: 0
 *       triggerPolarity: 0
 *       triggerPulseWidth: 20
 *       triggerThreshold0: 14368
 *       triggerThreshold1: 14436
 *       triggerThreshold10: 14446
 *       triggerThreshold11: 14440
 *       triggerThreshold12: 14581
 *       triggerThreshold13: 14596
 *       triggerThreshold14: 14566
 *       triggerThreshold15: 14400
 *       triggerThreshold2: 14406
 *       triggerThreshold3: 14554
 *       triggerThreshold4: 14366
 *       triggerThreshold5: 14435
 *       triggerThreshold6: 14461
 *       triggerThreshold7: 14457
 *       triggerThreshold8: 14549
 *       triggerThreshold9: 14443
 *       usePedestals: false
 *     }
 *     metrics: {} # ...
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * 
 * 
 */
class sbnd::PMTconfigurationExtractor
  : public sbnd::PMTconfigurationExtractorBase
{
  
    public:
  
  /// Constructor: no channel mapping, channel IDs won't be `finalize()`'d.
  PMTconfigurationExtractor() = default;
  
  /// Constructor: use channel mapping to `finalize()` channel IDs.
  PMTconfigurationExtractor
    (sbndDB::ISBNDChannelMap const& channelMappingService);
  
  
  // --- BEGIN -- Interface ----------------------------------------------------
  /// @name Interface
  /// @{
  
  /// Returns whether `pset` may contain the needed configuration.
  static bool mayHaveConfiguration(fhicl::ParameterSet const& pset)
    { return pset.has_key("configuration_documents"); }
  
  /// Returns whether the specified `key` of `pset` is a good configuration.
  static bool isGoodConfiguration
    (fhicl::ParameterSet const& pset, std::string const& key);
  
  
  /**
   * @brief Extracts all supported PMT configuration from `config`.
   * @param config a FHiCL parameter set with component configuration
   * @return an object with the supported PMT configuration
   * 
   * All PMT-related configuration that is known to this code is extracted and
   * returned.
   */
  sbn::PMTconfiguration extract(fhicl::ParameterSet const& config) const;
  
  
  /// Assigns unique channel IDs to the channel information.
  sbn::PMTconfiguration finalize(sbn::PMTconfiguration config) const;
  
  /// @}
  // --- END ---- Interface ----------------------------------------------------
  
    private:
  
  /// Regular expressions matching all names of supported PMT configurations.
  static std::vector<std::regex> const ConfigurationNames;
  
  /// Hardware PMT channel mapping to LArSoft's.
  sbndDB::ISBNDChannelMap const* fChannelMap = nullptr;
  
  /**
   * @brief Extracts PMT readout board configuration from `pset`.
   * @param pset information source (FHiCL configuration)
   * @param boardName name of the board we are given the configuration of
   * @return the requested configuration
   * 
   * 
   */
  sbn::V1730Configuration extractV1730configuration
    (fhicl::ParameterSet const& pset, std::string const& boardName)
    const;
  
  /**
   * @brief Returns the specified V1730 readout channel configuration.
   * @param boardPSet readout board configuration
   * @param channelNo number of channel on board (0-15)
   * @return the configuration of the specified channel
   */
  sbn::V1730channelConfiguration extractChannelConfiguration
    (fhicl::ParameterSet const& boardPSet, unsigned short int channelNo) const;
  
  /**
   * @brief Returns the specified PMT readout board configuration.
   * @param pset parameter set including `key`
   * @param key key of the PMT readout configuration candidate
   * @return the configuration, or an empty object if key does not represent one
   */
  std::optional<fhicl::ParameterSet> readBoardConfig
    (fhicl::ParameterSet const& pset, std::string const& key) const;
  
  
  /// Returns the fragment ID of the specified board as known by the database.
  static unsigned int readoutBoardDBfragmentID
    (sbn::V1730Configuration const& boardConfig);
  
}; // sbnd::PMTconfigurationExtractor


// -----------------------------------------------------------------------------
// --- Inline implementation
// -----------------------------------------------------------------------------
sbnd::PMTconfigurationExtractor::PMTconfigurationExtractor
  (sbndDB::ISBNDChannelMap const& channelMappingService)
  : fChannelMap(&channelMappingService)
  {}


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
namespace sbnd::details {
  
  template <typename ConfigMap>
  sbn::PMTconfiguration extractPMTreadoutConfigurationImpl
    (ConfigMap const& configMap, sbnd::PMTconfigurationExtractor extractor)
  {
    
    /*
    * Requirements: ConfigMap is a mapping type supporting iteration and whose
    * elements support structured assignment into a pair, the first element
    * being an identifier (process name, parameter set ID...) and the second
    * being a `fhicl::ParameterSet` with the configuration (or something
    * offering the same interface).
    * 
    * The plan is to look in all the FHiCL configuration fragments we can find
    * in the input config, and find all the useful configuration therein.
    * Given that there may be multiple input files, there may also be multiple
    * configurations for the same detector components.
    * In that case, we will extract parameters from each and every one of the
    * configurations, and throw an exception if they are not all consistent.
    * 
    * Consistency is tested only for the extracted parameters, not for the whole
    * FHiCL configuration fragment.
    */
    
    using Key_t = std::tuple_element_t<0U, typename ConfigMap::value_type>;
    
    std::optional<std::pair<Key_t, sbn::PMTconfiguration>> config;
    
    // look in the global configuration for all parameter sets which contain
    // `configuration_documents` as a (direct) name;
    for (auto const& [ id, pset ]: configMap) {
      if (!extractor.mayHaveConfiguration(pset)) continue;
      
      fhicl::ParameterSet const configDocs
        = extractor.convertConfigurationDocuments
          (pset, "configuration_documents", { std::regex{ "sbndpmt.*" } })
        ;
      
      sbn::PMTconfiguration candidateConfig = extractor.extract(configDocs);
      if (config) {
        if (config->second == candidateConfig) continue;
        mf::LogError log("extractPMTreadoutConfiguration");
        log << "Found two candidate configurations differring:"
          "\nFirst [" << config->first << "]:\n" << config->second
          << "\nSecond [" << id << "]:\n" << candidateConfig
          ;
        throw cet::exception("extractPMTreadoutConfiguration")
          << "extractPMTreadoutConfiguration() found inconsistent configurations.\n";
      } // if incompatible configurations
      
      config.emplace(std::move(id), std::move(candidateConfig));
    } // for all configuration documents
    
    if (!config) {
      throw cet::exception("extractPMTreadoutConfiguration")
        << "extractPMTreadoutConfiguration() could not find a suitable configuration.\n";
    }
    
    return extractor.finalize(std::move(config->second));
  } // extractPMTreadoutConfigurationImpl(ConfigMap)
  
  
} // namespace sbnd::details


// -----------------------------------------------------------------------------
template <typename RBegin, typename REnd>
bool sbnd::PMTconfigurationExtractorBase::matchKey
  (std::string const& key, RBegin rbegin, REnd rend)
{
  for (auto iRegex = rbegin; iRegex != rend; ++iRegex)
    if (std::regex_match(key, *iRegex)) return true;
  return false;
} // sbnd::PMTconfigurationExtractorBase::matchKey()


// -----------------------------------------------------------------------------
template <typename Principal>
sbn::PMTconfiguration sbnd::extractPMTreadoutConfiguration
  (Principal const& principal, sbnd::PMTconfigurationExtractor extractor)
{
  
  return details::extractPMTreadoutConfigurationImpl
    (util::readConfigurationFromArtPrincipal(principal), std::move(extractor));
  
} // sbnd::extractPMTreadoutConfiguration(Principal)


// ---------------------------------------------------------------------------


#endif // SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTCONFIGURATIONEXTRACTOR_H
