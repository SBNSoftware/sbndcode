/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTconfigurationExtractor.cxx
 * @brief  Utility to extract PMT readout configuration from data.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * @see    sbndcode/Decoders/DecoderTools/Dumpers/PMTconfigurationExtractor.h
 */


// SBND libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTconfigurationExtractor.h"

// C/C++ standard libraries
#include <algorithm> // std::find_if()
#include <memory> // std::unique_ptr<>



// -----------------------------------------------------------------------------
// ---  sbnd::PMTconfigurationExtractorBase
// -----------------------------------------------------------------------------
fhicl::ParameterSet
sbnd::PMTconfigurationExtractorBase::convertConfigurationDocuments(
  fhicl::ParameterSet const& container,
  std::string const& configListKey,
  std::initializer_list<std::regex const> components // TODO template with matcher functor
) {
  
  fhicl::ParameterSet const sourceConfig
    = container.get<fhicl::ParameterSet>(configListKey);
  
  fhicl::ParameterSet configDocs;
  for (auto const& key: sourceConfig.get_names()) {
    if (!sourceConfig.is_key_to_atom(key)) continue;
    
    if (!matchKey(key, components.begin(), components.end())) continue;
    
    std::string const psetStr = sourceConfig.get<std::string>(key);
    
    fhicl::ParameterSet pset;
    try {
      pset = fhicl::ParameterSet::make(psetStr);
    }
    catch (cet::exception& e) {
      throw cet::exception{ "convertConfigurationDocuments", "", e }
        << "Error parsing the content of key '" << configListKey << "." << key
        << "'; content was:\n" << psetStr << "\n";
    }
    
    configDocs.put(key, pset);
    
  } // for all main keys
  
  return configDocs;
} // convertConfigurationDocuments()



// -----------------------------------------------------------------------------
// ---  sbnd::PMTconfigurationExtractor
// -----------------------------------------------------------------------------
std::vector<std::regex> const
  sbnd::PMTconfigurationExtractor::ConfigurationNames
  { std::regex{ "sbndpmt.*" } }
  ;

// -----------------------------------------------------------------------------
bool sbnd::PMTconfigurationExtractor::isGoodConfiguration
  (fhicl::ParameterSet const& pset, std::string const& key)
{
  return matchKey(key, ConfigurationNames.begin(), ConfigurationNames.end());
} // sbnd::PMTconfigurationExtractor::isGoodConfiguration()


// -----------------------------------------------------------------------------
sbn::PMTconfiguration sbnd::PMTconfigurationExtractor::extract
  (fhicl::ParameterSet const& pset) const
{
  
  sbn::PMTconfiguration config;
  
  for (std::string const& key: pset.get_names()) {
    
    std::optional<fhicl::ParameterSet> boardConfig
      = readBoardConfig(pset, key);
    
    if (!boardConfig) continue;
    
    config.boards.push_back(extractV1730configuration(*boardConfig, key));
    
  } // for
  
  return config;
} // sbnd::PMTconfigurationExtractor::extract()


// -----------------------------------------------------------------------------
auto sbnd::PMTconfigurationExtractor::extractV1730configuration
  (fhicl::ParameterSet const& pset, std::string const& boardName) const
  -> sbn::V1730Configuration
{
  
  auto const& boardParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver");
  
  sbn::V1730Configuration rc; // readout config, for friends
  rc.boardName = boardName;
  
  rc.boardID = boardParams.get<unsigned int>("board_id");
  rc.fragmentID = boardParams.get<unsigned int>("fragment_id");
  
  rc.bufferLength = boardParams.get<int>("recordLength");
  rc.postTriggerFrac = boardParams.get<float>("postPercent") / 100.0f;
  
  rc.nChannels = boardParams.get<unsigned int>("nChannels");
  
  rc.useTimeTagForTimeStamp = boardParams.get("UseTimeTagForTimeStamp", false);
  
  for (unsigned short int channelNo = 0; channelNo < 16; ++channelNo)
    rc.channels.push_back(extractChannelConfiguration(boardParams, channelNo));
  
  return rc;
} // sbnd::PMTconfigurationExtractor::extractV1730configuration()


// -----------------------------------------------------------------------------
sbn::V1730channelConfiguration
sbnd::PMTconfigurationExtractor::extractChannelConfiguration
  (fhicl::ParameterSet const& boardPSet, unsigned short int channelNo) const
{
  std::string const ChannelStr = std::to_string(channelNo);
  
  sbn::V1730channelConfiguration channel;
  
  channel.channelNo = channelNo;
  
  channel.baseline = boardPSet.get<short signed int>
    ("BaselineCh" + std::to_string(channelNo + 1));
  channel.threshold = boardPSet.get<short signed int>
    ("triggerThreshold" + ChannelStr);
  channel.enabled = boardPSet.get<bool>("channelEnable" + ChannelStr);
  
  return channel;
} // sbnd::PMTconfigurationExtractor::extractChannelConfiguration()


// ---------------------------------------------------------------------------
sbn::PMTconfiguration sbnd::PMTconfigurationExtractor::finalize
  (sbn::PMTconfiguration config) const
{

  for (sbn::V1730Configuration& readoutBoardConfig: config.boards) {
    auto const fragmentID = readoutBoardDBfragmentID(readoutBoardConfig);
    if (!fChannelMap->hasPMTDigitizerID(fragmentID)) {
      mf::LogWarning("PMTconfigurationExtractor")
        << "No entry found in PMT channel mapping database for board '"
        << readoutBoardConfig.boardName << "' (fragment ID: "
        << readoutBoardConfig.fragmentID << " => " << std::hex << fragmentID
        << ")\n";
      continue;
    }
  
    sbndDB::DigitizerChannelChannelIDPairVec const& digitizerChannelVec
      = fChannelMap->getChannelIDPairVec(fragmentID);
    
    // finds the channel ID matching the specified channel number of this board
    auto const toChannelID = [&channelIDs=digitizerChannelVec]
      (short unsigned int channelNo)
      {
        auto const it = std::find_if(channelIDs.begin(), channelIDs.end(),
          [channelNo](auto const& p){ return std::get<0U>(p) == channelNo; });
        return (it != channelIDs.end())
          ? std::get<1U>(*it)
          : sbn::V1730channelConfiguration::NoChannelID
          ;
      };
    
    for (auto& channelInfo: readoutBoardConfig.channels)
      channelInfo.channelID = toChannelID(channelInfo.channelNo);
    
  } // for boards
  
  return config;
} // sbnd::PMTconfigurationExtractor::finalize()


// -----------------------------------------------------------------------------
std::optional<fhicl::ParameterSet>
sbnd::PMTconfigurationExtractor::readBoardConfig
  (fhicl::ParameterSet const& pset, std::string const& key) const
{
  static std::string const ExpectedFragmentType = "CAENV1730";
  
  std::optional<fhicl::ParameterSet> config;
  
  do { // fake loop for fast exit
    
    // it must be a parameter set
    if (!pset.has_key(key) || !pset.is_key_to_table(key)) break;
    
    auto boardPSet = pset.get<fhicl::ParameterSet>(key);
    
    // its "fragment_type" must be the expected one
    std::string fragmentType;
    if (!boardPSet.get_if_present("daq.fragment_receiver.fragment_type", fragmentType))
      break;
    if (fragmentType != ExpectedFragmentType) break;
    
    config.emplace(std::move(boardPSet)); // success
  } while (false);
  return config;
} // sbnd::PMTconfigurationExtractor::readBoardConfig()


// -----------------------------------------------------------------------------
// ---  free functions implementation
// -----------------------------------------------------------------------------
sbn::PMTconfiguration sbnd::extractPMTreadoutConfiguration
  (std::string const& srcFileName, sbnd::PMTconfigurationExtractor extractor)
{
  //
  // TFile::Open() call is needed to support non-local URL
  // (e.g. XRootD URL are not supported by TFile constructor).
  //
  return extractPMTreadoutConfiguration(*(
    std::unique_ptr<TFile>{ TFile::Open(srcFileName.c_str(), "READ") }
    ),
    std::move(extractor)
    );
} // sbnd::extractPMTreadoutConfiguration(string)


// -----------------------------------------------------------------------------
sbn::PMTconfiguration sbnd::extractPMTreadoutConfiguration
  (TFile& srcFile, sbnd::PMTconfigurationExtractor extractor)
{
  
  return details::extractPMTreadoutConfigurationImpl
    (util::readConfigurationFromArtFile(srcFile), std::move(extractor));
  
} // sbnd::extractPMTreadoutConfiguration(TFile)


// ---------------------------------------------------------------------------
unsigned int sbnd::PMTconfigurationExtractor::readoutBoardDBfragmentID
  (sbn::V1730Configuration const& boardConfig)
{
  return boardConfig.fragmentID & 0xFF; // secret recipe
} // sbnd::PMTconfigurationExtractor::readoutBoardDBfragmentID()


// -----------------------------------------------------------------------------

