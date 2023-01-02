/**
 * @file   sbndcode/Decode/DecoderTools/details/PMTDecoderUtils.h
 * @brief  Some helpers for PMT decoder tool.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_PMTDECODERUTILS_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_PMTDECODERUTILS_H


// SBND/SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "sbnobj/Common/PMT/Data/V1730channelConfiguration.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::lower_bound(), std::sort()
#include <vector>
#include <array>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <tuple>
#include <limits>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace daq::details {
  
  using namespace util::quantities::time_literals;
  using util::quantities::intervals::nanoseconds;
  
  // --- BEGIN -- Board information management ---------------------------------
  
  /// Information of the setup of a V1730 readout board.
  struct BoardSetup_t {
    
    ///< V1730B has 16 channels (can become a template parameter).
    static constexpr std::size_t NBoardChannels = 16U;
    
    /// Special value to mark the absence of fragment ID information.
    static constexpr unsigned int NoFragmentID
      = std::numeric_limits<unsigned int>::max();
    
    /// Special settings for one channel on the board.
    struct ChannelSetup_t {
      
      /// Associated off-line (LArSoft) channel ID.
      raw::Channel_t channelID = sbn::V1730channelConfiguration::NoChannelID;
      
      /// Whether the channel should be forcedly ignored or acquired.
      std::optional<bool> forcedSkip;
      
      /// Save the channel only when including the global trigger time.
      bool onGlobalOnly = false;
      
      std::uint16_t minSpan = 0; ///< Minimum acceptable waveform span.
      
      /// Category of this channel (will become instance name).
      std::string category;
      
      /// Returns if this channel is associated to a off-line channel number.
      bool hasChannel() const { return isChannel(channelID); }
      
      /// Whether this channel is requested to be skipped.
      bool mustSkip() const { return forcedSkip.value_or(false); }
      
      /// Whether this channel is requested to be saved.
      bool mustSave() const { return !forcedSkip.value_or(true); }
      
      /// Returns whether this is the default configuration.
      bool isDefault() const
        {
          // C++20: equality operator might be automatically generated
          // return *this == ChannelSetup_t{}; 
          return !hasChannel()
            && !forcedSkip && !onGlobalOnly && category.empty();
        }
      
      /// Returns if `channel` is a valid channel ID.
      static bool isChannel(raw::Channel_t channel)
        { return channel != sbn::V1730channelConfiguration::NoChannelID; }
      
    }; // ChannelSetup_t
    
    using AllChannelSetup_t = std::array<ChannelSetup_t, NBoardChannels>;
    
    std::string name; ///< Board name as specified in DAQ configuration.
    
    /// ID of the DAQ fragment associated to this board.
    unsigned int fragmentID = NoFragmentID;
    
    nanoseconds triggerDelay = 0_ns; ///< Delay from global trigger to TTT set.
    
    nanoseconds TTTresetDelay = 0_ns; ///< Delay from TTT reset issue to enact.
    
    /// Set of settings channel by channel.
    AllChannelSetup_t channelSettings;
    
    /// Returns whether this object contains a valid fragment ID.
    constexpr bool hasFragmentID() const { return fragmentID != NoFragmentID; }
    
  }; // struct BoardSetup_t

  std::ostream& operator<< (std::ostream&, BoardSetup_t::ChannelSetup_t const&);
  
  /// Derivative information of a V1730 readout board.
  struct BoardFacts_t {
    nanoseconds preTriggerTime; ///< How long from waveform start to PMT board trigger.
  }; // struct BoardFacts_t

  class BoardInfoLookup;

  std::ostream& operator<< (std::ostream&, BoardInfoLookup const&);

  // --- END -- Board information management -----------------------------------

  
  template <std::size_t KeyIndex = 0U>
  struct CompareWithKey;
  
  template <typename Coll, typename Key, typename KeyExtractor>
  typename Coll::value_type const* binarySearch
    (Coll const& coll, Key const& key, KeyExtractor&& extractKey);
    
  /// Finds the element in `coll` with the item `KeyIndex` equal to `key`.
  /// @return a pointer to the item found, or `nullptr` if none
  template <std::size_t KeyIndex = 0U, typename Coll, typename Key>
  typename Coll::value_type const* binarySearch
    (Coll const& coll, Key const& key);

  
} // namespace daq::details



// -----------------------------------------------------------------------------
/// Utility class for fast lookup of board data by fragment ID.
class daq::details::BoardInfoLookup {
  
    public:
  using FragmentID_t = unsigned int;
  
  /// Record of information about a readout board.
  struct BoardInfo_t {
    FragmentID_t fragmentID = std::numeric_limits<FragmentID_t>::max();
    BoardSetup_t const* setup = nullptr;
    sbn::V1730Configuration const* config = nullptr;
    BoardFacts_t facts;
    
    bool operator< (BoardInfo_t const& other) const
      { return fragmentID < other.fragmentID; }
  }; // struct BoardInfo_t

  using Database_t = std::vector<BoardInfo_t>;
  
  
  explicit BoardInfoLookup(Database_t&& database)
    : fDatabase(sortDatabase(std::move(database)))
    {}
  
  /// Returns the number of registered boards.
  unsigned int nBoardInfo() const { return fDatabase.size(); }
  
  /// Returns an object that can be iterated for all available information.
  auto const& allBoardInfo() const { return fDatabase; } // C++20: replace with std::span<>
  
  /// Returns a pointer to the information for board with `fragmentID`.
  BoardInfo_t const* findBoardInfo(FragmentID_t fragmentID) const;
  
  /// Returns a pointer to the setup information for board with `fragmentID`.
  BoardSetup_t const* findBoardSetup(FragmentID_t fragmentID) const
    {
      auto const* pInfo = findBoardInfo(fragmentID);
      return pInfo? pInfo->setup: nullptr;
    }
  
  /// Returns a pointer to the configuration of board with `fragmentID`.
  sbn::V1730Configuration const* findBoardConfig(FragmentID_t fragmentID) const
    {
      auto const* pInfo = findBoardInfo(fragmentID);
      return pInfo? pInfo->config: nullptr;
    }
  
  
    private:
  
  Database_t const fDatabase;
  
  /// Sorts the database structure in the argument and returns it.
  static Database_t sortDatabase(Database_t&& db)
    { std::sort(db.begin(), db.end()); return std::move(db); }
  
}; // class daq::details::BoardInfoLookup



// -----------------------------------------------------------------------------
// --- inline implementation
// -----------------------------------------------------------------------------
inline auto daq::details::BoardInfoLookup::findBoardInfo
  (FragmentID_t fragmentID) const -> BoardInfo_t const*
{
  return binarySearch(fDatabase, fragmentID, std::mem_fn(&BoardInfo_t::fragmentID));
}


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename Coll, typename Key, typename KeyExtractor>
typename Coll::value_type const* daq::details::binarySearch
  (Coll const& coll, Key const& key, KeyExtractor&& extractKey)
{
  using std::begin, std::end;
  auto const b = begin(coll);
  auto const e = end(coll);
  auto it = std::lower_bound(
    b, e, key,
    [&extractKey](auto const& a, Key const& key){ return extractKey(a) < key; }
    );
  return ((it != e) && (extractKey(*it) == key))? &*it: nullptr;
} // daq::details::binarySearch()


// -----------------------------------------------------------------------------
template <std::size_t KeyIndex /* = 0U */, typename Coll, typename Key>
typename Coll::value_type const* daq::details::binarySearch
  (Coll const& coll, Key const& key)
{
  return binarySearch
    (coll, key, [](auto const& a){ return std::get<KeyIndex>(a); });
} // daq::details::binarySearch()


// -----------------------------------------------------------------------------


#endif // SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_PMTDECODERUTILS_H
