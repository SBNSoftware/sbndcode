/**
 * @file    DaqDecoderSBNDPMT_module.cc
 * @brief   Produces `raw::OpDetWaveform` from V1730 artDAQ data fragments.
 * @authors Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date    Jan 28, 2023
 * 
 */


// SBND/SBN libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/OpDetWaveformMeta.h" 
#include "sbndcode/Decoders/DecoderTools/Dumpers/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker 
#include "sbndcode/Decoders/DecoderTools/details/PMTDecoderUtils.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/FragmentDumper.h"
#include "sbndcode/Decoders/ChannelMapping/ISBNDChannelMap.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/BinaryDumpUtils.h" // sbnd::ns::util::bin()
#include "sbndcode/Decoders/DecoderTools/Dumpers/ArtHandleTrackerManager.h"

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h" 
#include "sbnobj/Common/Trigger/BeamBits.h"
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/intervals_fhicl.h" // for nanoseconds in FHiCL
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

// artDAQ
#include "art/Framework/Core/ModuleMacros.h"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TTree.h"

// C/C++ standard libraries
#include <memory>
#include <ostream>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <optional>
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;
using detinfo::timescales::electronics_time;

//------------------------------------------------------------------------------
namespace sbnd { class DaqDecoderSBNDPMT; }
/**
 * @brief Produces `raw::OpDetWaveform` from V1730 artDAQ data fragments.
 * 
 * The module can read fragments from CAEN V1730 readout boards delivered by
 * artDAQ. It produces optical detector waveforms in LArSoft standards.
 * 
 * This decoder must support both a off-line mode (for storage and downstream
 * processing) and a on-line mode (for monitoring).
 * In particular, the on-line workflow is such that it may not be possible to
 * access the FHiCL configuration of the job and therefore the PMT configuration
 * data.
 * 
 * 
 * Configuration
 * --------------
 * 
 * The set of supported parameters can be seen on command line by running
 * `lar --print-description DaqDecoderSBNDPMT`.
 * 
 * Description of the configuration parameters:
 * * `FragmentsLabels` (list of input tags): a list of possible input data
 *     products. A valid input data product contain a list of fragments from
 *     V1730. All input tags are tried. If only one is available, its content
 *     is used for decoding; otherwise, an exception is thrown.
 * * `DiagnosticOutput` (flag, default: `false`): enables additional console
 *     output, including dumping of the fragments (that is huge output).
 * * `PMTconfigTag` (data product tag, optional): if specified, the pre-trigger
 *     buffer duration is read from there; although optional, it is strongly
 *     recommended that this information be provided, since it is essential for
 *     the correct timing of the PMT waveforms
 * * `BoardSetup` (list of board setup information): each entry specifies some
 *     information about a specific readout board; the boards are identified by
 *     their name; if a board is found in input that has no setup information,
 *     some time corrections are not applied 
 *     Each entry is in the form of a table:
 *     * `Name` (string, mandatory): the name of the board
 *        (e.g. `"sbndpmtwwtop01"`); this is use to match the setup
 *        information to a fragment ID in the PMT configuration.
 *     * `FragmentID` (integral, optional): if specified, allows the corrections
 *       using setup information to be applied even when no PMT configuration is
 *       provided (if neither PMT configuration nor setup information including
 *       `FragmentID` is available, no time correction is applied).
 *     * `SpecialChannels` (list of configurations): each entry described
 *       special features and settings of a specific channel, identified by
 *       its V1730B board channel number (`0` to `15`). The supported keys are:
 *         * `ChannelIndex` (integer, mandatory) the index of the channel in
 *           the board, from `0` to `15`; note that for this module to decode
 *           this channel index, the channel itself must be enabled on the
 *           readout board.
 *         * `Skip` (flag, optional): if specified, tells whether to always skip
 *           this channel or whether always store it. If not specified, a
 *           channel will be always stored unless it is not associated to any
 *           channel number (see `Channel` below). Note that if the channel is
 *           not enabled in the readout board, this option has no effect.
 *         * `Channel` (integral, optional): if specified, it overrides the
 *           channel number used for the waveforms on this channel; if not
 *           specified, the channel number is obtained from the channel mapping
 *           database. Note that insisting to save a waveform without an
 *           assigned channel number (from database or from this configuration)
 *           will trigger an exception.
 * * `RequireKnownBoards` (flag, default: `true`): if set, the readout boards
 *     in input must each have a setup configuration (`BoardSetup`) *and* must
 *     be present in the PMT DAQ configuration, or an exception is thrown;
 *     if not set, readout boards in input will be processed even when their
 *     setup or DAQ configuration is not known, at the cost of fewer corrections
 *     on the timestamps (which should then be considered unreliable).
 * * `RequireBoardConfig` (flag, default: `true`): if set, the readout boards
 *     which have a setup (`BoardSetup`) are required to be included in the DAQ
 *     configuration of the input file, or an exception is thrown; if not set,
 *     missing readout boards are unnoticed.
 * * `DataTrees` (list of strings, default: none): list of data trees to be
 *     produced; if none (default), then `TFileService` is not required.
 * * `SkipWaveforms` (flag, default: `false`) if set, waveforms won't be
 *     produced; this is intended as a debugging option for jobs where only the
 *     `DataTrees` are desired.
 * * `DropRawDataAfterUse` (flag, default: `true`): at the end of processing,
 *     the framework will be asked to remove the PMT data fragment from memory.
 *     Set this to `false` in the unlikely case where raw PMT fragments are
 *     still needed after decoding.
 * * `LogCategory` (string, default: `DaqDecoderSBNDPMT`): name of the message
 *     facility category where the output is sent.
 * 
 * 
 * Requirements
 * -------------
 * 
 * Services required include:
 * 
 * * `ISBNDChannelMap` for the association of fragments to LArSoft channel ID;
 * * `DetectorClocksService` for the correct decoding of the time stamps
 *   (always required, even when dumbed-down timestamp decoding is requested);
 * * `TFileService` only if the production of trees or plots is requested.
 * 
 * 
 * Output data products
 * ---------------------
 * 
 * * `std::vector<raw::OpDetWaveform>` (empty instance name): decoded waveforms
 *   from the "regular" PMT channels.
 *   Produced only if `SkipWaveforms` is `false`.
 * * `std::vector<raw::OpDetWaveform>` (non-empty instance name): decoded
 *   waveforms from the special PMT channels as configured in `BoardSetup`.
 *   Produced only if `SkipWaveforms` is `false`.
 * * `std::vector<sbn::OpDetWaveformMeta>`: for regular waveforms, metadata
 *   of each produced waveform (in the same order).
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>`: for regular
 *   waveforms, association with its metadata object (also in same order).
 *   Produced only if `SkipWaveforms` is `false`.
 * 
 * 
 * Waveform time stamp
 * --------------------
 * *Acquire from Fragment Metadata
 *  Fragment Metadata timestamp is stamped on CAEN1730 Boardreader
 *  
 * 
 * Data trees
 * -----------
 * 
 * The module supports the following ROOT trees production on demand:
 * 
 * * `PMTfragments`: data pertaining a single fragment; each entry is about a
 *   single fragment, and it includes full event ID, event time stamp (from
 *   _art_, i.e. the one assigned to the event by artDAQ), the ID of the
 *   fragment the entry is describing, and then the rest of the data, including:
 *     * `TTT` (64-bit integer): (Extended) Trigger Time Tag value, in readout
 *       board ticks (each worth 8 ns) from the last board reset;
 *       currently the value includes the 31-bit counter and in addition the
 *       overflow bit as MSB; the overflow bit is set by the readout board
 *       the first time the counter passes its limit (2^31) and wraps, and never
 *       cleared until the next board reset.
 *       While the tree has room for the 48-bit version (ETTT), the rest of the
 *       decoder does not yet support it.
 *     * `fragTime` (64-bit signed integer), `fragTimeSec` (32-bit signed
 *       integer): the timestamp of the PMT fragment, assigned by the board
 *       reader constructing the fragment.
 *     * `fragCount` (unsigned integer): the counter of the fragment within a
 *       readout board (its "event counter").
 *     * `waveformTime` (double): time assigned to the waveforms from this
 *       fragment (electronics time scale).
 *     * `waveformSize` (unsigned integer): number of ticks for the waveforms
 *       from this fragment.
 * 
 * 
 * Technical notes
 * ----------------
 * 
 * PMT configuration is optional, in the sense that it can be omitted; in that
 * case, some standard values will be used for it. This kind of setup may be
 * good for monitoring, but it does not meet the requirements for physics
 * analyses.
 * For a board to be served, an entry of that board must be present in the
 * module configuration (`BoardSetup`). It is an error for a fragment in input
 * not to have an entry for the corresponding board setup.
 * 
 * The module extracts the needed information and matches it into a
 * sort-of-database keyed by fragment ID, so that it can be quickly applied
 * when decoding a fragment. The matching is performed by board name.
 * 
 * 
 * ### Proto-waveforms
 * 
 * The module has grown complex enough that some decisions need to be taken much
 * later after a waveform is created, but using information available only
 * during the creation. For example, a waveform can be ultimately routed to a
 * different data product depending on which readout board it comes from,
 * but that information is lost by the time the waveform needs to be written.
 * Rather than attempting to recreate that information, it is saved when still
 * available and travels together with the waveform itself in a data structure
 * called "proto-waveform". The extra-information is eventually discarded when
 * putting the actual waveform into the _art_ event.
 * 
 * 
 * ### Processing pipeline
 * 
 * Processing happens according to the following structure (in `produce()`):
 * 1. pre-processing: currently nothing
 * 2. processing of each board data independently: at this level, all the
 *    buffers from the 16 channels of a single board are processed together
 *    (`processBoardFragments()`)
 *     1. the configuration and parameters specific to this board are fetched
 *     2. each data fragment is processed independently: at this level, data
 *        from all 16 channels _at a given time_ are processed together,
 *        producing up to 16 proto-waveforms
 * 3. post-processing of proto-waveforms:
 *     * sorting by time (as opposed as roughly by channel, as they come)
 * 4. conversion to data products and output
 * 
 * Glossary
 * ---------
 * 
 * * **setup**, **[PMT] configuration**: this is jargon specific to this module.
 *     Information about a readout board can come from two sources: the "setup"
 *     is information included in the `BoardSetup` configuration list of this
 *     module; the "PMT configuration" is information included in the DAQ
 *     configuration that is delivered via `PMTconfigTag`.
 * * **TAI** (International Atomic Time): a time standard defining a universal
 *     time with a precision higher than it will ever matter for SBND.
 * * **TTT**: trigger time tag, from the V1730 event record (31 bits); may be:
 * * **ETTT**: extended trigger time tag, from the V1730 event record (48 bits).
 * 
 */
class sbnd::DaqDecoderSBNDPMT: public art::EDProducer {
  
  // --- BEGIN -- some debugging tree declarations -----------------------------
  
  /// Enumerate the supported data trees.
  enum class DataTrees: std::size_t {
    Fragments, ///< Information about fragments
    N          ///< Counter.
  };
  using TreeNameList_t
    = std::array<std::string, static_cast<std::size_t>(DataTrees::N)>;
  static TreeNameList_t const TreeNames;
  
  /// Returns a string with all supported tree names.
  static std::string listTreeNames(std::string const& sep = "\n");
  
  // --- END ---- some debugging tree declarations -----------------------------
  
    public:
  
  // --- BEGIN -- public data types --------------------------------------------
  
  using nanoseconds = util::quantities::intervals::nanoseconds; ///< Alias.
  
  /// Data structure for trigger time.
  struct SplitTimestamp_t {
    static_assert(sizeof(int) >= 4U);
    struct Split_t {
        int seconds = std::numeric_limits<int>::min(); ///< The ongoing second.
        /// Nanoseconds from the start of current second.
        unsigned int nanoseconds = std::numeric_limits<unsigned int>::max();
    }; // Split_t
    
    /// Trigger time in nanoseconds from The Epoch.
    long long int time = std::numeric_limits<long long int>::min();
    /// Trigger time in nanoseconds from The Epoch (in components).
    Split_t split;
    
    constexpr SplitTimestamp_t() = default;
    constexpr SplitTimestamp_t(int sec, unsigned int ns);
    constexpr SplitTimestamp_t(long long int triggerTime);
  }; // SplitTimestamp_t
  
  // --- END ---- public data types --------------------------------------------
  
  
  // --- BEGIN -- FHiCL configuration ------------------------------------------
  
  /// Configuration of the V1730 readout board setup.
  struct BoardSetupConfig {
    
    struct ChannelSetupConfig {
      
      fhicl::Atom<unsigned short int> ChannelIndex {
        fhicl::Name("ChannelIndex"),
        fhicl::Comment
          ("index of the channel on the board these settings pertain")
        // mandatory
        };
      
      fhicl::OptionalAtom<bool> Skip {
        fhicl::Name("Skip"),
        fhicl::Comment(
          "set to true to force skipping this channel, false to force saving it"
          )
        };
      
      fhicl::Atom<raw::Channel_t> Channel {
        fhicl::Name("Channel"),
        fhicl::Comment("Off-line channel ID associated to this board channel"),
        sbn::V1730channelConfiguration::NoChannelID
        };
      
    }; // struct ChannelSetupConfig
    
    
    fhicl::Atom<std::string> Name {
      fhicl::Name("Name"),
      fhicl::Comment("board name, as specified in the DAQ configuration")
      };
    
    fhicl::OptionalAtom<unsigned int> FragmentID {
      fhicl::Name("FragmentID"),
      fhicl::Comment("ID of the fragments associated with the board")
      };
    
    fhicl::OptionalSequence<fhicl::Table<ChannelSetupConfig>> SpecialChannels {
      fhicl::Name("SpecialChannels"),
      fhicl::Comment("special settings for selected channels on the board")
      };
    
  }; // struct BoardSetupConfig
  
  
  /// Main module configuration.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> FragmentsLabels {
      Name("FragmentsLabels"),
      Comment("data product candidates with the PMT fragments from DAQ"),
      std::vector<art::InputTag>{ "daq:CAENV1730", "daq:ContainerCAENV1730" }
      };
    
    fhicl::Atom<bool> SurviveExceptions {
      Name("SurviveExceptions"),
      Comment
        ("when the decoding module throws an exception, print a message and move on"),
      true // default
      };
    
    fhicl::Atom<bool> DiagnosticOutput {
      Name("DiagnosticOutput"),
      Comment("enable additional console output"),
      false // default
      };
    
    fhicl::Atom<bool> PacketDump {
      Name("PacketDump"),
      Comment("enable dump of the whole V1730 data (huge)"),
      false // default
      };
    
    fhicl::Atom<bool> RequireKnownBoards {
      Name("RequireKnownBoards"),
      Comment
        ("all readout boards in input must be known (setup+PMT configuration)"),
      true
      };
    
    fhicl::Atom<bool> RequireBoardConfig {
      Name("RequireBoardConfig"),
      Comment
        ("all readout boards in setup must have a matching PMT configuration"),
      true
      };

    fhicl::OptionalAtom<art::InputTag> PMTconfigTag {
      Name("PMTconfigTag"),
      Comment("input tag for the PMT readout board configuration information")
      };
    
    fhicl::Sequence
      <fhicl::TableAs<daq::details::BoardSetup_t, BoardSetupConfig>>
    BoardSetup {
      Name("BoardSetup"),
      Comment("list of the setup settings for all relevant V1730 boards")
      };
    
    fhicl::Sequence<std::string> DataTrees {
      fhicl::Name("DataTrees"),
      fhicl::Comment
        ("produces the specified ROOT trees (" + listTreeNames(",") + ")"),
      std::vector<std::string>{} // default
      };
    
    fhicl::Atom<bool> DropRawDataAfterUse {
      Name("DropRawDataAfterUse"),
      Comment("drop PMT data fragments from memory after use"),
      true // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category for message stream"),
      "PMTDecoder" // default
      };
    
  }; // Config
  
  using Parameters = art::EDProducer::Table<Config>;

  static constexpr electronics_time NoTimestamp = std::numeric_limits<electronics_time>::min(); 
 
  // --- END ---- FHiCL configuration ------------------------------------------
  
  /// Constructor.
  explicit DaqDecoderSBNDPMT(Parameters const& params);
  
  /// On a new run: cache PMT configuration information.
  void beginRun(art::Run& run) override;
  
  /// Processes the event.
  void produce(art::Event& event) override;
  
  /// Prints a end-of-job message.
  void endJob() override;
  
  
    private:
  
  /// Type of setup of all channels in a readout board.
  using AllChannelSetup_t = daq::details::BoardSetup_t::AllChannelSetup_t;
  
  /// Collection of useful information from fragment data.
  struct FragmentInfo_t {
    artdaq::Fragment::fragment_id_t fragmentID
      = std::numeric_limits<artdaq::Fragment::fragment_id_t>::max();
    artdaq::Fragment::timestamp_t fragmentTimestamp;
    unsigned int eventCounter = 0U;
    std::uint32_t TTT = 0U;
    std::uint16_t enabledChannels = 0U;
    std::size_t nSamplesPerChannel = 0U;
    std::uint16_t const* data = nullptr;
  }; // FragmentInfo_t
  
  /// Information used in decoding from a board.
  struct NeededBoardInfo_t {
    static AllChannelSetup_t const DefaultChannelSetup; // default-initialized
    
    std::string const name;
    nanoseconds bufferLength;
    AllChannelSetup_t const* specialChannelSetup = nullptr;
    
    AllChannelSetup_t const& channelSetup() const
      { return specialChannelSetup? *specialChannelSetup: DefaultChannelSetup; }
    
  }; // NeededBoardInfo_t
  
  /// All the information collected about a waveform (with the waveform itself).
  struct ProtoWaveform_t {
    
    raw::OpDetWaveform waveform; ///< The complete waveform.
    
    /// Pointer to the settings for the channel this waveform belongs to.
    daq::details::BoardSetup_t::ChannelSetup_t const* channelSetup = nullptr;
    
    ///< Lowest sample in the original waveform.
    std::uint16_t minSample = std::numeric_limits<std::uint16_t>::max();
    ///< Highest sample in the original waveform.
    std::uint16_t maxSample =  std::numeric_limits<std::uint16_t>::max();
    
    /// Ordering: the same as the contained waveform.
    bool operator< (ProtoWaveform_t const& than) const
      { return waveform < than.waveform; }
    
  }; // struct ProtoWaveform_t

  /// Type of map: category to collection of waveforms under that category.
  using WaveformsByCategory_t
    = std::map<std::string, std::vector<raw::OpDetWaveform>>;
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  ///< List of candidate data products with artDAQ data fragments.
  std::vector<art::InputTag> const fInputTags;
  
  bool const fSurviveExceptions; ///< Whether to "ignore" errors.
  
  /// If true will spew endless messages to output.
  bool const fDiagnosticOutput;
  
  bool const fPacketDump; ///< Dump V1730 data.
  
  /// Whether info on all input boards is required.
  bool const fRequireKnownBoards;
  
  /// Whether setup info on all boards is required.
  bool const fRequireBoardConfig;

  /// Input tag of the PMT configuration.
  std::optional<art::InputTag> const fPMTconfigTag;
  
  /// All board setup settings.
  std::vector<daq::details::BoardSetup_t> const fBoardSetup;
  
  /// Clear fragment data product cache after use.
  bool const fDropRawDataAfterUse;
  
  std::string const fLogCategory; ///< Message facility category.
  
  // --- END ---- Configuration parameters -------------------------------------

  
  // --- BEGIN -- Services -----------------------------------------------------
  
  /// Interface to LArSoft configuration for detector timing.
  detinfo::DetectorTimings const fDetTimings;

  /// Fragment/channel mapping database.
  sbndDB::ISBNDChannelMap const& fChannelMap;

  // --- END ---- Services -----------------------------------------------------


  // --- BEGIN -- Cached values ------------------------------------------------
  
  /// Duration of the optical detector readout sampling tick (i.e. 2 ns; hush!).
  nanoseconds const fOpticalTick;
  
  bool const fSaveRegularWaveforms; ///< Whether regular waveforms are saved.
  
  // --- END ---- Cached values ------------------------------------------------
  
  
  // --- BEGIN -- Per-run data cache -------------------------------------------
  
  /// Find the information on a readout boards by fragment ID.
  std::optional<daq::details::BoardInfoLookup> fBoardInfoLookup;
  
  // --- END -- Per-run data cache ---------------------------------------------
  
  
  unsigned int fNFailures = 0U; ///< Number of event failures encountered.
  
  /// Instance name associated to the regular PMT waveforms.
  static std::string const RegularWaveformCategory;
  
  
  // --- BEGIN -- PMT readout configuration ------------------------------------
  
  /// Returns whether PMT configuration information is expected to be available.
  bool hasPMTconfiguration() const { return fPMTconfigTag.has_value(); }
  
  /// Updates the PMT configuration cache. How? Dunno. Placeholder.
  bool UpdatePMTConfiguration(sbn::PMTconfiguration const* PMTconfig);
  

  /**
   * @brief Returns a lookup object with board setup and configuration info.
   * @param PMTconfig the PMT configuration, if available
   * @return an object working like lookup table for all fragment information
   * 
   * This method merges the setup information from the module configuration with
   * the PMT configuration specified in the argument, and returns an object
   * that can look up all the information as a single record, with the
   * fragment ID as key. In addition, a few intermediate quantities ("facts",
   * see `BoardFacts_t`) are computed and stored in this object.
   * 
   * If a fragment ID is missing, it means that no PMT configuration was
   * provided and that the setup information did not include a fragment ID.
   * If some information (configuration or setup) is missing, the "facts"
   * depending on the missing information will have default values.
   */
  daq::details::BoardInfoLookup matchBoardConfigurationAndSetup
    (sbn::PMTconfiguration const* PMTconfig) const;
  
  /// Puts together all the needed information for a board.
  NeededBoardInfo_t fetchNeededBoardInfo(
    daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo,
    unsigned int fragmentID
    ) const;


  // --- END -- PMT readout configuration --------------------------------------


  // --- BEGIN -- Trees and their data ---------------------------------------
  
  /// Data structure for basic event information in simple ROOT trees.
  struct TreeData_EventID_t {
    unsigned int run;    ///< Run number.
    unsigned int subrun; ///< Subrun number.
    unsigned int event;  ///< Event number.
    SplitTimestamp_t timestamp; ///< Event timestamp (seconds from the epoch).
  }; // TreeData_EventID_t
  
  /// Structure collecting all data for a fragment ROOT tree.
  struct TreeFragment_t {
    struct Data_t: public TreeData_EventID_t {
      
      int fragmentID = 0; ///< ID of the fragment of this entry.
      
      unsigned int fragCount = 0U; ///< Counter of the fragment in the board.
      
      ///< Trigger time tag from the fragment.
      unsigned long int TriggerTimeTag = 0;
      
      SplitTimestamp_t fragTime; ///< PMT fragment time stamp.
      
      /// Time assigned to the waveforms.
      double waveformTime = std::numeric_limits<double>::lowest();
      
      unsigned int waveformSize = 0U; ///< Ticks in the waveforms.
      
    }; // Data_t
    
    Data_t data;
    TTree* tree = nullptr;
  }; // TreeFragment_t
  
  
  std::unique_ptr<TreeData_EventID_t> fEventInfo; ///< Event ID for trees.
  
  ///< Tree with fragment information.
  std::unique_ptr<TreeFragment_t> fTreeFragment;
  
  // --- END ---- Trees and their data -----------------------------------------
  
  
  // --- BEGIN -- Timestamps ---------------------------------------------------
  
  /// Returns the timestamp for the waveforms in the specified fragment.
  electronics_time fragmentWaveformTimestamp(
    FragmentInfo_t const& fragInfo,
    NeededBoardInfo_t const& boardInfo
    ) const;
  
  /**
   * @brief Returns the timestamp for the waveforms in the specified fragment.
   * This method assigns the timestamp from fragment metadata
   */
  electronics_time fragmentWaveformTimestampFromMetadata(
    FragmentInfo_t const& fragInfo,
    NeededBoardInfo_t const& boardInfo
    ) const;
 
  /// Returns the timestamp difference `a - b`.
  static long long int timestampDiff(std::uint64_t a, std::uint64_t b)
    { return static_cast<long long int>(a) - static_cast<long long int>(b); }
  
  /// Returns the fragment ID to be used with databases.
  static constexpr std::size_t effectivePMTboardFragmentID
    (artdaq::Fragment::fragment_id_t id)
    { return id & 0x0fff; }
  
  // --- END ---- Timestamps ---------------------------------------------------
  
  
  // --- BEGIN -- Input data management ----------------------------------------
  
  /// Tracks _art_ data products and removes their cached data on demand.
  // As often happens, caches need to be mutable. Also, not thread-safe.
  mutable util::ArtHandleTrackerManager<art::Event> fDataCacheRemover;
  
  /// Reads the fragments to be processed.
  std::vector<art::Handle<artdaq::Fragments>> readInputFragments(art::Event const& event) const;
  
  /// Throws an exception if `artdaqFragment` is not of type `CAEN1730`.
  void checkFragmentType(artdaq::Fragment const& artdaqFragment) const;
  
  /// Converts a fragment into a fragment collection
  /// (dispatcher based on fragment type).
  artdaq::FragmentPtrs makeFragmentCollection
    (artdaq::Fragment const& sourceFragment) const;

  /// Converts a plain fragment into a fragment collection.
  artdaq::FragmentPtrs makeFragmentCollectionFromFragment
    (artdaq::Fragment const& sourceFragment) const;

  /// Converts a container fragment into a fragment collection.
  artdaq::FragmentPtrs makeFragmentCollectionFromContainerFragment
    (artdaq::Fragment const& sourceFragment) const;

  /// Extracts waveforms from the specified fragments from a board.
  std::vector<ProtoWaveform_t> processBoardFragments
    ( artdaq::FragmentPtrs const& artdaqFragment);
  
  // --- END ---- Input data management ----------------------------------------
  
  
  // --- BEGIN -- Output waveforms ---------------------------------------------

  /// Sorts in place the specified waveforms in channel order, then in time.
  void sortWaveforms(std::vector<ProtoWaveform_t>& waveforms) const;
  
  /// Creates and returns waveform metadata.
  std::vector<sbn::OpDetWaveformMeta> createWaveformMetadata(std::vector<raw::OpDetWaveform> const& waveforms) const;
  
  /// Creates and returns waveform metadata associations.
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
  createWaveformMetadataAssociations(
    std::vector<raw::OpDetWaveform> const& waveforms,
    art::Event const& event, std::string const& instanceName = ""
    ) const;

  /**
   * @brief Converts the proto-waveforms into the final waveforms to be saved.
   * @param protoWaveforms the waveforms to be converted
   * @return a map: instance name to its collection of waveforms to save
   * 
   * Note that `protoWaveforms` is depleted in the process.
   */
  WaveformsByCategory_t prepareOutputWaveforms(std::vector<ProtoWaveform_t>&& protoWaveforms) const;

  /// Creates metadata and associations (if needed) for the regular waveforms.
  std::pair<
    std::vector<sbn::OpDetWaveformMeta>,
    art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
    >
  createRegularWaveformMetadata(
    WaveformsByCategory_t const& waveformProducts,
    art::Event const& event     ) const;

  // --- END ---- Output waveforms ---------------------------------------------
  
  using BoardID_t = short int; ///< Type used internally to represent board ID.
  
  /**
   * @brief Create waveforms and fills trees for the specified artDAQ fragment.
   * @param artdaqFragment the fragment to process
   * @param boardInfo board information needed, from configuration/setup
   * @param triggerTime absolute time of the trigger
   * @return collection of PMT waveforms from the fragment
   * 
   * This method fills the information for the PMT fragment tree
   * (`fillPMTfragmentTree()`) and creates PMT waveforms from the fragment data
   * (`createFragmentWaveforms()`).
   */
  std::vector<ProtoWaveform_t> processFragment(
    artdaq::Fragment const& artdaqFragment,
    NeededBoardInfo_t const& boardInfo
    );
  
  /**
   * @brief Creates `raw::OpDetWaveform` objects from the fragment data.
   * @param fragInfo information extracted from the fragment
   * @param channelSetup settings channel by channel
   * @param timeStamp timestamp of the waveforms in the fragment
   * @return collection of newly created `raw::OpDetWaveform`
   * 
   * All fragment information needed is enclosed in `fragInfo`
   * (`extractFragmentInfo()`). The timestamp can be obtained with a call to
   * `fragmentWaveformTimestamp()`.
   */
  std::vector<ProtoWaveform_t> createFragmentWaveforms(
    FragmentInfo_t const& fragInfo,
    AllChannelSetup_t const& channelSetup,
    electronics_time const timeStamp
    ) const;
  
  /// Extracts useful information from fragment data.
  FragmentInfo_t extractFragmentInfo
    (artdaq::Fragment const& artdaqFragment) const;
  
  /// Extracts the fragment ID (i.e. board ID) from the specified `fragment`.
  static BoardID_t extractFragmentBoardID(artdaq::Fragment const& fragment);
  
  /// Returns the board information for this fragment.
  NeededBoardInfo_t neededBoardInfo
    (artdaq::Fragment::fragment_id_t fragment_id) const;
  
  /// Throws an exception if the configuration of boards shows errors.
  void checkBoardSetup
    (std::vector<daq::details::BoardSetup_t> const& allBoardSetup) const;
  
  
  // --- BEGIN -- Tree-related methods -----------------------------------------
  
  /// Initializes all requested data trees.
  void initTrees(std::vector<std::string> const& treeNames);
  
  /// Initializes the event ID part of a tree.
  void initEventIDtree(TTree& tree, TreeData_EventID_t& data);
  
  /// Initializes the fragment data tree (`fTreeFragment`).
  void initFragmentsTree();

  /// Fills the base information of a tree data entry from an _art_ event.
  void fillTreeEventID
    (art::Event const& event, TreeData_EventID_t& treeData) const;

  /// Assigns the cached event information to the specified tree data.
  void assignEventInfo(TreeData_EventID_t& treeData) const;
  
  /// Fills the PMT fragment tree with the specified information
  /// (additional information needs to have been set already).
  void fillPMTfragmentTree(
    FragmentInfo_t const& fragInfo,
    electronics_time waveformTimestamp
    );
  
  /// Returns the name of the specified tree.
  static std::string const& treeName(DataTrees treeID);

  /// Static initialization.
  static TreeNameList_t initTreeNames();
  
  // --- END ---- Tree-related methods -----------------------------------------

  friend struct dumpChannel;
  friend std::ostream& operator<< (std::ostream&, ProtoWaveform_t const&);
  
}; // sbnd::DaqDecoderSBNDPMT


//------------------------------------------------------------------------------
// --- implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Moves the contend of `src` into the end of `dest`.
  template <typename T>
  std::vector<T>& appendTo(std::vector<T>& dest, std::vector<T>&& src) {
    if (dest.empty()) dest = std::move(src);
    else {
      dest.reserve(dest.size() + src.size());
      std::move(src.begin(), src.end(), std::back_inserter(dest));
    }
    src.clear();
    return dest;
  } // appendTo()
  
  /// Returns whether an element `value` is found in `coll`.
  template <typename Coll, typename T>
  bool contains(Coll const& coll, T const& value) { // C++20: remove
    auto const cend = end(coll);
    return std::find(begin(coll), cend, value) != cend;
  } // contains()
  
  /// Moves the content of `data` into a `std::unique_ptr`.
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& data)
    { return std::make_unique<T>(std::move(data)); }

} // local namespace


//------------------------------------------------------------------------------
namespace sbnd {

  /// Special function `fhicl::TableAs` uses to convert BoardSetupConfig.
  daq::details::BoardSetup_t convert
    (DaqDecoderSBNDPMT::BoardSetupConfig const& config)
  {
    
    daq::details::BoardSetup_t bs {
        config.Name()                                          // name
      , config.FragmentID()
          .value_or(daq::details::BoardSetup_t::NoFragmentID)  // fragmentID
      };
    
    // set the special configuration for the board channels that have one;
    // the others will stay with the default one (which implies that a channel
    // is saved in the standard data product, and only if it has a channel ID
    // entry from the database)
    for (
      DaqDecoderSBNDPMT::BoardSetupConfig::ChannelSetupConfig const& chConfig
      : config.SpecialChannels().value_or(
        std::vector<DaqDecoderSBNDPMT::BoardSetupConfig::ChannelSetupConfig>{}
        )
    ) {
      try {
        bs.channelSettings.at(chConfig.ChannelIndex()) = {
            chConfig.Channel()              // channelID
          , chConfig.Skip()                 // forcedSkip
          };
      }
      catch (std::out_of_range const&) {
        throw art::Exception(art::errors::Configuration)
          << "Configuration requested for invalid channel index "
          << chConfig.ChannelIndex() << " of the board '"
          << config.Name() << "'\n";
      }
    } // for
    
    return bs;
  } // convert(BoardSetupConfig)


  struct dumpChannel {
    DaqDecoderSBNDPMT::ProtoWaveform_t const& wf;
    dumpChannel(DaqDecoderSBNDPMT::ProtoWaveform_t const& wf): wf{ wf } {}
  };
  
  std::ostream& operator<< (std::ostream& out, dumpChannel const& d) {
    return out << (d.wf.channelSetup->category.empty()? std::dec: std::hex)
      << d.wf.waveform.ChannelNumber() << std::dec;
  }

} // namespace sbnd


//------------------------------------------------------------------------------
// --- sbnd::DaqDecoderSBNDPMT::SplitTimestamp_t
//------------------------------------------------------------------------------
constexpr sbnd::DaqDecoderSBNDPMT::SplitTimestamp_t::SplitTimestamp_t
  (int sec, unsigned int ns)
  : time { static_cast<long long int>(sec) * 1'000'000'000LL + ns }
  , split { sec, ns }
{}


//------------------------------------------------------------------------------
constexpr sbnd::DaqDecoderSBNDPMT::SplitTimestamp_t::SplitTimestamp_t
  (long long int triggerTime)
  : time { triggerTime }
  , split {
      static_cast<int>(time / 1'000'000'000), // seconds
      static_cast<unsigned int>(time % 1'000'000'000) // nanoseconds
    }
{}


//------------------------------------------------------------------------------
namespace sbnd {
  
  std::ostream& operator<<
    (std::ostream& out, DaqDecoderSBNDPMT::SplitTimestamp_t const& time)
  {
    out << time.split.seconds << '.'
      << std::setfill('0') << std::setw(9) << time.split.nanoseconds;
    return out;
  } // operator<< (std::ostream&, PMTDecoder::SplitTimestamp_t)
  
} // namespace sbnd


//------------------------------------------------------------------------------
// --- sbnd::DaqDecoderSBNDPMT
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// --- static definitions
//------------------------------------------------------------------------------
sbnd::DaqDecoderSBNDPMT::AllChannelSetup_t const
sbnd::DaqDecoderSBNDPMT::NeededBoardInfo_t::DefaultChannelSetup;

std::string const sbnd::DaqDecoderSBNDPMT::RegularWaveformCategory{ "" };


//------------------------------------------------------------------------------
sbnd::DaqDecoderSBNDPMT::TreeNameList_t const
sbnd::DaqDecoderSBNDPMT::TreeNames
  = sbnd::DaqDecoderSBNDPMT::initTreeNames();

auto sbnd::DaqDecoderSBNDPMT::initTreeNames() -> TreeNameList_t {
  TreeNameList_t names;
  names[static_cast<std::size_t>(DataTrees::Fragments)] = "PMTfragments";
  return names;
} // sbnd::DaqDecoderSBNDPMT::initTreeNames()


//------------------------------------------------------------------------------
std::string const& sbnd::DaqDecoderSBNDPMT::treeName(DataTrees treeID)
  { return TreeNames[static_cast<std::size_t>(treeID)]; }


//------------------------------------------------------------------------------
std::string sbnd::DaqDecoderSBNDPMT::listTreeNames
  (std::string const& sep /* = " " */)
{
  std::string l;
  for (std::string const& name: TreeNames) {
    if (!l.empty()) l += sep;
    l += '\'';
    l += name;
    l += '\'';
  } // for
  return l;
} // sbnd::DaqDecoderSBNDPMT::listTreeNames()


//------------------------------------------------------------------------------
// --- implementation
//------------------------------------------------------------------------------
sbnd::DaqDecoderSBNDPMT::DaqDecoderSBNDPMT(Parameters const& params)
  : art::EDProducer(params)
  , fInputTags{ params().FragmentsLabels() }
  , fSurviveExceptions{ params().SurviveExceptions() }
  , fDiagnosticOutput{ params().DiagnosticOutput() }
  , fPacketDump{ params().PacketDump() }
  , fRequireKnownBoards{ params().RequireKnownBoards() }
  , fRequireBoardConfig{ params().RequireBoardConfig() }
  , fPMTconfigTag{ params().PMTconfigTag() }
  , fBoardSetup{ params().BoardSetup() }
  , fDropRawDataAfterUse{ params().DropRawDataAfterUse() }
  , fLogCategory{ params().LogCategory() }
  , fDetTimings
    { art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
  , fChannelMap{ *(art::ServiceHandle<sbndDB::ISBNDChannelMap const>{}) }
  , fOpticalTick{ fDetTimings.OpticalClockPeriod() }
  , fSaveRegularWaveforms
    {  }
{
  //
  // configuration check
  //
  checkBoardSetup(fBoardSetup); // throws on error
  
  //
  // consumed data products declaration
  //
  for (art::InputTag const& inputTag: fInputTags)
    consumes<artdaq::Fragments>(inputTag);
  if (fPMTconfigTag) consumes<sbn::PMTconfiguration>(*fPMTconfigTag);
  //
  // produced data products declaration
  //
  
  // always write metadata for standard waveforms
  produces<std::vector<sbn::OpDetWaveformMeta>>();
  
  //
  // additional initialization
  //
  initTrees(params().DataTrees());
  
  
  //
  // configuration dump
  //
  mf::LogInfo log(fLogCategory);
  log << "Configuration:"
    << "\n * data from one of " << fInputTags.size() << " data products:";
  for (art::InputTag const& inputTag: fInputTags)
    log << " '" << inputTag.encode() << "'";
  log
    << "\n * boards with setup: " << fBoardSetup.size();
  if (fPMTconfigTag)
    log << "\n * PMT configuration from '" << fPMTconfigTag->encode() << "'";
  else 
    log << "\n * PMT configuration not used (and some corrections will be skipped)";
  if (fRequireKnownBoards) {
    log << "\n * all readout boards in input must be known (from `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name() << "`)"
      ;
  }
  else {
    log << "\n * readout boards with no information (from neither `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name()
      << "`) are processed at the best we can (skipping corrections)"
      ;
  }
  if (fRequireBoardConfig) {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` must appear in the PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }
  else {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` may lack a matching PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }

  //
  // sanity checks
  //
  
  if (fChannelMap.nPMTfragmentIDs() == 0) {
    throw cet::exception("DaqDecoderSBNDPMT")
      << "Channel mapping database does not report any PMT fragment ID!\n";
  } // if no PMT in database
  
  
} // sbnd::DaqDecoderSBNDPMT::DaqDecoderSBNDPMT()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::beginRun(art::Run& run) {
  
  //sbn::PMTconfiguration const* PMTconfig = fPMTconfigTag
  //  ? run.getPointerByLabel<sbn::PMTconfiguration>(*fPMTconfigTag): nullptr;
  sbn::PMTconfiguration const* PMTconfig = fPMTconfigTag
    ? run.getHandle<sbn::PMTconfiguration>(*fPMTconfigTag).product(): nullptr;
  
  UpdatePMTConfiguration(PMTconfig);

} // sbnd::DaqDecoderSBNDPMT::beginRun()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::produce(art::Event& event) {
 
  // ---------------------------------------------------------------------------
  // preparation
  //
  
  fDataCacheRemover.useEvent(event);
 
  //
  // event ID
  //
  
  // if needed, fill the record with the basic information of the event
  if (fEventInfo) fillTreeEventID(event, *fEventInfo);
 
  //LAN: event counter
  std::cout << std::endl;
  std::cout << "Processing Event: " << event.id().event() << std::endl;
  
  //
  // output data product initialization
  //
  std::vector<ProtoWaveform_t> protoWaveforms; // 
  
  // ---------------------------------------------------------------------------
  // pre-processing
  //
  
  
  // ---------------------------------------------------------------------------
  // processing
  //
  
  std::unordered_map<BoardID_t, unsigned int> boardCounts;
  bool duplicateBoards = false; 
  try { // catch-all

    auto fragmentHandles = readInputFragments(event);   

    for (auto const& handle : fragmentHandles) {

      for (auto const& frag : *handle) {   

	artdaq::FragmentPtrs const& fragmentCollection
        = makeFragmentCollection(frag);
      
        if (empty(fragmentCollection)) {
          mf::LogWarning("DaqDecoderSBNDPMT")
            << "Found a data fragment (ID=" << extractFragmentBoardID(frag)
            << ") containing no data.";
          continue;
        } // if no data
    
        BoardID_t const boardID = extractFragmentBoardID(*(fragmentCollection.front()));
        if (++boardCounts[boardID] > 1U) duplicateBoards = true;

        appendTo(
          protoWaveforms,
          processBoardFragments(fragmentCollection)
          );

      } // for fragment 
    } // for handle
    
  }
  catch (cet::exception const& e) {
    if (!fSurviveExceptions) throw;
    mf::LogError("DaqDecoderSBNDPMT")
      << "Error while attempting to decode PMT data:\n" << e.what() << '\n';
    protoWaveforms.clear();
    ++fNFailures;
  }
  catch (...) {
    if (!fSurviveExceptions) throw;
    mf::LogError("DaqDecoderSBNDPMT")
      << "Error while attempting to decode PMT data.\n";
    protoWaveforms.clear();
    ++fNFailures;
  }
 
  if (duplicateBoards) {
    mf::LogWarning log { "DaqDecoderSBNDPMT" };
    log << "found multiple data product entries for the same board:";
    for (auto [ boardID, count ]: boardCounts) {
      if (count < 2U) continue;
      log << " " << std::hex << boardID << std::dec << " (x" << count << ");";
    } // for
  } // if duplicate board fragments
  
  // we are done with the input: drop the caches
  // (if we were asked not to, no data is registered)
  fDataCacheRemover.removeCachedProducts();
  
  //
  // post-processing
  //
  /** LAN: COMMENT OUT SORTING WAVEFORMS
  sortWaveforms(protoWaveforms);
  END OF LAN**/ 
  
  std::map<std::string, std::vector<raw::OpDetWaveform>> waveformProducts = prepareOutputWaveforms(std::move(protoWaveforms));

  // create metadata for standard waveforms and, optionally, their associations
  auto [ waveformMeta, metaToWaveform ]
    = createRegularWaveformMetadata(waveformProducts, event);
  
  
  // put the data products into the event
  for (auto&& [ category, waveforms ]: waveformProducts) {
    mf::LogTrace(fLogCategory)
      << waveforms.size() << " PMT waveforms saved for "
      << (category.empty()? "standard": category) << " instance.";
    // the instance name is the category the waveforms belong to
    event.put(moveToUniquePtr(waveforms), category);
  }
  if (fSaveRegularWaveforms) event.put(moveToUniquePtr(metaToWaveform));
 
  event.put(moveToUniquePtr(waveformMeta)); // always save regular metadata
  
} // sbnd::DaqDecoderSBNDPMT::produce()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::endJob() {
  
  if (fNFailures > 0U) {
    mf::LogError(fLogCategory) << "Encountered errors on " << fNFailures
      << " events. Errors were ignored.";
  }
  
} // sbnd::DaqDecoderSBNDPMT::endJob()


//------------------------------------------------------------------------------
std::vector<art::Handle<artdaq::Fragments>> sbnd::DaqDecoderSBNDPMT::readInputFragments
  (art::Event const& event) const
{
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;

  for (art::InputTag const& inputTag: fInputTags) {
   
    mf::LogTrace(fLogCategory)
      << "DaqDecoderSBNDPMT trying data product: '" << inputTag.encode()
      << "'";
    
    art::Handle<artdaq::Fragments> thisHandle;

    if( !event.getByLabel<std::vector<artdaq::Fragment>>(inputTag, thisHandle) ) continue;
    if( !thisHandle.isValid() || thisHandle->empty() ) continue;  
 
    mf::LogTrace(fLogCategory)
      << "  => data product: '" << inputTag.encode() << "' is present and has "
      << thisHandle->size() << " entries";

    fragmentHandles.push_back(thisHandle); 
    
    if (fDropRawDataAfterUse) fDataCacheRemover.registerHandle(thisHandle);
  } // for

  if(fragmentHandles.empty()){
    cet::exception e { "DaqDecoderSBNDPMT" };
    e << "No suitable input data product found among:";
    for (art::InputTag const& inputTag: fInputTags)
      e << " '" << inputTag.encode() << "'";
    throw e << "\n";
  }

  return fragmentHandles;
} // sbnd::DaqDecoderSBNDPMT::readInputFragments()


//------------------------------------------------------------------------------
bool sbnd::DaqDecoderSBNDPMT::UpdatePMTConfiguration
  (sbn::PMTconfiguration const* PMTconfig)
{
  fBoardInfoLookup.emplace(matchBoardConfigurationAndSetup(PMTconfig));
  
  mf::LogDebug(fLogCategory)
    << "Board information as cached:\n" << *fBoardInfoLookup;
  
  return true;
} // sbnd::DaqDecoderSBNDPMT::UpdatePMTConfiguration()


auto sbnd::DaqDecoderSBNDPMT::matchBoardConfigurationAndSetup
  (sbn::PMTconfiguration const* PMTconfig) const
  -> daq::details::BoardInfoLookup
{
  /*
   * We need to support the case where no PMT configuration is known
   * (that is the standard situation in the online monitor).
   * The "strategy" is that in such cases we give up the correct time stamp
   * decoding; if the setup information contains a fragment ID, it may be
   * possible to do a little better, that is to use the setup information
   * (this is not possible without knowing the fragment ID that each bit of
   * setup information pertains).
   * 
   * So the cases for a board are:
   * * setup information is not present: encountering such a board will cause
   *   an exception to be thrown (implemented elsewhere)
   * * PMT configuration and setup present: full configuration
   *     * exception thrown if setup fragment ID is present and inconsistent
   * * PMT configuration not present: a general warning is printed;
   *     * boards with setup fragment ID information: add setup information
   *       to the "database" for the board: it will be used for partial
   *       timestamp reconstruction
   *     * boards without setup fragment ID information: board will not be
   *       added into the database; no specific correction will be performed;
   *       a warning is printed for each board
   * 
   */
  
  // dictionary of board configurations (if any)
  std::vector<std::pair<std::string, sbn::V1730Configuration const*>>
    configByName;
  if (PMTconfig) {
    if (!PMTconfig->boards.empty())
      configByName.reserve(PMTconfig->boards.size());
    for (sbn::V1730Configuration const& boardConfig: PMTconfig->boards)
      configByName.emplace_back(boardConfig.boardName, &boardConfig);
    std::sort(configByName.begin(), configByName.end()); // sorted by board name
  } // if we have configuration
  
  
  auto findPMTconfig = [this, &configByName]
    (std::string const& name) -> sbn::V1730Configuration const*
    {
      if (!hasPMTconfiguration()) return nullptr;
      auto const* ppBoardConfig
        = daq::details::binarySearch(configByName, name);
      if (!ppBoardConfig) {
        if (!fRequireBoardConfig) return nullptr;
        throw cet::exception("DaqDecoderSBNDPMT")
          << "No DAQ configuration found for PMT readout board '"
          << name << "'\n"
            "If this is expected, you may skip this check by setting "
            "DaqDecoderSBNDPMT module configuration `RequireBoardConfig`"
            " to `false`.\n"
          ;
      }
      return ppBoardConfig->second;
    }; // findPMTconfig()
  
  // the filling is driven by boards configured in the module
  // (which is the reason a setup entry is mandatory)
  daq::details::BoardInfoLookup::Database_t boardInfoByFragment;
  
  for (daq::details::BoardSetup_t const& boardSetup: fBoardSetup) {
    
    std::string const& boardName = boardSetup.name;
    
    sbn::V1730Configuration const* pBoardConfig = findPMTconfig(boardName);
    
    if (pBoardConfig) {
      // fragment ID from configuration and setup must match if both present
      if (boardSetup.hasFragmentID()
        && (boardSetup.fragmentID != pBoardConfig->fragmentID)
      ) {
        throw cet::exception("DaqDecoderSBNDPMT")
          << "Board '" << boardName << "' has fragment ID "
          << std::hex << pBoardConfig->fragmentID << std::dec
          << " but it is set up as "
          << std::hex << boardSetup.fragmentID << std::dec
          << "!\n";
      } // if fragment ID in setup
    }
    else {
      if (boardSetup.hasFragmentID()) {
        mf::LogPrint(fLogCategory)
          << "Board '" << boardName
          << "' has no configuration information;"
            " some time stamp corrections will be skipped.";
        // to avoid this, make a PMT configuration available from input file
      }
      else {
        mf::LogPrint(fLogCategory)
          << "Board '" << boardName
          << "' can't be associated to a fragment ID;"
            " its time stamp corrections will be skipped.";
        // to avoid this, add a `BoardSetup.FragmentID` entry for it in the
        // configuration of this module, or make a PMT configuration available
        continue; // no entry for this board at all
      }
    }
    
    unsigned int const fragmentID
      = pBoardConfig? pBoardConfig->fragmentID: boardSetup.fragmentID;
    assert(fragmentID != daq::details::BoardSetup_t::NoFragmentID);
    
    nanoseconds const preTriggerTime
      = pBoardConfig
      ? (pBoardConfig->bufferLength * (1.0f - pBoardConfig->postTriggerFrac))
        * fOpticalTick
      : nanoseconds{ 0.0 }
      ;
    
    daq::details::BoardFacts_t boardFacts {
      preTriggerTime  // ditto
      };
    
    boardInfoByFragment.push_back({
      fragmentID,               // fragmentID
      &boardSetup,              // setup
      pBoardConfig,             // config
      std::move(boardFacts)     // facts
      });
  } // for
  
  return daq::details::BoardInfoLookup{ std::move(boardInfoByFragment) };
  
} // sbnd::DaqDecoderSBNDPMT::matchBoardConfigurationAndSetup()


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::fetchNeededBoardInfo(
  daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo,
  unsigned int fragmentID
) const -> NeededBoardInfo_t {
  
  using util::quantities::intervals::microseconds;
  
  return NeededBoardInfo_t{
    // name
      ((boardInfo && boardInfo->config)
        ? boardInfo->config->boardName: ("<ID=" + std::to_string(fragmentID)))
    // bufferLength
    , ((boardInfo && boardInfo->config)
        ? boardInfo->config->bufferLength * fOpticalTick: nanoseconds{ 0.0 }
      )
    , ((boardInfo && boardInfo->setup)?
        &(boardInfo->setup->channelSettings): nullptr)
    };
      
} // sbnd::DaqDecoderSBNDPMT::fetchNeededBoardInfo()

//------------------------------------------------------------------------------
artdaq::FragmentPtrs sbnd::DaqDecoderSBNDPMT::makeFragmentCollection
  (artdaq::Fragment const& sourceFragment) const
{
  switch (sourceFragment.type()) {
    case sbndaq::FragmentType::CAENV1730:
      return makeFragmentCollectionFromFragment(sourceFragment);
    case artdaq::Fragment::ContainerFragmentType:
      return makeFragmentCollectionFromContainerFragment(sourceFragment);
    default:
      throw cet::exception("DaqDecoderSBNDPMT")
        << "Unexpected PMT data product fragment type: "
        << static_cast<int>(sourceFragment.type()) << " ('"
        << sbndaq::fragmentTypeToString
          (static_cast<sbndaq::FragmentType>(sourceFragment.type()))
        << "')\n";
  } // switch
} // sbnd::DaqDecoderSBNDPMT::makeFragmentCollection()


//------------------------------------------------------------------------------
artdaq::FragmentPtrs
sbnd::DaqDecoderSBNDPMT::makeFragmentCollectionFromFragment
  (artdaq::Fragment const& sourceFragment) const
{
  assert(sourceFragment.type() == sbndaq::FragmentType::CAENV1730);
  artdaq::FragmentPtrs fragColl;
  fragColl.push_back(std::make_unique<artdaq::Fragment>(sourceFragment));
  return fragColl;
} // sbnd::DaqDecoderSBNDPMT::makeFragmentCollectionFromFragment()


//------------------------------------------------------------------------------
artdaq::FragmentPtrs
sbnd::DaqDecoderSBNDPMT::makeFragmentCollectionFromContainerFragment
  (artdaq::Fragment const& sourceFragment) const
{
  assert(sourceFragment.type() == artdaq::Fragment::ContainerFragmentType);
  artdaq::ContainerFragment const containerFragment{ sourceFragment };
  
  if (containerFragment.block_count() == 0) return {};
    
  artdaq::FragmentPtrs fragColl;
  for (auto const iFrag: util::counter(containerFragment.block_count()))
    fragColl.push_back(containerFragment.at(iFrag));
  
  return fragColl;
} // sbnd::DaqDecoderSBNDPMT::makeFragmentCollectionFromContainerFragment()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::checkFragmentType
  (artdaq::Fragment const& artdaqFragment) const
{
  if (artdaqFragment.type() == sbndaq::FragmentType::CAENV1730) return;
  
  throw cet::exception("DaqDecoderSBNDPMT")
    << "Unexpected PMT fragment data type: '"
    << sbndaq::fragmentTypeToString
      (static_cast<sbndaq::FragmentType>(artdaqFragment.type()))
    << "'\n";
  
} // sbnd::DaqDecoderSBNDPMT::checkFragmentType


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::processBoardFragments(
  artdaq::FragmentPtrs const& artdaqFragments
) -> std::vector<ProtoWaveform_t> {
  
  if (artdaqFragments.empty()) return {};
  
  artdaq::Fragment const& referenceFragment = *(artdaqFragments.front());
  assert(&referenceFragment);
  
  checkFragmentType(referenceFragment);
  
  NeededBoardInfo_t const boardInfo
    = neededBoardInfo(artdaqFragments.front()->fragmentID());
  
  mf::LogTrace(fLogCategory)
    << " - " << boardInfo.name << ": " << artdaqFragments.size()
    << " fragments";
  
  std::vector<ProtoWaveform_t> waveforms;
  for (artdaq::FragmentPtr const& fragment: artdaqFragments)
    appendTo(waveforms, processFragment(*fragment, boardInfo));
  
  return { waveforms };
  
} // sbnd::DaqDecoderSBNDPMT::processBoardFragments()


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::processFragment(
  artdaq::Fragment const& artdaqFragment,
  NeededBoardInfo_t const& boardInfo
) -> std::vector<ProtoWaveform_t> {
  
  checkFragmentType(artdaqFragment);
  
  if (fPacketDump) {
    mf::LogVerbatim{ fLogCategory } << "PMT packet:"
      << "\n" << std::string(80, '-')
      << "\n" << sbndaq::dumpFragment(artdaqFragment)
      << "\n" << std::string(80, '-')
      ;
  } // if diagnostics
 
  //LAN: get fragment metadata 
  FragmentInfo_t const fragInfo = extractFragmentInfo(artdaqFragment);
  
  auto const timeStamp
    = fragmentWaveformTimestamp(fragInfo, boardInfo);

  if (fTreeFragment) fillPMTfragmentTree(fragInfo, timeStamp);
  
  return createFragmentWaveforms(fragInfo, boardInfo.channelSetup(), timeStamp);
} // sbnd::DaqDecoderSBNDPMT::processFragment()


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::createFragmentWaveforms(
  FragmentInfo_t const& fragInfo, AllChannelSetup_t const& channelSetup,
  electronics_time const timeStamp
) const -> std::vector<ProtoWaveform_t>
{

  std::vector<ProtoWaveform_t> protoWaveforms; // output collection
    
  std::optional<mf::LogVerbatim> diagOut;
  if (fDiagnosticOutput) diagOut.emplace(fLogCategory);
  
  sbndDB::DigitizerChannelChannelIDPairVec const& digitizerChannelVec
    = fChannelMap.getChannelIDPairVec
      (effectivePMTboardFragmentID(fragInfo.fragmentID))
    ;
  
  // allocate the vector outside the loop since we'll reuse it over and over
  std::vector<std::uint16_t> wvfm(fragInfo.nSamplesPerChannel);
  
  auto channelNumberToChannel
    = [&digitizerChannelVec](unsigned short int channelNumber) -> raw::Channel_t
    {
      for (auto const [ chNo, chID, _ ]: digitizerChannelVec) // too pythonic? 
        if (chNo == channelNumber) return chID;
      return sbn::V1730channelConfiguration::NoChannelID;
    };
  
  if (diagOut)
    (*diagOut) << "      " << digitizerChannelVec.size() << " channels:";
  
  std::size_t iNextChunk = 0;
  for (unsigned short int const channelNumber: util::counter(16U)) {
 
    if ((fragInfo.enabledChannels & (1 << channelNumber)) == 0) {
      if (diagOut){
        (*diagOut) << " " << channelNumber << " [disabled];";
      continue;}
    }
    
    std::size_t const iChunk = iNextChunk++;
    
    daq::details::BoardSetup_t::ChannelSetup_t const& thisChannelSetup
      = channelSetup.at(channelNumber); // this setup must be available
      
    if (thisChannelSetup.mustSkip()) {
      mf::LogTrace(fLogCategory)
        << "Channel number " << channelNumber << " of board 0x"
        << std::hex << fragInfo.fragmentID << std::dec
        << " is enabled but is requested to be skipped.";
      continue;
    }
    
    //
    // assign the offline channel ID
    //
 
    raw::Channel_t channel = thisChannelSetup.hasChannel()
      ? thisChannelSetup.channelID: channelNumberToChannel(channelNumber);
    
    assert(thisChannelSetup.isChannel(channel));
    std::cout << "Channel " << channel << " passes" << std::endl;   
    
    //
    // fill the waveform data
    //
    std::size_t const ch_offset = iChunk * fragInfo.nSamplesPerChannel;
    std::copy_n(fragInfo.data + ch_offset, wvfm.size(), wvfm.begin());
    
    //
    // create the proto-waveform
    //
    auto const [ itMin, itMax ] = std::minmax_element(wvfm.begin(), wvfm.end());
    protoWaveforms.push_back({ // create the waveform and its ancillary info
        raw::OpDetWaveform{ timeStamp.value(), channel, wvfm }  // waveform
      , &thisChannelSetup                                       // channelSetup
      , *itMin                                                  // minSample
      , *itMax                                                  // maxSample
      });
    
    mf::LogTrace log(fLogCategory);
    log << "PMT channel " << dumpChannel(protoWaveforms.back())
      << " has " << wvfm.size() << " samples (read from entry #" << iChunk
      << " in fragment data) starting at electronics time " << timeStamp;
    if (!protoWaveforms.back().channelSetup->category.empty()) {
      log << "; category: '" << protoWaveforms.back().channelSetup->category
        << "'";
    }
    
    //LAN
    std::cout << "PMT channel " << dumpChannel(protoWaveforms.back())
      << " has " << wvfm.size() << " samples (read from entry #" << iChunk
      << " in fragment data) starting at electronics time " << timeStamp
      << std::endl;
  } // for all channels in the board

  return protoWaveforms;
  
} // sbnd::DaqDecoderSBNDPMT::createFragmentWaveforms()

//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::prepareOutputWaveforms(
  std::vector<ProtoWaveform_t>&& protoWaveforms
) const -> WaveformsByCategory_t {
  
  // the output map will always include the regular waveforms...
  std::map<std::string, std::vector<raw::OpDetWaveform>> waveformProducts
    { { RegularWaveformCategory, std::vector<raw::OpDetWaveform>{} } };
  
  for (ProtoWaveform_t& waveform: protoWaveforms) {
    
    // select the destination data product for this waveform:
    auto const itOutputWaves
      = waveformProducts.find(waveform.channelSetup->category);
    if (itOutputWaves == waveformProducts.cend()) continue; // not to be saved
    
    // Set a new Timestamp
    waveform.waveform.SetTimeStamp(waveform.waveform.TimeStamp());
    itOutputWaves->second.push_back(std::move(waveform.waveform));
  } // for protowaveforms
  
  protoWaveforms.clear();
  
  return waveformProducts;
  
} // sbnd::DaqDecoderSBNDPMT::prepareOutputWaveforms()


//------------------------------------------------------------------------------
std::pair<
  std::vector<sbn::OpDetWaveformMeta>,
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
  >
sbnd::DaqDecoderSBNDPMT::createRegularWaveformMetadata(
  WaveformsByCategory_t const& waveformProducts,
  art::Event const& event
) const {

  // create metadata for standard waveforms and, optionally, their associations
  std::vector<raw::OpDetWaveform> const& regularWaveforms
    = waveformProducts.at(RegularWaveformCategory);
  
  std::vector<sbn::OpDetWaveformMeta> waveformMeta
    = createWaveformMetadata(regularWaveforms);
  
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> metaToWaveform;
  if (fSaveRegularWaveforms) {
    metaToWaveform
      = createWaveformMetadataAssociations(regularWaveforms, event);
  }
  
  return { std::move(waveformMeta), std::move(metaToWaveform) };
} // sbnd::DaqDecoderSBNDPMT::createRegularWaveformMetadata()

//------------------------------------------------------------------------------
std::vector<sbn::OpDetWaveformMeta>
sbnd::DaqDecoderSBNDPMT::createWaveformMetadata(
  std::vector<raw::OpDetWaveform> const& waveforms
) const {
  
  //apapadop
  /*
  sbn::OpDetWaveformMetaMaker const makeOpDetWaveformMeta
    { fOpticalTick, triggerInfo.relTriggerTime, triggerInfo.relBeamGateTime };
  */
  std::vector<sbn::OpDetWaveformMeta> waveformMeta;
  
  // apapadop
  /*
  for (auto const& waveform: waveforms)
    waveformMeta.push_back(makeOpDetWaveformMeta(waveform));
  */
  return waveformMeta;
} // sbnd::DaqDecoderSBNDPMT::createWaveformMetadata()


//------------------------------------------------------------------------------
art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
sbnd::DaqDecoderSBNDPMT::createWaveformMetadataAssociations(
  std::vector<raw::OpDetWaveform> const& waveforms,
  art::Event const& event, std::string const& instanceName /* = "" */
) const {
  
  art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr
    { event, instanceName };
  art::PtrMaker<sbn::OpDetWaveformMeta> const makeWaveMetaPtr
    { event, instanceName };
  
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> assns;
  for (std::size_t const iWaveform: util::counter(waveforms.size()))
    assns.addSingle(makeWaveformPtr(iWaveform), makeWaveMetaPtr(iWaveform));
  
  return assns;
} // sbnd::DaqDecoderSBNDPMT::createWaveformMetadataAssociations()



//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::extractFragmentBoardID
  (artdaq::Fragment const& fragment) -> BoardID_t
{
  return static_cast<BoardID_t>(fragment.fragmentID());
} // sbnd::DaqDecoderSBNDPMT::extractFragmentBoardID()


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::extractFragmentInfo
  (artdaq::Fragment const& artdaqFragment) const -> FragmentInfo_t
{
  //
  // fragment ID, timestamp and data begin
  //
  artdaq::Fragment::fragment_id_t const fragment_id
    = artdaqFragment.fragmentID();
  artdaq::Fragment::timestamp_t const fragmentTimestamp
    = artdaqFragment.timestamp();
  std::uint16_t const* data_begin = reinterpret_cast<std::uint16_t const*>
    (artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));

  //
  // event counter, trigger time tag, enabled channels
  //
  sbndaq::CAENV1730Fragment const fragment { artdaqFragment };
  sbndaq::CAENV1730FragmentMetadata const& metafrag = *(fragment.Metadata());
  sbndaq::CAENV1730EventHeader const& header = fragment.Event()->Header;
  
  unsigned int const eventCounter = header.eventCounter;
  
  unsigned int const TTT =  header.triggerTimeTag;
  
  std::uint16_t const enabledChannels = header.ChannelMask();

  //
  // samples per channel
  //
  // artDAQ size is measured in 4-byte words
  std::size_t const eventSize = header.eventSize * sizeof(std::uint32_t);
  
  constexpr std::size_t headerSize = sizeof(sbndaq::CAENV1730EventHeader);
  
  std::size_t const sampleDataSize = eventSize - headerSize;
  
  std::size_t const samplesInFragment = sampleDataSize / sizeof(std::uint16_t);
  
  std::size_t const nChannelsPerBoard = metafrag.nChannels;
  std::size_t const nSamplesPerChannel = samplesInFragment / nChannelsPerBoard;
  
  //
  // diagnostics dump
  //
  if (fDiagnosticOutput) {
    
    mf::LogVerbatim{ fLogCategory }
      << "----> PMT Fragment ID: " << std::hex << fragment_id << std::dec
        << ", size: " << eventSize << " B"
        << ", data size: " << samplesInFragment << " samples"
      << "\n    "
        << "  channels/board: " << nChannelsPerBoard
        << ", enabled: " << sbnd::ns::util::bin(enabledChannels)
        << ", samples/channel: " << nSamplesPerChannel
      << "\n    "
        << "  event counter: " << eventCounter
        << ", trigger time tag: " << TTT
        << ", time stamp: " << (fragmentTimestamp / 1'000'000'000UL)
          << "." << (fragmentTimestamp % 1'000'000'000UL) << " s"
      ;
  } // if diagnostics

  return { // C++20: write the member names explicitly
    fragment_id,
    fragmentTimestamp,
    eventCounter,
    TTT,
    enabledChannels,
    nSamplesPerChannel,
    data_begin
    };
  
} // sbnd::DaqDecoderSBNDPMT::extractFragmentInfo()


//------------------------------------------------------------------------------
electronics_time sbnd::DaqDecoderSBNDPMT::fragmentWaveformTimestamp(
  FragmentInfo_t const& fragInfo,
  NeededBoardInfo_t const& boardInfo
) const {
  
  // check availability of mapping for this board, otherwise can't do anything
  std::size_t const fragment_id = fragInfo.fragmentID;
  std::size_t const eff_fragment_id = effectivePMTboardFragmentID(fragment_id);
  
  if (!fChannelMap.hasPMTDigitizerID(eff_fragment_id)) {
    mf::LogError(fLogCategory)
      << "*** PMT could not find channel information for fragment: "
      << fragment_id
      ;
    return NoTimestamp;
  }
  
  return fragmentWaveformTimestampFromMetadata(fragInfo, boardInfo);
} // sbnd::DaqDecoderSBNDPMT::fragmentWaveformTimestamp()

//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::fragmentWaveformTimestampFromMetadata(
  FragmentInfo_t const&  fragInfo ,
  NeededBoardInfo_t const& /* boardInfo */
) const -> electronics_time {

  //LAN: Get the ns fraction starting after the current second in fragment timestamp
  auto const metadataTime{SplitTimestamp_t(fragInfo.fragmentTimestamp)};
  auto timestamp = util::quantities::nanosecond{static_cast<double>(metadataTime.split.nanoseconds)};
 
  std::cout << "Fragment Metadata Time = " << timestamp <<std::endl;

  return timestamp;
  
} // sbnd::DaqDecoderSBNDPMT::fragmentWaveformTimestampFromMetadata()


//------------------------------------------------------------------------------
auto sbnd::DaqDecoderSBNDPMT::neededBoardInfo
  (artdaq::Fragment::fragment_id_t fragment_id) const -> NeededBoardInfo_t
{
  
  assert(fBoardInfoLookup);

  daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo
    = fBoardInfoLookup->findBoardInfo(fragment_id);
  if (!boardInfo) {
    if (fRequireKnownBoards) {
      cet::exception e("DaqDecoderSBNDPMT");
      e << "Input fragment has ID " << fragment_id
        << " which has no associated board information (`BoardSetup`";
      if (!hasPMTconfiguration()) e << " + `.FragmentID`";
      throw e << ").\n";
    }
  }
  else {
    assert(boardInfo->fragmentID == fragment_id);
    assert(boardInfo->setup);
  }
  
  return fetchNeededBoardInfo(boardInfo, fragment_id);
} // sbnd::DaqDecoderSBNDPMT::neededBoardInfo()

//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::checkBoardSetup
  (std::vector<daq::details::BoardSetup_t> const& allBoardSetup) const
{
  //
  // duplicate special channel ID check
  //
  std::unordered_map<raw::Channel_t, unsigned int> assignedChannels;
  bool hasDuplicates = false;
  for (daq::details::BoardSetup_t const& setup: allBoardSetup) {
    for (daq::details::BoardSetup_t::ChannelSetup_t const& chSetup
      : setup.channelSettings
    ) {
      if (chSetup.mustSkip()) continue; // skipped channels do not matter anyway
      // special channels by definition have a non-empty category
      if (!chSetup.hasChannel() || chSetup.category.empty()) continue;
      if (++assignedChannels[chSetup.channelID] > 1) hasDuplicates = true;
    } // for all channels
  } // for all boards
  if (!hasDuplicates) return;
  
  art::Exception e { art::errors::Configuration };
  e << "Some channel ID are specified in multiple board special channel setup:";
  for (auto [ channelID, entries ]: assignedChannels) {
    if (entries <= 1) continue; // not a duplicate
    e << "\n - channel ID 0x" << std::hex << channelID << std::dec << ": "
      << entries << " times";
  } // for
  
  throw e;
} // sbnd::DaqDecoderSBNDPMT::checkBoardSetup()

//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::sortWaveforms
  (std::vector<ProtoWaveform_t>& waveforms) const
{
  auto byChannelThenTime = []
    (ProtoWaveform_t const& left, ProtoWaveform_t const& right)
    {
      return (left.waveform.ChannelNumber() != right.waveform.ChannelNumber())
        ? left.waveform.ChannelNumber() < right.waveform.ChannelNumber()
        : left.waveform.TimeStamp() < right.waveform.TimeStamp();
    };
  
  std::sort(waveforms.begin(), waveforms.end(), byChannelThenTime);

} // sbnd::DaqDecoderSBNDPMT::sortWaveforms()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::initTrees
  (std::vector<std::string> const& treeNames)
{
  
  auto findTree = [](std::string const& name)
    {
      return static_cast<DataTrees>(
          std::distance(TreeNames.begin(),
          std::find(TreeNames.begin(), TreeNames.end(), name))
        );
    };
  
  for (std::string const& name: treeNames) {
    switch (findTree(name)) {
      case DataTrees::Fragments: initFragmentsTree(); break;
      case DataTrees::N:
      default:
        throw cet::exception("DaqDecoderSBNDPMT")
          << "initTrees(): no data tree supported with name '" << name
          << "'.\n";
    } // switch
  } // for names
  
} // sbnd::DaqDecoderSBNDPMT::initTrees()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::initEventIDtree
  (TTree& tree, TreeData_EventID_t& data)
{
  
  if (!fEventInfo) fEventInfo = std::make_unique<TreeData_EventID_t>();
  
  tree.Branch("run", &data.run);
  tree.Branch("subrun", &data.subrun);
  tree.Branch("event", &data.event);
  tree.Branch("timestamp", &data.timestamp.time, "timestamp/L");
  
} // sbnd::DaqDecoderSBNDPMT::initEventIDtree()


//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::initFragmentsTree() {
  
  if (fTreeFragment) return;
  
  TTree* tree = art::ServiceHandle<art::TFileService>()
    ->make<TTree>("PMTfragments", "PMT fragment data");
  
  fTreeFragment = std::make_unique<TreeFragment_t>();
  fTreeFragment->tree = tree;
  auto& data = fTreeFragment->data;
  
  initEventIDtree(*tree, data);
  
  tree->Branch("fragmentID", &data.fragmentID);
  tree->Branch("fragCount", &data.fragCount);
  tree->Branch("fragTime", &data.fragTime.time, "fragTime/L"); // ROOT 6.24 can't detect 64-bit
  tree->Branch("fragTimeSec", &data.fragTime.split.seconds);
  tree->Branch("TTT", &data.TriggerTimeTag, "TTT/l"); // ROOT 6.24 can't detect 64-bit
  tree->Branch("waveformTime", &data.waveformTime);
  tree->Branch("waveformSize", &data.waveformSize);
  
} // sbnd::DaqDecoderSBNDPMT::initFragmentsTree()

//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::assignEventInfo
  (TreeData_EventID_t& treeData) const
{
  
  assert(fEventInfo);
  treeData = *fEventInfo; // nice slicing
  
} // sbnd::DaqDecoderSBNDPMT::assignEventInfo()

//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::fillTreeEventID
  (art::Event const& event, TreeData_EventID_t& treeData) const
{
  art::EventID const& id = event.id();
  treeData.run    = id.run();
  treeData.subrun = id.subRun();
  treeData.event  = id.event();
  
  art::Timestamp const& timestamp = event.time();
  treeData.timestamp
    = { static_cast<int>(timestamp.timeHigh()), timestamp.timeLow() };
  
} // sbnd::DaqDecoderSBNDPMT::fillTreeEventID()

//------------------------------------------------------------------------------
void sbnd::DaqDecoderSBNDPMT::fillPMTfragmentTree(
  FragmentInfo_t const& fragInfo,
  electronics_time waveformTimestamp
) {
  
  if (!fTreeFragment) return;
  
  fTreeFragment->data.fragmentID = fragInfo.fragmentID;
  fTreeFragment->data.fragCount = fragInfo.eventCounter;
  fTreeFragment->data.TriggerTimeTag = fragInfo.TTT;
  fTreeFragment->data.fragTime
    = { static_cast<long long int>(fragInfo.fragmentTimestamp) };
  fTreeFragment->data.waveformTime = waveformTimestamp.value();
  fTreeFragment->data.waveformSize = fragInfo.nSamplesPerChannel;
  assignEventInfo(fTreeFragment->data);
  fTreeFragment->tree->Fill();
  
} // sbnd::DaqDecoderSBNDPMT::fillPMTfragmentTree()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbnd::DaqDecoderSBNDPMT)


//------------------------------------------------------------------------------
