/**
 * @file sbndcode/Decoders/DecoderTools/Dumpers/OpDetWaveformMeta.h
 * @brief Derivative information from `raw::OpDetWaveform` data.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date Jan 28, 2023
 * 
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_OPDETWAVEFORMMETA_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_OPDETWAVEFORMMETA_H


// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Utilities/FlagSet.h"

// C/C++ standard libraries
#include <limits>


// -----------------------------------------------------------------------------
namespace sbn { struct OpDetWaveformMeta; }
/**
 * @brief Derivative information from `raw::OpDetWaveform` data.
 * 
 * This objects stores some of the information from `raw::OpDetWaveform`,
 * with the notable exception of the content of the waveform.
 * In particular, it reports the time range covered by one waveform, and it may
 * store whether the range includes the beam gate opening time.
 * 
 * Times are in the same scale as for `raw::OpDetWaveform`, that is the
 * @ref DetectorClocksElectronicsStartTime "electronics time scale".
 * The exact type of that time is `detinfo::timescales::electronics_time`.
 */
struct sbn::OpDetWaveformMeta {
  
  // this is the same type used in LArSoft tracking flags: ROOT dictionaries
  // are already provided in there
  
  using Flags_t = ::util::flags::FlagSet<32U>; ///< Type of flag interface.
  
  /// Magic value denoting the absence of time information.
  static constexpr double NoTime = std::numeric_limits<double>::lowest();
  static constexpr raw::Channel_t NoChannel
    = std::numeric_limits<raw::Channel_t>::max();
  
  /// Namespace for bits in the flags.
  struct bits {
    
    using Flag_t = Flags_t::Flag_t; ///< Type for a single flag.
    
    /// Whether this time interval includes the nominal beam gate opening.
    static constexpr Flag_t WithBeamGate { 0 };
    
    /// Whether this time interval includes the hardware trigger time.
    static constexpr Flag_t WithTrigger { 1 };
    
    /// The first flag available for future use.
    static constexpr Flag_t NFlags { 2 };
    
  }; // struct bits
  
  
  // --- BEGIN -- Data members -------------------------------------------------
  
  raw::Channel_t channel = NoChannel; ///< ID of the PMT channel.
  
  std::size_t nSamples = 0; ///< Number of samples in the waveform.
  
  /// Time of the first sample in the waveform [us]
  double startTime = NoTime;
  /// Time at the end of the last sample in the waveform [us]
  double endTime = NoTime;
  
  Flags_t flags; ///< All flags (may be set, unset or undefined); see `bits`.
  
  // --- END ---- Data members -------------------------------------------------
  
  
  // --- BEGIN -- raw::OpDetWaveform interface replica -------------------------
  /**
   * @name `raw::OpDetWaveform` interface replica
   * 
   * Partial mirror of `raw::OpDetWaveform` interface is provided here to
   * facilitate static polymorphism ("metaprogramming").
   * This version reflects the interface as found in LArSoft 9.36.00.
   * 
   * These methods are supposed to return a value equivalent to the one from
   * the original waveform.
   * 
   * @note At this time `TimeStamp()` returns directly `startTime`, implicitly
   * assuming that `raw::OpDetWaveform::TimeStamp()` time scale matches the
   * @ref DetectorClocksElectronicsStartTime "electronics time scale".
   * If this assumption is broken, an additional data member needs to be added.
   */
  /// @{
  
  /// Returns the channnel ID of the PMT waveform.
  raw::Channel_t ChannelNumber() const { return channel; }
  
  /// Returns the timestamp at the start of the waveform.
  raw::TimeStamp_t TimeStamp() const { return startTime; }
  
  /// @}
  // --- END ---- raw::OpDetWaveform interface replica -------------------------
  
  
  // --- BEGIN -- Short bit managing interface ---------------------------------
  
  /**
   * @name Bit query
   * 
   * Note that each bit may be in an undefined state.
   * The definition of a bit may be tested directly; for example:
   * `info.flags.isDefined(bits::WithBeamGate)` returns whether `info` has the
   * beam gate bit defined.
   * 
   */
  /// @{
  
  /// Returns whether the time interval includes for sure the beam opening time.
  bool withBeamGate() const; // inline
  
  /// Returns whether the time interval for sure does not include the beam
  /// opening time.
  bool withoutBeamGate() const; // inline
  
  /// Returns whether the time interval includes for sure the trigger time.
  bool withTrigger() const; // inline
  
  /// Returns whether the time interval for sure does not include the trigger
  /// time.
  bool withoutTrigger() const; // inline
  
  /// @}
  
  
}; // sbn::OpDetWaveformMeta


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
inline bool sbn::OpDetWaveformMeta::withBeamGate() const {
  return flags.isDefined(bits::WithBeamGate) && flags.isSet(bits::WithBeamGate);
}

inline bool sbn::OpDetWaveformMeta::withoutBeamGate() const {
  return
    flags.isDefined(bits::WithBeamGate) && flags.isUnset(bits::WithBeamGate);
}

inline bool sbn::OpDetWaveformMeta::withTrigger() const {
  return flags.isDefined(bits::WithTrigger) && flags.isSet(bits::WithTrigger);
}

inline bool sbn::OpDetWaveformMeta::withoutTrigger() const {
  return flags.isDefined(bits::WithTrigger) && flags.isUnset(bits::WithTrigger);
}


// -----------------------------------------------------------------------------


#endif // SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_OPDETWAVEFORMMETA_H
