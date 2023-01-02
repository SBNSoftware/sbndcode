/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/FragmentDumper.h
 * @brief  Utility to dump a artDAQ fragment on screen.
 * @date   Jan 28, 2023
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @see    sbndcode/Decoders/DecoderTools/Dumpers/FragmentDumper.cxx
 * 
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_FRAGMENTDUMPER_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_FRAGMENTDUMPER_H


// framework libraries
#include "artdaq-core/Data/Fragment.hh"

// C/C++ standard libraries
#include <iosfwd>


// -----------------------------------------------------------------------------
namespace sbndaq::details {
  
  struct DumpFragWrap { artdaq::Fragment const& frag; };
  
  std::ostream& operator<< (std::ostream& out, DumpFragWrap const& wrap);
  
} // namespace sbndaq::details
// -----------------------------------------------------------------------------


namespace sbndaq {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Dump a artDAQ fragment into an output stream.
   * @param frag fragment to be dumped
   * @return an object that can be inserted in an output stream
   * 
   * This function allows to dump a artDAQ fragment.
   * The actual dump is performed by
   * `sbndaq::details::dumpFragment(details::DumpFragWrap)`.
   * 
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * artdaq::Fragment const frag;
   * // ... initialized somehow ...
   * mf::LogTrace("DecoderTool") << sbndaq::dumpFragment(frag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * This is actually a helper function that just passes the task.
   */
  details::DumpFragWrap dumpFragment(artdaq::Fragment const& frag);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace sbndaq


#endif // SBNDCODE_DECODE_DECODERTOOLS_DUMPERS_FRAGMENTDUMPER_H
