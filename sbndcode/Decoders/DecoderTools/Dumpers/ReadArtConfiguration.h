/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/ReadArtConfiguration.h
 * @brief  Utilities to extract _art_ FHiCL configuration from different sources.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 23, 2023
 * @see    sbndcode/Decoders/DecoderTools/Dumpers/ReadArtConfiguration.cxx
 * 
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_READARTCONFIGURATION_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_READARTCONFIGURATION_H


// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetID.h"
#include "fhiclcpp/ParameterSetRegistry.h" // also defines ParameterSetID hash
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TFile.h"

// C/C++ standard libraries
#include <map>


// -----------------------------------------------------------------------------
namespace util {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Reads and returns the _art_ configuration stored in `sourceDir`.
   * @param file ROOT file where the configuration is stored
   * @return the full configuration
   * @throw cet::exception (category: `"readConfigurationFromArtFile"`) on error
   * @see `readConfigurationFromArtPrincipal()`
   * 
   * The configuration is expected to be stored by _art_ in the way it does
   * for _art_ ROOT files.
   * 
   * The configuration is returned as a map of parameter set ID to parameter
   * set.
   */
  std::map<fhicl::ParameterSetID, fhicl::ParameterSet>
    readConfigurationFromArtFile(TFile& file);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Reads and returns the complete _art_ configuration in the principal.
   * @tparam Principal type of framework principal class (e.g. `art::Event`)
   * @param principal the "principal" _art_ object to read the information from
   * @return the full configuration, as map: process name -> FHiCL parameter set
   * @throw cet::exception (category: `"readConfigurationFromArtPrincipal"`)
   *        on error
   * @see `readConfigurationFromArtFile()`
   * 
   * Compared to `readConfigurationFromArtFile()`, this function relays the same
   * information after it has been conveniently extracted by the framework.
   * A "principal" is framework jargon for `art::Event`, `art::SubRun` or
   * `art::Run` (all derived from `art::DataViewImpl`).
   * Therefore, this function can be called e.g. in `beginRun()` hook of a
   * module using its `art::Run` argument, or in a `analyze()` hook using
   * its `art::Event` argument.
   * 
   * This function is supposed to be more solid than
   * `readConfigurationFromArtFile()` because it relies less on internals of
   * how the framework works. , this function relays the same
   * information after it has been conveniently extracted by the framework.
   * Also, this function should also be compatible with `gallery::Event` too
   * (which is the reason why it is implemented as template instead of taking
   * a `art::DataViewImpl`, which is _art_-specific), making the name of this
   * function a misnomer.
   * 
   * The configuration is returned as a map of process names to parameter sets
   * (note that the key is different from the one returned by
   * `readConfigurationFromArtFile()`).
   */
  template <typename Principal>
  std::map<std::string, fhicl::ParameterSet>
  readConfigurationFromArtPrincipal(Principal const& principal);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util



// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename Principal>
std::map<std::string, fhicl::ParameterSet>
util::readConfigurationFromArtPrincipal(Principal const& principal) {
  
  std::map<std::string, fhicl::ParameterSet> configMap;
  
  for (art::ProcessConfiguration const& procConfig: principal.processHistory())
  {
    
    fhicl::ParameterSet config;
    if (!fhicl::ParameterSetRegistry::get(procConfig.parameterSetID(), config))
    {
      // this would be, as far as I understand, a logic error
      throw cet::exception("readConfigurationFromArtPrincipal")
        << "Configuration of process '" << procConfig.processName()
        << "' can't be found!\n";
    }
    
    configMap[procConfig.processName()] = std::move(config);
  } // for
  
  return configMap;
  
} // util::readConfigurationFromArtPrincipal()


// -----------------------------------------------------------------------------


#endif // SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_READARTCONFIGURATION_H
