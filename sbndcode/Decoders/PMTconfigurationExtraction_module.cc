/**
 * @file   PMTconfigurationExtraction_module.cc
 * @brief  Writes PMT configuration from FHiCL into data product.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 */

// SBND libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTconfigurationExtractor.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/SBNDChannelMap.h"
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace sbnd { class PMTconfigurationExtraction; }
/**
 * @brief Writes PMT configuration from FHiCL into a data product.
 * 
 * This module reads the configuration related to PMT from the FHiCL
 * configuration of the input runs and puts it into each run as a data product.
 * 
 * Input
 * ------
 * 
 * This module requires any input with _art_ run objects in it.
 * The format expected for that configuration is defined in
 * `sbnd::PMTconfigurationExtractor`, which is the utility used for the actual
 * extraction.
 * 
 * 
 * Output
 * -------
 * 
 * A data product of type `sbn::PMTconfiguration` is placed in each run.
 * Note that the module itself does not enforce any coherence in the
 * configuration.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * The following configuration parameters are supported:
 * 
 * * **AssignOfflineChannelIDs** (flag, default: `true`): when set, service
 *     `ISBNDChannelMap` is used to track the channel numbers back to the
 *     LArSoft PMT channel IDs, and this information is saved together with
 *     the channel information; if the service is not available, this flag
 *     should be set to `false`, in which case the channel ID will be marked
 *     as unknown (`sbn::V1730channelConfiguration::NoChannelID`).
 * * **RequireConsistentPMTconfigurations** (flag, default: `true`): requires
 *     that all input files contain compatible PMT configuration; while this is
 *     in general desired during decoding, when mixing files with different
 *     runs in the same process this check might fail.
 * * **Verbose** (flag, default: `false`): if set to `true`, it will print in
 *     full the configuration of the PMT the first time it is read and each time
 *     a different one is found.
 * * **LogCategory** (string, default: `PMTconfigurationExtraction`):
 *     category tag used for messages to message facility.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * This module does not support multithreading, and _art_ does not provide
 * multithreading for its functionality anyway: the only action is performed
 * at run and input file level, and the only concurrency in _art_ is currently
 * (_art_ 3.7) at event level.
 * 
 */
class sbnd::PMTconfigurationExtraction: public art::EDProducer {
  
  /// Current PMT configuration (may be still unassigned).
  std::optional<sbn::PMTconfiguration> fPMTconfig;
  
  /// Pointer to the online channel mapping service.
  sbndDB::SBNDChannelMap const* fChannelMap = nullptr;
  
  /// Whether PMT configuration inconsistency is fatal.
  bool fRequireConsistency = true;
  
  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.
  
    public:
  
  /// Configuration of the module.
  struct Config {
    
    fhicl::Atom<bool> AssignOfflineChannelIDs {
      fhicl::Name("AssignOfflineChannelIDs"),
      fhicl::Comment
        ("retrieves LArSoft channel ID of each PMT; requires SBNDChannelMap service"),
      true // default
      };
    
    fhicl::Atom<bool> RequireConsistentPMTconfigurations {
      fhicl::Name("RequireConsistentPMTconfigurations"),
      fhicl::Comment
        ("checks that all input files carry a compatible PMT configuration"),
      true // default
      };
    
    fhicl::Atom<bool> Verbose {
      fhicl::Name("Verbose"),
      fhicl::Comment("print in full each new PMT configuration read"),
      false // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "PMTconfigurationExtraction" // default
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  /// Constructor: just reads the configuration.
  PMTconfigurationExtraction(Parameters const& config);
  
  
  /// Action on new run: configuration is written.
  virtual void beginRun(art::Run& run) override;
  
  /// Mandatory method, unused.
  virtual void produce(art::Event&) override {}
  
  
    private:
  
  /// Throws an exception if the `newConfig` is not compatible with the current.
  bool checkConsistency
    (sbn::PMTconfiguration const& newConfig, std::string const& srcName) const;
  
  
}; // sbnd::PMTconfigurationExtraction


// -----------------------------------------------------------------------------
sbnd::PMTconfigurationExtraction::PMTconfigurationExtraction
  (Parameters const& config)
  : art::EDProducer(config)
  , fChannelMap(config().AssignOfflineChannelIDs()
      ? art::ServiceHandle<sbndDB::SBNDChannelMap const>{}.get()
      : nullptr
    )
  , fRequireConsistency(config().RequireConsistentPMTconfigurations())
  , fVerbose(config().Verbose())
  , fLogCategory(config().LogCategory())
{
  
  // no consummation here (except for FHiCL configuration)
  
  produces<sbn::PMTconfiguration, art::InRun>();
  
} // sbnd::PMTconfigurationExtraction::PMTconfigurationExtraction()


// -----------------------------------------------------------------------------
void sbnd::PMTconfigurationExtraction::beginRun(art::Run& run) {
  
  sbn::PMTconfiguration config = extractPMTreadoutConfiguration
    (run, sbnd::PMTconfigurationExtractor{ *fChannelMap });
  
  checkConsistency(config, "run " + std::to_string(run.run()));
  if (!fPMTconfig.has_value() && fVerbose)
    mf::LogInfo(fLogCategory) << "PMT readout:" << config;
  
  fPMTconfig = std::move(config);
  
  // put a copy of the current configuration
  run.put(std::make_unique<sbn::PMTconfiguration>(fPMTconfig.value()));
  
} // sbnd::PMTconfigurationExtraction::beginRun()


// -----------------------------------------------------------------------------
bool sbnd::PMTconfigurationExtraction::checkConsistency
  (sbn::PMTconfiguration const& config, std::string const& srcName) const
{
  if (!fPMTconfig.has_value() || (*fPMTconfig == config)) return true;
  
  // see debug information for more details
  if (fRequireConsistency) {
    throw cet::exception("PMTconfigurationExtraction")
      << "Configuration from " << srcName
      << " is incompatible with the previously found configuration.\n"
      ;
  } // if consistency is required
  
  mf::LogWarning(fLogCategory)
    << "Configuration from " << srcName
    << " is incompatible with the previously found configuration.\n"
    ;
  if (fVerbose) mf::LogVerbatim(fLogCategory) << "PMT readout:" << config;
  return false;
  
  
} // sbnd::PMTconfigurationExtraction::checkConsistency()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbnd::PMTconfigurationExtraction)


// -----------------------------------------------------------------------------
