/**
 * @file   DumpPMTconfiguration_module.cc
 * @brief  Dumps on console the content of `sbn::PMTconfiguration` data product.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 */

// SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <set>
#include <string>
#include <sstream> // std::ostringstream


//------------------------------------------------------------------------------
namespace sbn { class DumpPMTconfiguration; }

/**
 * @brief Dumps on console the content of `sbn::PMTconfiguration` data product.
 * 
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `sbn::PMTconfiguration`: configuration of PMT in `art::Run` data product.
 *     See e.g. `sbnd::PMTconfigurationExtraction` module.
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description DumpPMTconfiguration`.
 * 
 * * `PMTconfigurationTag` (data product input tag): the tag identifying the
 *     data product to dump; this data product must be in `art::Run`.
 * * `Verbosity` (integral, default: maximum): verbosity level used in the
 *     dump; see `sbn::PMTconfiguration::dump()` for details.
 * * `SkipDuplicateRuns` (flag, default: `true`): multiple files can contain
 *     information from the same run; with this flag set, only the first time
 *     a run is encountered its PMT configuration is dumped; otherwise, each
 *     time a run is opened by _art_, its configuration is printed.
 * * `OutputCategory` (string, default: `DumpPMTconfiguration`): name of the
 *     message facility output stream to dump the information into
 * 
 */
class sbn::DumpPMTconfiguration: public art::EDAnalyzer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> PMTconfigurationTag {
      Name("PMTconfigurationTag"),
      Comment("tag of PMT configuration data product (from art::Run)")
      };

    fhicl::Atom<unsigned int> Verbosity {
      Name("Verbosity"),
      Comment("verbosity level [default: maximum]"),
      sbn::PMTconfiguration::MaxDumpVerbosity // default
      };

    fhicl::Atom<bool> SkipDuplicateRuns {
      Name("SkipDuplicateRuns"),
      Comment("print only one PMT configuration from each run"),
      true // default
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the category used for the output"),
      "DumpPMTconfiguration"
      };

  }; // struct Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit DumpPMTconfiguration(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Dumps the data product.
  virtual void beginRun(art::Run const& run) override;
  
  /// Does nothing, but it is mandatory.
  virtual void analyze(art::Event const& event) override {}
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fPMTconfigurationTag; ///< Input PMT configuration tag.
  
  unsigned int const fVerbosity; ///< Verbosity level used for dumping.
  
  bool const fSkipDuplicateRuns; ///< Print only once from each run.

  /// Category used for message facility stream.
  std::string const fOutputCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  std::set<art::RunID> fRuns; ///< Set of runs already encountered.
  
}; // sbn::DumpPMTconfiguration


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- sbn::DumpPMTconfiguration
//------------------------------------------------------------------------------
sbn::DumpPMTconfiguration::DumpPMTconfiguration
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fPMTconfigurationTag(config().PMTconfigurationTag())
  , fVerbosity          (config().Verbosity())
  , fSkipDuplicateRuns  (config().SkipDuplicateRuns())
  , fOutputCategory     (config().OutputCategory())
{
  
  consumes<sbn::PMTconfiguration, art::InRun>(fPMTconfigurationTag);
  
} // sbn::DumpPMTconfiguration::DumpPMTconfiguration()


//------------------------------------------------------------------------------
void sbn::DumpPMTconfiguration::beginRun(art::Run const& run) {
  
  if (fSkipDuplicateRuns) {
    art::RunID const& id = run.id();
    if (fRuns.count(id)) {
      mf::LogTrace(fOutputCategory) << id << " has already been encountered.";
      return;
    }
    fRuns.insert(id);
  } // if skip duplicates
  
  auto const& config
    = run.getProduct<sbn::PMTconfiguration>(fPMTconfigurationTag);
  
  std::ostringstream sstr;
  config.dump(sstr, "  ", "", fVerbosity);
  
  mf::LogVerbatim(fOutputCategory) << run.id() << ": " << sstr.str();
  

} // sbn::DumpPMTconfiguration::beginRun()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpPMTconfiguration)


//------------------------------------------------------------------------------
