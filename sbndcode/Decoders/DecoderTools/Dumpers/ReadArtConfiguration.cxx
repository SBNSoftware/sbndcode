/**
 * @file   sbndcode/Decoders/DecoderTools/ReadArtConfiguration.cxx
 * @brief  Utilities to extract _art_ FHiCL configuration from different sources.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * @see    sbndcode/Decoders/DecoderTools/Dumpers/ReadArtConfiguration.h
 * 
 */


// library header
#include "sbndcode/Decoders/DecoderTools/Dumpers/ReadArtConfiguration.h"

// framework libraries
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "fhiclcpp/ParameterSet.h"


// -----------------------------------------------------------------------------
std::map<fhicl::ParameterSetID, fhicl::ParameterSet>
util::readConfigurationFromArtFile(TFile& file)
{
  /*
   * This code is ripped from `fhiclcpp/ParameterSetRegistry.h` and
   * `lardata/DetectorInfoServices/DetectorClocksServiceStandard_service.cc`
   * (LArSoft v 9.17.0).
   * 
   * The special wrapped defines hooks to support a "virtual file system"
   * within a ROOT file data base, in the way that art knows and I do not.
   * So we bite it and accept a dependency against art_root_io.
   */
  art::SQLite3Wrapper sqliteDB(&file, "RootFileDB");
  
  auto* db = static_cast<sqlite3*>(sqliteDB);
  
  auto const throwOnSQLiteNotOK = [&db](std::string const& msg = {})
    {
      if (db == nullptr) {
        throw cet::exception("readConfigurationFromArtFile")
          << "Can't open SQLite database.";
      }
      auto const errcode = sqlite3_errcode(db);
      // Caller's responsibility to make sure this really is an error
      // and not (say) SQLITE_ROW or SQLITE_DONE:
      if (errcode == SQLITE_OK) return;
      throw cet::exception("readConfigurationFromArtFile")
        << "SQLite3 error (code" << errcode << "): "
        << sqlite3_errstr(errcode) << (msg.empty() ? "" : (": " + msg))
        << "\n";
    };
  
  
  sqlite3_stmt* stmt = nullptr;
  sqlite3_prepare_v2
    (sqliteDB, "SELECT ID, PSetBlob from ParameterSets;", -1, &stmt, nullptr);
  throwOnSQLiteNotOK("[SELECT ID, PSetBlob from ParameterSets;]");
  
  std::map<fhicl::ParameterSetID, fhicl::ParameterSet> config;
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    
    // reinterpretation: `unsigned char*` -> `char*`
    std::string const psetIDstr
      = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
    auto pset = fhicl::ParameterSet::make
      (reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1)));
    config.emplace(fhicl::ParameterSetID{ psetIDstr }, std::move(pset));
    
  } // while
  sqlite3_finalize(stmt);
  throwOnSQLiteNotOK("[SELECT ID, PSetBlob from ParameterSets;]");
  
  return config;
} // util::readConfigurationFromArtFile()


// -----------------------------------------------------------------------------
