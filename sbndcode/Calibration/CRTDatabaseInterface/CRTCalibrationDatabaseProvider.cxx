// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Local
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabase.h"
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabaseProvider.h"

// Database interface helpers
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// C/C++ standard libraries
#include <string>
#include <vector>

//--------------------------------------------------------------------------------

sbndDB::CRTCalibrationDatabaseProvider::CRTCalibrationDatabaseProvider(
  const fhicl::ParameterSet& pset)
  : fVerbose{pset.get<bool>("Verbose", false)}
{
  fhicl::ParameterSet const databasepset{pset.get<fhicl::ParameterSet>("Database")};
  fDatabaseTag       = databasepset.get<std::string>("DatabaseTag");
  fDatabaseTimeStamp = databasepset.get<long>("DatabaseTimeStamp");
  fFEBTableName      = databasepset.get<std::string>("FEBTableName");
  fChannelTableName  = databasepset.get<std::string>("ChannelTableName");

  // if (fVerbose)
  //   mf::LogInfo(fLogCategory) << "Database tags for timing corrections:\n"
  //                             << "Cables corrections  " << fCRTCalibrationDatabaseTag << "\n";
}

// -------------------------------------------------------------------------------

uint64_t sbndDB::CRTCalibrationDatabaseProvider::RunToDatabaseTimestamp(uint32_t run) const
{
  uint64_t runNum = uint64_t(run);
  uint64_t timestamp = runNum + 1000000000;
  timestamp *= 1000000000;

  // if (fVerbose)
  //   mf::LogInfo(fLogCategory) << "Run " << runNum << " corrections from DB timestamp " << timestamp;

  return timestamp;
}

// -------------------------------------------------------------------------------

void sbndDB::CRTCalibrationDatabaseProvider::ReadCRTFEBCalibration(uint32_t run)
{
  lariov::DBFolder febTable(fFEBTableName, "", "", fDatabaseTag, true, false);

  febTable.UpdateData(fDatabaseTimeStamp);
  // bool ret = febTable.UpdateData(fDatabaseTimeStamp);
  // mf::LogDebug(fLogCategory) << fFEBTableName + " corrections" << (ret ? "" : " not")
  //                            << " updated for run " << run;
  // mf::LogTrace(fLogCategory)
  //   << "Fetched IoV [ " << db.CachedStart().DBStamp() << " ; " << db.CachedEnd().DBStamp()
  //   << " ] to cover t=" << RunToDatabaseTimestamp(run)
  //   << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp()
  //   << "]";

  std::vector<unsigned int> mac5List;
  if (int res = febTable.GetChannelList(mac5List); res != 0) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "GetChannelList() returned " << res << " on run " << run << " query in " << fFEBTableName << "\n";
  }

  if (mac5List.empty()) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "Got an empty channel list for run " << run << " in " << fFEBTableName << "\n";
  }

  for (auto mac5 : mac5List) {
    long moduleType = 0;
    int error = febTable.GetNamedChannelData(mac5, "type", moduleType);
    if (error)
      throw cet::exception("CRTCalibrationDatabaseProvider")
        << "Encountered error (code " << error
        << ") while trying to access 'type' on table " << fFEBTableName << "\n";
    fCRTFEBCalibrationData[mac5].moduleType = static_cast<int>(moduleType);
  }
}

// -----------------------------------------------------------------------------

void sbndDB::CRTCalibrationDatabaseProvider::readCRTCalibrationDatabase(const art::Run& run)
{
  fCRTFEBCalibrationData.clear();

  ReadCRTFEBCalibration(run.id().run());

  // if (fVerbose) {
  //   mf::LogInfo(fLogCategory) << "Dump information from database " << std::endl;
  //   mf::LogVerbatim(fLogCategory)
  //     << "channel, trigger cable delay, reset cable delay, laser corrections, muons corrections"
  //     << std::endl;
  //   for (auto const& [key, value] : fCRTCalibrationData) {
  //     mf::LogVerbatim(fLogCategory) << key << " " << value.breakoutBox << "," << std::endl;
  //   }
  // }
}
