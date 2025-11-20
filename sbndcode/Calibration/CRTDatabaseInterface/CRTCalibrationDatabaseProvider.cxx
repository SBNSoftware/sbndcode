// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Local
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabase.h"
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabaseProvider.h"

// C/C++ standard libraries
#include <string>
#include <vector>

//--------------------------------------------------------------------------------

sbndDB::CRTCalibrationDatabaseProvider::CRTCalibrationDatabaseProvider(const fhicl::ParameterSet& pset)
  : fVerbose(pset.get<bool>("Verbose", false))
{
  fhicl::ParameterSet const databasepset{pset.get<fhicl::ParameterSet>("Database")};
  fDatabaseTag       = databasepset.get<std::string>("DatabaseTag");
  fDatabaseTimeStamp = databasepset.get<long>("DatabaseTimeStamp");
  fFEBTableName      = databasepset.get<std::string>("FEBTableName");
  fChannelTableName  = databasepset.get<std::string>("ChannelTableName");

  if(fVerbose)
    mf::LogInfo("CRTCalibrationDatabaseProvider") << "CRT calibration database using tag:" << fDatabaseTag << "\n";
}

// -------------------------------------------------------------------------------

uint64_t sbndDB::CRTCalibrationDatabaseProvider::RunToDatabaseTimestamp(uint32_t run) const
{
  uint64_t runNum = uint64_t(run);
  uint64_t timestamp = runNum + 1000000000;
  timestamp *= 1000000000;

  if(fVerbose)
    mf::LogInfo("CRTCalibrationDatabaseProvider") << "Run " << runNum << " becomes timestamp: " << timestamp;

  return timestamp;
}

// -------------------------------------------------------------------------------

void sbndDB::CRTCalibrationDatabaseProvider::ReadCRTFEBCalibration(uint32_t run)
{
  lariov::DBFolder febTable(fFEBTableName, "", "", fDatabaseTag, true, false);

  febTable.UpdateData(fDatabaseTimeStamp);
  bool ret = febTable.UpdateData(fDatabaseTimeStamp);

  mf::LogDebug("CRTCalibrationDatabaseProvider") << fFEBTableName << (ret ? "" : " not")
                                                 << " updated for run " << run;

  mf::LogTrace("CRTCalibrationDatabaseProvider")
    << "Fetched IoV [ " << febTable.CachedStart().DBStamp() << " -> " << febTable.CachedEnd().DBStamp()
    << " ] to cover t=" << RunToDatabaseTimestamp(run)
    << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp()
    << "]";

  std::vector<unsigned int> mac5List;
  if (int res = febTable.GetChannelList(mac5List); res != 0) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "GetChannelList() returned " << res << " on run " << run << " query in " << fFEBTableName << "\n";
  }

  if (mac5List.empty()) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "Got an empty channel list for run " << run << " in " << fFEBTableName << "\n";
  }

  for(auto mac5 : mac5List)
    {
      long moduleType = 0;
      ReadElement(febTable, mac5, "type", moduleType);
      fCRTFEBCalibrationData[mac5].moduleType = static_cast<int>(moduleType);

      double t0CableLengthOffset = 0;
      ReadElement(febTable, mac5, "t0_timing_offset_cable_length", t0CableLengthOffset);
      fCRTFEBCalibrationData[mac5].t0CableLengthOffset = t0CableLengthOffset;

      double t0CalibratedOffset = 0;
      ReadElement(febTable, mac5, "t0_timing_offset_calibrated", t0CalibratedOffset);
      fCRTFEBCalibrationData[mac5].t0CalibratedOffset = t0CalibratedOffset;

      double t1CableLengthOffset = 0;
      ReadElement(febTable, mac5, "t1_timing_offset_cable_length", t1CableLengthOffset);
      fCRTFEBCalibrationData[mac5].t1CableLengthOffset = t1CableLengthOffset;

      double t1CalibratedOffset = 0;
      ReadElement(febTable, mac5, "t1_timing_offset_calibrated", t1CalibratedOffset);
      fCRTFEBCalibrationData[mac5].t1CalibratedOffset = t1CalibratedOffset;
    }
}

void sbndDB::CRTCalibrationDatabaseProvider::ReadCRTChannelCalibration(uint32_t run)
{
  lariov::DBFolder channelTable(fChannelTableName, "", "", fDatabaseTag, true, false);

  channelTable.UpdateData(fDatabaseTimeStamp);
  bool ret = channelTable.UpdateData(fDatabaseTimeStamp);

  mf::LogDebug("CRTCalibrationDatabaseProvider") << fChannelTableName << (ret ? "" : " not")
                                                 << " updated for run " << run;

  mf::LogTrace("CRTCalibrationDatabaseProvider")
    << "Fetched IoV [ " << channelTable.CachedStart().DBStamp() << " -> " << channelTable.CachedEnd().DBStamp()
    << " ] to cover t=" << RunToDatabaseTimestamp(run)
    << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp()
    << "]";

  std::vector<unsigned int> channelList;
  if (int res = channelTable.GetChannelList(channelList); res != 0) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "GetChannelList() returned " << res << " on run " << run << " query in " << fChannelTableName << "\n";
  }

  if (channelList.empty()) {
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "Got an empty channel list for run " << run << " in " << fChannelTableName << "\n";
  }

  for(auto ch : channelList)
    {
      long rawChannelNumber = 0;
      ReadElement(channelTable, ch, "raw_channel_number", rawChannelNumber);
      fCRTChannelCalibrationData[ch].rawChannelNumber = static_cast<int>(rawChannelNumber);

      long status = 0;
      ReadElement(channelTable, ch, "status", status);
      fCRTChannelCalibrationData[ch].status = static_cast<sbnd::crt::CRTChannelStatus>(status);

      long pedestal = 0;
      ReadElement(channelTable, ch, "pedestal", pedestal);
      fCRTChannelCalibrationData[ch].pedestal = static_cast<int>(pedestal);

      double gainFactor = 0;
      ReadElement(channelTable, ch, "gainfactor", gainFactor);
      fCRTChannelCalibrationData[ch].gainFactor = gainFactor;
    }
}

template <class T>
void sbndDB::CRTCalibrationDatabaseProvider::ReadElement(lariov::DBFolder &table, const int channel, const std::string &name, T &value)
{
  int error = table.GetNamedChannelData(channel, name, value);

  if(error)
    throw cet::exception("CRTCalibrationDatabaseProvider")
      << "Encountered error (code " << error
      << ") while trying to access '" << name << "' from table " << table.FolderName() << "\n";
}

// -----------------------------------------------------------------------------

void sbndDB::CRTCalibrationDatabaseProvider::readCRTCalibrationDatabase(const art::Run& run)
{
  fCRTFEBCalibrationData.clear();

  ReadCRTFEBCalibration(run.id().run());
  ReadCRTChannelCalibration(run.id().run());

  if(fVerbose)
    {
      mf::LogInfo("CRTCalibrationDatabaseProvider") << "Dump information from database " << std::endl;

      mf::LogVerbatim("CRTCalibrationDatabaseProvider")
        << "channel, type, t0_timing_offset_cable_length, t0_timing_offset_calibrated, t1_timing_offset_cable_length, t1_timing_offset_calibrated"
        << std::endl;

      for(auto const& [key, value] : fCRTFEBCalibrationData)
        mf::LogVerbatim("CRTCalibrationDatabaseProvider") << key << " " << value.moduleType << "," << value.t0CableLengthOffset << ","
                                                          << value.t0CalibratedOffset << "," << value.t1CableLengthOffset << ","
                                                          << value.t1CableLengthOffset << std::endl;

      mf::LogVerbatim("CRTCalibrationDatabaseProvider")
        << "channel, raw_channel_number, channel_status, pedestal, gain_factor"
        << std::endl;

      for(auto const& [key, value] : fCRTChannelCalibrationData)
        mf::LogVerbatim("CRTCalibrationDatabaseProvider") << key << " " << value.rawChannelNumber << "," << value.status << ","
                                                          << value.pedestal << "," << value.gainFactor << std::endl;
    }
}
