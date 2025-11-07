#ifndef SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASEPROVIDER_H
#define SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASEPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// Local
#include "sbndcode/Calibration/CRTDatabaseInterface/CRTCalibrationDatabase.h"

// C/C++ standard libraries
#include <string>
#include <map>
#include <stdint.h>


namespace sbndDB
{
  class CRTCalibrationDatabaseProvider;
}

class sbndDB::CRTCalibrationDatabaseProvider : public CRTCalibrationDatabase
{
public:
  CRTCalibrationDatabaseProvider(const fhicl::ParameterSet& pset);

  void readCRTCalibrationDatabase(const art::Run& run);

  int getModuleType(unsigned int febMAC5) const override {
    return getCRTFEBCalibrationOrDefault(febMAC5).moduleType;
  };

  double getT0CableLengthOffset(unsigned int febMAC5) const override {
    return getCRTFEBCalibrationOrDefault(febMAC5).t0CableLengthOffset;
  };

  double getT0CalibratedOffset(unsigned int febMAC5) const override {
    return getCRTFEBCalibrationOrDefault(febMAC5).t0CalibratedOffset;
  };

  double getT1CableLengthOffset(unsigned int febMAC5) const override {
    return getCRTFEBCalibrationOrDefault(febMAC5).t1CableLengthOffset;
  };

  double getT1CalibratedOffset(unsigned int febMAC5) const override {
    return getCRTFEBCalibrationOrDefault(febMAC5).t1CalibratedOffset;
  };

  int getRawChannelNumber(unsigned int channel) const override {
    return getCRTChannelCalibrationOrDefault(channel).rawChannelNumber;
  };

  sbnd::crt::CRTChannelStatus getChannelStatus(unsigned int channel) const override {
    return getCRTChannelCalibrationOrDefault(channel).status;
  };

  int getPedestal(unsigned int channel) const override {
    return getCRTChannelCalibrationOrDefault(channel).pedestal;
  };

  double getGainFactor(unsigned int channel) const override {
    return getCRTChannelCalibrationOrDefault(channel).gainFactor;
  };

private:
        
  bool fVerbose;

  std::string fDatabaseTag;  
  long        fDatabaseTimeStamp;
  std::string fFEBTableName;
  std::string fChannelTableName;

  struct CRTFEBCalibrationDB { 
    size_t moduleType          = 0;
    double t0CableLengthOffset = 0.;
    double t0CalibratedOffset  = 0.;
    double t1CableLengthOffset = 0.;
    double t1CalibratedOffset  = 0.;
  };

  struct CRTChannelCalibrationDB {
    size_t rawChannelNumber            = 0;
    sbnd::crt::CRTChannelStatus status = sbnd::crt::CRTChannelStatus::kGoodChannel;
    size_t pedestal                    = 0;
    double gainFactor                  = 0.;
  };
            
  const CRTFEBCalibrationDB CRTFEBDefaults         = {0, 0., 0., 0., 0.};
  const CRTChannelCalibrationDB CRTChannelDefaults = {0, sbnd::crt::CRTChannelStatus::kGoodChannel, 0, 0.};

  std::map<unsigned int, CRTFEBCalibrationDB> fCRTFEBCalibrationData;
  std::map<unsigned int, CRTChannelCalibrationDB> fCRTChannelCalibrationData;

  CRTFEBCalibrationDB const& getCRTFEBCalibrationOrDefault(unsigned int febMAC5) const
  {
    auto const it = fCRTFEBCalibrationData.find(febMAC5);
    return (it == fCRTFEBCalibrationData.end()) ? CRTFEBDefaults : it->second;
  }

  CRTChannelCalibrationDB const& getCRTChannelCalibrationOrDefault(unsigned int channel) const
  {
    auto const it = fCRTChannelCalibrationData.find(channel);
    return (it == fCRTChannelCalibrationData.end()) ? CRTChannelDefaults : it->second;
  }

  uint64_t RunToDatabaseTimestamp(uint32_t run) const;

  void ReadCRTFEBCalibration(uint32_t run);
  void ReadCRTChannelCalibration(uint32_t run);

  template <class T>
  void ReadElement(lariov::DBFolder &table, const int channel, const std::string &name, T &value);
};

#endif 
