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

private:
        
  bool fVerbose = false;

  std::string fDatabaseTag;  
  long        fDatabaseTimeStamp;
  std::string fFEBTableName;
  std::string fChannelTableName;

  struct CRTFEBCalibrationDB { 
    size_t moduleType=0;
  };
            
  const CRTFEBCalibrationDB CRTFEBDefaults = {0};

  std::map<unsigned int, CRTFEBCalibrationDB> fCRTFEBCalibrationData;
        
  CRTFEBCalibrationDB const& getCRTFEBCalibrationOrDefault (unsigned int febMAC5) const
  {
    auto const it = fCRTFEBCalibrationData.find(febMAC5);
    return (it == fCRTFEBCalibrationData.end()) ? CRTFEBDefaults: it->second;
  }

  uint64_t RunToDatabaseTimestamp(uint32_t run) const;

  void ReadCRTFEBCalibration(uint32_t run);
};

#endif 
