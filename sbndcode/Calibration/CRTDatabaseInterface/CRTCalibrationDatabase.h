#ifndef SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASE_H
#define SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASE_H

#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace sbndDB {

  class CRTCalibrationDatabase : lar::UncopiableClass {
  public:
    virtual ~CRTCalibrationDatabase() noexcept = default;
    virtual int getModuleType(unsigned int febMAC5) const = 0;
  }; // end class

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::CRTCalibrationDatabase, SHARED)

#endif
