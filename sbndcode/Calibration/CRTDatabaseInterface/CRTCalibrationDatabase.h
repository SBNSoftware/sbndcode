#ifndef SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASE_H
#define SBNDCODE_DATABASEINTERFACE_CRTCALIBRATIONDATABASE_H

#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include "sbnobj/SBND/CRT/CRTEnums.hh"

namespace sbndDB {

  class CRTCalibrationDatabase : lar::UncopiableClass {
  public:
    virtual ~CRTCalibrationDatabase() noexcept = default;
    virtual int getModuleType(unsigned int febMAC5) const = 0;
    virtual double getT0CableLengthOffset(unsigned int febMAC5) const = 0;
    virtual double getT0CalibratedOffset(unsigned int febMAC5) const = 0;
    virtual double getT1CableLengthOffset(unsigned int febMAC5) const = 0;
    virtual double getT1CalibratedOffset(unsigned int febMAC5) const = 0;
    virtual int getRawChannelNumber(unsigned int channel) const = 0;
    virtual sbnd::crt::CRTChannelStatus getChannelStatus(unsigned int channel) const = 0;
    virtual int getPedestal(unsigned int channel) const = 0;
    virtual double getGainFactor(unsigned int channel) const = 0;
  }; // end class

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::CRTCalibrationDatabase, SHARED)

#endif
