///////////////////////////////////////////////////////////////////////
///
/// Interface class between the calibration database and the PMT time corrections
///
///  A. Scarpelli
///
/// \mailto ascarpell@bnl.gov
///
////////////////////////////////////////////////////////////////////////
// Ported to SBND by Alejandro Sanchez-Castillo, Jan. 2025

#ifndef SBNDCCODE_DATABASEINTERFACE_PMTCALIBRATIONDATABASE_H
#define SBNDCCODE_DATABASEINTERFACE_PMTCALIBRATIONDATABASE_H


#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"

namespace sbndDB {

	class PMTCalibrationDatabase: lar::UncopiableClass
	{
		public: 
			
			virtual ~PMTCalibrationDatabase() noexcept = default;
			virtual int getBreakoutBox( unsigned int channelID ) const = 0;
			virtual int getCAENDigitizer( unsigned int channelID ) const = 0;
			virtual int getCAENDigitizerChannel( unsigned int channelID ) const = 0;
			virtual double getTotalTransitTime( unsigned int channelID ) const = 0;
			virtual double getSPEAmplitude( unsigned int channelID ) const = 0;
			virtual double getGaussFilterPower( unsigned int channelID ) const = 0;
			virtual double getGaussFilterWC( unsigned int channelID ) const = 0;
			virtual std::vector<double> getSER( unsigned int channelID ) const = 0;
	}; // end class

}// end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::PMTCalibrationDatabase, SHARED)

#endif