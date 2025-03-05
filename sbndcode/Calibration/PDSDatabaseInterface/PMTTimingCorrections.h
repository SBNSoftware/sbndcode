///////////////////////////////////////////////////////////////////////
///
/// \file   icaruscode/Timing/PMTTimingCorrections.h
///
/// \brief  Interface class between the calibration database and the PMT time corrections
///
/// \author A. Scarpelli
///
/// \mailto ascarpell@bnl.gov
///
////////////////////////////////////////////////////////////////////////

#ifndef ICARUSCODE_TIMING_PMTTIMINGCORRECTIONS_H
#define ICARUSCODE_TIMING_PMTTIMINGCORRECTIONS_H


#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"

namespace sbndDB {

	class PMTTimingCorrections: lar::UncopiableClass
	{
		public: 
			
			virtual ~PMTTimingCorrections() noexcept = default;
			
			virtual double getTriggerCableDelay( unsigned int channelID ) const = 0;
			
			virtual double getResetCableDelay( unsigned int channelID ) const = 0;

			virtual double getLaserCorrections( unsigned int channelID ) const = 0;

			virtual double getCosmicsCorrections( unsigned int channelID ) const = 0;

	}; // end class

}// end of namespace

DECLARE_ART_SERVICE_INTERFACE(sbndDB::PMTTimingCorrections, SHARED)

#endif