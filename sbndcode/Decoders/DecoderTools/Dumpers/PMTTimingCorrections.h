///////////////////////////////////////////////////////////////////////
///
/// \file   sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrections.h
///
/// \brief  Interface class between the calibration database and the PMT time corrections
///
/// \author Afroditi Papadopoulou
///
/// \mailto apapadopoulou@anl.gov
///
////////////////////////////////////////////////////////////////////////

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTTIMINGCORRECTIONS_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DUMPERS_PMTTIMINGCORRECTIONS_H


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
