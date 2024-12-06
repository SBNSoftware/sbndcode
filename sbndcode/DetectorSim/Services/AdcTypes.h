// AdcTypes.h

// David Adams
// December 2015
//
// Integer and floating type for Adc signals.

#ifndef AdcTypes_H
#define AdcTypes_H

#include <vector>
#include <utility>

typedef unsigned int AdcIndex;
typedef unsigned long AdcLongIndex;

typedef short AdcCount;
typedef std::vector<AdcCount> AdcCountVector;

typedef float AdcSignal;
typedef std::vector<AdcSignal> AdcSignalVector;
typedef std::vector<AdcSignalVector> AdcSignalVectorVector;

typedef float AdcPedestal;

typedef std::vector<bool> AdcFilterVector;

typedef short AdcFlag;
typedef std::vector<AdcFlag> AdcFlagVector;

const AdcFlag AdcGood         =  0; // ADC sample is fine
const AdcFlag AdcUnderflow    =  1; // ADC sample is underflow
const AdcFlag AdcOverflow     =  2; // ADC sample is overflow
const AdcFlag AdcStuck        =  8; // ADC sample is sticky unspecified
const AdcFlag AdcStuckOff     =  9; // ADC sample is sticky with low bits at 0
const AdcFlag AdcStuckOn      = 10; // ADC sample is sticky with low bits at 1
const AdcFlag AdcStuckPed     = 11; // ADC sample is sticky code near pedestal
const AdcFlag AdcStuckSig     = 12; // ADC sample is sticky code away from pedestal
const AdcFlag AdcMitigated    = 16; // ADC sample is mitigated unspecified
const AdcFlag AdcSetFixed     = 17; // ADC sample is mitigated with set to a fixed value
const AdcFlag AdcInterpolated = 18; // ADC sample is mitigated with interpolation
const AdcFlag AdcExtrapolated = 19; // ADC sample is mitigated with extrapolation

const AdcIndex AdcChannelStatusGood  = 0;
const AdcIndex AdcChannelStatusBad   = 1;
const AdcIndex AdcChannelStatusNoisy = 2;

typedef unsigned int AdcChannel;
typedef std::vector<AdcChannel> AdcChannelVector;

// A ROI includes all ticks from roi.first through roi.second.
typedef std::pair<AdcIndex, AdcIndex> AdcRoi;
typedef std::vector<AdcRoi> AdcRoiVector;

#endif
