///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceT1053.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author H. Greenlee 
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response.
/// Col3DCorrection - 3D path length correction for collection plane.
/// Ind3DCorrection - 3D path length correction for induction plane.
/// ColFieldRespAmp - Collection field response amplitude.
/// IndFieldRespAmp - Induction field response amplitude.
/// ShapeTimeConst  - Time constants for exponential shaping.
/// ColFilter       - Root parameterized collection plane filter function.
/// ColFilterParams - Collection filter function parameters.
/// IndFilter       - Root parameterized induction plane filter function.
/// IndFilterParams - Induction filter function parameters.
///
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGSERVICELARIAT_H
#define SIGNALSHAPINGSERVICELARIAT_H

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"


namespace util {
  class SignalShapingServiceT1053 {
  public:

    // Constructor, destructor.

    SignalShapingServiceT1053(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~SignalShapingServiceT1053();

    // Update configuration parameters.

    void reconfigure(const fhicl::ParameterSet& pset);

    // Accessors.

    const util::SignalShaping& SignalShaping(unsigned int channel) const;

    // Do convolution calcution (for simulation).

    template <class T> void Convolute(unsigned int channel, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).

    template <class T> void Deconvolute(unsigned int channel, std::vector<T>& func) const;

  private:

    // Private configuration methods.

    // Post-constructor initialization.

    void init() const{const_cast<SignalShapingServiceT1053*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireT1053.

    void SetFieldResponse();
    void SetElectResponse();

    // Calculate filter functions.

    void SetFilters();

    // Attributes.

    bool fInit;               ///< Initialization flag.

    // Fcl parameters.
    double fADCTicksPerPCAtLowestASICGainSetting; ///< Pulse area (in ADC*ticks) for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC

    double fASICGainInMVPerFC;                  ///< Cold electronics ASIC gain setting in mV/fC

    int fNFieldBins;         			///< number of bins for field response
    double fCol3DCorrection; 			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fInd3DCorrection;  			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fColFieldRespAmp;  			///< amplitude of response to field 
    double fIndUFieldRespAmp;  			///< amplitude of response to field in U plane
    double fIndVFieldRespAmp;  			///< amplitude of response to field in V plane
    std::vector<double> fShapeTimeConst;  	///< time constants for exponential shaping
    TF1* fColFilterFunc;      			///< Parameterized collection filter function.
    TF1* fIndUFilterFunc;      			///< Parameterized induction filter function for U plane.
    TF1* fIndVFilterFunc;      			///< Parameterized induction filter function for V plane

    
    bool fUseFunctionFieldShape;   		///< Flag that allows to use a parameterized field response instead of the hardcoded version
    bool fUseSimpleFieldShape;                 ///< Flag that turns on new field response shapes
    bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
    TF1* fColFieldFunc;      			///< Parameterized collection field shape function.
    TF1* fIndUFieldFunc;      			///< Parameterized induction field shape function for U plane.
    TF1* fIndVFieldFunc;      			///< Parameterized induction field shape function for V plane.
    
    TH1D *fFilterHist[3];    			///< Histogram used to hold the collection filter, hardcoded for the time being
    
    // Following attributes hold the convolution and deconvolution kernels

    util::SignalShaping fIndUSignalShaping;
    util::SignalShaping fIndVSignalShaping;
    util::SignalShaping fColSignalShaping;

    // Field response.

    std::vector<double> fIndUFieldResponse;
    std::vector<double> fIndVFieldResponse;
    std::vector<double> fColFieldResponse;

    // Electronics response.

    std::vector<double> fElectResponse;

    // Filters.

    std::vector<TComplex> fIndUFilter;
    std::vector<TComplex> fIndVFilter;
    std::vector<TComplex> fColFilter;
  };
}
//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceT1053::Convolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Convolute(func);
}


//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceT1053::Deconvolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Deconvolute(func);
}

DECLARE_ART_SERVICE(util::SignalShapingServiceT1053, LEGACY)
#endif
