// SBNDThermalNoiseServiceInFreq.h
// Andrew Scarff (University of Sheffield)
// July 2019

// Based upon SPhaseChannelNoiseService.h created by Jingbo Wang (UC Davis) for ProtoDUNE.
// This Service runs the thermal noise model, but implemeted in frequency, not time.
//
// The default parameters set in: sbndcode/DetectorSim/Services/noiseservice_sbnd.fcl
//

#ifndef SBNDThermalNoiseServiceInFreq_H
#define SBNDThermalNoiseServiceInFreq_H

#include "sbndcode/DetectorSim/Services/ChannelNoiseService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandomEngine.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art_root_io/TFileService.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"

#include "TH1F.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

#include <sstream>
#include <vector>
#include <iostream>


class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class SBNDThermalNoiseServiceInFreq : public ChannelNoiseService {

public:

  // Ctor.
  SBNDThermalNoiseServiceInFreq(fhicl::ParameterSet const& pset);

  // Ctor.
  SBNDThermalNoiseServiceInFreq(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~SBNDThermalNoiseServiceInFreq();

  // Add noise to a signal array.
  int addNoise(detinfo::DetectorClocksData const& clockData,
               Channel chan, AdcSignalVector& sigs) const override;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const override;

private:
 
  // General parameters
  unsigned int            fNoiseArrayPoints; ///< number of points in randomly generated noise array
  int                     fRandomSeed;       ///< Seed for random number service. If absent or zero, use SeedSvc.
  int                     fLogLevel;         ///< Log message level: 0=quiet, 1=init only, 2+=every event
  std::map< double, int > fShapingTimeOrder;
  double                  fSampleRate;
  double                  fNoiseWidth;       ///< exponential noise width (kHz)
  double                  fNoiseRand;        ///< fraction of random "wiggle" in noise in freq. spectrum
  double                  fLowCutoff;        ///< low frequency filter cutoff (kHz)
  
  //Declare noise engines.
  CLHEP::HepRandomEngine* m_pran;
  CLHEP::HepRandomEngine* fNoiseEngine;

  // Function to allow use of noise engine in ChannelNoiseService setup.
  void InitialiseProducerDeps(EngineCreator createEngine, fhicl::ParameterSet const& pset) override{
    
    CLHEP::HepRandomEngine& NoiseEngine(art::ServiceHandle<rndm::NuRandomService>{}->registerAndSeedEngine(
                                          createEngine("HepJamesRandom","noise"), "HepJamesRandom","noise",pset,"Seed"));
    fNoiseEngine = &NoiseEngine;
    return; 
  } 

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(SBNDThermalNoiseServiceInFreq, ChannelNoiseService, LEGACY)

#endif
