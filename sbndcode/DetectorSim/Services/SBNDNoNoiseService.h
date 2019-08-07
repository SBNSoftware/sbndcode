// SBNDNoNoiseService.h
// Andrew Scarff (University of Sheffield)
// July 2019

// Based upon SPhaseChannelNoiseService.h created by Jingbo Wang (UC Davis) for ProtoDUNE.9

// Implementation of a general TPC channel noise model with:
// (1) white noise
// (2) Inherent Gaussian noise in frequency
// (3) MicroBooNE noise in frequency
// (4) Coherent noise (exponential + Gaussian) in frequency 
//     (Note a: phase at each frequency bin is randamized at the moment. Will be updated soon
//      Note b: Currently, consecutive offline channels (configurable) are grouped together and 
//              the same coherent noise waveform is assigned to channels within the same group. )
//
// The default parameters are obtained from the ProtoDUNE-SP data (run 4096)
// fcl file: sbndcode/DetectorSim/Services/SBND_detsim_data_driven_noise.fcl
//

#ifndef SBNDNoNoiseService_H
#define SBNDNoNoiseService_H

#include "sbndcode/DetectorSim/Services/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class SBNDNoNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  SBNDNoNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  SBNDNoNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~SBNDNoNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
 

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(SBNDNoNoiseService, ChannelNoiseService, LEGACY)

#endif
