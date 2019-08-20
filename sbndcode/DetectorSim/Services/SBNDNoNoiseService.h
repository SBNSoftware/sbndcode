// SBNDNoNoiseService.h
// Andrew Scarff (University of Sheffield)
// July 2019

// Based upon SPhaseChannelNoiseService.h created by Jingbo Wang (UC Davis) for ProtoDUNE.

// This file is to be used if you don't want any noise added to the channels.
// fcl file: sbndcode/DetectorSim/Services/noiseservices_sbnd.fcl
//

#ifndef SBNDNoNoiseService_H
#define SBNDNoNoiseService_H

#include "sbndcode/DetectorSim/Services/ChannelNoiseService.h"

#include <vector>
#include <iostream>
#include <sstream>

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
