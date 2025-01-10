// ChannelNoiseService.h

// David Adams
// December 2015
//
// Interface for a service that adds noise to a TPC signal vector.

#ifndef ChannelNoiseService_H
#define ChannelNoiseService_H

#include <functional>
#include <vector>
#include <iostream>
#include "sbndcode/DetectorSim/Services/AdcTypes.h"
#include "fhiclcpp/ParameterSet.h"

namespace CLHEP { class HepRandomEngine; }
namespace detinfo { class DetectorClocksData; }

class ChannelNoiseService {

public:

  typedef unsigned int Channel;

  virtual ~ChannelNoiseService() =default;

  // Add noise to a signal vector sigs appropriate for channel chan.
  // Noise is added for all entries in the input vector.
  virtual int addNoise(detinfo::DetectorClocksData const&, Channel chan, AdcSignalVector& sigs) const =0;

  virtual void generateNoise(detinfo::DetectorClocksData const&){
    return;
  }

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;
  
  using EngineCreator = std::function<CLHEP::HepRandomEngine&(std::string const&, std::string const&)>;

  virtual void InitialiseProducerDeps(EngineCreator createEngine, fhicl::ParameterSet const& pset){
    return; 
  } 

};

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
DECLARE_ART_SERVICE_INTERFACE(ChannelNoiseService, LEGACY)
#endif

#endif
