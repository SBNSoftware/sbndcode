// SBNDNoNoiseService.cxx
// Andrew Scarff
// July 2019
// Based upon SPhaseChannelNoiseService.cxx developed by Jingbo Wang for ProtoDUNE.

#include "sbndcode/DetectorSim/Services/SBNDNoNoiseService.h"
#include <sstream>

using std::cout;
using std::ostream;
using std::endl;
using std::string;

//**********************************************************************

SBNDNoNoiseService::
SBNDNoNoiseService(fhicl::ParameterSet const& pset) { 
  const string myname = "SBNDNoNoiseService::ctor: ";  
}

//**********************************************************************

SBNDNoNoiseService::
SBNDNoNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: SBNDNoNoiseService(pset) { }

//**********************************************************************

SBNDNoNoiseService::~SBNDNoNoiseService() {
  const string myname = "SBNDNoNoiseService::dtor: ";
}

//**********************************************************************

int SBNDNoNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
  return 0;
}

//**********************************************************************

ostream& SBNDNoNoiseService::print(ostream& out, string prefix) const {
  out << prefix << "SBNDNoNoiseService: " << endl;
  
  return out;
}



//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(SBNDNoNoiseService, ChannelNoiseService)

//**********************************************************************
