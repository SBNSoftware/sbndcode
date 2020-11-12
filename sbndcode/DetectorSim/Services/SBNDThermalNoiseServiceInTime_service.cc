// SBNDThermalNoiseServiceInTime.cxx
// Andrew Scarff
// July 2019
// Based upon SPhaseChannelNoiseService.cxx developed by Jingbo Wang for ProtoDUNE.

#include "sbndcode/DetectorSim/Services/SBNDThermalNoiseServiceInTime.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;
using CLHEP::HepRandomEngine;

//**********************************************************************

SBNDThermalNoiseServiceInTime::
SBNDThermalNoiseServiceInTime(fhicl::ParameterSet const& pset): fRandomSeed(0), fLogLevel(1),
								m_pran(nullptr), fNoiseEngine(nullptr)
 {
   
  const string myname = "SBNDThermalNoiseServiceInTime::ctor: ";
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints",1);
  bool haveSeed      = pset.get_if_present<int>("RandomSeed", fRandomSeed);

  fShapingTimeOrder = { {0.5, 0 }, {1.0, 1}, {2.0, 2}, {3.0, 3} };

  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  int seed = fRandomSeed;
  
  string rname = "SBNDThermalNoiseServiceInTime";
  if ( haveSeed ) {
    if ( fLogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(seed);
  } else {
    if ( fLogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( fLogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( fLogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;
  if ( fLogLevel > 1 ) print() << endl;
}

//**********************************************************************

SBNDThermalNoiseServiceInTime::
SBNDThermalNoiseServiceInTime(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: SBNDThermalNoiseServiceInTime(pset) { }

//**********************************************************************

SBNDThermalNoiseServiceInTime::~SBNDThermalNoiseServiceInTime() {
  const string myname = "SBNDThermalNoiseServiceInTime::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int SBNDThermalNoiseServiceInTime::addNoise(detinfo::DetectorClocksData const&,
                                            Channel chan, AdcSignalVector& sigs) const {

  //Get services.
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::SignalShapingServiceSBND> sss;
  
  //Generate Noise:
  size_t view = (size_t)geo->View(chan);
  
  double noise_factor;
  auto tempNoiseVec = sss->GetNoiseFactVec();
  double shapingTime = 2.0; //sss->GetShapingTime(chan);
  double asicGain = sss->GetASICGain(chan);
  
  if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() ) {
    noise_factor = tempNoiseVec[view].at( fShapingTimeOrder.find( shapingTime )->second );
    noise_factor *= asicGain/4.7;
  }
  else {
    throw cet::exception("SBNDThermalNoiseServiceInTime_service")
      << "\033[93m"
      << "Shaping Time recieved from signalshapingservices_sbnd.fcl is not one of the allowed values"
      << std::endl
      << "Allowed values: 0.5, 1.0, 2.0, 3.0 us"
      << "\033[00m"
      << std::endl;
  }

  CLHEP::RandGaussQ rGauss(*fNoiseEngine, 0.0, noise_factor);
    

  //In this case fNoiseFact is a value in ADC counts
  //It is going to be the Noise RMS
  //loop over all bins in "noise" vector
  //and insert random noise value
  
  for (unsigned int i = 0; i < sigs.size(); i++){
    sigs.at(i) = rGauss.fire();
  }

  return 0;
}


//**********************************************************************

ostream& SBNDThermalNoiseServiceInTime::print(ostream& out, string prefix) const {
  out << prefix << "SBNDThermalNoiseServiceInTime: " << endl;
  
  out << prefix << "          LogLevel: " <<  fLogLevel << endl;
  out << prefix << "        RandomSeed: " <<  fRandomSeed << endl;
  out << prefix << "  NoiseArrayPoints: " << fNoiseArrayPoints << endl;
  
  return out;
}



//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(SBNDThermalNoiseServiceInTime, ChannelNoiseService)

//**********************************************************************
