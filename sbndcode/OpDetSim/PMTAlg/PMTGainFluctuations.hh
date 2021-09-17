///////////////////////////////////////////////////////////////////////
///
/// Interface class for PMTGainFluctuations1Dynode tool
///
////////////////////////////////////////////////////////////////////////

#ifndef SBND_PMTGainFluctuations_H
#define SBND_PMTGainFluctuations_H



namespace opdet {
  class PMTGainFluctuations;
}

//Base class
class opdet::PMTGainFluctuations {
public:
  //Constructor
  virtual ~PMTGainFluctuations() noexcept = default;

  //Returns fluctuated factor for SPR
  virtual double GainFluctuation(unsigned int npe, CLHEP::HepRandomEngine* eng) = 0;
};

#endif
