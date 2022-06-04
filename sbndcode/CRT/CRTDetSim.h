/**
 * \brief LArSoft plugin for SBND CRT detector simulation parameters
 *
 * \author Andy Mastbaum
 * \author Marco Del Tutto
 */

 ///////////////////////////////////////////////////////////////////////////////
/// Class: CRTDetSim
/// Module Type: producer
/// File: CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
///
/// Author: mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataalg/DetectorInfo/ElecClock.h"

#include "sbndcode/CRT/CRTSimulation/CRTDetSimAlg.h"

#include <string>

namespace sbnd {
namespace crt {

class CRTDetSim : public art::EDProducer {
public:
  explicit CRTDetSim(fhicl::ParameterSet const & p);

  CRTDetSim(CRTDetSim const &) = delete;
  CRTDetSim(CRTDetSim &&) = delete;
  CRTDetSim& operator = (CRTDetSim const &) = delete;
  CRTDetSim& operator = (CRTDetSim &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p) ;

  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

private:

  CLHEP::HepRandomEngine& fEngine; //!< Reference to art-managed random-number engine
  double fG4RefTime; //!< Stores the G4 reference time
  CRTDetSimAlg fDetAlg; //!< Instance of the CRT detector simulation algorithm
};

}  // namespace crt
}  // namespace sbnd
