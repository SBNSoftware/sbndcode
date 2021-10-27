#include <algorithm>

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace filt {

class SimEnergyDepFakeTriggerFilter : public art::EDFilter {
  public:
  explicit SimEnergyDepFakeTriggerFilter(fhicl::ParameterSet const& pset);
  virtual bool filter(art::Event& e) override;

  private:
  double fBeamTimeMin;   // Minimum time of beam window [us]
  double fBeamTimeMax;   // Maximum time of beam window [us]
  double fEnergyDeposit; // Minimum energy deposit in TPC for trigger [MeV]

  std::string fSimEnergyDepModuleName;
};

SimEnergyDepFakeTriggerFilter::SimEnergyDepFakeTriggerFilter(fhicl::ParameterSet const& pset)
    : EDFilter(pset)
    , fBeamTimeMin(pset.get<double>("BeamTimeMin"))
    , fBeamTimeMax(pset.get<double>("BeamTimeMax"))
    , fEnergyDeposit(pset.get<double>("EnergyDeposit"))
    , fSimEnergyDepModuleName(pset.get<std::string>("SimEnergyDepModuleName"))
{
}

bool SimEnergyDepFakeTriggerFilter::filter(art::Event& e)
{
  art::Handle<std::vector<sim::SimEnergyDeposit>> energyDeps;
  e.getByLabel(fSimEnergyDepModuleName, energyDeps);

  double energy(0);

  for (const sim::SimEnergyDeposit& energyDep : *energyDeps) {

    // Check particle time is within the beam time
    const double time = energyDep.Time() * 1e-3; // [ns] -> [us]
    if (time < fBeamTimeMin || time > fBeamTimeMax)
      continue;

    // Add up the energy deposit inside the TPC
    energy += energyDep.Energy(); // [MeV]
  }

  // If the energy deposit within the beam time is greater than some limit then trigger the event
  return energy > fEnergyDeposit;
}

DEFINE_ART_MODULE(SimEnergyDepFakeTriggerFilter)

}
