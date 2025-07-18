#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace filt {

class SimEnergyDepFakeTriggerFilter : public art::EDFilter {
  public:
  explicit SimEnergyDepFakeTriggerFilter(fhicl::ParameterSet const& pset);
  virtual bool filter(art::Event& e) override;

  private:
  const double fBeamTimeMin;   // Minimum time of beam window [us]
  const double fBeamTimeMax;   // Maximum time of beam window [us]
  const double fEnergyDeposit; // Minimum energy deposit in TPC for trigger [MeV]

  const std::string fSimEnergyDepModuleName;
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
  const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>&
      energyDeps(e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fSimEnergyDepModuleName));

  double energy(0);

  for (const sim::SimEnergyDeposit& energyDep : *energyDeps) {
    // Check particle time is within the beam time
    const double time(energyDep.Time() * 1e-3); // [ns] -> [us]
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
