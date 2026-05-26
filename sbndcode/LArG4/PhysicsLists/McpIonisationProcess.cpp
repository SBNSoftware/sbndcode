#include <Geant4/G4BraggModel.hh>
#include <Geant4/G4ICRU73QOModel.hh>
#include <Geant4/G4IonFluctuations.hh>
#include <Geant4/G4BetheBlochModel.hh>
#include "McpIonisationProcess.h"

McpIonisationProcess::McpIonisationProcess()
  : G4VEnergyLossProcess("mcpIonisation"),
    initialized(false) {
  SetProcessSubType(fIonisation);
}

G4bool McpIonisationProcess::IsApplicable(const G4ParticleDefinition &p) {
  return p.GetPDGCharge() != 0.0;
}

void McpIonisationProcess::InitialiseEnergyLossProcess(const G4ParticleDefinition *particle,
                                                       const G4ParticleDefinition *baseParticle) {
  std::cout << "initialized = " << initialized << std::endl;
  if (particle != nullptr) std::cout << "particle->GetParticleName() = " << particle->GetParticleName() << std::endl;
  if (baseParticle != nullptr)
    std::cout << "baseParticle.GetParticleName() = " << baseParticle->GetParticleName() << std::endl;
  if (initialized) return;
  initialized = true;

  G4VEmModel *emModel = EmModel(0);
  std::cout << "emModel = " << emModel << std::endl;
  if (!emModel) {
    if (particle->GetPDGCharge() > 0.0) SetEmModel(new G4BraggModel(particle));
    else SetEmModel(new G4ICRU73QOModel(particle));
  }
  //  EmModel(0)->SetLowEnergyLimit()                                                                                          
  AddEmModel(1, EmModel(0), new G4IonFluctuations());

  // high energy fluctuation model                                                                                           
  if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }

  // moderate energy model                                                                                                   
  if (!EmModel(1)) {
    auto *betheBloch = new G4BetheBlochModel(particle);
    G4DataVector shh;
    betheBloch->Initialise(particle, shh);
    SetEmModel(betheBloch);
  }
  //  EmModel(1)->SetLowEnergyLimit(elow);                                                                                     
  //  EmModel(1)->SetHighEnergyLimit(ehigh);                                                                                   
  AddEmModel(2, EmModel(1), FluctModel());


  for (int i = 0; i < 2; ++i) {
    G4VEmModel *model = EmModel(i);
    std::cout << "model = " << model->GetName() << std::endl;
  }
}
