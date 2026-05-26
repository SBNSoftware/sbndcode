#ifndef LARSIM_MCPIONISATIONPROCESS_H
#define LARSIM_MCPIONISATIONPROCESS_H

#include <Geant4/G4VEnergyLossProcess.hh>

class McpIonisationProcess : public G4VEnergyLossProcess {
  bool initialized;

 public:
  McpIonisationProcess();

  G4bool IsApplicable(const G4ParticleDefinition &p) override;

 protected:
  void InitialiseEnergyLossProcess(const G4ParticleDefinition *particle,
                                   const G4ParticleDefinition *baseParticle) override;
};


#endif //LARSIM_MCPIONISATIONPROCESS_H                                                                                       


