#include "SBND_MCP.hh"
//#include "art/Framework/Core/EDProducer.h"                                                                                                   
//#include "art/Framework/Core/ModuleMacros.h"                                                                                                 
//#include "art/Framework/Principal/Event.h"                                                                                                   
//#include "art/Framework/Principal/Handle.h"                                                                                                  
//#include "art/Framework/Principal/Run.h"                                                                                                     
//#include "art/Framework/Principal/SubRun.h"                                                                                                  
//#include "canvas/Utilities/InputTag.h"                                                                                                       
//#include "fhiclcpp/ParameterSet.h"                                                                                                           
////#include "messagefacility/MessageLogger/MessageLogger.h"                                                                                   
//#include "canvas/Utilities/Exception.h"                                                                                                      
//#include "canvas/Utilities/InputTag.h"                                                                                                       
//#include "art_root_io/TFileService.h"                                                                                                        
//#include "art/Framework/Services/Registry/ServiceRegistry.h"                                                                                 

#include <CLHEP/Units/SystemOfUnits.h>
#include <Geant4/G4EmStandardPhysics.hh>
#include <Geant4/G4EmExtraPhysics.hh>
#include <Geant4/G4DecayPhysics.hh>
#include <Geant4/G4HadronElasticPhysics.hh>
#include <Geant4/G4HadronPhysicsQGSP_BERT.hh>
#include <Geant4/G4StoppingPhysics.hh>
#include <Geant4/G4IonPhysics.hh>
#include <Geant4/G4BuilderType.hh>
#include <Geant4/G4PhysicsConstructorFactory.hh>
#include <G4ParallelWorldScoringProcess.hh>
#include <fhiclcpp/ParameterSet.h>
#include <art/Framework/Services/Registry/ServiceHandle.h>

#include "Geant4/G4PhysListStamper.hh" // defines macro for factory registration                                                               
#include "McpIonisationProcess.h"
#include "/exp/sbnd/app/users/sc5303/mCP/larsoft/larsoft_v09_66_00_e20prof/srcs/sbndcode/sbndcode/LArG4/MillichargedParticleService_service.cc"

class MillichargedParticle : public G4ParticleDefinition {
  /// The single instance of this particle definition                                                                                          
  static MillichargedParticle *instance;

  MillichargedParticle() = delete;

public:
  ///                                                                                                                                          
  static MillichargedParticle *MillichargedParticleDefinition(double massFraction, double charge) {
    // short circuit if                                                                                                                        
    if (instance) return instance;

    auto *muon = G4MuonMinus::Definition();
    auto *anInstance = new G4ParticleDefinition(
						"millicharged",
						muon->GetPDGMass() * massFraction,
						0,
						-charge * CLHEP::eplus,
						0,
						0,
						0,
						0,
						0,
						0,
						// todo make empty string or "millicharged"                                                                                            
						muon->GetParticleType(),
						0,
						0,
						5900052,
						true,
						0,
        nullptr
						);

    instance = reinterpret_cast<MillichargedParticle *>(anInstance);
    return instance;
  }
};

MillichargedParticle *MillichargedParticle::instance = nullptr;

class MillichargedPhysics : public G4VPhysicsConstructor {
  double massFraction, charge;


public:
  MillichargedPhysics() : G4VPhysicsConstructor("mcp_physics") {
    std::cout << "DEFAULT CONSTRUCTED MILLICHARGED PHYSICS" << std::endl;
  }

  MillichargedPhysics(double massFraction, double charge) :
    G4VPhysicsConstructor("mcp_physics"), massFraction(massFraction), charge(charge) {}

  void ConstructParticle() override {
    std::cout << "CONSTRUCTING PARTICLE" << std::endl;
    auto *mcp = MillichargedParticle::MillichargedParticleDefinition(massFraction, charge);
    mcp->DumpTable();
    G4MuonMinus::MuonMinusDefinition();
    G4Gamma::Definition();
    G4Electron::Definition();
    G4Positron::Definition();
  }

  void ConstructProcess() override {
    std::cout << "CONSTRUCTING PROCESS" << std::endl;

    G4PhysicsListHelper *listHelper = G4PhysicsListHelper::GetPhysicsListHelper();

    auto *table = G4ParticleTable::GetParticleTable();
    auto *iterator = table->GetIterator();
    iterator->reset();
    while ((*iterator)()) {
      auto *particle = iterator->value();
      particle->DumpTable();

      listHelper->RegisterProcess(new McpIonisationProcess, particle);

      auto *scoringProcess = new G4ParallelWorldScoringProcess("LArVoxelReadoutScoringProcess");
      scoringProcess->SetParallelWorld("LArVoxelReadoutGeometry");
      scoringProcess->SetProcessType(G4ProcessType::fParallel);
      scoringProcess->SetProcessSubType(491);
      std::cout << "LArVoxelParallelWorldScoringProcess->GetProcessType = "
                << scoringProcess->GetProcessType() << std::endl;
      std::cout << "LArVoxelParallelWorldScoringProcess->GetProcessSubType = "
                << scoringProcess->GetProcessSubType() << std::endl;
      listHelper->RegisterProcess(scoringProcess, particle);
    }
  }
};

G4_DECLARE_PHYSCONSTR_FACTORY(MillichargedPhysics);

G4_DECLARE_PHYSLIST_FACTORY(SBND_MCP);

SBND_MCP::SBND_MCP(int verbose) {
  if (verbose > 1) {
    G4cout << "<<< Geant4 Physics List simulation engine: SBND_MCP"
           << G4endl
           << G4endl;
  }

  defaultCutValue = 0.7 * CLHEP::mm;
  SetVerboseLevel(verbose);

  art::ServiceHandle<MillichargedParticleService> service;
  std::cout << "service->massFraction = " << service->massFraction << std::endl;
  std::cout << "service->charge = " << service->charge << std::endl;
  RegisterPhysics(new MillichargedPhysics(service->massFraction, service->charge));

}
