add_subdirectory(PhysicsLists)

set(
  MODULE_LIBRARIES
    art::Utilities
    cetlib_except::cetlib_except
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::Geom
    ROOT::MathCore
    CLHEP::CLHEP
    Geant4::G4digits_hits
    Geant4::G4event
    Geant4::G4geometry
    Geant4::G4global
    Geant4::G4graphics_reps
    Geant4::G4materials
    Geant4::G4intercoms
    Geant4::G4particles
    Geant4::G4physicslists
    Geant4::G4processes
    Geant4::G4track
    Geant4::G4run
    sbndcode_LArG4_PhysicsLists
)

cet_build_plugin(PhysListLoader art::module SOURCE PhysListLoader_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_source()
install_fhicl()

