add_subdirectory(PhysicsLists)

include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

art_make(
  MODULE_LIBRARIES
    ${ART_UTILITIES}
    cetlib_except
    ${ROOT_CORE}
    ${ROOT_HIST}
    ${ROOT_PHYSICS}
    ${ROOT_GEOM}
    ${ROOT_MATHCORE}
    ${CLHEP}
    ${G4DIGITS_HITS}
    ${G4EVENT}
    ${G4GEOMETRY}
    ${G4GLOBAL}
    ${G4GRAPHICS_REPS}
    ${G4MATERIALS}
    ${G4INTERCOMS}
    ${G4PARTICLES}
    ${G4PHYSICSLISTS}
    ${G4PROCESSES}
    ${G4TRACK}
    ${G4RUN}
    sbndcode_LArG4_PhysicsLists

  LIB_LIBRARIES
    ${ART_UTILITIES}
    cetlib_except
    ${ROOT_CORE}
    ${ROOT_HIST}
    ${ROOT_PHYSICS}
    ${ROOT_GEOM}
    ${ROOT_MATHCORE}
    ${CLHEP}
    ${G4DIGITS_HITS}
    ${G4EVENT}
    ${G4GEOMETRY}
    ${G4GLOBAL}
    ${G4GRAPHICS_REPS}
    ${G4MATERIALS}
    ${G4INTERCOMS}
    ${G4PARTICLES}
    ${G4PHYSICSLISTS}
    ${G4PROCESSES}
    ${G4TRACK}
    ${G4RUN}
    sbndcode_LArG4_PhysicsLists

  SERVICE_LIBRARIES
    ${ART_FRAMEWORK_CORE}
)

install_headers()
install_source()
install_fhicl()

