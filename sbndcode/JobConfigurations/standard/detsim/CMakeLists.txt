install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(detector_variations)
# LegacyLArG4 is not supported when using geant4 4.11
##add_subdirectory(legacy)

