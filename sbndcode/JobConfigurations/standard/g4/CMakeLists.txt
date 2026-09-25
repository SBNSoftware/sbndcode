install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(recomb_variations)
add_subdirectory(crt_filter)
add_subdirectory(optical_sim)

# LegacyLArG4 is not supported when using geant4 4.11
##add_subdirectory(legacy)
