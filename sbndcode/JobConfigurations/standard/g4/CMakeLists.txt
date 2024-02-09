install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(recomb_variations)
add_subdirectory(legacy)
add_subdirectory(2d_drift_sim)
