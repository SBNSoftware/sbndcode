install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(recomb_variations)
add_subdirectory(crt_filter)
add_subdirectory(optical_sim)

add_subdirectory(legacy)
