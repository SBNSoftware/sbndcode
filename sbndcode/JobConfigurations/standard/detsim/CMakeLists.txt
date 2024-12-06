install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(detector_variations)
add_subdirectory(legacy)

