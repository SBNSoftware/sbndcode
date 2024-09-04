add_subdirectory(other_flux_config)
add_subdirectory(to_be_removed)

install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )
