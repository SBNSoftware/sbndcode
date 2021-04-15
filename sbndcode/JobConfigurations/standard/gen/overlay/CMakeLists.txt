add_subdirectory(other_flux_config)

install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )
