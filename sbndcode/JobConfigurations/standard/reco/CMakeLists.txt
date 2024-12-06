add_subdirectory(deprecated)
add_subdirectory(config)

install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )
