add_subdirectory(flux_configb)
add_subdirectory(flux_configd)
add_subdirectory(flux_confige)
add_subdirectory(flux_configf)
add_subdirectory(flux_configh)
add_subdirectory(default)

install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )
