install_fhicl()
add_subdirectory(pmt_variations)

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )
