install_fhicl()
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

file(GLOB xml_file *.xml )
install_fw(LIST ${xml_file})
