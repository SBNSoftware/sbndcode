set( ${product}_script_dir "${product}/${version}/scripts" )
file(GLOB xml_file *.xml )
install_fw(LIST ${xml_file})
#install( FILES ${xml_file} DESTINATION ${${product}_script_dir} )
