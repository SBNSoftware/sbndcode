
# The xml files here are used to initialize pandora
# They need to be found via FW_SEARCH_PATH

# install xml files
set( ${product}_script_dir "${product}/${version}/scripts" )
file(GLOB xml_file *.xml )
install( FILES ${xml_file} DESTINATION ${${product}_script_dir} )

