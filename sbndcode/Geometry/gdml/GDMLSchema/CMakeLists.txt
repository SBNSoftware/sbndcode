# install schema files
file(GLOB gdml_schema *.xsd)
install( FILES ${gdml_schema} DESTINATION ${${product}_gdml_dir}/GDMLSchema )
