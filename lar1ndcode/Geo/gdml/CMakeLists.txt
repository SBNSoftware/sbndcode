
set( ${product}_gdml_dir "${product}/${version}/gdml" )

# install gdml files
file(GLOB gdml_file *.gdml *.C *.xml README )
file(GLOB gdml_bin *.pl genmake )
install( FILES ${gdml_file} DESTINATION ${${product}_gdml_dir} )
install( PROGRAMS ${gdml_bin} DESTINATION ${${product}_gdml_dir} )

# special note: 
# Preferred solution is to not install GDMLSchema files 
# and instead, usd the version found in geant4 at GDML_SCHEMA_DIR,
# but this is not working. So revert to installing the
# schema files:
add_subdirectory(GDMLSchema)
