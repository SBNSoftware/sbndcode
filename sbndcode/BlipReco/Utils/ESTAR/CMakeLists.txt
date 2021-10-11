file(GLOB root_file *.root )
install_fw(LIST ${root_file} )
#install( FILES ${root_file} DESTINATION "${product}/${version}/scripts" )
#install_fw( LIST ESTAREnergyLookupCurve.root )
