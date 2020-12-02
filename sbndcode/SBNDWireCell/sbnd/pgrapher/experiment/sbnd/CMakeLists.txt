# Wire-Cell configuration files, required at runtime.
FILE(GLOB cfg_file *.jsonnet *.fcl )
message(STATUS "Configuration files ${cfg_file}")
install_wp(LIST ${cfg_file} SUBDIRNAME pgrapher/experiment/sbnd)
