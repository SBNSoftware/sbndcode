add_subdirectory(sim)

file(GLOB jsonnet_file *.jsonnet )
install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/common)
