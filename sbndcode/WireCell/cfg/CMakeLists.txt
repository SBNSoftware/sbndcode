add_subdirectory(pgrapher)

file(GLOB jsonnet_file *.jsonnet)
install_wp(LIST ${jsonnet_file})
