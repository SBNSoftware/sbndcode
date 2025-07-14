file(GLOB jsonnet_file *.jsonnet )
install_wp(LIST ${jsonnet_file} SUBDIRNAME pgrapher/experiment/sbnd)
add_subdirectory(DNN_ROI)