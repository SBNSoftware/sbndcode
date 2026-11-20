file(GLOB MICHEL_MODEL_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.root")
install_fw(LIST ${MICHEL_MODEL_FILES})

file(GLOB xml_file *.xml )
install_fw(LIST ${xml_file})