add_executable(PlottingScript PlottingScript.cxx)
add_executable(PlottingScriptComplete PlottingScriptComplete.cxx)

target_link_libraries(PlottingScript ${ROOT_LIBRARIES})
target_link_libraries(PlottingScriptComplete ${ROOT_LIBRARIES})
