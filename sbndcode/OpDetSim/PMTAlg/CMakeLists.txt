set(
   TOOL_LIBRARIES
         nurandom::RandomUtils_NuRandomService_service
         CLHEP::CLHEP
)

cet_build_plugin(PMTGainFluctuations1Dynode art::tool SOURCE PMTGainFluctuations1Dynode_tool.cc LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
