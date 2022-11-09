set(
   TOOL_LIBRARIES
         nurandom::RandomUtils_NuRandomService_service
         art::Utilities
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         CLHEP::CLHEP
         cetlib::cetlib
         cetlib_except::cetlib_excep
)

cet_build_plugin(PMTGainFluctuations1Dynode art::tool SOURCE PMTGainFluctuations1Dynode_tool.cc LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
