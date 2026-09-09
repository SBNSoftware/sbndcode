set(
   MODULE_LIBRARIES
         sbndcode_OpDetSim
         larcore::Geometry_Geometry_service
         lardataobj::Simulation
         lardata::Utilities
         lardataobj::RawData
         sbndcode_Utilities_SignalShapingServiceSBND_service
         nurandom::RandomUtils_NuRandomService_service
         art::Framework_Core
         art::Framework_Principal
         art::Framework_Services_Optional_RandomNumberGenerator_service
         canvas::canvas
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         cetlib_except::cetlib_except
         CLHEP::CLHEP
         ${ROOT_BASIC_LIB_LIST}
)

cet_build_plugin(PMTPulseOscillation art::module SOURCE PMTPulseOscillation_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
add_subdirectory(job)
##FILE(GLOB fcl_files *.fcl)