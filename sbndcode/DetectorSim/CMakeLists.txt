
simple_plugin(SimWireSBND "module"
                           larcore_Geometry_Geometry_service
                           lardataobj_Simulation
                           lardata_Utilities
                           lardataobj_RawData
                           lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                           sbndcode_Utilities_SignalShapingServiceSBND_service
                           nurandom::RandomUtils_NuRandomService_service
                           art::Framework_Core
                           art::Framework_Principal
                           art_root_io::TFileService_service
                           art::Framework_Services_Optional_RandomNumberGenerator_service
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib cetlib_except
                           CLHEP::CLHEP
                           ${ROOT_BASIC_LIB_LIST}
			   
			   

)

add_subdirectory(Services)

install_headers()
install_fhicl()
install_source()

