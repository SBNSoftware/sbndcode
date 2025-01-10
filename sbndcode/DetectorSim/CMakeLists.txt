
cet_build_plugin( SimWireSBND
                           art::module LIBRARIES
                           larcore::Geometry_Geometry_service
                           lardataobj::Simulation
                           lardata::Utilities
                           lardataobj::RawData
                           lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                           larevt::CalibrationDBI_Providers
                           sbndcode_Utilities_SignalShapingServiceSBND_service
                           nurandom::RandomUtils_NuRandomService_service
                           art::Framework_Core
                           art::Framework_Principal
                           art_root_io::TFileService_service
                           art::Framework_Services_Optional_RandomNumberGenerator_service
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           CLHEP::CLHEP
                           ROOT::Core

)

add_subdirectory(Services)

install_headers()
install_fhicl()
install_source()

