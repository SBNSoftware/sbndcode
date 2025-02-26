
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

cet_build_plugin ( TPCRecomAna art::module
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  lardataobj::RawData
  lardataobj::RecoBase
  ROOT::Core
  ROOT::Tree
  art::Framework_Core
  art::Framework_Principal
  art::Utilities
  art_root_io::tfile_support
  art_root_io::TFileService_service
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger

  BASENAME_ONLY
              )

add_subdirectory(Services)

install_headers()
install_fhicl()
install_source()

