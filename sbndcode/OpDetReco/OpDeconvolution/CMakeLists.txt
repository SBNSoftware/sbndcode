art_make(
   LIB_LIBRARIES
         larcore_Geometry_Geometry_service
         lardataobj_RawData
         lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         nurandom::RandomUtils_NuRandomService_service
         art::Framework_Services_Optional_RandomNumberGenerator_service
         cetlib_except::cetlib_except
         messagefacility::MF_MessageLogger
         
         fhiclcpp::fhiclcpp
         cetlib::cetlib cetlib_except
         CLHEP::CLHEP
         ${ROOT_BASIC_LIB_LIST}
         ROOT::Core

   MODULE_LIBRARIES
         sbndcode_OpDetSim
         larcore_Geometry_Geometry_service
         lardataobj_Simulation
         lardata_Utilities
         lardataobj_RawData
         sbndcode_Utilities_SignalShapingServiceSBND_service
         nurandom::RandomUtils_NuRandomService_service
         art::Framework_Core
         art::Framework_Principal
         
         art::Framework_Services_Optional_RandomNumberGenerator_service
         canvas::canvas
         messagefacility::MF_MessageLogger
         
         fhiclcpp::fhiclcpp
         cetlib::cetlib cetlib_except
         CLHEP::CLHEP
         ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
##FILE(GLOB fcl_files *.fcl)

add_subdirectory(Alg)
add_subdirectory(reco_deco_sbnd)
