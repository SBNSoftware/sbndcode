art_make(
   TOOL_LIBRARIES
         art::Utilities
         larcore_Geometry_Geometry_service
         lardataobj_RawData
         lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         CLHEP::CLHEP
         ${ROOT_BASIC_LIB_LIST}
         art::Framework_Core
         art::Framework_Principal
         ${ART_FRAMEWORK_BASIC}
         art::Framework_Services_Registry
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         
         art::Persistency_Common
         art::Persistency_Provenance
         art::Utilities
         cetlib::cetlib
         cetlib_except::cetlib_except
         

)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
