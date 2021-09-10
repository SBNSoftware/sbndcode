art_make(
   TOOL_LIBRARIES
         art_Utilities
         larcore_Geometry_Geometry_service
         lardataobj_RawData
         lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         ${MF_MESSAGELOGGER}
         ${FHICLCPP}
         ${CETLIB}
         ${CLHEP}
         ${ROOT_BASIC_LIB_LIST}
         ${ART_FRAMEWORK_CORE}
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_FRAMEWORK_BASIC}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         ${ART_PERSISTENCY_COMMON}
         ${ART_PERSISTENCY_PROVENANCE}
         ${ART_UTILITIES}
         cetlib
         cetlib_except
         ${MF_UTILITIES}

)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
