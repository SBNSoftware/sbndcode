art_make(
   LIB_LIBRARIES
         larcore_Geometry_Geometry_service
         lardataobj_RawData
         lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         nurandom_RandomUtils_NuRandomService_service
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
	 ${ART_FRAMEWORK_CORE}
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_FRAMEWORK_BASIC}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         ${ART_PERSISTENCY_COMMON}
         ${ART_PERSISTENCY_PROVENANCE}
         ${ART_UTILITIES}
         cetlib_except
         ${MF_MESSAGELOGGER}
         ${MF_UTILITIES}
         ${FHICLCPP}
         cetlib cetlib_except
         ${CLHEP}
         ${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
