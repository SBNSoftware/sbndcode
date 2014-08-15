cet_enable_asserts()

set( lariat_util_lib_list LArFFT_service
                       LArProperties_service
			Geometry
			Geometry_service
			Utilities
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
    )


simple_plugin( SignalShapingServiceT1034  "service"
               ${lariat_util_lib_list}
	       BASENAME_ONLY
        )
             
 
install_headers()
install_fhicl()
install_source()



