simple_plugin( LArIATGeometryHelper service
	  	Geometry
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_IO_SOURCES}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_PERSISTENCY_COMMON}
		${ART_PERSISTENCY_PROVENANCE}
		${ART_UTILITIES}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
	        ${FHICLCPP}
		${CETLIB}
		${ROOT_BASIC_LIB_LIST}                
		BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

# testing
add_subdirectory(gdml)
