
simple_plugin(CalWireT1053 "module"
			lardata_Utilities
                       # larevt_Filters
			lardata_RawData
	                lardata_RecoBase
		        SignalShapingServiceT1053_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
