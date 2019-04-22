simple_plugin(SpaceChargeServiceSBND "service"
              sbndcode_SpaceCharge
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Persistency_Common canvas
	      art_Persistency_Provenance canvas
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
	      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art_Utilities canvas
	      ${MF_MESSAGELOGGER}
	      
  	      cetlib cetlib_except
   	      ${ROOT_BASIC_LIB_LIST}
)
install_headers()
install_fhicl()
install_source()
