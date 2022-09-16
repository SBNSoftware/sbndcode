art_make(MODULE_LIBRARIES
	larsim_Utils
        nusimdata_SimulationBase
	lardataobj_RecoBase
        nusimdata::SimulationBase
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILE_SUPPORT} 
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}

        ${ROOT_BASIC_LIB_LIST}

	sbndcode_GeoWrappers
	sbndcode_RecoUtils
	)

install_headers()
install_fhicl()
install_source()
