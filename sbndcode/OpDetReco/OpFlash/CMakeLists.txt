include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
art_make( 
    BASENAME_ONLY
    LIBRARY_NAME   sbndcode_OpDetReco_OpFlash
	LIB_LIBRARIES
        sbnd_Geometry
        sbndcode_OpDetSim
        larcore_Geometry_Geometry_service
        lardataobj_RecoBase
        larsim_Simulation 
        lardataobj_Simulation
        nusimdata_SimulationBase
        ${LARDATA_LIB}
        canvas
        cetlib
        cetlib_except
        ${ART_FRAMEWORK_ART}
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_UTILITIES}
        ${ART_FRAMEWORK_BASIC}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        art_Persistency_Common
        art_Persistency_Provenance
        art_Utilities
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_GEOM}
        ${ROOT_XMLIO}
        ${ROOT_GDML}
        ${ROOT_BASIC_LIB_LIST}
    MODULE_LIBRARIES
        sbndcode_OpDetReco_OpFlash_FlashFinder
        sbndcode_OpDetSim
        sbnd_Geometry
        larcore_Geometry_Geometry_service
        lardataobj_RecoBase
        larsim_Simulation 
        lardataobj_Simulation
        nusimdata_SimulationBase
        ${LARDATA_LIB}
        canvas
        cetlib
        cetlib_except
        ${ART_FRAMEWORK_ART}
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_BASIC}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        art_Persistency_Common
        art_Persistency_Provenance
        art_Utilities
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_GEOM}
        ${ROOT_XMLIO}
        ${ROOT_GDML}
        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(FlashFinder)
add_subdirectory(job)

