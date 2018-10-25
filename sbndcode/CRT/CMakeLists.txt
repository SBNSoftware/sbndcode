include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{IFDHX_FQ_DIR}/inc )

add_subdirectory(CRTProducts)

art_make( 
    EXCLUDE
        CRTChannelMapAlg.cxx
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTSimHitProducer_module.cc
        CRTTzeroProducer_module.cc
        CRTTrackProducer_module.cc
        CRTRecoAna_module.cc
    LIBRARY_NAME
        sbndcode_CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardataobj_RecoBase
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library( LIBRARY_NAME sbndcode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
    LIBRARIES 
        larcorealg_Geometry
        sbndcode_CRTData
        cetlib_except
        ${MF_MESSAGELOGGER}
        
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${CETLIB}
)

simple_plugin( CRTGeometryHelper service
        larcorealg_Geometry
        sbndcode_CRT
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        ${CETLIB}
)

simple_plugin( CRTDetSim module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTProducts
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}
)

simple_plugin( CRTSimHitProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTProducts
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}
)

simple_plugin(CRTTrackProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTProducts
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}
)

simple_plugin(CRTTzeroProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTProducts
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}
)

simple_plugin(CRTRecoAna module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTProducts
        sbndcode_RecoUtils
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}

        larcore_Geometry_Geometry_service
        larsim_Simulation lardataobj_Simulation
        larsim_MCCheater_BackTrackerService_service
        larsim_MCCheater_ParticleInventoryService_service
        lardata_Utilities
        larevt_Filters
        lardataobj_RawData
        lardataobj_RecoBase
        lardataobj_AnalysisBase
        lardata_RecoObjects
        larpandora_LArPandoraInterface
        nusimdata_SimulationBase
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        ${FHICLCPP}
        ${ROOT_GEOM}
        ${ROOT_XMLIO}
        ${ROOT_GDML}
        ${ROOT_SPECTRUM}
)


install_headers()
install_fhicl()
install_source()

