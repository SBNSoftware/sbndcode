include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

art_make(
    EXCLUDE
        CRTChannelMapAlg.cxx
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
    LIBRARY_NAME
        sbndcode_CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library( LIBRARY_NAME sbndcode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
    LIBRARIES larcorealg_Geometry
        sbndcode_CRTData
        cetlib_except
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
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
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CETLIB}
)

simple_plugin( CRTDetSim module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
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

install_headers()
install_fhicl()
install_source()

