cet_enable_asserts()

set( sbnd_util_lib_list lardata_Utilities_LArFFT_service
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          lardataalg_DetectorInfo
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          art_Persistency_Common
                          art_Persistency_Provenance 
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          art_Framework_Services_System_TriggerNamesService_service
                          art_Utilities
                          canvas
                          ${MF_MESSAGELOGGER}
                          cetlib cetlib_except
                          ${ROOT_BASIC_LIB_LIST}
                          ${ROOT_GEOM}
    )


simple_plugin( SignalShapingServiceT1053  "service"
               ${sbnd_util_lib_list}
        )

simple_plugin( SignalShapingServiceSBND  "service"
               ${sbnd_util_lib_list}
        )

install_headers()
install_fhicl()
install_source()
