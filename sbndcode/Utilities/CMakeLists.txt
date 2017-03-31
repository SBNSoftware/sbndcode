cet_enable_asserts()

set( sbnd_util_lib_list lardata_Utilities_LArFFT_service
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          art_Persistency_Common canvas_Persistency_Common
                          art_Persistency_Provenance canvas_Persistency_Provenance
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                          art_Utilities canvas_Utilities
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${CETLIB}
                          ${ROOT_BASIC_LIB_LIST}
                          ${ROOT_GEOM}
    )


simple_plugin( SignalShapingServiceT1053  "service"
               ${sbnd_util_lib_list}
        )

simple_plugin( SignalShapingServiceSBND  "service"
               ${sbnd_util_lib_list}
        )
             
simple_plugin( FileCatalogMetadataSBND  "service"
               ${sbnd_util_lib_list}
        )
install_headers()
install_fhicl()
install_source()



