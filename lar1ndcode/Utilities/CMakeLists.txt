cet_enable_asserts()

set( lar1nd_util_lib_list lardata_Utilities_LArFFT_service
                          lardata_DetectorInfoServices_LArPropertiesServiceStandard_service
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_PERSISTENCY_COMMON}
                          ${ART_PERSISTENCY_PROVENANCE}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                          ${ART_UTILITIES}
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${CETLIB}
                          ${ROOT_BASIC_LIB_LIST}
    )


simple_plugin( SignalShapingServiceT1053  "service"
               ${lar1nd_util_lib_list}
               BASENAME_ONLY
        )
             
 
simple_plugin( FileCatalogMetadataLAr1ND  "service"
               ${lar1nd_util_lib_list}
               BASENAME_ONLY
        )
install_headers()
install_fhicl()
install_source()



