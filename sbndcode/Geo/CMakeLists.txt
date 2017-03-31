art_make(
          LIBRARY_NAME sbnd_Geometry
          LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
          SERVICE_LIBRARIES larcore_Geometry
                        sbnd_Geometry
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        canvas_Utilities
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_CORE}
        )        

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)

