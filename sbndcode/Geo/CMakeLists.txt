art_make(
          LIBRARY_NAME sbnd_Geometry
          LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib
                        cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
          SERVICE_LIBRARIES larcore_Geometry
                        sbnd_Geometry
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        canvas
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib_except
                        ${ROOT_CORE}
        )        

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)

