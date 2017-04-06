art_make(
          LIBRARY_NAME sbnd_Geometry
          LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib
                        cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
          SERVICE_LIBRARIES
                        sbnd_Geometry
                        larcore_Geometry
                        art_Framework_Services_Registry
                        canvas
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib_except
                        ${ROOT_CORE}
        )


install_headers()
install_fhicl()
install_source()

# GDML detector descriptions
add_subdirectory(gdml)
