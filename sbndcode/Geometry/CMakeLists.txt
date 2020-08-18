art_make(
          LIB_LIBRARIES larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib
                        cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
          SERVICE_LIBRARIES
                        sbndcode_Geometry
                        larcorealg_Geometry
                        art_Framework_Services_Registry
                        canvas
                        ${FHICLCPP}
                        cetlib_except
                        ${ROOT_CORE}
        )


install_headers()
install_fhicl()
install_source()

# GDML detector descriptions
add_subdirectory(gdml)
# Wrappers for easier access to geometry properties
add_subdirectory(GeometryWrappers)
