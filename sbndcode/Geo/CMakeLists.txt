# simple_plugin( LArIATGeometryHelper service
#                   larcore_Geometry
#                 ${ART_FRAMEWORK_CORE}
#                 ${ART_FRAMEWORK_IO_SOURCES}
#                 ${ART_FRAMEWORK_PRINCIPAL}
#                 art_Persistency_Common canvas
#                 art_Persistency_Provenance canvas
#                 art_Utilities canvas
#                 ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                 ${MF_MESSAGELOGGER}
#                 ${MF_UTILITIES}
#                 ${FHICLCPP}
#                 cetlib cetlib_except
#                 ${ROOT_BASIC_LIB_LIST}                
# )



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

# testing
add_subdirectory(gdml)
