# simple_plugin( LArIATGeometryHelper service
#                   larcore_Geometry
#                 ${ART_FRAMEWORK_CORE}
#                 ${ART_FRAMEWORK_IO_SOURCES}
#                 ${ART_FRAMEWORK_PRINCIPAL}
#                 art_Persistency_Common canvas_Persistency_Common
#                 art_Persistency_Provenance canvas_Persistency_Provenance
#                 art_Utilities canvas_Utilities
#                 ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                 ${MF_MESSAGELOGGER}
#                 ${MF_UTILITIES}
#                 ${FHICLCPP}
#                 ${CETLIB}
#                 ${ROOT_BASIC_LIB_LIST}                
# )



art_make(
          LIBRARY_NAME sbnd_Geometry
          LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
          SERVICE_LIBRARIES larcore_Geometry
                          sbnd_Geometry
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
        )        


install_headers()
install_fhicl()
install_source()

# testing
add_subdirectory(gdml)
