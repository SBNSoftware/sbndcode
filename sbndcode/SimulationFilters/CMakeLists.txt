
art_make(
          MODULE_LIBRARIES
            larcore_Geometry_Geometry_service
            larcore_Geometry
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            art_Persistency_Common canvas_Persistency_Common
            art_Persistency_Provenance canvas_Persistency_Provenance
            art_Utilities canvas_Utilities
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${CLHEP}
            ${Boost_SYSTEM_LIBRARY}
            ${ROOT_BASIC_LIB_LIST}

#          MODULE_LIBRARIES larcore_Geometry
#                           larcore_Geometry_Geometry_service
#                           larsim_Simulation
#                           lardata_Utilities
#                           larevt_Filters
#                           lardataobj_RawData
#                           ${SIMULATIONBASE}
#                           ${LARRECO_LIB}
#                           ${LARDATA_LIB}
#                           ${ART_FRAMEWORK_CORE}
#                           ${ART_FRAMEWORK_PRINCIPAL}
#                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                           ${ART_PERSISTENCY_COMMON}
#                           ${ART_PERSISTENCY_PROVENANCE}
#                           ${ART_UTILITIES}
#                           ${MF_MESSAGELOGGER}
#                           ${MF_UTILITIES}
#                           ${FHICLCPP}
#                           ${CETLIB}
#                           ${ROOT_GEOM}
#                           ${ROOT_XMLIO}
#                           ${ROOT_GDML}
#                           ${ROOT_BASIC_LIB_LIST}
        )

#install_headers()
install_fhicl()
install_source()

