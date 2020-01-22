
art_make(
          MODULE_LIBRARIES
            larcore_Geometry_Geometry_service
            larcorealg_Geometry
            lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
            nusimdata_SimulationBase
            sbndcode_RecoUtils
            sbndcode_GeoWrappers
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            ${MF_MESSAGELOGGER}
            
            ${FHICLCPP}
            ${CETLIB}
            ${CLHEP}
            ${Boost_SYSTEM_LIBRARY}
            ${ROOT_BASIC_LIB_LIST}
            ${ROOT_GEOM}

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
#                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
#                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#                           ${ART_PERSISTENCY_COMMON}
#                           ${ART_PERSISTENCY_PROVENANCE}
#                           ${ART_UTILITIES}
#                           ${MF_MESSAGELOGGER}
#                           
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

