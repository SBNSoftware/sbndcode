
simple_plugin(SimWireT1053 "module"
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation
                           lardata_Utilities
                           lardata_RawData
                           lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                           lar1ndcode_Utilities_SignalShapingServiceT1053_service
                           larsim_RandomUtils_LArSeedService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_PERSISTENCY_PROVENANCE}
                           ${ART_UTILITIES}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
