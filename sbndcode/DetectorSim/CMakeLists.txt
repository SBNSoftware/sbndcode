
simple_plugin(SimWireT1053 "module"
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           lardata_Utilities
                           lardataobj_RawData
                           lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                           sbndcode_Utilities_SignalShapingServiceT1053_service
                           larsim_RandomUtils_LArSeedService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities canvas_Utilities
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
