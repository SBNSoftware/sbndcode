art_make(
#  LIB_LIBRRARIES

  MODULE_LIBRARIES
                           sbndcode_OpDetSim
                           larcore_Geometry_Geometry_service
                           lardataobj_Simulation
                           lardata_Utilities
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                           sbndcode_Utilities_SignalShapingServiceSBND_service
                           nutools_RandomUtils_NuRandomService_service
			   ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
)

# simple_plugin(SimPMTSBND "module"
#                            larcore_Geometry_Geometry_service
#                            lardataobj_Simulation
#                            lardata_Utilities
#                            lardataobj_RawData
#                            lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
#                            sbndcode_Utilities_SignalShapingServiceSBND_service
#                            nutools_RandomUtils_NuRandomService_service
# 			   ${ART_FRAMEWORK_CORE}
#                            ${ART_FRAMEWORK_PRINCIPAL}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                            canvas
#                            ${MF_MESSAGELOGGER}
#                            ${MF_UTILITIES}
#                            ${FHICLCPP}
#                            cetlib cetlib_except
#                            ${CLHEP}
#                            ${ROOT_BASIC_LIB_LIST}
# )
# 
# simple_plugin(SimArapucaSBND "module"
#                            larcore_Geometry_Geometry_service
#                            lardataobj_Simulation
#                            lardata_Utilities
#                            lardataobj_RawData
#                            lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
#                            sbndcode_Utilities_SignalShapingServiceSBND_service
#                            nutools_RandomUtils_NuRandomService_service
#                            ${ART_FRAMEWORK_CORE}
#                            ${ART_FRAMEWORK_PRINCIPAL}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                            canvas
#                            ${MF_MESSAGELOGGER}
#                            ${MF_UTILITIES}
#                            ${FHICLCPP}
#                            cetlib cetlib_except
#                            ${CLHEP}
#                            ${ROOT_BASIC_LIB_LIST}
# )

install_headers()
install_fhicl()
install_source()
